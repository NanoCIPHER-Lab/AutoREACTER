"""
AutoREACTER Command-Line Interface (CLI) module.

This module provides the main entry point for the AutoREACTER pipeline, orchestrating
the end-to-end workflow: functional group detection, reaction identification,
non-reactant selection, and simulation preparation (force field generation,
REACTER file building, and LAMMPS setup).

Classes:
    ErrorHandler: Tracks the completion state of each pipeline stage via a waterfall model.
    ARXCLI: The main CLI controller that ties together all pipeline components.
"""

from contextlib import contextmanager
import os
from pathlib import Path
import sys
import threading
from PIL import Image
from typing import Optional

from AutoREACTER.session import read_input
from AutoREACTER.input_parser import InputParser
from AutoREACTER.detectors.functional_groups_detector import FunctionalGroupsDetector
from AutoREACTER.detectors.reaction_detector import ReactionDetector
from AutoREACTER.detectors.non_monomer_detector import NonReactantsDetector
from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import PrepareReactions
from AutoREACTER.reaction_preparation.ff_wrapper.molecule_3d_preparation import Molecule3DPreparation
from AutoREACTER.reaction_preparation.ff_wrapper.ff_wrapper import FFWrapper
from AutoREACTER.reaction_preparation.ff_wrapper.REACTER_files_builder import REACTERFilesBuilder
from AutoREACTER.sim_setup.simulation_setup import SimulationSetupManager


class ErrorHandler:
    """
    Tracks pipeline stage completion using a waterfall model.

    The waterfall enforces a strict ordering: reactions must be selected
    before non-reactants, and both must be completed before processing.
    Each stage is represented as a boolean flag, initially set to False,
    and flipped to True once the corresponding step finishes successfully.
    """

    def waterfall_order(self):
        """
        Return the initial (all-false) waterfall state dictionary.

        Returns
        -------
        dict
            A dictionary with keys 'select_reactions', 'select_non_reactants',
            and 'process', each mapped to False.
        """
        order = {
            "select_reactions": False,
            "select_non_reactants": False,
            "process": False,
        }
        return order


class ARXCLI:
    """
    Main CLI controller for the AutoREACTER pipeline.

    This class accepts an input file path, initialises the session, and exposes
    methods to step through — or run end-to-end — the reaction-detection and
    simulation-preparation workflow.  Intermediate images (molecules, functional
    groups, reactions, non-reactants) are saved automatically to the session's
    image directory.

    Parameters
    ----------
    input : Path
        Path to the input configuration file consumed by ``read_input``.

    Attributes
    ----------
    session : Session
        The session object holding all parsed data and state.
    img_dir : Path
        Directory where debug/visualisation images are written.
    """

    # Shared error-handler instance across all ARXCLI objects
    EH = ErrorHandler()

    def __init__(self, input: Path) -> None:
        # Initialise the waterfall tracker (all stages False)
        self.error_handler = self.EH.waterfall_order()

        self.input = input
        abs_path = self.input.resolve()
        print(f"[OK] Read input from {abs_path}")

        # Parse the input file and create the session
        self.session = read_input(abs_path)
        self.img_dir = self.session.images_dir
        # with open(self.session.output_dir / "AutoREACTER.log", 'w') as f:
        #     f.write("--- Starting AutoREACTER Session ---\n")

        # Save an initial grid image of all monomers
        self._save_rdkit_img(
            InputParser().initial_molecules_image_grid(self.session),
            self.img_dir / "monomers.png",
        )

        # --- Pipeline stage flags ---
        self._fg_detected = False           # Functional groups have been detected
        self._reactions_detected = False    # Reactions have been detected
        self._reactions_selected = False    # User (or auto) has picked reactions
        self._non_reactants_detected = False
        self._non_reactants_selected = False

        # Bootstrap: detect functional groups and reactions right away
        self._ensure_reactions_detected()

    # ------------------------------------------------------------------
    # Public visualisation / inspection helpers
    # ------------------------------------------------------------------

    

    def show_molecules(self) -> Image:
        """
        Return a PIL/rdkit image grid of the initial molecules (monomers).

        Returns
        -------
        Image
            Grid image of all starting molecules.
        """
        return InputParser().initial_molecules_image_grid(self.session)

    def show_functional_groups(self) -> Image:
        """
        Return an image grid with functional groups highlighted on each molecule.

        Triggers functional-group detection if it hasn't run yet.

        Returns
        -------
        Image
            Grid image with functional groups annotated.
        """
        if not self._fg_detected:
            self._ensure_fg_detected()
        return FunctionalGroupsDetector().functional_group_highlighted_molecules_image_grid(self.session)

    def show_reactions(self) -> Image:
        """
        Return an image grid showing the detected reactions.

        Triggers reaction detection (and therefore FG detection) if needed.

        Returns
        -------
        Image
            Grid image of available reactions.
        """
        if not self._reactions_detected:
            self._ensure_reactions_detected()
        return ReactionDetector().available_reaction_image_grid(self.session)
    

    # ------------------------------------------------------------------
    # Selection steps
    # ------------------------------------------------------------------

    def select_reactions(self) -> None:
        """
        Interactively (or automatically) select which reactions to proceed with.

        If more than one reaction is found the user is prompted; otherwise
        selection is performed automatically.  Marks the 'select_reactions'
        waterfall stage as complete.
        """
        if not self._reactions_detected:
            self._ensure_reactions_detected()

        if not self._reactions_selected:
            if not self.session.reaction_instances:
                raise RuntimeError("No reactions detected. Cannot proceed.")
            else:
                selector = ReactionDetector()
                selector.reaction_selection(self.session)
                self._reactions_selected = True
                self.error_handler["select_reactions"] = True

    def show_non_reactants(self) ->  Optional['Image']:
        """
        Return an image visualising detected non-reactant species.

        This is the only place where the non-reactant visualization image is built
        and saved.
        """
        self._ensure_non_reactants_detected()

        detector = NonReactantsDetector()
        img = detector.non_reactants_to_visualization(self.session)

        if img is None:
            return None

        self._save_rdkit_img(
            img,
            self.img_dir / "non_reactants.png",
            is_non_reactant=True,
        )

        return img

    def select_non_reactants(self) -> None:
        """
        Interactively or automatically select non-reactant species.

        If non-reactants are found, the user is prompted to keep/discard/select them.
        If no non-reactants are found, the step is skipped.

        In all successful cases, mark the non-reactant selection stage complete.
        """
        self._ensure_non_reactants_detected()

        if self._non_reactants_selected:
            return

        if self.session.non_reactants and len(self.session.non_reactants) > 0:
            NonReactantsDetector().non_reactant_selection(self.session)

        self._non_reactants_selected = True
        self.error_handler["select_non_reactants"] = True


    def prepare_reactions(self) -> None:
        """
        Run the reaction template preparation stage.

        This includes template generation, edge/initiator detection, and
        template modification.  Saves visualisation images of the templates
        to disk.  Marks the 'process' waterfall stage as complete.
        """
        PrepareReactions(self.session).prepare_reactions(self.session)
        self.error_handler["process"] = True
        highlight_types = ["template", "edge", "initiators", "delete"]
        for highlight_type in highlight_types:
            img = PrepareReactions(self.session).reaction_templates_highlighted_image_grid(
                self.session, highlight_type=highlight_type
            )
            self._save_rdkit_img(
                img, self.img_dir / f"templates_{highlight_type}.png"
            )
        return None

    def show_reaction_templates(self, highlight_type: str = "template") -> Image:
        """
        Return an image grid visualising the reaction templates.

        Parameters
        ----------
        highlight_type : str
            The type of template highlighting to show. Options include:
            'template' (full template), 'edge' (reaction edge), 'initiators',
            and 'delete' (atoms to delete).

        Returns
        -------
        Image
            Grid image of reaction templates with specified highlights.
        """
        return PrepareReactions(self.session).reaction_templates_highlighted_image_grid(
            self.session, highlight_type=highlight_type
        )

    # ------------------------------------------------------------------
    # Full pipeline processing
    # ------------------------------------------------------------------

    def process(self):
        """
        Execute the back-half of the pipeline in one shot.

        This method runs reaction template preparation, 3D geometry setup,
        force-field generation via the LUNAR API, REACTER file building,
        and LAMMPS simulation writing — in that order.

        Raises
        ------
        RuntimeError
            If the waterfall stages 'select_reactions', 'select_non_reactants',
            or 'process' have not been completed beforehand.
        """
        # --- Guard: enforce waterfall ordering ---
        if not self.error_handler["select_reactions"]:
            raise RuntimeError("Reactions have not been selected.")

        if not self.error_handler["select_non_reactants"]:
            raise RuntimeError("Non-reactants have not been selected.")


        if not self.error_handler["process"]:
            raise RuntimeError("Processing has not been completed. Because not all selections were made.")
        
        # Run the entire processing pipeline, redirecting all print output to a log file in the output directory.
        with self._writer():
            # --- Run each stage in sequence ---
            Molecule3DPreparation(self.session).prepare_molecule_3d_geometry(self.session)

            FFWrapper(self.session).generate_force_field_files(self.session)

            REACTERFilesBuilder(self.session).molecule_template_preparation(self.session)

            SimulationSetupManager().setup_and_write_simulation(self.session)

    # ------------------------------------------------------------------
    # Internal helpers – lazy detection & image saving
    # ------------------------------------------------------------------

    def _ensure_fg_detected(self):
        """
        Run functional-group detection exactly once (idempotent).

        Sets ``self._fg_detected = True`` after the first call so subsequent
        invocations are no-ops.
        """
        if not self._fg_detected:
            FunctionalGroupsDetector().functional_groups_detector(self.session)
            self._fg_detected = True

    def _ensure_reactions_detected(self):
        """
        Ensure functional groups *and* reactions have been detected.

        This is the bootstrapping step: it first guarantees FG detection,
        saves a functional-group image, then runs reaction detection and
        saves a reaction image.  Both stages are idempotent.
        """
        self._ensure_fg_detected()

        # Save the functional-group visualisation for later inspection
        self._save_rdkit_img(
            FunctionalGroupsDetector().functional_group_highlighted_molecules_image_grid(self.session),
            self.img_dir / "functional_groups.png",
        )

        if not self._reactions_detected:
            ReactionDetector().reaction_detector(self.session)
            # Save the reaction visualisation
            self._save_rdkit_img(
                ReactionDetector().available_reaction_image_grid(self.session),
                self.img_dir / "reactions.png",
            )
            self._reactions_detected = True

    def _ensure_non_reactants_detected(self):
        """
        Ensure non-reactant species have been detected.

        This should only detect non-reactants. It should not generate or print
        visualization output, because show_non_reactants() handles that.
        """
        if not self._reactions_selected:
            self.select_reactions()

        if not self._non_reactants_detected:
            NonReactantsDetector().non_monomer_detector(self.session)
            self._non_reactants_detected = True

    def _save_rdkit_img(self, img, path: Path, is_non_reactant: bool = False):
        """
        Save an RDKit/PIL/IPython image to disk.

        Parameters
        ----------
        img
            Image returned by RDKit or IPython display.
        path : pathlib.Path
            Destination file path.
        is_non_reactant : bool
            If True, allow None images because no non-reactants may exist.
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)

        if img is None:
            if is_non_reactant:
                return
            raise ValueError("No image was generated. Cannot save molecule image.")

        # Case 1: PIL image
        if hasattr(img, "save"):
            img.save(path)
            return

        # Case 2: raw PNG/SVG bytes
        if isinstance(img, (bytes, bytearray)):
            path.write_bytes(img)
            return

        # Case 3: IPython.display.Image usually stores image bytes in .data
        data = getattr(img, "data", None)

        if isinstance(data, (bytes, bytearray)):
            path.write_bytes(data)
            return

        if isinstance(data, str):
            path.write_text(data, encoding="utf-8")
            return

        # Case 4: SVG string
        if isinstance(img, str):
            path.write_text(img, encoding="utf-8")
            return

        raise TypeError(
            f"Unsupported image type: {type(img)}. "
            "Expected PIL image, bytes, SVG string, or IPython display image."
        )
        
    # ------------------------------------------------------------------
    # Output management helpers
    # ------------------------------------------------------------------

    @contextmanager
    def _writer(self, filename="AutoREACTER.log"):
        """
        Helper function to tee standard output to both the terminal AND a text file.
        Catches Python prints AND system-level subprocess output.
        """
        self.session.output_dir.mkdir(parents=True, exist_ok=True)
        log_path = self.session.output_dir / filename

        # 1. Save the original OS-level terminal output
        original_stdout_fd = os.dup(1)

        # 2. Create an OS-level pipe (a temporary tunnel for our data)
        pipe_read_fd, pipe_write_fd = os.pipe()

        def tee_thread():
            """Background worker that reads the pipe and writes to both destinations."""
            with open(log_path, 'a') as log_file:
                while True:
                    # Read incoming data from the pipe
                    data = os.read(pipe_read_fd, 1024)
                    
                    # If the pipe is closed, stop the thread
                    if not data:
                        break
                        
                    # Write to the actual terminal
                    os.write(original_stdout_fd, data)
                    
                    # Write to the log file
                    log_file.write(data.decode('utf-8', errors='replace'))
                    log_file.flush()

        # 3. Start the background thread
        thread = threading.Thread(target=tee_thread)
        thread.start()

        # Flush Python's buffers before we switch the tracks
        sys.stdout.flush()

        try:
            # 4. Redirect all OS-level output to the write-end of our pipe
            os.dup2(pipe_write_fd, 1)
            yield
            
        finally:
            # Flush Python buffers one last time
            sys.stdout.flush()
            
            # 5. Restore the original terminal output
            os.dup2(original_stdout_fd, 1)
            
            # 6. Close the write end of the pipe (this tells the thread to stop)
            os.close(pipe_write_fd)
            
            # 7. Wait for the thread to finish processing the last bits of data
            thread.join()
            
            # 8. Clean up remaining file descriptors
            os.close(pipe_read_fd)
            os.close(original_stdout_fd)

    # ------------------------------------------------------------------
    # Magic Methods
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        """
        Provides a clean string representation for interactive environments (like Jupyter).
        Prevents the default '<...ARXCLI object at 0x...>' memory address from printing.
        """
        # Safely grab the filename if it exists
        if hasattr(self, 'input'):
            return ""
        return ""
