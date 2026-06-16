from pathlib import Path


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



class ARXCLI:
    def __init__(self, input: Path):
        self.input = input
        abs_path = self.input.resolve()
        print(f"[OK] Read input from {abs_path}")
        self.session = read_input(abs_path)
        self.img_dir = self.session.images_dir
        self._save_rdkit_img(
            InputParser().initial_molecules_image_grid(self.session),
            self.img_dir / "monomers.png"
        )
        # State the tracking of the session 
        self._fg_detected = False
        self._reactions_detected = False
        self._reactions_selected = False
        self._non_reactants_detected = False
        self._non_reactants_selected = False

        self._ensure_reactions_detected()


        
    def show_molecules(self):
        return InputParser().initial_molecules_image_grid(self.session)

    def show_functional_groups(self):
        self._ensure_fg_detected()
        return FunctionalGroupsDetector().functional_group_highlighted_molecules_image_grid(self.session)

    def show_reactions(self):
        self._ensure_reactions_detected()
        return ReactionDetector().available_reaction_image_grid(self.session)
    
    def select_reactions(self):
        self._ensure_reactions_detected()
        if not self._reactions_selected:
            # Only prompt if there is more than 1 reaction
            if self.session.reaction_instances and len(self.session.reaction_instances) > 1:
                print("[INFO] Multiple reactions found. Prompting for selection...")
                ReactionDetector().reaction_selection(self.session)
            else:
                print("[INFO] 1 or 0 reactions found. Auto-selecting...")
            
            self._reactions_selected = True
    
    def show_non_reactants(self):
        self._ensure_non_reactants_detected()
        return NonReactantsDetector().non_reactants_to_visualization(self.session)
    
    def select_non_reactants(self):
        self._ensure_non_reactants_detected()
        if not self._non_reactants_selected:
            if self.session.non_reactants and len(self.session.non_reactants) > 0:
                print("[INFO] Non-reactants found. Prompting for selection...")
                NonReactantsDetector().non_reactant_selection(self.session)
            else:
                print("[INFO] No non-reactants found. Skipping selection...")
            
            self._non_reactants_selected = True
    
    def process(self):
        """Runs the back-half of the pipeline all at once."""
        # Ensure the waterfall has reached the bottom before processing
        self.select_non_reactants()

        print("[INFO] Preparing reaction templates...")
        PrepareReactions(self.session).prepare_reactions(self.session)

        print("[INFO] Preparing 3D Geometries...")
        Molecule3DPreparation(self.session).prepare_molecule_3d_geometry(self.session)

        print("[INFO] Running LUNAR API force field generation...")
        FFWrapper(self.session).generate_force_field_files(self.session)

        print("[INFO] Building REACTER files...")
        REACTERFilesBuilder(self.session).molecule_template_preparation(self.session)

        print("[INFO] Writing LAMMPS simulation setup...")
        SimulationSetupManager().setup_and_write_simulation(
            self.session
        )

    def _ensure_fg_detected(self):
        if not self._fg_detected:
            FunctionalGroupsDetector().functional_groups_detector(self.session)
            self._fg_detected = True

    def _ensure_reactions_detected(self):
        self._ensure_fg_detected()
        # Save functional groups image
        self._save_rdkit_img(
            FunctionalGroupsDetector().functional_group_highlighted_molecules_image_grid(self.session),
            self.img_dir / "functional_groups.png"
        )
        if not self._reactions_detected:
            ReactionDetector().reaction_detector(self.session)
            # Save reactions image
            self._save_rdkit_img(
                ReactionDetector().reaction_highlighted_molecules_image_grid(self.session),
                self.img_dir / "reactions.png"
            )
            self._reactions_detected = True

    def _ensure_non_reactants_detected(self):
        if not self._reactions_selected:
            self.select_reactions()

        if not self._non_reactants_detected:
            print("[INFO] Detecting non-reactants...")
            NonReactantsDetector().non_monomer_detector(self.session)
            self._non_reactants_detected = True
            
        self._save_rdkit_img(
            NonReactantsDetector().non_reactants_to_visualization(self.session),
            self.img_dir / "non_reactants.png"
        )

    def _save_rdkit_img(self, img, path: Path):
        img.save(path)