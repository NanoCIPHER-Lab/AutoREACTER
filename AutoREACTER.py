import sys
import os
import time
from PIL import Image

from AutoREACTER.session import read_input
from AutoREACTER.cache import RunDirectoryManager
from AutoREACTER.detectors.functional_groups_detector import FunctionalGroupsDetector
from AutoREACTER.detectors.reaction_detector import ReactionDetector
from AutoREACTER.detectors.non_monomer_detector import NonReactantsDetector
from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import PrepareReactions
from AutoREACTER.reaction_preparation.lunar_client.molecule_3d_preparation import Molecule3DPreparation
from AutoREACTER.reaction_preparation.lunar_client.lunar_api_wrapper import LunarAPIWrapper
from AutoREACTER.reaction_preparation.lunar_client.REACTER_files_builder import REACTERFilesBuilder
from AutoREACTER.sim_setup.simulation_setup import SimulationSetupManager




def loading_message(message: str, duration: float = 3.0, interval: float = 0.5) -> None:
    """
    Print a loading message with animated dots to indicate progress.
    
    Args:
        message: The message to display
        duration: Total duration of the animation in seconds
        interval: Time between each dot in seconds
    """
    print(f"[INFO] {message}", end="", flush=True)
    steps = int(duration / interval)

    for _ in range(steps):
        print(".", end="", flush=True)
        time.sleep(interval)

    print()  # newline



def save_image(img: Image.Image, path: str, label: str = "Image") -> None:
    """
    Save a PIL Image to disk with robust error handling and user feedback.
    
    Args:
        img: PIL Image object to save
        path: Full path where the image should be saved
        label: Human-readable label used in log messages
    """
    if img is None:
        print(f"[WARN] {label}: Image is None, skipping save.")
        return
    try:
        if hasattr(img, "save"):
            img.save(path)
            print(f"\n[OK] {label} saved → {path}")
        else:
            print(f"[ERROR] {label}: Object has no .save() method")
    except Exception as e:
        print(f"[ERROR] Failed to save {label} to {path}: {e}")

def help_message() -> None:
    """
    Print usage instructions for the AutoREACTER script.
    This function is called when the user runs the script without arguments or with incorrect options.
    """
    print("Usage:")
    print("  python AutoREACTER.py -i <input.json>")
    print("\nOptions:")
    print("  -i, --input        Path to input JSON file")

def AutoREACTER(input_file: str) -> None:
    """
    Main function to run the AutoREACTER workflow.
    This function orchestrates the entire process from reading inputs, detecting reactions, preparing files, and setting up the simulation.
     Args:
        input_file: Path to the input JSON file containing simulation parameters and monomer information.
    
    Workflow Steps:
    1. Initialize Session: Read and validate inputs, set up staging and output directories.
    2. Functional Group Detection: Identify functional groups in the monomers and generate visualizations
    3. Reaction Discovery and Selection: Detect possible reactions based on functional groups and allow user selection.
    4. Non-monomer (Additive) Detection: Identify any non-reactant additives and allow user selection.
    5. Reaction Template Preparation: Prepare reaction templates and generate visualizations for review.
    6. 3D Geometry Preparation: Generate 3D geometries for molecules and reactions.
    7. Lunar API Processing: Send data to the Lunar API and retrieve results.
    8. Build REACTER Input Files: Create the necessary input files for REACTER based on the processed data.
    9. Final Simulation Setup and Output: Organize all outputs into the final directory structure and provide user feedback.
    """
    # === 1. Initialize Session ===
    session = read_input(input_file)

    # Generate initial visualization
    try:
        from AutoREACTER.input_parser import InputParser # Imported here just for the image generator
        img = InputParser().initial_molecules_image_grid(session.inputs)
        save_image(img, os.path.join(session.images_dir, "monomers.png"), "Monomers Grid")
    except Exception:
        print("[WARN] Failed to generate initial molecules image grid")

    # === 2. Functional Group Detection ===
    fg_detector = FunctionalGroupsDetector()
    functional_groups, functional_groups_imgs = fg_detector.functional_groups_detector(
        session.inputs.monomers
    )
    try:
        img = fg_detector.functional_group_highlighted_molecules_image_grid(functional_groups_imgs)
        save_image(img, os.path.join(session.images_dir, "functional_groups.png"), "Functional Groups")
    except Exception:
        pass  

    # === 3. Reaction Discovery and Selection ===
    reaction_detector = ReactionDetector()
    reaction_instances = reaction_detector.reaction_detector(functional_groups)
    try:
        img = reaction_detector.available_reaction_image_grid(reaction_instances)
        save_image(img, os.path.join(session.images_dir, "reactions.png"), "Available Reactions")
    except Exception:
        pass

    selected_reactions = reaction_detector.reaction_selection(reaction_instances)

    # === 4. Non-monomer (Additive) Detection ===
    non_monomer_detector = NonReactantsDetector()
    non_reactants = non_monomer_detector.non_monomer_detector(
        session.inputs, selected_reactions
    )
    try:
        img = non_monomer_detector.non_reactants_to_visualization(non_reactants)
        save_image(img, os.path.join(session.images_dir, "non_reactants.png"), "Non-Reactants")
    except Exception:
        pass

    updated_inputs = non_monomer_detector.non_reactant_selection(
        session.inputs, non_reactants
    )

    # === 5. Reaction Template Preparation ===
    prepare_reactions = PrepareReactions(session)
    prepared_reactions = prepare_reactions.prepare_reactions(selected_reactions)

    try:
        for highlight_type, filename in [
            ("template", "templates_all.png"),
            ("edge", "templates_edge.png"),
            ("initiators", "templates_initiators.png"),
            ("delete", "templates_delete.png")
        ]:
            img = prepare_reactions.reaction_templates_highlighted_image_grid(
                prepared_reactions, highlight_type=highlight_type
            )
            save_image(img, os.path.join(session.images_dir, filename))
    except Exception:
        print("[WARN] Failed to generate one or more reaction template visualizations")

    print(
        f"\n[INFO] Reaction preparation completed.\n"
        f"All generated images have been saved in the {session.images_dir} directory.\n"
        f"Please review the reaction templates before proceeding.\n"
    )

    # Interactive user confirmation
    ok_pass = input("Type 'ok' to continue: ").strip().lower()
    if ok_pass != "ok":
        print("[EXIT] Workflow stopped by user.")
        sys.exit(0)

    # === 6. 3D Geometry Preparation ===
    molecule_3d = Molecule3DPreparation(session)
    updated_inputs_3d, prepared_reactions_3d = molecule_3d.prepare_molecule_3d_geometry(
        updated_inputs, prepared_reactions
    )

    # === 7. Lunar API Processing ===
    lunar_api = LunarAPIWrapper(session)
    lunar_results = lunar_api.lunar_workflow(
        updated_inputs_3d, prepared_reactions_3d
    )

    print("\n")
    loading_message("Lunar API workflow completed. Proceeding to build REACTER files")
    time.sleep(1.5)
    print("\n")

    # === 8. Build REACTER Input Files ===
    builder = REACTERFilesBuilder(
        session=session,
        updated_inputs_with_3d_mols=updated_inputs_3d
    )

    reacter_files = builder.molecule_template_preparation(
        lunar_files=lunar_results,
        prepared_reactions_with_3d_mols=prepared_reactions_3d
    )

    # Move generated files to final output directory using RunDirectoryManager
    run_manager = RunDirectoryManager(session.output_dir.parent)
    reacter_files = run_manager.move_reacter_files(
        reacter_files,
        staging_dir=session.staging_dir,
        final_dir=session.output_dir
    )

    # === 9. Final simulation setup and output ===
    Simulation_setup_manager = SimulationSetupManager()
    updated_inputs_3d = Simulation_setup_manager.setup_and_write_simulation(
        setup=updated_inputs_3d,
        reacter_files=reacter_files,
        run_dir=session.output_dir
    )

    print("\n[INFO] AutoREACTER workflow completed successfully.\n")
    print(f"Final REACTER files are located in: {session.output_dir}")


if __name__ == "__main__":
    args = sys.argv[1:]

    if not args:
        help_message()
        sys.exit(1)

    if "-i" in args or "--input" in args:
        idx = args.index("-i") if "-i" in args else args.index("--input")
        if idx + 1 >= len(args):
            print("Error: Missing input file.")
            sys.exit(1)

        input_file = args[idx + 1]
        if not os.path.isfile(input_file):
            print(f"Error: Input file '{input_file}' does not exist.")
            sys.exit(1)

        AutoREACTER(input_file)
    else:
        help_message()

