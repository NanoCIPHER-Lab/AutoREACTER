import shutil
import sys
import os
import json
from pathlib import Path
import time

from PIL import Image

from AutoREACTER.initialization import Initialization
from AutoREACTER.input_parser import InputParser
from AutoREACTER.cache import GetCacheDir, RunDirectoryManager, RetentionCleanup
from AutoREACTER.detectors.functional_groups_detector import FunctionalGroupsDetector
from AutoREACTER.detectors.reaction_detector import ReactionDetector
from AutoREACTER.detectors.non_monomer_detector import NonReactantsDetector
from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import PrepareReactions
from AutoREACTER.reaction_preparation.lunar_client.molecule_3d_preparation import Molecule3DPreparation
from AutoREACTER.reaction_preparation.lunar_client.lunar_api_wrapper import LunarAPIWrapper
from AutoREACTER.reaction_preparation.lunar_client.REACTER_files_builder import REACTERFilesBuilder
from AutoREACTER.sim_setup.simulation_setup import SimulationSetupManager


def _move_image(src: Path, dest: Path) -> None:
    """
    Move all image files from source directory to a dedicated 'images' subdirectory
    within the destination run directory.
    
    Creates the images directory if it doesn't exist. Only PNG and JPEG files are moved.
    
    Args:
        src: Source directory containing image files
        dest: Destination run directory where images subdirectory will be created
    """
    images_dir = dest / "images"
    images_dir.mkdir(parents=True, exist_ok=True)

    for p in src.iterdir():
        if p.is_file() and p.suffix.lower() in [".png", ".jpg", ".jpeg"]:
            try:
                shutil.move(str(p), str(images_dir))
                print(f"[OK] Moved image: {p.name} → {images_dir}")
            except Exception as e:
                print(f"[ERROR] Failed to move image {p.name}: {e}")


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


def run_cleanup(mode: str = "skip") -> None:
    """Run cache cleanup with the specified retention mode.
    
    Args:
        mode: Cleanup mode - number of days (as string), 'all', or 'skip'
    """
    cache_manager = GetCacheDir()
    cleanup = RetentionCleanup(cache_manager.cache_base_dir)
    cleanup.run(mode=mode)


def help_message() -> None:
    """Display usage instructions and available command-line options."""
    print("Usage:")
    print("  python AutoREACTER.py -i <input.json>")
    print("  python AutoREACTER.py -c <days|all|skip>")
    print()
    print("Options:")
    print("  -i, --input        Path to input JSON file")
    print("  -c, --cleanup      Run cache cleanup")
    print()
    print("Cleanup modes:")
    print("  <days>   → delete runs older than N days (e.g., 7, 30)")
    print("  all      → delete all cached runs")
    print("  skip     → no cleanup (default)")


def AutoREACTER(input_file: str) -> None:
    """
    Main orchestrator for the AutoREACTER automated reaction workflow.
    
    This function coordinates the complete pipeline for polymer reaction discovery
    and REACTER input file generation. It performs the following major steps:
    
    1. Input parsing and validation
    2. Functional group detection with visualization
    3. Reaction discovery and user-guided selection
    4. Non-monomer (additive) detection and selection
    5. Reaction template preparation and visualization
    6. 3D molecular geometry optimization
    7. Lunar API integration for advanced processing
    8. Generation of final REACTER input files
    9. Simulation setup
    
    The workflow includes multiple visualization steps and an interactive
    confirmation point before proceeding with computationally intensive tasks.
    
    Args:
        input_file: Path to a JSON file containing monomer definitions,
                    reaction parameters, and simulation settings.
    """
    # Initialize logging, environment, and configuration
    Initialization()

    # Parse and validate user-provided input
    input_parser = InputParser()
    cache_manager = GetCacheDir()
    cache_dir = cache_manager.staging_dir

    # Create a dated directory for this specific run
    run_manager = RunDirectoryManager(cache_manager.cache_base_dir)
    output_dir = run_manager.make_dated_run_dir()

    print(f"[INFO] Staging directory: {cache_dir}")

    # Set up directory for storing all visualization images
    images_dir = cache_dir / "images"
    os.makedirs(images_dir, exist_ok=True)

    # Load input data from JSON
    with open(input_file, "r") as f:
        input_data = json.load(f)

    validated_inputs = input_parser.validate_inputs(input_data)

    # Generate and save initial visualization of input monomers
    try:
        img = input_parser.initial_molecules_image_grid(validated_inputs)
        save_image(img, os.path.join(images_dir, "monomers.png"), "Monomers Grid")
    except Exception:
        print("[WARN] Failed to generate initial molecules image grid")

    # === Functional Group Detection ===
    fg_detector = FunctionalGroupsDetector()
    functional_groups, functional_groups_imgs = fg_detector.functional_groups_detector(
        validated_inputs.monomers
    )

    try:
        img = fg_detector.functional_group_highlighted_molecules_image_grid(functional_groups_imgs)
        save_image(img, os.path.join(images_dir, "functional_groups.png"), "Functional Groups")
    except Exception:
        pass  # Visualization is optional

    # === Reaction Discovery and Selection ===
    reaction_detector = ReactionDetector()
    reaction_instances = reaction_detector.reaction_detector(functional_groups)

    try:
        img = reaction_detector.available_reaction_image_grid(reaction_instances)
        save_image(img, os.path.join(images_dir, "reactions.png"), "Available Reactions")
    except Exception:
        pass

    selected_reactions = reaction_detector.reaction_selection(reaction_instances)

    # === Non-monomer (Additive) Detection ===
    non_monomer_detector = NonReactantsDetector()
    non_reactants = non_monomer_detector.non_monomer_detector(
        validated_inputs, selected_reactions
    )

    try:
        img = non_monomer_detector.non_reactants_to_visualization(non_reactants)
        save_image(img, os.path.join(images_dir, "non_reactants.png"), "Non-Reactants")
    except Exception:
        pass

    updated_inputs = non_monomer_detector.non_reactant_selection(
        validated_inputs, non_reactants
    )

    # === Reaction Template Preparation ===
    prepare_reactions = PrepareReactions(cache_dir)
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
            save_image(img, os.path.join(images_dir, filename))
    except Exception:
        print("[WARN] Failed to generate one or more reaction template visualizations")

    print(
        f"\n[INFO] Reaction preparation completed.\n"
        f"All generated images have been saved in the {images_dir} directory.\n"
        f"Please review the reaction templates before proceeding.\n"
    )

    # Interactive user confirmation before heavy computation
    ok_pass = input("Type 'ok' to continue: ").strip().lower()
    if ok_pass != "ok":
        print("[EXIT] Workflow stopped by user.")
        sys.exit(0)

    # === 3D Geometry Preparation ===
    molecule_3d = Molecule3DPreparation(cache_dir)
    updated_inputs_3d, prepared_reactions_3d = molecule_3d.prepare_molecule_3d_geometry(
        updated_inputs, prepared_reactions
    )

    # === Lunar API Processing ===
    lunar_api = LunarAPIWrapper(cache_dir)
    lunar_results = lunar_api.lunar_workflow(
        updated_inputs_3d, prepared_reactions_3d
    )

    print("\n")
    loading_message("Lunar API workflow completed. Proceeding to build REACTER files")
    time.sleep(1.5)
    print("\n")

    # === Build REACTER Input Files ===
    builder = REACTERFilesBuilder(
        cache_dir=cache_dir,
        updated_inputs_with_3d_mols=updated_inputs_3d
    )

    reacter_files = builder.molecule_template_preparation(
        lunar_files=lunar_results,
        prepared_reactions_with_3d_mols=prepared_reactions_3d
    )

    # Move generated files to final output directory
    reacter_files = run_manager.move_reacter_files(
        reacter_files,
        staging_dir=cache_dir,
        final_dir=output_dir
    )

    _move_image(images_dir, output_dir)

    # Final simulation setup and output
    Simulation_setup_manager = SimulationSetupManager()
    updated_inputs_3d = Simulation_setup_manager.setup_and_write_simulation(
        setup=updated_inputs_3d,
        reacter_files=reacter_files,
        run_dir=output_dir
    )

    print("\n[INFO] AutoREACTER workflow completed successfully.\n")
    print(f"Final REACTER files are located in: {output_dir}")


if __name__ == "__main__":
    args = sys.argv[1:]

    if not args:
        help_message()
        sys.exit(1)

    # Handle input file mode
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

    # Handle cache cleanup mode
    elif "-c" in args or "--cleanup" in args:
        idx = args.index("-c") if "-c" in args else args.index("--cleanup")

        if idx + 1 >= len(args):
            print("Error: Missing cleanup mode.")
            sys.exit(1)

        mode = args[idx + 1]
        run_cleanup(mode)

    else:
        help_message()
