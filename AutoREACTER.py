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
    Move image files from source to a dedicated 'images' subdirectory in the destination.
    
    Creates the images directory if it doesn't exist and moves all PNG/JPG files.
    
    Args:
        src: Source directory containing images
        dest: Destination run directory
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


def loading_message(message, duration=3, interval=0.5):
    print(f"[INFO] {message}", end="", flush=True)
    steps = int(duration / interval)

    for _ in range(steps):
        print(".", end="", flush=True)
        time.sleep(interval)

    print()  # newline


def save_image(img: Image.Image, path: str, label: str = "Image") -> None:
    """
    Save a PIL Image to the specified path with error handling and user feedback.
    
    Args:
        img: PIL Image object to save
        path: Destination file path
        label: Human-readable label for logging messages
    """
    if img is None:
        print(f"[WARN] {label}: Image is None, skipping save.")
        return

    try:
        if hasattr(img, "save"):
            img.save(path)
            print(f"[OK] {label} saved → {path}\n")
        else:
            print(f"[ERROR] {label}: Object has no .save() method")
    except Exception as e:
        print(f"[ERROR] Failed to save {label} to {path}: {e}")

def run_cleanup(mode: str = "skip") -> None:
    """Run cache cleanup based on the specified mode."""
    cache_manager = GetCacheDir()
    cleanup = RetentionCleanup(cache_manager.cache_base_dir)
    cleanup.run(mode=mode)

def help_message() -> None:
    print("Usage:")
    print("  python AutoREACTER.py -i <input.json>")
    print("  python AutoREACTER.py -c <days|all|skip>")
    print()
    print("Options:")
    print("  -i, --input        Input JSON file")
    print("  -c, --cleanup      Cleanup cache")
    print()
    print("Cleanup modes:")
    print("  <days>   → delete runs older than given number of days (e.g., 7, 30)")
    print("  all      → delete all cached runs")
    print("  skip     → no cleanup")

def AutoREACTER(input_file: str) -> None:
    """
    Main workflow orchestrator for the AutoREACTER pipeline.
    
    This function coordinates the entire automated reaction discovery and preparation
    workflow for polymer chemistry, including:
      - Input validation and visualization
      - Functional group detection
      - Reaction discovery and selection
      - Non-monomer (additive) detection
      - 3D geometry preparation
      - Integration with the Lunar API
      - Generation of REACTER input files
    
    Args:
        input_file: Path to a JSON file containing monomer definitions and parameters.
    
    The function creates a dated output directory, generates multiple diagnostic images,
    and guides the user through an interactive review step before proceeding with
    computationally intensive 3D preparation and simulation file generation.
    """
    # Initialize the AutoREACTER environment and logging
    Initialization()

    # Parse and validate user input
    input_parser = InputParser()
    cache_manager = GetCacheDir()
    cache_dir = cache_manager.staging_dir

    # Set up run-specific output directory
    run_manager = RunDirectoryManager(cache_manager.cache_base_dir)
    output_dir = run_manager.make_dated_run_dir()

    print(f"[INFO] Staging directory: {cache_dir}")

    # Create directory for storing visualization images
    images_dir = cache_dir / "images"
    os.makedirs(images_dir, exist_ok=True)

    # Load and validate input data
    with open(input_file, "r") as f:
        input_data = json.load(f)

    validated_inputs = input_parser.validate_inputs(input_data)

    # Generate and save initial monomers visualization
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
        pass

    # === Reaction Discovery ===
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
            (None, "templates_all.png"),
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

    # Interactive user confirmation
    ok_pass = input("Type 'ok' to continue: ").strip().lower()
    if ok_pass != "ok":
        print("[EXIT] Workflow stopped by user.")
        sys.exit(0)

    # === 3D Geometry Preparation ===
    molecule_3d = Molecule3DPreparation(cache_dir)
    updated_inputs_3d, prepared_reactions_3d = molecule_3d.prepare_molecule_3d_geometry(
        updated_inputs, prepared_reactions
    )

    # === Lunar API Workflow ===
    lunar_api = LunarAPIWrapper(cache_dir)
    lunar_results = lunar_api.lunar_workflow(
        updated_inputs_3d, prepared_reactions_3d
    )

    print("\n")
    loading_message("Lunar API workflow completed. Proceeding to build REACTER files")
    time.sleep(1.5)
    print("\n")

    # === Build final REACTER files ===
    builder = REACTERFilesBuilder(
        cache_dir=cache_dir,
        updated_inputs_with_3d_mols=updated_inputs_3d
    )

    reacter_files = builder.molecule_template_preparation(
        lunar_files=lunar_results,
        prepared_reactions_with_3d_mols=prepared_reactions_3d
    )

    # Move final files to output directory
    reacter_files = run_manager.move_reacter_files(
        reacter_files,
        staging_dir=cache_dir,
        final_dir=output_dir
    )

    _move_image(images_dir, output_dir)

    print("\n[INFO] AutoREACTER workflow completed successfully.\n")
    print(f"Final REACTER files are located in: {output_dir}")

    # Debug prints (can be commented out in production)
    # print(f"Generated REACTER files: {reacter_files}")
    # print(f"Updated inputs with 3D molecules: {updated_inputs_3d}")
    # print(f"Prepared reactions with 3D molecules: {prepared_reactions_3d}")

    Simulation_setup_manager = SimulationSetupManager()
    _ = Simulation_setup_manager.populate_physical_parameters(updated_inputs_3d)


if __name__ == "__main__":
    args = sys.argv[1:]

    if not args:
        help_message()
        sys.exit(1)

    # INPUT MODE 
    if "-i" in args or "--input" in args:
        if "-i" in args:
            idx = args.index("-i")
        else:
            idx = args.index("--input")

        if idx + 1 >= len(args):
            print("Error: Missing input file.")
            sys.exit(1)

        input_file = args[idx + 1]

        if not os.path.isfile(input_file):
            print(f"Error: Input file '{input_file}' does not exist.")
            sys.exit(1)

        AutoREACTER(input_file)

    # CLEANUP MODE 
    elif "-c" in args or "--cleanup" in args:
        if "-c" in args:
            idx = args.index("-c")
        else:
            idx = args.index("--cleanup")

        if idx + 1 >= len(args):
            print("Error: Missing cleanup mode.")
            sys.exit(1)

        mode = args[idx + 1]
        run_cleanup(mode)

    else:
        help_message()