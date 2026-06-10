import json
import shutil
from pathlib import Path
from dataclasses import dataclass

# Import internal modules
from AutoREACTER.initialization import Initialization
from AutoREACTER.input_parser import InputParser
from AutoREACTER.cache import GetCacheDir

@dataclass
class ARXSession:
    """
    Holds the validated inputs and directory paths for a single AutoREACTER run.
    This acts as the 'state object' passed through the pipeline.
    """
    inputs: object
    staging_dir: Path
    output_dir: Path
    images_dir: Path

def _resolve_input_path(input_file_path: str) -> Path:
    """
    Resolves the input file path, ensuring it exists and is a valid JSON file.
    """
    input_path = Path(input_file_path).resolve()
    
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    
    if input_path.suffix.lower() != ".json" or not input_path.is_file():
        raise ValueError(f"Input file must be a JSON file: {input_path}")
    
    return input_path

def _clear_directory(path: Path):
    """
    Remove all contents of a directory (including subfolders) without deleting the root directory itself.
    """
    if not path.is_dir():
        return
    
    for item in path.iterdir():
        if item.is_file() or item.is_symlink():
            item.unlink()
        elif item.is_dir():
            shutil.rmtree(item) 

def read_input(input_file_path: str, clear_staging: bool = True) -> ARXSession:
    """
    Standard read function to initialize the AutoREACTER environment.
    
    1. Resolves the input file location.
    2. Sets up a temporary staging directory.
    3. Validates the JSON inputs.
    4. Creates a permanent output directory named after the simulation.
    """
    # 1. Resolve paths
    input_path = _resolve_input_path(input_file_path)

    # 2. Global Initialization
    Initialization()
    
    # 3. Setup temporary staging
    cache_manager = GetCacheDir(clear_staging=clear_staging)
    staging_dir = cache_manager.staging_dir

    # 4. Parse and Validate Inputs First (so we can get the simulation_name)
    input_parser = InputParser()
    with open(input_path, "r") as f:
        input_data = json.load(f)
    
    validated_inputs = input_parser.validate_inputs(input_data)

    # 5. Setup Output Directories based on Simulation Name
    base_output_dir = input_path.parent / "AutoREACTER_outputs"
    base_output_dir.mkdir(parents=True, exist_ok=True)

    # Create the specific folder for THIS simulation (sanitize to avoid path traversal)
    sim_name = str(validated_inputs.simulation_name)
    if Path(sim_name).name != sim_name or sim_name in {".", ".."}:
        raise ValueError(f"Invalid simulation_name for output directory: {sim_name!r}")
    output_dir = base_output_dir / sim_name

    if output_dir.exists():
        _clear_directory(output_dir)  # Only clear this specific run's old files
    output_dir.mkdir(parents=True, exist_ok=True)

    # Now create the images folder safely inside the simulation folder
    images_dir = Path((output_dir) / "images")
    images_dir.mkdir(parents=True, exist_ok=True)

    # 6. Return the State Object
    print(f"[INFO] Initialized AutoREACTER Session")
    print(f"[INFO] Simulation Name: {validated_inputs.simulation_name}")
    print(f"[INFO] Input File: {input_path}")
    print(f"[INFO] Temporary Staging: {staging_dir}")
    print(f"[INFO] Final Outputs will save to: {output_dir}")

    return ARXSession(
        inputs=validated_inputs,
        staging_dir=staging_dir,
        output_dir=output_dir,
        images_dir=images_dir
    )