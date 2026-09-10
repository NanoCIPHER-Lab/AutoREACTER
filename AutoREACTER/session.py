from __future__ import annotations  # 1. Must be the first line
from typing import TYPE_CHECKING
import json
import shutil
from dataclasses import dataclass, field
from pathlib import Path

# Import internal modules
from AutoREACTER.initialization import Initialization
from AutoREACTER.input_parser import InputParser, SimulationSetup
from AutoREACTER.detectors.reaction_detector import ReactionInstance
from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import ReactionMetadata
from AutoREACTER.reaction_preparation.ff_wrapper.ff_wrapper import FFFiles
from AutoREACTER.reaction_preparation.ff_wrapper.REACTER_files_builder import REACTERFiles
if TYPE_CHECKING:
    from AutoREACTER.detectors.functional_groups_detector import MonomerRole

@dataclass(slots=True)
class Session:
    """
    Runtime state container for one AutoREACTER workflow.

    A Session stores the validated input model, run directories, and all
    intermediate/final objects generated as the pipeline progresses.

    The object is intentionally passed through the full AutoREACTER pipeline so
    each stage can attach its own outputs without requiring large return-value
    chains.

    Pipeline state
    --------------
    inputs:
        Validated simulation setup parsed from the input JSON.

    staging_dir:
        Temporary working directory used for intermediate/cache files.

    output_dir:
        Final AutoREACTER output directory for this run.

    images_dir:
        Directory where molecule, functional-group, reaction, and template
        visualization images are saved.

    monomer_roles:
        MonomerRole objects after functional-group classification.

    reaction_instances:
        Detected reaction candidates before full reaction preparation.

    non_reactants:
        MonomerRole objects selected as non-reactive species.

    reaction_metadata:
        Prepared reaction objects. Each ReactionMetadata object stores RDKit
        mappings, template atoms, initiators, edge atoms, byproducts, generated
        REACTER template files, map files, and activity status.

    ff_files:
        Raw force-field and LUNAR output files before REACTER-specific file
        organization.

    reacter_files:
        Final REACTER file bundle. This should contain run-level files plus the
        final monomer and reaction metadata lists used by LAMMPS writers.
    """

    # Core run configuration
    inputs: SimulationSetup
    staging_dir: Path
    output_dir: Path
    images_dir: Path 

    # Pipeline-generated state
    monomer_roles: list["MonomerRole"] = field(default_factory=list)
    reaction_instances: list[ReactionInstance] = field(default_factory=list)
    non_reactants: list["MonomerRole"] = field(default_factory=list)
    reaction_metadata: list[ReactionMetadata] = field(default_factory=list)

    # File bundles generated later in the pipeline
    ff_files: FFFiles | None = None
    reacter_files: REACTERFiles | None = None

    # Runtime counters / sub-sessions attached during reaction preparation
    reaction_id_counter: int = 0
    reaction_progression_session: object | None = None


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

def _resolve_output_dir(
    raw_output_dir: str | None,
    input_path: Path,
    simulation_name: str,
) -> Path:
    """
    Resolve the final AutoREACTER output directory.

    Behavior:
    - If output_dir is missing, None, or empty, use the default location beside
      the input JSON:
          input_json_folder / AutoREACTER_outputs / simulation_name

    - If output_dir is relative, resolve it relative to the input JSON folder.

    - If output_dir is a Linux/WSL absolute path, use it directly.

    - If output_dir is a Windows-style path such as C:/Users/... or C:\\Users\\...,
      convert it to the WSL form /mnt/c/Users/....

    This function returns an absolute path. Directory creation/clearing happens
    in read_input().
    """
    if raw_output_dir is None or str(raw_output_dir).strip() == "":
        return (
            input_path.parent
            / "AutoREACTER_outputs"
            / simulation_name
        ).resolve()

    raw_output_dir = str(raw_output_dir).strip()

    # Windows path while running from WSL/Linux.
    if (
        len(raw_output_dir) >= 3
        and raw_output_dir[1] == ":"
        and raw_output_dir[2] in {"/", "\\"}
    ):
        drive = raw_output_dir[0].lower()
        rest = raw_output_dir[3:].replace("\\", "/")
        return Path(f"/mnt/{drive}/{rest}").resolve()

    output_dir = Path(raw_output_dir).expanduser()

    if not output_dir.is_absolute():
        output_dir = input_path.parent / output_dir

    return output_dir.resolve()

def read_input(input_file_path: str, clear_staging: bool = True) -> Session:
    """
    Standard read function to initialize the AutoREACTER environment.
    
    1. Resolves the input file location.
    2. Sets up a temporary staging directory.
    3. Validates the JSON inputs.
    4. Creates a permanent output directory named after the simulation.
    """
    from AutoREACTER.cache import GetCacheDir # Import here to avoid circular imports with Session
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

    # 5. Setup Output Directory
    sim_name = str(validated_inputs.simulation_name)

    if Path(sim_name).name != sim_name or sim_name in {".", ".."}:
        raise ValueError(f"Invalid simulation_name for output directory: {sim_name!r}")

    raw_output_dir = input_data.get("output_dir", None)

    output_dir = _resolve_output_dir(
        raw_output_dir=raw_output_dir,
        input_path=input_path,
        simulation_name=sim_name,
    )

    if output_dir.exists() and not output_dir.is_dir():
        raise ValueError(
            f"Resolved output_dir exists but is not a directory: {output_dir}"
        )

    if output_dir.exists():
        _clear_directory(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    images_dir = output_dir / "images"
    images_dir.mkdir(parents=True, exist_ok=True)

    # 6. Return the State Object
    print(f"[INFO] Initialized AutoREACTER Session")
    print(f"[INFO] Simulation Name: {validated_inputs.simulation_name}")
    print(f"[INFO] Input File: {input_path}")
    print(f"[INFO] Temporary Staging: {staging_dir}")
    print(f"[INFO] Final Outputs will save to: {output_dir}")

    return Session(
        inputs=validated_inputs,
        staging_dir=staging_dir,
        output_dir=output_dir,
        images_dir=images_dir,
    )

