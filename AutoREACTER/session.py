from __future__ import annotations  # 1. Must be the first line
from typing import TYPE_CHECKING
import json
import shutil
from pathlib import Path
from dataclasses import dataclass

# Import internal modules
from AutoREACTER.initialization import Initialization
from AutoREACTER.input_parser import InputParser, SimulationSetup
from AutoREACTER.detectors.reaction_detector import ReactionInstance
from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import ReactionMetadata
from AutoREACTER.reaction_preparation.ff_wrapper.ff_wrapper import FFFiles
from AutoREACTER.reaction_preparation.ff_wrapper.REACTER_files_builder import REACTERFiles
if TYPE_CHECKING:
    from AutoREACTER.detectors.functional_groups_detector import MonomerRole

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal


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

    Workflow options
    ----------------
    deep_search_enabled:
        If True, AutoREACTER performs deeper functional-group/reaction searching
        when supported by the detector logic.

    reaction_iteration_depth:
        Maximum reaction-progression depth. An integer limits the number of
        progression iterations. The string "all" can be used for unlimited or
        exhaustive progression if supported downstream.

    use_wildcard_atom_types:
        If True, generated templates may use wildcard atom typing behavior where
        supported.

    deduplicate_reaction_templates:
        If True, duplicate LAMMPS pre/post reaction template pairs are detected
        and disabled by setting ReactionMetadata.activity_stats = False.

    write_second_reaction_stage:
        If True, AutoREACTER writes the second reaction-stage LAMMPS input files.
        If False, only the first-stage reaction input is written.
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

    # Workflow options
    deep_search: bool = True
    reaction_iteration_depth: int | Literal["all"] = 5
    wildcards: bool = True
    deduplicate_reaction_templates: bool = True
    write_second_reaction_stage: bool = True

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

def _normalize_output_dir(raw_output_dir: str, input_path: Path) -> Path:
    """
    Normalize output_dir from JSON.

    Supports:
    - Linux/WSL absolute paths: /mnt/c/...
    - Windows paths: C:/Users/...
    - Relative paths: AutoREACTER_outputs/run_name
    """
    raw_output_dir = str(raw_output_dir).strip()

    # Handle Windows-style path when running in WSL/Linux.
    if len(raw_output_dir) >= 3 and raw_output_dir[1] == ":" and raw_output_dir[2] in {"/", "\\"}:
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

    if raw_output_dir is not None:
        output_dir = _normalize_output_dir(raw_output_dir, input_path)

    else:
        # Backward-compatible default behavior.
        output_dir = input_path.parent / "AutoREACTER_outputs" / sim_name

    if output_dir.exists():
        _clear_directory(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    images_dir = output_dir / "images"
    images_dir.mkdir(parents=True, exist_ok=True)

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

    return Session(
        inputs=validated_inputs,
        staging_dir=staging_dir,
        output_dir=output_dir,
        images_dir=images_dir
    )

