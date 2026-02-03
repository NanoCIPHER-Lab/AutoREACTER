"""
Main API wrapper for interacting with LUNAR (LAMMPS Unified Network for Atomistic Reactions).
Provides functions to run LUNAR atom typing, all2lmp conversion, and bond_react_merge utilities.
This module can handle arbitrary molecule files and manage temporary cache directories.

LUNAR is a tool for preparing molecular systems for LAMMPS molecular dynamics simulations,
particularly for reactive force fields like PCFF-IFF.
"""
import subprocess
import os
import sys
from pathlib import Path
import re
import time

"""
TODO:
# set the cache folder and import from cache.py later
# inside cache folder there should be other folders for each step of LUNAR and the simulation inputs also
# Clean the current cache folder before running any LUNAR steps
# At the end the files needs to be moved to the desired location
"""
# Define cache directories for LUNAR steps
# Base cache directory can be configured via the LUNAR_CACHE_DIR environment variable.
# If not set, it falls back to the original hardcoded path to preserve existing behavior.

# Import LUNAR location detection module
try:
    from locate_lunar import get_LUNAR_loc
except ImportError:
    from .locate_lunar import get_LUNAR_loc

# Get current working directory
cwd = os.getcwd()

# Locate LUNAR installation directory
LUNAR_LOCATION = get_LUNAR_loc(use_gui=False)

# Define paths to LUNAR utility scripts
atom_typing_py = os.path.join(LUNAR_LOCATION, "atom_typing.py")
all2lmp_py = os.path.join(LUNAR_LOCATION, "all2lmp.py")
bond_react_merge_py = os.path.join(LUNAR_LOCATION, "bond_react_merge.py")

import platform
import shutil

def _is_wsl() -> bool:
    """
    Detects whether the current Python process is running under Windows Subsystem for Linux (WSL).
    
    @returns
        True if running under WSL, False otherwise.
    """
    return ("microsoft" in platform.release().lower()) or ("WSL_INTEROP" in os.environ)

def normalize_path(p: str | Path) -> str:
    """
    Normalize a filesystem path for the current execution environment (native or WSL).
    
    Converts path separators and, when running under Windows Subsystem for Linux (WSL), translates Windows-style drive paths (e.g., `C:/path`) to WSL-style paths (`/mnt/c/path`). When not running under WSL, converts WSL-style paths (`/mnt/c/path` or `mnt/c/path`) to Windows-style drive paths (`C:\path`) and normalizes forward/backward slashes. Leading and trailing quotes are removed before conversion.
    
    Parameters:
        p (str | Path): The input path to normalize.
    
    Returns:
        str: The normalized path string appropriate for the current environment.
    """
    p = str(p).strip().strip('"').strip("'")
    # Convert all backslashes to forward slashes for internal consistency
    p = p.replace("\\", "/")

    is_wsl = _is_wsl()

    # 1. Convert Windows (C:/...) to WSL (/mnt/c/...)
    if is_wsl:
        m = re.match(r"^([A-Za-z]):/(.*)$", p)
        if m:
            drive = m.group(1).lower()
            rest = m.group(2)
            return f"/mnt/{drive}/{rest}"
        return p

    # 2. Convert WSL (/mnt/c/...) to Windows (C:/...)
    else:
        # Check for /mnt/c/ or mnt/c/
        m = re.match(r"^/?mnt/([A-Za-z])/(.*)$", p)
        if m:
            drive = m.group(1).upper()
            rest = m.group(2).replace("/", "\\")
            return f"{drive}:\\{rest}"
        
        # If it's already a Windows path, just fix backslashes
        return p.replace("/", "\\")

def move_merge_outputs(src_dir: Path, dst_dir: Path):
    """
    Move LUNAR bond_react_merge output files from a source directory into a destination directory.
    
    This creates the destination directory if it does not exist and relocates files matching these patterns:
    `*_merged.data`, `*_merged.lmpmol`, `force_field.data`, `log.lammps`, and `*.log`.
    
    Parameters:
        src_dir (Path | str): Directory containing generated merge outputs to move.
        dst_dir (Path | str): Destination directory where matching files will be placed; created if missing.
    """
    src_dir = Path(src_dir)
    dst_dir = Path(dst_dir)
    dst_dir.mkdir(parents=True, exist_ok=True)

    for pattern in ("*_merged.data", "*_merged.lmpmol", "force_field.data", "log.lammps", "*.log"):
        for f in src_dir.glob(pattern):
            shutil.move(str(f), str(dst_dir / f.name))

def loading_screen(name: str = "LUNAR") -> None:
    """
    Display an ASCII art banner with a loading spinner animation.
    
    This function provides visual feedback to users during long-running operations.
    
    Args:
        name (str): The name to display in the loading message. Defaults to "LUNAR".
    
    Returns:
        None
    """
    # ASCII art banner for LUNAR
    banner = r"""
    █████       █████  █████ ██████   █████   █████████   ███████████  
    ░░███       ░░███  ░░███ ░░██████ ░░███   ███░░░░░███ ░░███░░░░░███ 
    ░███        ░███   ░███  ░███░███ ░███  ░███    ░███  ░███    ░███ 
    ░███        ░███   ░███  ░███░░███░███  ░███████████  ░██████████  
    ░███        ░███   ░███  ░███ ░░██████  ░███░░░░░███  ░███░░░░░███ 
    ░███      █ ░███   ░███  ░███  ░░█████  ░███    ░███  ░███    ░███ 
    ███████████ ░░████████   █████  ░░█████ █████   █████ █████   █████
    ░░░░░░░░░░░   ░░░░░░░░   ░░░░░    ░░░░░ ░░░░░   ░░░░░ ░░░░░   ░░░░░ 
    """
    
    # Display the banner
    print(banner)
    print(f"Loading {name} ", end="", flush=True)

    # Loading animation characters
    animation = ["-", "\\", "|", "/"]
    
    # Run animation for 10 cycles
    for i in range(10):
        time.sleep(0.1)
        # Write backspace and next animation character
        sys.stdout.write("\b" + animation[i % len(animation)])
        sys.stdout.flush()
    
    print("\nReady!")


def get_ending_integer(s: str) -> int | None:
    """
    Extract and convert the trailing integer from a string.
    
    This function is useful for parsing reaction numbers from filenames
    (e.g., extracting "1" from "pre1" or "post1").
    
    Args:
        s (str): Input string that may end with an integer.
    
    Returns:
        int | None: The integer found at the end of the string, 
                   or None if no integer is found.
    
    Example:
        >>> get_ending_integer("pre1")
        1
        >>> get_ending_integer("data2")
        2
        >>> get_ending_integer("molecule")
        None
    """
    # Regular expression pattern to match one or more digits at the end of the string
    # r'\d+$' matches digits (\d+) at the end ($) of the string
    match = re.search(r'\d+$', s)
    
    if match:
        # Extract the matched substring and convert to integer
        return int(match.group())
    else:
        # No integer found at the end of the string
        return None
    

def run_LUNAR_atom_typing(molecule_files: dict[str, Path], cache_atom_typing: Path) -> dict[str, list[Path]]:
    """
    Run LUNAR atom typing for each input molecule and return paths to the generated typed files.
    
    Parameters:
        molecule_files (dict[str, Path]): Mapping of molecule names to input molecule files (expected topology files such as .mol).
        cache_atom_typing (Path): Directory where atom-typing outputs will be written.
    
    Returns:
        dict[str, list[Path]]: Mapping of molecule names to a two-element list of Paths:
            - [0]: Typed LAMMPS data file named "<stem>_typed.data"
            - [1]: Neighbor table analysis file named "<stem>_typed.nta"
    """
    molecule_files_typed: dict[str, list[Path]] = {}
    
    # Process each molecule file
    for name, file in molecule_files.items():
        # Construct and execute the atom typing command
        # Command format: python atom_typing.py -topo <file> -dir <cache> -ff PCFF-IFF -del-method mass -del-crit 0
        subprocess.run([
            sys.executable,          # Use the current Python interpreter
            atom_typing_py,          # Path to atom_typing.py script
            '-topo', str(file),      # Input topology file
            '-dir', cache_atom_typing,  # Output directory
            "-ff", "PCFF-IFF",       # Force field specification
            "-del-method", "mass",   # Deletion method (based on mass)
            "-del-crit", "0"         # Deletion criterion
        ])
        
        # Store paths to generated files
        # LUNAR generates two files: <stem>_typed.data and <stem>_typed.nta
        molecule_files_typed[name] = [
            Path(cache_atom_typing) / f"{Path(file).stem}_typed.data",
            Path(cache_atom_typing) / f"{Path(file).stem}_typed.nta"
        ]

    time.sleep(1)  # Ensure all files are written before returning
    for name, paths in molecule_files_typed.items():
        for path in paths:
            if path.is_file():
                print(f"[LUNAR bond_react_merge] Generated file for {name} at: {path}")
            else:
                raise FileNotFoundError(
                    f"Expected output file for {name} not found at {path}"
            )
    return molecule_files_typed


def run_LUNAR_all2lmp(molecule_files_typed: dict[str, list[Path]], cache_all2lmp: Path) -> dict[str, str]:
    """
    Convert typed molecule files and their corresponding NTA files into LAMMPS data filenames using LUNAR's all2lmp conversion.
    
    Parameters:
        molecule_files_typed (dict[str, list[Path]]): Mapping from molecule name to a two-element list of Paths [data_file, nta_file].
        cache_all2lmp (Path): Directory where converted LAMMPS files are written.
    
    Returns:
        dict[str, str]: Mapping from molecule name to the generated LAMMPS data filename (basename, not full path).
    """
    # Sometimes LUNAR can't identify the .frc file location automatically,
    # so we provide the full path to the force field file
    frc_file = os.path.join(LUNAR_LOCATION, "frc_files", "pcff.frc")
    
    molecule_files_all2lmp: dict[str, str] = {}
    
    # Process each typed molecule
    for name, (data_file, nta_file) in molecule_files_typed.items():
        # Construct and execute the all2lmp conversion command
        # Command format: python all2lmp.py -topo <data> -nta <nta> -frc <frc> -dir <cache>
        subprocess.run([
            sys.executable,          # Use the current Python interpreter
            all2lmp_py,              # Path to all2lmp.py script
            '-topo', str(data_file), # Typed data file
            '-nta', str(nta_file),   # NTA file
            '-frc', frc_file,        # Force field file
            '-dir', cache_all2lmp  # Output directory
        ])
        
        # Store the generated LAMMPS data filename
        # LUNAR generates: <stem>_IFF.data (e.g., data1_typed_IFF.data)
        molecule_files_all2lmp[name] = f"{data_file.stem}_IFF.data"

    time.sleep(1)  # Ensure all files are written before returning
    for name, path in molecule_files_all2lmp.items():
        # Construct full path to the generated file
        path = Path(cache_all2lmp) / path
        if Path(path).is_file():
            print(f"[LUNAR bond_react_merge] Generated file for {name} at: {path}")
        else:
            raise FileNotFoundError(
                f"Expected output file for {name} not found at {path}"
            )
    return molecule_files_all2lmp


def write_bond_react_merge_input(
    molecule_files_all2lmp: dict[str, str],
    cache_bond_react_merge: Path,
    cache_all2lmp: Path
) -> Path:
    """
    Create a merge_input.txt for Bond React Merge and prepare the merge cache directory.
    
    Generates a merge_input.txt in cache_bond_react_merge that lists each tagged LAMMPS file (from cache_all2lmp) with a short comment, normalizes paths written into the file, ensures the merge cache directory exists, and moves relevant all2lmp outputs into the merge cache.
    
    Parameters:
        molecule_files_all2lmp (dict[str, str]): Mapping from file tag (e.g., "data_1", "pre_1", "post_1") to the LAMMPS filename produced by the all2lmp step.
        cache_bond_react_merge (Path): Directory where merge_input.txt and merge results will be written (created if missing).
        cache_all2lmp (Path): Directory containing the LAMMPS files referenced by molecule_files_all2lmp.
    
    Returns:
        Path: The path to the cache_bond_react_merge directory containing the written merge_input.txt.
    """
    # Ensure the provided paths are normalized (resolving redundant separators, symlinks, etc.).
    cache_bond_react_merge_normalized = Path(normalize_path(cache_bond_react_merge))
    cache_all2lmp = Path(normalize_path(cache_all2lmp))

    # Create the output directory if it does not already exist.
    os.makedirs(cache_bond_react_merge_normalized, exist_ok=True)
    # Initialize the merge input file content with header comments and column titles.
    merge_files = f"""# anything following the "#" character will be ignored

{"# file-tag":<10}{"filename":<150}{"comment (required)"}
"""

    # Append one line per molecule file describing the tag, absolute path, and reason comment.
    for name, lammps_file in molecule_files_all2lmp.items():
        if name.startswith("data"):
            comment = "# This datafile will have all coeffs in it"
        else:
            rxn_number = get_ending_integer(name)
            comment = f"# for rxn{rxn_number}"

        full_file = normalize_path(Path(cache_all2lmp) / lammps_file)
        merge_files += f"{name:<10}{full_file:<150}{comment}\n"

    # Specify where the merge outputs should be written to within the merge_input file.
    merge_files += f"""
# Specify the parent_directory of where to write results (optional)
"""
    # Write the complete merge input text to the target file.
    with open(cache_bond_react_merge_normalized / "merge_input.txt", "w") as merge_input_file:
        merge_input_file.write(merge_files)

    move_merge_outputs(cache_all2lmp, cache_bond_react_merge)

    return cache_bond_react_merge



def run_bond_react_merge(merge_input_file_path: Path, 
                         molecule_files_all2lmp: dict[str, str]) -> dict[str, Path]:
    """
                         Run LUNAR's bond_react_merge using the merge_input.txt in the given directory and collect the generated output file paths.
                         
                         Parameters:
                             merge_input_file_path (Path): Directory containing merge_input.txt (working directory for the merge run).
                             molecule_files_all2lmp (dict[str, str]): Mapping of input tags (e.g., "data1", "pre1", "post1") to their LAMMPS data filenames.
                         
                         Returns:
                             dict[str, Path]: Mapping from each input tag and the special key "force_field_data" to the corresponding generated output Path.
                         
                         Raises:
                             FileNotFoundError: If any expected output file is missing after the merge run.
                         """
    molecule_files_bond_react_merge: dict[str, Path] = {}
    merge_input_file_path = Path(normalize_path(merge_input_file_path))
    # Construct and execute the bond_react_merge command
    # Command format: python bond_react_merge.py -files infile:<input_file> -atomstyle full
    subprocess.run(
        [
            sys.executable,                     # Use the current Python interpreter
            bond_react_merge_py,                # Path to bond_react_merge.py script
            "-files", f"infile:merge_input.txt", # Input file specification
            "-atomstyle", "full",               # LAMMPS atom style
        ],
        cwd=merge_input_file_path,              # Run in the directory with input file
        check=True                              # Raise exception if command fails
    )
    
    # Map input names to expected output filenames
    for name in molecule_files_all2lmp.keys():
        if name.startswith("data"):
            # Data files generate merged .data files
            molecule_files_bond_react_merge[name] = (
                os.path.join(merge_input_file_path, f"{name}_typed_IFF_merged.data")
            )
        elif name.startswith("pre"):
            # Pre-reaction files generate .lmpmol files
            rxn_number = get_ending_integer(name)
            molecule_files_bond_react_merge[name] = (
                os.path.join(merge_input_file_path, f"pre_{rxn_number}_typed_IFF_merged.lmpmol")
            )
        elif name.startswith("post"):
            # Post-reaction files generate .lmpmol files
            rxn_number = get_ending_integer(name)
            molecule_files_bond_react_merge[name] = (
                os.path.join(merge_input_file_path, f"post_{rxn_number}_typed_IFF_merged.lmpmol")
            )
    
    # Add the force field data file (generated by bond_react_merge)
    molecule_files_bond_react_merge["force_field_data"] = (
        os.path.join(merge_input_file_path, "force_field.data")
    )
    
    # Verify that all expected output files were created
    for name, path in molecule_files_bond_react_merge.items():
        if os.path.isfile(path):
            print(f"[LUNAR bond_react_merge] Generated file for {name} at: {path}")
        else:
            raise FileNotFoundError(
                f"Expected output file for {name} not found at {path}"
            )
    
    return molecule_files_bond_react_merge


def lunar_workflow(molecule_files: dict[str, Path], cache_dir: Path) -> dict[str, Path]:
    """
    Execute the full LUNAR preparation workflow and return generated output paths.
    
    Runs atom typing, converts typed files to LAMMPS format, writes the bond_react_merge input, and runs bond_react_merge. Creates and uses a LUNAR cache directory under `cache_dir/lunar` unless overridden by the `LUNAR_CACHE_DIR` environment variable; subdirectories `atom_typing`, `all2lmp`, and `bond_react_merge` are created and used for intermediate and final files.
    
    Parameters:
        molecule_files (dict[str, Path]): Mapping of molecule identifiers to their input molecule file paths.
        cache_dir (Path): Base directory under which the LUNAR cache directory will be created (unless `LUNAR_CACHE_DIR` is set).
    
    Returns:
        dict[str, Path]: Mapping of molecule identifiers and special keys (e.g., force field file) to their final output file paths produced by the bond_react_merge step.
    """
    lunar_cache = os.path.join(cache_dir, "lunar")
    os.makedirs(lunar_cache, exist_ok=True)
    LUNAR_CACHE_DIR = os.environ.get(
        "LUNAR_CACHE_DIR",
        lunar_cache
    )
    cache_atom_typing = os.path.join(LUNAR_CACHE_DIR, "atom_typing")
    cache_all2lmp = os.path.join(LUNAR_CACHE_DIR, "all2lmp")
    cache_bond_react_merge = os.path.join(LUNAR_CACHE_DIR, "bond_react_merge")

    for p in (LUNAR_CACHE_DIR, cache_atom_typing, cache_all2lmp, cache_bond_react_merge):
        os.makedirs(p, exist_ok=True)
    # Ensure cache directory exists (redundant but safe)
    # Step 0: Show loading screen
    loading_screen("LUNAR Workflow")
    
    # Step 1: Run atom typing
    typed_files = run_LUNAR_atom_typing(molecule_files, cache_atom_typing)
    
    # Step 2: Run all2lmp conversion
    all2lmp_files = run_LUNAR_all2lmp(typed_files, cache_all2lmp)
    
    # Step 3: Write bond_react_merge input file
    merge_input_path = write_bond_react_merge_input(all2lmp_files, cache_bond_react_merge, cache_all2lmp)
    
    # Step 4: Run bond_react_merge
    lunar_files = run_bond_react_merge(merge_input_path, all2lmp_files)
    
    return lunar_files
    
if __name__ == "__main__":
    """
    Example usage of the LUNAR API wrapper.
    
    This block demonstrates how to use the lunar_workflow function with
    sample molecule files. In practice, users should modify the molecule_files
    dictionary to point to their actual .mol files.
    """
    # Import Path here to avoid circular imports if this module is imported elsewhere
    from pathlib import Path
    try:
        from molecule_3d_preparation import prepare_3d_molecule
    except ImportError:
        from lunar_api_wrapper.molecule_3d_preparation import prepare_3d_molecule
    from rdkit import Chem
    from rdkit.Chem import AllChem
    # Example usage: construct reactants, run a SMARTS-based reaction, prepare 3D structures,
    # and write optimized .mol files to the cache directory.

    # Use the local `cache` variable (TODO: replace with config/cache module)
    cache = r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\lunar"

    # Example SMILES for two reactants 
    reactant_smiles1 = "C1=CC=C(C(=C1)C(=O)O)O"
    reactant_smiles2 = "OCCC(O)=O"

    # SMARTS describing esterification example
    reaction_smarts = (
        "[O;!$(OC=*):1]-[H:3].[CX3:2](=[O:5])[OX2H1:4]>>[OX2:1]-[CX3:2](=[O:5]).[O:4]-[H:3]"
    )

    # Load reactants from SMILES and add explicit hydrogens (recommended for 3D embedding)
    reactant1 = Chem.MolFromSmiles(reactant_smiles1)
    reactant1 = Chem.AddHs(reactant1)
    reactant2 = Chem.MolFromSmiles(reactant_smiles2)
    reactant2 = Chem.AddHs(reactant2)

    # Create an RDKit reaction from SMARTS
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)

    # Combine reactants into a single RDKit Mol for visualization/complex handling
    combined_reactants = Chem.CombineMols(reactant1, reactant2)

    # Run the reaction to produce product sets. RunReactants returns a tuple of tuples,
    # where each outer tuple is one possible product set. We'll take the first generated set.
    products = rxn.RunReactants((reactant1, reactant2))
    if not products:
        raise RuntimeError("Reaction produced no products for the provided reactants and SMARTS.")
    products = products[0]  # Take the first product set

    # Combine produced product fragments into a single RDKit Mol
    combined_products = Chem.CombineMols(*products)

    # Create a mapping of names -> RDKit Mol objects to be prepared (embedded/optimized)
    molecule_dict = {
        "data1": reactant1,
        "data2": reactant2,
        "pre1": combined_reactants,   # combined reactant complex
        "post1": combined_products,   # combined product complex
    }

    print("Cache directory:", cache)
    prepared_molecules = prepare_3d_molecule(cache_dir=cache, molecule_dict=molecule_dict)
    print("Prepared molecules and saved paths:", prepared_molecules)
    
    # Define sample molecule files for testing
    # Note: These paths are specific to the original developer's system
    molecule_files: dict[str, Path] = {
        "data_1": Path(r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\lunar\data1.mol"),
        "data_2": Path(r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\lunar\data2.mol"),
        "pre_1":  Path(r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\lunar\pre1.mol"),
        "post_1": Path(r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\lunar\post1.mol"),
    }
    
    # Execute the complete LUNAR workflow
    lunar_files = lunar_workflow(molecule_files)
    
    # Output summary of generated files
    print("\n" + "=" * 80)
    print("LUNAR Workflow Complete!")
    print("Generated files:")
    for name, path in lunar_files.items():
        print(f"  {name}: {path}")
    print("=" * 80)
    import json
    with open(Path(test_cache) / "lunar_workflow_summary.txt", "w") as summary_file:
        summary_file.write(json.dumps(lunar_files, indent=4, default=str))