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
test_cache = r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\lunar"
test_cache_atom_typing = f"{test_cache}\\atom_typing"
test_cache_all2lmp = f"{test_cache}\\all2lmp"
test_cache_bond_react_merge = f"{test_cache}\\bond_react_merge"

for p in (test_cache, test_cache_atom_typing, test_cache_all2lmp, test_cache_bond_react_merge):
    os.makedirs(p, exist_ok=True)
# Ensure cache directory exists (redundant but safe)

# Import LUNAR location detection module
from locate_LUNAR import get_LUNAR_loc

# Get current working directory
cwd = os.getcwd()

# Locate LUNAR installation directory
LUNAR_LOCATION = get_LUNAR_loc(use_gui=False)

# Define paths to LUNAR utility scripts
atom_typing_py = os.path.join(LUNAR_LOCATION, "atom_typing.py")
all2lmp_py = os.path.join(LUNAR_LOCATION, "all2lmp.py")
bond_react_merge_py = os.path.join(LUNAR_LOCATION, "bond_react_merge.py")


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
    

def run_LUNAR_atom_typing(molecule_files: dict[str, Path]) -> dict[str, list[Path]]:
    """
    Run LUNAR atom typing on a set of molecule files.
    
    This function assigns atom types to molecules using the PCFF-IFF force field
    and generates typed data files and NTA (Neighbor Table Analysis) files.
    
    Args:
        molecule_files (dict[str, Path]): Dictionary mapping molecule names to 
                                         their file paths. Files should be in 
                                         .mol format.
    
    Returns:
        dict[str, list[Path]]: Dictionary mapping molecule names to a list of 
                              two Path objects:
                              - [0]: Typed data file (.data)
                              - [1]: NTA file (.nta)
    
    Note:
        Uses the following LUNAR atom_typing.py command:
        python3 atom_typing.py -topo <file.mol> -dir <cache_dir> -ff PCFF-IFF 
               -del-method mass -del-crit 0
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
            '-dir', test_cache_atom_typing,  # Output directory
            "-ff", "PCFF-IFF",       # Force field specification
            "-del-method", "mass",   # Deletion method (based on mass)
            "-del-crit", "0"         # Deletion criterion
        ])
        
        # Store paths to generated files
        # LUNAR generates two files: <stem>_typed.data and <stem>_typed.nta
        molecule_files_typed[name] = [
            Path(test_cache_atom_typing) / f"{file.stem}_typed.data",
            Path(test_cache_atom_typing) / f"{file.stem}_typed.nta"
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


def run_LUNAR_all2lmp(molecule_files_typed: dict[str, list[Path]]) -> dict[str, str]:
    """
    Convert typed molecule files to LAMMPS data format using LUNAR's all2lmp utility.
    
    This function takes the typed data and NTA files and converts them to 
    LAMMPS-readable data files with force field parameters.
    
    Args:
        molecule_files_typed (dict[str, list[Path]]): Dictionary mapping molecule 
                                                     names to a list of two Paths:
                                                     [data_file, nta_file]
    
    Returns:
        dict[str, str]: Dictionary mapping molecule names to the generated 
                       LAMMPS data filenames (without full path).
    
    Note:
        Uses the following LUNAR all2lmp.py command:
        python3 all2lmp.py -topo <data_file> -nta <nta_file> -frc <frc_file> 
               -dir <cache_dir>
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
            '-dir', test_cache_all2lmp  # Output directory
        ])
        
        # Store the generated LAMMPS data filename
        # LUNAR generates: <stem>_IFF.data (e.g., data1_typed_IFF.data)
        molecule_files_all2lmp[name] = f"{data_file.stem}_IFF.data"

    time.sleep(1)  # Ensure all files are written before returning
    for name, path in molecule_files_all2lmp.items():
        # Construct full path to the generated file
        path = Path(test_cache_all2lmp) / path
        if Path(path).is_file():
            print(f"[LUNAR bond_react_merge] Generated file for {name} at: {path}")
        else:
            raise FileNotFoundError(
                f"Expected output file for {name} not found at {path}"
            )
    return molecule_files_all2lmp


def write_bond_react_merge_input(molecule_files_all2lmp: dict[str, str]) -> Path:
    """
    Generate input file for LUNAR's bond_react_merge utility.
    
    This function creates a configuration file that specifies which LAMMPS data
    files to merge and how to process them for reactive simulations.
    
    Args:
        molecule_files_all2lmp (dict[str, str]): Dictionary mapping molecule 
                                                names to LAMMPS data filenames.
    
    Returns:
        Path: Path to the directory containing the generated merge input file.
    
    Note:
        The generated file follows LUNAR's bond_react_merge input format with
        sections for file tags, filenames, and comments.
    """
    # Start building the merge input file content
    # Header with format explanation
    merge_files: str = f"""# anything following the "#" character will be ignored 
 
# Specify a desired path to append to the front of each filename (optional however if not present 
# and the files are a path from LUNAR on your computer each file below must have that path 
# specified in front of the filename) 

path = "{test_cache_all2lmp}" 

{"# file-tag":<15}{"filename":<40}{"comment (required)"}  
"""
    
    # Add entries for each molecule file
    for name, lammps_file in molecule_files_all2lmp.items():
        # Determine appropriate comment based on file type
        if name.startswith("data"):
            # Data files contain all force field coefficients
            comment = "# This datafile will have all coeffs in it"
        else:
            # Pre- and post-reaction files are associated with specific reactions
            rxn_number = get_ending_integer(name)
            comment = f"# for rxn{rxn_number}"
        
        # Add formatted line to the merge file
        merge_files += f"{name:<15}{lammps_file:<40}{comment} file\n"
    
    # Add footer with output directory specification
    merge_files += f""" 
# Specify the parent_directory of where to write results (optional) 
parent_directory = "{test_cache_bond_react_merge}"
"""
    
    # Ensure the output directory exists
    try:
        os.makedirs(test_cache_bond_react_merge, exist_ok=True)
    except FileExistsError:
        # Directory already exists, continue without error
        pass

    # Write the merge input file
    merge_input_path = os.path.join(test_cache_bond_react_merge, "merge_input.txt")
    with open(merge_input_path, "w") as f:
        f.write(merge_files)
    
    time.sleep(1)  # Ensure the merge input file is fully written before returning
    # Return the directory containing the merge input file
    return Path(test_cache_bond_react_merge)

def run_bond_react_merge(merge_input_file_path: Path, 
                         molecule_files_all2lmp: dict[str, str]) -> dict[str, Path]:
    """
    Execute LUNAR's bond_react_merge utility to combine and process LAMMPS data files.
    
    This function merges multiple LAMMPS data files into a unified format suitable
    for reactive molecular dynamics simulations, generating both merged data files
    and a force field parameter file.
    
    Args:
        merge_input_file_path (Path): Directory containing the merge input file.
        molecule_files_all2lmp (dict[str, str]): Dictionary mapping molecule names 
                                                to LAMMPS data filenames.
    
    Returns:
        dict[str, Path]: Dictionary mapping molecule names and special keys to 
                        their output file paths:
                        - "dataX": Merged data files for data molecules
                        - "preX": LAMMPS molecule files for pre-reaction states
                        - "postX": LAMMPS molecule files for post-reaction states
                        - "force_field_data": Combined force field parameter file
    
    Raises:
        FileNotFoundError: If expected output files are not generated.
    
    Note:
        Uses the following LUNAR bond_react_merge.py command:
        python3 bond_react_merge.py -files infile:merge_input.txt -atomstyle full
    """
    molecule_files_bond_react_merge: dict[str, Path] = {}
    
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
                merge_input_file_path / f"{name}_typed_IFF_merged.data"
            )
        elif name.startswith("pre"):
            # Pre-reaction files generate .lmpmol files
            rxn_number = get_ending_integer(name)
            molecule_files_bond_react_merge[name] = (
                merge_input_file_path / f"pre{rxn_number}_typed_IFF_merged.lmpmol"
            )
        elif name.startswith("post"):
            # Post-reaction files generate .lmpmol files
            rxn_number = get_ending_integer(name)
            molecule_files_bond_react_merge[name] = (
                merge_input_file_path / f"post{rxn_number}_typed_IFF_merged.lmpmol"
            )
    
    # Add the force field data file (generated by bond_react_merge)
    molecule_files_bond_react_merge["force_field_data"] = (
        merge_input_file_path / "force_field.data"
    )
    
    # Verify that all expected output files were created
    for name, path in molecule_files_bond_react_merge.items():
        if path.is_file():
            print(f"[LUNAR bond_react_merge] Generated file for {name} at: {path}")
        else:
            raise FileNotFoundError(
                f"Expected output file for {name} not found at {path}"
            )
    
    return molecule_files_bond_react_merge


def lunar_workflow(molecule_files: dict[str, Path]) -> dict[str, Path]:
    """
    Execute the complete LUNAR workflow for preparing molecules for LAMMPS simulations.
    
    This function orchestrates the entire process:
    1. Atom typing using PCFF-IFF force field
    2. Conversion to LAMMPS data format
    3. Generation of all topological requirements 
    4. Merging of files for reactive simulations
    
    Args:
        molecule_files (dict[str, Path]): Dictionary mapping molecule names to 
                                         their .mol file paths.
    
    Returns:
        dict[str, Path]: Dictionary mapping molecule names and special keys to 
                        their final output file paths after the complete workflow.
    
    Workflow Steps:
        0. Display loading screen
        1. Run LUNAR atom typing
        2. Run all2lmp conversion
        3. Write bond_react_merge input file
        4. Run bond_react_merge utility
    """
    # Step 0: Show loading screen
    loading_screen("LUNAR Workflow")
    
    # Step 1: Run atom typing
    typed_files = run_LUNAR_atom_typing(molecule_files)
    
    # Step 2: Run all2lmp conversion
    all2lmp_files = run_LUNAR_all2lmp(typed_files)
    
    # Step 3: Write bond_react_merge input file
    merge_input_path = write_bond_react_merge_input(all2lmp_files)
    
    # Step 4: Run bond_react_merge
    final_files = run_bond_react_merge(merge_input_path, all2lmp_files)
    
    return final_files
    
if __name__ == "__main__":
    """
    Example usage of the LUNAR API wrapper.
    
    This block demonstrates how to use the lunar_workflow function with
    sample molecule files. In practice, users should modify the molecule_files
    dictionary to point to their actual .mol files.
    """
    # Import Path here to avoid circular imports if this module is imported elsewhere
    from pathlib import Path
    
    # Define sample molecule files for testing
    # Note: These paths are specific to the original developer's system
    molecule_files: dict[str, Path] = {
        "data1": Path(r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\lunar\data1.mol"),
        "data2": Path(r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\lunar\data2.mol"),
        "pre1":  Path(r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\lunar\pre1.mol"),
        "post1": Path(r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\lunar\post1.mol"),
    }
    
    # Execute the complete LUNAR workflow
    final_files = lunar_workflow(molecule_files)
    
    # Output summary of generated files
    print("\n" + "=" * 80)
    print("LUNAR Workflow Complete!")
    print("Generated files:")
    for name, path in final_files.items():
        print(f"  {name}: {path}")
    print("=" * 80)