"""Main API wrapper for interacting with LUNAR.
Provides functions to run LUNAR atom typing, all2lmp conversion, and bond_react_merge utilities.
This module can handle arbitrary molecule files and manage temporary cache directories."""
import subprocess
import os, sys
import subprocess
from pathlib import Path
import re
import time

test_cache = r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\lunar"
test_cache_atom_typing = f"{test_cache}\\atom_typing"
test_cache_all2lmp = f"{test_cache}\\all2lmp"
test_cache_bond_react_merge = f"{test_cache}\\bond_react_merge"
try:
    os.makedirs(test_cache)
    os.makedirs(test_cache_atom_typing)
    os.makedirs(test_cache_all2lmp) 
    os.makedirs(test_cache_bond_react_merge)
except FileExistsError:
    pass
os.makedirs(test_cache, exist_ok=True)
print(f"Test cache folder: {test_cache}")


from locate_LUNAR import get_LUNAR_loc
cwd = os.getcwd()
LUNAR_LOCATION = get_LUNAR_loc(use_gui=False)

atom_typing_py = os.path.join(LUNAR_LOCATION, "atom_typing.py")
all2lmp_py = os.path.join(LUNAR_LOCATION, "all2lmp.py")
bond_react_merge_py = os.path.join(LUNAR_LOCATION, "bond_react_merge.py")

def loading_screen(name="LUNAR") -> None:
        """Show ASCII art banner with loading spinner."""
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
        print(banner)
        print(f"Loading {name} ", end="", flush=True)

        animation = ["-", "\\", "|", "/"]
        for i in range(10):
            time.sleep(0.1)
            sys.stdout.write("\b" + animation[i % len(animation)])
            sys.stdout.flush()
        print("\nReady!")


def get_ending_integer(s):
    """
    Extracts and converts the trailing integer from a string.
    Returns the integer if found, otherwise None.
    """
    # The pattern r'\d+$' matches one or more digits (\d+) at the end of the string ($)
    match = re.search(r'\d+$', s)
    if match:
        # Get the matched substring and convert it to an integer
        return int(match.group())
    else:
        return None
def run_LUNAR_atom_typing(molecule_files: dict[str, Path]) -> dict[str, Path]:
    molecule_files_typed: dict[str, Path] = {}
    for name, file in molecule_files.items():
    #    "python3   atom_typing.py    -topo   test1.mol   -dir   testing_directory    -ff   PCFF-IFF   -del-method   mass   -del-crit   0 "
        subprocess.run([
            sys.executable, 
            atom_typing_py, 
            '-topo', str(file),  
            '-dir', test_cache_atom_typing, 
            "-ff", "PCFF-IFF", 
            "-del-method", "mass", 
            "-del-crit", "0"
            ])
        molecule_files_typed[name] = [Path(test_cache_atom_typing) / f"{file.stem}_typed.data", Path(test_cache_atom_typing) / f"{file.stem}_typed.nta"]
    return molecule_files_typed

def run_LUNAR_all2lmp(molecule_files_typed: dict[str, Path]) -> dict[str, Path]:
    # sometime LUNAR can't identify the .frc file location automatically so we give the full path
    frc_file = os.path.join(LUNAR_LOCATION, "frc_files", "pcff.frc")
    molecule_files_all2lmp: dict[str, Path] = {}
    for name, (data_file, nta_file) in molecule_files_typed.items():
    #    "python3   all2lmp.py    -topo   test1.data   -nta  test1.nta   -class  2  -frc  pcff.frc   -dir  testing_directory "
        subprocess.run([
            sys.executable, 
            all2lmp_py, 
            '-topo', str(data_file), 
            '-nta', str(nta_file), 
            '-frc', frc_file, 
            '-dir', test_cache_all2lmp
            ])
        molecule_files_all2lmp[name] = f"{data_file.stem}_IFF.data"
    return molecule_files_all2lmp

def write_bond_react_merge_input(molecule_files_all2lmp: dict[str, Path]) -> Path:
    # Example of bond_react_merge usage
    merge_files:str = f"""# anything following the “#” character will be ignored 
 
# Specify a desired path to append to the front of each filename (optional however if not present 
# and the files are a path from LUNAR on your computer each file below must have that path 
# specified in front of the filename) 

path = "{test_cache_all2lmp}" 

{"# file-tag":<15}{"filename":<40}{"comment (required)"}  
"""
    for name, lammps_file in molecule_files_all2lmp.items():
        if name.startswith("data"):
            comment = "# This datafile will have all coeffs in it"
        else:
            rxn_number = get_ending_integer(name)
            comment = f"# for rxn{rxn_number}" 
        merge_files += f"{name:<15}{lammps_file:<40}{comment} file\n"
    merge_files +=f""" 
# Specify the parent_directory of where to write results (optional) 
parent_directory = "{test_cache_bond_react_merge}"
"""
    try:
        os.makedirs(test_cache_bond_react_merge, exist_ok=True)
    except FileExistsError:
        pass
    with open(os.path.join(test_cache_bond_react_merge, "merge_input.txt"), "w") as f:
        f.write(merge_files)
    return test_cache_bond_react_merge

def run_bond_react_merge(merge_input_file_path: Path, molecule_files_all2lmp: dict[str, Path]) -> None:
    # To run bond_react_merge, uncomment the following lines:

    molecule_files_bond_react_merge: dict[str, Path] = {}
    # python3   bond_react_merge.py    -files  infile:merge_files.txt   -atomstyle full   -class 2
    subprocess.run(
        [
            sys.executable,
            bond_react_merge_py,
            "-files",
            f"infile:merge_input.txt",
            "-atomstyle",
            "full",
        ],
        cwd=merge_input_file_path,
        check=True,
    )
    for name in molecule_files_all2lmp.keys():
        if name.startswith("data"):
            molecule_files_bond_react_merge[name] = Path(merge_input_file_path) / f"{name}_typed_IFF_merged.data"
        if name.startswith("pre"):
            rxn_number = get_ending_integer(name)
            molecule_files_bond_react_merge[name] = Path(merge_input_file_path) / f"pre{rxn_number}_typed_IFF_merged.lmpmol"
        if name.startswith("post"):
            rxn_number = get_ending_integer(name)
            molecule_files_bond_react_merge[name] = Path(merge_input_file_path) / f"post{rxn_number}_typed_IFF_merged.lmpmol"
    molecule_files_bond_react_merge["force_field_data"] = Path(merge_input_file_path) / "force_field.data"
    for name, path in molecule_files_bond_react_merge.items():
        if path.is_file():
            print(f"[LUNAR bond_react_merge] Generated file for {name} at: {path}")
            pass
        else:
            raise FileNotFoundError(f"Expected output file for {name} not found at {path}")
    return molecule_files_bond_react_merge

def lunar_workflow(molecule_files: dict[str, Path]) -> dict[str, Path]:
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
    from pathlib import Path

    molecule_files: dict[str, Path] = {
        "data1": Path(r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\lunar\data1.mol"),
        "data2": Path(r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\lunar\data2.mol"),
        "pre1":  Path(r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\lunar\pre1.mol"),
        "post1": Path(r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\lunar\post1.mol"),
    }
    final_files = lunar_workflow(molecule_files)





        



    
    