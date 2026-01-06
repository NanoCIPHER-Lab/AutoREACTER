"""Main API wrapper for interacting with LUNAR.
Provides functions to run LUNAR atom typing, all2lmp conversion, and bond_react_merge utilities.
This module can handle arbitrary molecule files and manage temporary cache directories."""
import subprocess
import os, sys
import subprocess
from pathlib import Path
import shutil
import datetime
import time
from rdkit import Chem
from rdkit.Chem import AllChem
import re

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

def run_LUNAR_atom_typing(input_file: str, output_file: str):
    """Run LUNAR atom typing on the input file and save to output file."""
    command = [LUNAR_LOCATION, '-i', input_file, '-o', output_file]
    subprocess.run(command, check=True)

def run_LUNAR_all2lmp(input_file: str, output_file: str):
    """Run LUNAR all2lmp conversion on the input file and save to output file."""
    command = [LUNAR_LOCATION, '-i', input_file, '-o', output_file, '--all2lmp']
    subprocess.run(command, check=True)

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
    
if __name__ == "__main__":
    from pathlib import Path

    molecule_files: dict[str, Path] = {
        "data1": Path(r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\lunar\data1.mol"),
        "data2": Path(r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\lunar\data2.mol"),
        "pre1":  Path(r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\lunar\pre1.mol"),
        "post1": Path(r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\lunar\post1.mol"),
    }

    atom_typing_py = os.path.join(LUNAR_LOCATION, "atom_typing.py")
    all2lmp_py = os.path.join(LUNAR_LOCATION, "all2lmp.py")

    bond_react_merge_py = os.path.join(LUNAR_LOCATION, "bond_react_merge.py")
    molecule_files_typed: dict[str, Path] = {}
    for name, file in molecule_files.items():
        print(f"Processing {name} from {file}")
        "python3   atom_typing.py    -topo   test1.mol   -dir   testing_directory    -ff   PCFF-IFF   -del-method   mass   -del-crit   0 "
        subprocess.run([sys.executable, atom_typing_py, '-topo', str(file),  '-dir', test_cache_atom_typing, "-ff", "PCFF-IFF", "-del-method", "mass", "-del-crit", "0"])
        molecule_files_typed[name] = [Path(test_cache_atom_typing) / f"{file.stem}_typed.data", Path(test_cache_atom_typing) / f"{file.stem}_typed.nta"]

    # sometime LUNAR can't identify the .frc file location automatically so we give the full path
    frc_file = os.path.join(LUNAR_LOCATION, "frc_files", "pcff.frc")
    molecule_files_all2lmp: dict[str, Path] = {}
    for name, (data_file, nta_file) in molecule_files_typed.items():
        print(f"Converting {name} typed files to LAMMPS format")
        "python3   all2lmp.py    -topo   test1.data   -nta  test1.nta   -class  2  -frc  pcff.frc   -dir  testing_directory "
        subprocess.run([sys.executable, all2lmp_py, '-topo', str(data_file), '-nta', str(nta_file), '-frc', frc_file, '-dir', test_cache_all2lmp])
        molecule_files_all2lmp[name] = f"{data_file.stem}_IFF.data"

    # Example of bond_react_merge usage
    sub_run_bond_react_merge = []
    merge_files:str = f"""# anything following the “#” character will be ignored 
 
# Specify a desired path to append to the front of each filename (optional however if not present 
# and the files are a path from LUNAR on your computer each file below must have that path 
# specified in front of the filename) 

path = f"{test_cache_all2lmp}" 

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
{test_cache_bond_react_merge}
"""
    try:
        os.makedirs(test_cache_bond_react_merge, exist_ok=True)
    except FileExistsError:
        pass
    with open(os.path.join(test_cache_bond_react_merge, "merge_input.txt"), "w") as f:
        f.write(merge_files)

        



    
    