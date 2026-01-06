"""Main API wrapper for interacting with LUNAR."""
import subprocess
import os, sys
import subprocess
from pathlib import Path
import shutil
import datetime
import time
from rdkit import Chem
from rdkit.Chem import AllChem

test_cache = r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\lunar"
test_cache_atom_typing = r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\lunar\atom_typing"
test_cache_all2lmp = r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\lunar\all2lmp"
test_cache_bond_react_merge = r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\lunar\bond_react_merge"

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
    for name, file in molecule_files.items():
        print(f"Processing {name} from {file}")
        subprocess.run([sys.executable, atom_typing_py, '-topo', str(file),  '-dir', test_cache_atom_typing, "-ff", "PCFF-IFF", "-del-method", "mass", "-del-crit", "0"])
