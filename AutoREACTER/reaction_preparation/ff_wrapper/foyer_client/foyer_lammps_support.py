import os
import re
from pathlib import Path
from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional
from lammps import lammps

if TYPE_CHECKING:
    from AutoREACTER.reaction_preparation.lunar_client.lunar_api_wrapper import AtomTypingResult

class GetForceFieldDataFile:
    def __init__(self, atom_typing_results: list['AtomTypingResult'], output_dir: Optional[Path] = None):
        self.atom_typing_results = atom_typing_results
        self.output_dir = output_dir
        self.force_field_data_file = self._parse_abstract_force_field_data_file(atom_typing_results)


    def _types_counter(self, atom_typing_results: list['AtomTypingResult']) -> dict[str, int]:
        """Count the number of unique occurrences of each atom type across LUNAR-processed data files."""
        unique_types = {
            "atom_types": set(),
            "bond_types": set(),
            "angle_types": set(),
            "dihedral_types": set(),
            "improper_types": set()
        }
        
        keys = [
            "Atom Type Labels", 
            "Bond Type Labels", 
            "Angle Type Labels", 
            "Dihedral Type Labels", 
            "Improper Type Labels"
        ]

        for result in atom_typing_results:
            data_file = result.typed_data_file
            with open(data_file, 'r') as f:
                content = f.read()
                for key in keys:
                    pattern = rf"{key}:\s*(\d+)\s*([^\n]+)"
                    match = re.search(pattern, content)
                    
                    if match:
                        dict_key = key.split()[0].lower() + "_types"
                        types_line = match.group(2)
                        unique_types[dict_key].update(types_line.split())

        # Convert the sets to integer counts to satisfy dict[str, int]
        return {k: len(v) for k, v in unique_types.items()}

    def _generate_empty_data_file(self, output_dir: Optional[Path]) -> str:
        """Helper to create a minimal LAMMPS data file to establish the simulation box."""
        filename = "empty.data"
        minimal_empty_data = """LAMMPS minimal empty data file

0 atoms
0 bonds
0 angles
0 dihedrals
0 impropers

0 atom types
0 bond types
0 angle types
0 dihedral types
0 improper types

-50.0 50.0 xlo xhi
-50.0 50.0 ylo yhi
-50.0 50.0 zlo zhi
"""
        if output_dir:
            filename = output_dir / "empty.data"
        else:
            filename = "empty.data"
        with open(filename, "w") as f:
            f.write(minimal_empty_data)
        
        return filename


    def _parse_abstract_force_field_data_file(
        self, atom_typing_results: list['AtomTypingResult']
    ) -> Path:
        """
        Dynamically counts types, allocates space in LAMMPS, appends all molecule 
        files, and outputs the combined force_field.data file.
        """
        
        # 1. Get dynamic counts from all files
        types_count = self._types_counter(atom_typing_results)
        
        # 2. Create the dummy file required by LAMMPS to initialize the box
        empty_data_path = self._generate_empty_data_file(self.output_dir)

        # 3. Initialize LAMMPS (suppress output log)
        lmp = lammps(cmdargs=['-log', 'none'])
        lmp.command("atom_style full")

        # 4. Construct read_data command with dynamic values from our counter
        read_empty_cmd = (
            f"read_data {empty_data_path} "
            f"extra/atom/types {types_count.get('atom_types', 0)} "
            f"extra/bond/types {types_count.get('bond_types', 0)} "
            f"extra/angle/types {types_count.get('angle_types', 0)} "
            f"extra/dihedral/types {types_count.get('dihedral_types', 0)} "
            f"extra/improper/types {types_count.get('improper_types', 0)} "
            "extra/bond/per/atom 50 "
            "extra/angle/per/atom 50 "
            "extra/dihedral/per/atom 50 "
            "extra/improper/per/atom 50 "
            "extra/special/per/atom 50"
        )
        
        # Allocate the box
        lmp.command(read_empty_cmd)

        # 5. Append every file into the LAMMPS system
        for result in atom_typing_results:
            # "add append" safely merges the new atoms into the existing system
            lmp.command(f"read_data {result.typed_data_file} add append")

        # 6. Write out the final combined file
        if self.output_dir:
            output_file = self.output_dir / "force_field.data"
        else:
            output_file = Path("force_field.data").resolve()
        lmp.command(f"write_data {output_file.as_posix()}")

        # Clean up the dummy file
        if os.path.exists(empty_data_path):
            os.remove(empty_data_path)

        return output_file
    

class DATAFile2LAMMPSMolecule:
    def __init__(self, data_file: Path, units: str = "real", atom_style: str = "full"):
        """
        Initializes the LAMMPS instance and loads the data file.
        Units and atom_style are parameterized to match the data file.
        """
        self.data_file = data_file
        self.units = units
        self.atom_style = atom_style
        
        # Keep the active LAMMPS instance as an attribute
        self.lmp = self._load_data_file(self.data_file)

    def _load_data_file(self, data_file: Path) -> lammps:
        """Loads a LAMMPS data file into the LAMMPS instance."""
        lmp = lammps(cmdargs=['-log', 'none', '-echo', 'none'])
        lmp.command(f"units {self.units}")
        lmp.command(f"atom_style {self.atom_style}")
        lmp.command(f"read_data {data_file.as_posix()}")
        return lmp

    def extract_and_write(self, molecule_id: int = 1, output_file: Path = None) -> Path:
        """
        Groups a specific molecule by its ID and writes it to a LAMMPS molecule file.
        If no output_file is specified, it automatically creates one by changing 
        the input file's extension from .data to .molecule.
        
        Returns the Path object of the generated molecule file.
        """
        # Automatically handle the .data -> .molecule conversion if no path is given
        if output_file is None:
            output_file = self.data_file.with_suffix('.molecule')
            
        group_name = f"single_mol_{molecule_id}"
        
        try:
            # Group the specific molecule ID
            self.lmp.command(f"group {group_name} molecule {molecule_id}")
            
            # Write out to native LAMMPS molecule format
            self.lmp.command(f"write_molecule {group_name} {output_file.as_posix()}")
            
            # Return the path so it can be passed directly to your dataclass
            return output_file
            
        except Exception as e:
            print(f"An error occurred while writing the molecule: {e}")
            raise e






            