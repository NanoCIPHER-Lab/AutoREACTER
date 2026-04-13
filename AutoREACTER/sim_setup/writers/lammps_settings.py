from dataclasses import dataclass, fields

from AutoREACTER.reaction_preparation.lunar_client.REACTER_files_builder import REACTERFiles

@dataclass(slots=True)
class LammpsSettings:
    units: str 
    dimension: str 
    boundary: str 
    atom_style: str 
    bond_style: str 
    angle_style: str
    dihedral_style: str
    improper_style: str
    special_bonds: str
    pair_style: str 
    kspace_style: str
    pair_modify: str 
    neighbor: str
    neigh_modify: str 

class LammpsInitialSettings:
    def __init__(self, reacter_files: REACTERFiles):
        self.reacter_files = reacter_files
        self.lunar_in_file_location = reacter_files.in_file


    def _defults(self) -> LammpsSettings:
        return LammpsSettings(
            units="real",
            dimension="3",
            boundary="p p p",
            atom_style="full",
            bond_style="harmonic",
            angle_style="harmonic",
            dihedral_style="opls",
            improper_style="cvff",
            special_bonds="lj/coul 0 0 1",
            pair_style="lj/class2/coul/long 8.5",
            kspace_style="pppm 1e-4",
            pair_modify="tail yes mix sixthpower",
            neighbor=None,
            neigh_modify=None
        )

    def get_LUNAR_lammps_settings(self) -> None:
        settings = self._defults()
        try:
            with open(self.lunar_in_file_location, 'r') as f:
                lines = f.readlines()
        except Exception as e:
            print(f"[WARNING] Could not read LAMMPS input file at {self.lunar_in_file_location}: {e}")
            print("[WARNING] Proceeding with default LAMMPS settings.")
        
        for field in fields(settings):
            keyword = field.name.replace('_', ' ')
            for line in lines:
                clean_line = line.split('#')[0].strip()
                if clean_line.startswith(keyword):
                    value = clean_line[len(keyword):].strip()
                    if value:
                        setattr(settings, field.name, value)
        return settings
        
        
