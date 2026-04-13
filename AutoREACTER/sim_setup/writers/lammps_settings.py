from dataclasses import dataclass
from pathlib import Paths

from cairo import Path
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

class InitialSetupWriter:
    def __init__(self, reacter_files: REACTERFiles):
        self.reacter_files = reacter_files
        lunar_in_file_location = reacter_files.in_file

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
            special_bonds="lj/coul 0.0 0.0 0.5",
            pair_style="lj/cut/coul/long 10.0 10.0",
            kspace_style="pppm 1e-4",
            pair_modify="mix arithmetic",
            neighbor="2.0 bin",
            neigh_modify="delay 5 every 1"
        )

    def _get_LUNAR_lammps_settings(self, output_path: Path) -> None:
        with open(output_path, 'w') as f:
            pass

