from dataclasses import dataclass
from pathlib import Paths
from AutoREACTER.reaction_preparation.lunar_client.REACTER_files_builder import REACTERFiles

@dataclass(slots=True)
class LammpsSettings:
    units: str = "real"
    dimension: str = "3"
    boundary: str = "p p p"
    atom_style: str = "full"
    newton: str = "on"
    bond_style: str = "class2"
    angle_style: str = "class2"
    dihedral_style: str = "class2"
    improper_style: str = "class2"
    special_bonds: str = "lj/coul 0 0 1"
    pair_style: str = "lj/class2/coul/long 8.5"
    kspace_style: str = "pppm 1.0e-4"
    pair_modify: str = "tail yes mix sixthpower"
    neighbor: str = "2.0 bin"
    neigh_modify: str = "delay 0 every 1 check yes one 5000 page 100000"

class InitialSetupWriter:
    def __init__(self, reacter_files: REACTERFiles):
        self.reacter_files = reacter_files
