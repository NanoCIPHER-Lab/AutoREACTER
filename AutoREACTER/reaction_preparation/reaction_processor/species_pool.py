from __future__ import annotations

from dataclasses import dataclass, field
from rdkit import Chem
from typing import Optional

@dataclass(slots=True)
class PoolSpecies:
    """
    Represents a species in the reaction pool.
    Attributes:
        species_id (str): The unique identifier of the species.
        is_monomer (bool): Indicates whether the species is a monomer.
        smiles (Optional[str]): The SMILES representation of the species.
        mol (Optional[Chem.Mol]): The RDKit Mol object of the species.
        template_idxes (Optional[list[int]]): The list of atom indices used as a template.
    """
    species_id: str
    is_monomer: bool
    smiles: Optional[str] = None
    mol: Optional[Chem.Mol] = None
    template_idxes: Optional[list[int]] = None


def _populate_mols(pool: list[PoolSpecies]) -> list[PoolSpecies]:
    """
    Populate the `mol` and `template_idxes` attributes of each monomer species in the pool based on its SMILES string.
    Args:
        pool (list[PoolSpecies]): The list of species in the pool.

    Returns:
        list[PoolSpecies]: The updated list of species with populated `mol` and `template_idxes` attributes.
    """
    for species in pool:
        if species.smiles and species.mol is None and species.is_monomer:
            species.mol = Chem.MolFromSmiles(species.smiles)
            atoms = []
            for atom in species.mol.GetAtoms():
                atoms.append(atom.GetIdx())
            species.template_idxes = atoms
    return pool
