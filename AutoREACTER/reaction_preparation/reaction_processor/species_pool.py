from __future__ import annotations

from dataclasses import dataclass, field
from rdkit import Chem
from typing import Optional

@dataclass(slots=True)
class PoolSpecies:
    """
    Represents a species in the reaction pool.
    """
    species_id: str
    is_monomer: bool
    smiles: Optional[str] = None
    mol: Optional[Chem.Mol] = None
    template_idxes: Optional[list[int]] = None


def _populate_mols(pool: list[PoolSpecies]) -> list[PoolSpecies]:
    """
    Populate the `mol` attribute of each monomer species in the pool based on its SMILES string.
    Args:
        pool (list[PoolSpecies]): The list of species in the pool.

    Returns:
        list[PoolSpecies]: The updated list of species with populated `mol` attributes.
    """
    for species in pool:
        if species.smiles and species.mol is None and species.is_monomer:
            species.mol = Chem.MolFromSmiles(species.smiles)
            atoms = []
            for atom in species.mol.GetAtoms():
                atoms.append(atom.GetIdx())
            species.template_idxes = atoms



    return pool
