"""
3D Molecule Preparation and Optimization Utilities

This module provides utilities to:
- Separate disconnected fragments in a molecule by translating one fragment in 3D.
- Embed and optimize 3D geometries for molecules (using RDKit).
- Prepare a dictionary of molecules by creating optimized 3D .mol files saved to a cache directory.

Notes / TODO:
- The `cache` path is currently hard-coded and should be obtained from a configuration/cache module
  (e.g., a `cache.py`) in a production setting.
- This code currently assumes simple cases (single mapping, one reactant/product complex per entry).
  It does not handle multiple alternative reaction products or complex multi-reactant reaction topologies.
"""

import os
from pathlib import Path
from pathlib import Path
from typing import Tuple
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol
from AutoREACTER.input_parser import SimulationSetup
from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import ReactionMetadata

class Molecule3DPreparationError(Exception):
    """Custom exception for errors during 3D molecule preparation."""

class FragmentSeparationError(Molecule3DPreparationError):
    """Raised when there is an issue with separating fragments in 3D."""

class OptimizationError(Molecule3DPreparationError):
    """Raised when there is an issue during the optimization process."""

class Molecule3DPreparation:
    """Utility class for preparing and optimizing 3D molecule geometries."""

    def __init__(self, cache_dir: Path):
        self.cache_dir = Path (cache_dir / "3D_molecules")
        os.makedirs(self.cache_dir, exist_ok=True)
        self.molecule_3d_path = Path(self.cache_dir / "molecules_3Dmol")
        os.makedirs(self.molecule_3d_path, exist_ok=True)
        self.full_templates_path = Path(self.cache_dir / "full_templates_3Dmol")
        os.makedirs(self.full_templates_path, exist_ok=True)

    @property
    def cache(self):
        return self.cache_dir

    def prepare_3d_molecules(self, updated_inputs: SimulationSetup, prepared_reactions: list[ReactionMetadata]) -> tuple[SimulationSetup, list[ReactionMetadata]]:
        """
        Prepare and optimize a dictionary of molecules by embedding and minimizing 3D geometries.

        For each (name, mol) entry in `molecule_dict`, `optimization` is called and the
        resulting file path (string) replaces the Mol object in the returned dictionary.

        Parameters
        ----------
        cache_dir : str
            Directory to save optimized .mol files to.
        molecule_dict : dict of str -> rdkit.Chem.rdchem.Mol
            Mapping of molecule names to RDKit Mol objects to be prepared.

        Returns
        -------
        dict of str -> str
            Mapping of molecule names to file paths of the saved optimized .mol files.
        """
        monomer_entries = updated_inputs.monomers
        for entry in monomer_entries:
            rdkit_mol = entry.rdkit_mol
            if rdkit_mol is None:
                raise Molecule3DPreparationError(f"RDKit Mol object is missing for molecule {entry.data_id}. Cannot prepare 3D geometry without it.")
            molecule_name = entry.data_id
            try:
                optimized_path, mol_block = self._optimization(molecule_name, rdkit_mol, self.molecule_3d_path)
                entry.molecule_3Dmol_path = optimized_path
                entry.molecule_3Dmol_file = mol_block
            except Exception as e:
                raise OptimizationError(f"Error optimizing molecule {molecule_name}: {str(e)}") from e
        
        for reaction in prepared_reactions:
            combined_reactants = reaction.combined_reactants
            combined_products = reaction.combined_products

            if combined_reactants is not None:
                try:
                    optimized_reactants_path, reactants_mol_block = self._optimization(f"{reaction.reaction_id}_pre", combined_reactants, self.full_templates_path)
                    reaction.reactants_3Dmol_path =  optimized_reactants_path
                    reaction.reactants_3Dmol_file = reactants_mol_block
                except Exception as e:
                    raise OptimizationError(f"Error optimizing reactant complex for reaction {reaction.reaction_id}: {str(e)}") from e
            
            if combined_products is not None:
                try:
                    optimized_products_path, products_mol_block = self._optimization(f"{reaction.reaction_id}_post", combined_products, self.full_templates_path)
                    reaction.products_3Dmol_path =  optimized_products_path
                    reaction.products_3Dmol_file = products_mol_block
                except Exception as e:
                    raise OptimizationError(f"Error optimizing product complex for reaction {reaction.reaction_id}: {str(e)}") from e


    def _separate_fragments_3d(self, mol: Mol, shift: Tuple[float, float, float] = (8.0, 0.0, 0.0)) -> Mol:
        # Get tuple of atom index lists for each fragment (asMols=False returns tuples of atom indices)
        frags = Chem.GetMolFrags(mol, asMols=False)

        # If there's only one fragment, nothing to do
        if len(frags) < 2:
            return mol

        # For now we only support exactly two fragments. Raise if more are present to avoid silent errors.
        if len(frags) > 2:
            raise FragmentSeparationError(f"Expected 2 fragments for separation, but found {len(frags)}. Please check the molecule structure.")
        # Get the conformer (3D coordinates) for the molecule
        conf = mol.GetConformer()
        shift_arr = np.array(shift, dtype=float)

        # Translate atoms in the second fragment (fragment index 1)
        for atom_idx in frags[1]:
            pos = np.array(conf.GetAtomPosition(atom_idx))
            conf.SetAtomPosition(atom_idx, pos + shift_arr)

        return mol


    def _optimization(self, molecule_name: str, mol: Mol, cache_dir: Path, separate_fragments: bool = False) -> Tuple[Path, str]:
        
        # Record number of explicit atoms before doing anything (sanity check later)
        n_atoms_start = mol.GetNumAtoms(onlyExplicit=True)

        # Update RDKit's property cache (e.g., valence information). strict=False avoids exceptions
        # for molecules that are slightly invalid but salvageable.
        mol.UpdatePropertyCache(strict=False)

        # Sanitize the molecule. We skip SANITIZE_PROPERTIES because writing or reading various
        # properties can be fragile; we keep the core sanitization operations.
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
        )

        # Embed the molecule into 3D space using ETKDG (a robust, knowledge-based method).
        # Embedding produces a conformer on the molecule that we can then optimize.
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())

        # If the molecule name indicates reactant/product complexes (e.g., starts with 'pre'/'post'),
        # separate disconnected fragments in 3D to avoid overlapping coordinates.
        if separate_fragments:
            mol = self._separate_fragments_3d(mol, shift=(8.0, 0.0, 0.0))

        # Optimize geometry using the MMFF force field (works well for many organic molecules).
        # This modifies the conformation in-place.
        AllChem.MMFFOptimizeMolecule(mol)

        # Integrity check: make sure explicit atom count didn't change during operations.
        # A mismatch usually indicates an unintended removal or addition of atoms.
        if n_atoms_start != mol.GetNumAtoms(onlyExplicit=True):
            raise OptimizationError(f"Atom count mismatch for {molecule_name}: started with {n_atoms_start} explicit atoms but ended with {mol.GetNumAtoms(onlyExplicit=True)}. This may indicate an issue during sanitization or embedding.")

        # Create output directory if it doesn't exist
        os.makedirs(cache_dir, exist_ok=True)

        # Save optimized geometry as a .mol file
        output_path = os.path.join(cache_dir, f"{molecule_name}.mol")
        print(f"Saving optimized {molecule_name} to {output_path}")
        Chem.MolToMolFile(mol, output_path)
        mol_block = Chem.MolToMolBlock(mol)

        return output_path, mol_block