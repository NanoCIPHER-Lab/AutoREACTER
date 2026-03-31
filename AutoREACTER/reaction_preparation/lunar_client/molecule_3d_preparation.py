"""
3D Molecule Preparation and Optimization Utilities

This module provides utilities for preparing 3D molecular geometries using RDKit.
It handles:
- Separation of disconnected molecular fragments in 3D space
- 3D coordinate generation (embedding) and geometry optimization
- Saving optimized structures as .mol files to a cache directory

The main class `Molecule3DPreparation` is used to prepare both individual monomers
and combined reactant/product complexes for simulations.
"""

# Standard library imports
import os
from pathlib import Path
from typing import Tuple

# Third-party imports
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import Descriptors

# Local imports from the AutoREACTER package
from AutoREACTER.input_parser import SimulationSetup
from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import ReactionMetadata


class Molecule3DPreparationError(Exception):
    """Custom exception for errors during 3D molecule preparation."""


class FragmentSeparationError(Molecule3DPreparationError):
    """Raised when there is an issue with separating molecular fragments in 3D."""


class OptimizationError(Molecule3DPreparationError):
    """Raised when an error occurs during 3D embedding or geometry optimization."""


class Molecule3DPreparation:
    """Utility class for preparing and optimizing 3D molecule geometries.
    
    This class handles the generation of 3D coordinates for monomers and
    reaction complexes, performs geometry optimization, and saves the results
    to disk for later use in simulations.
    
    Attributes:
        cache_dir: Base directory where 3D molecule files will be stored.
        molecule_3d_path: Directory for individual monomer 3D structures.
        full_templates_path: Directory for combined reactant/product complexes.
    """

    def __init__(self, cache_dir: Path):
        """Initialize the 3D preparation utility and create required directories.

        Args:
            cache_dir: Base directory where 3D molecule files will be stored.
        """
        self.cache_dir = Path(cache_dir / "3D_molecules")
        os.makedirs(self.cache_dir, exist_ok=True)

        self.molecule_3d_path = Path(self.cache_dir / "molecules_3Dmol")
        os.makedirs(self.molecule_3d_path, exist_ok=True)

        self.full_templates_path = Path(self.cache_dir / "full_templates_3Dmol")
        os.makedirs(self.full_templates_path, exist_ok=True)

    @property
    def cache(self) -> Path:
        """Return the main cache directory for 3D molecules.
        
        Returns:
            Path object pointing to the cache directory.
        """
        return self.cache_dir

    def prepare_molecule_3d_geometry(
        self,
        updated_inputs: SimulationSetup,
        prepared_reactions: list[ReactionMetadata],
    ) -> tuple[SimulationSetup, list[ReactionMetadata]]:
        """Prepare 3D geometries for all monomers and reaction complexes.

        This method processes individual monomer molecules and the combined
        reactant/product complexes for each reaction, optimizing their
        3D structures and attaching the resulting file paths and mol blocks
        to the respective data objects.

        Args:
            updated_inputs: Simulation setup containing monomer information.
            prepared_reactions: List of reaction metadata objects.

        Returns:
            Updated simulation inputs and reaction metadata with 3D information added.

        Raises:
            OptimizationError: If any molecule fails to embed or optimize.
            Molecule3DPreparationError: If RDKit Mol object is missing for a monomer.
        """
        # Process individual monomer molecules
        monomer_entries = updated_inputs.monomers
        for entry in monomer_entries:
            if entry.rdkit_mol is None:
                raise Molecule3DPreparationError(
                    f"RDKit Mol object is missing for molecule {entry.data_id}. "
                    "Cannot prepare 3D geometry without it."
                )

            try:
                optimized_path, mol_block = self._optimization(
                    molecule_name=entry.data_id,
                    mol=entry.rdkit_mol,
                    cache_dir=self.molecule_3d_path,
                )
                entry.molecule_3Dmol_path = optimized_path
                entry.molecule_3Dmol_file = mol_block
            except Exception as e:
                raise OptimizationError(
                    f"Error optimizing molecule {entry.data_id}: {str(e)}"
                ) from e

        # Process combined reactant and product complexes for each reaction
        for reaction in prepared_reactions:
            if reaction.combined_reactants is not None:
                try:
                    optimized_path, mol_block = self._optimization(
                        molecule_name=f"{reaction.reaction_id}_pre",
                        mol=reaction.combined_reactants,
                        cache_dir=self.full_templates_path,
                        separate_fragments=True,
                    )
                    reaction.reactants_3Dmol_path = optimized_path
                    reaction.reactants_3Dmol_file = mol_block
                except Exception as e:
                    raise OptimizationError(
                        f"Error optimizing reactant complex for reaction {reaction.reaction_id}: {str(e)}"
                    ) from e

            if reaction.combined_products is not None:
                try:
                    optimized_path, mol_block = self._optimization(
                        molecule_name=f"{reaction.reaction_id}_post",
                        mol=reaction.combined_products,
                        cache_dir=self.full_templates_path,
                        separate_fragments=True,
                    )
                    reaction.products_3Dmol_path = optimized_path
                    reaction.products_3Dmol_file = mol_block
                except Exception as e:
                    raise OptimizationError(
                        f"Error optimizing product complex for reaction {reaction.reaction_id}: {str(e)}"
                    ) from e

        return updated_inputs, prepared_reactions

    def _separate_fragments_3d(
        self, mol: Mol,
    ) -> Mol:
        """Separate disconnected fragments by translating one fragment in 3D space.

        Useful for reactant or product complexes that consist of multiple
        disconnected molecules to prevent overlapping coordinates. The separation
        distance is calculated based on the molecular weight of the complex.

        Args:
            mol: RDKit molecule containing one or more fragments.

        Returns:
            The modified molecule with separated fragments.

        Raises:
            FragmentSeparationError: If more than two fragments are detected.
        """
        # Get molecular fragments (connected components)
        frags = Chem.GetMolFrags(mol, asMols=False)
        
        # No separation needed if only one fragment is present
        if len(frags) < 2:
            return mol
        
        # Check if we have exactly two fragments
        if len(frags) > 2:
            raise FragmentSeparationError(
                f"Expected 2 fragments for separation, but found {len(frags)}. "
                "Please check the molecule structure."
            )
        
        # Calculate separation distance based on molecular weight
        mw = Descriptors.MolWt(mol)
        mw_magnitude = round(mw/100)
        shift = ((mw_magnitude+1)*8.0, 0.0, 0.0)
    
        # Get the molecule's conformation
        conf = mol.GetConformer()
        shift_arr = np.array(shift, dtype=float)

        # Translate all atoms in the second fragment
        for atom_idx in frags[1]:
            pos = np.array(conf.GetAtomPosition(atom_idx))
            conf.SetAtomPosition(atom_idx, pos + shift_arr)

        return mol

    def _optimization(
        self,
        molecule_name: str,
        mol: Mol,
        cache_dir: Path,
        separate_fragments: bool = False,
    ) -> Tuple[Path, str]:
        """Embed a molecule in 3D, optionally separate fragments, optimize geometry,
        and save the result as a .mol file.

        The optimization process includes:
        1. Sanitization and property cache update
        2. 3D coordinate embedding using ETKDG method
        3. Optional fragment separation for multi-molecule complexes
        4. Geometry optimization using MMFF force field
        5. File saving and mol block generation

        Args:
            molecule_name: Name used for the output file.
            mol: RDKit molecule to be optimized.
            cache_dir: Directory where the .mol file will be saved.
            separate_fragments: Whether to separate disconnected fragments before optimization.

        Returns:
            Tuple containing:
                - Path to the saved .mol file
                - Mol block string of the optimized structure

        Raises:
            OptimizationError: If atom count changes or optimization fails.
        """
        # Record initial atom count for integrity check
        n_atoms_start = mol.GetNumAtoms(onlyExplicit=True)

        # Update property cache and sanitize molecule
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES,
        )

        # Generate initial 3D coordinates using ETKDG method
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        if result == -1:
            raise OptimizationError(f"Failed to embed molecule {molecule_name} in 3D.")

        # Separate fragments if this is a complex (reactants or products)
        if separate_fragments:
            mol = self._separate_fragments_3d(mol)

        # Perform geometry optimization using MMFF force field
        ff_result = AllChem.MMFFOptimizeMolecule(mol)
        if ff_result == -1:
            raise OptimizationError(f"MMFF optimization failed for {molecule_name}.")
        if ff_result == 1:
            print(f"Warning: MMFF optimization did not converge for {molecule_name}.")

        # Verify atom count integrity (no atoms lost during optimization)
        if n_atoms_start != mol.GetNumAtoms(onlyExplicit=True):
            raise OptimizationError(
                f"Atom count mismatch for {molecule_name}: started with {n_atoms_start} "
                f"explicit atoms but ended with {mol.GetNumAtoms(onlyExplicit=True)}."
            )

        # Ensure cache directory exists
        os.makedirs(cache_dir, exist_ok=True)

        # Save optimized structure to file
        output_path = Path(cache_dir) / f"{molecule_name}.mol"
        print(f"Saving optimized {molecule_name} to {output_path}")

        Chem.MolToMolFile(mol, str(output_path))
        mol_block = Chem.MolToMolBlock(mol)

        return output_path, mol_block
