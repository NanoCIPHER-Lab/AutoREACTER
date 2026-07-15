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

# Third-party imports
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import Descriptors

# Local imports from the AutoREACTER package
from AutoREACTER.input_parser import SimulationSetup
from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import ReactionMetadata
from typing import TYPE_CHECKING

from AutoREACTER.session import Session

if TYPE_CHECKING:
    from AutoREACTER.session import Session


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

    def __init__(self, session: "Session"):
        """Initialize the 3D preparation utility using the shared AutoREACTER session.

        In AutoREACTER, session.staging_dir is the working/cache directory.
        """
        self.session = session
        self.inputs = session.inputs

        self.cache_dir = Path(session.staging_dir) / "3D_molecules"
        os.makedirs(self.cache_dir, exist_ok=True)

        self.molecule_3d_path = self.cache_dir / "molecules_3Dmol"
        os.makedirs(self.molecule_3d_path, exist_ok=True)

        self.full_templates_path = self.cache_dir / "full_templates_3Dmol"
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
        session: Session,
    ) -> None:
        """Prepare 3D geometries for all monomers and reaction complexes.

        This method processes individual monomer molecules and the combined
        reactant/product complexes for each reaction, optimizing their
        3D structures and attaching the resulting file paths and mol blocks
        to the respective data objects.

        Args:
            session: The current Session containing validated inputs and reaction metadata.

        Returns:
            None. This method updates the session's inputs and reaction metadata in place.

        Raises:
            OptimizationError: If any molecule fails to embed or optimize.
            Molecule3DPreparationError: If RDKit Mol object is missing for a monomer.
        """
        # Process individual monomer molecules
        updated_inputs, prepared_reactions = session.inputs, session.reaction_metadata
        monomer_entries = updated_inputs.monomers
        for entry in monomer_entries:
            if not entry.status:
                continue  # Skip monomers that ommitted during selections
            if entry.rdkit_mol is None:
                raise Molecule3DPreparationError(
                    f"RDKit Mol object is missing for molecule {entry.name}. "
                    "Cannot prepare 3D geometry without it."
                )

            try:
                mol = Chem.AddHs(entry.rdkit_mol)
                optimized_path = self._optimization(
                    molecule_name=entry.name,
                    mol=mol,
                    cache_dir=self.molecule_3d_path,
                )
                entry.molecule_3Dmol_path = optimized_path
            except Exception as e:
                raise OptimizationError(
                    f"Error optimizing molecule {entry.name}: {str(e)}"
                ) from e

        # Process combined reactant and product complexes for each reaction
        for reaction in prepared_reactions:
            if not reaction.activity_stats:
                continue  # Skip reactions that were filtered out during duplicate filtering
            if reaction.reactant_combined_RDmol is not None:
                try:
                    optimized_path = self._optimization(
                        molecule_name=f"pre{reaction.reaction_id}",
                        mol=reaction.reactant_combined_RDmol,
                        cache_dir=self.full_templates_path,
                        separate_fragments=True,
                    )
                    reaction.reactant_combined_3Dmol_path = optimized_path
                except Exception as e:
                    raise OptimizationError(
                        f"Error optimizing reactant complex for reaction {reaction.reaction_id}: {str(e)}"
                    ) from e

            if reaction.product_combined_RDmol is not None:
                try:
                    optimized_path = self._optimization(
                        molecule_name=f"post{reaction.reaction_id}",
                        mol=reaction.product_combined_RDmol,
                        cache_dir=self.full_templates_path,
                        separate_fragments=True,
                    )
                    reaction.product_combined_3Dmol_path = optimized_path
                except Exception as e:
                    raise OptimizationError(
                        f"Error optimizing product complex for reaction {reaction.reaction_id}: {str(e)}"
                    ) from e

        return None

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
        shift = ((mw_magnitude+1)*4.0, 0.0, 0.0)
    
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
    ) -> Path:
        """Repair, embed, optimize, and save a molecule without adding hydrogens."""

        # Work on a copy so the ReactionMetadata molecule and its indexing remain
        # unchanged outside this 3D preparation step.
        mol = Chem.Mol(mol)
        n_atoms_start = mol.GetNumAtoms(onlyExplicit=True)

        try:
            mol = self._repair_reaction_molecule_for_3d(mol)
        except Exception as error:
            raise OptimizationError(
                f"Failed to repair molecule {molecule_name} before 3D embedding: "
                f"{error}"
            ) from error

        # Remove any old or partial conformers before embedding.
        mol.RemoveAllConformers()

        params = AllChem.ETKDGv3()
        params.randomSeed = 0xF00D

        embed_result = AllChem.EmbedMolecule(mol, params)

        if embed_result == -1:
            raise OptimizationError(
                f"Failed to embed molecule {molecule_name} in 3D."
            )

        if separate_fragments:
            mol = self._separate_fragments_3d(mol)

        # MMFF generally works for these carbon radicals, but keep a UFF fallback
        # for structures for which MMFF lacks parameters.
        if AllChem.MMFFHasAllMoleculeParams(mol):
            ff_result = AllChem.MMFFOptimizeMolecule(
                mol,
                maxIters=1000,
            )
            force_field_name = "MMFF"
        elif AllChem.UFFHasAllMoleculeParams(mol):
            ff_result = AllChem.UFFOptimizeMolecule(
                mol,
                maxIters=1000,
            )
            force_field_name = "UFF"
        else:
            ff_result = None
            force_field_name = None
            print(
                f"Warning: no MMFF or UFF parameters are available for "
                f"{molecule_name}. Saving the embedded geometry without "
                "force-field optimization."
            )

        if ff_result == -1:
            raise OptimizationError(
                f"{force_field_name} optimization failed for {molecule_name}."
            )

        if ff_result == 1:
            print(
                f"Warning: {force_field_name} optimization did not converge "
                f"for {molecule_name}."
            )

        n_atoms_end = mol.GetNumAtoms(onlyExplicit=True)

        if n_atoms_start != n_atoms_end:
            raise OptimizationError(
                f"Atom count mismatch for {molecule_name}: started with "
                f"{n_atoms_start} explicit atoms but ended with {n_atoms_end}."
            )

        os.makedirs(cache_dir, exist_ok=True)

        output_path = Path(cache_dir) / f"{molecule_name}.mol"
        print(f"Saving optimized {molecule_name} to {output_path}")

        Chem.MolToMolFile(
            mol,
            str(output_path),
            includeStereo=True,
            kekulize=False,
        )

        return output_path

    def _repair_reaction_molecule_for_3d(self, mol: Mol) -> Mol:
        """Repair RDKit reaction-product valence bookkeeping before 3D embedding.

        RunReactants may produce atoms that have both:
        1. explicit hydrogen atoms as neighbors, and
        2. a nonzero NumExplicitHs property inherited from the product SMARTS.

        That double-counts hydrogens and can produce apparent carbon valences of
        five or six. This method changes atom properties only; it does not add or
        remove atoms.

        Neutral, non-aromatic carbon atoms with bond valence three are retained as
        carbon-centered radicals instead of being given an implicit hydrogen.
        """
        repaired = Chem.RWMol(Chem.Mol(mol))
        repaired.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(repaired)

        for atom in repaired.GetAtoms():
            explicit_h_neighbors = sum(
                neighbor.GetAtomicNum() == 1
                for neighbor in atom.GetNeighbors()
            )

            # The hydrogen atoms already exist as real graph atoms. Clear only the
            # duplicate SMARTS hydrogen-count property.
            if explicit_h_neighbors > 0:
                if atom.GetNumExplicitHs() > 0:
                    atom.SetNumExplicitHs(0)

                # Do not let RDKit add another implicit hydrogen during sanitization.
                atom.SetNoImplicit(True)

            if atom.GetAtomicNum() != 6:
                continue

            bond_valence = sum(
                bond.GetValenceContrib(atom)
                for bond in atom.GetBonds()
            )

            # Preserve neutral carbon-centered radical chain ends.
            if (
                not atom.GetIsAromatic()
                and atom.GetFormalCharge() == 0
                and abs(bond_valence - 3.0) < 1.0e-6
            ):
                atom.SetNoImplicit(True)
                atom.SetNumRadicalElectrons(1)

            # Radical-radical coupling gives each carbon its fourth valence.
            elif bond_valence >= 4.0:
                atom.SetNumRadicalElectrons(0)

        repaired_mol = repaired.GetMol()
        repaired_mol.ClearComputedProps()
        repaired_mol.UpdatePropertyCache(strict=False)

        # Full sanitization is now safe because the duplicate hydrogen counts
        # and radical valences have been repaired.
        Chem.SanitizeMol(repaired_mol)

        return repaired_mol