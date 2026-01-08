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
from typing import Dict, Tuple
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol

# TODO: Move this to a configuration module (cache.py) or an environment variable.
# This is the directory where optimized .mol files will be written.
# The cache directory can be configured via the LUNAR_CACHE_DIR environment variable.
# If not set, it falls back to a repository-relative "cache/lunar" directory.
CACHE_ENV_VAR = "LUNAR_CACHE_DIR"
cache = os.environ.get(
    CACHE_ENV_VAR,
    os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
        "cache",
        "lunar",
    ),
)


def separate_fragments_3d(mol: Mol, shift: Tuple[float, float, float] = (8.0, 0.0, 0.0)) -> Mol:
    """
    If `mol` contains multiple disconnected fragments, translate the second fragment by `shift` in 3D.

    This is useful when you have a combined molecule (e.g., two reactants placed in the same RDKit Mol
    object) that are topologically disconnected and you want to ensure their 3D coordinates don't overlap.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The molecule that may contain multiple disconnected fragments. Must have a conformer
        (3D coordinates) for translational changes to apply.
    shift : tuple of float, optional
        A 3-tuple specifying the (x, y, z) translation to apply to the second fragment.
        Default is (8.0, 0.0, 0.0) which shifts the fragment +8 A along x.

    Returns
    -------
    rdkit.Chem.rdchem.Mol
        The same RDKit molecule object with the second fragment translated in-place.

    Raises
    ------
    ValueError
        If more than two disconnected fragments are found. The function is currently designed
        to handle at most two fragments.
    """
    # Get tuple of atom index lists for each fragment (asMols=False returns tuples of atom indices)
    frags = Chem.GetMolFrags(mol, asMols=False)

    # If there's only one fragment, nothing to do
    if len(frags) < 2:
        return mol

    # For now we only support exactly two fragments. Raise if more are present to avoid silent errors.
    if len(frags) > 2:
        raise ValueError("More than 2 fragments found")

    # Get the conformer (3D coordinates) for the molecule
    conf = mol.GetConformer()
    shift_arr = np.array(shift, dtype=float)

    # Translate atoms in the second fragment (fragment index 1)
    for atom_idx in frags[1]:
        pos = np.array(conf.GetAtomPosition(atom_idx))
        conf.SetAtomPosition(atom_idx, pos + shift_arr)

    return mol


def optimization(molecule_name: str, mol: Mol, cache_dir: str) -> str:
    """
    Embed and optimize a molecule's 3D geometry and save the result to a .mol file in `cache_dir`.

    The function:
    - Records the initial explicit atom count.
    - Updates RDKit internal properties and sanitizes the molecule (except 'properties' sanitization).
    - Embeds the molecule into 3D using RDKit's ETKDG algorithm.
    - If the molecule name begins with "pre" or "post" it attempts to separate fragments in 3D
      (useful for reactant/product complexes where there are disconnected fragments).
    - Runs an MMFF optimization to locally relax the geometry.
    - Verifies the explicit atom count hasn't changed (sanity check).
    - Writes the resulting structure to `{cache_dir}/{molecule_name}.mol`.

    Parameters
    ----------
    molecule_name : str
        A name for the molecule that will be used to create the output filename.
    mol : rdkit.Chem.rdchem.Mol
        The molecule to embed/optimize. It may or may not already have H atoms added.
    cache_dir : str
        Directory to save the optimized .mol file to.

    Returns
    -------
    str
        Path to the saved .mol file.

    Raises
    ------
    ValueError
        If the explicit atom count changes during processing, indicating something went wrong.
    """
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
    if molecule_name.startswith("pre") or molecule_name.startswith("post"):
        mol = separate_fragments_3d(mol, shift=(8.0, 0.0, 0.0))

    # Optimize geometry using the MMFF force field (works well for many organic molecules).
    # This modifies the conformation in-place.
    AllChem.MMFFOptimizeMolecule(mol)

    # Integrity check: make sure explicit atom count didn't change during operations.
    # A mismatch usually indicates an unintended removal or addition of atoms.
    if n_atoms_start != mol.GetNumAtoms(onlyExplicit=True):
        raise ValueError("Explicit atom count changed!")

    # Create output directory if it doesn't exist
    os.makedirs(cache_dir, exist_ok=True)

    # Save optimized geometry as a .mol file
    output_path = os.path.join(cache_dir, f"{molecule_name}.mol")
    print(f"Saving optimized {molecule_name} to {output_path}")
    Chem.MolToMolFile(mol, output_path)

    return output_path


def prepare_3d_molecule(cache_dir: str, molecule_dict: Dict[str, Mol]) -> Dict[str, str]:
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
    for mol_name, mol in molecule_dict.items():
        output_path = optimization(mol_name, mol, cache_dir)
        molecule_dict[mol_name] = output_path
    return molecule_dict


if __name__ == "__main__":
    # Example usage: construct reactants, run a SMARTS-based reaction, prepare 3D structures,
    # and write optimized .mol files to the cache directory.

    # Use the local `cache` variable (TODO: replace with config/cache module)
    cache = r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\lunar"

    # Example SMILES for two reactants 
    reactant_smiles1 = "C1=CC=C(C(=C1)C(=O)O)O"
    reactant_smiles2 = "OCCC(O)=O"

    # SMARTS describing esterification example
    reaction_smarts = (
        "[O;!$(OC=*):1]-[H:3].[CX3:2](=[O:5])[OX2H1:4]>>[OX2:1]-[CX3:2](=[O:5]).[O:4]-[H:3]"
    )

    # Load reactants from SMILES and add explicit hydrogens (recommended for 3D embedding)
    reactant1 = Chem.MolFromSmiles(reactant_smiles1)
    reactant1 = Chem.AddHs(reactant1)
    reactant2 = Chem.MolFromSmiles(reactant_smiles2)
    reactant2 = Chem.AddHs(reactant2)

    # Create an RDKit reaction from SMARTS
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)

    # Combine reactants into a single RDKit Mol for visualization/complex handling
    combined_reactants = Chem.CombineMols(reactant1, reactant2)

    # Run the reaction to produce product sets. RunReactants returns a tuple of tuples,
    # where each outer tuple is one possible product set. We'll take the first generated set.
    products = rxn.RunReactants((reactant1, reactant2))
    if not products:
        raise RuntimeError("Reaction produced no products for the provided reactants and SMARTS.")
    products = products[0]  # Take the first product set

    # Combine produced product fragments into a single RDKit Mol
    combined_products = Chem.CombineMols(*products)

    # Create a mapping of names -> RDKit Mol objects to be prepared (embedded/optimized)
    molecule_dict = {
        "data1": reactant1,
        "data2": reactant2,
        "pre1": combined_reactants,   # combined reactant complex
        "post1": combined_products,   # combined product complex
    }

    print("Cache directory:", cache)
    prepared_molecules = prepare_3d_molecule(cache_dir=cache, molecule_dict=molecule_dict)
    print("Prepared molecules and saved paths:", prepared_molecules)
