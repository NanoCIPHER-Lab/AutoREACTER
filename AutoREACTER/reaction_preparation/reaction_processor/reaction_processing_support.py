from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops
from rdkit.Chem.rdChemReactions import ChemicalReaction

from AutoREACTER.reaction_preparation.reaction_processor.utils import compare_set

if TYPE_CHECKING:
    from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import ReactionMetadata


class MappingError(Exception):
    """Raised when atom mapping between reactants and products fails or is inconsistent."""


class SMARTSParsingError(Exception):
    """Raised when parsing SMILES or reaction SMARTS fails."""


def detect_duplicates(reaction_metadata_list: list["ReactionMetadata"]) -> list["ReactionMetadata"]:
    """
    Filters duplicate reactions based on reactant and product molecules.
    
    Args:
        reaction_metadata_list: List of reaction metadata to filter
        
    Returns:
        List of unique reactions; duplicates marked with activity_stats=False
    """
    unique_metadata: list["ReactionMetadata"] = []
    
    for reaction in reaction_metadata_list:
        # Compare current reaction's reactants and products against unique reactions collected so far
        reactants = reaction.reactant_combined_RDmol
        products = reaction.product_combined_RDmol

        # Keep reaction if it's unique, otherwise mark as duplicate
        if compare_set(unique_metadata, reactants, products):
            unique_metadata.append(reaction)
        else:
            reaction.activity_stats = False
            
    return unique_metadata


# --- CORE REACTION LOGIC ---

def assign_first_shell_and_initiators(  reactant_combined: Chem.Mol, 
                                        product_combined: Chem.Mol, 
                                        reversed_mapping_dict: dict[int, int]) -> tuple[list[int], list[int]]:
    """
    Identifies atoms in the first coordination shell (atoms with map numbers < 999) and initiator atoms.
    Initiators are atoms with map numbers 1 or 2 (typically the reactive centers).
    
    Args:
        reactant_combined: Combined reactant molecule
        product_combined: Combined product molecule
        reversed_mapping_dict: Product idx -> Reactant idx mapping
        
    Returns:
        Tuple of (first_shell atom indices, initiator atom indices)
        
    Raises:
        ValueError: If exactly 2 initiators are not found
    """
    first_shell = []
    initiator_idxs = []

    for p_atom in product_combined.GetAtoms():
        # Only process atoms with valid map numbers (< 999 indicates non-byproduct)
        if p_atom.GetAtomMapNum() < 999:
            p_idx = p_atom.GetIdx()

            if p_idx not in reversed_mapping_dict:
                raise ValueError(f"Mapping error: product atom {p_idx} not found in mapping_dict")

            r_idx = reversed_mapping_dict[p_idx]
            atom = reactant_combined.GetAtomWithIdx(r_idx)
            atom.SetAtomMapNum(p_atom.GetAtomMapNum())

            first_shell.append(r_idx)

            # Initiators are atoms with map numbers 1 or 2
            if p_atom.GetAtomMapNum() in [1, 2]:
                initiator_idxs.append(r_idx)
    
    if len(initiator_idxs) != 2:
        raise ValueError(f"Expected 2 initiators, got {len(initiator_idxs)}: {initiator_idxs}")

    return first_shell, initiator_idxs

def detect_byproducts(  product_combined: Chem.Mol, 
                        reversed_mapping_dict: dict[int, int], 
                        delete_atoms: bool) -> list[int]:
    """
    Identifies byproduct atoms (smallest molecular fragment) and maps them back to reactant space.
    
    Args:
        product_combined: Combined product molecule
        reversed_mapping_dict: Product idx -> Reactant idx mapping
        delete_atoms: Whether to perform byproduct detection
        
    Returns:
        List of reactant indices corresponding to byproduct atoms
    """
    if not delete_atoms:
        return []

    # Get tuples of original atom indices for each fragment
    frags_indices = rdmolops.GetMolFrags(product_combined)
    
    # Find the tuple with the smallest number of atoms
    smallest_frag_indices = min(frags_indices, key=len)

    byproduct_reactant_indices = []

    # Map byproduct product indices back to reactant indices
    for p_idx in smallest_frag_indices:
        if p_idx in reversed_mapping_dict:
            byproduct_reactant_indices.append(reversed_mapping_dict[p_idx])

    return byproduct_reactant_indices

def validate_mapping(df: pd.DataFrame, reactant: Chem.Mol, product: Chem.Mol) -> None:
    """
    Validates atom mapping consistency: checks for required columns, duplicates, bounds, and completeness.
    
    Args:
        df: Dataframe containing reactant_idx and product_idx columns
        reactant: Reactant molecule
        product: Product molecule
        
    Raises:
        MappingError: If any validation check fails
    """
    # Ensure dataframe exists and has required columns
    if df is None or df.empty:
        raise MappingError("Mapping validation failed: empty dataframe")

    # Check for required columns
    required_cols = {"reactant_idx", "product_idx"}
    if not required_cols.issubset(df.columns):
        raise MappingError(f"Mapping validation error: required columns {required_cols} not found in dataframe.")

    # Extract indices and perform validation checks
    r_idxs = df["reactant_idx"].dropna().tolist()
    p_idxs = df["product_idx"].dropna().tolist()

    # Atom counts must match
    if len(r_idxs) != len(p_idxs):
        raise MappingError(f"Mapping validation error: mismatch in atom counts between reactant and product.")

    # No duplicate mappings (1-to-1 mapping required)
    if len(set(r_idxs)) != len(r_idxs):
        raise MappingError(f"Mapping validation error: duplicate indices found in reactant mapping.")
    if len(set(p_idxs)) != len(p_idxs):
        raise MappingError(f"Mapping validation error: duplicate indices found in product mapping.")

    # Indices must be within molecule bounds
    if any(idx >= reactant.GetNumAtoms() for idx in r_idxs):
        raise MappingError(f"Mapping validation error: reactant index out of bounds.")
    if any(idx >= product.GetNumAtoms() for idx in p_idxs):
        raise MappingError(f"Mapping validation error: product index out of bounds.")

    # All atoms must be mapped (complete mapping)
    if len(r_idxs) != reactant.GetNumAtoms():
        raise MappingError(f"Mapping validation error: incomplete mapping for reactant.")
    if len(p_idxs) != product.GetNumAtoms():
        raise MappingError(f"Mapping validation error: incomplete mapping for product.")

# --- ATOM MAPPING ---

def assign_atom_map_numbers_and_set_isotopes(r1: Chem.Mol, r2: Chem.Mol) -> None:
    """
    Assigns unique map numbers and isotopes to reactant atoms for tracking through reaction.
    Isotopes survive RDKit's reaction engine, allowing atom identity recovery post-reaction.
    
    Args:
        r1: First reactant molecule
        r2: Second reactant molecule
    """
    # Assign map numbers 1001+ to first reactant atoms
    for atom in r1.GetAtoms():
        idx = 1001 + atom.GetIdx()
        atom.SetAtomMapNum(idx)
        atom.SetIsotope(idx)  # Isotope survives the reaction

    # Assign map numbers 2001+ to second reactant atoms
    for atom in r2.GetAtoms():
        idx = 2001 + atom.GetIdx()
        atom.SetAtomMapNum(idx)
        atom.SetIsotope(idx) # Isotope survives the reaction
        
def reassign_atom_map_numbers_by_isotope(mol: Chem.Mol) -> None:
    """
    Restores atom map numbers from isotope values after reaction.
    RDKit's reaction engine preserves isotopes, allowing recovery of original atom identities.
    
    Args:
        mol: Product molecule with isotope information
    """
    for atom in mol.GetAtoms():
        surviving_idx = atom.GetIsotope()
        if surviving_idx != 0:
            atom.SetAtomMapNum(surviving_idx)  # Restore original map number
            atom.SetIsotope(0)                 # Clear isotope to restore normal chemistry

def build_atom_index_mapping(   reactant_combined: Chem.Mol, 
                                product_combined: Chem.Mol) -> tuple[dict[int, int], pd.DataFrame]:
    """
    Builds bidirectional atom index mapping between reactants and products using map numbers.
    
    Args:
        reactant_combined: Combined reactant molecule
        product_combined: Combined product molecule
        
    Returns:
        Tuple of (mapping dict: reactant_idx -> product_idx, dataframe with mapping)
    """
    mapping_dict = {}
    
    # Pre-index product atoms by map number for O(1) lookup
    product_map = {
        atom.GetAtomMapNum(): atom.GetIdx()
        for atom in product_combined.GetAtoms()
        if atom.GetAtomMapNum() != 0
    }
    
    rows = []
    for r_atom in reactant_combined.GetAtoms():
        r_map_num = r_atom.GetAtomMapNum()

        # Match reactant atom to product atom via map number
        if r_map_num != 0 and r_map_num in product_map:
            r_idx = r_atom.GetIdx()
            p_idx = product_map[r_map_num]

            mapping_dict[r_idx] = p_idx
            rows.append({
                "reactant_idx": r_idx,
                "product_idx": p_idx
            })

    df = pd.DataFrame(rows)
    return mapping_dict, df

def reveal_template_map_numbers(mol: Chem.Mol) -> None:
    """
    Restores map numbers from RDKit's internal 'old_mapno' property for visualization.
    RDKit stores original map numbers in this property after reaction execution.
    
    Args:
        mol: Product molecule
    """
    for atom in mol.GetAtoms():
        if atom.HasProp('old_mapno'):
            map_num = atom.GetIntProp('old_mapno')
            atom.SetAtomMapNum(map_num)

def clear_isotopes(mol_1: Chem.Mol, mol_2: Chem.Mol) -> None:
    """
    Clears isotope values from molecules to restore normal chemistry after using isotopes for atom tracking.
    
    Args:
        mol_1: First molecule to clear
        mol_2: Second molecule to clear
    """
    for atom in mol_1.GetAtoms():
        atom.SetIsotope(0)
    for atom in mol_2.GetAtoms():
        atom.SetIsotope(0)

# --- BUILDERS ---

def build_reaction(rxn_smarts: str) -> ChemicalReaction:
    """Builds RDKit ChemicalReaction object from SMARTS string."""
    return AllChem.ReactionFromSmarts(rxn_smarts)

def build_reactants(reactant_smiles_1: str, reactant_smiles_2: str) -> tuple[Chem.Mol, Chem.Mol]:
    """
    Builds reactant molecules from SMILES strings with explicit hydrogens added.
    
    Args:
        reactant_smiles_1: SMILES string for first reactant
        reactant_smiles_2: SMILES string for second reactant
        
    Returns:
        Tuple of (reactant1 molecule, reactant2 molecule) with explicit hydrogens
    """
    mol_reactant_1 = Chem.MolFromSmiles(reactant_smiles_1)
    if mol_reactant_1 is None:
        raise SMARTSParsingError(f"Failed to parse first reactant SMILES: {reactant_smiles_1!r}")
    mol_reactant_1 = Chem.AddHs(mol_reactant_1)
    
    mol_reactant_2 = Chem.MolFromSmiles(reactant_smiles_2)
    if mol_reactant_2 is None:
        raise SMARTSParsingError(f"Failed to parse second reactant SMILES: {reactant_smiles_2!r}")
    mol_reactant_2 = Chem.AddHs(mol_reactant_2)
    
    return mol_reactant_1, mol_reactant_2

def build_reaction_tuple(same_reactants: bool, mol_reactant_1: Chem.Mol, mol_reactant_2: Chem.Mol) -> list:
    """
    Builds list of reactant pairs to process. If reactants are identical, returns single pair.
    Otherwise returns both orderings to account for reaction directionality.
    
    Args:
        same_reactants: Whether both reactants are identical
        mol_reactant_1: First reactant molecule
        mol_reactant_2: Second reactant molecule
        
    Returns:
        List of reactant pairs [[r1, r2], ...] to process
    """
    if same_reactants:
        return [[mol_reactant_1, mol_reactant_1]]

    return [[mol_reactant_1, mol_reactant_2], [mol_reactant_2, mol_reactant_1]]

