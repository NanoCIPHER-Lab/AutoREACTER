from rdkit import Chem
import pandas as pd
from typing import List
from pathlib import Path
import os

from AutoREACTER.detectors.reaction_detector import ReactionInstance

def prepare_paths( cache, subdir):
    """
    Create and prepare directory structure for CSV output files.
        
    Args:
        cache (str or Path): Base cache directory path
            
    Returns:
        Path: Path object pointing to the CSV cache directory
            
    Example:
        >>> csv_dir = prepare_paths("/path/to/cache")
        >>> print(csv_dir)  # /path/to/cache/csv
    """
    csv_cache = Path(cache) / subdir
    os.makedirs(csv_cache, exist_ok=True)
    return csv_cache

def add_dict_as_new_columns(df_existing, data_dict, titles=("template_reactant_idx", "template_product_idx")):
    df_existing[titles[0]] = pd.Series(list(data_dict.keys())).astype("Int64")
    df_existing[titles[1]] = pd.Series(list(data_dict.values())).astype("Int64")
    return df_existing


def add_column_safe(df, list_data, column_name):
    df[column_name] = pd.Series(list_data).astype("Int64")
    return df


def extract_unique_references(detected_reactions: List[ReactionInstance]) -> list[str]:
    seen = set()
    refs = []

    for metadata in detected_reactions:
        # Prefer the correct 'references' attribute on ReactionInstance, with
        # a backward-compatible fallback to 'reference' if present.
        ref = getattr(metadata, "references", None)
        if ref is None:
            ref = getattr(metadata, "reference", {}) or {}
        else:
            if ref is None:
                ref = {}

        for v in ref.values():
            if isinstance(v, str):
                if v not in seen:
                    seen.add(v)
                    refs.append(v)

            elif isinstance(v, (list, tuple)):
                for u in v:
                    if isinstance(u, str) and u not in seen:
                        seen.add(u)
                        refs.append(u)

            elif isinstance(v, dict):
                for u in v.values():
                    if isinstance(u, str) and u not in seen:
                        seen.add(u)
                        refs.append(u)

    return refs

def compare_set(reaction_metadata_list, _react2, _prod2):
    """
    Compare reactant + product pair against existing reactions.
    Ignores atom mapping.

    Returns False if duplicate pair found.
    """

    # normalize new pair
    react2 = Chem.Mol(_react2)
    prod2 = Chem.Mol(_prod2)

    for atom in react2.GetAtoms():
        atom.SetAtomMapNum(0)
    for atom in prod2.GetAtoms():
        atom.SetAtomMapNum(0)

    smi_r2 = Chem.MolToSmiles(react2, canonical=True)
    smi_p2 = Chem.MolToSmiles(prod2, canonical=True)

    for metadata in reaction_metadata_list:
        react1 = Chem.Mol(metadata.reactant_combined_RDmol)
        prod1 = Chem.Mol(metadata.product_combined_RDmol)

        for atom in react1.GetAtoms():
            atom.SetAtomMapNum(0)
        for atom in prod1.GetAtoms():
            atom.SetAtomMapNum(0)

        smi_r1 = Chem.MolToSmiles(react1, canonical=True)
        smi_p1 = Chem.MolToSmiles(prod1, canonical=True)

        if smi_r1 == smi_r2 and smi_p1 == smi_p2:
            return False

    return True


def compare_rdkit_molecules_canonical(data_smiles_list, mol_smi_2):
    """
    Compares two RDKit molecule SMILES strings to determine if they
    represent the same chemical structure using canonical SMILES.

    Args:
        data_smiles_list (list): List of SMILES strings to compare against.
        mol_smi_2 (str): SMILES string of the second molecule.

    Returns:
        bool: True if molecules are chemically identical, False otherwise.
    """
    if not mol_smi_2:
        return data_smiles_list, False
    mol2 = Chem.MolFromSmiles(mol_smi_2)
    if mol2 is None:
        return data_smiles_list, False
        
    for mol_smi_1 in data_smiles_list:
        mol1 = Chem.MolFromSmiles(mol_smi_1)
        
        # Handle cases where SMILES might be invalid
        if mol1 is None or mol2 is None:
            return data_smiles_list, False  # or raise a ValueError
        # Generate canonical SMILES and compare them
        canonical_smi_1 = Chem.MolToSmiles(mol1, canonical=True)
        canonical_smi_2 = Chem.MolToSmiles(mol2, canonical=True)

        if canonical_smi_1 == canonical_smi_2:
            return data_smiles_list, True
    
    data_smiles_list.append(mol_smi_2)
    return data_smiles_list, False



def prep_for_3d_molecule_generation(data_smiles_list, molecule_dict_csv_path_dict):
    """
    Prepares a dictionary of RDKit molecule objects for 3D generation.
    Combines molecules from a list of SMILES strings and a dictionary of CSV paths.
    Args:
        data_smiles_list (list): List of SMILES strings to convert to RDKit molecules.
        molecule_dict_csv_path_dict (dict): Dictionary containing reactant and product molecules from CSV paths.
    Returns:
        dict: A dictionary where keys are identifiers (e.g., 'data_1', 'pre_1', 'post_1') and values are RDKit molecule objects.
    """
    molecule_dict = {}
    for i, smiles in enumerate(data_smiles_list):
        key = f"data_{i+1}"
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        molecule = Chem.AddHs(molecule)
        molecule_dict[key] = molecule
    for i, (key, value) in enumerate(molecule_dict_csv_path_dict.items()):
        reactant = value.get("reactant")
        product = value.get("product")
        if reactant and product:
            molecule_dict[f"pre_{i+1}"] = reactant
            molecule_dict[f"post_{i+1}"] = product
    return molecule_dict
