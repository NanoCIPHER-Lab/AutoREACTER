"""
Reaction Mapping and Analysis Module

This module provides tools for processing chemical reactions, mapping atoms between
reactants and products, and generating comprehensive analysis data for polymer chemistry.
It uses RDKit for molecular operations and pandas for data management.

Key features:
- Atom mapping between reactants and products
- Reaction validation and consistency checking
- CSV export of mapping data
- Visualization of reaction grids
- Support for polymerization reactions

Author: Janitha Mahanthe
Date: 1/21/2026
Version: 0.1.0
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops
import pandas as pd
from pathlib import Path
import os
import itertools
try:
    # Case 1: correct when reaction_template_pipeline is a proper package
    from .util import compare_products

except (ImportError, ModuleNotFoundError):
    try:
        # Case 2: fully-qualified import under your installed namespace
        from reaction_lammps_mupt.reaction_template_builder.reaction_template_pipeline.util import (
            compare_products,
        )

    except (ImportError, ModuleNotFoundError):
        try:
            # Case 3: if reaction_template_pipeline is installed as a top-level package
            from reaction_template_pipeline.util import compare_products

        except (ImportError, ModuleNotFoundError):
            # Case 4: running as a loose script from the same directory
            from util import compare_products


def is_number_in_set(set_of_tuples, reactant):
    """
    Find a tuple in the set that contains atom indices matching specific atom map numbers.
    
    This function searches for atoms in the reactant molecule that have atom map numbers
    101 or 102, then checks if those atom indices exist in any tuple within the provided set.
    
    Args:
        set_of_tuples (set of tuples): Set of atom index tuples to search through
        reactant (Chem.Mol): RDKit molecule object to search for mapped atoms
        
    Returns:
        tuple: The first tuple found containing a matching atom index
        
    Raises:
        ValueError: If no matching atom is found in the provided set of tuples
        
    Example:
        >>> matches = {(0, 1, 2), (3, 4, 5)}
        >>> result = is_number_in_set(matches, reactant_mol)
    """
    for atom in reactant.GetAtoms():
        # Check for specific atom map numbers (101 or 102)
        if atom.GetAtomMapNum() == 101 or atom.GetAtomMapNum() == 102:
            x = atom.GetIdx()
            for t in set_of_tuples:
                if x in t:
                    print(f"Found matching atom index {x} in tuple {t}")
                    return t
    raise ValueError("No matching atom found in the provided set of tuples.")


def smart_mapping(reactant, smarts_template, match_tuple):
    """
    Apply atom mapping from SMARTS template to reactant molecule based on match indices.
    
    This function transfers atom map numbers from a SMARTS template to corresponding
    atoms in the reactant molecule using the provided match tuple that maps SMARTS
    atom indices to reactant atom indices.
    
    Args:
        reactant (Chem.Mol): The reactant molecule to apply mapping to
        smarts_template (Chem.Mol): SMARTS pattern molecule with atom map numbers
        match_tuple (tuple): Tuple mapping SMARTS atom indices to reactant atom indices
        
    Returns:
        None: Modifies the reactant molecule in-place by setting atom map numbers
        
    Note:
        Only atoms with non-zero map numbers in the SMARTS template are processed.
        The function does nothing if match_tuple is empty.
    """
    if not match_tuple:
        return
    
    # Create mapping from SMARTS atom indices to their map numbers
    SMARTS_index_to_Map_Number = {}
    for smarts_atom in smarts_template.GetAtoms():
        map_num = smarts_atom.GetAtomMapNum()
        if map_num != 0:  # Only process atoms with explicit mapping
            SMARTS_index_to_Map_Number[smarts_atom.GetIdx()] = map_num
    
    # Apply mapping to corresponding atoms in reactant
    for smarts_pos, map_num in SMARTS_index_to_Map_Number.items():
        if smarts_pos < len(match_tuple):
            atom_index_in_mol = match_tuple[smarts_pos]
            atom = reactant.GetAtomWithIdx(atom_index_in_mol)
            atom.SetAtomMapNum(map_num)


def is_consecutive(num_list):
    """
    Check if a list of numbers forms a consecutive sequence without gaps.
    
    Args:
        num_list (list): List of integers to check
        
    Returns:
        bool: True if the numbers are consecutive, False otherwise
        
    Examples:
        >>> is_consecutive([1, 2, 3, 4])  # Returns True
        >>> is_consecutive([1, 3, 4])     # Returns False
        >>> is_consecutive([])            # Returns False
    """
    if not num_list:
        return False
    
    # Check for consecutive sequence: unique numbers and max-min+1 equals length
    return (len(set(num_list)) == len(num_list) and 
            max(num_list) - min(num_list) + 1 == len(num_list))


def prepare_paths(cache):
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
    csv_cache = Path(cache) / "csv" / "reactant_product_mapping"
    os.makedirs(csv_cache, exist_ok=True)
    return csv_cache


def build_reaction(rxn_smarts):
    """
    Create an RDKit reaction object from SMARTS string.
    
    Args:
        rxn_smarts (str): Reaction SMARTS string
        
    Returns:
        Chem.rdChemReactions.ChemicalReaction: RDKit reaction object
        
    Example:
        >>> rxn = build_reaction("[C:1]=[O:2]>>[C:1]-[O:2]")
    """
    return AllChem.ReactionFromSmarts(rxn_smarts)


def build_reactants(reactant_smiles_1, reactant_smiles_2):
    """
    Create RDKit molecule objects from SMILES strings with explicit hydrogens.
    
    Args:
        reactant_smiles_1 (str): SMILES string for first reactant
        reactant_smiles_2 (str): SMILES string for second reactant
        
    Returns:
        tuple: (mol_reactant_1, mol_reactant_2) - Both molecules with explicit hydrogens
        
    Example:
        >>> mol1, mol2 = build_reactants("CCO", "CC=O")
    """
    mol_reactant_1 = Chem.MolFromSmiles(reactant_smiles_1)
    mol_reactant_1 = Chem.AddHs(mol_reactant_1)  # Add explicit hydrogens
    mol_reactant_2 = Chem.MolFromSmiles(reactant_smiles_2)
    mol_reactant_2 = Chem.AddHs(mol_reactant_2)  # Add explicit hydrogens
    return mol_reactant_1, mol_reactant_2


def process_reactions(rxn, csv_cache, reaction_tuple, key=None, 
                      molecule_and_csv_path_dict=None, delete_atoms=False):
    """
    Process chemical reactions, map atoms between reactants and products, and save data.
    
    This is the main processing function that:
    1. Runs reactions on reactant pairs
    2. Maps atoms between reactants and products
    3. Validates mapping consistency
    4. Saves mapping data to CSV files
    5. Identifies reaction features (first shell, initiators, byproducts)
    
    Args:
        rxn (Chem.rdChemReactions.ChemicalReaction): RDKit reaction object
        csv_cache (Path): Directory path for CSV output files
        reaction_tuple (list): List of reactant pairs to process
        key (str or int, optional): Identifier for the reaction set
        molecule_and_csv_path_dict (dict, optional): Dictionary to store results
        delete_atoms (bool, optional): Whether to identify byproduct atoms for deletion
        
    Returns:
        tuple: (mols, output, molecule_and_csv_path_dict)
            - mols: List of reactant and product molecules for visualization
            - output: String output (currently empty, reserved for future use)
            - molecule_and_csv_path_dict: Dictionary containing all processed data
            
    Raises:
        ValueError: If mapping validation fails at any step
        
    Note:
        This function performs extensive validation including:
        - Atom count consistency between reactants and products
        - Consecutive indexing requirements
        - Initiator count validation (must be exactly 2)
    """
    if molecule_and_csv_path_dict is None:
        molecule_and_csv_path_dict = {}
    
    total_products = 0
    mols = []

    # Process each reactant pair in the reaction tuple
    for j, pair in enumerate(reaction_tuple):
        r1, r2 = Chem.Mol(pair[0]), Chem.Mol(pair[1])
        products = rxn.RunReactants((r1, r2))
        
        # Get substructure matches for each reactant against reaction templates
        matches_1 = r1.GetSubstructMatches(rxn.GetReactants()[0])
        matches_2 = r2.GetSubstructMatches(rxn.GetReactants()[1])

        # Process each product set generated by the reaction
        for product_set in products:
            # Clean up any existing DataFrame
            if "df" in locals():
                del df
            
            
            # Initialize data structures
            df = pd.DataFrame(columns=["reactant_idx", "product_idx"])
            first_shell = []      # Atoms in first coordination shell
            initiator_idxs = []   # Initiator atom indices
            mapping_dict = {}     # Additional mapping dictionary

            # Apply atom mapping from reaction properties
            num_total_atoms = 0
            for i, product in enumerate(product_set):
                if i > 0:
                    num_total_atoms += product_set[i - 1].GetNumAtoms()
                for atom in product.GetAtoms():
                    if atom.HasProp("react_idx"):
                        r_idx = atom.GetIntProp("react_idx")
                        a_idx = atom.GetIntProp("react_atom_idx")
                        r = (r1, r2)[r_idx]
                        r_atom = r.GetAtomWithIdx(a_idx)
                        map_num = atom.GetIdx() + 1 + num_total_atoms + 100
                        r_atom.SetAtomMapNum(map_num)
                        atom.SetAtomMapNum(map_num)

            # Combine reactants and products for easier processing
            reactant_combined = Chem.CombineMols(r1, r2)
            product_combined = Chem.CombineMols(*product_set)

            # Create mapping between reactant and product atoms based on map numbers
            for r_atom in reactant_combined.GetAtoms():
                for p_atom in product_combined.GetAtoms():
                    if r_atom.GetAtomMapNum() == p_atom.GetAtomMapNum():
                        new_row = pd.DataFrame([{"reactant_idx": r_atom.GetIdx(), 
                                                "product_idx": p_atom.GetIdx()}])
                        df = pd.concat([df, new_row], ignore_index=True)
                        break

            # Validate mapping completeness and consistency
            num_reactant_atoms = reactant_combined.GetNumAtoms()
            num_product_atoms = product_combined.GetNumAtoms()
            reactant_mapped = df["reactant_idx"].notna().sum()
            product_mapped = df["product_idx"].notna().sum()

            # Validation 1: Mapped atom counts must match between columns
            if reactant_mapped != product_mapped:
                raise ValueError(
                    "Mismatch in mapped atom counts between columns: "
                    f"reactant_idx mapped={reactant_mapped}, product_idx mapped={product_mapped}"
                )

            # Validation 2: All reactant atoms must be mapped
            if reactant_mapped != num_reactant_atoms:
                df.to_csv(csv_cache / f"debug_mapping_{key}_{total_products}.csv", index=False)
                raise ValueError(
                    "Mapping does not cover all reactant atoms: "
                    f"reactant atoms={num_reactant_atoms}, mapped={reactant_mapped}"
                )

            # Validation 3: All product atoms must be mapped
            if product_mapped != num_product_atoms:
                raise ValueError(
                    "Mapping does not cover all product atoms: "
                    f"product atoms={num_product_atoms}, mapped={product_mapped}"
                )

            # Validation 4: Mapping indices must be consecutive
            if (is_consecutive(df["reactant_idx"].tolist()) is False or 
                is_consecutive(df["product_idx"].tolist()) is False):
                raise ValueError(
                    "Mapping indices are not consecutive: "
                    f"reactant indices={df['reactant_idx'].tolist()}, "
                    f"product indices={df['product_idx'].tolist()}"
                )

            # Apply smart mapping to both reactants
            sub_set1 = is_number_in_set(matches_1, r1)
            smart_mapping(reactant=r1, smarts_template=rxn.GetReactants()[0], 
                         match_tuple=sub_set1 if sub_set1 else ())
            

            sub_set2 = is_number_in_set(matches_2, r2)
            smart_mapping(reactant=r2, smarts_template=rxn.GetReactants()[1], 
                         match_tuple=sub_set2 if sub_set2 else ())

            # Recombine reactants after mapping
            reactant_combined = Chem.CombineMols(r1, r2)
            
            # Identify first shell atoms and initiators
            for atom in reactant_combined.GetAtoms():
                if atom.GetAtomMapNum() <= 99:
                    first_shell.append(atom.GetIdx())
                # we could get the initiators by their map numbers
                if atom.GetAtomMapNum() in [1, 2]:
                    initiator_idxs.append(atom.GetIdx())

            # Identify byproduct atoms if deletion is requested
            byproduct_product_idxs = []
            byproduct_reactant_idxs = []

            if delete_atoms:
                frags = rdmolops.GetMolFrags(product_combined, asMols=True)
                smallest_mol = min(frags, key=lambda m: m.GetNumAtoms())

                # 1) get byproduct atom indices in PRODUCT space
                for atom in smallest_mol.GetAtoms():
                    byproduct_map_number = atom.GetAtomMapNum()
                    for p_atom in product_combined.GetAtoms():
                        if p_atom.GetAtomMapNum() == byproduct_map_number:
                            byproduct_product_idxs.append(p_atom.GetIdx())
                            break

                # 2) convert PRODUCT indices -> REACTANT indices using df mapping
                byproduct_reactant_idxs = (
                    df.loc[df["product_idx"].isin(byproduct_product_idxs), "reactant_idx"]
                    .astype(int)
                    .tolist()
                )

            # if you want the column to store REACTANT indices now:
            byproduct_indexs = byproduct_reactant_idxs

            # Add any additional mappings from mapping_dict
            for r_idx, p_idx in mapping_dict.items():
                new_row = pd.DataFrame([{"reactant_index": r_idx, "product_idx": p_idx}])
                df = pd.concat([df, new_row], ignore_index=True)

            # Create comprehensive DataFrame with all analysis columns
            first_shell_column = pd.Series(first_shell, name="first_shell")
            initiator_idxs_column = pd.Series(initiator_idxs, name="initiators")
            
            # Validation: Must have exactly 2 initiators
            if len(initiator_idxs) != 2:
                raise ValueError(f"Expected 2 initiators, got {len(initiator_idxs)}: {initiator_idxs}")
            
            by_product_indexs_column = pd.Series(byproduct_indexs, name="byproduct_indices")
            
            # Combine all data into final DataFrame
            df_combined = pd.concat([df, first_shell_column, initiator_idxs_column,
                                     by_product_indexs_column], axis=1).astype(pd.Int64Dtype())
            
            # Check for duplicate products
            if not compare_products(molecule_and_csv_path_dict, product_combined):
                continue

            # Update counters and store results
            total_products += 1
            mols.append(reactant_combined)
            mols.append(product_combined)

            # Organize results in the output dictionary
            for i in itertools.count(1):
                if i not in molecule_and_csv_path_dict: # Check if i IS a key
                    sub_dict = molecule_and_csv_path_dict[i] = {} # Initialize if not found
                    dict_key = i
                    break

            # Save to CSV file
            df_combined.to_csv(csv_cache / f"reaction_{dict_key}.csv", index=False)
            print(f"Saved reaction {dict_key} to CSV")

            sub_dict["reactant"] = reactant_combined
            sub_dict["product"] = product_combined
            sub_dict["csv_path"] = csv_cache / f"reaction_{dict_key}.csv"
            sub_dict["delete_atoms"] = delete_atoms
            sub_dict["reaction_dataframe"] = df_combined

    return mols, molecule_and_csv_path_dict


def save_grid_image(mols, cache, key=None):
    from pathlib import Path
    from rdkit.Chem import Draw
    import base64

    out_dir = Path(cache) / "grid_images"
    out_dir.mkdir(parents=True, exist_ok=True)

    out_path = out_dir / (f"reaction_grid_{key}.png" if key is not None else "reaction_grid.png")

    # Ask RDKit for a raster image (PNG-like)
    img = Draw.MolsToGridImage(mols, useSVG=False, subImgSize=(900, 900))

    # Case 1: PIL.Image.Image (has .save)
    if hasattr(img, "save"):
        img.save(str(out_path))
        return out_path

    # Case 2: IPython/RDKit display object (often has .data)
    data = getattr(img, "data", None)
    if data is None:
        raise TypeError(f"Unexpected image type {type(img)}; cannot save.")

    # data can be bytes OR base64 string depending on wrapper
    if isinstance(data, bytes):
        png_bytes = data
    elif isinstance(data, str):
        # try base64 decode; if it fails, treat as raw text
        try:
            png_bytes = base64.b64decode(data)
        except Exception:
            png_bytes = data.encode("utf-8")
    else:
        raise TypeError(f"Unexpected img.data type {type(data)}; cannot save.")

    with open(out_path, "wb") as f:
        f.write(png_bytes)

    return out_path



def reaction_tuples(same_reactants, mol_reactant_1, mol_reactant_2):
    """
    Generate reaction tuples based on whether reactants are the same or different.
    
    Args:
        same_reactants (bool): Whether both reactants are the same molecule
        mol_reactant_1 (Chem.Mol): First reactant molecule
        mol_reactant_2 (Chem.Mol): Second reactant molecule (can be None)
        
    Returns:
        list: List of reactant pairs for reaction processing
        
    Note:
        For same reactants, returns both [A,B] and [B,A] to account for symmetry.
        For different reactants, returns only [A,B].
    """
    if same_reactants:
        mol_reactant_2 = mol_reactant_1

    return [[mol_reactant_1, mol_reactant_2], [mol_reactant_2, mol_reactant_1]]



def process_reaction_dict(detected_reactions, cache):
    """
    Process a dictionary of detected reactions and generate mapping data.
    
    This is a high-level function that processes multiple reactions defined in
    a dictionary structure, handling all steps from SMILES to final CSV output.
    
    Args:
        detected_reactions (dict): Dictionary containing reaction definitions
        cache (str or Path): Base cache directory for output files
        
    Returns:
        tuple: (molecule_and_csv_path_dict, detected_reactions)
            - molecule_and_csv_path_dict: Dictionary with all processed reaction data
            - detected_reactions: The original input dictionary (unchanged)
            
    Example:
        >>> reactions_dict = {
        ...     1: {"reaction": "[C:1]=[O:2]>>[C:1]-[O:2]", 
        ...         "monomer_1": {"smiles": "CCO"}, ...}
        ... }
        >>> results, original = process_reaction_dict(reactions_dict, "/path/to/cache")
    """
    molecule_and_csv_path_dict = {}
    csv_cache = prepare_paths(cache)

    # TODO: Add function to compare products if products are identicals drop the duplicates
    

    # Process each reaction in the dictionary
    for key in detected_reactions:
        reaction_dict = detected_reactions[key]
        rxn_smarts = reaction_dict["reaction"]
        reactant_smiles_1 = reaction_dict["monomer_1"]["smiles"]
        same_reactants = reaction_dict["same_reactants"]
        
        # Handle same vs different reactants
        if same_reactants:
            reactant_smiles_2 = reactant_smiles_1
        else:
            reactant_smiles_2 = reaction_dict["monomer_2"]["smiles"]
        
        delete_atoms = reaction_dict.get("delete_atom", False)
        
        # Build reaction and reactants
        rxn = build_reaction(rxn_smarts)
        mol_reactant_1, mol_reactant_2 = build_reactants(reactant_smiles_1, reactant_smiles_2)
        reaction_tuple = reaction_tuples(same_reactants, mol_reactant_1, mol_reactant_2)

        # Process the reaction
        mols, molecule_and_csv_path_dict = process_reactions(
            rxn, csv_cache, reaction_tuple, key, molecule_and_csv_path_dict, 
            delete_atoms=delete_atoms
        )
        
        # Save outputs
        save_grid_image(mols, cache, key)
    
    return molecule_and_csv_path_dict, detected_reactions


def run_all(cache, rxn_smarts, reactant_smiles_1, reactant_smiles_2):
    """
    Run a complete reaction processing pipeline for a single reaction.
    
    This is a convenience function that handles all steps for processing
    a single reaction from SMILES strings to final outputs.
    
    Args:
        cache (str or Path): Base cache directory for output files
        rxn_smarts (str): Reaction SMARTS string
        reactant_smiles_1 (str): SMILES string for first reactant
        reactant_smiles_2 (str): SMILES string for second reactant
        
    Returns:
        dict: Dictionary containing all processed reaction data
        
    Example:
        >>> results = run_all("/path/to/cache", 
        ...                   "[C:1]=[O:2]>>[C:1]-[O:2]", 
        ...                   "CCO", "CC=O")
    """
    molecule_and_csv_path_dict = {}
    csv_cache = prepare_paths(cache)
    
    # Build reaction and reactants
    rxn = build_reaction(rxn_smarts)
    mol_reactant_1, mol_reactant_2 = build_reactants(reactant_smiles_1, reactant_smiles_2)
    reaction_tuple = [[mol_reactant_1, mol_reactant_2]]
    delete_atoms = True  # Default to identifying byproducts
    
    # Process the reaction
    mols, molecule_and_csv_path_dict = process_reactions(
        rxn, csv_cache, reaction_tuple, None, molecule_and_csv_path_dict, delete_atoms
    )
    
    # Save outputs
    save_grid_image(mols, cache, None)

    return molecule_and_csv_path_dict


if __name__ == "__main__":
    """
    Main execution block with example reaction data.
    
    This demonstrates the usage of the module with example polyesterification
    reactions between hydroxy carboxylic acids.
    """
    
    # Example reaction database for polyesterification reactions
    detected_reactions = {
        1: {
            "reaction_name": "Hydroxy Carboxylic Acid Polycondensation(Polyesterification)",
            "same_reactants": True,
            "reactant_1": "hydroxy_carboxylic_acid",
            "product": "polyester_chain",
            "delete_atom": True,
            "reaction": "[O;!$(OC=*):1]-[H:3].[CX3:2](=[O:5])[OX2H1:4]>>[OX2:1]-[CX3:2](=[O:5]).[O:4]-[H:3]",
            "reference": {
                "smarts": "https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329",
                "reaction_and_mechanism": [
                    "https://pubs.acs.org/doi/10.1021/ed048pA734.1",
                    "https://pubs.acs.org/doi/10.1021/ed073pA312",
                ],
            },
            "monomer_1": {
                "smiles": "O=C(O)c1cc(O)c(Cl)c(C(=O)O)c1",
                "functionality_type": "di_different",
                "functional_group_name": "hydroxy_carboxylic_acid",
                "functional_group_smarts_1": "[OX2H1;!$(OC=*):1]",
                "functional_count_1": 1,
                "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
                "functional_count_2": 2,
            },
        },
        2: {
            "reaction_name": "Hydroxy Carboxylic Acid Polycondensation(Polyesterification)",
            "same_reactants": True,
            "reactant_1": "hydroxy_carboxylic_acid",
            "product": "polyester_chain",
            "delete_atom": True,
            "reaction": "[O;!$(OC=*):1]-[H:3].[CX3:2](=[O:5])[OX2H1:4]>>[OX2:1]-[CX3:2](=[O:5]).[O:4]-[H:3]",
            "reference": {
                "smarts": "https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329",
                "reaction_and_mechanism": [
                    "https://pubs.acs.org/doi/10.1021/ed048pA734.1",
                    "https://pubs.acs.org/doi/10.1021/ed073pA312",
                ],
            },
            "monomer_1": {
                "smiles": "O=C(O)CCCC(O)CCCO",
                "functionality_type": "di_different",
                "functional_group_name": "hydroxy_carboxylic_acid",
                "functional_group_smarts_1": "[OX2H1;!$(OC=*):1]",
                "functional_count_1": 2,
                "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
                "functional_count_2": 1,
            },
        },
        3: {
            "reaction_name": "Hydroxy Carboxylic and Hydroxy Carboxylic Polycondensation(Polyesterification)",
            "same_reactants": False,
            "reactant_1": "hydroxy_carboxylic_acid",
            "reactant_2": "hydroxy_carboxylic_acid",
            "product": "polyester_chain",
            "delete_atom": True,
            "reaction": "[O;!$(OC=*):1]-[H:3].[CX3:2](=[O:5])[OX2H1:4]>>[OX2:1]-[CX3:2](=[O:5]).[O:4]-[H:3]",
            "reference": {
                "smarts": "https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329",
                "reaction_and_mechanism": [
                    "https://pubs.acs.org/doi/10.1021/ed048pA734.1",
                    "https://pubs.acs.org/doi/10.1021/ed073pA312",
                ],
            },
            "monomer_1": {
                "smiles": "O=C(O)c1cc(O)c(Cl)c(C(=O)O)c1",
                "functionality_type": "di_different",
                "functional_group_name": "hydroxy_carboxylic_acid",
                "functional_group_smarts_1": "[OX2H1;!$(OC=*):1]",
                "functional_count_1": 1,
                "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
                "functional_count_2": 2,
            },
            "monomer_2": {
                "smiles": "O=C(O)CCCC(O)CCCO",
                "functionality_type": "di_different",
                "functional_group_name": "hydroxy_carboxylic_acid",
                "functional_group_smarts_1": "[OX2H1;!$(OC=*):1]",
                "functional_count_1": 2,
                "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
                "functional_count_2": 1,
            },
        },
    }

    # Configuration parameters
    cache = "C:\\Users\\Janitha\\Documents\\GitHub\\reaction_lammps_mupt\\cache\\00_cache"
    rxn_smarts = "[O;!$(OC=*):1]-[H:3].[CX3:2](=[O:5])[OX2H1:4]>>[OX2:1]-[CX3:2](=[O:5]).[O:4]-[H:3]"
    reactant_smiles_1 = "O=C(O)c1cc(O)c(Cl)c(C(=O)O)c1"
    reactant_smiles_2 = "O=C(O)CCCC(O)CCCO"

    # Execute the reaction processing pipeline
    print("Starting reaction processing pipeline...")
    
    # Process a single reaction using run_all
    # run_all(cache, rxn_smarts, reactant_smiles_1, reactant_smiles_2)
    
    # Process all reactions in the dictionary
    molecule_dict_csv_path_dict, detected_reactions = process_reaction_dict(detected_reactions, cache)
    
    # Display results
    import pprint
    print("\n" + "="*80)
    print("REACTION PROCESSING RESULTS")
    print("="*80)
    pprint.pprint(molecule_dict_csv_path_dict)
    
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Total reactions processed: {len(molecule_dict_csv_path_dict)}")
    
    # Count total products across all reactions
    print("\nProcessing complete! CSV files and images have been saved to the cache directory.")
