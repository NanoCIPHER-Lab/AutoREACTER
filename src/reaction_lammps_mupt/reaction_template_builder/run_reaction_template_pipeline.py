"""
Reaction Template Pipeline Module
---------------------------------
This module provides functionality to process chemical reaction dictionaries, 
perform atom mapping, identify reaction templates via atom walking, 
and handle deduplication of reaction fragments.
"""

from reaction_template_pipeline.map_reactant_atoms import process_reaction_dict
from reaction_template_pipeline.util import format_detected_reactions_dict, prep_for_3d_molecule_generation
from reaction_template_pipeline.walker import reaction_atom_walker
from reaction_template_pipeline.compare_rdkit_fragments import compare_rdkit_fragments
from lunar_client.molecule_3d_preparation import prepare_3d_molecule
from lunar_client.lunar_api_wrapper import lunar_workflow
import pandas as pd
from pathlib import Path
import os

def is_continuous(d):
    """
    Checks if the integer keys of a dictionary form a continuous sequence starting from 1.

    Args:
        d (dict): The dictionary to check.

    Returns:
        bool: True if keys are [1, 2, ..., len(d)], False otherwise.
    """
    keys = sorted(d.keys())
    if not keys: 
        return True
    return keys == list(range(1, len(keys) + 1))

def molecule_dict_csv_path_dict_rearrange(molecule_dict_csv_path_dict):
    """
    Rearranges the reaction dictionary keys to be continuous and renames 
    associated CSV files on the filesystem to match the new keys.

    This is typically called after duplicates are removed to ensure the 
    output file numbering remains sequential.

    Args:
        molecule_dict_csv_path_dict (dict): Dictionary containing reaction data and file paths.

    Returns:
        dict: A new dictionary with normalized continuous keys (1 to N).
    """
    if is_continuous(molecule_dict_csv_path_dict):
        return molecule_dict_csv_path_dict
        
    print("Rearranging molecule_dict_csv_path_dict keys to be continuous...")
    new_dict = {}
    
    # Sort keys to ensure we process them in the original numerical order
    sorted_old_keys = sorted(molecule_dict_csv_path_dict.keys())
    
    for new_key_idx, old_key in enumerate(sorted_old_keys):
        new_key = new_key_idx + 1  # Start new keys from 1
        reaction_data = molecule_dict_csv_path_dict[old_key]
        csv_save_path = reaction_data.get("csv_path")
        
        if new_key != old_key:
            # Construct the new filename based on the new index
            new_csv_save_path = os.path.join(
                os.path.dirname(csv_save_path),
                f"reaction_{new_key}.csv"
            )
            
            # Check if destination exists to avoid crashing or accidental overwrites
            if os.path.exists(new_csv_save_path):
                os.remove(new_csv_save_path)
                
            # Rename the physical file on the disk
            os.rename(csv_save_path, new_csv_save_path)
            # Update the path in the metadata
            reaction_data["csv_path"] = new_csv_save_path
            
        new_dict[new_key] = reaction_data

    return new_dict

def add_dict_as_new_columns(df_existing, data_dict, titles=["template_reactant_idx", "template_product_idx"]):
    """
    Adds dictionary keys and values as new columns to an existing DataFrame.

    Args:
        df_existing (pd.DataFrame): The DataFrame to modify.
        data_dict (dict): Dictionary where keys and values will become column data.
        titles (list): List of two strings for the new column names.

    Returns:
        pd.DataFrame: The modified DataFrame with new columns.
    """
    # Use pd.Series to ensure data starts at the top and matches indices
    # Use .astype("Int64") to keep integers as whole numbers despite potential <NA> values
    df_existing[titles[0]] = pd.Series(list(data_dict.keys())).astype("Int64")
    df_existing[titles[1]] = pd.Series(list(data_dict.values())).astype("Int64")
    
    return df_existing

def add_column_safe(df, list_data, column_name):
    """
    Safely adds a list as a new column to a DataFrame, handling potential 
    length mismatches by using Series alignment.

    Args:
        df (pd.DataFrame): Target DataFrame.
        list_data (list): Data to be added.
        column_name (str): Name of the new column.

    Returns:
        pd.DataFrame: Modified DataFrame.
    """
    # This adds the list as a column starting from the top
    # .astype("Int64") prevents integer values (e.g. 18) from becoming floats (18.0)
    df[column_name] = pd.Series(list_data).astype("Int64")
    return df

def run_reaction_template_pipeline(detected_reactions_dict, cache):
    """
    Main execution pipeline for mapping reactions, identifying templates, 
    detecting duplicates, and saving results to CSV.

    Args:
        detected_reactions_dict (dict): Raw dictionary of detected reactions.
        cache (str): Path to the cache directory for processing.

    Returns:
        tuple: (updated_molecule_dict, formatted_summary_dict)
    TODO: Handle duplicates by re-indexing and renaming files accordingly.
    """
    duplicated = False
    processed_dict = {}
    
    # Initial processing and atom mapping
    molecule_dict_csv_path_dict, detected_reactions = process_reaction_dict(detected_reactions_dict, cache)
    formatted_dict = format_detected_reactions_dict(detected_reactions)
    
    for key, reaction in molecule_dict_csv_path_dict.items():
        combined_reactant_molecule_object = reaction.get("reactant")
        combined_product_molecule_object = reaction.get("product")
        reaction_dataframe = reaction.get("reaction_dataframe")
        csv_save_path = reaction.get("csv_path")
        
        # Create a mapping dictionary from the dataframe
        fully_mapped_dict = reaction_dataframe.set_index("reactant_idx")["product_idx"].to_dict()
        first_shell = reaction_dataframe['first_shell'].dropna().tolist()
        
        # Perform atom walking to determine the reaction template and edge atoms
        template_mapped_dict, edge_atoms = reaction_atom_walker(
            combined_reactant_molecule_object,
            first_shell,
            fully_mapped_dict
        )
        
        # # Check for structural duplicates using RDKit fragments
        # is_duplicate, processed_dict = compare_rdkit_fragments(
        #     processed_dict,
        #     combined_reactant_molecule_object,
        #     combined_product_molecule_object,
        #     template_mapped_dict
        # )
        
        # if is_duplicate:
        #     duplicated = True
        #     print(f"Duplicate reaction found for reaction ID {key}, skipping saving CSV.")
        #     continue

        # Update the dataframe with template mapping indices
        reaction_dataframe = add_dict_as_new_columns(
            reaction_dataframe,
            template_mapped_dict,
            titles = ["template_reactant_idx", "template_product_idx"]
        )
        
        # Add edge atoms information
        reaction_dataframe = add_column_safe(
            reaction_dataframe,
            edge_atoms,
            "edge_atoms"
        )
        
        # Store the updated dataframe back in the reaction object and save to disk
        reaction["reaction_dataframe"] = reaction_dataframe.copy() # To avoid SettingWithCopyWarning
        reaction_dataframe.to_csv(csv_save_path, index=False)
        
    print("Formatted Detected Reactions Summary:")
    print(molecule_dict_csv_path_dict)
    
    # If duplicates were found and skipped, re-index the dictionary and files to be continuous
    # if duplicated:
    #     molecule_dict_csv_path_dict = molecule_dict_csv_path_dict_rearrange(molecule_dict_csv_path_dict)
    #     del molecule_dict_csv_path_dict[key] # Remove the last duplicate entry
    
    return molecule_dict_csv_path_dict , formatted_dict

def execute_pipeline(detected_reactions, retain_smiles, cache):
    """
    Executes the full reaction template pipeline, 3D molecule generation, 
    and LUNAR workflow.
    """
    # 1. Run the main reaction template pipeline
    molecule_dict_csv_path_dict, formatted_dict = run_reaction_template_pipeline(
        detected_reactions, 
        cache
    )
    
    # 2. Re-format and extract unique SMILES
    formatted_dict, data_smiles_list = format_detected_reactions_dict(
        detected_reactions, 
        retain_smiles
    )
    print("Unique SMILES List:", data_smiles_list)

    # 3. Prepare dictionary for 3D generation
    molecule_dict = prep_for_3d_molecule_generation(
        data_smiles_list, 
        molecule_dict_csv_path_dict
    )
    
    # 4. Generate 3D structures
    cache_mol = Path(cache) / "mol_files"
    prepared_molecules = prepare_3d_molecule(
        cache_dir=cache_mol, 
        molecule_dict=molecule_dict
    )
    print("Prepared 3D Molecules:", prepared_molecules)
    return
    # 5. Execute final LUNAR workflow
    final_files = lunar_workflow(molecule_files=prepared_molecules, cache_dir=cache)
    print("Final LUNAR Workflow Files:", final_files)
    
    return final_files, formatted_dict


if __name__ == "__main__":
    # Example configuration for testing the pipeline
    detected_reactions_dict = {
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
    
    # Configuration for retained molecules and cache path
    non_monomer_molecules_to_retain = ["CCO"]
    
    # Execute the pipeline
    cache_path = r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\00_cache"

    # --- Execution ---
    execute_pipeline(
        detected_reactions=detected_reactions_dict,
        retain_smiles=non_monomer_molecules_to_retain,
        cache=cache_path
    )