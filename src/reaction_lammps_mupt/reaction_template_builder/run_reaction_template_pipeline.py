"""
Reaction Template Pipeline Module
---------------------------------
This module provides a comprehensive pipeline to process chemical reaction data. 
It handles atom mapping, reaction template identification via atom walking, 
structural deduplication, 3D molecule generation, and integration with the 
LUNAR workflow for molecular template preparation.

Main features:
1. Atom mapping and template extraction.
2. Sequential indexing and file management for reaction data.
3. Integration with RDKit for fragment comparison.
4. 3D structure generation and workflow automation.

TODO:
- Implement duplicate reaction detection and handling.
"""

# --- preferred: package-relative imports; fallback: absolute imports for script/notebook ---

try:
    # Case 1 (correct when imported as part of the package)
    from .reaction_template_pipeline.map_reactant_atoms import process_reaction_dict
    from .reaction_template_pipeline.util import (
        format_detected_reactions_dict,
        prep_for_3d_molecule_generation,
    )
    from .reaction_template_pipeline.walker import reaction_atom_walker
    from .reaction_template_pipeline.compare_rdkit_fragments import compare_rdkit_fragments

    from .lunar_client.molecule_3d_preparation import prepare_3d_molecule
    from .lunar_client.lunar_api_wrapper import lunar_workflow
    from .lunar_client.molecule_template_preparation import molecule_template_preparation

except (ImportError, ModuleNotFoundError):
    try:
        # Case 2 (works if running as a loose script from inside reaction_template_builder/)
        from reaction_template_pipeline.map_reactant_atoms import process_reaction_dict
        from reaction_template_pipeline.util import (
            format_detected_reactions_dict,
            prep_for_3d_molecule_generation,
        )
        from reaction_template_pipeline.walker import reaction_atom_walker
        from reaction_template_pipeline.compare_rdkit_fragments import compare_rdkit_fragments

        from lunar_client.molecule_3d_preparation import prepare_3d_molecule
        from lunar_client.lunar_api_wrapper import lunar_workflow
        from lunar_client.molecule_template_preparation import molecule_template_preparation

    except (ImportError, ModuleNotFoundError):
        # Case 3 (fully-qualified import from repo root / installed package)
        from reaction_lammps_mupt.reaction_template_builder.reaction_template_pipeline.map_reactant_atoms import (
            process_reaction_dict,
        )
        from reaction_lammps_mupt.reaction_template_builder.reaction_template_pipeline.util import (
            format_detected_reactions_dict,
            prep_for_3d_molecule_generation,
        )
        from reaction_lammps_mupt.reaction_template_builder.reaction_template_pipeline.walker import (
            reaction_atom_walker,
        )
        from reaction_lammps_mupt.reaction_template_builder.reaction_template_pipeline.compare_rdkit_fragments import (
            compare_rdkit_fragments,
        )

        from reaction_lammps_mupt.reaction_template_builder.lunar_client.molecule_3d_preparation import (
            prepare_3d_molecule,
        )
        from reaction_lammps_mupt.reaction_template_builder.lunar_client.lunar_api_wrapper import (
            lunar_workflow,
        )
        from reaction_lammps_mupt.reaction_template_builder.lunar_client.molecule_template_preparation import (
            molecule_template_preparation,
        )


# standard libs / third-party
import pandas as pd
from pathlib import Path
import os

def extract_unique_references(detected_reactions: dict) -> list[str]:
    """
    Collect ordered, deduplicated reference URLs from a mapping of detected reactions.
    
    Parameters:
        detected_reactions (dict): Mapping of reaction identifiers to reaction dictionaries. Each reaction may contain a "reference" entry whose values can be a string, a list/tuple of strings, or a nested dict of strings.
    
    Returns:
        list[str]: Ordered list of unique reference URLs, preserving the first-seen occurrence for each URL.
    """
    seen = set()
    refs: list[str] = []

    for rxn in detected_reactions.values():
        ref = rxn.get("reference") or {}

        # single URL fields
        for v in ref.values():
            if isinstance(v, str):
                if v not in seen:
                    seen.add(v)
                    refs.append(v)

            # list-of-URLs fields
            elif isinstance(v, (list, tuple)):
                for u in v:
                    if isinstance(u, str) and u not in seen:
                        seen.add(u)
                        refs.append(u)

            # (optional) nested dict handling, if you ever add that later
            elif isinstance(v, dict):
                for u in v.values():
                    if isinstance(u, str) and u not in seen:
                        seen.add(u)
                        refs.append(u)

    return refs


def save_grid_image(mols, cache, key=None):
    """
    Save a grid image of the provided RDKit molecule objects to the cache's grid_images directory.
    
    Creates a PNG file named "reaction_grid.png" or "reaction_grid_{key}.png" inside a "grid_images" subdirectory of `cache`, rendering `mols` as a grid image and returning the file path.
    
    Parameters:
        mols: Iterable of RDKit molecule objects to render into a grid.
        cache (str or pathlib.Path): Directory where the "grid_images" subdirectory will be created.
        key (str, optional): Optional identifier appended to the filename to avoid collisions.
    
    Returns:
        pathlib.Path: Path to the written PNG image file.
    
    Raises:
        TypeError: If the RDKit drawing result is an unexpected object type or its data is neither bytes nor a decodable base64/string that can be written as PNG.
    """
    from pathlib import Path
    from rdkit.Chem import Draw
    import base64

    out_dir = Path(cache) / "grid_images"
    out_dir.mkdir(parents=True, exist_ok=True)

    out_path = out_dir / (f"reaction_grid_{key}.png" if key is not None else "reaction_grid.png")

    # Ask RDKit for a raster image (PNG-like)
    img = Draw.MolsToGridImage(mols, useSVG=False)

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

def is_continuous(d):
    """
    Checks if the integer keys of a dictionary form a continuous sequence starting from 1.

    This is used to verify if reaction IDs are sequential after potential 
    filtering or deduplication steps.

    Args:
        d (dict): The dictionary to check, where keys are expected to be integers.

    Returns:
        bool: True if keys are exactly [1, 2, ..., len(d)], False otherwise.
    """
    keys = sorted(d.keys())
    if not keys: 
        return True
    # Compare sorted keys against a generated range of the same length
    return keys == list(range(1, len(keys) + 1))

def molecule_dict_csv_path_dict_rearrange(molecule_dict_csv_path_dict):
    """
    Normalize reaction keys to a continuous 1..N sequence and rename their CSV files accordingly.
    
    Parameters:
        molecule_dict_csv_path_dict (dict): Mapping of integer reaction keys to metadata dicts.
            Each metadata dict is expected to include a "csv_path" entry pointing to the reaction's CSV file.
    
    Returns:
        dict: New mapping with consecutive integer keys starting at 1 and metadata dictionaries whose
            "csv_path" values have been updated to match any renamed files.
    """
    # Skip if already continuous to save processing time
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
            # Construct the new filename based on the new index (e.g., reaction_2.csv)
            new_csv_save_path = os.path.join(
                os.path.dirname(csv_save_path),
                f"reaction_{new_key}.csv"
            )
            
            # Check if destination exists to avoid crashing or accidental overwrites
            if os.path.exists(new_csv_save_path):
                os.remove(new_csv_save_path)
                
            # Rename the physical file on the disk to match the new key
            os.rename(csv_save_path, new_csv_save_path)
            
            # Update the path in the metadata dictionary
            reaction_data["csv_path"] = new_csv_save_path
            
        new_dict[new_key] = reaction_data

    return new_dict

def add_dict_as_new_columns(df_existing, data_dict, titles=["template_reactant_idx", "template_product_idx"]):
    """
    Adds dictionary keys and values as new columns to an existing DataFrame.

    This is specifically used to map reactant indices to product indices 
    within the reaction template.

    Args:
        df_existing (pd.DataFrame): The DataFrame to modify.
        data_dict (dict): Dictionary where keys and values will become column data.
        titles (list): List of two strings for the new column names.

    Returns:
        pd.DataFrame: The modified DataFrame with new columns added.
    """
    # Convert dict keys and values to Series to ensure alignment with the DataFrame
    # Use .astype("Int64") to allow for potential Null/NaN values while keeping integers
    df_existing[titles[0]] = pd.Series(list(data_dict.keys())).astype("Int64")
    df_existing[titles[1]] = pd.Series(list(data_dict.values())).astype("Int64")
    
    return df_existing

def add_column_safe(df, list_data, column_name):
    """
    Safely adds a list as a new column to a DataFrame, handling potential 
    length mismatches by using Series alignment.

    Args:
        df (pd.DataFrame): Target DataFrame.
        list_data (list): Data to be added to the column.
        column_name (str): Name of the new column.

    Returns:
        pd.DataFrame: Modified DataFrame with the new column.
    """
    # Creating a Series from the list ensures it starts from the top (index 0)
    # and fills missing rows with NaN if the list is shorter than the DataFrame
    df[column_name] = pd.Series(list_data).astype("Int64")
    return df

def run_reaction_template_pipeline(detected_reactions_dict, cache):
    """
    Orchestrates atom mapping, reaction-template extraction, and persistence of processed reaction data to CSV.
    
    Parameters:
        detected_reactions_dict (dict): Raw detected-reactions input; keys map to reaction entries containing reactant/product objects and mapping data.
        cache (str): Filesystem path to a cache directory used for intermediate files and output CSVs.
    
    Returns:
        tuple: (molecule_dict_csv_path_dict, formatted_dict)
            molecule_dict_csv_path_dict (dict): Updated mapping of reaction keys to reaction records (each record includes the updated DataFrame and its `csv_path`).
            formatted_dict (dict): A human-readable summary of the detected reactions after initial processing.
    """
    # Initial processing: atom mapping and basic dictionary formatting
    molecule_dict_csv_path_dict, detected_reactions = process_reaction_dict(detected_reactions_dict, cache)
    formatted_dict = format_detected_reactions_dict(detected_reactions)
    
    # Iterate through each detected reaction to perform template analysis
    for key, reaction in molecule_dict_csv_path_dict.items():
        combined_reactant_molecule_object = reaction.get("reactant")
        combined_product_molecule_object = reaction.get("product")
        reaction_dataframe = reaction.get("reaction_dataframe")
        csv_save_path = reaction.get("csv_path")
        
        # Create a mapping dictionary {reactant_idx: product_idx} from the dataframe
        fully_mapped_dict = reaction_dataframe.set_index("reactant_idx")["product_idx"].to_dict()
        
        # Extract the 'first shell' (atoms directly involved in bond changes)
        first_shell = reaction_dataframe['first_shell'].dropna().tolist()
        
        # Perform atom walking to determine the reaction template and boundary (edge) atoms
        template_mapped_dict, edge_atoms = reaction_atom_walker(
            combined_reactant_molecule_object,
            first_shell,
            fully_mapped_dict
        )
        
        # Note: Duplicate checking logic is currently commented out but available for future use
        # is_duplicate, processed_dict = compare_rdkit_fragments(...)

        # Update the dataframe with specific template mapping indices
        reaction_dataframe = add_dict_as_new_columns(
            reaction_dataframe,
            template_mapped_dict,
            titles = ["template_reactant_idx", "template_product_idx"]
        )
        
        # Append edge atoms information (crucial for building polymer chains)
        reaction_dataframe = add_column_safe(
            reaction_dataframe,
            edge_atoms,
            "edge_atoms"
        )
        
        # Store the updated dataframe back in the reaction object
        # Using .copy() to avoid SettingWithCopyWarning in pandas
        reaction["reaction_dataframe"] = reaction_dataframe.copy() 
        
        # Persist the processed reaction data to a CSV file
        reaction_dataframe.to_csv(csv_save_path, index=False)

    # For debugging: print the final molecule dictionary with CSV paths
    # print("Formatted Detected Reactions Summary:")
    # print(molecule_dict_csv_path_dict)
    return molecule_dict_csv_path_dict, formatted_dict


def execute_pipeline(detected_reactions, retain_smiles, cache):
    """
    Orchestrates the reaction template pipeline, producing 3D structures, running the LUNAR workflow, and writing molecule template files.
    
    Parameters:
        detected_reactions (dict): Detected reaction definitions and associated metadata.
        retain_smiles (list): SMILES strings for non-monomer molecules that must be retained for 3D generation.
        cache (str or Path): Directory path used for intermediate files, outputs, and cache storage.
    
    Returns:
        tuple: 
            formatted_dict (dict): Reformatted representation of detected reactions prepared for downstream use.
            molecule_template_files (dict): Mapping from molecule identifier to the generated molecule template file path.
    """
    # 1. Run the core reaction template mapping and walking pipeline
    molecule_dict_csv_path_dict, formatted_dict = run_reaction_template_pipeline(
        detected_reactions, 
        cache
    )
    
    # 2. Re-format data and extract a list of unique SMILES for 3D generation
    formatted_dict, data_smiles_list = format_detected_reactions_dict(
        detected_reactions, 
        retain_smiles
    )
    print("Unique SMILES List:", data_smiles_list)

    # 3. Prepare the molecule dictionary structure required for 3D generation
    molecule_dict = prep_for_3d_molecule_generation(
        data_smiles_list, 
        molecule_dict_csv_path_dict
    )
    
    # 4. Generate 3D structures (MOL files) for all involved molecules
    cache_mol = Path(cache) / "mol_files"
    prepared_molecules = prepare_3d_molecule(
        cache_dir=cache_mol, 
        molecule_dict=molecule_dict
    )
    # For debugging: print the prepared 3D molecule file paths
    # print("Prepared 3D Molecules:", prepared_molecules)

    # 5. Execute the LUNAR workflow (e.g., force field assignment, energy minimization)
    lunar_out_loc_dict = lunar_workflow(molecule_files=prepared_molecules, cache_dir=cache)
    
    # 6. Final Step: Generate molecule template files combining 3D data and reaction maps
    molecule_template_files = molecule_template_preparation(
            molecule_dict_csv_path_dict,
            lunar_out_loc_dict,
            cache
        )   
    
    # Log the location of generated template files
    for name, path in molecule_template_files.items():
        print(f"Molecule Template File for {name}: {path}")

    # Extract and save unique references to a text file for user reference
    references = extract_unique_references(detected_reactions)

    # Print and save references to a text file
    print("\nReferences used in detected reactions:")
    for ref in references:
        print(ref)
    
    with open(Path(cache) / "molecule_template_files.txt", "w") as f:
        for name, path in molecule_template_files.items():
            f.write(f"{name}: {path}\n")
        f.write("\nReferences:\n")
        for ref in references:
            f.write(f"{ref}\n")

    return formatted_dict, molecule_template_files
    # If duplicates were found and skipped, re-index the dictionary and files to be continuous
    if duplicated:
        molecule_dict_csv_path_dict = molecule_dict_csv_path_dict_rearrange(molecule_dict_csv_path_dict)


if __name__ == "__main__":
    # --- Example configuration for testing the pipeline ---
    # This dictionary defines several polycondensation reactions with monomer details
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
    
    # Define molecules that should be retained in the system (e.g., solvents or catalysts)
    non_monomer_molecules_to_retain = ["CCO"]
    
    # Path to the cache directory for storing intermediate files
    cache_path = r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\00_cache"

    # --- Execution Start ---
    execute_pipeline(
        detected_reactions=detected_reactions_dict,
        retain_smiles=non_monomer_molecules_to_retain,
        cache=cache_path
    )