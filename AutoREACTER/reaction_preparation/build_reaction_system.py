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

"""
TODO: Preserve stable reaction_id for each reaction (never renumber).
Maintain a seq_id mapping for the current filtered set (continuous 1..N) for UI/export.
Keep on-disk filenames stable as reaction_{reaction_id}.csv.
If sequential filenames are needed, generate them only in an export step, and store both reaction_id and seq_id in the CSV/dataframe.
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

from AutoREACTER.detectors.reaction_detector import ReactionInstance, ReactionTemplate
from AutoREACTER.detectors.functional_groups_detector import FunctionalGroupInfo, MonomerRole
from AutoREACTER.input_parser import SimulationSetup, MonomerEntry

# standard libs / third-party
import pandas as pd
from pathlib import Path
import os
from dataclasses import dataclass, field
from typing import List, Any, Optional, Dict
from rdkit import Chem



@dataclass (slots=True)
class ReactionMetadata:
    """
    A dataclass to hold all relevant metadata and file paths for a single reaction instance.
    Attributes:
    - reaction_id: Unique identifier for the reaction.
    - reactant_combined_mol: RDKit molecule object representing the combined reactants.
    - product_combined_mol: RDKit molecule object representing the combined products.
    - reaction_smarts: Optional SMARTS string representing the reaction.
    - reactant_smarts: Optional SMARTS string for the reactants.
    - product_smiles: Optional SMILES string for the products.
    - csv_path: Path to the CSV file containing processing reaction data.
    - mol_3d_path: Optional path to the generated 3D molecule file.
    - reaction_dataframe: Optional pandas DataFrame holding detailed reaction mapping and analysis data.
    - reactant_to_product_mapping: Dictionary mapping reactant atom indices to product atom indices.
    - template_reactant_to_product_mapping: Dictionary mapping reactant atom indices to product atom indices specifically for the reaction template.
    - edge_atoms: List of atom indices that are considered 'edge atoms' connecting the reaction template to the rest of the molecule.
    - delete_atom: Boolean flag indicating whether an atom is deleted in the reaction (e.g., for polycondensation).
    - delete_atom_idx: Optional index of the atom that is deleted in the reaction, if applicable.
    - activity_stats: Boolean flag indicating whether to calculate and include activity statistics for the reaction.
    """
    reaction_id: int
    reactant_combined_mol: Chem.Mol
    product_combined_mol: Chem.Mol
    reaction_smarts: Optional[str] = None
    reactant_smarts: Optional[str] = None
    product_smiles: Optional[str] = None
    csv_path: Path = None
    mol_3d_path: Optional[Path] = None
    reaction_dataframe: Optional[pd.DataFrame] = None
    reactant_to_product_mapping: Dict[int, int]
    template_reactant_to_product_mapping: Dict[int, int] 
    edge_atoms: List[int] 
    delete_atom: bool = True
    delete_atom_idx: Optional[int] = None
    activity_stats: bool = True
    

 
@dataclass (slots=True)
class PipelineOutput:
    """
    Final structured output of the ReactionTemplatePipeline.
    Attributes:
    - template_files: Dictionary mapping reaction identifiers to their corresponding molecule template file paths.
    - all_reactions: Optional list of ReactionMetadata instances containing detailed information about each processed reaction
    """
    template_files: Dict[str, Path]
    all_reactions: Optional[List[ReactionMetadata]] = None

# The main class encapsulating the reaction template pipeline
class ReactionTemplatePipeline:
    """
    A class to encapsulate the entire reaction template pipeline.

    This class provides methods to process detected reactions, identify templates,
    and prepare molecule files for 3D generation and LUNAR workflow execution.
    It serves as a structured way to organize the complex sequence of steps involved
    in the reaction template preparation process.
    """
    def __init__(self, detected_reactions_dict: dict, cache: str, non_monomer_molecules_to_retain: list):
        """
        Initializes the ReactionTemplatePipeline with the given detected reactions and cache directory.

        Args:
            detected_reactions_dict (dict): A dictionary containing detected reaction data.
            cache (str): Path to the cache directory for processing.
            non_monomer_molecules_to_retain (list): List of SMILES strings for non-monomer molecules to keep.
        """
        self.detected_reactions_dict = detected_reactions_dict
        self.cache = cache
        self.non_monomer_molecules_to_retain = non_monomer_molecules_to_retain
        self.formatted_dict, self.molecule_template_files, self.molecule_images = self.execute_pipeline(self.detected_reactions_dict, self.non_monomer_molecules_to_retain, self.cache)

    def _add_column_safe(self, df, list_data, column_name):
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
    
    def _add_dict_as_new_columns(self, df_existing, data_dict, titles=["template_reactant_idx", "template_product_idx"]):
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
    
    def run_reaction_template_pipeline(self, detected_reactions, cache):
        """
        Main execution pipeline for mapping reactions, identifying templates, 
        and saving results to CSV.

        Steps:
        1. Maps atoms between reactants and products.
        2. Uses an atom walker to identify the 'reaction center' (template).
        3. Identifies 'edge atoms' that connect the template to the rest of the molecule.
        4. Updates reaction dataframes with mapping indices and saves them.

        Args:
            detected_reactions_dict (dict): Raw dictionary of detected reactions.
            cache (str): Path to the cache directory for processing.

        Returns:
            tuple: (updated_molecule_dict, formatted_summary_dict, molecule_images)
        """
        # Initial processing: atom mapping and basic dictionary formatting
        molecule_dict_csv_path_dict, detected_reactions, molecule_images = process_reaction_dict(self.detected_reactions_dict, self.cache)
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
            reaction_dataframe = self._add_dict_as_new_columns(
                reaction_dataframe,
                template_mapped_dict,
                titles = ["template_reactant_idx", "template_product_idx"]
            )
            
            # Append edge atoms information (crucial for building polymer chains)
            reaction_dataframe = self._add_column_safe(
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
        return molecule_dict_csv_path_dict, formatted_dict, molecule_images


    def execute_pipeline(self, detected_reactions, non_monomer_molecules_to_retain, cache):
        """
        Orchestrates the full reaction template pipeline, including 3D generation 
        and LUNAR workflow execution.

        Args:
            detected_reactions (dict): Input reaction definitions.
            non_monomer_molecules_to_retain (list): List of SMILES strings for non-monomer molecules to keep.
            cache (str): Directory path for temporary files and outputs.

        Returns:
            tuple: (formatted_dict, molecule_template_files)
        """
        # 1. Run the core reaction template mapping and walking pipeline
        molecule_dict_csv_path_dict, formatted_dict, molecule_images = self.run_reaction_template_pipeline(
            detected_reactions, 
            cache
        )
        
        # 2. Re-format data and extract a list of unique SMILES for 3D generation
        formatted_dict, data_smiles_list = format_detected_reactions_dict(
            detected_reactions, 
            non_monomer_molecules_to_retain
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
        
        with open(Path(cache) / "molecule_template_files.txt", "w") as f:
            for name, path in molecule_template_files.items():
                f.write(f"{name}: {path}\n")
            

        return formatted_dict, molecule_template_files, molecule_images
        # TODO: If duplicates were found and skipped, re-index the dictionary and files to be continuous
        # if duplicated:
        #     molecule_dict_csv_path_dict = molecule_dict_csv_path_dict_rearrange(molecule_dict_csv_path_dict)


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
    pipeline = ReactionTemplatePipeline(
        detected_reactions_dict=detected_reactions_dict, 
        cache=cache_path, 
        non_monomer_molecules_to_retain=non_monomer_molecules_to_retain
    )
    print(pipeline.formatted_dict)
    print(pipeline.molecule_template_files)
