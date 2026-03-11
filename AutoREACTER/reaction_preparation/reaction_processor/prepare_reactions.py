
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional
import pandas as pd
from rdkit import Chem

from AutoREACTER.detectors.reaction_detector import ReactionInstance, ReactionTemplate
from AutoREACTER.detectors.functional_groups_detector import FunctionalGroupInfo, MonomerRole

@dataclass (slots=True)
class ReactionMetadata:
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
    
class PrepareReactions:
    def __init__(self, reaction_templates: List[ReactionInstance]):
        self.reaction_templates = reaction_templates

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
    
    def prepare_reactions(self, reactions: List[ReactionInstance]) -> List[ReactionMetadata]:
        for reaction in reactions:
            # Placeholder for the actual preparation logic
            # This is where you would implement the logic to process each reaction,
            # generate the combined mols, create the DataFrame, and populate the metadata.
            pass
