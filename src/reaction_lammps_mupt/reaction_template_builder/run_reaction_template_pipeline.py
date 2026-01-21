from reaction_template_pipeline.map_reactant_atoms import process_reaction_dict
from reaction_template_pipeline.util import format_detected_reactions_dict
from reaction_template_pipeline.walker import reaction_atom_walker

import pandas as pd

import pandas as pd

def add_dict_as_new_columns(df_existing, data_dict, titles=["template_reactant_idx", "template_product_idx"]):
    """
    Adds dictionary keys and values as new columns starting from row 0.
    """
    # Use pd.Series to ensure data starts at the top and matches indices
    # Use .astype("Int64") to keep integers as whole numbers despite <NA>
    df_existing[titles[0]] = pd.Series(list(data_dict.keys())).astype("Int64")
    df_existing[titles[1]] = pd.Series(list(data_dict.values())).astype("Int64")
    
    return df_existing

def add_column_safe(df, list_data, column_name):
    # This adds the list as a column starting from the top
    # .astype("Int64") prevents 18 from becoming 18.0
    df[column_name] = pd.Series(list_data).astype("Int64")
    return df

def run_reaction_template_pipeline(detected_reactions_dict, cache):
    molecule_dict_csv_path_dict, detected_reactions = process_reaction_dict(detected_reactions_dict, cache)
    formatted_dict = format_detected_reactions_dict(detected_reactions)
    for key, reaction in molecule_dict_csv_path_dict.items():
        combined_reactant_molecule_object = reaction.get("reactant")
        reaction_dataframe = reaction.get("reaction_dataframe")
        csv_save_path = reaction.get("csv_path")
        fully_mapped_dict = reaction_dataframe.set_index("reactant_idx")["product_idx"].to_dict()
        first_shell = reaction_dataframe['first_shell'].dropna().tolist()
        
        template_mapped_dict, edge_atoms = reaction_atom_walker(
            combined_reactant_molecule_object,
            first_shell,
            fully_mapped_dict
        )
        
        # Now updates columns in place rather than appending rows
        reaction_dataframe = add_dict_as_new_columns(
            reaction_dataframe,
            template_mapped_dict,
            titles = ["template_reactant_idx", "template_product_idx"]
        )
        
        reaction_dataframe = add_column_safe(
            reaction_dataframe,
            edge_atoms,
            "edge_atoms"
        )
        
        reaction["reaction_dataframe"] = reaction_dataframe
        reaction_dataframe.to_csv(csv_save_path, index=False)
    return molecule_dict_csv_path_dict

if __name__ == "__main__":
    
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
    Non_monomer_molecules_to_retain = ["CCO"]
    cache = "C:\\Users\\Janitha\\Documents\\GitHub\\reaction_lammps_mupt\\cache\\00_cache"
    molecule_dict_csv_path_dict = run_reaction_template_pipeline(detected_reactions_dict, cache)
    import pprint
    pprint.pprint(molecule_dict_csv_path_dict)

    