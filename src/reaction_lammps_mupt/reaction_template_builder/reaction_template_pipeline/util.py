from rdkit import Chem
from rdkit.Chem import AllChem


def get_all_possible_reactions(detected_reactions):
    """Get all possible reactions from detected reactions.
    
    Args:
        detected_reactions (dict): Detected reactions dictionary.
    Returns:
        list: List of all possible reactions.
    """
    for key in detected_reactions:
        reaction_dict = detected_reactions[key]
        reactant_1 = reaction_dict["monomer_1"]["smiles"]
        same_reactants = reaction_dict["same_reactants"]
        if same_reactants:
            reactant_2 = reactant_1
        else:
            reactant_2 = reaction_dict["monomer_2"]["smiles"]
        rxn_smarts = reaction_dict["reaction"]
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
        reactant_mol_1 = Chem.MolFromSmiles(reactant_1)
        reactant_mol_1 = AllChem.AddHs(reactant_mol_1)
        reactant_mol_2 = Chem.MolFromSmiles(reactant_2)
        reactant_mol_2 = AllChem.AddHs(reactant_mol_2)
        products = rxn.RunReactants((reactant_mol_1, reactant_mol_2))
        if not products:
            raise ValueError("No products generated from the reaction.")
        number_of_products = len(products)
        if not same_reactants:
            products = rxn.RunReactants((reactant_mol_2, reactant_mol_1))
            number_of_products += len(products)
        reaction_dict["number_of_products"] = number_of_products
    return detected_reactions
 
    

if __name__ == "__main__":
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
    Non_monomer_molecules_to_retain = ["CCO"]
    import json
    print(json.dumps(get_all_possible_reactions(detected_reactions), indent=2))
