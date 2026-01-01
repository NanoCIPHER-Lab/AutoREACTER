import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions

def reaction_selector(selected_reactions_dict):
    for key, reaction in selected_reactions_dict.items():
        print(f"Processing Reaction ID: {key}")
        print(f"Reaction Name: {reaction['reaction_name']}")
        print(f"Same Reactants: {reaction['same_reactants']}")
        print(f"Reactant 1: {reaction['reactant_1']}")
        if not reaction['same_reactants']:
            print(f"Reactant 2: {reaction['reactant_2']}")
        print(f"Product: {reaction['product']}")
        print(f"Delete Atom: {reaction['delete_atom']}")
        print(f"Reaction SMARTS: {reaction['reaction']}")
        print("Monomer 1 Details:")
        for monomer_key, monomer_value in reaction['monomer_1'].items():
            print(f"  {monomer_key}: {monomer_value}")
        if not reaction['same_reactants']:
            print("Monomer 2 Details:")
            for monomer_key, monomer_value in reaction['monomer_2'].items():
                print(f"  {monomer_key}: {monomer_value}")
        print("-" * 40)        
        break  # Remove this break to process all reactions

if __name__ == "__main__":
     # Example monomer dictionary
    selected_reactions_dict = {
        30: {
            "reaction_name": "Hydroxy Carboxylic Acid Polycondensation(Polyesterification)",
            "same_reactants": True,
            "reactant_1": "hydroxy_carboxylic_acid",
            "product": "polyester_chain",
            "delete_atom": True,
            "reaction": "[OX2H1;!$(OC=*):1].[CX3:2](=[O])[OX2H1]>>[OX2:1]-[CX3:2](=[O]).O",
            "monomer_1": {
                "smiles": "C1=CC=C(C(=C1)C(=O)O)O",
                "functionality_type": "di_different",
                "functional_group_name": "hydroxy_carboxylic_acid",
                "functional_group_smarts_1": "[OX2H1;!$(OC=*):1]",
                "functional_count_1": 1,
                "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
                "functional_count_2": 1
            }
        },
        31: {
            "reaction_name": "Hydroxy Carboxylic and Hydroxy Carboxylic Polycondensation(Polyesterification)",
            "same_reactants": False,
            "reactant_1": "hydroxy_carboxylic_acid",
            "reactant_2": "hydroxy_carboxylic_acid",
            "product": "polyester_chain",
            "delete_atom": True,
            "reaction": "[OX2H1;!$(OC=*):1].[CX3:2](=[O])[OX2H1]>>[OX2:1]-[CX3:2](=[O]).O",
            "monomer_1": {
                "smiles": "OCCC(O)=O",
                "functionality_type": "di_different",
                "functional_group_name": "hydroxy_carboxylic_acid",
                "functional_group_smarts_1": "[OX2H1;!$(OC=*):1]",
                "functional_count_1": 1,
                "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
                "functional_count_2": 1
            },
            "monomer_2": {
                "smiles": "C1=CC=C(C(=C1)C(=O)O)O",
                "functionality_type": "di_different",
                "functional_group_name": "hydroxy_carboxylic_acid",
                "functional_group_smarts_1": "[OX2H1;!$(OC=*):1]",
                "functional_count_1": 1,
                "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
                "functional_count_2": 1
            }
        },
        33: {
            "reaction_name": "Diol and Di-Carboxylic Acid Polycondensation(Polyesterification)",
            "same_reactants": False,
            "reactant_1": "diol",
            "reactant_2": "di_carboxylic_acid",
            "product": "polyester_chain",
            "delete_atom": True,
            "reaction": "[CX3:2](=[O])[OX2H1,Cl,Br].[O,S;X2;H1;!$([O,S]C=*):3]>>[CX3:2](=[O])-[O,S;X2;!$([O,S]C=*):3]",
            "monomer_1": {
                "smiles": "OCCO",
                "functionality_type": "di_identical",
                "functional_group_name": "diol",
                "functional_group_smarts_1": "[O,S;X2;H1;!$([O,S]C=*):3]",
                "functional_count_1": 2
            },
            "monomer_2": {
                "smiles": "O=C(O)CCCCC(=O)O",
                "functionality_type": "di_identical",
                "functional_group_name": "di_carboxylic_acid",
                "functional_group_smarts_1": "[CX3:1](=[O])[OX2H1]",
                "functional_count_1": 2
            }
        },
        34: {
            "reaction_name": "Amino Acid Polycondensation (Polyamidation)",
            "same_reactants": True,
            "reactant_1": "amino_acid",
            "product": "polyamide_chain",
            "delete_atom": True,
            "reaction": "[NX3;H2,H1;!$(OC=*):1].[CX3:2](=[O])[OX2H1]>>[NX3:1]-[CX3:2](=[O]).O",
            "monomer_1": {
                "smiles": "NCC(=O)O",
                "functionality_type": "di_different",
                "functional_group_name": "amino_acid",
                "functional_group_smarts_1": "[NH2;!$(NC=O)]",
                "functional_count_1": 1,
                "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
                "functional_count_2": 1
            }
        }
    }
    reaction_selector(selected_reactions_dict)
