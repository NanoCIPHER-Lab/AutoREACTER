from os import path
import rdkit
import pathlib
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import Draw
path = pathlib.Path(__file__).parent.resolve()
from reaction_template_pipeline.map_reactant_atoms import map_reactant_atoms


def reaction_selector(selected_reactions_dict):
    """Still there is no method to verify the all possible reactions 
    can be implemented correctly. So, for now, we will just process with one or two reactions.
    TODO: Add verification method for all reactions."""
    for key, reaction in selected_reactions_dict.items():
        reaction_smarts = reaction["reaction"]
        reactant_smiles1 = reaction["monomer_1"]["smiles"]
        same_reactants = reaction["same_reactants"]
        delete_atom = reaction["delete_atom"]
        reactant_smiles2 = reaction["monomer_2"]["smiles"] if same_reactants is False else reactant_smiles1
        print(f"Processing Reaction ID: {key}")
        print(f"Reaction Name: {reaction_smarts}")
        print(f"Same Reactants: {reaction['same_reactants']}")
        print(f"Reactant 1: {reaction['reactant_1']}")
        if reaction["same_reactants"] is False:
            print(f"Reactant 2: {reaction['reactant_2']}")
        print(f"Product: {reaction['product']}")
        print(f"Delete Atom: {reaction['delete_atom']}")
        print(f"Reaction SMARTS: {reaction['reaction']}")
        print("Monomer 1 Details:")
        for monomer_key, monomer_value in reaction['monomer_1'].items():
            print(f"  {monomer_key}: {monomer_value}")
        if reaction["same_reactants"] is False:
            print("Monomer 2 Details:")
            for monomer_key, monomer_value in reaction["monomer_2"].items():
                print(f"  {monomer_key}: {monomer_value}")
        print("-" * 40)
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        if rxn is None:
            raise ValueError(f"Invalid reaction SMARTS: {reaction_smarts}")
        reactant1 = Chem.MolFromSmiles(reactant_smiles1)
        reactant2 = Chem.MolFromSmiles(reactant_smiles2) 
        if reactant1 is None:
            raise ValueError(f"Invalid SMILES string: {reactant_smiles1}")
        if reactant2 is None:
            raise ValueError(f"Invalid SMILES string: {reactant_smiles2}")
        if reactant1 is None:
            raise ValueError(f"Invalid SMILES string for reactant 1: {reactant_smiles1}")
        if reactant2 is None:
            raise ValueError(f"Invalid SMILES string for reactant 2: {reactant_smiles2}")
        reactant1 = Chem.AddHs(reactant1)
        reactant2 = Chem.AddHs(reactant2)
        combined_reactants, combined_products, byproduct_map_numbers = map_reactant_atoms(reactant1, reactant2, rxn , delete_atom)
        ## for debugging purpose, save the image
        Draw.MolsToGridImage([combined_reactants, combined_products], molsPerRow=2, subImgSize=(1800, 1800)).save(path / f"reaction_{key}.png")
        print("Delete Atom Map Numbers:", byproduct_map_numbers)

        if not same_reactants:
            print("Processing second reaction case:")
            combined_reactants, combined_products, byproduct_map_numbers = map_reactant_atoms(reactant2, reactant1, rxn, delete_atom)
            ## for debugging purpose, save the image
            Draw.MolsToGridImage([combined_reactants, combined_products], molsPerRow=2, subImgSize=(1800, 1800)).save(path / f"reaction_{key}_same_reaction.png")
            print("Delete Atom Map Numbers (swapped reactants):", byproduct_map_numbers)
        
        print("=" * 80)     
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
    selected_reactions_dict = {
        31: {
            "reaction_name": "Hydroxy Carboxylic and Hydroxy Carboxylic Polycondensation(Polyesterification)",
            "same_reactants": False,
            "reactant_1": "hydroxy_carboxylic_acid",
            "reactant_2": "hydroxy_carboxylic_acid",
            "product": "polyester_chain",
            "delete_atom": True,
            "reaction": "[O;!$(OC=*):1]-[H:3].[CX3:2](=[O:5])[OX2H1:4]>>[OX2:1]-[CX3:2](=[O:5]).[O:4]-[H:3]",
            "monomer_1": {
                "smiles": "C1=CC=C(C(=C1)C(=O)O)O",
                "functionality_type": "di_different",
                "functional_group_name": "hydroxy_carboxylic_acid",
                "functional_group_smarts_1": "[OX2H1;!$(OC=*):1]",
                "functional_count_1": 1,
                "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
                "functional_count_2": 1
            },
            "monomer_2": {
                "smiles": "OCCC(O)=O",
                "functionality_type": "di_different",
                "functional_group_name": "hydroxy_carboxylic_acid",
                "functional_group_smarts_1": "[OX2H1;!$(OC=*):1]",
                "functional_count_1": 1,
                "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
                "functional_count_2": 1
            }
        }
    }
    reaction_selector(selected_reactions_dict)
