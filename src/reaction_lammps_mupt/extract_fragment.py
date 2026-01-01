import os
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops
from rdkit.Chem.AllChem import MolFragmentToSmiles
from pathlib import Path


def create_the_reaction(reaction_smarts, reactant_smiles_1, reactant_smiles_2=None):
    # Create reactant molecules from SMILES strings
    reaction = AllChem.ReactionFromSmarts(reaction_smarts)
    if reaction is None:
        reactant_smiles_2 = reactant_smiles_1  # Handle case with only one reactant
    # 3. Create reactant molecules
    reactant1 = Chem.MolFromSmiles(reactant_smiles_1)
    reactant2 = Chem.MolFromSmiles(reactant_smiles_2) 
    reactants = (f"{reactant_smiles_1}.{reactant_smiles_2}")
    # 4. Run the reaction
    
    products = reaction.RunReactants((reactant1, reactant2))
    if not products:
        print("No products were formed from the reaction.")
    # 5. Process and display results
    products_list = []
    for result_set in products:
        product_string = ""
        for j,product in enumerate(result_set):
            product_string += (Chem.MolToSmiles(product))
            if j < len(result_set) - 1:
                product_string += "."
        products_list.append(product_string)
    # for debugging purposes
    # print("Products:", products_list)
    return reactants, products_list
    

def template_fragment_extractor(
        reaction_smarts = "[OX2H1;!$(OC=*):1].[CX3:2](=[O])[OX2H1]>>[OX2:1]-[CX3:2](=[O]).O",
        reactants ="C1=CC=C(C(=C1)C(=O)O)O.OCCC(O)=O", 
        products_list = 'O=C(CCO)Oc1ccccc1C(=O)O.O'):
    monomers = Chem.MolFromSmiles(reactants)
    products_list = Chem.MolFromSmiles(products_list)
    img =Draw.MolToImage(monomers)
    monomers = Chem.AddHs(monomers)
    products_list = Chem.AddHs(products_list)
    img_2 =Draw.MolToImage(products_list)
    img_2.save(f"highlighted_molecule_products_45454.png")
    substructure = Chem.MolFromSmarts("[OX2H1;!$(OC=*):1].[CX3:2](=[O])[OX2H1]")
    substructure = Chem.MolFromSmarts("[CX3:2](=[O])[OX2H1]")
    matches = monomers.GetSubstructMatches(substructure)
    print("Matches found:", matches)
    atom_indices = [index for match in matches for index in match]
    print("Atom indices to highlight:", atom_indices)
    img = Draw.MolToImage(monomers, highlightAtoms=atom_indices)
    img.save(f"testing_highlighted.png")
    monomers = Chem.MolFromSmiles(reactants)
    print(monomers)
    pass




if __name__ == "__main__":
    template_fragment_extractor()
    # create_the_reaction(
    #     "[OX2H1;!$(OC=*):1].[CX3:2](=[O])[OX2H1]>>[OX2:1]-[CX3:2](=[O]).O",
    #     "C1=CC=C(C(=C1)C(=O)O)O",
    #     "OCCC(O)=O"
    # )
    selected_reactions = {
  32: {
    "reaction_name": "Hydroxy Carboxylic and Hydroxy Carboxylic Polycondensation(Polyesterification)",
    "same_reactants": False,
    "reactant_1": "hydroxy_carboxylic_acid",
    "reactant_2": "hydroxy_carboxylic_acid",
    "product": "polyester_chain",
    "delete_atom": True,
    "reaction": "[OX2H1;!$(OC=*):1].[CX3:2](=[O])[OX2H1]>>[OX2:1]-[CX3:2](=[O]).O",
    "monomer_1": {
      "smiles": "C1=CC=C(C(=C1)C(=O)O)O",
      "functionality_type": "di_different",
      "functional_group_name": "hydroxy_carboxylic_acid",
      "functional_hroup_smarts_1": "[OX2H1;!$(OC=*):1]",
      "functional_count_1": 1,
      "functional_hroup_smarts_2": "[CX3:2](=[O])[OX2H1]",
      "functional_count_2": 1
    },
    "monomer_2": {
      "smiles": "OCCC(O)=O",
      "functionality_type": "di_different",
      "functional_group_name": "hydroxy_carboxylic_acid",
      "functional_hroup_smarts_1": "[OX2H1;!$(OC=*):1]",
      "functional_count_1": 1,
      "functional_hroup_smarts_2": "[CX3:2](=[O])[OX2H1]",
      "functional_count_2": 1
    }
  }
}
