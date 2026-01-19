from os import path
import rdkit
import pathlib
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import Draw
path = pathlib.Path(__file__).parent.resolve()
from reaction_template_pipeline.map_reactant_atoms import map_reactant_atoms, map_product_atoms
from reaction_template_pipeline.walker import reaction_atom_walker
if __name__ == "__main__":
    detected_reactions = {
    1: {
        "reaction_name": "Hydroxy Carboxylic Acid Polycondensation (Polyesterification)",
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
        "smiles": "C1=CC=C(C(=C1)C(=O)O)O",
        "functionality_type": "di_different",
        "functional_group_name": "hydroxy_carboxylic_acid",
        "functional_group_smarts_1": "[OX2H1;!$(OC=*):1]",
        "functional_count_1": 1,
        "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
        "functional_count_2": 1,
        },
    },

    2: {
        "reaction_name": "Hydroxy Carboxylic Acid Polycondensation (Polyesterification)",
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
        "smiles": "OCCC(O)=O",
        "functionality_type": "di_different",
        "functional_group_name": "hydroxy_carboxylic_acid",
        "functional_group_smarts_1": "[OX2H1;!$(OC=*):1]",
        "functional_count_1": 1,
        "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
        "functional_count_2": 1,
        },
    },

    3: {
        "reaction_name": "Hydroxy Carboxylic Acid + Hydroxy Carboxylic Acid Polycondensation (Polyesterification)",
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
        "smiles": "C1=CC=C(C(=C1)C(=O)O)O",
        "functionality_type": "di_different",
        "functional_group_name": "hydroxy_carboxylic_acid",
        "functional_group_smarts_1": "[OX2H1;!$(OC=*):1]",
        "functional_count_1": 1,
        "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
        "functional_count_2": 1,
        },
        "monomer_2": {
        "smiles": "OCCC(O)=O",
        "functionality_type": "di_different",
        "functional_group_name": "hydroxy_carboxylic_acid",
        "functional_group_smarts_1": "[OX2H1;!$(OC=*):1]",
        "functional_count_1": 1,
        "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
        "functional_count_2": 1,
        },
    },
    }
    {
  1: {
    "reaction_name": "Hydroxy Carboxylic Acid Polycondensation (Polyesterification)",
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
      "smiles": "C1=CC=C(C(=C1)C(=O)O)O",
      "functionality_type": "di_different",
      "functional_group_name": "hydroxy_carboxylic_acid",
      "functional_group_smarts_1": "[OX2H1;!$(OC=*):1]",
      "functional_count_1": 1,
      "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
      "functional_count_2": 1,
    },
  },

  2: {
    "reaction_name": "Hydroxy Carboxylic Acid Polycondensation (Polyesterification)",
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
      "smiles": "OCCC(O)=O",
      "functionality_type": "di_different",
      "functional_group_name": "hydroxy_carboxylic_acid",
      "functional_group_smarts_1": "[OX2H1;!$(OC=*):1]",
      "functional_count_1": 1,
      "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
      "functional_count_2": 1,
    },
  },

  3: {
    "reaction_name": "Hydroxy Carboxylic Acid + Hydroxy Carboxylic Acid Polycondensation (Polyesterification)",
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
      "smiles": "C1=CC=C(C(=C1)C(=O)O)O",
      "functionality_type": "di_different",
      "functional_group_name": "hydroxy_carboxylic_acid",
      "functional_group_smarts_1": "[OX2H1;!$(OC=*):1]",
      "functional_count_1": 1,
      "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
      "functional_count_2": 1,
    },
    "monomer_2": {
      "smiles": "OCCC(O)=O",
      "functionality_type": "di_different",
      "functional_group_name": "hydroxy_carboxylic_acid",
      "functional_group_smarts_1": "[OX2H1;!$(OC=*):1]",
      "functional_count_1": 1,
      "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
      "functional_count_2": 1,
    },
  },
}
Non_monomer_molecules_to_retain = ["CCO"]
