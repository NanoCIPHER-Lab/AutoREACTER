import sys
import os
try:
    import rdkit
except ModuleNotFoundError:
    message = (
            "ERROR: RDKit is not installed in this environment.\n\n"
            "👉 To install RDKit, run one of the following:\n"
            "   conda install -c conda-forge rdkit   (recommended)\n"
            "   pip install rdkit-pypi              (alternative)\n\n"
            "Exiting program..."
        )
    sys.exit(message)
# all of above will go to the main.py file
# Purpose of this file is to parse inputs from user
# and feed to main.py file 
# when inputting monomers make sure no to duplicate smiles strings
monomers = {1: 'CCO', 2: 'CCC'}


# relative Import 
# this will go to the main.py file and import the lunar handler class
lunar = LUNAR.LunarHandler(monomer_dictionary=monomers)
    
    
class InputParser:
    def __init__(self):
        pass
    def validate_smiles_rdkit(self, monomer_number: int) -> str:
        """
        Prompt user for a SMILES string and validate with RDKit.

        Args:
            monomer_number (int): The index of the current monomer being requested.

        Returns:
            str: A validated SMILES string.
        """

        while True:
            smiles_string = input(
                f"Please enter the SMILES string for monomer {monomer_number}: "
            )
            mol = Chem.MolFromSmiles(smiles_string)
            if mol is None:
                print("Invalid SMILES string. Please try again.")
            else:
                return smiles_string
            

# Each every data will be stored in one dictionary
# will have series of Temparatures 
# Will have series of Number of monomers/ atoms (most probably will)
# Density
# Simulation Name
# May be will have a GUI
# Example
inputs = {
    "simulation_name": "MySimulation",
    "temperature": [300, 400, 500],

    "density": 0.8,
    "monomers": {
        1: "CCO",
        2: "CCC",
    },
    "Number of monomers": {
        1: 1000,
        2: 1000,
    }, # or
    "stoichiometric_ratio": {
        1: 1,
        2: 1,
    },
    "Number of total atoms": [10000, 100000, 1000000],
}


def input(inputs=inputs):
    return inputs