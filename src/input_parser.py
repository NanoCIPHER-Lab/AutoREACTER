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


from rdkit import Chem
    
    
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
            
    def validate_smiles_rdkit(self, inputs: dict) -> None:
        """
        Validate all SMILES strings in the provided inputs dictionary using RDKit.

        Args:
            inputs (dict): A dictionary containing monomer SMILES strings under the
                "monomers" key, mapping monomer indices to SMILES strings.

        Returns:
            None

        Raises:
            ValueError: If any SMILES string in the inputs is invalid.
        """
        def validate_smiles_rdkit(smiles_monomer: str) -> bool:
                mol = Chem.MolFromSmiles(smiles_monomer)
                if mol is None:
                    raise ValueError("Invalid SMILES string. Terminating program.")
                else:
                    return True

        for monomer_number in inputs["monomers"]:
            smiles_string = inputs["monomers"][monomer_number]
            validate_smiles_rdkit(smiles_string)
        return
    
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

def input(inputs):
    parser = InputParser()
    parser.validate_smiles_rdkit(inputs)
    return inputs

if __name__ == "__main__":
    user_inputs = input(inputs)
    print(user_inputs)
"""
TODO: Proceed with the rest of the program using user_inputs
    in here we have to process the inputs dictionary
    e.g.,
"Number of monomers": {
        1: 1000,
        2: 1000,
    }, # or
    "stoichiometric_ratio": {
        1: 1,
        2: 1,
    },
    as for given inputs we have to calculate the number of monomers based on total atoms and stoichiometric ratio
    also in here we will be able to calculate the box size and inital density based on the number of monomers and their molar masses can be obtained from rdkit
    at the end we will have to give properly formatted inputs to the main.py file
"""