import importlib
import sys
from initializaion import initialize
from input_parser import InputParser





if __name__ == "__main__":
    inputs =  inputs = {
    "simulation_name": "MySimulation",
    "temperature": [300, 400, 500],
    "density": 0.8,
    "monomers": {
        1:  "C=CCO",                           # vinyl monomer
        2:  "CCCCCCC=C",
        3:  "CC=C1CC2CC1C=C2",                 # cyclic olefin (norbornene)
        4:  "C1=CC2CCC1C2",                    # cyclic olefin (norbornadiene)
        5:  "C1CCC(=O)OCC1",                   # lactone (ε-caprolactone)
        6:  "OCCC(O)=O",                       # hydroxy carboxylic acid (lactic acid)
        7:  "C1=CC=C(C(=C1)C(=O)O)O",          # hydroxy carboxylic acid (p-hydroxybenzoic acid)
        8:  "O=C(O)CCCCC(=O)O",                # di-carboxylic acid (adipic acid)
        9:  "O=C(Cl)c1ccc(cc1)C(=O)Cl",        # di-acid chloride (terephthaloyl chloride)
        10: "OCCO",                            # diol (ethylene glycol)
        11: "O=C1OC(=O)C=C1",                  # cyclic anhydride (maleic anhydride)
        12: "C1CO1",                           # epoxide (ethylene oxide)
        13: "O=C1CCCCN1",                      # lactam (caprolactam)
        14: "NCC(=O)O",                        # amino acid (glycine)
        15: "NCCCC(=O)O",                      # amino acid (aminocaproic acid)
        16: "NCCN",                            # diamine (ethylenediamine)
        17: "NCCCCCCN",                        # diamine (hexamethylenediamine)
        18: "O=C1OC(=O)c2cc(C3OC(=O)C=CC3=O)ccc2C1=O",  # dianhydride (pyromellitic dianhydride)
        19: "O=C=N-C-N=C=O",                   # di-isocyanate (generic)
        20: "C1=CC(=CC=C1C(C2=CC=CC=C2)(O)CC3CO3)CC4CO4" # diepoxide (BADGE)
        },
    "number_of_monomers": {
        1:  500,   # simple vinyl
        2:  400,
        3:  300,   # norbornene
        4:  300,   # norbornadiene
        5:  250,   # lactone
        6:  300,   # lactic acid
        7:  200,   # aromatic hydroxy acid
        8:  300,   # adipic acid
        9:  150,   # acid chloride (reactive, keep lower)
        10: 400,   # ethylene glycol
        11: 200,   # maleic anhydride
        12: 350,   # epoxide (small)
        13: 250,   # caprolactam
        14: 300,   # glycine
        15: 200,   # aminocaproic acid
        16: 300,   # ethylenediamine
        17: 250,   # hexamethylenediamine
        18: 100,   # pyromellitic dianhydride (bulky)
        19: 150,   # di-isocyanate
        20: 150,   # BADGE diepoxide
        },
    }

    initialize()
    inputs = InputParser(inputs)
    print(inputs.validated_inputs)


# first Check RDKit is imported correctly
# Then Check LUNAR is imported correctly
# There parts a will very strict - will close the program if not imported correctly

# run User Input Parser
# will have class InputParser

# then will go to the reaction detector part 
# if we failed to generate reactions we say "No reactions detected with given monomers"
# we generate scripts without automatic reactions we request for templates from user and the map file -> we say this to lunar and The INPUT Script Generator

# then will come back to main.py file 

# Next part will be to go to the LUNAR for parameterization

# then will come back to main.py file 

# Then with all the data will go to INPUT Script Generator

# and finally will result in OUTPUT files as a ready Folder

# MAIN.PY name will changed in the future to another name 