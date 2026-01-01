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
    sys.exit(message) # change this to not to close while MuPT, But our packafe will close

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

# MAIN.PY name will changeed in the future to another name 