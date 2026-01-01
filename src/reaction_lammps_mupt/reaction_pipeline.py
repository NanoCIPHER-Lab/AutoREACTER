import os
import sys
import json

try:
    # Case 1 — Running as part of a package
    from .functional_groups_detector import functional_groups_detection
    from .reaction_detector import reaction_selector

except ImportError:
    # Case 2 — Running standalone: fix sys.path automatically
    CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
    PARENT_DIR = os.path.dirname(CURRENT_DIR)

    if PARENT_DIR not in sys.path:
        sys.path.insert(0, PARENT_DIR)

    from REACTION_DETECTOR.functional_groups_detector import functional_groups_detection
    from REACTION_DETECTOR.reaction_detector import reaction_selector




def reaction_pipeline(monomer_dictionary):
    """
    Main pipeline function to detect and select polymerization reactions for given monomers.

    Args:
        monomer_dictionary (dict): A dictionary where keys are monomer names (e.g., integers) and
                                  values are their SMILES strings representing chemical structures.

    Returns:
        dict: A dictionary containing selected polymerization reactions for each monomer,
              keyed by monomer name.
    """
    # Step 1: Detect functional groups in the monomers using the detection module
    monomer_with_functional_groups = functional_groups_detection(monomer_dictionary)

    # Step 2: Select appropriate polymerization reactions based on detected functional groups
    selected_reactions = reaction_selector(monomer_with_functional_groups)

    # Placeholder: Add reaction cutter functionality here in future iterations
    # (e.g., to refine or segment selected reactions)

    return selected_reactions


if __name__ == "__main__":
    # Example monomer dictionary for testing the pipeline
    # Keys are integer identifiers; values are SMILES strings for various monomers
    monomer_dictionary = {
        1:  "C=CCO",                 # Allyl alcohol (vinyl alcohol-like monomer) - verified
        2:  "CCCCCCC=C",             # 1-Octene (linear olefin monomer)
        3:  "CC=C1CC2CC1C=C2",       # Norbornene (cyclic olefin monomer) - verified
        4:  "C1=CC2CCC1C2",          # Bicyclo[2.2.1]hept-2-ene derivative (cyclic olefin, e.g., norbornadiene-like) - verified
        5:  "C1CCC(=O)OCC1",         # ε-Caprolactone (lactone monomer) - verified
        6:  "OCCC(O)=O",             # Lactic acid (hydroxy carboxylic acid monomer) - verified
        7:  "C1=CC=C(C(=C1)C(=O)O)O", # p-Hydroxybenzoic acid (hydroxy carboxylic acid monomer) - verified
        8:  "O=C(O)CCCCC(=O)O",      # Adipic acid (di-carboxylic acid monomer) - verified
        9:  "O=C(Cl)c1ccc(cc1)C(=O)Cl", # Terephthaloyl chloride (di-carboxylic acid halide monomer) - verified
        10: "OCCO",                  # Ethylene glycol (diol monomer) - verified
        11: "O=C1OC(=O)C=C1",        # Maleic anhydride (cyclic anhydride monomer) - verified
        12: "C1CO1",                 # Ethylene oxide (epoxide monomer) - verified
        13: "O=C1CCCCN1",            # ε-Caprolactam (lactam monomer) - verified
        14: "NCC(=O)O",              # Glycine (amino acid monomer) - verified
        15: "NCCCC(=O)O",            # 6-Aminocaproic acid (amino acid monomer) - verified
        16: "NCCN",                  # Ethylenediamine (di-amine monomer)
        17: "NCCCCCCN",              # Hexamethylenediamine (primary di-amine monomer)
        18: "O=C1OC(=O)c2cc(C3OC(=O)C=CC3=O)ccc2C1=O",  # Pyromellitic dianhydride (di-cyclic anhydride monomer)
        19: "O=C=NC6=C(N=C=O)C=CC=C6", # Toluene diisocyanate (di-isocyanate monomer; corrected SMILES for TDI core)
        20: "C1=CC(=CC=C1C(C2=CC=CC=C2)(O)CC3CO3)CC4CO4" # Bisphenol A diglycidyl ether (di-epoxide monomer)
    }
    # Run the pipeline on the example monomers
    selected_reactions = reaction_pipeline(monomer_dictionary)
    # Print the results in a formatted JSON structure for readability
    print("Selected Reactions:")
    print(json.dumps(selected_reactions, indent=2))
# End of reaction_pipeline.py