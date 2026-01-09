from rdkit import Chem
import json

# Dictionary defining various types of monomers, including their functionality types (e.g., 'vinyl', 'mono', 'di_different', 'di_identical'),
# SMARTS patterns for substructure matching, and group names for identification.
# Additional functional groups can be added based on references like "J. Chem. Inf. Model. 2023, 63, 5539−5548".
# Note: Potential for monomers with mixed groups like COCl and COOH is unaddressed.
# Can be many more functional groups added here

"""
Monomer Functionality Detection Module
--------------------------------------
This module identifies and categorizes polymer monomers based on their chemical 
functionality using SMARTS pattern matching. 

Functionality Elaboration Logic:
1. Vinyl/Mono: Requires at least one reactive site ([C]=[C] or ring structure).
2. Di_Identical: Strictly requires >= 2 occurrences of the SAME functional group 
   (e.g., two -OH groups for a diol). If only one is found, it is treated as 
   non-functional for step-growth to prevent chain termination errors.
3. Di_Different: Requires at least one occurrence of two DISTINCT functional 
   groups (e.g., one -NH2 and one -COOH for an amino acid).

This categorization ensures that only monomers capable of propagation are 
passed to the polymerization reaction engine.
"""

monomer_types = {
    "vinyl_monomer": {
        "functionality_type": "vinyl",
        "smarts_1": "[C]=[C;D1]",
        "group_name": "vinyl"
    },
    "cyclic_olefin_monomer": {
        "functionality_type": "vinyl",
        "smarts_1": "[CX3;R:1]=[CX3;R:2]",
        "group_name": "cyclic_olefin"
    },
    "lactone_monomer": {
        "functionality_type": "mono",
        "smarts_1": "[CX3;R:1](=[OX1])[OX2;R:2]",
        "group_name": "lactone"
    },
    "hydroxy_carboxylic_acid_monomer": {
        "functionality_type": "di_different",
        "smarts_1": "[OX2H1;!$(OC=*):1]",
        "smarts_2": "[CX3:2](=[O])[OX2H1]",
        "group_name": "hydroxy_carboxylic_acid"
    },
    "di_carboxylic_acid_monomer": {
        "functionality_type": "di_identical",
        "smarts_1": "[CX3:1](=[O])[OX2H1]",
        "group_name": "di_carboxylic_acid"
    },
    "di_carboxylic_acidhalide_monomer": {
        "functionality_type": "di_identical",
        "smarts_1": "[CX3:1](=[O])[Cl,Br]",
        "group_name": "di_carboxylic_acidhalide"
    },
    "diol_monomer": {
        "functionality_type": "di_identical",
        "smarts_1": "[O,S;X2;H1;!$([O,S]C=*):3]",
        "group_name": "diol"
    },
    "cyclic_anhydride_monomer": {
        "functionality_type": "mono",
        "smarts_1": "[C,c;R:1][CX3,c;R](=[OX1])[OX2,o;R][CX3,c;R](=[OX1])[C,c;R:2]",
        "group_name": "cyclic_anhydride"
    },
    "epoxide_monomer": {
        "functionality_type": "mono",
        "smarts_1": "[CX4;R:3]1[OX2;R:4][CX4;R:5]1",
        "group_name": "epoxide"
    },
    "lactam_monomer": {
        "functionality_type": "mono",
        "smarts_1": "[CX3;R:1](=[OX1])[NX3;R:2]",
        "group_name": "lactam"
    },
    "amino_acid_monomer": {
        "functionality_type": "di_different",
        "smarts_1": "[NH2;!$(NC=O)]",
        "smarts_2": "[CX3:2](=[O])[OX2H1]",
        "group_name": "amino_acid"
    },
    "di_amine_monomer": {
        "functionality_type": "di_identical",
        "smarts_1": "[N&X3;H2,H1;!$(NC=*):3]",
        "group_name": "di_amine"
    },
    "primery_di_amine_monomer": {
        "functionality_type": "di_identical",
        "smarts_1": "[C,c:6][NX3;H2;!$(N[C,S]=*)]",
        "group_name": "di_primery_amine"
    },
    "di_cyclic_anhydride_monomer": {
        "functionality_type": "di_identical",
        "smarts_1": "[CX3,c;R:1](=[OX1])[OX2,o;R][CX3,c;R:2](=[OX1])",
        "group_name": "di_cyclic_anhydride"
    },
    "di_isocyanate_monomer": {
        "functionality_type": "di_identical",
        "smarts_1": "[NX2:1]=[CX2]=[OX1,SX1:2]",
        "group_name": "di_isocyanate"
    },
    "di_epoxide_monomer": {
        "functionality_type": "di_identical",
        "smarts_1": "[CX4;H2,H1,H0;R:1]1[OX2;R:2][CX4;H1,H0;R:3]1",
        "group_name": "di_epoxide"
    }
    # need to add more functional groups here from "J. Chem. Inf. Model. 2023, 63, 5539−5548"
}
# is there monomers with both COCl and COOH groups?

def detect_monomer_functionality(smiles: str, functionality_type: str, smarts_1: str, smarts_2: str = None) -> tuple:
    """
    Detect the functionality type of a monomer based on its SMILES string using RDKit substructure matching with SMARTS patterns.
    
    Args:
        smiles (str): SMILES representation of the monomer.
        functionality_type (str): Type of functionality (e.g., 'mono', 'di_identical', 'di_different', 'vinyl').
        smarts_1 (str): Primary SMARTS pattern for matching.
        smarts_2 (str, optional): Secondary SMARTS pattern for 'di_different' types. Defaults to None.
    
    Returns:
        tuple: (functionality_count: int, count_1: int or None, count_2: int or None)
            - functionality_count: 0 (no match), 1 (single match), or 2 (dual matches).
            - count_1: Number of matches for smarts_1.
            - count_2: Number of matches for smarts_2 (if applicable).
    
    Notes:
        For 'di_identical', requires at least 2 matches of smarts_1.
        Duplicate definition present in original code; this is the first instance.
    """
    # Convert SMILES to RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0, None, None
        
    # Create pattern from primary SMARTS and check for initial match
    patt1 = Chem.MolFromSmarts(smarts_1)
    matches1 = mol.GetSubstructMatches(patt1)
    count1 = len(matches1)
    
    # has to fix di identical case
    if smarts_2:
        # Handle 'di_different' types with two distinct patterns
        patt2 = Chem.MolFromSmarts(smarts_2)
        matches2 = mol.GetSubstructMatches(patt2)
        count2 = len(matches2)
        if count1 >= 1 and count2 >= 1:
            # Count all occurrences of each pattern
            return 2, count1, count2
        else:
            return 0, None, None
    else:
        if count1 >= 1:
            # Count occurrences of primary pattern
            if functionality_type == "di_identical":
                # Check if the match occurs more than once for di_identical functionality
                if count1 >= 2:
                    return 2, count1, None
                else:
                    return 0, count1, None
            # Vinyl and mono types return 1 if at least one match is found
            return 1, count1, None
        else:
            return 0, None, None

def functional_groups_detection(monomer_dictionary: dict) -> dict:
    """
    Detect functional groups in a dictionary of monomers and categorize them based on predefined types.
    
    Args:
        monomer_dictionary (dict): Dictionary with monomer indices as keys and SMILES strings as values.
    
    Returns:
        dict: Selected monomers with detected functionalities, structured as:
            {index: {'smiles': str, functional_group_index: {details of functionality}}}
    
    Notes:
        Iterates over monomer_types and uses detect_monomer_functionality to match patterns.
        Prints detection results for each matching monomer.
        Code structure has repetitive if-blocks for each functionality_type; loops noted as needing fixes in original.
        Debugging print statement commented out.
    """
    selected_monomers = {}

    # from here code is still messed up 
    # need to fix the loops
    
    # Iterate over each monomer in the input dictionary
    for indexm, smiles in monomer_dictionary.items():
        functional_group_index = 0
        # Iterate over each predefined functional group type
        for functional_group in monomer_types.values():
            f_type = functional_group["functionality_type"]
            s1 = functional_group["smarts_1"]
            s2 = functional_group.get("smarts_2") # Get smarts_2 if it exists
            
            # Use the core detection logic for all types
            functionality, count1, count2 = detect_monomer_functionality(smiles, f_type, s1, s2)
            
            # Check for matches based on functionality requirements
            is_match = False
            if f_type in ["vinyl", "mono"] and functionality == 1:
                is_match = True
            elif f_type in ["di_different", "di_identical"] and functionality == 2:
                is_match = True
                
            if is_match:
                functional_group_index += 1
                print(f"Monomer {indexm} ({smiles}) has functionality: {functional_group['group_name']}")
                
                if indexm not in selected_monomers:
                    selected_monomers[indexm] = {"smiles": smiles}
                
                # Construct result entry
                res_entry = {
                    "functionality_type": f_type,
                    "functional_group_name": functional_group["group_name"],
                    "functional_group_smarts_1": s1,
                    "functional_count_1": count1
                }
                
                if s2:
                    res_entry["functional_group_smarts_2"] = s2
                    res_entry["functional_count_2"] = count2
                    
                selected_monomers[indexm][functional_group_index] = res_entry

    # for debugging purposes
    print(json.dumps(selected_monomers, indent=4))
    return selected_monomers


if __name__ == "__main__":
    monomer_dictionary = {
    1: 'C=CCCCO', 
    2: 'CCC',
    3: 'C=Cc1cccc(C=C)c1',
    4: 'C=Cc1cc(C(=O)O)cc(C(=O)O)c1',
    5: 'Nc1cc(C(=O)O)cc(C(=O)O)c1',}

    # Example monomer dictionary for testing
    monomer_dictionary = {
        1:  "C=CCO",                 # vinyl monomer # veryfied
        2:  "CCCCCCC=C",           
        3:  " CC=C1CC2CC1C=C2",           # cyclic olefin monomer (norbornene) # veryfied
        4: "C1=CC2CCC1C2",           # cyclic olefin monomer (norbornadiene) # veryfied
        5:  "C1CCC(=O)OCC1",           # lactone monomer (ε-caprolactone) # veryfied
        6:  "OCCC(O)=O",           # hydroxy carboxylic acid (lactic acid) # veryfied
        7:  "C1=CC=C(C(=C1)C(=O)O)O", # hydroxy carboxylic acid (p-hydroxybenzoic acid) # veryfied
        8:  "O=C(O)CCCCC(=O)O",      # di-carboxylic acid (adipic acid) # veryfied
        9:  "O=C(Cl)c1ccc(cc1)C(=O)Cl", # di-carboxylic acid halide (terephthaloyl chloride) # veryfied
        10: "OCCO",                  # diol (ethylene glycol) # veryfied
        11: "O=C1OC(=O)C=C1",        # cyclic anhydride (maleic anhydride) # veryfied
        12: "C1CO1",                 # epoxide monomer (ethylene oxide) # veryfied
        13: "O=C1CCCCN1",            # lactam monomer (caprolactam) # veryfied
        14: "NCC(=O)O",              # amino acid (glycine) # veryfied
        15: "NCCCC(=O)O",            # amino acid (aminocaproic acid) # veryfied
        16: "NCCN",                  # di-amine (ethylenediamine)
        17: "NCCCCCCN",              # primary di-amine (hexamethylenediamine)
        18: "O=C1OC(=O)c2cc(C3OC(=O)C=CC3=O)ccc2C1=O",  # di-cyclic anhydride (pyromellitic dianhydride)
        19: "O=C=N-C-N=C=O",                # di-isocyanate (simplified; represents TDI core)
        20: "C1=CC(=CC=C1C(C2=CC=CC=C2)(O)CC3CO3)CC4CO4" # di-epoxide (bisphenol A diglycidyl ether)
    }
    monomer_dictionary = {1: "OC(F)c1cc(C(O)Cl)cc(C(O)Br)c1",
                          2: "O=C(O)c1c(F)c(C(=O)O)c(Br)c(C(=O)O)c1Cl"}
    # Run functional groups detection
    print(json.dumps(functional_groups_detection(monomer_dictionary), indent=4))

        