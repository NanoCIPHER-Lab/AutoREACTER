"""
* Monomer Functionality Detection Module
--------------------------------------
Purpose
This module takes a set of monomer SMILES strings and determines which monomers
are chemically “eligible” to participate in polymerization reactions. It does
this by scanning each monomer for known reactive functional groups using RDKit
SMARTS substructure matching. The output is a structured dictionary describing
which functional group(s) were found in each monomer and how many reactive sites
of each type exist.
* What “functionality” means here
- A “functional group” is a chemical motif that can react (e.g., -OH, -COOH,
  epoxide rings, vinyl double bonds, etc.).
- A “reactive site count” is the number of times that motif appears in the
  monomer. For example, ethylene glycol has two -OH sites; adipic acid has two
  -COOH sites.
- A “functionality type” is a rule for deciding whether the monomer is allowed
  to proceed to the polymerization engine. This module is intentionally a
  FILTER + ANNOTATOR:
  (1) Filter: remove monomers that would cause chain termination or invalid
      step-growth logic.
  (2) Annotator: retain counts so later steps can handle multi-site reactions.
* Core approach
1) Convert SMILES → RDKit Mol
2) For each defined functional group class (stored in monomer_types):
   - Build RDKit pattern(s) from SMARTS
   - Find ALL matches using substructure matching
   - Count how many matches exist
   - Decide whether the monomer “qualifies” under the rules below
   - If qualified, record:
       • functionality_type
       • functional_group_name
       • SMARTS pattern(s)
       • site count(s)
* Functionality Elaboration Logic (rules enforced)
1) Vinyl and Mono (presence-based eligibility)
   - Goal: identify monomers that can undergo chain-growth / ring-opening style
     Reactions where a single reactive motif is sufficient to participate.
   - Pass condition:
       • at least one match to the SMARTS pattern (count >= 1)
   - Notes:
       • If a monomer has multiple vinyl sites (e.g., divinyl crosslinker),
         it still “passes” here because it is reactive. The multiplicity is
         tracked in the site count and is intended to be exploited later by
         the reaction engine (crosslinking / branching logic).
       • This stage does NOT decide how many reactions will happen; it only
         verifies that the reaction is possible and records the capacity.
* 2) Di_Identical (step-growth safety rule: prevent termination)
   - Goal: enforce correct step-growth behavior for monomers expected to link
     through identical end-groups (A-A type monomers such as diols, diacids,
     diamines, diisocyanates, diepoxides).
   - Pass condition:
       • at least two occurrences of the SAME functional group (count >= 2)
   - Fail behavior:
       • if only one match is found (count == 1), the monomer is treated as
         non-functional for this category (i.e., rejected for step-growth).
   - Why this matters:
       • A monomer with only one reactive end group behaves like a chain
         stopper in step-growth polymerization, which can silently break
         propagation and produce wrong molecular weight distributions.
       • This module intentionally blocks those “A1” species from entering
         The step-growth reaction engine under a “di_identical” assumption.
* 3) Di_Different (heterobifunctional eligibility)
   - Goal: identify monomers that carry two distinct reactive groups (A-B type)
     such as amino acids (NH2 + COOH), hydroxy acids (OH + COOH), etc.
   - Pass condition:
       • at least one match to SMARTS_1 (count1 >= 1)
       • AND at least one match to SMARTS_2 (count2 >= 1)
   - Notes:
       • The module records BOTH counts. This is important because a monomer
         may have uneven stoichiometry (e.g., OH=1 but COOH=2), which can enable
         branching or different pairing possibilities later.
       • This stage still does not enumerate pairings; it only guarantees both
         chemical capabilities exist and preserve the information required for
         Later enumeration.
* What this module guarantees downstream
- Every monomer returned by this detector has at least one valid reactive motif
  according to the category rules above.
- Step-growth “di_identical” monomers have at least two identical reactive sites,
  reducing the risk of chain termination artifacts.
- “di_different” monomers contain both required group types.
- Reactive site multiplicities are explicitly captured (counts), enabling later
  logic to generate many possible reactions for monomers with multiple sites.
* What this module does NOT attempt to do (by design)
    - It does not enumerate all possible site-to-site reaction pairings.
    - It does not resolve conflicts between overlapping SMARTS matches.
* Output concept (high-level)
For each monomer index, the detector returns:
- "smiles": the original monomer SMILES
- One or more detected functional group entries (indexed 1..N), each including:
    • functionality_type (vinyl / mono / di_identical / di_different)
    • functional_group_name (e.g., diol, di_carboxylic_acid, epoxide, vinyl)
    • functional_group_smarts_1 (+ smarts_2 if applicable)
    • functional_count_1 (+ count_2 if applicable)
* This output is the “functional inventory” that the polymerization engine uses
later to determine how many reactions are possible and how to build reaction
templates safely.
"""
from rdkit import Chem
import json

# Dictionary defining various types of monomers, including their functionality types (e.g., 'vinyl', 'mono', 'di_different', 'di_identical'),
# SMARTS patterns for substructure matching, and group names for identification.
# Additional functional groups can be added based on references like "J. Chem. Inf. Model. 2023, 63, 5539−5548".
# Note: Potential for monomers with mixed groups like COCl and COOH is unaddressed.
# Can be many more functional groups added here
monomer_types = {
    "hydroxy_carboxylic_acid_monomer": {
        "functionality_type": "di_different",
        "smarts_1": "[OX2H1;!$(OC=*):1]",
        "smarts_2": "[CX3:2](=[O])[OX2H1]",
        "group_name": "hydroxy_carboxylic_acid"
    },
    "hydroxy_acid_halides_monomer": {
        "functionality_type": "di_different",
        "smarts_1": "[OX2H1;!$(OC=*):1]",
        "smarts_2": "[CX3:2](=[O])[Cl,Br,I]",
        "group_name": "hydroxy_acid_halide",
        "comments": "Hydroxy acid halides are highly reactive and less commonly used monomers for polyesterification compared to hydroxy carboxylic acids."
    },  
    "diol_monomer": {
        "functionality_type": "di_identical",
        "smarts_1": "[O,S;X2;H1;!$([O,S]C=*):3]",
        "group_name": "diol",
        "comments": None,
     },
    "di_carboxylic_acid_halide_monomer": {
        "functionality_type": "di_identical",
        "smarts_1": "[CX3:1](=[O])[Cl,Br,I]",
        "group_name": "di_carboxylic_acid_halide",
        "comments": None,
    },
    "di_carboxylic_acid_monomer": {
        "functionality_type": "di_identical",
        "smarts_1": "[CX3:1](=[O])[OX2H1]",
        "group_name": "di_carboxylic_acid",
        "comments": None,
    },
    # "vinyl_monomer": {
    #     "functionality_type": "vinyl",
    #     "smarts_1": "[C]=[C;D1]",
    #     "group_name": "vinyl"
    # },
    # "cyclic_olefin_monomer": {
    #     "functionality_type": "vinyl",
    #     "smarts_1": "[CX3;R:1]=[CX3;R:2]",
    #     "group_name": "cyclic_olefin"
    # },
    # "lactone_monomer": {
    #     "functionality_type": "mono",
    #     "smarts_1": "[CX3;R:1](=[OX1])[OX2;R:2]",
    #     "group_name": "lactone"
    # },
    # "cyclic_anhydride_monomer": {
    #     "functionality_type": "mono",
    #     "smarts_1": "[C,c;R:1][CX3,c;R](=[OX1])[OX2,o;R][CX3,c;R](=[OX1])[C,c;R:2]",
    #     "group_name": "cyclic_anhydride"
    # },
    # "epoxide_monomer": {
    #     "functionality_type": "mono",
    #     "smarts_1": "[CX4;R:3]1[OX2;R:4][CX4;R:5]1",
    #     "group_name": "epoxide"
    # },
    # "lactam_monomer": {
    #     "functionality_type": "mono",
    #     "smarts_1": "[CX3;R:1](=[OX1])[NX3;R:2]",
    #     "group_name": "lactam"
    # },
    # "amino_acid_monomer": {
    #     "functionality_type": "di_different",
    #     "smarts_1": "[NH2;!$(NC=O)]",
    #     "smarts_2": "[CX3:2](=[O])[OX2H1]",
    #     "group_name": "amino_acid"
    # },
    # "di_amine_monomer": {
    #     "functionality_type": "di_identical",
    #     "smarts_1": "[N&X3;H2,H1;!$(NC=*):3]",
    #     "group_name": "di_amine"
    # },
    # "primery_di_amine_monomer": {
    #     "functionality_type": "di_identical",
    #     "smarts_1": "[C,c:6][NX3;H2;!$(N[C,S]=*)]",
    #     "group_name": "di_primery_amine"
    # },
    # "di_cyclic_anhydride_monomer": {
    #     "functionality_type": "di_identical",
    #     "smarts_1": "[CX3,c;R:1](=[OX1])[OX2,o;R][CX3,c;R:2](=[OX1])",
    #     "group_name": "di_cyclic_anhydride"
    # },
    # "di_isocyanate_monomer": {
    #     "functionality_type": "di_identical",
    #     "smarts_1": "[NX2:1]=[CX2]=[OX1,SX1:2]",
    #     "group_name": "di_isocyanate"
    # },
    # "di_epoxide_monomer": {
    #     "functionality_type": "di_identical",
    #     "smarts_1": "[CX4;H2,H1,H0;R:1]1[OX2;R:2][CX4;H1,H0;R:3]1",
    #     "group_name": "di_epoxide"
    # }
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
        Duplicate definition present in original code; this is the second instance.
    """
    # Convert SMILES to RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0, None, None

    # Create pattern from primary SMARTS and check for initial match
    patt1 = Chem.MolFromSmarts(smarts_1)
    if patt1 is None:
        return 0, None, None

    if smarts_2:
        # Handle 'di_different' types with two distinct patterns
        patt2 = Chem.MolFromSmarts(smarts_2)
        if patt2 is None:
            return 0, None, None

        # Count all occurrences of each pattern
        functional_count_1 = len(mol.GetSubstructMatches(patt1))
        functional_count_2 = len(mol.GetSubstructMatches(patt2))

        if functional_count_1 >= 1 and functional_count_2 >= 1:
            return 2, functional_count_1, functional_count_2
        else:
            return 0, functional_count_1, functional_count_2
        
    else:
        # Count occurrences of primary pattern
        functional_count_1 = len(mol.GetSubstructMatches(patt1))

        if functionality_type == "di_identical":
            # Check if the match occurs more than once for di_identical functionality
            if functional_count_1 >= 2:
                return 2, functional_count_1, None
            else:
                return 0, functional_count_1, None
        
        # For vinyl or mono
        if functional_count_1 >= 1:
            return 1, functional_count_1, None
        else:
            return 0, 0, None


def functional_groups_detector(monomer_dictionary: dict) -> dict:
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
            # Handle 'vinyl' functionality (single match required)
            if functional_group["functionality_type"] == "vinyl":
                functionality, functional_count_1, functional_count_2 = detect_monomer_functionality(
                    smiles,
                    functional_group["functionality_type"],
                    functional_group["smarts_1"],
                )
                if functionality >= 1:
                    functional_group_index += 1
                    if functional_group["comments"]:
                        print(f"Note: Monomer {indexm} ({smiles}) - {functional_group['comments']}")
                    print(f"Monomer {indexm} ({smiles}) has functionality: {functional_group['group_name']}")
                    if indexm not in selected_monomers:
                        selected_monomers[indexm] = {"smiles": smiles}
                    selected_monomers[indexm][functional_group_index] = {
                        "functionality_type": "vinyl",
                        "functional_group_name": functional_group["group_name"],
                        "functional_group_smarts_1": functional_group["smarts_1"],
                        "functional_count_1": functional_count_1
                    }

            # Handle 'mono' functionality (single match required)
            if functional_group["functionality_type"] == "mono":
                functionality, functional_count_1, functional_count_2 = detect_monomer_functionality(
                    smiles,
                    functional_group["functionality_type"],
                    functional_group["smarts_1"],
                )
                if functionality >= 1:
                    if functional_group["comments"]:
                        print(f"Note: Monomer {indexm} ({smiles}) - {functional_group['comments']}")
                    functional_group_index += 1
                    print(f"Monomer {indexm} ({smiles}) has functionality: {functional_group['group_name']}")
                    if indexm not in selected_monomers:
                        selected_monomers[indexm] = {"smiles": smiles}
                    selected_monomers[indexm][functional_group_index] = {
                        "functionality_type": "mono",
                        "functional_group_name": functional_group["group_name"],
                        "functional_group_smarts_1": functional_group["smarts_1"],
                        "functional_count_1": functional_count_1
                    }
                

            # Handle 'di_different' functionality (matches for both smarts_1 and smarts_2 required)
            if functional_group["functionality_type"] == "di_different":
                functionality, functional_count_1, functional_count_2 = detect_monomer_functionality(
                    smiles,
                    functional_group["functionality_type"],
                    functional_group["smarts_1"],
                    functional_group["smarts_2"],
                )
                if functionality >= 2:
                    if functional_group["comments"]:
                        print(f"Note: Monomer {indexm} ({smiles}) - {functional_group['comments']}")
                    functional_group_index += 1
                    print(f"Monomer {indexm} ({smiles}) has functionality: {functional_group['group_name']}")
                    if indexm not in selected_monomers:
                        selected_monomers[indexm] = {"smiles": smiles}
                    selected_monomers[indexm][functional_group_index] = {
                        "functionality_type": "di_different",
                        "functional_group_name": functional_group["group_name"],
                        "functional_group_smarts_1": functional_group["smarts_1"],
                        "functional_count_1": functional_count_1,
                        "functional_group_smarts_2": functional_group["smarts_2"],
                        "functional_count_2": functional_count_2
                    }

            # Handle 'di_identical' functionality (at least 2 matches of smarts_1 required)
            if functional_group["functionality_type"] == "di_identical":
                functionality, functional_count_1, functional_count_2 = detect_monomer_functionality(
                    smiles,
                    functional_group["functionality_type"],
                    functional_group["smarts_1"],
                )

                if functionality >= 2:
                    # Print comments if available
                    if functional_group["comments"]:
                        print(f"Note: Monomer {indexm} ({smiles}) - {functional_group['comments']}")
                    functional_group_index += 1
                    print(f"Monomer {indexm} ({smiles}) has functionality: {functional_group['group_name']}")
                    if indexm not in selected_monomers:
                        selected_monomers[indexm] = {"smiles": smiles}
                    selected_monomers[indexm][functional_group_index] = {
                        "functionality_type": "di_identical",
                        "functional_group_name": functional_group["group_name"],
                        "functional_group_smarts_1": functional_group["smarts_1"],
                        "functional_count_1": functional_count_1,
                    }

    # for debugging purposes
    # print(json.dumps(selected_monomers, indent=4))
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
    # Run functional groups detection
    print(json.dumps(functional_groups_detector(monomer_dictionary), indent=4))

        