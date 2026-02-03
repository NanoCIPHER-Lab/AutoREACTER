import json
from typing import Dict, Any
# Dictionary defining various polymerization reactions.
# Each reaction includes metadata such as whether reactants are the same,
# the types of reactants, the product type, whether to delete atoms,
# and the SMARTS string for the reaction pattern using RDKit syntax.

reactions = {
    # "Vinyl Addition Polymerization": {
    #     "same_reactants": True,
    #     "reactant_1": "vinyl",
    #     "product": "polyvinyl_chain",
    #     "delete_atom": False,
    #     "reaction" : "[CH2:1]=[CH;H1,H0;!R:2].[CH2:3]=[CH;H1,H0;!R:4]>>[CH2:1]-[CH:2]-[CH2:3]-[CH:4]" 
    # },
    # "Cyclic Olefin Addition Polymerization": {
    #     "same_reactants": True,
    #     "reactant_1": "cyclic_olefin",   
    #     "product": "polycyclic_chain",
    #     "delete_atom": False,
    #     "reaction": "[CX3;R:1]=[CX3;R:2].[CX3;R:3]=[CX3;R:4]>>[CX4:1]-[CX4:2]-[CX4:3]-[CX4:4]"
    # },
    # "Vinyl Copolymerization": {
    #     "same_reactants": False,
    #     "reactant_1": "vinyl",
    #     "reactant_2": "vinyl",
    #     "product": "copolyvinyl_chain",
    #     "delete_atom": False,
    #     "reaction" : "[CH2:1]=[CH;H1,H0;!R:2].[CH2:3]=[CH;H1,H0;!R:4]>>[CH2:1]-[CH:2]-[CH2:3]-[CH:4]"
    # },
    # "Cyclic Olefin and Vinyl Copolymerization": {
    #     "same_reactants": False,
    #     "reactant_1": "vinyl",
    #     "reactant_2": "cyclic_olefin",
    #     "product": "copolycyclicvinyl_chain",
    #     "delete_atom": False,
    #     "reaction": "[CH2:1]=[CH;H1,H0;!R:2].[CX3;R:3]=[CX3;R:4]>>[CX4:1]-[CX4:2]-[CH2:3]-[CH:4]" # has to give reactants as in the order 
    # },
    # "Cyclic Olefin Copolymerization": {
    #     "same_reactants": False,
    #     "reactant_1": "cyclic_olefin",
    #     "reactant_2": "cyclic_olefin",
    #     "product": "copolycyclic_chain",
    #     "delete_atom": False,
    #     "reaction": "[CX3;R:1]=[CX3;R:2].[CX3;R:3]=[CX3;R:4]>>[CX4:1]-[CX4:2]-[CX4:3]-[CX4:4]"
    # },
    # "Lactone Ring-Opening Polyesterification": { # does not work yet
    #     "reactant_1": "lactone",
    #     "reactant_2": "initiator",
    #     "product": "polyester_chain"
    # },
    "Hydroxy Carboxylic Acid Polycondensation(Polyesterification)": {
        "same_reactants": True,
        "reactant_1": "hydroxy_carboxylic_acid",
        "product": "polyester_chain",
        "delete_atom": True,
        "reaction": "[O;!$(OC=*):1]-[H:3].[CX3:2](=[O:5])[OX2H1:4]>>[OX2:1]-[CX3:2](=[O:5]).[O:4]-[H:3]",
        "reference": {"smarts": "https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329",
                      "reaction_and_mechanism": ["https://pubs.acs.org/doi/10.1021/ed048pA734.1", 
                                                 "https://pubs.acs.org/doi/10.1021/ed073pA312"]}
    },
    "Hydroxy Carboxylic and Hydroxy Carboxylic Polycondensation(Polyesterification)": {
        "same_reactants": False,
        "reactant_1": "hydroxy_carboxylic_acid",
        "reactant_2": "hydroxy_carboxylic_acid",
        "product": "polyester_chain",
        "delete_atom": True,
        "reaction": "[O;!$(OC=*):1]-[H:3].[CX3:2](=[O:5])[OX2H1:4]>>[OX2:1]-[CX3:2](=[O:5]).[O:4]-[H:3]",
        "reference": {"smarts": "https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329",
                      "reaction_and_mechanism": ["https://pubs.acs.org/doi/10.1021/ed048pA734.1", 
                                                 "https://pubs.acs.org/doi/10.1021/ed073pA312"]}
    },
    # "Diol and Di-Carboxylic Acid Polycondensation(Polyesterification)": {
    #     "same_reactants": False,
    #     "reactant_1": "diol",
    #     "reactant_2": "di_carboxylic_acid",
    #     "product": "polyester_chain",   
    #     "delete_atom": True,
    #     "reaction": "[CX3:2](=[O])[OX2H1,Cl,Br].[O,S;X2;H1;!$([O,S]C=*):3]>>[CX3:2](=[O])-[O,S;X2;!$([O,S]C=*):3]"
    #     # "[CX3:2](=[O])[OX2H1,Cl,Br:1].[O,S;X2;H1;!$([O,S]C=*):3]>>[CX3:2](=[O])-[O,S;X2;!$([O,S]C=*):3].[OX2H1,Cl,Br:1]"
    #     # by product not working properly for now 
    # },
    # "Cyclic Anhydride and Epoxide Polyesterification": {
    #     "same_reactants": False,
    #     "reactant_1": "cyclic_anhydride",
    #     "reactant_2": "diol",
    #     "product": "polyester_chain",
    #     "delete_atom": False # Not sure about this
    # },
    # "Cyclic Anhydride and Epoxide Polyetherification": { # Not sure untill tested
    #     "same_reactants": False,
    #     "reactant_1": "cyclic_anhydride",
    #     "reactant_2": "epoxide",   
    #     "product": "polyether_chain",
    #     "delete_atom": False # Not sure about this
    # },
    # "Epoxide Ring-Opening Polyetherification": { # Do not have a clear idea about this reaction yet
    #     "same_reactants": False,
    #     "reactant_1": "epoxide",
    #     "reactant_2": "initiator",
    #     "product": "polyether_chain",
    #     "delete_atom": False
    # },
    # "Hindered Phenol Polyetherification": { # Same as above
    #     "same_reactants": True,
    #     "reactant_1": "hindered_phenol",
    #     "product": "polyether_chain",
    #     "delete_atom": False
    # },
    # "Hindered Phenol Hindered Phenol Polyetherification": { # Same as above
    #     "same_reactants": False,
    #     "reactant_1": "hindered_phenol",
    #     "reactant_2": "hindered_phenol",
    #     "product": "polyether_chain",
    #     "delete_atom": False
    # },
    # "Bis(p-halogenatedaryl)sulfone Diol (without thiol) polycondensation": { # Same as above
    #     "same_reactants": False,
    #     "reactant_1": "bis(p-halogenatedaryl)sulfone",
    #     "reactant_2": "diol",
    #     "product": "polyether_chain",
    #     "delete_atom": True
    # },
    # "Bis(bis(p-fluoroaryl)ketone Diol (without thiol) polycondensation": { # Same as above
    #     "same_reactants": False,
    #     "reactant_1": "(bis(p-fluoroaryl)ketone",
    #     "reactant_2": "diol",
    #     "product": "polyether_chain",
    #     "delete_atom": True
    # },
    # "Bis(bis(p-fluoroaryl)ketone Diol (without thiol) polycondensation": { # Same as above
    #     "same_reactants": False,
    #     "reactant_1": "(bis(p-fluoroaryl)ketone",
    #     "reactant_2": "diol",
    #     "product": "polyether_chain",
    #     "delete_atom": True
    # },
    # "Lactam Ring-Opening Polyamidation": { # does not work yet
    #     "same_reactants": True,
    #     "reactant_1": "lactam",
    #     "product": "polyamide_chain",
    #     "delete_atom": False
    # },
    # "Amino Acid Polycondensation (Polyamidation)": {
    #     "same_reactants": True,
    #     "reactant_1": "amino_acid",
    #     "product": "polyamide_chain",
    #     "delete_atom": True,
    #     "reaction": "[NX3;H2,H1;!$(OC=*):1].[CX3:2](=[O])[OX2H1]>>[NX3:1]-[CX3:2](=[O]).O"
    # },
    # "Amino Acid and Amino Acid Polycondensation (Polyamidation)": {
    #     "same_reactants": False,
    #     "reactant_1": "amino_acid",
    #     "reactant_2": "amino_acid",
    #     "product": "polyamide_chain",
    #     "delete_atom": True,
    #     "reaction": "[NX3;H2,H1;!$(OC=*):1].[CX3:2](=[O])[OX2H1]>>[NX3:1]-[CX3:2](=[O]).O"
    # },
    # "Di-Amine and Di-Carboxylic Acid Polycondensation (Polyamidation)": {
    #     "same_reactants": False,
    #     "reactant_1": "di_amine",
    #     "reactant_2": "di_carboxylic_acid",
    #     "product": "polyamide_chain",
    #     "delete_atom": True,
    #     "reaction": "[CX3:2](=[O])[OX2H1,Cl,Br:1].[N&X3;H2,H1;!$(NC=*):3]>>[CX3:2](=[O])-[NX3;!$(NC=*):3].[OX2H1,Cl,Br:1]" # same product issue as above in polyesterification
    # },
    # "Di-cyclic Anhydride and Di-Primery ammine Polycondensation (Polyimidation)": { # Do not have a clear idea about this reaction yet
    #     "same_reactants": False,
    #     "reactant_1": "di_cyclic_anhydride",
    #     "reactant_2": "di_isocyanate",
    #     "product": "polyurethane_chain",
    #     "delete_atom": True
    # },
#     "Di-Isocyanate and Diol Polyurethane Formation": {
#         "same_reactants": False,
#         "reactant_1": "di_isocyanate",
#         "reactant_2": "diol",
#         "product": "polyurethane_chain",
#         "delete_atom": True,
#         "reaction": "[NX3;H1,H0;!$(N[C,S]=*):1].[CX4;H2,H1;!$([CX4](=O)):2]>>[NX3:1]-[CX4:2](=O)"
#     },
#     "Di-Epoxide and Di-Isocyanate Polyamination": {
#         "same_reactants": False,
#         "reactant_1": "di_epoxide",
#         "reactant_2": "di_isocyanate",
#         "product": "polyamine_chain",
#         "delete_atom": False,
#         "reaction": "[NX2:3]=[CX2:4]=[OX1,SX1:5].[OX2,SX2;H1;!$([O,S]C=*):6]>>[NX3:3][CX3:4](=[OX1,SX1:5])[OX2,SX2;!$([O,S]C=*):6]"
#     }
# }
}


def reaction_detector(monomer_dictionary: dict) -> dict:
    """
    Detects possible polymerization reactions based on the functional groups in the provided monomer dictionary.

    This function iterates through predefined reactions and matches them against the functional groups
    in the monomers. It supports homopolymerization (same reactants) and copolymerization (different reactants).
    For each match, it records the monomer details and reaction index.

    Args:
        monomer_dictionary (dict): A dictionary where keys are monomer indices or IDs,
                                  and values are dictionaries containing monomer details,
                                  including 'smiles' and functional groups with 'functional_group_name'.

    Returns:
        dict: A dictionary of detected reactions, structured as:
              {reaction_name: {<reaction metadata keys>, rx_index: {monomer_1: {...}, monomer_2: {... (if applicable)}}}}

    Note:
        This logic matches functional_group_name against reaction reactant types.
        It does not validate the reaction SMARTS against the exact monomer structures (RDKit reaction execution is later).
    """
    detected_reactions = {}

    # Helper: get all functional-group indices in a monomer that match a target functional_group_name
    def _matching_fg_indices(monomer_entry: dict, target_group_name: str) -> list:
        fg_hits = []
        for fg_index in monomer_entry:
            if isinstance(fg_index, int):
                if monomer_entry[fg_index].get("functional_group_name") == target_group_name:
                    fg_hits.append(fg_index)
        return fg_hits

    # Iterate over each predefined reaction
    for reaction_name, reaction_info in reactions.items():
        rx_index = 0  # Counter for reaction instances for this reaction_name
        seen_pairs = set()  # Used to avoid duplicates for symmetric cases

        # Always prepare a clean metadata dict per reaction_name (prevents mutating global "reactions")
        if reaction_name not in detected_reactions:
            detected_reactions[reaction_name] = {}
            for k, v in reaction_info.items():
                detected_reactions[reaction_name][k] = v  # copy metadata only (same_reactants/reactants/product/delete_atom/reaction)

        reactant_1_name = reaction_info.get("reactant_1")
        reactant_2_name = reaction_info.get("reactant_2")  # may be None
        same_reactants = reaction_info.get("same_reactants", False)

        # Case 1: True homopolymerization definition (single reactant only; no reactant_2)
        if same_reactants and reactant_2_name is None:
            for i in monomer_dictionary:
                fg_hits_i = _matching_fg_indices(monomer_dictionary[i], reactant_1_name)
                for fg_index_i in fg_hits_i:
                    # Each functional group hit is a valid candidate "instance"
                    rx_index += 1
                    detected_reactions[reaction_name][rx_index] = {}
                    detected_reactions[reaction_name][rx_index]["monomer_1"] = {}
                    detected_reactions[reaction_name][rx_index]["monomer_1"]["smiles"] = monomer_dictionary[i]["smiles"]
                    detected_reactions[reaction_name][rx_index]["monomer_1"].update(monomer_dictionary[i][fg_index_i])  # attach fg metadata

        # Case 2: Two-reactant reaction definition (reactant_2 exists)
        else:
            # Build candidate lists for reactant_1 and reactant_2 across all monomers
            reactant_1_candidates = []  # list of tuples: (monomer_id, fg_index)
            reactant_2_candidates = []  # list of tuples: (monomer_id, fg_index)

            for i in monomer_dictionary:
                fg_hits_i_r1 = _matching_fg_indices(monomer_dictionary[i], reactant_1_name)
                for fg_index_i in fg_hits_i_r1:
                    reactant_1_candidates.append((i, fg_index_i))

                if reactant_2_name is not None:
                    fg_hits_i_r2 = _matching_fg_indices(monomer_dictionary[i], reactant_2_name)
                    for fg_index_i in fg_hits_i_r2:
                        reactant_2_candidates.append((i, fg_index_i))

            # If reactant_2 is missing for some reason, nothing to pair
            if reactant_2_name is None:
                continue

            # Pair candidates (cartesian product) so we don't miss any valid combinations
            for (i, fg_index_i) in reactant_1_candidates:
                for (j, fg_index_j) in reactant_2_candidates:
                    # Do not pair the same monomer with itself (keeps your current behavior)
                    if i == j:
                        continue

                    # Duplicate control:
                    # If reactant_1 and reactant_2 are the SAME type, reaction is symmetric in practice.
                    # So we treat (i,fg_i,j,fg_j) same as (j,fg_j,i,fg_i) and store only one.
                    if reactant_1_name == reactant_2_name:
                        ordered = tuple(sorted([(i, fg_index_i), (j, fg_index_j)]))
                        pair_key = (reaction_name, ordered[0], ordered[1])
                    else:
                        # If types differ, keep directionality (reactant_1 -> monomer_1, reactant_2 -> monomer_2)
                        pair_key = (reaction_name, (i, fg_index_i), (j, fg_index_j))

                    if pair_key in seen_pairs:
                        continue
                    seen_pairs.add(pair_key)

                    rx_index += 1
                    detected_reactions[reaction_name][rx_index] = {}
                    detected_reactions[reaction_name][rx_index]["monomer_1"] = {}
                    detected_reactions[reaction_name][rx_index]["monomer_1"]["smiles"] = monomer_dictionary[i]["smiles"]
                    detected_reactions[reaction_name][rx_index]["monomer_1"].update(monomer_dictionary[i][fg_index_i])

                    detected_reactions[reaction_name][rx_index]["monomer_2"] = {}
                    detected_reactions[reaction_name][rx_index]["monomer_2"]["smiles"] = monomer_dictionary[j]["smiles"]
                    detected_reactions[reaction_name][rx_index]["monomer_2"].update(monomer_dictionary[j][fg_index_j])

    # Clean up: remove reaction_name entries with no detected instances
    cleaned_detected_reactions = {}
    for reaction_name, rx_block in detected_reactions.items():
        has_any_instance = False
        for k in rx_block:
            if isinstance(k, int):
                has_any_instance = True
                break
        if has_any_instance:
            cleaned_detected_reactions[reaction_name] = rx_block

    return cleaned_detected_reactions


def reaction_arranger(monomer_dictionary: dict) -> dict:
    """
    Arranges and prints detected reactions from the monomer dictionary.

    This function processes the output of reaction_detector, organizes reactions by index,
    and prints user-friendly messages about identified homomonomers or comonomers.
    It restructures the data for easier selection in subsequent steps.

    Args:
        monomer_dictionary (dict): The same dictionary as input to reaction_detector.

    Returns:
        dict: An arranged dictionary of reactions, structured as:
              {reaction_name_index: {reaction_name: str, monomer_1: {...}, monomer_2: {... (if applicable)}, ...}}
    """
    detected_reactions = reaction_detector(monomer_dictionary)
    arranged_reactions = {}
    reaction_name_index = 0

    for reaction_name, reaction_info in detected_reactions.items():
        for rx_index in reaction_info:
            if isinstance(rx_index, int):
                monomer_1 = reaction_info[rx_index]["monomer_1"]

                try:
                    monomer_2 = reaction_info[rx_index]["monomer_2"]
                except KeyError:
                    monomer_2 = None

                # Homomonomer case (single monomer in entry)
                if monomer_2 is None:
                    reaction_name_index += 1
                    monomer_1_smiles = monomer_1["smiles"]
                    print(f"{reaction_name_index}. {reaction_name} homomonomer identified: {monomer_1_smiles}")

                    arranged_reactions[reaction_name_index] = {}
                    arranged_reactions[reaction_name_index]["reaction_name"] = reaction_name

                    # Copy metadata keys (non-int keys only)
                    for rx_index_meta, rx_info in reaction_info.items():
                        if not isinstance(rx_index_meta, int):
                            arranged_reactions[reaction_name_index][rx_index_meta] = rx_info

                    arranged_reactions[reaction_name_index]["monomer_1"] = monomer_1

                # Comonomer / two-monomer case
                else:
                    reaction_name_index += 1
                    monomer_1_smiles = monomer_1["smiles"]
                    monomer_2_smiles = monomer_2["smiles"]
                    print(f"{reaction_name_index}. {reaction_name} comonomers identified: {monomer_1_smiles} and {monomer_2_smiles}.")

                    arranged_reactions[reaction_name_index] = {}
                    arranged_reactions[reaction_name_index]["reaction_name"] = reaction_name

                    # Copy metadata keys (non-int keys only)
                    for rx_index_meta, rx_info in reaction_info.items():
                        if not isinstance(rx_index_meta, int):
                            arranged_reactions[reaction_name_index][rx_index_meta] = rx_info

                    arranged_reactions[reaction_name_index]["monomer_1"] = monomer_1
                    arranged_reactions[reaction_name_index]["monomer_2"] = monomer_2

    return arranged_reactions


def reaction_selector(monomer_dictionary: Dict[str, Any]) -> Dict[int, Any]:
    """
    Allows the user to select specific reactions from the arranged list.

    This function uses reaction_arranger to get the list of detected and arranged reactions,
    then prompts the user for indices to select via input. It filters and returns only the selected reactions.

    Rules:
    - Re-asks if any token is not a positive integer.
    - Re-asks if any index is not present in arranged_reactions.
    - Re-asks if the input is empty.

    Args:
        monomer_dictionary: Input dictionary containing monomers + metadata.

    Returns:
        A dictionary containing only the user-selected reactions, keyed by their indices.
        If no reactions are detected, raises a ValueError.
    """
    arranged_reactions = reaction_arranger(monomer_dictionary)
    
    if not arranged_reactions:
        raise ValueError("No possible reactions detected from the given monomers.")

    valid_indices = sorted(arranged_reactions.keys())

    while True:
        raw = input(
            f"Enter reaction indices (comma-separated). Valid options: {valid_indices}\n> "
        ).strip()

        if not raw:
            print("Empty input. Please enter at least one index.")
            continue

        parts = [p.strip() for p in raw.split(",") if p.strip()]
        bad_tokens = [p for p in parts if not p.isdigit()]

        if bad_tokens:
            print(f"Invalid token(s): {bad_tokens}. Use only integers like: 1,2,3")
            continue

        selected_indices: list[int] = []
        for p in parts:
            # safe because since checked isdigit()
            selected_indices.append(int(p))

        out_of_range: list[int] = []
        for idx in selected_indices:
            if idx not in arranged_reactions:
                out_of_range.append(idx)

        if out_of_range:
            print(f"Out of range / not available: {out_of_range}. Valid: {valid_indices}")
            continue

        # If we got here, everything is valid
        selected_reactions: Dict[int, Any] = {}
        for idx in selected_indices:
            selected_reactions[idx] = arranged_reactions[idx]

        return selected_reactions




if __name__ == "__main__":
    # Example monomer dictionary for testing
    monomer_dictionary = {
        1: {
            "smiles": "C=CCO",
            1: {
                "functionality_type": "vinyl",
                "functional_group_name": "vinyl",
                "functional_hroup_smarts_1": "[C]=[C;D1]",
                "functional_count_1": 1
            }
        },

        2: {
            "smiles": "CCCCCCC=C",
            1: {
                "functionality_type": "vinyl",
                "functional_group_name": "vinyl",
                "functional_hroup_smarts_1": "[C]=[C;D1]",
                "functional_count_1": 1
            }
        },

        3: {
            "smiles": " CC=C1CC2CC1C=C2",
            1: {
                "functionality_type": "mono",
                "functional_group_name": "cyclic_olefin",
                "functional_hroup_smarts_1": "[CX3;R:1]=[CX3;R:2]",
                "functional_count_1": 1
            }
        },

        4: {
            "smiles": "C1=CC2CCC1C2",
            1: {
                "functionality_type": "mono",
                "functional_group_name": "cyclic_olefin",
                "functional_hroup_smarts_1": "[CX3;R:1]=[CX3;R:2]",
                "functional_count_1": 1
            }
        },

        5: {
            "smiles": "C1CCC(=O)OCC1",
            1: {
                "functionality_type": "mono",
                "functional_group_name": "lactone",
                "functional_hroup_smarts_1": "[CX3;R:1](=[OX1])[OX2;R:2]",
                "functional_count_1": 1
            }
        },

        6: {
            "smiles": "OCCC(O)=O",
            1: {
                "functionality_type": "di_different",
                "functional_group_name": "hydroxy_carboxylic_acid",
                "functional_hroup_smarts_1": "[OX2H1;!$(OC=*):1]",
                "functional_count_1": 1,
                "functional_hroup_smarts_2": "[CX3:2](=[O])[OX2H1]",
                "functional_count_2": 1
            }
        },

        7: {
            "smiles": "C1=CC=C(C(=C1)C(=O)O)O",
            1: {
                "functionality_type": "di_different",
                "functional_group_name": "hydroxy_carboxylic_acid",
                "functional_hroup_smarts_1": "[OX2H1;!$(OC=*):1]",
                "functional_count_1": 1,
                "functional_hroup_smarts_2": "[CX3:2](=[O])[OX2H1]",
                "functional_count_2": 1
            }
        },

        8: {
            "smiles": "O=C(O)CCCCC(=O)O",
            1: {
                "functionality_type": "di_identical",
                "functional_group_name": "di_carboxylic_acid",
                "functional_hroup_smarts_1": "[CX3:1](=[O])[OX2H1]",
                "functional_count_1": 2
            }
        },

        9: {
            "smiles": "O=C(Cl)c1ccc(cc1)C(=O)Cl",
            1: {
                "functionality_type": "di_identical",
                "functional_group_name": "di_carboxylic_acidhalide",
                "functional_hroup_smarts_1": "[CX3:1](=[O])[Cl,Br]",
                "functional_count_1": 2
            }
        },

        10: {
            "smiles": "OCCO",
            1: {
                "functionality_type": "di_identical",
                "functional_group_name": "diol",
                "functional_hroup_smarts_1": "[O,S;X2;H1;!$([O,S]C=*):3]",
                "functional_count_1": 2
            }
        },

        11: {
            "smiles": "O=C1OC(=O)C=C1",
            1: {
                "functionality_type": "mono",
                "functional_group_name": "cyclic_olefin",
                "functional_hroup_smarts_1": "[CX3;R:1]=[CX3;R:2]",
                "functional_count_1": 1
            },
            2: {
                "functionality_type": "mono",
                "functional_group_name": "lactone",
                "functional_hroup_smarts_1": "[CX3;R:1](=[OX1])[OX2;R:2]",
                "functional_count_1": 2
            },
            3: {
                "functionality_type": "mono",
                "functional_group_name": "cyclic_anhydride",
                "functional_hroup_smarts_1":
                    "[C,c;R:1][CX3,c;R](=[OX1])[OX2,o;R]"
                    "[CX3,c;R](=[OX1])[C,c;R:2]",
                "functional_count_1": 1
            }
        },

        12: {
            "smiles": "C1CO1",
            1: {
                "functionality_type": "mono",
                "functional_group_name": "epoxide",
                "functional_hroup_smarts_1": "[CX4;R:3]1[OX2;R:4][CX4;R:5]1",
                "functional_count_1": 1
            }
        },

        13: {
            "smiles": "O=C1CCCCN1",
            1: {
                "functionality_type": "mono",
                "functional_group_name": "lactam",
                "functional_hroup_smarts_1": "[CX3;R:1](=[OX1])[NX3;R:2]",
                "functional_count_1": 1
            }
        },

        14: {
            "smiles": "NCC(=O)O",
            1: {
                "functionality_type": "di_different",
                "functional_group_name": "amino_acid",
                "functional_hroup_smarts_1": "[NH2;!$(NC=O)]",
                "functional_count_1": 1,
                "functional_hroup_smarts_2": "[CX3:2](=[O])[OX2H1]",
                "functional_count_2": 1
            }
        },

        15: {
            "smiles": "NCCCC(=O)O",
            1: {
                "functionality_type": "di_different",
                "functional_group_name": "amino_acid",
                "functional_hroup_smarts_1": "[NH2;!$(NC=O)]",
                "functional_count_1": 1,
                "functional_hroup_smarts_2": "[CX3:2](=[O])[OX2H1]",
                "functional_count_2": 1
            }
        },

        16: {
            "smiles": "NCCN",
            1: {
                "functionality_type": "di_identical",
                "functional_group_name": "di_amine",
                "functional_hroup_smarts_1": "[N&X3;H2,H1;!$(NC=*):3]",
                "functional_count_1": 2
            },
            2: {
                "functionality_type": "di_identical",
                "functional_group_name": "di_primery_amine",
                "functional_hroup_smarts_1": "[C,c:6][NX3;H2;!$(N[C,S]=*)]",
                "functional_count_1": 2
            }
        },

        17: {
            "smiles": "NCCCCCCN",
            1: {
                "functionality_type": "di_identical",
                "functional_group_name": "di_amine",
                "functional_hroup_smarts_1": "[N&X3;H2,H1;!$(NC=*):3]",
                "functional_count_1": 2
            },
            2: {
                "functionality_type": "di_identical",
                "functional_group_name": "di_primery_amine",
                "functional_hroup_smarts_1": "[C,c:6][NX3;H2;!$(N[C,S]=*)]",
                "functional_count_1": 2
            }
        },

        18: {
            "smiles": "O=C1OC(=O)c2cc(C3OC(=O)C=CC3=O)ccc2C1=O",
            1: {
                "functionality_type": "mono",
                "functional_group_name": "cyclic_olefin",
                "functional_hroup_smarts_1": "[CX3;R:1]=[CX3;R:2]",
                "functional_count_1": 1
            },
            2: {
                "functionality_type": "mono",
                "functional_group_name": "lactone",
                "functional_hroup_smarts_1": "[CX3;R:1](=[OX1])[OX2;R:2]",
                "functional_count_1": 3
            },
            3: {
                "functionality_type": "mono",
                "functional_group_name": "cyclic_anhydride",
                "functional_hroup_smarts_1":
                    "[C,c;R:1][CX3,c;R](=[OX1])[OX2,o;R]"
                    "[CX3,c;R](=[OX1])[C,c;R:2]",
                "functional_count_1": 1
            }
        },

        19: {
            "smiles": "O=C=N-C-N=C=O",
            1: {
                "functionality_type": "di_identical",
                "functional_group_name": "di_isocyanate",
                "functional_hroup_smarts_1": "[NX2:1]=[CX2]=[OX1,SX1:2]",
                "functional_count_1": 2
            }
        },

        20: {
            "smiles": "C1=CC(=CC=C1C(C2=CC=CC=C2)(O)CC3CO3)CC4CO4",
            1: {
                "functionality_type": "mono",
                "functional_group_name": "epoxide",
                "functional_hroup_smarts_1": "[CX4;R:3]1[OX2;R:4][CX4;R:5]1",
                "functional_count_1": 2
            },
            2: {
                "functionality_type": "di_identical",
                "functional_group_name": "di_epoxide",
                "functional_hroup_smarts_1":
                    "[CX4;H2,H1,H0;R:1]1[OX2;R:2][CX4;H1,H0;R:3]1",
                "functional_count_1": 2
            }
        }
    }
    import json
    print(json.dumps(reaction_selector(monomer_dictionary), indent=4))

