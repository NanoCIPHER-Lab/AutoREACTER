import json
# Dictionary defining various polymerization reactions.
# Each reaction includes metadata such as whether reactants are the same,
# the types of reactants, the product type, whether to delete atoms,
# and the SMARTS string for the reaction pattern using RDKit syntax.
reactions = {
    "Vinyl Addition Polymerization": {
        "same_reactants": True,
        "reactant_1": "vinyl",
        "product": "polyvinyl_chain",
        "delete_atom": False,
        "reaction" : "[CH2:1]=[CH;H1,H0;!R:2].[CH2:3]=[CH;H1,H0;!R:4]>>[CH2:1]-[CH:2]-[CH2:3]-[CH:4]" 
    },
    "Cyclic Olefin Addition Polymerization": {
        "same_reactants": True,
        "reactant_1": "cyclic_olefin",   
        "product": "polycyclic_chain",
        "delete_atom": False,
        "reaction": "[CX3;R:1]=[CX3;R:2].[CX3;R:3]=[CX3;R:4]>>[CX4:1]-[CX4:2]-[CX4:3]-[CX4:4]"
    },
    "Vinyl Copolymerization": {
        "same_reactants": False,
        "reactant_1": "vinyl",
        "reactant_2": "vinyl",
        "product": "copolyvinyl_chain",
        "delete_atom": False,
        "reaction" : "[CH2:1]=[CH;H1,H0;!R:2].[CH2:3]=[CH;H1,H0;!R:4]>>[CH2:1]-[CH:2]-[CH2:3]-[CH:4]"
    },
    "Cyclic Olefin and Vinyl Copolymerization": {
        "same_reactants": False,
        "reactant_1": "vinyl",
        "reactant_2": "cyclic_olefin",
        "product": "copolycyclicvinyl_chain",
        "delete_atom": False,
        "reaction": "[CH2:1]=[CH;H1,H0;!R:2].[CX3;R:3]=[CX3;R:4]>>[CX4:1]-[CX4:2]-[CH2:3]-[CH:4]" # has to give reactants as in the order 
    },
    "Cyclic Olefin Copolymerization": {
        "same_reactants": False,
        "reactant_1": "cyclic_olefin",
        "reactant_2": "cyclic_olefin",
        "product": "copolycyclic_chain",
        "delete_atom": False,
        "reaction": "[CX3;R:1]=[CX3;R:2].[CX3;R:3]=[CX3;R:4]>>[CX4:1]-[CX4:2]-[CX4:3]-[CX4:4]"
    },
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
        "reaction": "[OX2H1;!$(OC=*):1].[CX3:2](=[O])[OX2H1]>>[OX2:1]-[CX3:2](=[O]).O"
    },
    "Hydroxy Carboxylic and Hydroxy Carboxylic Polycondensation(Polyesterification)": {
        "same_reactants": False,
        "reactant_1": "hydroxy_carboxylic_acid",
        "reactant_2": "hydroxy_carboxylic_acid",
        "product": "polyester_chain",
        "delete_atom": True,
        "reaction": "[OX2H1;!$(OC=*):1].[CX3:2](=[O])[OX2H1]>>[OX2:1]-[CX3:2](=[O]).O"
    },
    "Diol and Di-Carboxylic Acid Polycondensation(Polyesterification)": {
        "same_reactants": False,
        "reactant_1": "diol",
        "reactant_2": "di_carboxylic_acid",
        "product": "polyester_chain",   
        "delete_atom": True,
        "reaction": "[CX3:2](=[O])[OX2H1,Cl,Br].[O,S;X2;H1;!$([O,S]C=*):3]>>[CX3:2](=[O])-[O,S;X2;!$([O,S]C=*):3]"
        # "[CX3:2](=[O])[OX2H1,Cl,Br:1].[O,S;X2;H1;!$([O,S]C=*):3]>>[CX3:2](=[O])-[O,S;X2;!$([O,S]C=*):3].[OX2H1,Cl,Br:1]"
        # by product not working properly for now 
    },
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
    "Amino Acid Polycondensation (Polyamidation)": {
        "same_reactants": True,
        "reactant_1": "amino_acid",
        "product": "polyamide_chain",
        "delete_atom": True,
        "reaction": "[NX3;H2,H1;!$(OC=*):1].[CX3:2](=[O])[OX2H1]>>[NX3:1]-[CX3:2](=[O]).O"
    },
    "Amino Acid and Amino Acid Polycondensation (Polyamidation)": {
        "same_reactants": False,
        "reactant_1": "amino_acid",
        "reactant_2": "amino_acid",
        "product": "polyamide_chain",
        "delete_atom": True,
        "reaction": "[NX3;H2,H1;!$(OC=*):1].[CX3:2](=[O])[OX2H1]>>[NX3:1]-[CX3:2](=[O]).O"
    },
    "Di-Amine and Di-Carboxylic Acid Polycondensation (Polyamidation)": {
        "same_reactants": False,
        "reactant_1": "di_amine",
        "reactant_2": "di_carboxylic_acid",
        "product": "polyamide_chain",
        "delete_atom": True,
        "reaction": "[CX3:2](=[O])[OX2H1,Cl,Br:1].[N&X3;H2,H1;!$(NC=*):3]>>[CX3:2](=[O])-[NX3;!$(NC=*):3].[OX2H1,Cl,Br:1]" # same product issue as above in polyesterification
    },
    # "Di-cyclic Anhydride and Di-Primery ammine Polycondensation (Polyimidation)": { # Do not have a clear idea about this reaction yet
    #     "same_reactants": False,
    #     "reactant_1": "di_cyclic_anhydride",
    #     "reactant_2": "di_isocyanate",
    #     "product": "polyurethane_chain",
    #     "delete_atom": True
    # },
    "Di-Isocyanate and Diol Polyurethane Formation": {
        "same_reactants": False,
        "reactant_1": "di_isocyanate",
        "reactant_2": "diol",
        "product": "polyurethane_chain",
        "delete_atom": True,
        "reaction": "[NX3;H1,H0;!$(N[C,S]=*):1].[CX4;H2,H1;!$([CX4](=O)):2]>>[NX3:1]-[CX4:2](=O)"
    },
    "Di-Epoxide and Di-Isocyanate Polyamination": {
        "same_reactants": False,
        "reactant_1": "di_epoxide",
        "reactant_2": "di_isocyanate",
        "product": "polyamine_chain",
        "delete_atom": False,
        "reaction": "[NX2:3]=[CX2:4]=[OX1,SX1:5].[OX2,SX2;H1;!$([O,S]C=*):6]>>[NX3:3][CX3:4](=[OX1,SX1:5])[OX2,SX2;!$([O,S]C=*):6]"
    }
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
              {reaction_name: {rx_index: {monomer_1: {...}, monomer_2: {... (if applicable)}, ...}, ...}}

    Note:
        This is a placeholder implementation. The actual logic checks functional group names
        against reaction reactant types but does not perform structural validation with RDKit.
    """
    detected_reactions = {}
    
    for reaction_name, reaction_info in reactions.items():
        rx_index = 0
        
        # Determine reactant targets
        r1_target = reaction_info["reactant_1"]
        r2_target = reaction_info.get("reactant_2")
        
        if reaction_info["same_reactants"]:
            # Homopolymerization Case: Single monomer species reacts with itself
            for i, monomer_i in monomer_dictionary.items():
                for fg_idx, fg_data in monomer_i.items():
                    if isinstance(fg_idx, int) and fg_data["functional_group_name"] == r1_target:
                        rx_index += 1
                        if reaction_name not in detected_reactions:
                            detected_reactions[reaction_name] = reaction_info.copy()
                        
                        detected_reactions[reaction_name][rx_index] = {
                            "monomer_1": {"smiles": monomer_i["smiles"], **fg_data}
                        }
        else:
            # Copolymerization Case: Monomer i reacts with Monomer j
            for i, monomer_i in monomer_dictionary.items():
                for fg_idx_i, fg_data_i in monomer_i.items():
                    if isinstance(fg_idx_i, int) and fg_data_i["functional_group_name"] == r1_target:
                        for j, monomer_j in monomer_dictionary.items():
                            # For copolymerization, we allow same monomer ID if it has multiple functional sites
                            # or just different monomer IDs.
                            for fg_idx_j, fg_data_j in monomer_j.items():
                                if isinstance(fg_idx_j, int) and fg_data_j["functional_group_name"] == r2_target:
                                    # Prevent redundant pairs (A+B is same as B+A) if reactant types are identical
                                    if r1_target == r2_target and i > j:
                                        continue
                                    
                                    rx_index += 1
                                    if reaction_name not in detected_reactions:
                                        detected_reactions[reaction_name] = reaction_info.copy()
                                    
                                    detected_reactions[reaction_name][rx_index] = {
                                        "monomer_1": {"smiles": monomer_i["smiles"], **fg_data_i},
                                        "monomer_2": {"smiles": monomer_j["smiles"], **fg_data_j}
                                    }

    return detected_reactions

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

    Note:
        This function calls reaction_detector internally and uses print statements for user feedback.
    """
    detected_reactions = reaction_detector(monomer_dictionary)
    arranged_reactions = {}
    reaction_name_index = 0
    
    for reaction_name, reaction_info in detected_reactions.items():
        # Metadata keys (strings)
        metadata = {k: v for k, v in reaction_info.items() if not isinstance(k, int)}
        
        for rx_idx, rx_data in reaction_info.items():
            if isinstance(rx_idx, int):
                reaction_name_index += 1
                monomer_1 = rx_data["monomer_1"]
                monomer_2 = rx_data.get("monomer_2")
                
                entry = {"reaction_name": reaction_name, **metadata, "monomer_1": monomer_1}
                
                if monomer_2 is None:
                    print(f"{reaction_name_index}. {reaction_name} homomonomer identified: {monomer_1['smiles']}")
                else:
                    entry["monomer_2"] = monomer_2
                    print(f"{reaction_name_index}. {reaction_name} comonomers identified: {monomer_1['smiles']} and {monomer_2['smiles']}.")
                
                arranged_reactions[reaction_name_index] = entry
                
    return arranged_reactions

def reaction_selector(monomer_dictionary: dict) -> dict:
    """
    Allows the user to select specific reactions from the arranged list.

    This function uses reaction_arranger to get the list of detected and arranged reactions,
    then prompts the user for indices to select via input. It filters and returns only the selected reactions.

    Args:
        monomer_dictionary (dict): The same dictionary as input to previous functions.

    Returns:
        dict: A dictionary containing only the user-selected reactions, keyed by their indices.

    Note:
        This is interactive and relies on user input. For non-interactive use, consider modifying to accept indices as a parameter.
    """
    arranged_reactions = reaction_arranger(monomer_dictionary)
    selected_reactions = {}
    
    if not arranged_reactions:
        print("No reactions detected.")
        return {}

    user_input = input("Enter the reaction indices you want to select (comma-separated): ")
    selected_indices = [idx.strip() for idx in user_input.split(",") if idx.strip().isdigit()]
    
    for idx_str in selected_indices:
        idx = int(idx_str)
        if idx in arranged_reactions:
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

