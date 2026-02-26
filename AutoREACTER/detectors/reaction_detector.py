
# Dictionary defining various polymerization reactions.
# Each reaction includes metadata such as whether reactants are the same,
# the types of reactants, the product type, whether to delete atoms,
# and the SMARTS string for the reaction pattern using RDKit syntax.

"""
TODO: Missing Polymerization Mechanisms
- Polycarbonates: Diols + Phosgene or Diphenyl Carbonate
- Polyureas: Diamines + Diisocyanates
- Aromatic Polyethers (PEEK/Sulfones): Activated dihalides + Bisphenols
- Aromatic Polyimides: Dianhydrides + Diamines
- Polybenzimidazoles (PBI): Tetraamines + Dicarboxylates
- Phenol-Formaldehyde (Bakelite): Phenol + Formaldehyde
- Polysiloxanes (Silicones): Hydrolysis/Condensation of Dichlorosilanes
- Polysulfides: Dihalides + Sodium Sulfide
- Thiol-Ene Click Polymerizations
- Ring-Opening Metathesis Polymerization (ROMP)
- Dendritic Polymers: Random Hyperbranched and Dendrimers
- Enzymatic Polymerizations: In Vivo / In Vitro biocatalysis
- Cycloaddition (Four-Center) Reactions
- Spiro Polymers
- Pseudopolyrotaxanes and Polyrotaxanes
- Polymerization in Supercritical Carbon Dioxide
"""

import json
from typing import Dict, Any, Tuple, Optional
from dataclasses import dataclass

try:
    from reactions_library import ReactionLibrary
except (ImportError, ModuleNotFoundError):
    from .reactions_library import ReactionLibrary

try:
    from .functional_groups_detector import FunctionalGroupInfo, MonomerRole
except (ImportError, ModuleNotFoundError):
    from functional_groups_detector import FunctionalGroupInfo, MonomerRole


@dataclass(slots=True, frozen=True)
class FunctionalGroupInfo:
    functionality_type: str
    fg_name: str
    fg_smarts_1: str
    fg_count_1: int
    fg_smarts_2: Optional[str] = None
    fg_count_2: Optional[int] = None

@dataclass(slots=True, frozen=True)
class MonomerRole:
    smiles: str
    name: str
    functionalities: Tuple[FunctionalGroupInfo, ...]

@dataclass(slots=True, frozen=True)
class ReactionTemplate:
    reaction_name: str
    reactant_1: str                 # functional group name
    reactant_2: Optional[str]       # functional group name or None
    reaction_smarts: str
    same_reactants: bool
    delete_atom: bool
    references: dict

@dataclass(slots=True, frozen=True)
class ReactionInstance:
    reaction_name: str
    reaction_smarts: str
    delete_atom: bool
    references: dict

    monomer_1: MonomerRole
    fg_1: FunctionalGroupInfo

    monomer_2: Optional[MonomerRole] = None
    fg_2: Optional[FunctionalGroupInfo] = None


class ReactionDetector:
    def __init__(self):
        self.reactions = ReactionLibrary().reactions

    def reaction_detector(self, monomer_dictionary: dict) -> dict:
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
        self.detected_reactions = {}

        # Helper: get all functional-group indices in a monomer that match a target functional_group_name
        def _matching_fg_indices(monomer_entry: dict, target_group_name: str) -> list:
            fg_hits = []
            for fg_index in monomer_entry:
                if isinstance(fg_index, int):
                    if monomer_entry[fg_index].get("functional_group_name") == target_group_name:
                        fg_hits.append(fg_index)
            return fg_hits

        # Iterate over each predefined reaction
        for reaction_name, reaction_info in self.reactions.items():
            rx_index = 0  # Counter for reaction instances for this reaction_name
            seen_pairs = set()  # Used to avoid duplicates for symmetric cases

            # Always prepare a clean metadata dict per reaction_name (prevents mutating global "reactions")
            if reaction_name not in self.detected_reactions:
                self.detected_reactions[reaction_name] = {}
                for k, v in reaction_info.items():
                    self.detected_reactions[reaction_name][k] = v  # copy metadata only (same_reactants/reactants/product/delete_atom/reaction)

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
                        self.detected_reactions[reaction_name][rx_index] = {}
                        self.detected_reactions[reaction_name][rx_index]["monomer_1"] = {}
                        self.detected_reactions[reaction_name][rx_index]["monomer_1"]["smiles"] = monomer_dictionary[i]["smiles"]
                        self.detected_reactions[reaction_name][rx_index]["monomer_1"].update(monomer_dictionary[i][fg_index_i])  # attach fg metadata

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
                        self.detected_reactions[reaction_name][rx_index] = {}
                        self.detected_reactions[reaction_name][rx_index]["monomer_1"] = {}
                        self.detected_reactions[reaction_name][rx_index]["monomer_1"]["smiles"] = monomer_dictionary[i]["smiles"]
                        self.detected_reactions[reaction_name][rx_index]["monomer_1"].update(monomer_dictionary[i][fg_index_i])

                        self.detected_reactions[reaction_name][rx_index]["monomer_2"] = {}
                        self.detected_reactions[reaction_name][rx_index]["monomer_2"]["smiles"] = monomer_dictionary[j]["smiles"]
                        self.detected_reactions[reaction_name][rx_index]["monomer_2"].update(monomer_dictionary[j][fg_index_j])

        # Clean up: remove reaction_name entries with no detected instances
        self.cleaned_detected_reactions = {}
        for reaction_name, rx_block in self.detected_reactions.items():
            has_any_instance = False
            for k in rx_block:
                if isinstance(k, int):
                    has_any_instance = True
                    break
            if has_any_instance:
                self.cleaned_detected_reactions[reaction_name] = rx_block

        return self.cleaned_detected_reactions


    def reaction_arranger(self, monomer_dictionary: dict) -> dict:
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
        detected_reactions = self.reaction_detector(monomer_dictionary)
        self.arranged_reactions = {}
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
                        # Optional comment display
                        try:
                            comment = reaction_info.get("comments", None)
                        except Exception:
                            comment = None
                        if comment:
                            print(f"   Note: {comment}")
                    
                        self.arranged_reactions[reaction_name_index] = {}
                        self.arranged_reactions[reaction_name_index]["reaction_name"] = reaction_name

                        # Copy metadata keys (non-int keys only)
                        for rx_index_meta, rx_info in reaction_info.items():
                            if not isinstance(rx_index_meta, int):
                                self.arranged_reactions[reaction_name_index][rx_index_meta] = rx_info

                        self.arranged_reactions[reaction_name_index]["monomer_1"] = monomer_1

                    # Comonomer / two-monomer case
                    else:
                        reaction_name_index += 1
                        monomer_1_smiles = monomer_1["smiles"]
                        monomer_2_smiles = monomer_2["smiles"]
                        print(f"{reaction_name_index}. {reaction_name} comonomers identified: {monomer_1_smiles} and {monomer_2_smiles}.")

                        # Optional comment display
                        try:
                            comment = reaction_info.get("comments", None)
                        except Exception:
                            comment = None
                        if comment:
                            print(f"   Note: {comment}")

                        self.arranged_reactions[reaction_name_index] = {}
                        self.arranged_reactions[reaction_name_index]["reaction_name"] = reaction_name

                        # Copy metadata keys (non-int keys only)
                        for rx_index_meta, rx_info in reaction_info.items():
                            if not isinstance(rx_index_meta, int):
                                self.arranged_reactions[reaction_name_index][rx_index_meta] = rx_info

                        self.arranged_reactions[reaction_name_index]["monomer_1"] = monomer_1
                        self.arranged_reactions[reaction_name_index]["monomer_2"] = monomer_2

        return self.arranged_reactions


    def reaction_selector(self, monomer_dictionary: Dict[str, Any]) -> Dict[int, Any]:
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
        arranged_reactions = self.reaction_arranger(monomer_dictionary)
        
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
    detector = ReactionDetector()
    # Example monomer dictionary for testing
    monomer_dictionary = {
        1: {
            "smiles": "C=CCO",
            1: {
                "functionality_type": "vinyl",
                "functional_group_name": "vinyl",
                "functional_group_smarts_1": "[C]=[C;D1]",
                "functional_count_1": 1
            }
        },

        2: {
            "smiles": "CCCCCCC=C",
            1: {
                "functionality_type": "vinyl",
                "functional_group_name": "vinyl",
                "functional_group_smarts_1": "[C]=[C;D1]",
                "functional_count_1": 1
            }
        },

        3: {
            "smiles": " CC=C1CC2CC1C=C2",
            1: {
                "functionality_type": "mono",
                "functional_group_name": "cyclic_olefin",
                "functional_group_smarts_1": "[CX3;R:1]=[CX3;R:2]",
                "functional_count_1": 1
            }
        },

        4: {
            "smiles": "C1=CC2CCC1C2",
            1: {
                "functionality_type": "mono",
                "functional_group_name": "cyclic_olefin",
                "functional_group_smarts_1": "[CX3;R:1]=[CX3;R:2]",
                "functional_count_1": 1
            }
        },

        5: {
            "smiles": "C1CCC(=O)OCC1",
            1: {
                "functionality_type": "mono",
                "functional_group_name": "lactone",
                "functional_group_smarts_1": "[CX3;R:1](=[OX1])[OX2;R:2]",
                "functional_count_1": 1
            }
        },

        6: {
            "smiles": "OCCC(O)=O",
            1: {
                "functionality_type": "di_different",
                "functional_group_name": "hydroxy_carboxylic_acid",
                "functional_group_smarts_1": "[OX2H1;!$(OC=*):1]",
                "functional_count_1": 1,
                "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
                "functional_count_2": 1
            }
        },

        7: {
            "smiles": "C1=CC=C(C(=C1)C(=O)O)O",
            1: {
                "functionality_type": "di_different",
                "functional_group_name": "hydroxy_carboxylic_acid",
                "functional_group_smarts_1": "[OX2H1;!$(OC=*):1]",
                "functional_count_1": 1,
                "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
                "functional_count_2": 1
            }
        },

        8: {
            "smiles": "O=C(O)CCCCC(=O)O",
            1: {
                "functionality_type": "di_identical",
                "functional_group_name": "di_carboxylic_acid",
                "functional_group_smarts_1": "[CX3:1](=[O])[OX2H1]",
                "functional_count_1": 2
            }
        },

        9: {
            "smiles": "O=C(Cl)c1ccc(cc1)C(=O)Cl",
            1: {
                "functionality_type": "di_identical",
                "functional_group_name": "di_carboxylic_acidhalide",
                "functional_group_smarts_1": "[CX3:1](=[O])[Cl,Br]",
                "functional_count_1": 2
            }
        },

        10: {
            "smiles": "OCCO",
            1: {
                "functionality_type": "di_identical",
                "functional_group_name": "diol",
                "functional_group_smarts_1": "[O,S;X2;H1;!$([O,S]C=*):3]",
                "functional_count_1": 2
            }
        },

        11: {
            "smiles": "O=C1OC(=O)C=C1",
            1: {
                "functionality_type": "mono",
                "functional_group_name": "cyclic_olefin",
                "functional_group_smarts_1": "[CX3;R:1]=[CX3;R:2]",
                "functional_count_1": 1
            },
            2: {
                "functionality_type": "mono",
                "functional_group_name": "lactone",
                "functional_group_smarts_1": "[CX3;R:1](=[OX1])[OX2;R:2]",
                "functional_count_1": 2
            },
            3: {
                "functionality_type": "mono",
                "functional_group_name": "cyclic_anhydride",
                "functional_group_smarts_1":
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
                "functional_group_smarts_1": "[CX4;R:3]1[OX2;R:4][CX4;R:5]1",
                "functional_count_1": 1
            }
        },

        13: {
            "smiles": "O=C1CCCCN1",
            1: {
                "functionality_type": "mono",
                "functional_group_name": "lactam",
                "functional_group_smarts_1": "[CX3;R:1](=[OX1])[NX3;R:2]",
                "functional_count_1": 1
            }
        },

        14: {
            "smiles": "NCC(=O)O",
            1: {
                "functionality_type": "di_different",
                "functional_group_name": "amino_acid",
                "functional_group_smarts_1": "[NH2;!$(NC=O)]",
                "functional_count_1": 1,
                "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
                "functional_count_2": 1
            }
        },

        15: {
            "smiles": "NCCCC(=O)O",
            1: {
                "functionality_type": "di_different",
                "functional_group_name": "amino_acid",
                "functional_group_smarts_1": "[NH2;!$(NC=O)]",
                "functional_count_1": 1,
                "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
                "functional_count_2": 1
            }
        },

        16: {
            "smiles": "NCCN",
            1: {
                "functionality_type": "di_identical",
                "functional_group_name": "di_amine",
                "functional_group_smarts_1": "[N&X3;H2,H1;!$(NC=*):3]",
                "functional_count_1": 2
            },
            2: {
                "functionality_type": "di_identical",
                "functional_group_name": "di_primery_amine",
                "functional_group_smarts_1": "[C,c:6][NX3;H2;!$(N[C,S]=*)]",
                "functional_count_1": 2
            }
        },

        17: {
            "smiles": "NCCCCCCN",
            1: {
                "functionality_type": "di_identical",
                "functional_group_name": "di_amine",
                "functional_group_smarts_1": "[N&X3;H2,H1;!$(NC=*):3]",
                "functional_count_1": 2
            },
            2: {
                "functionality_type": "di_identical",
                "functional_group_name": "di_primery_amine",
                "functional_group_smarts_1": "[C,c:6][NX3;H2;!$(N[C,S]=*)]",
                "functional_count_1": 2
            }
        },

        18: {
            "smiles": "O=C1OC(=O)c2cc(C3OC(=O)C=CC3=O)ccc2C1=O",
            1: {
                "functionality_type": "mono",
                "functional_group_name": "cyclic_olefin",
                "functional_group_smarts_1": "[CX3;R:1]=[CX3;R:2]",
                "functional_count_1": 1
            },
            2: {
                "functionality_type": "mono",
                "functional_group_name": "lactone",
                "functional_group_smarts_1": "[CX3;R:1](=[OX1])[OX2;R:2]",
                "functional_count_1": 3
            },
            3: {
                "functionality_type": "mono",
                "functional_group_name": "cyclic_anhydride",
                "functional_group_smarts_1":
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
                "functional_group_smarts_1": "[NX2:1]=[CX2]=[OX1,SX1:2]",
                "functional_count_1": 2
            }
        },

        20: {
            "smiles": "C1=CC(=CC=C1C(C2=CC=CC=C2)(O)CC3CO3)CC4CO4",
            1: {
                "functionality_type": "mono",
                "functional_group_name": "epoxide",
                "functional_group_smarts_1": "[CX4;R:3]1[OX2;R:4][CX4;R:5]1",
                "functional_count_1": 2
            },
            2: {
                "functionality_type": "di_identical",
                "functional_group_name": "di_epoxide",
                "functional_group_smarts_1":
                    "[CX4;H2,H1,H0;R:1]1[OX2;R:2][CX4;H1,H0;R:3]1",
                "functional_count_1": 2
            }
        }
    }
    import json
    detector = ReactionDetector()
    detected = detector.reaction_detector(monomer_dictionary)
    arranged = detector.reaction_arranger(monomer_dictionary)
    print("\nDetected Reactions:")
    print(json.dumps(detected, indent=2))
    # print("\nArranged Reactions:")
    # print(json.dumps(arranged, indent=2))
