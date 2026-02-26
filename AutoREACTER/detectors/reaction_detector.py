
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

import pathlib
from typing import Dict, Any, Tuple, Optional
from dataclasses import dataclass
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import Draw, rdChemReactions

try:
    from reactions_library import ReactionLibrary
except (ImportError, ModuleNotFoundError):
    from .reactions_library import ReactionLibrary

try:
    from .functional_groups_detector import FunctionalGroupInfo, MonomerRole
except (ImportError, ModuleNotFoundError):
    from functional_groups_detector import FunctionalGroupInfo, MonomerRole

class SMARTSERROR(Exception):
    """Error from developer-defined SMARTS patterns not producing expected products.
    Please contact the developer to resolve this issue. Raise a GitHub issue if you encounter this error.
    at <https://github.com/NanoCIPHER-Lab/AutoREACTER/issues>"""
    pass

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
    functional_group_1: FunctionalGroupInfo
    monomer_2: Optional[MonomerRole] = None
    functional_group_2: Optional[FunctionalGroupInfo] = None

class ReactionDetector:
    def __init__(self):
        self.reactions = ReactionLibrary().reactions
    
    def _matching_fgs(self, monomer_entry: MonomerRole, target_group_name: str) -> list:
            fg_hits = []
            for fg in monomer_entry.functionalities:
                if fg.fg_name == target_group_name:
                    fg_hits.append(fg)
            return fg_hits

    def _seen_pair_key(
            self,
            reaction_name: str,
            monomer_role_1: MonomerRole,
            fg_1: FunctionalGroupInfo,
            monomer_role_2: Optional[MonomerRole] = None,
            fg_2: Optional[FunctionalGroupInfo] = None,
        ) -> Tuple:

        if monomer_role_2 is None or fg_2 is None:
            return (reaction_name, monomer_role_1, fg_1)

        # symmetric case
        ordered = tuple(sorted(
            [(monomer_role_1, fg_1), (monomer_role_2, fg_2)],
            key=lambda x: (id(x[0]), id(x[1]))
        ))
        return (reaction_name, ordered[0], ordered[1])

    def reaction_detector(self, monomer_roles: list[MonomerRole]) -> list[ReactionInstance]:
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

        # create a new list to hold reaction instances
        reaction_instances = []
        seen_pairs = set()  # to track unique reactant pairs and avoid duplicates

        # Iterate over each predefined reaction
        for reaction_name, reaction_info in self.reactions.items():
            reactant_1_name = reaction_info.get("reactant_1")
            reactant_2_name = reaction_info.get("reactant_2")  # may be None
            same_reactants = reaction_info.get("same_reactants", False)
            
            # HOMO CASE
            if same_reactants and reactant_2_name is None:
                for monomer_role in monomer_roles: 
                    fg_hits = self._matching_fgs(monomer_role, reactant_1_name)
                    if not fg_hits:
                        continue
                    for fg in fg_hits:
                        pair_key = self._seen_pair_key(
                            reaction_name,
                            monomer_role,
                            fg
                        )

                        if pair_key in seen_pairs:
                            continue

                        seen_pairs.add(pair_key)
                        reaction_instances.append(
                            ReactionInstance(
                                reaction_name=reaction_name,
                                reaction_smarts=reaction_info["reaction"],
                                delete_atom=reaction_info["delete_atom"],
                                references=reaction_info["reference"],
                                monomer_1=monomer_role,
                                functional_group_1=fg,
                                monomer_2=None,
                                functional_group_2=None
                            )
                        )
            # COMONOMER CASE
            else:
                for monomer_role_i in monomer_roles:
                    fg_hits_i = self._matching_fgs(monomer_role_i, reactant_1_name)
                    if not fg_hits_i:
                        continue
                    for fg_i in fg_hits_i:
                        for monomer_role_j in monomer_roles:
                            if monomer_role_i == monomer_role_j:
                                continue  # skip same monomer
                            fg_hits_j = self._matching_fgs(monomer_role_j, reactant_2_name)
                            if not fg_hits_j:
                                continue
                            for fg_j in fg_hits_j:
                                # Create a unique pair identifier to avoid duplicates
                                pair_key = self._seen_pair_key(
                                            reaction_name,
                                            monomer_role_i,
                                            fg_i,
                                            monomer_role_j,
                                            fg_j
                                        )
                                if pair_key not in seen_pairs:
                                    seen_pairs.add(pair_key)
                                    reaction_instances.append(
                                        ReactionInstance(
                                            reaction_name=reaction_name,
                                            reaction_smarts=reaction_info["reaction"],
                                            delete_atom=reaction_info["delete_atom"],
                                            references=reaction_info["reference"],
                                            monomer_1=monomer_role_i,
                                            functional_group_1=fg_i,
                                            monomer_2=monomer_role_j,
                                            functional_group_2=fg_j
                                        )
                                        )

        return reaction_instances

    def create_reaction_image(self, reactant_a_smiles, reactant_b_smiles, reaction_smarts, reaction_name):

        reactant_a = Chem.MolFromSmiles(reactant_a_smiles)
        reactant_b = Chem.MolFromSmiles(reactant_b_smiles)
        reactant_a = Chem.AddHs(reactant_a)
        reactant_b = Chem.AddHs(reactant_b)
        # debug prints
        # print(Chem.MolToSmiles(reactant_a))
        # print(Chem.MolToSmiles(reactant_b))
        rxn_engine = rdChemReactions.ReactionFromSmarts(reaction_smarts)
        products_sets = rxn_engine.RunReactants((reactant_a, reactant_b))
        # debug print
        # print(f"Products sets for {reaction_smarts}")

        if not products_sets:
            products_sets = rxn_engine.RunReactants((reactant_b, reactant_a))
            if not products_sets:
                print(f"No products generated for {reaction_smarts} with either reactant order.")
                raise SMARTSERROR(f"Reaction SMARTS '{reaction_smarts}' did not produce any products for reactants '{reactant_a_smiles}' and '{reactant_b_smiles}' in either order.")


        display_rxn = rdChemReactions.ChemicalReaction()
        display_rxn.AddReactantTemplate(reactant_a)
        display_rxn.AddReactantTemplate(reactant_b)

        for prod in products_sets[0]:
            try:
                Chem.SanitizeMol(prod)
            except:
                pass
            display_rxn.AddProductTemplate(prod)

        return Draw.ReactionToImage(display_rxn, subImgSize=(400, 400))


    def create_reaction_image_grid(self, selected_reaction_instances):

        img_list = []

        for rxn in selected_reaction_instances:

            reactant_a = rxn.monomer_1.smiles

            if rxn.monomer_2:
                reactant_b = rxn.monomer_2.smiles
            else:
                reactant_b = rxn.monomer_1.smiles

            img = self.create_reaction_image(
                reactant_a,
                reactant_b,
                rxn.reaction_smarts,
                rxn.reaction_name
            )

            if img:
                img_list.append(img)

        if not img_list:
            return None

        single_width, single_height = img_list[0].size
        total_height = single_height * len(img_list)
        total_width = single_width + 80

        new_img = Image.new("RGB", (total_width, total_height), (255, 255, 255))
        draw = ImageDraw.Draw(new_img)

        try:
            font = ImageFont.truetype("arial.ttf", 40)
        except:
            font = ImageFont.load_default()

        y_offset = 0
        for i, img in enumerate(img_list, start=1):
            draw.text((10, y_offset + 10), f"{i}.", fill=(0, 0, 0), font=font)
            new_img.paste(img, (70, y_offset))
            y_offset += single_height

        return new_img

    def reaction_selection(self, reaction_instances: list[ReactionInstance]) -> list[ReactionInstance]:
        print("Detected Reactions:")
        for i, rx in enumerate(reaction_instances, start=1):
            print(f"{i}. {rx.reaction_name}")
            print(f"   Monomer 1: {rx.monomer_1.smiles}")
            print(f"   FG 1: {rx.functional_group_1.fg_name}")

            if rx.monomer_2:
                print(f"   Monomer 2: {rx.monomer_2.smiles}")
                print(f"   FG 2: {rx.functional_group_2.fg_name}")

            print()
        if not reaction_instances:
            print("No reaction instances available.")
            return []

        while True:
            raw = input("Enter reaction numbers (comma-separated): ").strip()

            if not raw:
                print("Empty input. Try again.")
                continue

            tokens = [t.strip() for t in raw.split(",")]

            if not all(t.isdigit() for t in tokens):
                print("Invalid input. Use integers like: 1,2,3")
                continue

            indices = sorted(set(int(t) for t in tokens))

            if any(i < 1 or i > len(reaction_instances) for i in indices):
                print("One or more indices out of range.")
                continue

            break
        selected = [reaction_instances[i-1] for i in indices]
        print(f"Selected {len(selected)} reactions.")
        for rx in selected:
            print(f"{rx.reaction_name}: {rx.monomer_1.smiles} ({rx.functional_group_1.fg_name})", end="")
            if rx.monomer_2:
                print(f" + {rx.monomer_2.smiles} ({rx.functional_group_2.fg_name})")
            else:
                print()
        return selected


if __name__ == "__main__":
    detector = ReactionDetector()
    # Example monomer dictionary for testing
    def convert_dict_to_monomer_roles(monomer_dictionary: dict) -> list[MonomerRole]:
        monomer_roles = []

        for _, entry in monomer_dictionary.items():

            smiles = entry["smiles"]

            functionalities = []

            for key, value in entry.items():
                if not isinstance(key, int):
                    continue

                functionalities.append(
                    FunctionalGroupInfo(
                        functionality_type=value["functionality_type"],
                        fg_name=value["functional_group_name"],
                        fg_smarts_1=value.get("functional_group_smarts_1"),
                        fg_count_1=value.get("functional_count_1", 0),
                        fg_smarts_2=value.get("functional_group_smarts_2"),
                        fg_count_2=value.get("functional_count_2", 0),
                    )
                )

            monomer_roles.append(
                MonomerRole(
                    smiles=smiles,
                    name=smiles,  # or replace with actual name if you have one
                    functionalities=tuple(functionalities),
                )
            )

        return monomer_roles
    
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

    detector = ReactionDetector()

    # Convert dictionary to object model
    monomer_roles = convert_dict_to_monomer_roles(monomer_dictionary)

    # Run detector
    reaction_instances = detector.reaction_detector(monomer_roles)

    import os
    from pathlib import Path
    autoreacter_dir = pathlib.Path(__file__).parent.parent.parent.resolve()
    save_dir = autoreacter_dir / "cache" / "_cache_test" / "reaction_visualization"
    os.makedirs(save_dir, exist_ok=True)
    
    file_path =  save_dir / f"reaction_instances.png"
    img = detector.create_reaction_image_grid(reaction_instances)
    if img:
        img.save(file_path)
        print(f"Reaction grid image saved to: {file_path}")
    else:
        print("No valid reaction images generated.")
    
    # Optional: Allow user to select reactions
    selected_reactions = detector.reaction_selection(reaction_instances)