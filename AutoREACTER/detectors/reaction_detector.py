"""
Polymerization Reaction Detector and Visualizer.

This module provides tools to identify potential polymerization reactions between 
monomers based on their functional groups using SMARTS patterns and RDKit. 
It includes functionality to detect homo-polymerization and co-polymerization, 
and generates visual grids of the resulting chemical reactions.
"""

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
import os
from typing import Dict, Any, Tuple, Optional, List, Set
from dataclasses import dataclass
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import Draw, rdChemReactions

# Attempt to import internal library components
try:
    from reactions_library import ReactionLibrary
except (ImportError, ModuleNotFoundError):
    from .reactions_library import ReactionLibrary

try:
    from .functional_groups_detector import FunctionalGroupInfo, MonomerRole
except (ImportError, ModuleNotFoundError):
    from functional_groups_detector import FunctionalGroupInfo, MonomerRole

class SMARTSerror(Exception):
    """Error from developer-defined SMARTS patterns not producing expected products.
    Please contact the developer to resolve this issue. 
    Raise a GitHub issue if you encounter this error.
    at <https://github.com/NanoCIPHER-Lab/AutoREACTER/issues>"""
    pass

class EmptyReactionListError(Exception):
    """Raised when no reaction instances are detected, but the user attempts to select reactions.
    This should be prevented by the reaction_selection method, but this error serves as a safeguard."""
    pass

@dataclass(slots=True, frozen=True)
class FunctionalGroupInfo:
    """
    Stores metadata about a specific functional group detected in a molecule.
    
    Attributes:
        functionality_type: Category of the group (e.g., 'vinyl', 'diol').
        fg_name: Human-readable name of the functional group.
        fg_smarts_1: Primary SMARTS pattern for detection.
        fg_count_1: Number of occurrences of the primary pattern.
        fg_smarts_2: Secondary SMARTS pattern (optional).
        fg_count_2: Number of occurrences of the secondary pattern (optional).
    """
    functionality_type: str
    fg_name: str
    fg_smarts_1: str
    fg_count_1: int
    fg_smarts_2: Optional[str] = None
    fg_count_2: Optional[int] = None

@dataclass(slots=True, frozen=True)
class MonomerRole:
    """
    Represents a monomer and its associated functional groups.
    
    Attributes:
        smiles: SMILES string of the monomer.
        name: Identifier or name of the monomer.
        functionalities: A tuple of FunctionalGroupInfo objects present in this monomer.
    """
    smiles: str
    name: str
    functionalities: Tuple[FunctionalGroupInfo, ...]

@dataclass(slots=True, frozen=True)
class ReactionTemplate:
    """
    Defines a generic template for a polymerization reaction.
    """
    reaction_name: str
    reactant_1: str                 # functional group name
    reactant_2: Optional[str]       # functional group name or None
    reaction_smarts: str
    same_reactants: bool
    delete_atom: bool
    references: dict

@dataclass(slots=True, frozen=True)
class ReactionInstance:
    """
    Represents a specific instance of a reaction between identified monomers.
    
    Attributes:
        reaction_name: Name of the polymerization type.
        reaction_smarts: The SMARTS string used to perform the reaction.
        delete_atom: Flag indicating if specific atoms should be removed post-reaction.
        references: Dictionary containing literature references.
        monomer_1: The first monomer involved.
        functional_group_1: The specific FG on monomer 1 participating.
        monomer_2: The second monomer (None for homo-polymerization).
        functional_group_2: The specific FG on monomer 2 participating.
    """
    reaction_name: str
    reaction_smarts: str
    delete_atom: bool
    references: dict
    monomer_1: MonomerRole
    functional_group_1: FunctionalGroupInfo
    monomer_2: Optional[MonomerRole] = None
    functional_group_2: Optional[FunctionalGroupInfo] = None

class ReactionDetector:
    """
    Detects valid chemical reactions between a list of monomers based on a 
    predefined library of polymerization mechanisms.
    """
    def __init__(self):
        """Initializes the detector by loading the reaction library."""
        self.reactions = ReactionLibrary().reactions
    
    def _matching_fgs(self, monomer_entry: MonomerRole, target_group_name: str) -> List[FunctionalGroupInfo]:
        """
        Filters functional groups in a monomer that match a target name.
        
        Args:
            monomer_entry: The monomer to search.
            target_group_name: The name of the FG to look for.
            
        Returns:
            A list of matching FunctionalGroupInfo objects.
        """
        return [fg for fg in monomer_entry.functionalities if fg.fg_name == target_group_name]

    def _seen_pair_key(
            self,
            reaction_name: str,
            monomer_role_1: MonomerRole,
            fg_1: FunctionalGroupInfo,
            monomer_role_2: Optional[MonomerRole] = None,
            fg_2: Optional[FunctionalGroupInfo] = None,
        ) -> Tuple:
        """
        Generates a unique, order-independent key for a reaction pair to avoid duplicates.
        
        Args:
            reaction_name: Name of the reaction.
            monomer_role_1: First monomer.
            fg_1: First functional group.
            monomer_role_2: Second monomer (optional).
            fg_2: Second functional group (optional).
            
        Returns:
            A tuple representing the unique state of this reaction instance.
        """
        if monomer_role_2 is None or fg_2 is None:
            return (reaction_name, monomer_role_1, fg_1)

        # Sort reactants by memory ID to ensure (A, B) and (B, A) produce the same key
        ordered = tuple(sorted(
            [(monomer_role_1, fg_1), (monomer_role_2, fg_2)],
            key=lambda x: (id(x[0]), id(x[1]))
        ))
        return (reaction_name, ordered[0], ordered[1])

    def reaction_detector(self, monomer_roles: List[MonomerRole]) -> List[ReactionInstance]:
        """
        Scans a list of monomers to find all possible polymerization reactions.
        
        Args:
            monomer_roles: List of monomers with their detected functional groups.
            
        Returns:
            A list of ReactionInstance objects representing valid matches.
        """
        reaction_instances = []
        seen_pairs: Set[Tuple] = set()

        for reaction_name, reaction_info in self.reactions.items():
            reactant_1_name = reaction_info.get("reactant_1")
            reactant_2_name = reaction_info.get("reactant_2")
            same_reactants = reaction_info.get("same_reactants", False)
            
            # CASE 1: HOMO-POLYMERIZATION (e.g., A + A)
            if same_reactants and reactant_2_name is None:
                for monomer_role in monomer_roles: 
                    fg_hits = self._matching_fgs(monomer_role, reactant_1_name)
                    for fg in fg_hits:
                        pair_key = self._seen_pair_key(reaction_name, monomer_role, fg)
                        if pair_key not in seen_pairs:
                            seen_pairs.add(pair_key)
                            reaction_instances.append(
                                ReactionInstance(
                                    reaction_name=reaction_name,
                                    reaction_smarts=reaction_info["reaction"],
                                    delete_atom=reaction_info["delete_atom"],
                                    references=reaction_info["reference"],
                                    monomer_1=monomer_role,
                                    functional_group_1=fg
                                )
                            )
            # CASE 2: CO-POLYMERIZATION (e.g., A + B)
            else:
                for monomer_role_i in monomer_roles:
                    fg_hits_i = self._matching_fgs(monomer_role_i, reactant_1_name)
                    if not fg_hits_i: continue
                    
                    for fg_i in fg_hits_i:
                        for monomer_role_j in monomer_roles:
                            # Prevent a monomer reacting with itself in a co-monomer definition
                            if monomer_role_i == monomer_role_j:
                                continue 
                            
                            fg_hits_j = self._matching_fgs(monomer_role_j, reactant_2_name)
                            for fg_j in fg_hits_j:
                                pair_key = self._seen_pair_key(
                                    reaction_name, monomer_role_i, fg_i, monomer_role_j, fg_j
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

    def create_reaction_image(self, reactant_a_smiles: str, reactant_b_smiles: str, 
                              reaction_smarts: str, reaction_name: str) -> Image.Image:
        """
        Simulates a reaction using RDKit and generates a 2D image of the transformation.
        
        Args:
            reactant_a_smiles: SMILES for the first reactant.
            reactant_b_smiles: SMILES for the second reactant.
            reaction_smarts: SMARTS pattern defining the reaction.
            reaction_name: Name of the reaction for metadata.
            
        Returns:
            A PIL Image object of the reaction.
            
        Raises:
            SMARTSerror: If the reaction SMARTS fails to produce products.
        """
        reactant_a = Chem.AddHs(Chem.MolFromSmiles(reactant_a_smiles))
        reactant_b = Chem.AddHs(Chem.MolFromSmiles(reactant_b_smiles))
        
        rxn_engine = rdChemReactions.ReactionFromSmarts(reaction_smarts)
        
        # Try reacting A + B
        products_sets = rxn_engine.RunReactants((reactant_a, reactant_b))

        # Fallback: Try reacting B + A if the first order failed
        if not products_sets:
            products_sets = rxn_engine.RunReactants((reactant_b, reactant_a))
            if not products_sets:
                raise SMARTSerror(f"Reaction SMARTS '{reaction_smarts}' failed for {reactant_a_smiles} + {reactant_b_smiles}")

        # Prepare the reaction for visualization
        display_rxn = rdChemReactions.ChemicalReaction()
        display_rxn.AddReactantTemplate(reactant_a)
        display_rxn.AddReactantTemplate(reactant_b)

        # Sanitize and add the first set of products found
        for prod in products_sets[0]:
            try:
                Chem.SanitizeMol(prod)
            except:
                pass
            display_rxn.AddProductTemplate(prod)

        return Draw.ReactionToImage(display_rxn, subImgSize=(400, 400))

    def create_reaction_image_grid(self, selected_reaction_instances: List[ReactionInstance]) -> Optional[Image.Image]:
        """
        Creates a vertical grid image containing all selected reaction visualizations.
        
        Args:
            selected_reaction_instances: List of reactions to visualize.
            
        Returns:
            A combined PIL Image or None if no images were generated.
        """
        img_list = []

        for rxn in selected_reaction_instances:
            reactant_a = rxn.monomer_1.smiles
            reactant_b = rxn.monomer_2.smiles if rxn.monomer_2 else rxn.monomer_1.smiles

            try:
                img = self.create_reaction_image(
                    reactant_a, reactant_b, rxn.reaction_smarts, rxn.reaction_name
                )
                if img: img_list.append(img)
            except SMARTSerror as e:
                print(f"Skipping visualization: {e}")

        if not img_list:
            return None

        # Calculate dimensions for the vertical stack
        single_width, single_height = img_list[0].size
        total_height = single_height * len(img_list)
        total_width = single_width + 80 # Extra space for numbering

        new_img = Image.new("RGB", (total_width, total_height), (255, 255, 255))
        draw = ImageDraw.Draw(new_img)

        try:
            font = ImageFont.truetype("arial.ttf", 40)
        except:
            font = ImageFont.load_default()

        # Paste each reaction image into the grid
        y_offset = 0
        for i, img in enumerate(img_list, start=1):
            draw.text((10, y_offset + 10), f"{i}.", fill=(0, 0, 0), font=font)
            new_img.paste(img, (70, y_offset))
            y_offset += single_height

        return new_img

    def reaction_selection(self, reaction_instances: List[ReactionInstance]) -> List[ReactionInstance]:
        """
        Interactive CLI to allow users to select specific reactions from the detected list.
        
        Args:
            reaction_instances: List of all detected reactions.
            
        Returns:
            A list of user-selected ReactionInstance objects.
        """
        if not reaction_instances:
            raise EmptyReactionListError("No reaction instances detected. Cannot proceed with selection.")
        
        if len(reaction_instances) == 1:
            print("Only one reaction detected. Automatically selecting it.")
            return reaction_instances

        print("Detected Reactions:")
        for i, rx in enumerate(reaction_instances, start=1):
            print(f"{i}. {rx.reaction_name}")
            print(f"   Monomer 1: {rx.monomer_1.smiles} ({rx.functional_group_1.fg_name})")
            if rx.monomer_2:
                print(f"   Monomer 2: {rx.monomer_2.smiles} ({rx.functional_group_2.fg_name})")
            print()

        while True:
            raw = input("Enter reaction numbers (comma-separated, e.g., 1,3): ").strip()
            if not raw:
                continue

            try:
                tokens = [int(t.strip()) for t in raw.split(",")]
                if any(i < 1 or i > len(reaction_instances) for i in tokens):
                    print("Out of range. Try again.")
                    continue
                indices = sorted(set(tokens))
                break
            except ValueError:
                print("Invalid input. Please use numbers.")

        selected = [reaction_instances[i-1] for i in indices]
        print(f"Selected {len(selected)} reactions.")
        return selected


if __name__ == "__main__":
    # --- Test Setup ---
    def convert_dict_to_monomer_roles(monomer_dictionary: dict) -> List[MonomerRole]:
        """Helper to convert raw dictionary data into MonomerRole objects."""
        monomer_roles = []
        for _, entry in monomer_dictionary.items():
            smiles = entry["smiles"]
            functionalities = []
            for key, value in entry.items():
                if not isinstance(key, int): continue
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
            monomer_roles.append(MonomerRole(smiles=smiles, name=smiles, functionalities=tuple(functionalities)))
        return monomer_roles
    
    # Mock monomer data for testing detection logic
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
    monomer_roles = convert_dict_to_monomer_roles(monomer_dictionary)
    reaction_instances = detector.reaction_detector(monomer_roles)

    # --- Save Results ---
    autoreacter_dir = pathlib.Path(__file__).parent.parent.parent.resolve()
    save_dir = autoreacter_dir / "cache" / "_cache_test" / "reaction_visualization"
    os.makedirs(save_dir, exist_ok=True)
    
    file_path = save_dir / "reaction_instances.png"
    img = detector.create_reaction_image_grid(reaction_instances)
    
    if img:
        img.save(file_path)
        print(f"Reaction grid image saved to: {file_path}")
    
    # Optional interactive selection
    selected_reactions = detector.reaction_selection(reaction_instances)

    print(f"Total detected reactions: {len(selected_reactions)}")
    for i, rx in enumerate(selected_reactions, start=1):
        print(f"{i}. {rx.reaction_name} - {rx.monomer_1.smiles} ({rx.functional_group_1.fg_name})", end="")
        if rx.monomer_2:
            print(f" + {rx.monomer_2.smiles} ({rx.functional_group_2.fg_name})")
        else:
            print()
    print(selected_reactions)