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
- CF2=CF2 (Tetrafluoroethylene) Polymerization : Seprate case from typical CH2=CHY vinyl polymerization due to unique reactivity and applications. 
"""

import pathlib
import os
from typing import Dict, Any, Tuple, Optional, List, Set
from dataclasses import dataclass
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import Draw, rdChemReactions
from AutoREACTER.detectors.reactions_library import ReactionLibrary
from AutoREACTER.detectors.functional_groups_detector import FunctionalGroupInfo, MonomerRole

class InfiniteReactionLoopError(Exception):
    """Raised when the reaction pooling process exceeds a reasonable number of iterations, indicating a potential infinite loop."""
    pass

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


@dataclass(slots=True)
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
    same_reactants: bool

    monomer_1: MonomerRole
    functional_group_1: FunctionalGroupInfo
    monomer_1_indexes: Tuple[int, ...] = ()

    monomer_2: Optional[MonomerRole] = None
    functional_group_2: Optional[FunctionalGroupInfo] = None
    monomer_2_indexes: Tuple[int, ...] = ()


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
        Creates a unique key for a given reaction instance to track seen pairs and avoid duplicates.
        For co-polymerization, the order of monomers does not matter (A + B is the same as B + A).
        Args:
            reaction_name: Name of the reaction type.
            monomer_role_1: The first monomer role.
            fg_1: The functional group on the first monomer.
            monomer_role_2: The second monomer role (optional).
            fg_2: The functional group on the second monomer (optional).
        Returns:
            A tuple representing the unique key for the reaction instance.
    """

        if monomer_role_2 is None or fg_2 is None:
            return (
                reaction_name,
                monomer_role_1.smiles,
                fg_1.fg_name
            )

        pair1 = (monomer_role_1.smiles, fg_1.fg_name)
        pair2 = (monomer_role_2.smiles, fg_2.fg_name)

        ordered = tuple(sorted([pair1, pair2]))

        return (reaction_name, ordered)

    def reaction_detector(self, monomer_roles: List[MonomerRole]) -> List[ReactionInstance]:
        """
        Scans a list of monomers to find all possible polymerization reactions.
        
        Args:
            monomer_roles: List of monomers with their detected functional groups.
            reactants_pool: Pool of reactants for possible reactions.

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
                                        same_reactants=same_reactants,
                                        monomer_1=monomer_role,
                                        functional_group_1=fg,
                                        monomer_1_indexes=tuple(monomer_role.indexes),
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
                                            same_reactants=same_reactants,
                                            monomer_1=monomer_role_i,
                                            functional_group_1=fg_i,
                                            monomer_1_indexes=tuple(monomer_role_i.indexes),
                                            monomer_2=monomer_role_j,
                                            functional_group_2=fg_j,
                                            monomer_2_indexes=tuple(monomer_role_j.indexes),
                                            )
                                        )
                    
            # case 1.1: If same reacants has two functional groups (e.g., A + A with FG1 and FG2)
            if not same_reactants and reactant_2_name is not None:
                for monomer_role in monomer_roles:

                    fg_hits_1 = self._matching_fgs(monomer_role, reactant_1_name)
                    fg_hits_2 = self._matching_fgs(monomer_role, reactant_2_name)

                    for fg_1 in fg_hits_1:
                        for fg_2 in fg_hits_2:

                            # Skip identical FG objects
                            if fg_1 == fg_2:
                                continue

                            pair_key = self._seen_pair_key(
                                reaction_name, monomer_role, fg_1, monomer_role, fg_2
                            )

                            if pair_key not in seen_pairs:
                                seen_pairs.add(pair_key)

                                reaction_instances.append(
                                    ReactionInstance(
                                        reaction_name=reaction_name,
                                        reaction_smarts=reaction_info["reaction"],
                                        delete_atom=reaction_info["delete_atom"],
                                        references=reaction_info["reference"],
                                        same_reactants=same_reactants,
                                        monomer_1=monomer_role,
                                        functional_group_1=fg_1,
                                        monomer_1_indexes=tuple(monomer_role.indexes),
                                        monomer_2=monomer_role,
                                        functional_group_2=fg_2,
                                        monomer_2_indexes=tuple(monomer_role.indexes),
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

    def available_reaction_image_grid(self, selected_reaction_instances: List[ReactionInstance]) -> Optional[Image.Image]:
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


