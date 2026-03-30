
"""
Module for preparing chemical reactions for analysis, including atom mapping between reactants and products,
reaction metadata extraction, and visualization utilities using RDKit.

This module processes reaction SMARTS, applies atom mappings, identifies key reaction features (e.g., first shell,
initiators, byproducts), and generates metadata and visualizations for downstream tasks like reaction prediction
or mechanistic analysis.
"""

# WARNING:
# When modifying this file for dataframe or any other indexing varibales use idx and idxs, do not use index or indicies or similar.

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional
import itertools
import os

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdmolops
from PIL.Image import Image
import pandas as pd

from AutoREACTER.detectors.reaction_detector import ReactionInstance
from AutoREACTER.reaction_preparation.reaction_processor.utils import (
    add_dict_as_new_columns, add_column_safe, compare_set, prepare_paths
)
from AutoREACTER.reaction_preparation.reaction_processor.walker import reaction_atom_walker



class MappingError(Exception):
    """Custom exception raised when atom mapping between reactants and products fails or is inconsistent."""


class SMARTSParsingError(Exception):
    """Custom exception raised when parsing SMILES or reaction SMARTS fails."""


@dataclass(slots=True)
class ReactionMetadata:
    reaction_id: int
    reactant_combined_mol: Chem.Mol
    product_combined_mol: Chem.Mol
    reactant_to_product_mapping: Dict[int, int]
    template_reactant_to_product_mapping: Optional[Dict[int, int]] = None
    edge_atoms: Optional[List[int]] = None

    first_shell: Optional[List[int]] = None
    initiators: Optional[List[int]] = None
    byproduct_indices: Optional[List[int]] = None

    reaction_smarts: Optional[str] = None
    reactant_smarts: Optional[str] = None
    product_smiles: Optional[str] = None
    csv_path: Optional[Path] = None
    mol_3d_path: Optional[Path] = None
    reaction_dataframe: Optional[pd.DataFrame] = None
    delete_atom: bool = True
    delete_atom_idx: Optional[int] = None
    activity_stats: bool = True

class PrepareReactions:
    """
    Class for preparing reaction instances by processing SMILES/SMARTS, performing atom mapping,
    extracting reaction features, and generating metadata and visualizations.
    """

    def __init__(self, cache: Path):
        """
        Initialize the reaction preparer with a cache directory for storing CSV outputs.

        Args:
            cache: Path to the base cache directory for storing processed reaction data.
        """
        self.cache = Path(cache)
        self.csv_cache = prepare_paths(self.cache, "csv_cache")

    def prepare_reactions(self, reaction_instances: list[ReactionInstance]) -> list[ReactionMetadata]:
        # Initial processing: atom mapping and basic dictionary formatting
        reactions_metadata = self.process_reaction_instances(reaction_instances)

        # Detect duplicates after processing all reactions to ensure comprehensive comparison
        reaction_metadata = self._detect_duplicates(reactions_metadata)
        
        # Iterate through each detected reaction to perform template analysis
        for reaction in reactions_metadata:
            if not reaction.activity_stats:
                continue  # Skip reactions marked as duplicates
            combined_reactant_molecule_object = reaction.reactant_combined_mol
            reaction_dataframe = reaction.reaction_dataframe
            csv_save_path = reaction.csv_path

            # Create a mapping dictionary {reactant_idx: product_idx} from the dataframe
            fully_mapped_dict = reaction_dataframe.set_index("reactant_idx")["product_idx"].to_dict()
            
            # Extract the 'first shell' (atoms directly involved in bond changes)
            first_shell = reaction_dataframe["first_shell"].dropna().tolist()
            
            # Perform atom walking to determine the reaction template and boundary (edge) atoms
            template_mapped_dict, edge_atoms = reaction_atom_walker(
                combined_reactant_molecule_object,
                first_shell,
                fully_mapped_dict
            )

            # Update the dataframe with specific template mapping indices
            reaction_dataframe = add_dict_as_new_columns(
                reaction_dataframe,
                template_mapped_dict,
                titles = ["template_reactant_idx", "template_product_idx"]
            )
            
            # Append edge atoms information (crucial for building polymer chains)
            reaction_dataframe = add_column_safe(
                reaction_dataframe,
                edge_atoms,
                "edge_atoms"
            )
            
            # Store the updated dataframe back in the reaction object
            # Using .copy() to avoid SettingWithCopyWarning in pandas
            reaction.reaction_dataframe = reaction_dataframe.copy() 
            
            # Persist the processed reaction data to a CSV file
            reaction_dataframe.to_csv(csv_save_path, index=False)
            reaction.edge_atoms = edge_atoms
            reaction.template_reactant_to_product_mapping = template_mapped_dict

        return reaction_metadata

    def _detect_duplicates(self, reaction_metadata_list: list[ReactionMetadata]) -> list[ReactionMetadata]:
        unique_metadata: list[ReactionMetadata] = []
        for reaction in reaction_metadata_list:
            reactants = reaction.reactant_combined_mol
            products = reaction.product_combined_mol

            if compare_set(unique_metadata, reactants, products):
                unique_metadata.append(reaction)
            else:
                reaction.activity_stats = False  # mark duplicate reactions with False for activity_stats

        return unique_metadata
    
    def process_reaction_instances(self, detected_reactions: list[ReactionInstance]) -> list[ReactionMetadata]:
        csv_cache = self.csv_cache
        reaction_metadata = []

        # Process each reaction in the detected reactions class
        for reaction in detected_reactions:
            rxn_smarts = reaction.reaction_smarts
            reactant_smiles_1 = reaction.monomer_1.smiles
            same_reactants = reaction.same_reactants

            # Handle same vs different reactants
            if same_reactants:
                reactant_smiles_2 = reactant_smiles_1
            else:
                reactant_smiles_2 = reaction.monomer_2.smiles

            delete_atoms = reaction.delete_atom

            # Build reaction and reactants
            rxn = self._build_reaction(rxn_smarts)
            mol_reactant_1, mol_reactant_2 = self._build_reactants(reactant_smiles_1, reactant_smiles_2)
            reaction_tuple = self._build_reaction_tuple(same_reactants, mol_reactant_1, mol_reactant_2)

            # Process the reaction
            reaction_metadata = self.process_reaction_products( 
                                   rxn, 
                                   csv_cache, 
                                   reaction_tuple, 
                                   delete_atoms,
                                   reaction_metadata
                                   )
        
        return reaction_metadata
    
    def process_reaction_products(self, 
                                   rxn: Chem.rdChemReactions.ChemicalReaction, 
                                   csv_cache: Path, 
                                   reaction_tuple: list, 
                                   delete_atoms: bool = True,
                                    reaction_metadata: Optional[list[ReactionMetadata]] = None
                                   ) -> list[ReactionMetadata]:
        
        # Process each reactant pair in the reaction tuple
        for pair in reaction_tuple:
            r1, r2 = Chem.Mol(pair[0]), Chem.Mol(pair[1])
            products = rxn.RunReactants((r1, r2))
            
            if not products:
                continue  # Skip if no products are generated
            
            self._assign_atom_map_numbers(r1, r2) # Assign unique map numbers to reactant atoms to track them through the reaction

            # Process each product set generated by the reaction
            for product_set in products:
                # loop over ALL match combinations

                # IMPORTANT: copy reactants for each combination
                r1_copy = Chem.Mol(r1)
                r2_copy = Chem.Mol(r2)  


                df = pd.DataFrame(columns=["reactant_idx", "product_idx"])

                # Step 2: combine for reactant→product mapping
                reactant_combined = Chem.CombineMols(r1_copy, r2_copy)
                product_combined = Chem.CombineMols(*product_set)

                # Restore Atom Map Numbers from Isotopes in the products
                self._reassign_atom_map_numbers_by_isotope(product_combined)
            
                mapping_dict, df = self._build_atom_index_mapping(reactant_combined, product_combined)  # Build mapping dictionary and dataframe of reactant to product atom indices based on map numbers

                self._reveal_template_map_numbers(product_combined) # Make map numbers visible for debugging and visualization

                # Validate the mapping to ensure consistency before proceeding
                if not self._validate_mapping(df, reactant_combined, product_combined):
                    raise MappingError("Atom mapping validation failed: inconsistent mapping between reactants and products.")
                
                first_shell, initiator_idxs = self._assign_first_shell_and_initiators(reactant_combined, 
                                                                                      product_combined, 
                                                                                      mapping_dict)


                byproduct_reactant_idxs = self._detect_byproducts(product_combined, mapping_dict, delete_atoms)

                # Build dataframe
                df_combined = pd.concat([
                    df,
                    pd.Series(first_shell, name="first_shell"),
                    pd.Series(initiator_idxs, name="initiators"),
                    pd.Series(byproduct_reactant_idxs, name="byproduct_idx")
                ], axis=1).astype(pd.Int64Dtype())

                total_products = len(reaction_metadata) + 1

                self._clear_isotopes(reactant_combined, product_combined) # Clear isotopes to restore normal chemistry before saving

                df_combined.to_csv(csv_cache / f"reaction_{total_products}.csv", index=False) 

                reaction_metadata.append(
                    ReactionMetadata(
                                reaction_id = total_products,
                                reactant_combined_mol = reactant_combined,
                                product_combined_mol = product_combined,
                                reactant_to_product_mapping = mapping_dict,
                                first_shell = first_shell,
                                initiators = initiator_idxs,
                                byproduct_indices = byproduct_reactant_idxs,
                                csv_path = csv_cache / f"reaction_{total_products}.csv",
                                reaction_dataframe = df_combined,
                                delete_atom = delete_atoms,
                                delete_atom_idx = byproduct_reactant_idxs[0] if byproduct_reactant_idxs else None,
                                activity_stats = True
                            )
                    )

        return reaction_metadata
    
    def _assign_first_shell_and_initiators(self, reactant_combined, product_combined, mapping_dict):
        first_shell = []
        initiator_idxs = []
        reversed_mapping_dict = {}
        for react_idx, prod_idx in mapping_dict.items():
            reversed_mapping_dict[prod_idx] = react_idx
        

        for p_atom in product_combined.GetAtoms():
            if p_atom.GetAtomMapNum() < 999:

                p_idx = p_atom.GetIdx()

                if p_idx not in reversed_mapping_dict:
                    raise ValueError(f"Mapping error: product atom {p_idx} not found in mapping_dict")

                r_idx = reversed_mapping_dict[p_idx]

                atom = reactant_combined.GetAtomWithIdx(r_idx)
                atom.SetAtomMapNum(p_atom.GetAtomMapNum())

                first_shell.append(r_idx)

                if p_atom.GetAtomMapNum() in [1, 2]:
                    initiator_idxs.append(r_idx)
        if len(initiator_idxs) != 2:
            raise ValueError(f"Expected 2 initiators, got {len(initiator_idxs)}: {initiator_idxs}") 

        return first_shell, initiator_idxs
    
    def _reassign_atom_map_numbers_by_isotope(self, mol: Chem.Mol) -> None:
        """
        surviving_idx = atom.GetIsotope() inside the atom's "Isotope" property because 
        RDKit's reaction engine wouldn't delete them there. Here, we are opening up that 
        safe box to retrieve the ID.""" 
        for atom in mol.GetAtoms():
            surviving_idx = atom.GetIsotope()
            if surviving_idx != 0:
                atom.SetAtomMapNum(surviving_idx) # Restore map number
                atom.SetIsotope(0)                # Clear isotope to restore normal chemistry

    def _build_atom_index_mapping(self, reactant_combined: Chem.Mol, product_combined: Chem.Mol) -> tuple[dict[int, int], pd.DataFrame]:
        mapping_dict = {}
        # Pre-index product atoms by map number (avoids nested loop)
        product_map = {
            atom.GetAtomMapNum(): atom.GetIdx()
            for atom in product_combined.GetAtoms()
            if atom.GetAtomMapNum() != 0
        }
        rows = []
        for r_atom in reactant_combined.GetAtoms():
            r_map_num = r_atom.GetAtomMapNum()

            if r_map_num != 0 and r_map_num in product_map:
                r_idx = r_atom.GetIdx()
                p_idx = product_map[r_map_num]

                mapping_dict[r_idx] = p_idx
                rows.append({
                    "reactant_idx": r_idx,
                    "product_idx": p_idx
                })

        df = pd.DataFrame(rows)

        return mapping_dict, df
    
    def _clear_isotopes(self, mol_1: Chem.Mol, mol_2: Chem.Mol) -> None:
        '''Clear isotopes from a molecule to restore normal chemistry after using isotopes to store custom IDs'''
        for atom in mol_1.GetAtoms():
            atom.SetIsotope(0)
        for atom in mol_2.GetAtoms():
            atom.SetIsotope(0)

    def _reveal_template_map_numbers(self, mol : Chem.Mol) -> None:
        '''Make map numbers assigned to products of reaction visible for display'''
        for atom in mol.GetAtoms():
            if atom.HasProp('old_mapno'): # RDKit sets certain "magic" properties post-reaction, including map numbers and reactant atom indices
                map_num = atom.GetIntProp('old_mapno')
                atom.SetAtomMapNum(map_num)

    def _assign_atom_map_numbers(self, r1: Chem.Mol, r2: Chem.Mol) -> None:
        for atom in r1.GetAtoms():
            idx = 1001 + atom.GetIdx()
            atom.SetAtomMapNum(idx)
            atom.SetIsotope(idx)  # The isotope is what will survive the reaction!

        for atom in r2.GetAtoms():
            idx = 2001 + atom.GetIdx()
            atom.SetAtomMapNum(idx)
            atom.SetIsotope(idx) # The isotope is what will survive the reaction!

    def _validate_mapping(self, df: pd.DataFrame, reactant: Chem.Mol, product: Chem.Mol) -> bool:
        if df["reactant_idx"].notna().sum() != df["product_idx"].notna().sum():
            return False
        if df["reactant_idx"].notna().sum() != reactant.GetNumAtoms():
            return False
        if df["product_idx"].notna().sum() != product.GetNumAtoms():
            return False
        if not self._is_consecutive(df["reactant_idx"].tolist()):
            return False
        if not self._is_consecutive(df["product_idx"].tolist()):
            return False
        return True
    

    def _detect_byproducts(self, product_combined: Chem.Mol, mapping_dict: dict[int, int], delete_atoms: bool) -> list[int]:
        if not delete_atoms:
            return []

        frags = rdmolops.GetMolFrags(product_combined, asMols=True)
        smallest = min(frags, key=lambda m: m.GetNumAtoms())

        byproduct = []

        # reverse mapping: product_idx → reactant_idx
        reversed_mapping_dict = {}
        for react_idx, prod_idx in mapping_dict.items():
            reversed_mapping_dict[prod_idx] = react_idx

        for atom in smallest.GetAtoms():
            p_idx = atom.GetIdx()
            if p_idx in reversed_mapping_dict:
                byproduct.append(reversed_mapping_dict[p_idx])

        return byproduct

    def _is_consecutive(self, num_list: list[int]) -> bool:
        if not num_list:
            return False

        return (
            len(set(num_list)) == len(num_list)
            and max(num_list) - min(num_list) + 1 == len(num_list)
    )
    
    def _build_reaction_tuple(self, same_reactants, mol_reactant_1, mol_reactant_2):
        if same_reactants:
            return [[mol_reactant_1, mol_reactant_1]]

        return [[mol_reactant_1, mol_reactant_2], [mol_reactant_2, mol_reactant_1]]
    
    
    def _build_reaction(self, rxn_smarts):
        return AllChem.ReactionFromSmarts(rxn_smarts)

    def _build_reactants(self, reactant_smiles_1, reactant_smiles_2):
        mol_reactant_1 = Chem.MolFromSmiles(reactant_smiles_1)
        mol_reactant_1 = Chem.AddHs(mol_reactant_1)  # Add explicit hydrogens
        mol_reactant_2 = Chem.MolFromSmiles(reactant_smiles_2)
        mol_reactant_2 = Chem.AddHs(mol_reactant_2)  # Add explicit hydrogens
        return mol_reactant_1, mol_reactant_2

    def reaction_templates_highlighted_image_grid(
        self,
        metadata_list: List[ReactionMetadata],
        highlight_type: str = "template",
    ) -> Image:
        mols = []
        highlight_lists = []
        highlight_colors = []

        for metadata in metadata_list:
            reactant = Chem.RWMol(metadata.reactant_combined_mol)
            product = Chem.RWMol(metadata.product_combined_mol)

            # Clear atom maps for clean visualization
            for atom in reactant.GetAtoms():
                atom.SetAtomMapNum(0)
            for atom in product.GetAtoms():
                atom.SetAtomMapNum(0)

            df = metadata.reaction_dataframe
            atoms: List[int] = []
            color_map: Dict[int, tuple] = {}

            if highlight_type == "template":
                atoms = list((metadata.template_reactant_to_product_mapping or {}).keys())
                for a in atoms:
                    color_map[a] = (0.2, 0.6, 1.0)  # blue

            elif highlight_type == "edge":
                atoms = metadata.edge_atoms or []
                for a in atoms:
                    color_map[a] = (1.0, 0.4, 0.0)  # orange

            elif highlight_type == "initiators":
                atoms = df["initiators"].dropna().astype(int).tolist() if df is not None else []
                for a in atoms:
                    color_map[a] = (0.0, 0.8, 0.2)  # green

            elif highlight_type == "delete":
                if metadata.delete_atom and metadata.byproduct_indices:
                    atoms = metadata.byproduct_indices
                    for a in atoms:
                        color_map[a] = (1.0, 0.0, 0.0)  # red

            mols.extend([reactant, product])
            highlight_lists.append(atoms)
            highlight_lists.append([])  # no highlights on products
            highlight_colors.append(color_map)
            highlight_colors.append({})

        img = Draw.MolsToGridImage(
            mols,
            molsPerRow=2,
            highlightAtomLists=highlight_lists,
            highlightAtomColors=highlight_colors,
            subImgSize=(400, 400),
            useSVG=False,
        )
        return img