
"""
Module for preparing chemical reactions for analysis, including atom mapping between reactants and products,
reaction metadata extraction, and visualization utilities using RDKit.

This module processes reaction SMARTS, applies atom mappings, identifies key reaction features (e.g., first shell,
initiators, byproducts), and generates metadata and visualizations for downstream analysis. It also includes validation checks to ensure mapping 
consistency and completeness.
"""


# WARNING:
# When modifying this file for dataframe or any other indexing variables use idx and idxs, do not use index or indices or similar.

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

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
    """
    Stores comprehensive metadata for a single reaction including molecular structures, 
    atom mappings, and analysis results.
    reaction_id: Unique identifier for the reaction instance
    reactant_combined_mol: RDKit molecule object representing combined reactants
    product_combined_mol: RDKit molecule object representing combined products
    reactant_to_product_mapping: Dictionary mapping reactant atom indices to product atom indices
    product_to_reactant_mapping: Dictionary mapping product atom indices back to reactant atom indices
    template_reactant_to_product_mapping: Optional dictionary mapping reactant indices in the template to product indices
    edge_atoms: Optional list of reactant atom indices that are at the edge of the local environment
    first_shell: Optional list of reactant atom indices in the first coordination shell (reaction center)
    initiators: Optional list of reactant atom indices that are initiators (map numbers 1 or 2)
    byproduct_indices: Optional list of reactant atom indices corresponding to detected byproducts
    reaction_smarts: Optional string of the reaction SMARTS pattern
    reactant_smarts: Optional string of the combined reactant SMILES
    product_smarts: Optional string of the combined product SMILES
    csv_path: Optional Path to the CSV file containing atom mappings and analysis
    reaction_dataframe: Optional pandas DataFrame containing detailed mapping and analysis results
    delete_atom: Boolean indicating whether the reaction involves a delete atom (byproduct)
    delete_atom_idx: Optional integer index of the reactant atom that corresponds to the byproduct
    reactant_combined_3Dmol_path: Optional Path to the 3D structure file for the combined reactants
    product_combined_3Dmol_path: Optional Path to the 3D structure file for the combined products
    activity_stats: Boolean indicating whether this reaction should be included in activity statistics (e.g., not a duplicate)
    """
    reaction_id: int
    reactant_combined_RDmol: Chem.Mol
    product_combined_RDmol: Chem.Mol
    reactant_to_product_mapping: Dict[int, int]
    product_to_reactant_mapping: Dict[int, int]
    template_reactant_to_product_mapping: Optional[Dict[int, int]] = None
    edge_atoms: Optional[List[int]] = None
    first_shell: Optional[List[int]] = None
    initiators: Optional[List[int]] = None
    byproduct_indices: Optional[List[int]] = None
    reaction_smarts: Optional[str] = None
    reactant_smiles: Optional[str] = None
    product_smiles: Optional[str] = None
    csv_path: Optional[Path] = None
    reaction_dataframe: Optional[pd.DataFrame] = None
    delete_atom: bool = True
    delete_atom_idx: Optional[int] = None
    reactant_combined_3Dmol_path: Optional[Path] = None
    product_combined_3Dmol_path: Optional[Path] = None
    activity_stats: bool = True


class PrepareReactions:
    """Processes chemical reactions: builds atom mappings, identifies reaction centers, and detects byproducts."""

    def __init__(self, cache: Path):
        """Initialize with cache directory for storing reaction CSVs."""
        self.cache = Path(cache)
        self.csv_cache = prepare_paths(self.cache, "csv_cache")

    # --- PUBLIC ---

    def prepare_reactions(self, reaction_instances: list[ReactionInstance]) -> list[ReactionMetadata]:
        """
        Main pipeline: processes reaction instances, detects duplicates, and enriches metadata with template mappings.
        
        Args:
            reaction_instances: List of detected reaction instances to process
            
        Returns:
            List of processed ReactionMetadata objects with template mappings and edge atoms
        """
        reactions_metadata = self._process_reaction_instances(reaction_instances)
        reaction_metadata = self._detect_duplicates(reactions_metadata)
        
        for reaction in reactions_metadata:
            # Skip reactions marked as duplicates
            if not reaction.activity_stats:
                continue  # Skip reactions marked as duplicates
            
            combined_reactant_molecule_object = reaction.reactant_combined_RDmol
            reaction_dataframe = reaction.reaction_dataframe
            csv_save_path = reaction.csv_path
            
            # Build full atom mapping dictionary from dataframe
            fully_mapped_dict = reaction_dataframe.set_index("reactant_idx")["product_idx"].to_dict()
            first_shell = reaction_dataframe["first_shell"].dropna().tolist()
            
            # Generate template mapping by walking reaction graph
            template_mapped_dict, edge_atoms = reaction_atom_walker(
                combined_reactant_molecule_object,
                first_shell,
                fully_mapped_dict
            )
            
            # Add template mapping and edge atoms to dataframe
            reaction_dataframe = add_dict_as_new_columns(
                reaction_dataframe,
                template_mapped_dict,
                titles=["template_reactant_idx", "template_product_idx"]
            )

            # Add edge atoms as a new column in the dataframe
            reaction_dataframe = add_column_safe(
                reaction_dataframe,
                edge_atoms,
                "edge_atoms"
            )
            
            # Save updated dataframe back to CSV and update metadata
            reaction.reaction_dataframe = reaction_dataframe.copy()
            # Save the updated dataframe with template mappings and edge atoms to CSV
            reaction_dataframe.to_csv(csv_save_path, index=False)
            # Update metadata with template mapping and edge atoms
            reaction.edge_atoms = edge_atoms
            # Store the template mapping in the metadata for later use
            reaction.template_reactant_to_product_mapping = template_mapped_dict

        return reaction_metadata
    
    # --- PIPELINE STEPS (PRIVATE) ---
    
    def _process_reaction_instances(self, detected_reactions: list[ReactionInstance]) -> list[ReactionMetadata]:
        """
        Converts ReactionInstance objects into ReactionMetadata by building molecules and running reactions.
        
        Args:
            detected_reactions: List of detected reaction instances
            
        Returns:
            List of ReactionMetadata objects with atom mappings
        """
        csv_cache = self.csv_cache
        reaction_metadata = []

        for reaction in detected_reactions:
            rxn_smarts = reaction.reaction_smarts
            reactant_smiles_1 = reaction.monomer_1.smiles
            same_reactants = reaction.same_reactants
            
            # Handle case where both reactants are identical
            if same_reactants:
                reactant_smiles_2 = reactant_smiles_1
            else:
                reactant_smiles_2 = reaction.monomer_2.smiles

            delete_atoms = reaction.delete_atom

            # Build reaction and reactant molecules
            rxn = self._build_reaction(rxn_smarts)
            # This function also runs the reaction and builds metadata for each product set, including atom mappings and byproduct detection
            mol_reactant_1, mol_reactant_2 = self._build_reactants(reactant_smiles_1, reactant_smiles_2)
            # Build reaction tuple based on whether reactants are the same or different. If same, only one ordering is needed. If different, both orderings are processed to account for reaction directionality.
            reaction_tuple = self._build_reaction_tuple(same_reactants, mol_reactant_1, mol_reactant_2)

            # Process products and build metadata
            reaction_metadata = self._process_reaction_products(
                rxn,
                csv_cache,
                reaction_tuple,
                delete_atoms,
                reaction_metadata
            )
        
        return reaction_metadata

    def _detect_duplicates(self, reaction_metadata_list: list[ReactionMetadata]) -> list[ReactionMetadata]:
        """
        Filters duplicate reactions based on reactant and product molecules.
        
        Args:
            reaction_metadata_list: List of reaction metadata to filter
            
        Returns:
            List of unique reactions; duplicates marked with activity_stats=False
        """
        unique_metadata: list[ReactionMetadata] = []
        
        for reaction in reaction_metadata_list:
            # Compare current reaction's reactants and products against unique reactions collected so far
            reactants = reaction.reactant_combined_RDmol
            products = reaction.product_combined_RDmol

            # Keep reaction if it's unique, otherwise mark as duplicate
            if compare_set(unique_metadata, reactants, products):
                unique_metadata.append(reaction)
            else:
                reaction.activity_stats = False
                
        return unique_metadata
    
    def _process_reaction_products(self, 
                                   rxn: Chem.rdChemReactions.ChemicalReaction, 
                                   csv_cache: Path, 
                                   reaction_tuple: list, 
                                   delete_atoms: bool = True,
                                   reaction_metadata: Optional[list[ReactionMetadata]] = None
                                   ) -> list[ReactionMetadata]:
        """
        Runs reactions on reactant pairs and builds metadata for each product set.
        
        Args:
            rxn: RDKit ChemicalReaction object
            csv_cache: Path to cache directory for saving CSVs
            reaction_tuple: List of reactant pairs to process
            delete_atoms: Whether to detect and track byproducts
            reaction_metadata: Accumulator list for metadata objects
            
        Returns:
            Updated list of ReactionMetadata objects
        """
        if reaction_metadata is None:
            reaction_metadata = []
        
        for pair in reaction_tuple:
            r1, r2 = Chem.Mol(pair[0]), Chem.Mol(pair[1])
            
            # Assign unique map numbers and isotopes to track atoms through reaction
            self._assign_atom_map_numbers_and_set_isotopes(r1, r2)
            # Run the reaction to get products
            products = rxn.RunReactants((r1, r2))
            
            # If no products are generated, skip to the next reactant pair
            if not products:
                continue
            
            # Process each product set generated by the reaction
            for product_set in products:
                df = pd.DataFrame(columns=["reactant_idx", "product_idx"])

                # Combine molecules for mapping
                reactant_combined = Chem.CombineMols(r1, r2)
                product_combined = Chem.CombineMols(*product_set)

                # Restore atom map numbers from isotopes (which survive the reaction)
                self._reassign_atom_map_numbers_by_isotope(product_combined)
                
                # Build bidirectional atom index mappings
                mapping_dict, df = self._build_atom_index_mapping(reactant_combined, product_combined)
                reverse_mapping = {v: k for k, v in mapping_dict.items()}

                # Restore map numbers for visualization
                self._reveal_template_map_numbers(product_combined)

                # Validate mapping consistency
                self._validate_mapping(df, reactant_combined, product_combined)

                # Identify atoms involved in reaction center and initiators
                first_shell, initiator_idxs = self._assign_first_shell_and_initiators(
                    reactant_combined,
                    product_combined,
                    reverse_mapping
                )

                # Detect byproducts (smallest fragments)
                byproduct_reactant_idxs = self._detect_byproducts(product_combined, reverse_mapping, delete_atoms)

                # Combine all mapping data into single dataframe
                df_combined = pd.concat([
                    df,
                    pd.Series(first_shell, name="first_shell"),
                    pd.Series(initiator_idxs, name="initiators"),
                    pd.Series(byproduct_reactant_idxs, name="byproduct_idx")
                ], axis=1).astype(pd.Int64Dtype())

                total_products = len(reaction_metadata) + 1

                # Clear isotopes before saving to restore normal chemistry
                self._clear_isotopes(reactant_combined, product_combined)

                # Save mapping dataframe to CSV
                df_combined.to_csv(csv_cache / f"reaction_{total_products}.csv", index=False)

                # Create and store metadata object
                reaction_metadata.append(
                    ReactionMetadata(
                        reaction_id=total_products,
                        reactant_combined_RDmol=reactant_combined,
                        product_combined_RDmol=product_combined,
                        reactant_to_product_mapping=mapping_dict,
                        product_to_reactant_mapping=reverse_mapping,
                        first_shell=first_shell,
                        initiators=initiator_idxs,
                        byproduct_indices=byproduct_reactant_idxs,
                        csv_path=csv_cache / f"reaction_{total_products}.csv",
                        reaction_dataframe=df_combined,
                        delete_atom=delete_atoms,
                        delete_atom_idx=byproduct_reactant_idxs[0] if byproduct_reactant_idxs else None,
                        activity_stats=True
                    )
                )

        return reaction_metadata
    
    # --- CORE REACTION LOGIC ---
    
    def _assign_first_shell_and_initiators(self, 
                                           reactant_combined: Chem.Mol, 
                                           product_combined: Chem.Mol, 
                                           reversed_mapping_dict: dict[int, int]) -> tuple[list[int], list[int]]:
        """
        Identifies atoms in the first coordination shell (atoms with map numbers < 999) and initiator atoms.
        Initiators are atoms with map numbers 1 or 2 (typically the reactive centers).
        
        Args:
            reactant_combined: Combined reactant molecule
            product_combined: Combined product molecule
            reversed_mapping_dict: Product idx -> Reactant idx mapping
            
        Returns:
            Tuple of (first_shell atom indices, initiator atom indices)
            
        Raises:
            ValueError: If exactly 2 initiators are not found
        """
        first_shell = []
        initiator_idxs = []

        for p_atom in product_combined.GetAtoms():
            # Only process atoms with valid map numbers (< 999 indicates non-byproduct)
            if p_atom.GetAtomMapNum() < 999:
                p_idx = p_atom.GetIdx()

                if p_idx not in reversed_mapping_dict:
                    raise ValueError(f"Mapping error: product atom {p_idx} not found in mapping_dict")

                r_idx = reversed_mapping_dict[p_idx]
                atom = reactant_combined.GetAtomWithIdx(r_idx)
                atom.SetAtomMapNum(p_atom.GetAtomMapNum())

                first_shell.append(r_idx)

                # Initiators are atoms with map numbers 1 or 2
                if p_atom.GetAtomMapNum() in [1, 2]:
                    initiator_idxs.append(r_idx)
        
        if len(initiator_idxs) != 2:
            raise ValueError(f"Expected 2 initiators, got {len(initiator_idxs)}: {initiator_idxs}")

        return first_shell, initiator_idxs
    
    def _detect_byproducts(self, 
                          product_combined: Chem.Mol, 
                          reversed_mapping_dict: dict[int, int], 
                          delete_atoms: bool) -> list[int]:
        """
        Identifies byproduct atoms (smallest molecular fragment) and maps them back to reactant space.
        
        Args:
            product_combined: Combined product molecule
            reversed_mapping_dict: Product idx -> Reactant idx mapping
            delete_atoms: Whether to perform byproduct detection
            
        Returns:
            List of reactant indices corresponding to byproduct atoms
        """
        # If delete_atoms is False, we skip byproduct detection and return an empty list
        if not delete_atoms:
            return []

        # Fragment molecule and find smallest fragment
        frags = rdmolops.GetMolFrags(product_combined, asMols=True)
        smallest = min(frags, key=lambda m: m.GetNumAtoms())

        byproduct_reactant_indices = []

        # Map byproduct atoms back to reactant indices
        for atom in smallest.GetAtoms():
            p_idx = atom.GetIdx()
            if p_idx in reversed_mapping_dict:
                byproduct_reactant_indices.append(reversed_mapping_dict[p_idx])

        return byproduct_reactant_indices
    
    def _validate_mapping(self, df: pd.DataFrame, reactant: Chem.Mol, product: Chem.Mol) -> None:
        """
        Validates atom mapping consistency: checks for required columns, duplicates, bounds, and completeness.
        
        Args:
            df: Dataframe containing reactant_idx and product_idx columns
            reactant: Reactant molecule
            product: Product molecule
            
        Raises:
            MappingError: If any validation check fails
        """
        # Ensure dataframe exists and has required columns
        if df is None or df.empty:
            raise MappingError("Mapping validation failed: empty dataframe")

        # Check for required columns
        required_cols = {"reactant_idx", "product_idx"}
        if not required_cols.issubset(df.columns):
            raise MappingError(f"Mapping validation error: required columns {required_cols} not found in dataframe.")

        # Extract indices and perform validation checks
        r_idxs = df["reactant_idx"].dropna().tolist()
        p_idxs = df["product_idx"].dropna().tolist()

        # Atom counts must match
        if len(r_idxs) != len(p_idxs):
            raise MappingError(f"Mapping validation error: mismatch in atom counts between reactant and product.")

        # No duplicate mappings (1-to-1 mapping required)
        if len(set(r_idxs)) != len(r_idxs):
            raise MappingError(f"Mapping validation error: duplicate indices found in reactant mapping.")
        if len(set(p_idxs)) != len(p_idxs):
            raise MappingError(f"Mapping validation error: duplicate indices found in product mapping.")

        # Indices must be within molecule bounds
        if any(idx >= reactant.GetNumAtoms() for idx in r_idxs):
            raise MappingError(f"Mapping validation error: reactant index out of bounds.")
        if any(idx >= product.GetNumAtoms() for idx in p_idxs):
            raise MappingError(f"Mapping validation error: product index out of bounds.")

        # All atoms must be mapped (complete mapping)
        if len(r_idxs) != reactant.GetNumAtoms():
            raise MappingError(f"Mapping validation error: incomplete mapping for reactant.")
        if len(p_idxs) != product.GetNumAtoms():
            raise MappingError(f"Mapping validation error: incomplete mapping for product.")

    # --- ATOM MAPPING ---
    
    def _assign_atom_map_numbers_and_set_isotopes(self, r1: Chem.Mol, r2: Chem.Mol) -> None:
        """
        Assigns unique map numbers and isotopes to reactant atoms for tracking through reaction.
        Isotopes survive RDKit's reaction engine, allowing atom identity recovery post-reaction.
        
        Args:
            r1: First reactant molecule
            r2: Second reactant molecule
        """
        # Assign map numbers 1001+ to first reactant atoms
        for atom in r1.GetAtoms():
            idx = 1001 + atom.GetIdx()
            atom.SetAtomMapNum(idx)
            atom.SetIsotope(idx)  # Isotope survives the reaction

        # Assign map numbers 2001+ to second reactant atoms
        for atom in r2.GetAtoms():
            idx = 2001 + atom.GetIdx()
            atom.SetAtomMapNum(idx)
            atom.SetIsotope(idx) # Isotope survives the reaction
			
    def _reassign_atom_map_numbers_by_isotope(self, mol: Chem.Mol) -> None:
        """
        Restores atom map numbers from isotope values after reaction.
        RDKit's reaction engine preserves isotopes, allowing recovery of original atom identities.
        
        Args:
            mol: Product molecule with isotope information
        """
        for atom in mol.GetAtoms():
            surviving_idx = atom.GetIsotope()
            if surviving_idx != 0:
                atom.SetAtomMapNum(surviving_idx)  # Restore original map number
                atom.SetIsotope(0)                 # Clear isotope to restore normal chemistry

    def _build_atom_index_mapping(self, 
                                  reactant_combined: Chem.Mol, 
                                  product_combined: Chem.Mol) -> tuple[dict[int, int], pd.DataFrame]:
        """
        Builds bidirectional atom index mapping between reactants and products using map numbers.
        
        Args:
            reactant_combined: Combined reactant molecule
            product_combined: Combined product molecule
            
        Returns:
            Tuple of (mapping dict: reactant_idx -> product_idx, dataframe with mapping)
        """
        mapping_dict = {}
        
        # Pre-index product atoms by map number for O(1) lookup
        product_map = {
            atom.GetAtomMapNum(): atom.GetIdx()
            for atom in product_combined.GetAtoms()
            if atom.GetAtomMapNum() != 0
        }
        
        rows = []
        for r_atom in reactant_combined.GetAtoms():
            r_map_num = r_atom.GetAtomMapNum()

            # Match reactant atom to product atom via map number
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

    def _reveal_template_map_numbers(self, mol: Chem.Mol) -> None:
        """
        Restores map numbers from RDKit's internal 'old_mapno' property for visualization.
        RDKit stores original map numbers in this property after reaction execution.
        
        Args:
            mol: Product molecule
        """
        for atom in mol.GetAtoms():
            if atom.HasProp('old_mapno'):
                map_num = atom.GetIntProp('old_mapno')
                atom.SetAtomMapNum(map_num)

    def _clear_isotopes(self, mol_1: Chem.Mol, mol_2: Chem.Mol) -> None:
        """
        Clears isotope values from molecules to restore normal chemistry after using isotopes for atom tracking.
        
        Args:
            mol_1: First molecule to clear
            mol_2: Second molecule to clear
        """
        for atom in mol_1.GetAtoms():
            atom.SetIsotope(0)
        for atom in mol_2.GetAtoms():
            atom.SetIsotope(0)

    # --- BUILDERS ---
    
    def _build_reaction(self, rxn_smarts: str) -> Chem.rdChemReactions.ChemicalReaction:
        """Builds RDKit ChemicalReaction object from SMARTS string."""
        return AllChem.ReactionFromSmarts(rxn_smarts)
    
    def _build_reactants(self, reactant_smiles_1: str, reactant_smiles_2: str) -> tuple[Chem.Mol, Chem.Mol]:
        """
        Builds reactant molecules from SMILES strings with explicit hydrogens added.
        
        Args:
            reactant_smiles_1: SMILES string for first reactant
            reactant_smiles_2: SMILES string for second reactant
            
        Returns:
            Tuple of (reactant1 molecule, reactant2 molecule) with explicit hydrogens
        """
        mol_reactant_1 = Chem.MolFromSmiles(reactant_smiles_1)
        if mol_reactant_1 is None:
            raise SMARTSParsingError(f"Failed to parse first reactant SMILES: {reactant_smiles_1!r}")
        mol_reactant_1 = Chem.AddHs(mol_reactant_1)
        
        mol_reactant_2 = Chem.MolFromSmiles(reactant_smiles_2)
        if mol_reactant_2 is None:
            raise SMARTSParsingError(f"Failed to parse second reactant SMILES: {reactant_smiles_2!r}")
        mol_reactant_2 = Chem.AddHs(mol_reactant_2)
        
        return mol_reactant_1, mol_reactant_2
    
    def _build_reaction_tuple(self, same_reactants: bool, mol_reactant_1: Chem.Mol, mol_reactant_2: Chem.Mol) -> list:
        """
        Builds list of reactant pairs to process. If reactants are identical, returns single pair.
        Otherwise returns both orderings to account for reaction directionality.
        
        Args:
            same_reactants: Whether both reactants are identical
            mol_reactant_1: First reactant molecule
            mol_reactant_2: Second reactant molecule
            
        Returns:
            List of reactant pairs [[r1, r2], ...] to process
        """
        if same_reactants:
            return [[mol_reactant_1, mol_reactant_1]]

        return [[mol_reactant_1, mol_reactant_2], [mol_reactant_2, mol_reactant_1]]
    
    # --- HELPERS ---
    
    def _is_consecutive(self, num_list: list[int]) -> bool:
        """
        Checks if list contains consecutive integers with no duplicates.
        
        Args:
            num_list: List of integers to check
            
        Returns:
            True if list is consecutive and has no duplicates, False otherwise
        """
        if not num_list:
            return False

        return (
            len(set(num_list)) == len(num_list)
            and max(num_list) - min(num_list) + 1 == len(num_list)
        )

    # --- VISUALIZATION ---
    
    def reaction_templates_highlighted_image_grid(
        self,
        metadata_list: List[ReactionMetadata],
        highlight_type: str = "template",
    ) -> Image:
        """
        Generates grid image of reactions with highlighted atoms based on type.
        
        Args:
            metadata_list: List of reaction metadata to visualize
            highlight_type: Type of atoms to highlight - "template", "edge", "initiators", or "delete"
            
        Returns:
            PIL Image containing 2-column grid of reactant-product pairs with highlighted atoms
        """
        mols = []
        highlight_lists = []
        highlight_colors = []
        names = []

        for metadata in metadata_list:
            reactant = Chem.RWMol(metadata.reactant_combined_RDmol)
            product = Chem.RWMol(metadata.product_combined_RDmol)
            names.extend([f"pre_{metadata.reaction_id}", f"post_{metadata.reaction_id}"])

            # Clear atom maps for clean visualization
            for atom in reactant.GetAtoms():
                atom.SetAtomMapNum(0)
            for atom in product.GetAtoms():
                atom.SetAtomMapNum(0)

            df = metadata.reaction_dataframe
            atoms: List[int] = []
            color_map: Dict[int, tuple] = {}

            # Select atoms to highlight based on type
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
            
            # Build atom mappings for highlighting
            forward_map = metadata.reactant_to_product_mapping
            reactant_atoms = atoms

            # Map reactant atoms to product atoms
            product_atoms = []
            for r_idx in reactant_atoms:
                if r_idx in forward_map:
                    product_atoms.append(forward_map[r_idx])

            # Build color maps for reactant and product
            reactant_color_map = {a: color_map[a] for a in reactant_atoms}
            product_color_map = {p: color_map[r] for r, p in forward_map.items() if r in reactant_atoms}

            highlight_lists.append(reactant_atoms)
            highlight_lists.append(product_atoms)
            highlight_colors.append(reactant_color_map)
            highlight_colors.append(product_color_map)

        # Generate grid image with 2 molecules per row
        img = Draw.MolsToGridImage(
            mols,
            legends=names,
            molsPerRow=2,
            highlightAtomLists=highlight_lists,
            highlightAtomColors=highlight_colors,
            subImgSize=(400, 400),
            useSVG=False,
        )
        return img
