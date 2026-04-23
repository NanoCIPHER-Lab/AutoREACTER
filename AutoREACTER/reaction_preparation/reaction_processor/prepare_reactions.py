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
from rdkit.Chem import Draw
from rdkit.Chem.rdChemReactions import ChemicalReaction
from PIL.Image import Image
import pandas as pd

from AutoREACTER.detectors.reaction_detector import ReactionInstance
from AutoREACTER.reaction_preparation.reaction_processor.utils import (
    add_dict_as_new_columns, add_column_safe, prepare_paths
)
from AutoREACTER.reaction_preparation.reaction_processor.walker import reaction_atom_walker
from AutoREACTER.reaction_preparation.reaction_processor.reaction_processing_support import (
    assign_first_shell_and_initiators,
    assign_atom_map_numbers_and_set_isotopes,
    build_atom_index_mapping,
    build_reaction,
    build_reactants,
    build_reaction_tuple,
    clear_isotopes,
    detect_byproducts,
    detect_duplicates,
    reassign_atom_map_numbers_by_isotope,
    reveal_template_map_numbers,
    validate_mapping,
)

from AutoREACTER.reaction_preparation.reaction_processor.reaction_propagation import ReactionPropagation
from typing import TYPE_CHECKING, Dict, List, Optional

if TYPE_CHECKING:
    from AutoREACTER.detectors.functional_groups_detector import MonomerRole

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
        self.reaction_propagation = ReactionPropagation(self.csv_cache) 

    # --- PUBLIC ---

    def prepare_reactions(
        self, 
        reaction_instances: list[ReactionInstance], 
        monomer_roles: list['MonomerRole'] = None
    ) -> list[ReactionMetadata]:
        """
        Main pipeline: processes reaction instances, detects duplicates, and enriches metadata with template mappings.
        
        Args:
            reaction_instances: List of detected reaction instances to process
            
        Returns:
            List of processed ReactionMetadata objects with template mappings and edge atoms
        """
        reactions_metadata = self._process_reaction_instances(reaction_instances)
        reaction_metadata = detect_duplicates(reactions_metadata)
        
        for reaction in reaction_metadata:
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
        
        if monomer_roles is not None:
            reaction_metadata = self.reaction_propagation.run_propagation_loop(
                reaction_metadata, 
                monomer_roles
            )
            
        return reaction_metadata

    # --- PRIVATE ---


    def _process_reaction_instances(self, detected_reactions: list[ReactionInstance]) -> list[ReactionMetadata]:
        """
        Converts ReactionInstance objects into ReactionMetadata by building molecules and running reactions.
        
        Args:
            detected_reactions: List of detected reaction instances
            
        Returns:
            List of ReactionMetadata objects with atom mappings
        """
        csv_cache = self.csv_cache  # Ensure csv_cache is available for processing
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
            rxn = build_reaction(rxn_smarts)
            # This function also runs the reaction and builds metadata for each product set, including atom mappings and byproduct detection
            mol_reactant_1, mol_reactant_2 = build_reactants(reactant_smiles_1, reactant_smiles_2)
            # Build reaction tuple based on whether reactants are the same or different. If same, only one ordering is needed. If different, both orderings are processed to account for reaction directionality.
            reaction_tuple = build_reaction_tuple(same_reactants, mol_reactant_1, mol_reactant_2)

            # Process products and build metadata
            reaction_metadata = self._process_reaction_products(
                rxn,
                csv_cache,
                reaction_tuple,
                delete_atoms,
                reaction_metadata
            )
        
        return reaction_metadata
    
    def _process_reaction_products(self, 
                                   rxn: ChemicalReaction, 
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
            assign_atom_map_numbers_and_set_isotopes(r1, r2)
            # Run the reaction to get products
            products = rxn.RunReactants((r1, r2))
            
            # If no products are generated, skip to the next reactant pair
            if not products:
                continue
            
            # Process each product set generated by the reaction
            for product_set in products:
                df = pd.DataFrame(columns=["reactant_idx", "product_idx"])

                # Combine molecules for mapping safely
                reactant_combined = Chem.CombineMols(r1, r2)

                if len(product_set) == 1:
                    product_combined = product_set[0]
                else:
                    product_combined = product_set[0]
                    for p in product_set[1:]:
                        product_combined = Chem.CombineMols(product_combined, p)

                # Restore atom map numbers from isotopes (which survive the reaction)
                reassign_atom_map_numbers_by_isotope(product_combined)
                
                # Build bidirectional atom index mappings
                mapping_dict, df = build_atom_index_mapping(reactant_combined, product_combined)
                reverse_mapping = {v: k for k, v in mapping_dict.items()}

                # Restore map numbers for visualization
                reveal_template_map_numbers(product_combined)

                # Validate mapping consistency
                validate_mapping(df, reactant_combined, product_combined)

                # Identify atoms involved in reaction center and initiators
                first_shell, initiator_idxs = assign_first_shell_and_initiators(
                    reactant_combined,
                    product_combined,
                    reverse_mapping
                )

                # Detect byproducts (smallest fragments)
                byproduct_reactant_idxs = detect_byproducts(product_combined, reverse_mapping, delete_atoms)

                # Combine all mapping data into single dataframe
                df_combined = pd.concat([
                    df,
                    pd.Series(first_shell, name="first_shell"),
                    pd.Series(initiator_idxs, name="initiators"),
                    pd.Series(byproduct_reactant_idxs, name="byproduct_idx")
                ], axis=1).astype(pd.Int64Dtype())

                total_products = len(reaction_metadata) + 1

                # Clear isotopes before saving to restore normal chemistry
                clear_isotopes(reactant_combined, product_combined)

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
            
    def _process_reaction_instances(self, detected_reactions: list[ReactionInstance]) -> list[ReactionMetadata]:
        """
        Converts ReactionInstance objects into ReactionMetadata by building molecules and running reactions.
        
        Args:
            detected_reactions: List of detected reaction instances
            csv_cache: Path to the CSV cache file
        Returns:
            List of ReactionMetadata objects with atom mappings
        """
        csv_cache = self.csv_cache  # Ensure csv_cache is available for processing
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
            rxn = build_reaction(rxn_smarts)
            # This function also runs the reaction and builds metadata for each product set, including atom mappings and byproduct detection
            mol_reactant_1, mol_reactant_2 = build_reactants(reactant_smiles_1, reactant_smiles_2)
            # Build reaction tuple based on whether reactants are the same or different. If same, only one ordering is needed. If different, both orderings are processed to account for reaction directionality.
            reaction_tuple = build_reaction_tuple(same_reactants, mol_reactant_1, mol_reactant_2)

            # Process products and build metadata
            reaction_metadata = self._process_reaction_products(
                rxn,
                csv_cache,
                reaction_tuple,
                delete_atoms,
                reaction_metadata
            )
        
        return reaction_metadata
    
    def _process_reaction_products(self, 
                                   rxn: ChemicalReaction, 
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
            assign_atom_map_numbers_and_set_isotopes(r1, r2)
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
                if len(product_set) == 1:
                    product_combined = product_set[0]
                else:
                    product_combined = product_set[0]
                    for p in product_set[1:]:
                        product_combined = Chem.CombineMols(product_combined, p)

                # Restore atom map numbers from isotopes (which survive the reaction)
                reassign_atom_map_numbers_by_isotope(product_combined)
                
                # Build bidirectional atom index mappings
                mapping_dict, df = build_atom_index_mapping(reactant_combined, product_combined)
                reverse_mapping = {v: k for k, v in mapping_dict.items()}

                # Restore map numbers for visualization
                reveal_template_map_numbers(product_combined)

                # Validate mapping consistency
                validate_mapping(df, reactant_combined, product_combined)

                # Identify atoms involved in reaction center and initiators
                first_shell, initiator_idxs = assign_first_shell_and_initiators(
                    reactant_combined,
                    product_combined,
                    reverse_mapping
                )

                # Detect byproducts (smallest fragments)
                byproduct_reactant_idxs = detect_byproducts(product_combined, reverse_mapping, delete_atoms)

                # Combine all mapping data into single dataframe
                df_combined = pd.concat([
                    df,
                    pd.Series(first_shell, name="first_shell"),
                    pd.Series(initiator_idxs, name="initiators"),
                    pd.Series(byproduct_reactant_idxs, name="byproduct_idx")
                ], axis=1).astype(pd.Int64Dtype())

                total_products = len(reaction_metadata) + 1

                # Clear isotopes before saving to restore normal chemistry
                clear_isotopes(reactant_combined, product_combined)

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