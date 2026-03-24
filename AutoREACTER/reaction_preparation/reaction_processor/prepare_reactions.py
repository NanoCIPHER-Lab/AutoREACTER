"""
Module for preparing chemical reactions for analysis, including atom mapping between reactants and products,
reaction metadata extraction, and visualization utilities using RDKit.

This module processes reaction SMARTS, applies atom mappings, identifies key reaction features (e.g., first shell,
initiators, byproducts), and generates metadata and visualizations for downstream tasks like reaction prediction
or mechanistic analysis.
"""

# WARNING:
# When modifying this file for dataframe or any other indexing varibales use idx, do not use index or indicies or similar.

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
from AutoREACTER.reaction_preparation.reaction_processor.atom_mapping import smart_mapping
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
    Dataclass holding comprehensive metadata for a processed reaction instance.

    Attributes:
        reaction_id: Unique identifier for the reaction instance.
        reactant_combined_mol: RDKit Mol object combining all reactant molecules.
        product_combined_mol: RDKit Mol object combining all product molecules.
        reactant_to_product_mapping: Dictionary mapping reactant atom indices to product atom indices.
        template_reactant_to_product_mapping: Dictionary mapping template-matched reactant atoms to products.
        edge_atoms: List of atom indices on the 'edge' of the reaction (e.g., involved in changes).
        first_shell: List of atoms in the first coordination shell around reaction centers (optional).
        initiators: List of initiator atom indices (typically 2 per reaction).
        byproduct_indices: List of byproduct atom indices (optional).
        reaction_smarts: Original reaction SMARTS string (optional).
        reactant_smarts: Reactant SMARTS string (optional).
        product_smiles: Product SMILES string (optional).
        csv_path: Path to CSV file storing reaction mapping data (optional).
        mol_3d_path: Path to 3D molecule file (optional).
        reaction_dataframe: Pandas DataFrame with detailed atom mapping and features.
        delete_atom: Flag indicating if byproduct atoms should be deleted.
        delete_atom_idx: Specific index of atom to delete (optional).
        activity_stats: Flag for including activity statistics (optional).
    """
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
        """
        Main execution pipeline for mapping reactions, identifying templates, 
        and saving results to CSV.

        Steps:
        1. Maps atoms between reactants and products.
        2. Uses an atom walker to identify the 'reaction center' (template).
        3. Identifies 'edge atoms' that connect the template to the rest of the molecule.
        4. Updates reaction dataframes with mapping indices and saves them.

        Args:
            reaction_instances (list[ReactionInstance]): List of reaction instances to process.
            cache (str): Path to the cache directory for processing.

        Returns:
            tuple: (updated_molecule_dict, formatted_summary_dict)
        """
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
        """
        Process a dictionary of detected reactions and generate mapping data.
        
        This is a high-level function that processes multiple reactions defined in
        a dictionary structure, handling all steps from SMILES to final CSV output.
        
        Args:
            detected_reactions (list[ReactionInstance]): List of detected reaction instances
            cache (Path): Base cache directory for output files
            
        Returns:
            tuple: (molecule_and_csv_path_dict, detected_reactions)
                - molecule_and_csv_path_dict: Dictionary with all processed reaction data
                - detected_reactions: The original input dictionary (unchanged)
                
        Example:
            >>> reactions_dict = {
            ...     1: {"reaction": "[C:1]=[O:2]>>[C:1]-[O:2]", 
            ...         "monomer_1": {"smiles": "CCO"}, ...}
            ... }
            >>> results, original = process_reaction_instances(reactions_dict, "/path/to/cache")
        """
        csv_cache = self.csv_cache
        all_metadata = []

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
                                   )
            all_metadata.extend(reaction_metadata)
        
        return all_metadata
    
    def process_reaction_products(self, 
                                   rxn: Chem.rdChemReactions.ChemicalReaction, 
                                   csv_cache: Path, 
                                   reaction_tuple: list, 
                                   delete_atoms: bool = True
                                   ) -> list[ReactionMetadata]:
        """
        Process chemical reactions, map atoms between reactants and products, and save data.
        
        This is the main processing function that:
        1. Runs reactions on reactant pairs
        2. Maps atoms between reactants and products
        3. Validates mapping consistency
        4. Saves mapping data to CSV files
        5. Identifies reaction features (first shell, initiators, byproducts)
        
        Args:
            rxn (Chem.rdChemReactions.ChemicalReaction): RDKit reaction object
            csv_cache (Path): Directory path for CSV output files
            reaction_tuple (list): List of reactant pairs to process
            key (str or int, optional): Identifier for the reaction set
            molecule_and_csv_path_dict (dict, optional): Dictionary to store results
            delete_atoms (bool, optional): Whether to identify byproduct atoms for deletion
            
        Returns:
            tuple: (mols, output, molecule_and_csv_path_dict)
                - mols: List of reactant and product molecules for visualization
                - output: String output (currently empty, reserved for future use)
                - molecule_and_csv_path_dict: Dictionary containing all processed data
                
        Raises:
            ValueError: If mapping validation fails at any step
            
        Note:
            This function performs extensive validation including:
            - Atom count consistency between reactants and products
            - Consecutive indexing requirements
            - Initiator count validation (must be exactly 2)
        """
        reaction_metadata = []
        
        # Process each reactant pair in the reaction tuple
        for pair in reaction_tuple:
            r1, r2 = Chem.Mol(pair[0]), Chem.Mol(pair[1])
            products = rxn.RunReactants((r1, r2))
            
            if not products:
                continue  # Skip if no products are generated

            # Get substructure matches for each reactant against reaction templates
            matches_1 = r1.GetSubstructMatches(rxn.GetReactants()[0])
            matches_2 = r2.GetSubstructMatches(rxn.GetReactants()[1])
            
            # Process each product set generated by the reaction
            for product_set in products:
                # loop over ALL match combinations
                for match1 in matches_1 if matches_1 else [()]:
                    for match2 in matches_2 if matches_2 else [()]:

                        # IMPORTANT: copy reactants for each combination
                        r1_copy = Chem.Mol(r1)
                        r2_copy = Chem.Mol(r2)

                        # Clean up any existing DataFrame
                        if "df" in locals():
                            del df

                        df = pd.DataFrame(columns=["reactant_idx", "product_idx"])
                        first_shell = []
                        initiator_idxs = []
                        mapping_dict = {}


                        # Apply atom mapping from reaction properties
                        num_total_atoms = 0
                        for i, product in enumerate(product_set):
                            if i > 0:
                                num_total_atoms += product_set[i - 1].GetNumAtoms()
                            for atom in product.GetAtoms():
                                if atom.HasProp("react_idx"):
                                    r_idx = atom.GetIntProp("react_idx")
                                    a_idx = atom.GetIntProp("react_atom_idx")
                                    r = (r1_copy, r2_copy)[r_idx]
                                    r_atom = r.GetAtomWithIdx(a_idx)
                                    map_num = atom.GetIdx() + num_total_atoms + 101
                                    r_atom.SetAtomMapNum(map_num)
                                    atom.SetAtomMapNum(map_num)

                        # Combine reactants and products
                        reactant_combined = Chem.CombineMols(r1_copy, r2_copy)
                        product_combined = Chem.CombineMols(*product_set)

                        # Create mapping dictionary and dataframe of reactant to product atom indices
                        for r_atom in reactant_combined.GetAtoms():
                            for p_atom in product_combined.GetAtoms():
                                if r_atom.GetAtomMapNum() == p_atom.GetAtomMapNum():
                                    r_atom_idx = r_atom.GetIdx()
                                    p_atom_idx = p_atom.GetIdx()
                                    mapping_dict[r_atom_idx] = p_atom_idx
                                    new_row = pd.DataFrame([{
                                        "reactant_idx": r_atom_idx,
                                        "product_idx": p_atom_idx
                                    }])
                                    df = pd.concat([df, new_row], ignore_index=True)
                                    break

                        if not self._validate_mapping(df, reactant_combined, product_combined):
                            raise MappingError("Atom mapping validation failed: inconsistent mapping between reactants and products.")

                        # Identify first shell + initiators
                        for atom in reactant_combined.GetAtoms():
                            if atom.GetAtomMapNum() <= 99:
                                first_shell.append(atom.GetIdx())
                            if atom.GetAtomMapNum() in [1, 2]:
                                initiator_idxs.append(atom.GetIdx())

                        if len(initiator_idxs) != 2:
                            raise ValueError(f"Expected 2 initiators, got {len(initiator_idxs)}: {initiator_idxs}") 
                        
                        """
                        This step enforces SMARTS-based atom mapping onto the reactants using selected substructure matches. 
                        It assigns atom map numbers from the SMARTS template to corresponding atoms in the reactants, ensuring 
                        consistent labeling of reaction centers.

                        For each reactant, one matching substructure is selected and used to apply mapping. This helps standardize 
                        atom indices so that downstream processes, such as the atom walker, can identify reaction templates and edge atoms.
                        """
                        # Apply SMARTS-based mapping for THIS combination
                        smart_mapping(
                            reactant=r1_copy,
                            smarts_template=rxn.GetReactants()[0],
                            match_tuple=match1 if match1 else (),
                        )

                        smart_mapping(
                            reactant=r2_copy,
                            smarts_template=rxn.GetReactants()[1],
                            match_tuple=match2 if match2 else (),
                        )

                        # Byproduct detection
                        byproduct_reactant_idxs = []

                        if delete_atoms:
                            byproduct_reactant_idxs = self._detect_byproducts(product_combined, mapping_dict, delete_atoms)

                        # Build dataframe
                        df_combined = pd.concat([
                            df,
                            pd.Series(first_shell, name="first_shell"),
                            pd.Series(initiator_idxs, name="initiators"),
                            pd.Series(byproduct_reactant_idxs, name="byproduct_idx")
                        ], axis=1).astype(pd.Int64Dtype())

                        total_products = len(reaction_metadata) + 1 if reaction_metadata else 1

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
        reverse_map = {v: k for k, v in mapping_dict.items()}

        for atom in smallest.GetAtoms():
            p_idx = atom.GetIdx()
            if p_idx in reverse_map:
                byproduct.append(reverse_map[p_idx])

        return byproduct

    def _is_consecutive(self, num_list: list[int]) -> bool:
        """
        Check if a list of numbers forms a consecutive sequence without gaps.
        
        Args:
            num_list (list): List of integers to check
            
        Returns:
            bool: True if the numbers are consecutive, False otherwise
            
        Examples:
            >>> is_consecutive([1, 2, 3, 4])  # Returns True
            >>> is_consecutive([1, 3, 4])     # Returns False
            >>> is_consecutive([])            # Returns False
        """
        if not num_list:
            return False

        return (
            len(set(num_list)) == len(num_list)
            and max(num_list) - min(num_list) + 1 == len(num_list)
    )
    
    def _is_number_in_set(self, set_of_tuples, reactant):
        """
        Find a tuple in the set that contains atom indices matching specific atom map numbers.
        
        This function searches for atoms in the reactant molecule that have atom map numbers
        101 or 102, then checks if those atom indices exist in any tuple within the provided set.
        
        Args:
            set_of_tuples (set of tuples): Set of atom index tuples to search through
            reactant (Chem.Mol): RDKit molecule object to search for mapped atoms
            
        Returns:
            tuple: The first tuple found containing a matching atom index
            
        Raises:
            ValueError: If no matching atom is found in the provided set of tuples
            
        Example:
            >>> matches = {(0, 1, 2), (3, 4, 5)}
            >>> result = self._is_number_in_set(matches, reactant_mol)
        """
        for atom in reactant.GetAtoms():
            # Check for specific atom map numbers (101 or 102)
            if atom.GetAtomMapNum() == 101 or atom.GetAtomMapNum() == 102:
                x = atom.GetIdx()
                for t in set_of_tuples:
                    if x in t:
                        print(f"Found matching atom index {x} in tuple {t}")
                        return t
        raise ValueError("No matching atom found in the provided set of tuples.")
    
    def _build_reaction_tuple(self, same_reactants, mol_reactant_1, mol_reactant_2):
        """
        Generate reactant pair orderings for reaction processing.

        - For same reactants: only one unique ordering [A, A].
        - For different reactants: return both [A, B] and [B, A] to avoid relying
        on reactant ordering and to support order-dependent reaction definitions.
        """
        if same_reactants:
            return [[mol_reactant_1, mol_reactant_1]]

        return [[mol_reactant_1, mol_reactant_2], [mol_reactant_2, mol_reactant_1]]
    
    
    def _build_reaction(self, rxn_smarts):
        """
        Create an RDKit reaction object from SMARTS string.
        
        Args:
            rxn_smarts (str): Reaction SMARTS string
            
        Returns:
            Chem.rdChemReactions.ChemicalReaction: RDKit reaction object
            
        Example:
            >>> rxn = build_reaction("[C:1]=[O:2]>>[C:1]-[O:2]")
        """
        return AllChem.ReactionFromSmarts(rxn_smarts)


    def _build_reactants(self, reactant_smiles_1, reactant_smiles_2):
        """
        Create RDKit molecule objects from SMILES strings with explicit hydrogens.
        
        Args:
            reactant_smiles_1 (str): SMILES string for first reactant
            reactant_smiles_2 (str): SMILES string for second reactant
            
        Returns:
            tuple: (mol_reactant_1, mol_reactant_2) - Both molecules with explicit hydrogens
            
        Example:
            >>> mol1, mol2 = build_reactants("CCO", "CC=O")
        """
        mol_reactant_1 = Chem.MolFromSmiles(reactant_smiles_1)
        mol_reactant_1 = Chem.AddHs(mol_reactant_1)  # Add explicit hydrogens
        mol_reactant_2 = Chem.MolFromSmiles(reactant_smiles_2)
        mol_reactant_2 = Chem.AddHs(mol_reactant_2)  # Add explicit hydrogens
        return mol_reactant_1, mol_reactant_2

    def _is_initiator_in_set(self, set_of_tuples, reactant):
        """
        Find a tuple in the set that contains atom indices matching specific atom map numbers.
        
        This function searches for atoms in the reactant molecule that have atom map numbers
        101 or 102, then checks if those atom indices exist in any tuple within the provided set.
        
        Args:
            set_of_tuples (set of tuples): Set of atom index tuples to search through
            reactant (Chem.Mol): RDKit molecule object to search for mapped atoms
            
        Returns:
            tuple: The first tuple found containing a matching atom index
            
        Raises:
            ValueError: If no matching atom is found in the provided set of tuples
            
        Example:
            >>> matches = {(0, 1, 2), (3, 4, 5)}
            >>> result = is_initiator_in_set(matches, reactant_mol)
        """
        for atom in reactant.GetAtoms():
            # Check for specific atom map numbers (101 or 102)
            if atom.GetAtomMapNum() == 101 or atom.GetAtomMapNum() == 102:
                x = atom.GetIdx()
                for t in set_of_tuples:
                    if x in t:
                        print(f"Found matching atom index {x} in tuple {t}")
                        return t
        raise ValueError("No matching atom found in the provided set of tuples.")
    

    def reaction_templates_highlighted_image_grid(
        self,
        metadata_list: List[ReactionMetadata],
        highlight_type: str = "template",
    ) -> Image:
        """
        Generate a grid image of reactant-product pairs with highlighted atoms.

        Clears atom maps before drawing. Supports multiple highlight types.

        Args:
            metadata_list: List of ReactionMetadata.
            highlight_type: Type of atoms to highlight:
                - "template": Template-matched atoms (blue).
                - "edge": Edge atoms (orange).
                - "initiators": Initiator atoms (green).
                - "delete": Byproduct atoms if delete_atom=True (red).

        Returns:
            PIL Image of the grid.
        """
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
    

