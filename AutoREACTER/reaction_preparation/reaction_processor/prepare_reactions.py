
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdChemReactions, rdmolops
from rdkit.Chem.rdchem import Mol
from PIL.Image import Image
import itertools
import pandas as pd

from AutoREACTER.detectors.reaction_detector import ReactionInstance, ReactionTemplate
from AutoREACTER.detectors.functional_groups_detector import FunctionalGroupInfo, MonomerRole
from AutoREACTER.reaction_preparation.reaction_processor.atom_mapping import smart_mapping
from AutoREACTER.reaction_preparation.reaction_processor.utils import compare_products
from AutoREACTER.reaction_preparation.reaction_processor.walker import reaction_atom_walker

class MappingError(Exception):
    """Custom exception for errors in atom mapping."""
    pass

class SMARTSParsingError(Exception):
    """Custom exception for errors in parsing reaction SMARTS."""
    pass

@dataclass(slots=True)
class ReactionMetadata:
    reaction_id: int
    reactant_combined_mol: Chem.Mol
    product_combined_mol: Chem.Mol
    reactant_to_product_mapping: Dict[int, int]
    template_reactant_to_product_mapping: Dict[int, int]
    edge_atoms: List[int]

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
    def __init__(self, cache: Path):
        self.cache = self._prepare_paths(cache)

    def _add_column_safe(self, df, list_data, column_name):
        """
        Safely adds a list as a new column to a DataFrame, handling potential 
        length mismatches by using Series alignment.

        Args:
            df (pd.DataFrame): Target DataFrame.
            list_data (list): Data to be added to the column.
            column_name (str): Name of the new column.

        Returns:
            pd.DataFrame: Modified DataFrame with the new column.
        """
        # Creating a Series from the list ensures it starts from the top (index 0)
        # and fills missing rows with NaN if the list is shorter than the DataFrame
        df[column_name] = pd.Series(list_data).astype("Int64")
        return df
    
    def _add_dict_as_new_columns(self, df_existing, data_dict, titles=["template_reactant_idx", "template_product_idx"]):
        # Convert dict keys and values to Series to ensure alignment with the DataFrame
        # Use .astype("Int64") to allow for potential Null/NaN values while keeping integers
        df_existing[titles[0]] = pd.Series(list(data_dict.keys())).astype("Int64")
        df_existing[titles[1]] = pd.Series(list(data_dict.values())).astype("Int64")
    
        return df_existing
    
    def _reaction_tuples(self, same_reactants, mol_reactant_1, mol_reactant_2):
        
        if same_reactants:
            return [[mol_reactant_1, mol_reactant_1]]

        return [[mol_reactant_1, mol_reactant_2], [mol_reactant_2, mol_reactant_1]]
    
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
                        return t
        raise ValueError("No matching atom found in the provided set of tuples")
    
    def _build_reactants(self, reactant_smiles_1, reactant_smiles_2):
        """
        Create RDKit molecule objects from SMILES strings with explicit hydrogens.
        
        Args:
            reactant_smiles_1 (str): SMILES string for first reactant
            reactant_smiles_2 (str): SMILES string for second reactant
            
        Returns:
            tuple: (mol_reactant_1, mol_reactant_2) - Both molecules with explicit hydrogens
            
        Example:
            >>> mol1, mol2 = self._build_reactants("CCO", "CC=O")
        """
        mol_reactant_1 = Chem.MolFromSmiles(reactant_smiles_1)
        mol = Chem.MolFromSmiles(reactant_smiles_1)
        if mol is None:
            raise SMARTSParsingError(f"Invalid SMILES: {reactant_smiles_1}")
        mol_reactant_1 = Chem.AddHs(mol_reactant_1)  # Add explicit hydrogens
        mol_reactant_2 = Chem.MolFromSmiles(reactant_smiles_2)
        if mol_reactant_2 is None:
            raise SMARTSParsingError(f"Invalid SMILES: {reactant_smiles_2}")
        mol_reactant_2 = Chem.AddHs(mol_reactant_2)  # Add explicit hydrogens
        return mol_reactant_1, mol_reactant_2

    def _prepare_paths(self, cache):
        csv_cache = Path(cache) / "atom_mapping_btw_reactants_products_in_csv"
        csv_cache.mkdir(parents=True, exist_ok=True)
        return csv_cache
    
    def _build_reaction(self, rxn_smarts):
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
        if rxn is None:
            raise SMARTSParsingError(f"Invalid reaction SMARTS: {rxn_smarts}")
        return rxn

    def _is_consecutive(self, num_list):
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
        
        # Check for consecutive sequence: unique numbers and max-min+1 equals length
        return (len(set(num_list)) == len(num_list) and 
                max(num_list) - min(num_list) + 1 == len(num_list))

    
    def _process_reactions(self, rxn, csv_cache, reaction_tuple, key=None, 
                      molecule_and_csv_path_dict=None, delete_atoms=False):
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
            mols, molecule_and_csv_path_dict:
                - mols: List of reactant and product molecules for visualization
                - molecule_and_csv_path_dict: Dictionary containing all processed data
                
        Raises:
            ValueError: If mapping validation fails at any step
            
        Note:
            This function performs extensive validation including:
            - Atom count consistency between reactants and products
            - Consecutive indexing requirements
            - Initiator count validation (must be exactly 2)
        """
        if molecule_and_csv_path_dict is None:
            molecule_and_csv_path_dict = {}
        
        total_products = 0
        mols = []

        # Process each reactant pair in the reaction tuple
        for j, pair in enumerate(reaction_tuple):
            r1, r2 = Chem.Mol(pair[0]), Chem.Mol(pair[1])
            products = rxn.RunReactants((r1, r2))
            if not products:
                continue  # Skip if no products are generated
            
            # Get substructure matches for each reactant against reaction templates
            matches_1 = r1.GetSubstructMatches(rxn.GetReactants()[0])
            matches_2 = r2.GetSubstructMatches(rxn.GetReactants()[1])

            # Process each product set generated by the reaction
            for product_set in products:

                # Initialize data structures
                df = pd.DataFrame(columns=["reactant_idx", "product_idx"])
                first_shell = []      # Atoms in first coordination shell
                initiator_idxs = []   # Initiator atom indices

                # Apply atom mapping from reaction properties
                num_total_atoms = 0
                for i, product in enumerate(product_set):
                    if i > 0:
                        num_total_atoms += product_set[i - 1].GetNumAtoms()
                    for atom in product.GetAtoms():
                        if atom.HasProp("react_idx"):
                            r_idx = atom.GetIntProp("react_idx")
                            a_idx = atom.GetIntProp("react_atom_idx")
                            r = (r1, r2)[r_idx]
                            r_atom = r.GetAtomWithIdx(a_idx)
                            map_num = atom.GetIdx() + 1 + num_total_atoms + 100
                            r_atom.SetAtomMapNum(map_num)
                            atom.SetAtomMapNum(map_num)

                # Combine reactants and products for easier processing
                reactant_combined = Chem.CombineMols(r1, r2)
                product_combined = Chem.CombineMols(*product_set)

                # Create mapping between reactant and product atoms based on map numbers
                for r_atom in reactant_combined.GetAtoms():
                    for p_atom in product_combined.GetAtoms():
                        if r_atom.GetAtomMapNum() == p_atom.GetAtomMapNum():
                            new_row = pd.DataFrame([{"reactant_idx": r_atom.GetIdx(), 
                                                    "product_idx": p_atom.GetIdx()}])
                            df = pd.concat([df, new_row], ignore_index=True)
                            break

                # Validate mapping completeness and consistency
                num_reactant_atoms = reactant_combined.GetNumAtoms()
                num_product_atoms = product_combined.GetNumAtoms()
                reactant_mapped = df["reactant_idx"].notna().sum()
                product_mapped = df["product_idx"].notna().sum()

                # Validation 1: Mapped atom counts must match between columns
                if reactant_mapped != product_mapped:
                    raise MappingError(
                        "Mismatch in mapped atom counts between columns: "
                        f"reactant_idx mapped={reactant_mapped}, product_idx mapped={product_mapped}"
                    )

                # Validation 2: All reactant atoms must be mapped
                if reactant_mapped != num_reactant_atoms:
                    df.to_csv(csv_cache / f"debug_mapping_{key}_{total_products}.csv", index=False)
                    raise MappingError(
                        "Mapping does not cover all reactant atoms: "
                        f"reactant atoms={num_reactant_atoms}, mapped={reactant_mapped}"
                    )

                # Validation 3: All product atoms must be mapped
                if product_mapped != num_product_atoms:
                    raise MappingError(
                        "Mapping does not cover all product atoms: "
                        f"product atoms={num_product_atoms}, mapped={product_mapped}"
                    )

                # Validation 4: Mapping indices must be consecutive
                if (self._is_consecutive(df["reactant_idx"].tolist()) is False or 
                    self._is_consecutive(df["product_idx"].tolist()) is False):
                    raise MappingError(
                        "Mapping indices are not consecutive: "
                        f"reactant indices={df['reactant_idx'].tolist()}, "
                        f"product indices={df['product_idx'].tolist()}"
                    )

                # Apply smart mapping to both reactants
                try:
                    sub_set1 = self._is_number_in_set(matches_1, r1)
                except ValueError:
                    sub_set1 = ()
                smart_mapping(
                    reactant=r1,
                    smarts_template=rxn.GetReactants()[0],
                    match_tuple=sub_set1 if sub_set1 else (),
                )
                
                try:
                    sub_set2 = self._is_number_in_set(matches_2, r2)
                except ValueError:
                    sub_set2 = ()
                smart_mapping(
                    reactant=r2,
                    smarts_template=rxn.GetReactants()[1],
                    match_tuple=sub_set2 if sub_set2 else (),
                )

                # Recombine reactants after mapping
                reactant_combined = Chem.CombineMols(r1, r2)
                
                # Identify first shell atoms and initiators
                for atom in reactant_combined.GetAtoms():
                    map_num = atom.GetAtomMapNum()
                    idx = atom.GetIdx()
                    if 0 < map_num <= 99:
                        first_shell.append(idx)

                    if 0 < map_num <= 2:
                        initiator_idxs.append(idx)
                
                # Identify byproduct atoms if deletion is requested
                byproduct_product_idxs = []
                byproduct_reactant_idxs = []

                if delete_atoms:
                    frags = rdmolops.GetMolFrags(product_combined, asMols=True)
                    smallest_mol = min(frags, key=lambda m: m.GetNumAtoms())

                    # 1) get byproduct atom indices in PRODUCT space
                    for atom in smallest_mol.GetAtoms():
                        byproduct_map_number = atom.GetAtomMapNum()
                        for p_atom in product_combined.GetAtoms():
                            if p_atom.GetAtomMapNum() == byproduct_map_number:
                                byproduct_product_idxs.append(p_atom.GetIdx())
                                break

                    # 2) convert PRODUCT indices -> REACTANT indices using df mapping
                    byproduct_reactant_idxs = (
                        df.loc[df["product_idx"].isin(byproduct_product_idxs), "reactant_idx"]
                        .astype(int)
                        .tolist()
                    )

                # if you want the column to store REACTANT indices now:
                byproduct_indexs = byproduct_reactant_idxs

                # Create comprehensive DataFrame with all analysis columns
                first_shell_column = pd.Series(first_shell, name="first_shell")
                initiator_idxs_column = pd.Series(initiator_idxs, name="initiators")
                
                # Validation: Must have exactly 2 initiators
                if len(initiator_idxs) != 2:
                    raise ValueError(f"Expected 2 initiators, got {len(initiator_idxs)}: {initiator_idxs}")
                
                by_product_indexs_column = pd.Series(byproduct_indexs, name="byproduct_indices")
                
                # Combine all data into final DataFrame
                df_combined = pd.concat([df, first_shell_column, initiator_idxs_column,
                                        by_product_indexs_column], axis=1).astype(pd.Int64Dtype())
                
                # Check for duplicate products
                if not compare_products(molecule_and_csv_path_dict, product_combined):
                    continue

                # Update counters and store results
                total_products += 1
                mols.append(reactant_combined)
                mols.append(product_combined)

                # Organize results in the output dictionary
                for i in itertools.count(1):
                    if i not in molecule_and_csv_path_dict: # Check if i IS a key
                        sub_dict = molecule_and_csv_path_dict[i] = {} # Initialize if not found
                        dict_key = i
                        break

                # Save to CSV file
                df_combined.to_csv(csv_cache / f"reaction_{dict_key}.csv", index=False)
                print(f"Saved reaction {dict_key} to CSV")

                sub_dict["reactant"] = reactant_combined
                sub_dict["product"] = product_combined
                sub_dict["csv_path"] = csv_cache / f"reaction_{dict_key}.csv"
                sub_dict["delete_atoms"] = delete_atoms
                sub_dict["reaction_dataframe"] = df_combined


        return mols, molecule_and_csv_path_dict
    
    def prepare_reactions(self, reactions: List[ReactionInstance],) -> List[ReactionMetadata]:

        csv_cache = self._prepare_paths(self.cache)
        metadata_list = []

        for reaction in reactions:

            rxn = self._build_reaction(reaction.reaction_smarts)

            mol_r1, mol_r2 = self._build_reactants(
                reaction.monomer_1.smiles,
                reaction.monomer_2.smiles
            )

            reaction_tuple = self._reaction_tuples(
                reaction.same_reactants,
                mol_r1,
                mol_r2
            )

            mols, result_dict = self._process_reactions(
                rxn,
                csv_cache,
                reaction_tuple,
                key=reaction.reaction_id
            )

            for rid, data in result_dict.items():

                df = data["reaction_dataframe"]

                reactant_mol = data["reactant"]

                # FULL atom mapping
                fully_mapped_dict = dict(zip(
                    df["reactant_idx"].dropna(),
                    df["product_idx"].dropna()
                ))

                # first shell atoms
                first_shell = df["first_shell"].dropna().astype(int).tolist()

                # WALKER
                template_mapped_dict, edge_atoms = reaction_atom_walker(
                    reactant_mol,
                    first_shell,
                    fully_mapped_dict
                )

                # update dataframe
                df = self._add_dict_as_new_columns(
                    df,
                    template_mapped_dict,
                    titles=["template_reactant_idx", "template_product_idx"]
                )

                df = self._add_column_safe(
                    df,
                    edge_atoms,
                    "edge_atoms"
                )

                # store updated dataframe
                data["reaction_dataframe"] = df.copy()

                # build metadata
                metadata = ReactionMetadata(
                    reaction_id=rid,
                    reactant_combined_mol=data["reactant"],
                    product_combined_mol=data["product"],
                    reactant_to_product_mapping=fully_mapped_dict,
                    template_reactant_to_product_mapping=template_mapped_dict,
                    edge_atoms=edge_atoms,
                    reaction_smarts=reaction.reaction_smarts,
                    csv_path=data["csv_path"],
                    reaction_dataframe=df
                )

                metadata_list.append(metadata)

        return metadata_list
    


    def reaction_templates_highligted_image_grid(
        self,
        metadata_list: List[ReactionMetadata],
        highlight_type: str = "template"
    ) -> Image:
        """
        Create a grid image highlighting atoms in reactant molecules.

        highlight_type options:
            "template"
            "edge"
            "initiators"
        """

        mols = []
        highlight_lists = []

        for metadata in metadata_list:

            reactant = metadata.reactant_combined_mol
            product = metadata.product_combined_mol
            df = metadata.reaction_dataframe

            if highlight_type == "template":
                atoms = list(metadata.template_reactant_to_product_mapping.keys())

            elif highlight_type == "edge":
                atoms = metadata.edge_atoms

            elif highlight_type == "initiators":
                atoms = df["initiators"].dropna().astype(int).tolist()

            else:
                atoms = []

            mols.extend([reactant, product])
            highlight_lists.append(atoms)
            highlight_lists.append([])  # don't highlight product

        img = Draw.MolsToGridImage(
            mols,
            molsPerRow=2,
            highlightAtomLists=highlight_lists,
            subImgSize=(400, 400),
            useSVG=False
        )

        return img

            
