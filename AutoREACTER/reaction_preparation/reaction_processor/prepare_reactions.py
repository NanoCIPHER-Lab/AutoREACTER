"""
Module for preparing chemical reactions for analysis, including atom mapping between reactants and products,
reaction metadata extraction, and visualization utilities using RDKit.

This module processes reaction SMARTS, applies atom mappings, identifies key reaction features (e.g., first shell,
initiators, byproducts), and generates metadata and visualizations for downstream tasks like reaction prediction
or mechanistic analysis.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdmolops
from PIL.Image import Image
import pandas as pd

from AutoREACTER.detectors.reaction_detector import ReactionInstance
from AutoREACTER.reaction_preparation.reaction_processor.atom_mapping import smart_mapping
from AutoREACTER.reaction_preparation.reaction_processor.utils import compare_products
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
    template_reactant_to_product_mapping: Dict[int, int]
    edge_atoms: List[int]

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
        self.csv_cache = self._prepare_paths(self.cache)

    def _add_column_safe(self, df: pd.DataFrame, list_data: List, column_name: str) -> pd.DataFrame:
        """
        Safely add a new integer column to a DataFrame, handling NaN values with nullable Int64 dtype.

        Args:
            df: Input DataFrame.
            list_data: List of data to add as a new column.
            column_name: Name of the new column.

        Returns:
            DataFrame with the new column added.
        """
        df[column_name] = pd.Series(list_data).astype("Int64")
        return df

    def _add_dict_as_new_columns(
        self,
        df_existing: pd.DataFrame,
        data_dict: Dict[int, int],
        titles: List[str] = ["template_reactant_idx", "template_product_idx"]
    ) -> pd.DataFrame:
        """
        Add two new columns from a dictionary keys and values, using nullable Int64 dtype.

        Args:
            df_existing: Input DataFrame.
            data_dict: Dictionary where keys go to titles[0] and values to titles[1].
            titles: List of two column names for keys and values.

        Returns:
            DataFrame with the two new columns added.
        """
        df_existing[titles[0]] = pd.Series(list(data_dict.keys())).astype("Int64")
        df_existing[titles[1]] = pd.Series(list(data_dict.values())).astype("Int64")
        return df_existing

    def _reaction_tuples(
        self,
        same_reactants: bool,
        mol_reactant_1: Chem.Mol,
        mol_reactant_2: Chem.Mol
    ) -> List[List[Chem.Mol]]:
        """
        Generate reactant pairs for reaction processing, handling homo-polymerization (same reactants).

        Args:
            same_reactants: If True, use the same molecule twice (homo-polymerization).
            mol_reactant_1: First reactant molecule.
            mol_reactant_2: Second reactant molecule.

        Returns:
            List of reactant pairs [[r1, r2], ...].
        """
        if same_reactants:
            return [[mol_reactant_1, mol_reactant_1]]
        return [[mol_reactant_1, mol_reactant_2], [mol_reactant_2, mol_reactant_1]]

    def _is_number_in_set(self, set_of_tuples: List[tuple], reactant: Chem.Mol) -> tuple:
        """
        Find the substructure match tuple containing an atom with map number >= 100.

        Args:
            set_of_tuples: List of match tuples from GetSubstructMatches.
            reactant: RDKit Mol of the reactant.

        Returns:
            The matching tuple.

        Raises:
            ValueError: If no matching atom/tuple found.
        """
        for atom in reactant.GetAtoms():
            if atom.GetAtomMapNum() >= 100:
                x = atom.GetIdx()
                for t in set_of_tuples:
                    if x in t:
                        return t
        raise ValueError("No matching atom found in the provided set of tuples")

    def _build_reactants(self, reactant_smiles_1: str, reactant_smiles_2: str) -> tuple[Chem.Mol, Chem.Mol]:
        """
        Parse SMILES strings into RDKit Mol objects, add hydrogens, and validate.

        Args:
            reactant_smiles_1: SMILES for first reactant.
            reactant_smiles_2: SMILES for second reactant.

        Returns:
            Tuple of (mol_reactant_1, mol_reactant_2).

        Raises:
            SMARTSParsingError: If SMILES parsing fails.
        """
        mol_reactant_1 = Chem.MolFromSmiles(reactant_smiles_1)
        if mol_reactant_1 is None:
            raise SMARTSParsingError(f"Invalid SMILES: {reactant_smiles_1}")
        mol_reactant_1 = Chem.AddHs(mol_reactant_1)

        mol_reactant_2 = Chem.MolFromSmiles(reactant_smiles_2)
        if mol_reactant_2 is None:
            raise SMARTSParsingError(f"Invalid SMILES: {reactant_smiles_2}")
        mol_reactant_2 = Chem.AddHs(mol_reactant_2)

        return mol_reactant_1, mol_reactant_2

    def _prepare_paths(self, cache: Path) -> Path:
        """
        Create the CSV cache directory for storing atom mapping data.

        Args:
            cache: Base cache Path.

        Returns:
            Path to the CSV cache subdirectory.
        """
        csv_cache = Path(cache) / "atom_mapping_btw_reactants_products_in_csv"
        csv_cache.mkdir(parents=True, exist_ok=True)
        return csv_cache

    def _build_reaction(self, rxn_smarts: str) -> AllChem.ChemicalReaction:
        """
        Parse reaction SMARTS into an RDKit ChemicalReaction object.

        Args:
            rxn_smarts: Reaction SMARTS string.

        Returns:
            RDKit ChemicalReaction object.

        Raises:
            SMARTSParsingError: If parsing fails.
        """
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
        if rxn is None:
            raise SMARTSParsingError(f"Invalid reaction SMARTS: {rxn_smarts}")
        return rxn

    def _is_consecutive(self, num_list: List[int]) -> bool:
        """
        Check if a list of numbers is consecutive without duplicates or gaps.

        Args:
            num_list: List of integers.

        Returns:
            True if consecutive and unique.
        """
        if not num_list:
            return False
        return (len(set(num_list)) == len(num_list) and
                max(num_list) - min(num_list) + 1 == len(num_list))

    def _process_reactions(
        self,
        rxn: AllChem.ChemicalReaction,
        csv_cache: Path,
        reaction_tuple: List[List[Chem.Mol]],
        key: Optional[str] = None,
        molecule_and_csv_path_dict: Optional[Dict] = None,
        delete_atom: bool = False
    ) -> Dict:
        """
        Process chemical reactions, map atoms between reactants and products, and save data.

        This is the core processing function that:
        1. Runs reactions on reactant pairs using RDKit.RunReactants.
        2. Applies atom mapping from reaction properties and validates consistency.
        3. Identifies reaction centers, first shell, initiators, and byproducts.
        4. Saves detailed mapping DataFrames to CSV files.
        5. Stores results in a dictionary keyed by local reaction IDs.

        Args:
            rxn: RDKit ChemicalReaction object.
            csv_cache: Directory to save CSV files.
            reaction_tuple: List of reactant pairs to process.
            key: Identifier for the reaction set (e.g., reaction index).
            molecule_and_csv_path_dict: Dictionary to store results (updated in-place).
            delete_atom: If True, identify and mark byproduct atoms for deletion.

        Returns:
            Dictionary of processed reaction data, including mols, DataFrames, and paths.
        """
        if molecule_and_csv_path_dict is None:
            molecule_and_csv_path_dict = {}

        total_products = 0

        # Process each reactant pair in the reaction tuple
        for j, pair in enumerate(reaction_tuple):
            r1, r2 = Chem.Mol(pair[0]), Chem.Mol(pair[1])
            products = rxn.RunReactants((r1, r2))
            if not products:
                continue

            # Get substructure matches for each reactant against reaction templates
            matches_1 = r1.GetSubstructMatches(rxn.GetReactants()[0])
            matches_2 = r2.GetSubstructMatches(rxn.GetReactants()[1])

            # Process each product set generated by the reaction
            for product_set in products:
                if "df" in locals():
                    del df

                # Initialize data structures for mapping
                df = pd.DataFrame(columns=["reactant_idx", "product_idx"])
                first_shell = []
                initiator_idxs = []
                mapping_dict = {}

                # Apply atom mapping from reaction properties (react_idx and react_atom_idx props)
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

                # Combine reactants and products for unified processing
                reactant_combined = Chem.CombineMols(r1, r2)
                product_combined = Chem.CombineMols(*product_set)
                if not compare_products(molecule_and_csv_path_dict, product_combined):
                    continue

                # Build mapping DataFrame based on shared atom map numbers
                for r_atom in reactant_combined.GetAtoms():
                    for p_atom in product_combined.GetAtoms():
                        if r_atom.GetAtomMapNum() == p_atom.GetAtomMapNum():
                            new_row = pd.DataFrame([{
                                "reactant_idx": r_atom.GetIdx(),
                                "product_idx": p_atom.GetIdx()
                            }])
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

                # Validation 2: All reactant atoms must be mapped (debug CSV if fails)
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

                # Validation 4: Mapping indices must be consecutive (no gaps/duplicates)
                if (not self._is_consecutive(df["reactant_idx"].tolist()) or
                    not self._is_consecutive(df["product_idx"].tolist())):
                    raise MappingError(
                        "Mapping indices are not consecutive: "
                        f"reactant indices={df['reactant_idx'].tolist()}, "
                        f"product indices={df['product_idx'].tolist()}"
                    )

                # Apply smart atom mapping to reactants using reaction templates
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

                # Recombine reactants after smart mapping
                reactant_combined = Chem.CombineMols(r1, r2)

                # Identify first shell atoms (map num <=99) and initiators (map num 1 or 2)
                for atom in reactant_combined.GetAtoms():
                    if atom.GetAtomMapNum() <= 99:
                        first_shell.append(atom.GetIdx())
                    if atom.GetAtomMapNum() in [1, 2]:
                        initiator_idxs.append(atom.GetIdx())

                # Identify byproduct atoms if deletion is enabled
                byproduct_product_idxs = []
                byproduct_reactant_idxs = []

                if delete_atom:
                    frags = rdmolops.GetMolFrags(product_combined, asMols=True)
                    smallest_mol = min(frags, key=lambda m: m.GetNumAtoms())

                    # Map byproduct atoms from smallest fragment back to product_combined indices
                    for atom in smallest_mol.GetAtoms():
                        byproduct_map_number = atom.GetAtomMapNum()
                        for p_atom in product_combined.GetAtoms():
                            if p_atom.GetAtomMapNum() == byproduct_map_number:
                                byproduct_product_idxs.append(p_atom.GetIdx())
                                break

                    # Convert product byproduct indices to reactant indices using df
                    byproduct_reactant_idxs = (
                        df.loc[df["product_idx"].isin(byproduct_product_idxs), "reactant_idx"]
                        .astype(int)
                        .tolist()
                    )

                byproduct_indexs = byproduct_reactant_idxs  # Store reactant indices

                # Append any additional manual mappings
                for r_idx, p_idx in mapping_dict.items():
                    new_row = pd.DataFrame([{"reactant_idx": r_idx, "product_idx": p_idx}])
                    df = pd.concat([df, new_row], ignore_index=True)

                # Create columns for first_shell, initiators, and byproducts
                first_shell_column = pd.Series(first_shell, name="first_shell")
                initiator_idxs_column = pd.Series(initiator_idxs, name="initiators")

                # Validation: Exactly 2 initiators expected
                if len(initiator_idxs) != 2:
                    raise ValueError(f"Expected 2 initiators, got {len(initiator_idxs)}: {initiator_idxs}")

                by_product_indexs_column = pd.Series(byproduct_indexs, name="byproduct_indices")

                # Combine into final DataFrame with nullable integer dtype
                df_combined = pd.concat(
                    [df, first_shell_column, initiator_idxs_column, by_product_indexs_column],
                    axis=1
                ).astype(pd.Int64Dtype())

                # Update counters and store in output dict
                total_products += 1
                local_id = len(molecule_and_csv_path_dict) + 1
                dict_key = f"{key}_{local_id}"

                sub_dict = molecule_and_csv_path_dict.setdefault(local_id, {})

                # Save DataFrame to CSV
                csv_file = csv_cache / f"reaction_{dict_key}.csv"
                df_combined.to_csv(csv_file, index=False)
                print(f"Saved reaction {dict_key} to CSV")

                sub_dict["reactant"] = reactant_combined
                sub_dict["product"] = product_combined
                sub_dict["csv_path"] = csv_file
                sub_dict["delete_atom"] = delete_atom
                sub_dict["reaction_dataframe"] = df_combined

        return molecule_and_csv_path_dict

    def _extract_reactant_info(self, reaction: ReactionInstance) -> tuple[str, str, bool]:
        """
        Extract SMILES for reactants from a ReactionInstance, handling homo- vs. co-polymerization.

        Args:
            reaction: ReactionInstance with monomer_1 and optional monomer_2.

        Returns:
            Tuple of (smiles_1, smiles_2, same_reactants_flag).

        Raises:
            ValueError: If monomer_1 is missing.
        """
        if reaction.monomer_1 is None:
            raise ValueError(f"Reaction {reaction.reaction_name} has no monomer_1")

        smiles_1 = reaction.monomer_1.smiles

        if reaction.monomer_2 is None:
            return smiles_1, smiles_1, True

        smiles_2 = reaction.monomer_2.smiles
        return smiles_1, smiles_2, False

    def prepare_reactions(self, reactions: List[ReactionInstance]) -> List[ReactionMetadata]:
        """
        Main public method to process a list of ReactionInstances into ReactionMetadata.

        For each reaction:
        1. Builds reaction object and reactant mols.
        2. Processes pairs via _process_reactions.
        3. Extracts mappings, edge atoms via reaction_atom_walker.
        4. Builds and augments DataFrames.
        5. Constructs ReactionMetadata objects.

        Args:
            reactions: List of ReactionInstance objects.

        Returns:
            List of ReactionMetadata objects.
        """
        csv_cache = self.csv_cache
        metadata_list = []

        for rxn_idx, reaction in enumerate(reactions, start=1):
            rxn = self._build_reaction(reaction.reaction_smarts)

            reactant_smiles_1, reactant_smiles_2, same_reactants = \
                self._extract_reactant_info(reaction)

            mol_r1, mol_r2 = self._build_reactants(
                reactant_smiles_1,
                reactant_smiles_2
            )

            reaction_tuple = self._reaction_tuples(
                same_reactants,
                mol_r1,
                mol_r2
            )

            result_dict = self._process_reactions(
                rxn,
                csv_cache,
                reaction_tuple,
                key=rxn_idx,
                delete_atom=reaction.delete_atom
            )

            for rid, data in result_dict.items():
                df = data["reaction_dataframe"]
                reactant_mol = data["reactant"]

                # Extract core reactant-to-product mapping (drop NaNs)
                mapping_df = df.dropna(subset=["reactant_idx", "product_idx"])
                fully_mapped_dict = dict(zip(
                    mapping_df["reactant_idx"].astype(int),
                    mapping_df["product_idx"].astype(int)
                ))

                first_shell = df["first_shell"].dropna().astype(int).tolist()

                # Walk the reaction graph to get template mappings and edge atoms
                template_mapped_dict, edge_atoms = reaction_atom_walker(
                    reactant_mol,
                    first_shell,
                    fully_mapped_dict
                )

                # Augment DataFrame with template mappings and edge atoms
                df = self._add_dict_as_new_columns(
                    df,
                    template_mapped_dict,
                    titles=["template_reactant_idx", "template_product_idx"]
                )

                df = self._add_column_safe(df, edge_atoms, "edge_atoms")
                data["reaction_dataframe"] = df.copy()

                # Build metadata object
                metadata = ReactionMetadata(
                    reaction_id=rid,
                    reactant_combined_mol=data["reactant"],
                    product_combined_mol=data["product"],
                    reactant_to_product_mapping=fully_mapped_dict,
                    template_reactant_to_product_mapping=template_mapped_dict,
                    edge_atoms=edge_atoms,
                    initiators=df["initiators"].dropna().astype(int).tolist(),
                    byproduct_indices=df["byproduct_indices"].dropna().astype(int).tolist(),
                    reaction_smarts=reaction.reaction_smarts,
                    csv_path=data["csv_path"],
                    reaction_dataframe=df,
                    delete_atom=data.get("delete_atom", True),
                )

                metadata_list.append(metadata)

        return metadata_list

    def reaction_templates_highlighted_image_grid(
        self,
        metadata_list: List[ReactionMetadata],
        highlight_type: str = "template"
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
            reactant = metadata.reactant_combined_mol
            product = metadata.product_combined_mol
            # Clear atom maps for clean visualization
            for atom in reactant.GetAtoms():
                atom.SetAtomMapNum(0)
            for atom in product.GetAtoms():
                atom.SetAtomMapNum(0)
            df = metadata.reaction_dataframe

            atoms = []
            color_map = {}

            if highlight_type == "template":
                atoms = list(metadata.template_reactant_to_product_mapping.keys())
                for a in atoms:
                    color_map[a] = (0.2, 0.6, 1.0)  # blue

            elif highlight_type == "edge":
                atoms = metadata.edge_atoms
                for a in atoms:
                    color_map[a] = (1.0, 0.4, 0.0)  # orange

            elif highlight_type == "initiators":
                atoms = df["initiators"].dropna().astype(int).tolist()
                for a in atoms:
                    color_map[a] = (0.0, 0.8, 0.2)  # green

            elif highlight_type == "delete":
                if not metadata.delete_atom:
                    atoms = []
                else:
                    atoms = metadata.byproduct_indices
                    for a in atoms:
                        color_map[a] = (1.0, 0.0, 0.0)  # red

            else:
                atoms = []

            mols.extend([reactant, product])
            highlight_lists.append(atoms)
            highlight_lists.append([])  # No highlights on products
            highlight_colors.append(color_map)
            highlight_colors.append({})

        # Draw grid: 2 columns (reactant | product)
        img = Draw.MolsToGridImage(
            mols,
            molsPerRow=2,
            highlightAtomLists=highlight_lists,
            highlightAtomColors=highlight_colors,
            subImgSize=(400, 400),
            useSVG=False
        )

        return img
