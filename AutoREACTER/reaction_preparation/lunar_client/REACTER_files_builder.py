"""
LAMMPS Molecule File Processor
==============================

This module processes LAMMPS molecule files for bond-react simulations.
It handles parsing, modification, and formatting of molecular structure data
including atoms, bonds, angles, dihedrals, and impropers.
Based on specified atom indices (template), it filters and reindexes the relevant data.
The indexes obtained from the walker module are used to identify which atoms to for REACTER

The workflow converts atom indices from template molecules to product molecules
using a mapping dictionary (reactant_to_product).

Dependencies:
    - pandas: For DataFrame operations
    - os: For file path operations
    - re: For regular expression operations
"""

import os
import pandas as pd
import datetime
from pathlib import Path
from dataclasses import dataclass
from typing import Optional
from AutoREACTER.input_parser import SimulationSetup
from AutoREACTER.reaction_preparation.lunar_client.lunar_api_wrapper import LunarFiles
from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import ReactionMetadata 
from AutoREACTER.input_parser import SimulationSetup
from AutoREACTER.reaction_preparation.lunar_client.modifiers_molecule_files import (
    modify_types, modify_charges, modify_coords,
    modify_bonds, modify_angles, modify_dihedrals, modify_impropers,
)
# from AutoREACTER.reaction_preparation.lunar_client.modifiers_data_files import (
#     modify_atoms_data,
#     modify_bonds_data, modify_angles_data,
#     modify_dihedrals_data, modify_impropers_data,
)
now = datetime.datetime.now() 
import logging
logger = logging.getLogger(__name__)

@dataclass(slots=True)
class DataFiles:
    """Simple container for paired LAMMPS data and molecule files."""
    data_file: Path             # Main *.data file for LAMMPS
    lmp_molecule_file: Path     # Associated *.lmpmol molecule template

@dataclass(slots=True)
class MoleculeFile:
    """Wrapper associating a molecule ID with its generated data files."""
    id: str
    molecule_files: Optional[DataFiles]

@dataclass(slots=True)
class TemplateFile:
    """
    Container for pre- and post-reaction template file pairs.
    """
    reaction_id: Optional[int]
    map_file: Optional[Path]
    pre_reaction_file: Optional[DataFiles]
    post_reaction_file: Optional[DataFiles]

@dataclass(slots=True)
class REACTERFiles:
    """
    Complete collection of output files from the LUNAR workflow.
    """
    force_field_data: Path             # force_field.data (FF parameters)
    in_file: Path                      # in.fix_bond_react.script (LAMMPS input)
    molecule_files: list[MoleculeFile]
    template_files: list[TemplateFile]

class REACTERFilesBuilder:
    def __init__(self, cache_dir: Path, updated_inputs_with_3d_mols: SimulationSetup):
        self.cache_dir = Path(cache_dir) / "lunar" / "REACTER_files"
        os.makedirs(self.cache_dir, exist_ok=True)
        self.force_field = updated_inputs_with_3d_mols.force_field
            

    def _get_ending_integer(self, s: str) -> int | None:
        """
        Extract and convert the trailing integer from a string.
        
        This function is useful for parsing reaction numbers from filenames
        (e.g., extracting "1" from "pre1" or "post1").
        
        Args:
            s (str): Input string that may end with an integer.
        
        Returns:
            int | None: The integer found at the end of the string, 
                    or None if no integer is found.
        
        Example:
            >>> get_ending_integer("pre1")
            1
            >>> get_ending_integer("data2")
            2
            >>> get_ending_integer("molecule")
            None
        """
        # Regular expression pattern to match one or more digits at the end of the string
        # r'\d+$' matches digits (\d+) at the end ($) of the string
        match = re.search(r'\d+$', s)
        
        if match:
            # Extract the matched substring and convert to integer
            return int(match.group())
        else:
            # No integer found at the end of the string
            return None
    
    def _ensure_dir(self, p: str, remove_blocking_file: bool = False) -> str:
        if os.path.exists(p) and not os.path.isdir(p):
            if not remove_blocking_file:
                raise FileExistsError(f"Path exists and is a file (expected dir): {p}")
            logger.warning("ensure_dir: removing blocking file: %s", p)
            os.remove(p)
        os.makedirs(p, exist_ok=True)
        return p

    def _col_int_list(self, col: str, df: pd.DataFrame) -> list[int]:
        if col not in df.columns:
            return []
        s = df[col].dropna()
        if s.empty:
            return []
        # preserve order while removing duplicates
        seen = set()
        vals = []
        for x in s.tolist():
            try:
                fx = float(x)
            except (TypeError, ValueError):
                raise ValueError(f"{col} contains non-numeric value: {x!r}")
            if not fx.is_integer():
                raise ValueError(f"{col} contains non-integer value: {x!r}")
            ix = int(fx)
            if ix not in seen:
                seen.add(ix)
                vals.append(ix)
        return vals

    def _load_molecule_file(self, molecule_file_path: Path)-> tuple[list[str], int, int, int, int, int, int, int]:
        """
        Load and parse a molecule configuration file to extract structural data sections.
        
        This function reads a molecule file and identifies the starting indices for various
        structural sections (Types, Charges, Coords, Bonds, Angles, Dihedrals, Impropers).
        The file is expected to have a specific format with section headers followed by data.
        
        Args:
            molecule_file_path (Path): Path to the molecule configuration file to be loaded.
        
        Returns:
            tuple: A tuple containing:
                - lines (list): All lines from the file with whitespace stripped
                - type_start_index (int): Starting line index for the Types section
                - charge_start_index (int): Starting line index for the Charges section
                - coord_start_index (int): Starting line index for the Coords section
                - bond_start_index (int): Starting line index for the Bonds section
                - angle_start_index (int): Starting line index for the Angles section
                - dihedral_start_index (int): Starting line index for the Dihedrals section
                - improper_start_index (int): Starting line index for the Impropers section
        
        Raises:
            ValueError: If the file is not found or essential sections (Types, Charges, Coords)
                    are missing from the file.
        
        """
        # Open and read the molecule file
        with open(molecule_file_path, 'r') as file:
            if file is None:
                raise ValueError(f"File {molecule_file_path} not found.")
            # Strip whitespace from each line for cleaner processing
            lines = [line.strip() for line in file]
        
        # Initialize section indices to None (will be set when sections are found)
        type_start_index = None
        charge_start_index = None
        coord_start_index = None
        bond_start_index = None
        angle_start_index = None
        dihedral_start_index = None
        improper_start_index = None
        
        # Iterate through lines to find section headers and set their starting indices
        for index, line in enumerate(lines):
            if line.strip() == "Types":
                type_start_index = index + 2
            if line.strip() == "Charges":
                charge_start_index = index + 2
            if line.strip() == "Coords":
                coord_start_index = index + 2
            if line.strip() == "Bonds":
                bond_start_index = index + 2
            if line.strip() == "Angles":
                angle_start_index = index + 2
            if line.strip() == "Dihedrals":
                dihedral_start_index = index + 2
            if line.strip() == "Impropers":
                improper_start_index = index + 2
        
        # Validate that essential sections were found
        if type_start_index is None or charge_start_index is None or coord_start_index is None:
            raise ValueError("Essential sections (Types, Charges, Coords) not found in the file.")
        
        return lines, type_start_index, charge_start_index, coord_start_index, bond_start_index, angle_start_index, dihedral_start_index, improper_start_index


    


    def _molecule_file_format(self, number_of_types: int, number_of_bonds: int, number_of_angles: int, 
                            number_of_dihedrals: int, number_of_impropers: int, types_section: str, charge_section: str, coord_section: str, 
                            bond_section: str, angle_section: str, dihedral_section: str, improper_section: str)-> str:
        """
        Format LUNAR molecule data into LAMMPS REACTER molecule file format.
        
        Constructs a complete molecule file string following LAMMPS format
        with header and all sections (types, charges, coordinates, bonds, angles,
        dihedrals, and impropers).
        
        Parameters:
        -----------
        number_of_types : int
            Total number of atom types in the molecule.
        number_of_bonds : int
            Total number of bonds in the molecule.
        number_of_angles : int
            Total number of angles in the molecule.
        number_of_dihedrals : int
            Total number of dihedrals in the molecule.
        number_of_impropers : int
            Total number of impropers in the molecule.
        types_section : str
            Formatted string containing atom type information.
        charge_section : str
            Formatted string containing atomic charges.
        coord_section : str
            Formatted string containing atomic coordinates.
        bond_section : str
            Formatted string containing bond information.
        angle_section : str
            Formatted string containing angle information.
        dihedral_section : str
            Formatted string containing dihedral information.
        improper_section : str
            Formatted string containing improper dihedral information.
        
        Returns:
        --------
        str
            Complete formatted molecule file content ready for output.
        
        Notes:
        ------
        - Header includes version information and atom typing method
        - Class 2 force field format is used
        - All sections are clearly labeled with headers
        """
        # Create molecule file with header and count information
        modified_molecule_file = f"""HEADER,  pre1.mol  read w/ mol2lmp > atom_typing v1.11 / 14 April 2025 w/PCFF-IFF atom-types > all2lmp: v1.22 / 14 April 2025  Class: 2 > bond_react_merge: v1.17 / 14 April 2025 molecule file (filetag: pre1)

    {number_of_types} atoms
    {number_of_bonds} bonds
    {number_of_angles} angles
    {number_of_dihedrals} dihedrals
    {number_of_impropers} impropers

    Types

    """
        # Append all data sections in order
        modified_molecule_file += types_section
        modified_molecule_file += "\nCharges\n\n"
        modified_molecule_file += charge_section
        modified_molecule_file += "\nCoords\n\n"
        modified_molecule_file += coord_section
        modified_molecule_file += "\nBonds\n\n"
        modified_molecule_file += bond_section
        modified_molecule_file += "\nAngles\n\n"
        modified_molecule_file += angle_section
        modified_molecule_file += "\nDihedrals\n\n"
        modified_molecule_file += dihedral_section
        modified_molecule_file += "\nImpropers\n\n"
        modified_molecule_file += improper_section
        
        return modified_molecule_file


    def _molecule_file_preparation(self, template_file_path: Path, template_indexes: list[int]) -> str:
        """
        Orchestrate the complete molecule file processing pipeline.
        
        This function coordinates all steps needed to convert a template molecule
        file to a product molecule file by parsing sections, modifying atom indices,
        and reformatting the data.
        
        Parameters:
        -----------
        template_file_path : Path
            File path to the input molecule file to be processed.
        template_indexes : list[int]
            List of atom indices to use for mapping (1-based indexing).
            These are the atoms that will be tracked through the reaction.
        
        Returns:
        --------
        str
            Complete formatted molecule file content ready for output.
        
        Workflow:
        ---------
        1. Load and parse the input molecule file
        2. Extract and modify types section with template indices
        3. Modify charges based on new atom mappings
        4. Modify coordinates based on new atom mappings
        5. Modify bonds with new atom indices
        6. Modify angles with new atom indices
        7. Modify dihedrals with new atom indices
        8. Modify impropers with new atom indices
        9. Format all sections into complete molecule file
        
        Notes:
        ------
        - Depends on external functions: load_molecule_file, modify_* functions
        - Uses type_df returned by modify_types for all subsequent mappings
        """
        # Parse molecule file and get starting indices for each section
        lines, \
            type_start_index, \
            charge_start_index, \
                coord_start_index, \
                    bond_start_index, \
                        angle_start_index, \
                            dihedral_start_index, \
                                improper_start_index = self._load_molecule_file(template_file_path)
        
        # Process types section - returns DataFrame for mapping and formatted string
        df_types, types_section, number_of_types, index_change_dict = modify_types(lines, template_indexes, type_start_index)
        
        # Process charges section using atom mapping from types
        charge_section = modify_charges(lines, df_types, charge_start_index)
        
        # Process coordinates section using atom mapping from types
        coord_section = modify_coords(lines, df_types, coord_start_index)
        
        # Process bonds section and get count
        bond_section, number_of_bonds = modify_bonds(lines, df_types, bond_start_index, legacy_mode=False)
        
        # Process angles section and get count
        angle_section, number_of_angles = modify_angles(lines, df_types, angle_start_index, legacy_mode=False)
        
        # Process dihedrals section and get count
        dihedral_section, number_of_dihedrals = modify_dihedrals(lines, df_types, dihedral_start_index, legacy_mode=False)
        
        # Process impropers section and get count
        improper_section, number_of_impropers = modify_impropers(lines, df_types, improper_start_index)
        
        # Assemble all sections into final molecule file format
        modified_molecule_file = self._molecule_file_format(number_of_types, number_of_bonds, number_of_angles, 
                                                    number_of_dihedrals, number_of_impropers, types_section, charge_section, coord_section, 
                                                    bond_section, angle_section, dihedral_section, improper_section)
        
        return modified_molecule_file, index_change_dict

    def _map_file_write(self, reactant_to_product, initiator_atoms, edge_atoms, delete_ids):
        """
        Construct a textual ".map" file describing a superimposition mapping
        between reactant and product atom indices.

        The textual format constructed here is a simple, human-readable structure
        expected downstream by other tools in the workflow. The important details:

        - reactant_to_product: dict mapping reactant_index (0-based) -> product_index (0-based)
        These indices refer to template atom indices after reindexing/filtering; the
        entries will be written as 1-based values in the map file.
        - initiator atoms: list of two atom indices (0-based) that define the initiator pair.
        - edge_atoms: list of atom indices (0-based) that are considered edge atoms.
        - delete_ids: list of atom indices (0-based) to mark as deleted; an empty or falsy
        value disables the DeleteIDs section.

        Returns:
        A single string containing the contents of the map file.
        """
        # Header giving a nominal description
        map_file = "this is a nominal superimpose file\n\n"

        # Write counts: number of edge IDs and number of equivalences
        # (len(reactant_to_product) is the number of mapping pairs)
        map_file += f"""{len(edge_atoms)} edgeIDs
    {len(reactant_to_product)} equivalences"""

        # Optionally include a deleteIDs comment line that indicates how many delete IDs
        if delete_ids:
            map_file += f"\n#{len(delete_ids)} deleteIDs\n"
        else:
            map_file += "\n"

        # The map format expects 1-based indices for initiators; the code callers supply
        # 0-based indices, so add 1 for output. The initiator list is expected to contain
        # exactly two atoms; we write them on separate lines.
        # The blank lines before and after are also part of the format and are
        # preserved to avoid changing function signature.
        if len(initiator_atoms) != 2:
            raise ValueError(
                f"Expected exactly 2 initiator atoms (0-based) after filtering, got {len(initiator_atoms)}: {initiator_atoms}"
            )
        map_file += f"\nInitiatorIDs\n\n{initiator_atoms[0]+1}\n{initiator_atoms[1]+1}\n\nEdgeIDs\n\n"
        # Edge atom indices are output as 1-based, one per line
        for atom in edge_atoms:
            map_file += f"{atom+1}\n"

        # Equivalences section: write reactant_index  product_index pairs.
        map_file += "\nEquivalences\n\n"
        for key, value in reactant_to_product.items():
            # key, value are 0-based internally; write as 1-based
            map_file += f"{(key + 1):<5} {value + 1}\n"

        # If delete IDs are present, output them prefixed with '#' on their own lines.
        if delete_ids:
            map_file += "\n#DeleteIDs\n\n"
            for atom in delete_ids:
                map_file += f"#{atom+1}\n"

        return map_file


    def _build_bond_react_templates(self,
        file_dict: dict,
        reactant_to_product: dict,
        initiator_atoms: list,
        edge_atoms: list,
        delete_ids: list,
    ):
        """
        Process pairs of 'pre' and 'post' molecule files and generate reindexed
        template files plus a .map file for each reaction.

        Parameters
        ----------
        file_dict : dict
        Mapping of filenames to file paths for the current reaction context.
        Expected keys look like "pre_{num}" and "post_{num}".
        reactant_to_product : dict
        Mapping from reactant template atom indices (0-based) -> product template atom indices (0-based).
        These indices refer to the full original molecule file indices before filtering/reindexing.
        initiator_atoms : list[int]
        List (typically of length 2) of atom indices (0-based) in the reactant template that act as initiators.
        edge_atoms : list[int]
        List of atom indices (0-based) that represent the "edge" atoms in the template.
        delete_ids : list[int]
        List of atom indices (0-based) to be considered deleted (byproducts) if present.
        cache_dir_reactor : str
        Directory path where reactor-related files should be cached. Output files are written directly
        under this directory.

        Returns
        -------
        molecule_file_dict : dict
        A dictionary mapping generated logical names to filesystem paths, e.g.:
            {
            "pre_{num}": "<path to reindexed pre file>",
            "post_{num}": "<path to reindexed post file>",
            "map_file_{num}": "<path to generated .map file>",
            ...
            }

        Notes
        -----
        - This function relies on an external helper `molecule_file_preparation` which is
        expected to return a tuple (modified_contents, index_change_dict). The index_change_dict
        should map original 1-based indices to new 1-based indices after filtering.
        - Conversions between 0-based and 1-based indexing are handled carefully and
        documented inline.
        """
        template_indexes_reactant = []
        template_indexes_product = []


        # Build 1-based index lists for filtering: molecule_file_preparation expects
        # 1-based indices (hence +1 conversion from reactant_to_product keys/values).
        for key, value in reactant_to_product.items():
            template_indexes_reactant.append(int(key + 1))
            template_indexes_product.append(int(value + 1))

        # Iterate the provided file dictionary and only process entries that start with "pre_"
        # so that we can find their matching "post_{num}" counterparts.
        for name, pre_path in file_dict.items():
            if not name.startswith("pre_"):
                continue

            # Extract the trailing integer to form the matching post key
            num = self._get_ending_integer(name)
            post_key = f"post_{num}"
            post_path = file_dict.get(post_key)
            
            if post_path is None:
                # If no matching post file is present, this is a critical error for building templates.
                raise ValueError(f"Corresponding post-reaction file for pre_{num} not found.")

            # Reindex and filter the molecule files so they only include atoms relevant to the templates.
            # Each call returns new file contents and an index_change_dict mapping original 1-based indices
            # -> new 1-based indices.
            if not pre_path or not os.path.isfile(pre_path):
                raise FileNotFoundError(f"Missing pre file: {pre_path}")
            if not post_path or not os.path.isfile(post_path):
                raise FileNotFoundError(f"Missing post file: {post_path}")
            
            pre_modified, index_change_dict_reactant = self._molecule_file_preparation(
                pre_path,
                template_indexes_reactant
            )

            post_modified, index_change_dict_product =  self._molecule_file_preparation(
                post_path,
                template_indexes_product
            )

            # Write the modified molecule files to the cache directory's parent (consistent with previous code)
            os.makedirs(self.cache_dir, exist_ok=True)

            pre_out = os.path.join(self.cache_dir, f"template_pre_{num}.molecule")
            with open(pre_out, "w") as f:
                f.write(pre_modified)

            post_out = os.path.join(self.cache_dir, f"template_post_{num}.molecule")
            with open(post_out, "w") as f:
                f.write(post_modified)

            # Record the generated file paths in the return dictionary


            # Build equivalence mapping in 0-based template index space.
            # reactant_to_product maps original 0-based full-file indices. We need to map those
            # to indices in the reindexed template files. The index_change_dicts map original
            # 1-based indices -> new 1-based indices, so we convert accordingly.
            template_reactant_to_template_product = {}
            missing_pairs = []  # (r_full0, p_full0, r_missing, p_missing)
            
            for r_full0, p_full0 in reactant_to_product.items():
                r_key = r_full0 + 1  # full -> 1-based for index_change_dict keys
                p_key = p_full0 + 1
            
                r_missing = r_key not in index_change_dict_reactant
                p_missing = p_key not in index_change_dict_product
            
                if r_missing or p_missing:
                    missing_pairs.append((r_full0, p_full0, r_missing, p_missing))
                    continue
            
                r_tmp = index_change_dict_reactant[r_key]  # 1-based template index
                p_tmp = index_change_dict_product[p_key]
                # Convert back to 0-based template indices in the map that will be passed to map_file_write
                template_reactant_to_template_product[r_tmp - 1] = p_tmp - 1
            
            if missing_pairs:
                N = 20
                preview = missing_pairs[:N]
                r_miss_count = sum(1 for *_, rm, _ in missing_pairs if rm)
                p_miss_count = sum(1 for *_, _, pm in missing_pairs if pm)
            
                raise ValueError(
                    "Template mapping indices missing after filtering. "
                    "reactant_to_product contains atoms removed during template trimming.\n"
                    f"Missing reactant mappings: {r_miss_count}; missing product mappings: {p_miss_count}; "
                    f"total missing pairs: {len(missing_pairs)}.\n"
                    f"First {min(N, len(missing_pairs))} missing (r_full0, p_full0, r_missing, p_missing): {preview}"
                )
            
            # Recompute initiator, delete, and edge atom lists for the reindexed template.
            # Only include atoms that survived the filtering process (i.e., have an entry in the index_change_dict).
            initiator_atoms_t = [
                index_change_dict_reactant[i + 1] - 1
                for i in initiator_atoms
                if (i + 1) in index_change_dict_reactant
            ]
            delete_ids_t = [
                index_change_dict_reactant[i + 1] - 1
                for i in delete_ids
                if (i + 1) in index_change_dict_reactant
            ]
            edge_atoms_t = [
                index_change_dict_reactant[i + 1] - 1
                for i in edge_atoms
                if (i + 1) in index_change_dict_reactant
            ]


            # Construct the .map file contents using the reindexed template-space mappings
            map_file = self._map_file_write(
                template_reactant_to_template_product,
                initiator_atoms_t,
                edge_atoms_t,
                delete_ids_t
            )

            map_path = os.path.join(self.cache_dir, f"RXN_{num}.map")
            with open(map_path, "w") as f:
                f.write(map_file)

            

        return pre_out, post_out, map_path


    def molecule_template_preparation(self, lunar_files: LunarFiles, prepared_reactions_with_3d_mols: list[ReactionMetadata], updated_inputs_with_3d_mols: SimulationSetup) -> REACTERFiles:
        """
        Top-level orchestrator that loops over reactions and prepares template
        files and mappings for each reaction.

        Parameters
        ----------
        lunar_files : LunarFiles
        Container for all input files and paths needed for the preparation process.
        prepared_reactions_with_3d_mols : list[ReactionMetadata]
        List of ReactionMetadata objects, each containing information about a reaction,
            including the reaction DataFrame and any relevant flags (e.g., delete_atoms).
        updated_inputs_with_3d_mols : SimulationSetup

        Returns
        -------
        REACTERFiles
        A complete collection of output files generated from the preparation process,
        """
        template_files = []
        pre_and_post_files = lunar_files.template_files
        
        # Iterate each reaction entry and build templates for that single reaction only.
        for rxn in pre_and_post_files:
            id = rxn.reaction_id
            print(f"Preparing templates and map file for reaction ID: {id}")
            # The reaction DataFrame is expected under the key "reaction_dataframe"
            # Find the specific metadata object in the list that matches the current ID
            current_rxn_metadata = next(m for m in prepared_reactions_with_3d_mols if m.reaction_id == id)
            df: pd.DataFrame = current_rxn_metadata.reaction_dataframe
            # Optional indicator whether byproducts listed should be treated as delete IDs
            delete_atoms = getattr(current_rxn_metadata, "delete_atoms", False)

            # Extract integer lists from DataFrame columns using helper _col_int_list
            byproducts = self._col_int_list("byproduct_indices", df = df)
            initiators = self._col_int_list("initiators", df = df)
            edge_atoms = self._col_int_list("edge_atoms", df = df)

            # Build a template-level mapping (if present in DataFrame). Keys/values are expected
            # to be 0-based indices in the original CSV context.
            template_map = {}
            if {"template_reactant_idx", "template_product_idx"}.issubset(df.columns):
                tmp = df[["template_reactant_idx", "template_product_idx"]].dropna()
                for r_idx, p_idx in tmp.itertuples(index=False):
                    template_map[int(r_idx)] = int(p_idx)

            # If delete_atoms flag is set, delete_ids come from byproducts; otherwise none.
            delete_ids = byproducts if delete_atoms else []

            # The current reaction's pre and post molecule files are expected to be found in lunar_files.template_files with keys like "pre_{id}" and "post_{id}".
            current_rxn_files = {
                f"pre_{id}": rxn.pre_reaction_file.lmp_molecule_file,
                f"post_{id}": rxn.post_reaction_file.lmp_molecule_file
            }

            # Build the per-reaction template files and mappings.
            pre_out, post_out, map_path = self._build_bond_react_templates(
                file_dict = current_rxn_files,
                reactant_to_product = template_map,
                initiator_atoms = initiators,
                edge_atoms = edge_atoms,
                delete_ids = delete_ids
            )
            template_files.append(TemplateFile(
                reaction_id = id,
                map_file = Path(map_path),
                pre_reaction_file = DataFiles(
                    data_file = None,  # Placeholder, as the current workflow focuses on molecule files
                    lmp_molecule_file = Path(pre_out)
                ),
                post_reaction_file = DataFiles(
                    data_file = None,  # Placeholder
                    lmp_molecule_file = Path(post_out)
                )
            ))

        return REACTERFiles(template_files=template_files)

