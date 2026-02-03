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

import shutil
import pandas as pd
import re
import os
import json
import logging
logger = logging.getLogger(__name__)

def get_ending_integer(s: str) -> int | None:
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
    
def ensure_dir(p: str, remove_blocking_file: bool = False) -> str:
    if os.path.exists(p) and not os.path.isdir(p):
        if not remove_blocking_file:
            raise FileExistsError(f"Path exists and is a file (expected dir): {p}")
        logger.warning("ensure_dir: removing blocking file: %s", p)
        os.remove(p)
    os.makedirs(p, exist_ok=True)
    return p

def _col_int_list(col: str, df: pd.DataFrame) -> list[int]:
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

def load_molecule_file(molecule_file):
    """
    Load and parse a molecule configuration file to extract structural data sections.
    
    This function reads a molecule file and identifies the starting indices for various
    structural sections (Types, Charges, Coords, Bonds, Angles, Dihedrals, Impropers).
    The file is expected to have a specific format with section headers followed by data.
    
    Args:
        molecule_file (str): Path to the molecule configuration file to be loaded.
    
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
    with open(molecule_file, 'r') as file:
        if file is None:
            raise ValueError(f"File {molecule_file} not found.")
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


def modify_types(lines, template_indexes, type_start_index):
    """
    Extract, filter, and reindex atom type information from the molecule file.
    
    This function processes the Types section of a molecule file, keeping only atoms
    whose indices are in the template_indexes list, and reassigns their indices sequentially.
    The modified data is saved to a file and formatted as a string section.
    
    Args:
        lines (list): All lines from the molecule file.
        template_indexes (list): List of original atom indices to keep in the output.
        type_start_index (int): Starting line index of the Types section.
    
    Returns:
        tuple: A tuple containing:
            - df (pd.DataFrame): DataFrame with filtered and reindexed atom type data
            - types_section (str): Formatted string representation of the types section
            - number_of_types (int): Total number of atom types in the filtered data
    
    Note:
        - Atoms not in template_indexes are excluded from the output
        - Atom indices are reassigned sequentially starting from 1
    """
    # Sort template indices for consistent processing
    template_indexes.sort()
    
    # Initialize list to store parsed atom type data
    data = []
    
    # Parse the Types section from the input lines
    for line in lines[type_start_index:]:
        # Stop processing when an empty line is encountered
        if not line.strip():
            break
        
        # Split the line into individual components
        parts = line.split()
        
        # Validate that the line has sufficient data fields
        if len(parts) >= 4:
            # Extract atom properties from the parsed line
            atom_index = int(parts[0])
            atom_type = int(parts[1])
            hash_value = parts[2]
            atom_type_real = parts[3]
            
            # Append parsed atom data to the list
            data.append({
                "atom_index": atom_index,
                "atom_type": atom_type,
                "hash": hash_value,
                "atom_type_real": atom_type_real
            })

    # Create a DataFrame from the parsed data
    df = pd.DataFrame(data)
    
    # Initialize a new column for reindexed atom indices
    df["new_atom_index"] = None 
    
    # Map original atom indices to new indices for atoms in template_indexes
    for value in df["atom_index"]:
        if value in template_indexes:
            new_type = value
            df.loc[df["atom_index"] == value, "new_atom_index"] = new_type
    
    # Identify rows where atoms are not in the template (to be removed)
    indices_to_drop = []        
    for index, row in df.iterrows():
        if row["new_atom_index"] is None:
            indices_to_drop.append(index)
    
    # Remove rows for atoms not in the template and reset the index
    df = df.drop(indices_to_drop).reset_index(drop=True)        
    
    # Reassign new atom indices sequentially starting from 1
    for i, value in enumerate(df["new_atom_index"]):
        df.at[i, "new_atom_index"] = i + 1
    
    # Format the types section as a string with proper spacing
    types_section = ""
    for i, row in df.iterrows():
        types_section += f"{row['new_atom_index']:>4}{row['atom_type']:>4}{row['hash']:>4}{row['atom_type_real']:>4}\n"
    
    # Calculate the total number of atom types
    number_of_types = len(df)

    index_change_dict = {}
    for i, row in df.iterrows():
        if row['new_atom_index']:
            index_change_dict[row['atom_index']] = row['new_atom_index']
    
    return df, types_section, number_of_types, index_change_dict



def modify_charges(lines, type_df, charge_start_index):
    """
    Extract and reindex atomic charge information based on filtered atom types.
    
    This function processes the Charges section of a molecule file, keeping only charges
    for atoms that exist in the filtered type_df DataFrame, and reassigns their indices
    to match the new atom indices from the types modification.
    
    Args:
        lines (list): All lines from the molecule file.
        type_df (pd.DataFrame): DataFrame containing filtered atom types with new indices.
        charge_start_index (int): Starting line index of the Charges section.
    
    Returns:
        str: Formatted string representation of the charges section with proper spacing.
    
    Note:
        - Only charges for atoms present in type_df are retained
        - Atom indices are updated to match the new indices from type_df
    """
    # Initialize list to store parsed charge data
    charge_data = []
    
    # Parse the Charges section from the input lines
    for line in lines[charge_start_index:]:
        # Stop processing when an empty line is encountered
        if not line.strip():
            break

        # Split the line into individual components
        parts = line.split()
        
        # Validate that the line has sufficient data fields
        if len(parts) >= 3:
            # Extract charge properties from the parsed line
            atom_index = int(parts[0])
            charge_value = float(parts[1])
            hash_value = parts[2]
            atom_type_real = parts[3]
            
            # Append parsed charge data to the list
            charge_data.append({
                "atom_index": atom_index,
                "charge_value": charge_value,
                "hash": hash_value,
                "atom_type_real": atom_type_real
            })

    # Create a DataFrame from the parsed charge data
    charge_df = pd.DataFrame(charge_data)
    
    # Initialize a new column for reindexed atom indices
    charge_df["new_atom_index"] = None 
    
    # Map original atom indices to new indices using the type_df mapping
    for value in charge_df["atom_index"]:
        if value in type_df["atom_index"].values:
            # Look up the new index from the type_df
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            charge_df.loc[charge_df["atom_index"] == value, "new_atom_index"] = new_index
    
    # Identify rows where atoms are not in the filtered type data (to be removed)
    indices_to_drop = []        
    for index, row in charge_df.iterrows():
        if row["new_atom_index"] is None:
            indices_to_drop.append(index)
    
    # Remove rows for atoms not in the filtered data and reset the index
    charge_df = charge_df.drop(indices_to_drop).reset_index(drop=True)        
    
    # Reassign new atom indices sequentially starting from 1
    for i, value in enumerate(charge_df["new_atom_index"]):
        charge_df.at[i, "new_atom_index"] = i + 1
    
    # Format the charges section as a string with proper spacing
    charge_section = ""
    for i, row in charge_df.iterrows():
        charge_section += f"{row['new_atom_index']:>4}{row['charge_value']:>12.6f}{row['hash']:>4}{row['atom_type_real']:>4}\n"
    
    return charge_section


def modify_coords(lines, type_df, coord_start_index):
    """
    Extract and reindex atomic coordinate information based on filtered atom types.
    
    This function processes the Coords section of a molecule file, keeping only coordinates
    for atoms that exist in the filtered type_df DataFrame, and reassigns their indices
    to match the new atom indices from the types modification.
    
    Args:
        lines (list): All lines from the molecule file.
        type_df (pd.DataFrame): DataFrame containing filtered atom types with new indices.
        coord_start_index (int): Starting line index of the Coords section.
    
    Returns:
        str: Formatted string representation of the coordinates section with proper spacing.
    
    Note:
        - Only coordinates for atoms present in type_df are retained
        - Atom indices are updated to match the new indices from type_df
        - Coordinates are formatted to 6 decimal places
    """
    # Initialize list to store parsed coordinate data
    coord_data = []
    
    # Parse the Coords section from the input lines
    for line in lines[coord_start_index:]:
        # Stop processing when an empty line is encountered
        if not line.strip():
            break

        # Split the line into individual components
        parts = line.split()
        
        # Validate that the line has sufficient data fields
        if len(parts) >= 5:
            # Extract coordinate properties from the parsed line
            atom_index = int(parts[0])
            x = float(parts[1])
            y = float(parts[2])
            z = float(parts[3])
            hash_value = parts[4]
            atom_type_real = parts[5]
            
            # Append parsed coordinate data to the list
            coord_data.append({
                "atom_index": atom_index,
                "x": x,
                "y": y,
                "z": z,
                "hash": hash_value,
                "atom_type_real": atom_type_real
            })

    # Create a DataFrame from the parsed coordinate data
    coord_df = pd.DataFrame(coord_data)
    
    # Initialize a new column for reindexed atom indices
    coord_df["new_atom_index"] = None 
    
    # Map original atom indices to new indices using the type_df mapping
    for value in coord_df["atom_index"]:
        if value in type_df["atom_index"].values:
            # Look up the new index from the type_df
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            coord_df.loc[coord_df["atom_index"] == value, "new_atom_index"] = new_index
    
    # Identify rows where atoms are not in the filtered type data (to be removed)
    indices_to_drop = []        
    for index, row in coord_df.iterrows():
        if row["new_atom_index"] is None:
            indices_to_drop.append(index)
    
    # Remove rows for atoms not in the filtered data and reset the index
    coord_df = coord_df.drop(indices_to_drop).reset_index(drop=True)        
    
    # Reassign new atom indices sequentially starting from 1
    for i, value in enumerate(coord_df["new_atom_index"]):
        coord_df.at[i, "new_atom_index"] = i + 1
    
    # Format the coordinates section as a string with proper spacing
    coord_section = ""
    for i, row in coord_df.iterrows():
        coord_section += f"{row['new_atom_index']:>4}{row['x']:>12.6f}{row['y']:>12.6f}{row['z']:>12.6f}{row['hash']:>4}{row['atom_type_real']:>4}\n"
    
    return coord_section


def modify_bonds(lines, type_df, bond_start_index):
    """
    Extract and reindex bond information based on filtered atom types.
    
    This function processes the Bonds section of a molecule file, keeping only bonds
    where both atoms exist in the filtered type_df DataFrame. Bond indices and atom
    indices are reassigned to match the new atom indices from the types modification.
    
    Args:
        lines (list): All lines from the molecule file.
        type_df (pd.DataFrame): DataFrame containing filtered atom types with new indices.
        bond_start_index (int): Starting line index of the Bonds section.
    
    Returns:
        tuple: A tuple containing:
            - bond_section (str): Formatted string representation of the bonds section
            - number_of_bonds (int): Total number of bonds in the filtered data
    
    Note:
        - Only bonds where both atoms are present in type_df are retained
        - Bond indices are reassigned sequentially starting from 1
        - Atom indices are updated to match the new indices from type_df
    """
    # Initialize list to store parsed bond data
    bond_data = []
    
    # Parse the Bonds section from the input lines
    for line in lines[bond_start_index:]:
        # Stop processing when an empty line is encountered
        if not line.strip():
            break

        # Split the line into individual components
        parts = line.split()
        
        # Validate that the line has sufficient data fields
        if len(parts) >= 7:
            # Extract bond properties from the parsed line
            bond_index = int(parts[0])
            bond_type = int(parts[1])
            atom1_index = int(parts[2])
            atom2_index = int(parts[3])
            hash_value = parts[4]
            atom1_type_real = parts[5]
            atom2_type_real = parts[6]
            
            # Append parsed bond data to the list
            bond_data.append({
                "bond_index": bond_index,
                "bond_type": bond_type,
                "atom1_index": atom1_index,
                "atom2_index": atom2_index,
                "hash": hash_value,
                "atom1_type_real": atom1_type_real,
                "atom2_type_real": atom2_type_real
            })

    # Create a DataFrame from the parsed bond data
    bond_df = pd.DataFrame(bond_data)
    
    # Initialize new columns for reindexed atom indices
    bond_df["new_atom1_index"] = None 
    bond_df["new_atom2_index"] = None 
    
    # Map original atom1 indices to new indices using the type_df mapping
    for value in bond_df["atom1_index"]:
        if value in type_df["atom_index"].values:
            # Look up the new index from the type_df
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            bond_df.loc[bond_df["atom1_index"] == value, "new_atom1_index"] = new_index
    
    # Map original atom2 indices to new indices using the type_df mapping
    for value in bond_df["atom2_index"]:
        if value in type_df["atom_index"].values:
            # Look up the new index from the type_df
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            bond_df.loc[bond_df["atom2_index"] == value, "new_atom2_index"] = new_index
    
    # Identify rows where either atom is not in the filtered type data (to be removed)
    indices_to_drop = []        
    for index, row in bond_df.iterrows():
        if row["new_atom1_index"] is None or row["new_atom2_index"] is None:
            indices_to_drop.append(index)
    
    # Remove rows for bonds with missing atoms and reset the index
    bond_df = bond_df.drop(indices_to_drop).reset_index(drop=True)        
    
    # Reassign bond indices sequentially starting from 1
    for i, value in enumerate(bond_df.index):
        bond_df.at[i, "bond_index"] = i + 1
    
    # Format the bonds section as a string with proper spacing
    bond_section = ""
    for i, row in bond_df.iterrows():
        bond_section += f"{row['bond_index']:>4}{row['bond_type']:>4}{row['new_atom1_index']:>4}{row['new_atom2_index']:>4}{row['hash']:>4}{row['atom1_type_real']:>4}{row['atom2_type_real']:>4}\n"
    
    # Calculate the total number of bonds
    number_of_bonds = len(bond_df)
    
    return bond_section, number_of_bonds

def modify_angles(lines, type_df, angle_start_index):
    """
    Extract and reindex angle information based on filtered atom types.
    
    This function processes the Angles section of a molecule file, keeping only angles
    where all three atoms exist in the filtered type_df DataFrame. Angle indices and atom
    indices are reassigned to match the new atom indices from the types modification.
    
    Args:
        lines (list): All lines from the molecule file.
        type_df (pd.DataFrame): DataFrame containing filtered atom types with new indices.
        angle_start_index (int): Starting line index of the Angles section.
    
    Returns:
        tuple: A tuple containing:
            - angle_section (str): Formatted string representation of the angles section
            - number_of_angles (int): Total number of angles in the filtered data
    
    Note:
        - Only angles where all three atoms are present in type_df are retained
        - Angle indices are reassigned sequentially starting from 1
        - Atom indices are updated to match the new indices from type_df
    """
    # Initialize list to store parsed angle data
    angle_data = []
    
    # Parse the Angles section from the input lines
    for line in lines[angle_start_index:]:
        # Stop processing when an empty line is encountered
        if not line.strip():
            break
        
        # Split the line into individual components
        parts = line.split()
        
        # Validate that the line has sufficient data fields
        if len(parts) >= 9:
            # Extract angle properties from the parsed line
            angle_index = int(parts[0])
            angle_type = int(parts[1])
            atom1_index = int(parts[2])
            atom2_index = int(parts[3])
            atom3_index = int(parts[4])
            hash_value = parts[5]
            atom1_type_real = parts[6]
            atom2_type_real = parts[7]
            atom3_type_real = parts[8]
            
            # Append parsed angle data to the list
            angle_data.append({
                "angle_index": angle_index,
                "angle_type": angle_type,
                "atom1_index": atom1_index,
                "atom2_index": atom2_index,
                "atom3_index": atom3_index,
                "hash": hash_value,
                "atom1_type_real": atom1_type_real,
                "atom2_type_real": atom2_type_real,
                "atom3_type_real": atom3_type_real
            })
    
    # Create a DataFrame from the parsed angle data
    angle_df = pd.DataFrame(angle_data)
    
    # Initialize new columns for reindexed atom indices
    angle_df["new_atom1_index"] = None 
    angle_df["new_atom2_index"] = None 
    angle_df["new_atom3_index"] = None 
    
    # Map original atom1 indices to new indices using the type_df mapping
    for value in angle_df["atom1_index"]:
        if value in type_df["atom_index"].values:
            # Look up the new index from the type_df
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            angle_df.loc[angle_df["atom1_index"] == value, "new_atom1_index"] = new_index
    
    # Map original atom2 indices to new indices using the type_df mapping
    for value in angle_df["atom2_index"]:
        if value in type_df["atom_index"].values:
            # Look up the new index from the type_df
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            angle_df.loc[angle_df["atom2_index"] == value, "new_atom2_index"] = new_index
    
    # Map original atom3 indices to new indices using the type_df mapping
    for value in angle_df["atom3_index"]:
        if value in type_df["atom_index"].values:
            # Look up the new index from the type_df
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            angle_df.loc[angle_df["atom3_index"] == value, "new_atom3_index"] = new_index
    
    # Identify rows where any atom is not in the filtered type data (to be removed)
    indices_to_drop = []        
    for index, row in angle_df.iterrows():
        if (row["new_atom1_index"] is None or 
            row["new_atom2_index"] is None or 
            row["new_atom3_index"] is None):
            indices_to_drop.append(index)
    
    # Remove rows for angles with missing atoms and reset the index
    angle_df = angle_df.drop(indices_to_drop).reset_index(drop=True)        
    
    # Reassign angle indices sequentially starting from 1
    for i, value in enumerate(angle_df.index):
        angle_df.at[i, "angle_index"] = i + 1
    
    # Format the angles section as a string with proper spacing
    angle_section = ""
    for i, row in angle_df.iterrows():
        angle_section += f"{row['angle_index']:>4}{row['angle_type']:>4}{row['new_atom1_index']:>4}{row['new_atom2_index']:>4}{row['new_atom3_index']:>4}{row['hash']:>4}{row['atom1_type_real']:>4}{row['atom2_type_real']:>4}{row['atom3_type_real']:>4}\n"
    
    # Calculate the total number of angles
    number_of_angles = len(angle_df)
    
    return angle_section, number_of_angles

def modify_dihedrals(lines, type_df, dihedral_start_index):
    """
    Extract and reindex dihedral information based on filtered atom types.
    
    This function processes the Dihedrals section of a molecule file, keeping only dihedrals
    where all four atoms exist in the filtered type_df DataFrame. Dihedral indices and atom
    indices are reassigned to match the new atom indices from the types modification.
    
    Args:
        lines (list): All lines from the molecule file.
        type_df (pd.DataFrame): DataFrame containing filtered atom types with new indices.
        dihedral_start_index (int): Starting line index of the Dihedrals section.
    
    Returns:
        tuple: A tuple containing:
            - dihedral_section (str): Formatted string representation of the dihedrals section
            - number_of_dihedrals (int): Total number of dihedrals in the filtered data
    
    Note:
        - Only dihedrals where all four atoms are present in type_df are retained
        - Dihedral indices are reassigned sequentially starting from 1
        - Atom indices are updated to match the new indices from type_df
    """
    # Initialize list to store parsed dihedral data
    angle_data = []
    
    # Parse the Dihedrals section from the input lines
    for line in lines[dihedral_start_index:]:
        # Stop processing when an empty line is encountered
        if not line.strip():
            break
        
        # Split the line into individual components
        parts = line.split()
        
        # Validate that the line has sufficient data fields
        if len(parts) >= 11:
            # Extract dihedral properties from the parsed line
            dihedral_index = int(parts[0])
            dihedral_type = int(parts[1])
            atom1_index = int(parts[2])
            atom2_index = int(parts[3])
            atom3_index = int(parts[4])
            atom4_index = int(parts[5])
            hash_value = parts[6]
            atom1_type_real = parts[7]
            atom2_type_real = parts[8]
            atom3_type_real = parts[9]
            atom4_type_real = parts[10]
            
            # Append parsed dihedral data to the list
            angle_data.append({
                "dihedral_index": dihedral_index,
                "dihedral_type": dihedral_type,
                "atom1_index": atom1_index,
                "atom2_index": atom2_index,
                "atom3_index": atom3_index,
                "atom4_index": atom4_index,
                "hash": hash_value,
                "atom1_type_real": atom1_type_real,
                "atom2_type_real": atom2_type_real,
                "atom3_type_real": atom3_type_real,
                "atom4_type_real": atom4_type_real
            })
    
    # Create a DataFrame from the parsed dihedral data
    dihedral_df = pd.DataFrame(angle_data)
    
    # Initialize new columns for reindexed atom indices
    dihedral_df["new_atom1_index"] = None 
    dihedral_df["new_atom2_index"] = None 
    dihedral_df["new_atom3_index"] = None 
    dihedral_df["new_atom4_index"] = None
    
    # Map original atom1 indices to new indices using the type_df mapping
    for value in dihedral_df["atom1_index"]:
        if value in type_df["atom_index"].values:
            # Look up the new index from the type_df
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            dihedral_df.loc[dihedral_df["atom1_index"] == value, "new_atom1_index"] = new_index
    
    # Map original atom2 indices to new indices using the type_df mapping
    for value in dihedral_df["atom2_index"]:
        if value in type_df["atom_index"].values:
            # Look up the new index from the type_df
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            dihedral_df.loc[dihedral_df["atom2_index"] == value, "new_atom2_index"] = new_index
    
    # Map original atom3 indices to new indices using the type_df mapping
    for value in dihedral_df["atom3_index"]:
        if value in type_df["atom_index"].values:
            # Look up the new index from the type_df
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            dihedral_df.loc[dihedral_df["atom3_index"] == value, "new_atom3_index"] = new_index
    
    # Map original atom4 indices to new indices using the type_df mapping
    for value in dihedral_df["atom4_index"]:
        if value in type_df["atom_index"].values:
            # Look up the new index from the type_df
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            dihedral_df.loc[dihedral_df["atom4_index"] == value, "new_atom4_index"] = new_index
    
    # Identify rows where any atom is not in the filtered type data (to be removed)
    indices_to_drop = []        
    for index, row in dihedral_df.iterrows():
        if (row["new_atom1_index"] is None or 
            row["new_atom2_index"] is None or 
            row["new_atom3_index"] is None or
            row["new_atom4_index"] is None):
            indices_to_drop.append(index)
    
    # Remove rows for dihedrals with missing atoms and reset the index
    dihedral_df = dihedral_df.drop(indices_to_drop).reset_index(drop=True)        
    
    # Reassign dihedral indices sequentially starting from 1
    for i, value in enumerate(dihedral_df.index):
        dihedral_df.at[i, "dihedral_index"] = i + 1
    
    # Format the dihedrals section as a string with proper spacing
    dihedral_section = ""
    for i, row in dihedral_df.iterrows():
        dihedral_section += f"{row['dihedral_index']:>4}{row['dihedral_type']:>4}{row['new_atom1_index']:>4}{row['new_atom2_index']:>4}{row['new_atom3_index']:>4}{row['new_atom4_index']:>4}{row['hash']:>4}{row['atom1_type_real']:>4}{row['atom2_type_real']:>4}{row['atom3_type_real']:>4}{row['atom4_type_real']:>4}\n"
    
    # Calculate the total number of dihedrals
    number_of_dihedrals = len(dihedral_df)
    
    return dihedral_section, number_of_dihedrals

def modify_impropers(lines, type_df, improper_start_index):
    """
    Parse and modify improper dihedral data from molecule file lines.
    
    This function extracts improper dihedral information from the input lines,
    maps old atom indices to new atom indices using the type dataframe, and
    returns a formatted improper section string along with the count.
    
    Parameters:
    -----------
    lines : list[str]
        All lines from the molecule file to be processed.
    type_df : pd.DataFrame
        DataFrame containing atom index mapping with columns:
        - 'atom_index': Original atom indices
        - 'new_atom_index': Mapped atom indices for the product molecule
    improper_start_index : int
        Line number where the Impropers section begins in the file.
    
    Returns:
    --------
    tuple[str, int]
        - improper_section (str): Formatted string of improper data ready for output
        - number_of_impropers (int): Total count of valid impropers after filtering
    
    Raises:
    -------
    IndexError: If a line has fewer than 11 parts (malformed data)
    
    Notes:
    ------
    - Impropers with any missing atom index mappings are dropped
    - Improper indices are re-indexed starting from 1
    """
    # Initialize list to store parsed improper dihedral data
    angle_data = []
    
    # Iterate through lines starting from improper section
    for line in lines[improper_start_index:]:
        # Break if empty line encountered (marks end of section)
        if not line.strip():
            break
        
        # Split line into parts
        parts = line.split()
        
        # Validate that line has minimum required fields (11 parts)
        if len(parts) >= 11:
            # Parse improper data from line parts
            improper_index = int(parts[0])
            improper_type = int(parts[1])
            atom1_index = int(parts[2])
            atom2_index = int(parts[3])
            atom3_index = int(parts[4])
            atom4_index = int(parts[5])
            hash_value = parts[6]
            atom1_type_real = parts[7]
            atom2_type_real = parts[8]
            atom3_type_real = parts[9]
            atom4_type_real = parts[10]
            
            # Store parsed data as dictionary
            angle_data.append({
                "improper_index": improper_index,
                "improper_type": improper_type,
                "atom1_index": atom1_index,
                "atom2_index": atom2_index,
                "atom3_index": atom3_index,
                "atom4_index": atom4_index,
                "hash": hash_value,
                "atom1_type_real": atom1_type_real,
                "atom2_type_real": atom2_type_real,
                "atom3_type_real": atom3_type_real,
                "atom4_type_real": atom4_type_real
            })
    
    # Convert list of dictionaries to DataFrame
    improper_df = pd.DataFrame(angle_data)
    
    # Initialize columns for new atom indices (will be populated with mapped values)
    improper_df["new_atom1_index"] = None
    improper_df["new_atom2_index"] = None
    improper_df["new_atom3_index"] = None
    improper_df["new_atom4_index"] = None
    
    # Map atom1 indices from template to product molecule
    for value in improper_df["atom1_index"]:
        if value in type_df["atom_index"].values:
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            improper_df.loc[improper_df["atom1_index"] == value, "new_atom1_index"] = new_index
    
    # Map atom2 indices from template to product molecule
    for value in improper_df["atom2_index"]:
        if value in type_df["atom_index"].values:
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            improper_df.loc[improper_df["atom2_index"] == value, "new_atom2_index"] = new_index
    
    # Map atom3 indices from template to product molecule
    for value in improper_df["atom3_index"]:
        if value in type_df["atom_index"].values:
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            improper_df.loc[improper_df["atom3_index"] == value, "new_atom3_index"] = new_index
    
    # Map atom4 indices from template to product molecule
    for value in improper_df["atom4_index"]:
        if value in type_df["atom_index"].values:
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            improper_df.loc[improper_df["atom4_index"] == value, "new_atom4_index"] = new_index
    
    # Identify impropers with incomplete mappings (any None value indicates missing atom)
    indices_to_drop = []
    for index, row in improper_df.iterrows():
        if (row["new_atom1_index"] is None or 
            row["new_atom2_index"] is None or 
            row["new_atom3_index"] is None or
            row["new_atom4_index"] is None):
            indices_to_drop.append(index)
    
    # Remove impropers with incomplete mappings and reset index
    improper_df = improper_df.drop(indices_to_drop).reset_index(drop=True)
    
    # Re-index impropers sequentially starting from 1
    for i, value in enumerate(improper_df.index):
        improper_df.at[i, "improper_index"] = i + 1
    
    # Build formatted improper section string with proper spacing
    improper_section = ""
    for i, row in improper_df.iterrows():
        improper_section += f"{row['improper_index']:>4}{row['improper_type']:>4}{row['new_atom1_index']:>4}{row['new_atom2_index']:>4}{row['new_atom3_index']:>4}{row['new_atom4_index']:>4}{row['hash']:>4}{row['atom1_type_real']:>4}{row['atom2_type_real']:>4}{row['atom3_type_real']:>4}{row['atom4_type_real']:>4}\n"
    
    # Get final count of impropers
    number_of_impropers = len(improper_df)
    
    return improper_section, number_of_impropers


def molecule_file_format(number_of_types, number_of_bonds, number_of_angles, 
                         number_of_dihedrals, number_of_impropers, types_section, charge_section, coord_section, 
                         bond_section, angle_section, dihedral_section, improper_section):
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


def molecule_file_preparation(test_molecule_file, template_indexes):
    """
    Orchestrate the complete molecule file processing pipeline.
    
    This function coordinates all steps needed to convert a template molecule
    file to a product molecule file by parsing sections, modifying atom indices,
    and reformatting the data.
    
    Parameters:
    -----------
    test_molecule_file : str
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
    lines, type_start_index, charge_start_index, coord_start_index, bond_start_index, angle_start_index, dihedral_start_index, improper_start_index = load_molecule_file(test_molecule_file)
    
    # Process types section - returns DataFrame for mapping and formatted string
    df_types, types_section, number_of_types, index_change_dict = modify_types(lines, template_indexes, type_start_index)
    
    # Process charges section using atom mapping from types
    charge_section = modify_charges(lines, df_types, charge_start_index)
    
    # Process coordinates section using atom mapping from types
    coord_section = modify_coords(lines, df_types, coord_start_index)
    
    # Process bonds section and get count
    bond_section, number_of_bonds = modify_bonds(lines, df_types, bond_start_index)
    
    # Process angles section and get count
    angle_section, number_of_angles = modify_angles(lines, df_types, angle_start_index)
    
    # Process dihedrals section and get count
    dihedral_section, number_of_dihedrals = modify_dihedrals(lines, df_types, dihedral_start_index)
    
    # Process impropers section and get count
    improper_section, number_of_impropers = modify_impropers(lines, df_types, improper_start_index)
    
    # Assemble all sections into final molecule file format
    modified_molecule_file = molecule_file_format(number_of_types, number_of_bonds, number_of_angles, 
                                                 number_of_dihedrals, number_of_impropers, types_section, charge_section, coord_section, 
                                                 bond_section, angle_section, dihedral_section, improper_section)
    
    return modified_molecule_file, index_change_dict

def map_file_write(reactant_to_product, initiator_atoms, edge_atoms, delete_ids):
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


def build_bond_react_templates(
    file_dict,
    reactant_to_product,
    initiator_atoms,
    edge_atoms,
    delete_ids,
    cache_dir_reactor
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
      Directory path where reactor-related files should be cached. The function will write output files
      to the directory that contains this path (i.e., os.path.dirname(cache_dir_reactor)).

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
    molecule_file_dict = {}

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
        num = get_ending_integer(name)
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
        
        pre_modified, index_change_dict_reactant = molecule_file_preparation(
            pre_path,
            template_indexes_reactant
        )

        post_modified, index_change_dict_product = molecule_file_preparation(
            post_path,
            template_indexes_product
        )

        # Write the modified molecule files to the cache directory's parent (consistent with previous code)
        os.makedirs(cache_dir_reactor, exist_ok=True)

        pre_out = os.path.join(cache_dir_reactor, f"template_pre_{num}.molecule")
        with open(pre_out, "w") as f:
            f.write(pre_modified)

        post_out = os.path.join(cache_dir_reactor, f"template_post_{num}.molecule")
        with open(post_out, "w") as f:
            f.write(post_modified)

        # Record the generated file paths in the return dictionary
        molecule_file_dict[f"pre_{num}"] = pre_out
        molecule_file_dict[f"post_{num}"] = post_out

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
        map_file = map_file_write(
            template_reactant_to_template_product,
            initiator_atoms_t,
            edge_atoms_t,
            delete_ids_t
        )

        map_path = os.path.join(cache_dir_reactor, f"RXN_{num}.map")
        with open(map_path, "w") as f:
            f.write(map_file)

        molecule_file_dict[f"map_file_{num}"] = map_path

    return molecule_file_dict


def molecule_template_preparation(molecule_dict_csv_path_dict, lunar_out_loc_dict, cache_path):
    """
    Top-level orchestrator that loops over reactions and prepares template
    files and mappings for each reaction.

    Parameters
    ----------
    molecule_dict_csv_path_dict : dict
      Mapping reaction id -> dict containing at least:
        - "reaction_dataframe": a pandas DataFrame describing the reaction and template indices
        Optionally it may include "delete_atoms" boolean to indicate byproducts should be removed.
    lunar_out_loc_dict : dict
      Mapping names like "pre_{rxn_id}" and "post_{rxn_id}" to file paths for molecule files.
    cache_path : str
      Base cache directory where reactor files will be stored. This function will create
      a subdirectory "reactor_files" under cache_path.

    Returns
    -------
    total_molecule_file_dict : dict
      Combined dictionary of generated file paths for all processed reactions.
    """
    reactor_dir = ensure_dir(os.path.join(cache_path, "reactor_files"))
    
    total_molecule_file_dict = {}

    # Iterate each reaction entry and build templates for that single reaction only.
    for rxn_id, rxn in molecule_dict_csv_path_dict.items():
        print(f"Preparing templates and map file for reaction ID: {rxn_id}")
        # The reaction DataFrame is expected under the key "reaction_dataframe"
        df: pd.DataFrame = rxn["reaction_dataframe"]
        # Optional indicator whether byproducts listed should be treated as delete IDs
        delete_atoms: bool = bool(rxn.get("delete_atoms", False))

        # Extract integer lists from DataFrame columns using helper _col_int_list
        byproducts = _col_int_list("byproduct_indices", df = df)
        initiators = _col_int_list("initiators", df = df)
        edge_atoms = _col_int_list("edge_atoms", df = df)

        # Build a template-level mapping (if present in DataFrame). Keys/values are expected
        # to be 0-based indices in the original CSV context.
        template_map = {}
        if {"template_reactant_idx", "template_product_idx"}.issubset(df.columns):
            tmp = df[["template_reactant_idx", "template_product_idx"]].dropna()
            for r_idx, p_idx in tmp.itertuples(index=False):
                template_map[int(r_idx)] = int(p_idx)

        # If delete_atoms flag is set, delete_ids come from byproducts; otherwise none.
        delete_ids = byproducts if delete_atoms else []

        # Target only the specific files for this rxn_id so we avoid redundant looping.
        current_rxn_files = {
            f"pre_{rxn_id}": lunar_out_loc_dict.get(f"pre_{rxn_id}"),
            f"post_{rxn_id}": lunar_out_loc_dict.get(f"post_{rxn_id}")
        }

        # Build the per-reaction template files and mappings.
        molecule_file_dict = build_bond_react_templates(
            file_dict = current_rxn_files,
            reactant_to_product = template_map,
            initiator_atoms = initiators,
            edge_atoms = edge_atoms,
            delete_ids = delete_ids,
            cache_dir_reactor = reactor_dir
        )
        total_molecule_file_dict.update(molecule_file_dict)

    return total_molecule_file_dict


if __name__ == "__main__":
    """
    Main execution block for testing the template preparation functions.
    
    This demonstrates how to use the functions with sample file paths
    and generates output for multiple reactions.
    """
    # Sample file dictionary with molecule file paths
    molecule_file_dict = {
        "data_1": r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\00_cache\lunar\bond_react_merge\data_1_typed_IFF_merged.data",
        "data_2": r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\00_cache\lunar\bond_react_merge\data_2_typed_IFF_merged.data",
        "data_3": r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\00_cache\lunar\bond_react_merge\data_3_typed_IFF_merged.data",

        "pre_1":  r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\00_cache\lunar\bond_react_merge\pre_1_typed_IFF_merged.lmpmol",
        "post_1": r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\00_cache\lunar\bond_react_merge\post_1_typed_IFF_merged.lmpmol",

        "pre_2":  r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\00_cache\lunar\bond_react_merge\pre_2_typed_IFF_merged.lmpmol",
        "post_2": r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\00_cache\lunar\bond_react_merge\post_2_typed_IFF_merged.lmpmol",

        "pre_3":  r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\00_cache\lunar\bond_react_merge\pre_3_typed_IFF_merged.lmpmol",
        "post_3": r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\00_cache\lunar\bond_react_merge\post_3_typed_IFF_merged.lmpmol",

        "pre_4":  r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\00_cache\lunar\bond_react_merge\pre_4_typed_IFF_merged.lmpmol",
        "post_4": r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\00_cache\lunar\bond_react_merge\post_4_typed_IFF_merged.lmpmol",

        "force_field_data": r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\00_cache\lunar\bond_react_merge\force_field.data",
    }

    # CSV file paths containing reaction mapping data
    csv_dir = r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\00_cache\csv\reactant_product_mapping"
    csv_paths = {
        1: os.path.join(csv_dir, "reaction_1.csv"),
        2: os.path.join(csv_dir, "reaction_2.csv"),
        3: os.path.join(csv_dir, "reaction_3.csv"),
        4: os.path.join(csv_dir, "reaction_4.csv"),
    }

    all_outputs = {}

    # Process each reaction CSV file
    for rxn_id, csv_path in csv_paths.items():
        df = pd.read_csv(csv_path)

        # Extract atom lists from CSV
        byproducts = _col_int_list("byproduct_indices", df=df)
        initiators = _col_int_list("initiators", df=df)
        edge_atoms = _col_int_list("edge_atoms", df=df)

        # Build template mapping
        template_map = {}
        if {"template_reactant_idx", "template_product_idx"}.issubset(df.columns):
            tmp = df[["template_reactant_idx", "template_product_idx"]].dropna()
            for r_idx, p_idx in tmp.itertuples(index=False):
                template_map[int(r_idx)] = int(p_idx)

        # Determine if atoms should be deleted
        delete_atoms = False
        if "delete_atoms" in df.columns:
            delete_atoms = bool(df["delete_atoms"].dropna().iloc[0]) if not df["delete_atoms"].dropna().empty else False

        delete_ids = byproducts if delete_atoms else []

        # Generate bond/react templates
        out = build_bond_react_templates(
            file_dict=molecule_file_dict,
            reactant_to_product=template_map,
            initiator_atoms=initiators,
            edge_atoms=edge_atoms,
            delete_ids=delete_ids,
            cache_dir_reactor=os.path.join(os.path.dirname(csv_dir), "reactor_files"),
        )

        all_outputs[rxn_id] = out

    # Save results to file
    save_path = r"C:\Users\Janitha\Documents\GitHub\reaction_lammps_mupt\cache\00_cache\lunar"
    os.makedirs(save_path, exist_ok=True)

    with open(os.path.join(save_path, "molecule_file_dict.txt"), "w") as f:
        f.write(json.dumps(all_outputs, indent=4))

    # Print results for verification
    print(json.dumps(all_outputs, indent=4))
