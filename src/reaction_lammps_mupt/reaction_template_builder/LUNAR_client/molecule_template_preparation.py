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

import pandas as pd
import re
import os

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
    
    Note:
        The function assumes that data for each section starts 2 lines after the section header.
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
        - Output is saved to a hardcoded file path (should be parameterized)
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
    
    return df, types_section, number_of_types



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
    df_types, types_section, number_of_types = modify_types(lines, template_indexes, type_start_index)
    
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
    
    return modified_molecule_file


def load_molecule_files(file_dict, reactant_to_product):
    """
    Process multiple molecule files for bond-react simulation setup.
    
    Iterates through pre- and post-reaction molecule files, processes each one
    using the appropriate atom index mapping, and saves the results as template
    molecule files.
    
    Parameters:
    -----------
    file_dict : dict[str, str]
        Dictionary mapping file identifiers to file paths:
        - Keys: 'pre1', 'post1', etc.
        - Values: Full file paths to molecule files
    reactant_to_product : dict[int, int]
        Mapping of reactant atom indices to product atom indices.
        - Keys: Reactant atom indices (0-based in dictionary, converted to 1-based)
        - Values: Product atom indices (0-based in dictionary, converted to 1-based)
    
    Returns:
    --------
    None
        Writes template molecule files to disk in the same directory as input files.
    
    Output Files:
    --------
    - template_pre{n}.molecule: Processed pre-reaction molecule file
    - template_post{n}.molecule: Processed post-reaction molecule file
    
    Raises:
    -------
    ValueError: If expected molecule data file is not found in file_dict.
    
    TODO:
    -----
    Currently processes files with pre/post naming. Future enhancement could:
    - Support arbitrary file patterns
    - Process multiple reaction pathways
    - Validate file consistency before processing
    """
    # Initialize lists to store atom index mappings for reactants and products
    template_indexes_reactant = []
    template_indexes_product = []
    
    # Build 1-based index lists from reactant_to_product mapping
    # LAMMPS uses 1-based indexing while dictionary uses 0-based
    for key, value in reactant_to_product.items():
        # Convert reactant indices to 1-based indexing
        template_indexes_reactant.append(int(key + 1))
        # Convert product indices to 1-based indexing
        template_indexes_product.append(int(value + 1))
    
    # Iterate through all files in the file dictionary
    for name, filepath in file_dict.items():
        # Determine if file is pre-reaction or post-reaction and extract number
        if name.startswith('pre'):
            # Extract file number from name (e.g., 'pre1' -> 1)
            num = get_ending_integer(name)
            # Retrieve molecule file path from dictionary
            molecule_file = file_dict.get(f'pre{num}')
            # Store file identifier for output naming
            file_name = f'pre{num}'
        elif name.startswith('post'):
            # Extract file number from name (e.g., 'post1' -> 1)
            num = get_ending_integer(name)
            # Retrieve molecule file path from dictionary
            molecule_file = file_dict[f'post{num}']
            # Store file identifier for output naming
            file_name = f'post{num}'
        else:
            # Skip files that don't match pre/post pattern
            continue
        
        # Validate that molecule file path was found
        if molecule_file is None:
            raise ValueError(f"Data file for {name} not found.")
        
        # Process molecule file with appropriate index mapping
        # Use reactant indices for 'pre' files, product indices for 'post' files
        modified_molecule = molecule_file_preparation(
            filepath,
            template_indexes_reactant if name.startswith('pre') else template_indexes_product
        )
        
        # Get directory path from input file path
        save_path = os.path.dirname(filepath)
        
        # Write processed molecule file to disk
        with open(os.path.join(save_path, f"template_{file_name}.molecule"), "w") as f:
            f.write(modified_molecule)
            file_path = os.path.join(save_path, f"template_{file_name}.molecule")
        
        # Update the dictionary with the modified molecule template file path
        molecule_file_dict[name] = file_path
    return molecule_file_dict


if __name__ == "__main__":
    """
    Main execution block for bond-react molecule file processing.
    
    Sets up file paths and atom index mappings, then processes molecule files
    for a bond-react LAMMPS simulation.
    """
    
    # Dictionary mapping file identifiers to their full file paths
    molecule_file_dict = {
        # Data files for system definition
        "data1": "C:\\Users\\Janitha\\Documents\\GitHub\\reaction_lammps_mupt\\cache\\lunar\\bond_react_merge\\data1_typed_IFF_merged.data",
        "data2": "C:\\Users\\Janitha\\Documents\\GitHub\\reaction_lammps_mupt\\cache\\lunar\\bond_react_merge\\data2_typed_IFF_merged.data",
        # Pre-reaction template molecule file
        "pre1": "C:\\Users\\Janitha\\Documents\\GitHub\\reaction_lammps_mupt\\cache\\lunar\\bond_react_merge\\pre1_typed_IFF_merged.lmpmol",
        # Post-reaction template molecule file
        "post1": "C:\\Users\\Janitha\\Documents\\GitHub\\reaction_lammps_mupt\\cache\\lunar\\bond_react_merge\\post1_typed_IFF_merged.lmpmol",
        # Force field parameters
        "force_field_data": "C:\\Users\\Janitha\\Documents\\GitHub\\reaction_lammps_mupt\\cache\\lunar\\bond_react_merge\\force_field.data"
    }

    # Mapping of reactant atom indices to product atom indices
    # Used for tracking atoms through the chemical reaction
    # Keys/Values are 0-based indices (converted to 1-based in processing)
    reactant_to_product = {
        1: 6,      
        2: 4,      
        3: 3,      
        4: 5,      
        5: 8,      
        6: 9,      
        9: 0,      
        12: 7,     
        15: 26,    
        16: 21,    
        17: 18,    
        18: 17,    
        19: 1,     
        20: 25,    
        21: 2, 
        23: 22,    
        24: 23,    
        25: 19,    
        26: 20,    
        27: 27,    
    }
    
    # Initialize lists for storing converted 1-based indices
    template_indexes_reactant = []
    template_indexes_product = []
    
    # Process all molecule files and generate templates
    molecule_file_dict = load_molecule_files(molecule_file_dict, reactant_to_product)    
    print("Molecule files processed and templates generated.")
    import json
    print(json.dumps(molecule_file_dict, indent=4))