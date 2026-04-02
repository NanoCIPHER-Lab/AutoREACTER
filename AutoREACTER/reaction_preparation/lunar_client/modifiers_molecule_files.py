import pandas as pd

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
    legacy_mode = True
    
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
            try:
                atom_index = int(parts[0])
                atom_type = int(parts[1])
            except (ValueError, IndexError): # to faciltate the new lammps atom typings with string types
                atom_index = int(parts[0])
                atom_type = str(parts[1])
                legacy_mode = False
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
    
    return df, types_section, number_of_types, index_change_dict, legacy_mode



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
        if len(parts) >= 4:
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
        if len(parts) >= 6:
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


def modify_bonds(lines, type_df, bond_start_index, legacy_mode):
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
            if legacy_mode:
                bond_type = int(parts[1])
            else:
                bond_type = parts[1]
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

def modify_angles(lines, type_df, angle_start_index, legacy_mode):
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
            angle_type = int(parts[1]) if legacy_mode else parts[1]
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

def modify_dihedrals(lines, type_df, dihedral_start_index, legacy_mode):
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
            dihedral_type = int(parts[1]) if legacy_mode else parts[1]
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

def modify_impropers(lines, type_df, improper_start_index, legacy_mode):
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
            if legacy_mode:
                improper_type = int(parts[1])
            else:
                improper_type = parts[1]
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