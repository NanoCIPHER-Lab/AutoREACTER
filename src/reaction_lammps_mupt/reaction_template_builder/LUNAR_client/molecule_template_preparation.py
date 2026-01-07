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
    with open(molecule_file, 'r') as file:
        if file is None:
            raise ValueError(f"File {molecule_file} not found.")
        lines = [line.strip() for line in file]
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
    if type_start_index is None or charge_start_index is None or coord_start_index is None:
        raise ValueError("Essential sections (Types, Charges, Coords) not found in the file.")
    return lines, type_start_index, charge_start_index, coord_start_index, bond_start_index, angle_start_index, dihedral_start_index, improper_start_index

def modify_types(lines, template_indexes, type_start_index):
    template_indexes.sort()
    data = []
    for line in lines[type_start_index:]:
        if not line.strip(): # More Pythonic way to check for empty strings
            break
        
        parts = line.split()
        
        # Check if the line has enough parts to avoid an IndexError
        if len(parts) >= 4:
            atom_index = int(parts[0])
            atom_type = int(parts[1])
            hash_value = parts[2]
            atom_type_real = parts[3]
            
            # Using the already converted variables directly
            data.append({
                "atom_index": atom_index,
                "atom_type": atom_type,
                "hash": hash_value,
                "atom_type_real": atom_type_real
            })

    df = pd.DataFrame(data)
    df["new_atom_index"] = None 
    for value in df["atom_index"]:
        if value in template_indexes:
            new_type = value
            df.loc[df["atom_index"] == value, "new_atom_index"] = new_type
    indices_to_drop = []        
    for index, row in df.iterrows():
        if row["new_atom_index"] is None:
            indices_to_drop.append(index)
    df = df.drop(indices_to_drop).reset_index(drop=True)        
    for i,value in enumerate(df["new_atom_index"]):
        df.at[i, "new_atom_index"] = i + 1
    
    with open("C:\\Users\\Janitha\\Documents\\GitHub\\reaction_lammps_mupt\\cache\\lunar\\bond_react_merge\\modified_types.txt", "w") as f:
        f.write(df.to_string(index=False))
    types_section = ""
    for i, row in df.iterrows():
        types_section += f"{row['new_atom_index']:>4}{row['atom_type']:>4}{row['hash']:>4}{row['atom_type_real']:>4}\n"
    number_of_types = len(df)
    return df, types_section, number_of_types

def modify_charges(lines, type_df, charge_start_index):
    charge_data = []
    for line in lines[charge_start_index:]:
        if not line.strip(): # More Pythonic way to check for empty strings
            break

        parts = line.split()
        
        # Check if the line has enough parts to avoid an IndexError
        if len(parts) >= 3:
            atom_index = int(parts[0])
            charge_value = float(parts[1])
            hash_value = parts[2]
            atom_type_real = parts[3]
            
            # Using the already converted variables directly
            charge_data.append({
                "atom_index": atom_index,
                "charge_value": charge_value,
                "hash": hash_value,
                "atom_type_real": atom_type_real
            })

    charge_df = pd.DataFrame(charge_data)
    charge_df["new_atom_index"] = None 
    for value in charge_df["atom_index"]:
        if value in type_df["atom_index"].values:
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            charge_df.loc[charge_df["atom_index"] == value, "new_atom_index"] = new_index
    indices_to_drop = []        
    for index, row in charge_df.iterrows():
        if row["new_atom_index"] is None:
            indices_to_drop.append(index)
    charge_df = charge_df.drop(indices_to_drop).reset_index(drop=True)        
    for i,value in enumerate(charge_df["new_atom_index"]):
        charge_df.at[i, "new_atom_index"] = i + 1
    with open("C:\\Users\\Janitha\\Documents\\GitHub\\reaction_lammps_mupt\\cache\\lunar\\bond_react_merge\\modified_charges.txt", "w") as f:
        f.write(charge_df.to_string(index=False))
    charge_section = ""
    for i, row in charge_df.iterrows():
        charge_section += f"{row['new_atom_index']:>4}{row['charge_value']:>12.6f}{row['hash']:>4}{row['atom_type_real']:>4}\n"
    return  charge_section

def modify_coords(lines, type_df, coord_start_index):
    coord_data = []
    for line in lines[coord_start_index:]:
        if not line.strip(): # More Pythonic way to check for empty strings
            break

        parts = line.split()
        
        # Check if the line has enough parts to avoid an IndexError
        if len(parts) >= 5:
            atom_index = int(parts[0])
            x = float(parts[1])
            y = float(parts[2])
            z = float(parts[3])
            hash_value = parts[4]
            atom_type_real = parts[5]
            
            # Using the already converted variables directly
            coord_data.append({
                "atom_index": atom_index,
                "x": x,
                "y": y,
                "z": z,
                "hash": hash_value,
                "atom_type_real": atom_type_real
            })

    coord_df = pd.DataFrame(coord_data)
    coord_df["new_atom_index"] = None 
    for value in coord_df["atom_index"]:
        if value in type_df["atom_index"].values:
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            coord_df.loc[coord_df["atom_index"] == value, "new_atom_index"] = new_index
    indices_to_drop = []        
    for index, row in coord_df.iterrows():
        if row["new_atom_index"] is None:
            indices_to_drop.append(index)
    coord_df = coord_df.drop(indices_to_drop).reset_index(drop=True)        
    for i,value in enumerate(coord_df["new_atom_index"]):
        coord_df.at[i, "new_atom_index"] = i + 1
    with open("C:\\Users\\Janitha\\Documents\\GitHub\\reaction_lammps_mupt\\cache\\lunar\\bond_react_merge\\modified_coords.txt", "w") as f:
        f.write(coord_df.to_string(index=False))
    coord_section = ""
    for i, row in coord_df.iterrows():
        coord_section += f"{row['new_atom_index']:>4}{row['x']:>12.6f}{row['y']:>12.6f}{row['z']:>12.6f}{row['hash']:>4}{row['atom_type_real']:>4}\n"
    return coord_section

def modify_bonds(lines, type_df, bond_start_index):
    bond_data = []
    for line in lines[bond_start_index:]:
        if not line.strip(): # More Pythonic way to check for empty strings
            break

        parts = line.split()
        
        # Check if the line has enough parts to avoid an IndexError
        if len(parts) >= 7:
            bond_index = int(parts[0])
            bond_type = int(parts[1])
            atom1_index = int(parts[2])
            atom2_index = int(parts[3])
            hash_value = parts[4]
            atom1_type_real = parts[5]
            atom2_type_real = parts[6]
            
            # Using the already converted variables directly
            bond_data.append({
                "bond_index": bond_index,
                "bond_type": bond_type,
                "atom1_index": atom1_index,
                "atom2_index": atom2_index,
                "hash": hash_value,
                "atom1_type_real": atom1_type_real,
                "atom2_type_real": atom2_type_real
            })

    bond_df = pd.DataFrame(bond_data)
    bond_df["new_atom1_index"] = None 
    bond_df["new_atom2_index"] = None 
    for value in bond_df["atom1_index"]:
        if value in type_df["atom_index"].values:
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            bond_df.loc[bond_df["atom1_index"] == value, "new_atom1_index"] = new_index
    for value in bond_df["atom2_index"]:
        if value in type_df["atom_index"].values:
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            bond_df.loc[bond_df["atom2_index"] == value, "new_atom2_index"] = new_index
    indices_to_drop = []        
    for index, row in bond_df.iterrows():
        if row["new_atom1_index"] is None or row["new_atom2_index"] is None:
            indices_to_drop.append(index)
    bond_df = bond_df.drop(indices_to_drop).reset_index(drop=True)        
    for i,value in enumerate(bond_df.index):
        bond_df.at[i, "bond_index"] = i + 1
    with open("C:\\Users\\Janitha\\Documents\\GitHub\\reaction_lammps_mupt\\cache\\lunar\\bond_react_merge\\modified_bonds.txt", "w") as f:
        f.write(bond_df.to_string(index=False))
    bond_section = ""
    for i, row in bond_df.iterrows():
        bond_section += f"{row['bond_index']:>4}{row['bond_type']:>4}{row['new_atom1_index']:>4}{row['new_atom2_index']:>4}{row['hash']:>4}{row['atom1_type_real']:>4}{row['atom2_type_real']:>4}\n"
    number_of_bonds = len(bond_df)
    return bond_section, number_of_bonds

def modify_angles(lines, type_df, angle_start_index):
    angle_data = []
    for line in lines[angle_start_index:]:
        if not line.strip(): # More Pythonic way to check for empty strings
            break
        parts = line.split()
        # Check if the line has enough parts to avoid an IndexError
        if len(parts) >= 9:
            angle_index = int(parts[0])
            angle_type = int(parts[1])
            atom1_index = int(parts[2])
            atom2_index = int(parts[3])
            atom3_index = int(parts[4])
            hash_value = parts[5]
            atom1_type_real = parts[6]
            atom2_type_real = parts[7]
            atom3_type_real = parts[8]
            # Using the already converted variables directly
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
    angle_df = pd.DataFrame(angle_data)
    angle_df["new_atom1_index"] = None 
    angle_df["new_atom2_index"] = None 
    angle_df["new_atom3_index"] = None 
    for value in angle_df["atom1_index"]:
        if value in type_df["atom_index"].values:
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            angle_df.loc[angle_df["atom1_index"] == value, "new_atom1_index"] = new_index
    for value in angle_df["atom2_index"]:
        if value in type_df["atom_index"].values:
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            angle_df.loc[angle_df["atom2_index"] == value, "new_atom2_index"] = new_index
    for value in angle_df["atom3_index"]:
        if value in type_df["atom_index"].values:
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            angle_df.loc[angle_df["atom3_index"] == value, "new_atom3_index"] = new_index
    indices_to_drop = []        
    for index, row in angle_df.iterrows():
        if (row["new_atom1_index"] is None or 
            row["new_atom2_index"] is None or 
            row["new_atom3_index"] is None):
            indices_to_drop.append(index)
    angle_df = angle_df.drop(indices_to_drop).reset_index(drop=True)        
    for i,value in enumerate(angle_df.index):
        angle_df.at[i, "angle_index"] = i + 1
    with open("C:\\Users\\Janitha\\Documents\\GitHub\\reaction_lammps_mupt\\cache\\lunar\\bond_react_merge\\modified_angles.txt", "w") as f:
        f.write(angle_df.to_string(index=False))
    angle_section = ""
    for i, row in angle_df.iterrows():
        angle_section += f"{row['angle_index']:>4}{row['angle_type']:>4}{row['new_atom1_index']:>4}{row['new_atom2_index']:>4}{row['new_atom3_index']:>4}{row['hash']:>4}{row['atom1_type_real']:>4}{row['atom2_type_real']:>4}{row['atom3_type_real']:>4}\n"
    number_of_angles = len(angle_df)
    return angle_section, number_of_angles

def modify_dihedrals(lines, type_df, dihedral_start_index):
    angle_data = []
    for line in lines[dihedral_start_index:]:
        if not line.strip(): # More Pythonic way to check for empty strings
            break
        parts = line.split()
        # Check if the line has enough parts to avoid an IndexError
        if len(parts) >= 11:
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
            # Using the already converted variables directly
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
    dihedral_df = pd.DataFrame(angle_data)
    dihedral_df["new_atom1_index"] = None 
    dihedral_df["new_atom2_index"] = None 
    dihedral_df["new_atom3_index"] = None 
    dihedral_df["new_atom4_index"] = None
    for value in dihedral_df["atom1_index"]:
        if value in type_df["atom_index"].values:
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            dihedral_df.loc[dihedral_df["atom1_index"] == value, "new_atom1_index"] = new_index
    for value in dihedral_df["atom2_index"]:
        if value in type_df["atom_index"].values:
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            dihedral_df.loc[dihedral_df["atom2_index"] == value, "new_atom2_index"] = new_index
    for value in dihedral_df["atom3_index"]:
        if value in type_df["atom_index"].values:
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            dihedral_df.loc[dihedral_df["atom3_index"] == value, "new_atom3_index"] = new_index
    for value in dihedral_df["atom4_index"]:
        if value in type_df["atom_index"].values:
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            dihedral_df.loc[dihedral_df["atom4_index"] == value, "new_atom4_index"] = new_index
    indices_to_drop = []        
    for index, row in dihedral_df.iterrows():
        if (row["new_atom1_index"] is None or 
            row["new_atom2_index"] is None or 
            row["new_atom3_index"] is None or
            row["new_atom4_index"] is None):
            indices_to_drop.append(index)
    dihedral_df = dihedral_df.drop(indices_to_drop).reset_index(drop=True)        
    for i,value in enumerate(dihedral_df.index):
        dihedral_df.at[i, "dihedral_index"] = i + 1
    with open("C:\\Users\\Janitha\\Documents\\GitHub\\reaction_lammps_mupt\\cache\\lunar\\bond_react_merge\\modified_dihedrals.txt", "w") as f:
        f.write(dihedral_df.to_string(index=False))
    dihedral_section = ""
    for i, row in dihedral_df.iterrows():
        dihedral_section += f"{row['dihedral_index']:>4}{row['dihedral_type']:>4}{row['new_atom1_index']:>4}{row['new_atom2_index']:>4}{row['new_atom3_index']:>4}{row['new_atom4_index']:>4}{row['hash']:>4}{row['atom1_type_real']:>4}{row['atom2_type_real']:>4}{row['atom3_type_real']:>4}{row['atom4_type_real']:>4}\n"
    number_of_dihedrals = len(dihedral_df)
    return dihedral_section, number_of_dihedrals

def modify_impropers(lines, type_df, improper_start_index):
    angle_data = []
    for line in lines[improper_start_index:]:
        if not line.strip(): # More Pythonic way to check for empty strings
            break
        parts = line.split()
        # Check if the line has enough parts to avoid an IndexError
        if len(parts) >= 11:
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
            # Using the already converted variables directly
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
    improper_df = pd.DataFrame(angle_data)
    improper_df["new_atom1_index"] = None 
    improper_df["new_atom2_index"] = None 
    improper_df["new_atom3_index"] = None 
    improper_df["new_atom4_index"] = None
    for value in improper_df["atom1_index"]:
        if value in type_df["atom_index"].values:
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            improper_df.loc[improper_df["atom1_index"] == value, "new_atom1_index"] = new_index
    for value in improper_df["atom2_index"]:
        if value in type_df["atom_index"].values:
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            improper_df.loc[improper_df["atom2_index"] == value, "new_atom2_index"] = new_index
    for value in improper_df["atom3_index"]:
        if value in type_df["atom_index"].values:
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            improper_df.loc[improper_df["atom3_index"] == value, "new_atom3_index"] = new_index
    for value in improper_df["atom4_index"]:
        if value in type_df["atom_index"].values:
            new_index = type_df.loc[type_df["atom_index"] == value, "new_atom_index"].values[0]
            improper_df.loc[improper_df["atom4_index"] == value, "new_atom4_index"] = new_index
    indices_to_drop = []        
    for index, row in improper_df.iterrows():
        if (row["new_atom1_index"] is None or 
            row["new_atom2_index"] is None or 
            row["new_atom3_index"] is None or
            row["new_atom4_index"] is None):
            indices_to_drop.append(index)
    improper_df = improper_df.drop(indices_to_drop).reset_index(drop=True)        
    for i,value in enumerate(improper_df.index):
        improper_df.at[i, "improper_index"] = i + 1
    with open("C:\\Users\\Janitha\\Documents\\GitHub\\reaction_lammps_mupt\\cache\\lunar\\bond_react_merge\\modified_impropers.txt", "w") as f:
        f.write(improper_df.to_string(index=False))
    improper_section = ""
    for i, row in improper_df.iterrows():
        improper_section += f"{row['improper_index']:>4}{row['improper_type']:>4}{row['new_atom1_index']:>4}{row['new_atom2_index']:>4}{row['new_atom3_index']:>4}{row['new_atom4_index']:>4}{row['hash']:>4}{row['atom1_type_real']:>4}{row['atom2_type_real']:>4}{row['atom3_type_real']:>4}{row['atom4_type_real']:>4}\n"
    number_of_impropers = len(improper_df)
    return improper_section, number_of_impropers

def molecule_file_format(number_of_types, number_of_bonds, number_of_angles, 
                         number_of_dihedrals, number_of_impropers, types_section, charge_section, coord_section, 
                         bond_section, angle_section, dihedral_section, improper_section):
    modified_molecule_file = f"""HEADER,  pre1.mol  read w/ mol2lmp > atom_typing v1.11 / 14 April 2025 w/PCFF-IFF atom-types > all2lmp: v1.22 / 14 April 2025  Class: 2 > bond_react_merge: v1.17 / 14 April 2025 molecule file (filetag: pre1)

{number_of_types} atoms
{number_of_bonds} bonds
{number_of_angles} angles
{number_of_dihedrals} dihedrals
{number_of_impropers} impropers

Types

"""
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
    lines, type_start_index, charge_start_index, coord_start_index, bond_start_index, angle_start_index, dihedral_start_index, improper_start_index = load_molecule_file(test_molecule_file)
    df_types, types_section, number_of_types = modify_types(lines, template_indexes, type_start_index)
    charge_section = modify_charges(lines, df_types, charge_start_index)
    coord_section = modify_coords(lines, df_types, coord_start_index)
    bond_section, number_of_bonds = modify_bonds(lines, df_types, bond_start_index)
    angle_section, number_of_angles = modify_angles(lines, df_types, angle_start_index)
    dihedral_section, number_of_dihedrals = modify_dihedrals(lines, df_types, dihedral_start_index)
    improper_section, number_of_impropers = modify_impropers(lines, df_types, improper_start_index)
    modified_molecule_file = molecule_file_format(number_of_types, number_of_bonds, number_of_angles, 
                                                 number_of_dihedrals, number_of_impropers, types_section, charge_section, coord_section, 
                                                 bond_section, angle_section, dihedral_section, improper_section)
    return modified_molecule_file

def load_molecule_files(file_dict, reactant_to_product):
    """
    TODO: Make these changes for multiple reactions by looping through reactant_to_product dictionary
    1. For each reactant and product file, read the molecule file.
    """
    template_indexes_reactant = []
    template_indexes_product = []
    for key, value in reactant_to_product.items(): 
        template_indexes_reactant.append(int(key+1)) # Convert to 1-based indexing which used in LAMMPS files in mol files they are 0-based
        template_indexes_product.append(int(value+1))  # Convert to 1-based indexing which used in LAMMPS files in mol files they are 0-based
    for name, filepath in file_dict.items():
        if name.startswith('pre'):
            num = get_ending_integer(name)
            molecule_file = file_dict.get(f'pre{num}')
            file_name = f'pre{num}'
        elif name.startswith('post'):
            num = get_ending_integer(name)
            molecule_file = file_dict[f'post{num}']
            file_name = f'post{num}'
        else:
            continue
        if molecule_file is None:
            raise ValueError(f"Data file for {name} not found.")
        modified_molecule = molecule_file_preparation(filepath, template_indexes_reactant if name.startswith('pre') else template_indexes_product)
        save_path = os.path.dirname(filepath)
        with open(os.path.join(save_path, f"template_{file_name}.molecule"), "w") as f:
            f.write(modified_molecule)
            
if __name__ == "__main__":

    molecule_file_dict = {   
    "data1": "C:\\Users\\Janitha\\Documents\\GitHub\\reaction_lammps_mupt\\cache\\lunar\\bond_react_merge\\data1_typed_IFF_merged.data",
    "data2": "C:\\Users\\Janitha\\Documents\\GitHub\\reaction_lammps_mupt\\cache\\lunar\\bond_react_merge\\data2_typed_IFF_merged.data",
    "pre1": "C:\\Users\\Janitha\\Documents\\GitHub\\reaction_lammps_mupt\\cache\\lunar\\bond_react_merge\\pre1_typed_IFF_merged.lmpmol",
    "post1": "C:\\Users\\Janitha\\Documents\\GitHub\\reaction_lammps_mupt\\cache\\lunar\\bond_react_merge\\post1_typed_IFF_merged.lmpmol",
    "force_field_data": "C:\\Users\\Janitha\\Documents\\GitHub\\reaction_lammps_mupt\\cache\\lunar\\bond_react_merge\\force_field.data"
}

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
    template_indexes_reactant = []
    template_indexes_product = []
    for key, value in reactant_to_product.items(): 
        template_indexes_reactant.append(int(key+1)) # Convert to 1-based indexing which used in LAMMPS files in mol files they are 0-based
        template_indexes_product.append(int(value+1))  # Convert to 1-based indexing which used in LAMMPS files in mol files they are 0-based
    print(template_indexes_reactant)
    print(template_indexes_product)
    load_molecule_files(molecule_file_dict, reactant_to_product)

