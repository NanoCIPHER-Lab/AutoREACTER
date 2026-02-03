"""
Chemical Fragment Extraction and Comparison Utility
This module provides functions to extract specific fragments from RDKit molecules,
cap open valences with placeholder atoms (Francium), and compare these fragments
against a history of processed structures to identify unique chemical transformations.
"""

from rdkit import Chem

import copy

def dict_keys_to_list(input_dict):
    """
    Converts a dictionary of atom mappings into two separate lists of indices.

    Args:
        input_dict (dict): A dictionary where keys represent reactant atom indices 
                           and values represent product atom indices.

    Returns:
        tuple: (reactant_indices, product_indices) as lists of integers.
    """
    # Convert keys and values to integers to ensure consistent indexing
    reactant_indices = [int(k) for k in input_dict.keys()]
    product_indices = [int(v) for v in input_dict.values()]
    return reactant_indices, product_indices


def extract_fragment_by_indices(mol, atom_indices_to_keep):
    """
    Creates a new molecule containing only the atoms specified by the provided indices.
    All other atoms are removed.

    Args:
        mol (rdkit.Chem.rdchem.Mol): The source RDKit molecule.
        atom_indices_to_keep (list): List of atom indices to retain in the fragment.

    Returns:
        rdkit.Chem.rdchem.Mol: The extracted molecular fragment, or None if input is invalid.
    """
    if mol is None:
        return None

    # 1. Convert to an RWMol (Read-Write Molecule) object to allow structural editing
    rwmol = Chem.RWMol(mol)

    # 2. Identify atoms to remove
    # We find the difference between all atoms in the molecule and the ones we want to keep
    all_indices = [a.GetIdx() for a in mol.GetAtoms()]
    atom_indices_to_remove = [idx for idx in all_indices if idx not in atom_indices_to_keep]

    # 3. Sort indices in reverse (descending) order
    # Crucial: Removing atoms shifts the indices of subsequent atoms. 
    # Removing from highest index to lowest prevents this shifting issue.
    sorted_indices = sorted(atom_indices_to_remove, reverse=True)

    # 4. Iterate and remove atoms from the RWMol
    for idx in sorted_indices:
        rwmol.RemoveAtom(idx) 

    # Convert back to a standard Molecule object
    new_mol = rwmol.GetMol()
    
    # 5. Attempt Sanitization
    # Fragments often have "broken" valences. We try to sanitize to refresh 
    # molecular properties, but wrap it in a try-except to prevent valence errors
    # from crashing the execution.
    try:
        Chem.SanitizeMol(new_mol)
    except Exception:
        # If sanitization fails (e.g., due to invalid valences), we proceed with the raw fragment
        pass
    
    return new_mol

def cap_open_valences_with_fr(new_mol_fragment, francium_atomic_num=87):
    """
    Identifies atoms with unsatisfied valences and attaches Francium (Fr) atoms 
    as placeholders. This is useful for maintaining structural context in fragments.

    Args:
        new_mol_fragment (rdkit.Chem.rdchem.Mol): The fragment to cap.
        francium_atomic_num (int): The atomic number to use for capping (default 87 for Fr).

    Returns:
        dict: A dictionary containing the RDKit object, SMILES string, and InChI string.
    """
    # Handle empty or null fragments
    if new_mol_fragment is None or new_mol_fragment.GetNumAtoms() == 0:
        return {"object": None, "smiles": "", "inchi": ""}

    rw_mol = Chem.RWMol(new_mol_fragment)
    
    # Use a static list of atoms to avoid iterator invalidation during modification
    atoms = list(rw_mol.GetAtoms())
    
    for atom in atoms:
        atom_idx = atom.GetIdx()
        try:
            # Determine the expected valence for the atom type
            default_valence = Chem.GetPeriodicTable().GetDefaultValence(atom.GetAtomicNum())
            
            # Calculate the current valence based on existing bonds
            current_valence = sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds())
            
            # Handle cases where default valence might be a tuple (multiple oxidation states)
            base_valence = default_valence[0] if isinstance(default_valence, tuple) else default_valence
            
            # Calculate how many bonds are "missing"
            open_valences = max(0, int(base_valence) - int(current_valence))
            
            # Add Francium atoms for each open valence
            if open_valences > 0:
                for _ in range(open_valences):
                    fr_atom = Chem.Atom(francium_atomic_num)
                    fr_idx = rw_mol.AddAtom(fr_atom)
                    # Connect the placeholder atom with a single bond
                    rw_mol.AddBond(atom_idx, fr_idx, Chem.BondType.SINGLE)
        except Exception:
            # Skip atoms where valence cannot be determined (e.g., certain metals)
            continue

    # Finalize the molecule after capping
    capped_mol = rw_mol.GetMol()
    try:
        Chem.SanitizeMol(capped_mol)
    except Exception as e:
        # Some fragments (e.g., highly unusual or intentionally "broken" structures)
        # may fail RDKit sanitization. We keep the unsanitized molecule to preserve
        # behavior, but log the issue for easier debugging.
        print(f"Warning: RDKit sanitization failed for capped fragment: {e}")

    return {
        "object": capped_mol,
        "smiles": Chem.MolToSmiles(capped_mol),
        "inchi": Chem.MolToInchi(capped_mol),
    }

def compare_fragments(mol1_info, mol2_info):
    """
    Compares two molecular info dictionaries to determine if they represent 
    the same chemical structure.

    Args:
        mol1_info (dict): Dictionary containing 'smiles' and/or 'inchi'.
        mol2_info (dict): Dictionary containing 'smiles' and/or 'inchi'.

    Returns:
        bool: True if molecules match by SMILES or InChI, False otherwise.
    """
    if not mol1_info or not mol2_info:
        return False
    
    # Check SMILES identity
    smiles1, smiles2 = mol1_info.get("smiles"), mol2_info.get("smiles")
    if smiles1 and smiles2 and smiles1 == smiles2:
        return True

    # Check InChI identity (more robust for tautomers/stereoisomers in some cases)
    inchi1, inchi2 = mol1_info.get("inchi"), mol2_info.get("inchi")
    if inchi1 and inchi2 and inchi1 == inchi2:
        return True
    return False

def compare_rdkit_fragments(processed_dict, combined_reactant_mol, combined_product_mol, template_mapped_dict):
    """
    Main logic to extract fragments from a reaction and check if this specific 
    transformation has been encountered before.

    Args:
        processed_dict (dict): A history of previously seen fragments.
        combined_reactant_mol (rdkit.Chem.rdchem.Mol): The full reactant molecule.
        combined_product_mol (rdkit.Chem.rdchem.Mol): The full product molecule.
        template_mapped_dict (dict): Mapping of atom indices involved in the reaction.

    Returns:
        tuple: (bool, updated_processed_dict). True if the fragment pair was already known.
    """
    # Create deep copies to avoid modifying the original molecules in memory
    react_mol = copy.deepcopy(combined_reactant_mol)
    prod_mol = copy.deepcopy(combined_product_mol)
    
    # Get the indices of the atoms involved in the reaction center
    react_indices, prod_indices = dict_keys_to_list(template_mapped_dict)
    
    # Process Reactant Fragment
    raw_react_frag = extract_fragment_by_indices(react_mol, react_indices)
    try:
        react_smiles = Chem.MolToSmiles(raw_react_frag) if raw_react_frag else ""
    except Exception:
        react_smiles = ""
    try:
        react_inchi = Chem.MolToInchi(raw_react_frag) if raw_react_frag else ""
    except Exception:
        react_inchi = ""
    raw_react_info = {"smiles": react_smiles, "inchi": react_inchi}
    capped_react_info = cap_open_valences_with_fr(raw_react_frag)
    
    # Process Product Fragment
    raw_prod_frag = extract_fragment_by_indices(prod_mol, prod_indices)
    try:
        prod_smiles = Chem.MolToSmiles(raw_prod_frag) if raw_prod_frag else ""
    except Exception:
        prod_smiles = ""
    try:
        prod_inchi = Chem.MolToInchi(raw_prod_frag) if raw_prod_frag else ""
    except Exception:
        prod_inchi = ""
    raw_prod_info = {"smiles": prod_smiles, "inchi": prod_inchi}
    capped_prod_info = cap_open_valences_with_fr(raw_prod_frag)

    # Compare against history
    for proc_id, history in processed_dict.items():
        hist_react = history['reactant_info']
        hist_prod = history['product_info']
        
        # Check if current reactant matches history (either raw or capped)
        react_match = compare_fragments(hist_react, raw_react_info) or \
                      compare_fragments(hist_react, capped_react_info)
        
        # Check if current product matches history (either raw or capped)
        prod_match = compare_fragments(hist_prod, raw_prod_info) or \
                     compare_fragments(hist_prod, capped_prod_info)

        # If both reactant and product fragments match an entry, it's a duplicate
        if react_match and prod_match:
            return True, processed_dict
        
    # If it's a new transformation, add it to the history
    new_id = len(processed_dict) + 1
    processed_dict[new_id] = {
        "reactant_info": capped_react_info,
        "product_info": capped_prod_info
    }
    
    return False, processed_dict
