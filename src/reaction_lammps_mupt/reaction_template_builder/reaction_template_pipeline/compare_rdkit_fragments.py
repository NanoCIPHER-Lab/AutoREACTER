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
    Create two lists of integer atom indices from a mapping of reactant-to-product atom indices.
    
    Parameters:
        input_dict (dict): Mapping where keys are reactant atom indices and values are product atom indices.
    
    Returns:
        tuple: (reactant_indices, product_indices) — two lists of integers corresponding to the reactant and product atom indices.
    """
    # Convert keys and values to integers to ensure consistent indexing
    reactant_indices = [int(k) for k in input_dict.keys()]
    product_indices = [int(v) for v in input_dict.values()]
    return reactant_indices, product_indices


def extract_fragment_by_indices(mol, atom_indices_to_keep):
    """
    Create a molecule containing only the atoms at the specified indices.
    
    Parameters:
        mol (rdkit.Chem.rdchem.Mol): Source RDKit molecule to fragment.
        atom_indices_to_keep (list[int]): Atom indices to retain in the resulting fragment.
    
    Returns:
        rdkit.Chem.rdchem.Mol or None: The extracted fragment as an RDKit Mol, or `None` if `mol` is `None`. A sanitization attempt is made; if sanitization fails the raw (unsanitized) fragment is returned.
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
    Attach Francium placeholder atoms to any atoms in the fragment that have unsatisfied valences.
    
    Parameters:
        new_mol_fragment (rdkit.Chem.rdchem.Mol): RDKit molecule fragment whose open valences should be capped.
        francium_atomic_num (int): Atomic number to use for placeholder atoms (default 87 for Fr).
    
    Returns:
        dict: {
            "object": rdkit.Chem.rdchem.Mol or None — the capped molecule object (None if input was None or empty),
            "smiles": str — SMILES of the capped molecule (empty string if input was None or empty),
            "inchi": str — InChI of the capped molecule (empty string if input was None or empty)
        }
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
        except:
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
    Determine whether two molecular info dictionaries represent the same chemical structure.
    
    Parameters:
        mol1_info (dict): Mapping that may contain "smiles" and/or "inchi" string entries.
        mol2_info (dict): Mapping that may contain "smiles" and/or "inchi" string entries.
    
    Returns:
        True if the SMILES or InChI strings are equal, False otherwise.
    """
    if not mol1_info or not mol2_info:
        return False
    
    # Check SMILES identity
    if mol1_info.get("smiles") == mol2_info.get("smiles"):
        return True

    # Check InChI identity (more robust for tautomers/stereoisomers in some cases)
    if mol1_info.get("inchi") == mol2_info.get("inchi"):
        return True
        
    return False

def compare_rdkit_fragments(processed_dict, combined_reactant_mol, combined_product_mol, template_mapped_dict):
    """
    Extract reactant and product fragments for a mapped reaction, cap open valences, and determine whether the transformation has been seen before.
    
    Processes the mapped atom indices from template_mapped_dict to build raw and francium-capped fragment representations for both reactant and product, compares those representations against entries in processed_dict, and appends a new history entry when no match is found.
    
    Parameters:
        processed_dict (dict): History mapping IDs to entries with keys "reactant_info" and "product_info".
        combined_reactant_mol (rdkit.Chem.rdchem.Mol): Full reactant molecule to extract the reactant fragment from.
        combined_product_mol (rdkit.Chem.rdchem.Mol): Full product molecule to extract the product fragment from.
        template_mapped_dict (dict): Mapping of atom indices involved in the reaction (reactant index -> product index).
    
    Returns:
        tuple: `True` if the reactant/product fragment pair already exists in processed_dict, `False` otherwise; the second element is the (possibly updated) processed_dict.
    """
    # Create deep copies to avoid modifying the original molecules in memory
    react_mol = copy.deepcopy(combined_reactant_mol)
    prod_mol = copy.deepcopy(combined_product_mol)
    
    # Get the indices of the atoms involved in the reaction center
    react_indices, prod_indices = dict_keys_to_list(template_mapped_dict)
    
    # Process Reactant Fragment
    raw_react_frag = extract_fragment_by_indices(react_mol, react_indices)
    raw_react_info = {
        "smiles": Chem.MolToSmiles(raw_react_frag) if raw_react_frag else "", 
        "inchi": Chem.MolToInchi(raw_react_frag) if raw_react_frag else ""
    }
    capped_react_info = cap_open_valences_with_fr(raw_react_frag)
    
    # Process Product Fragment
    raw_prod_frag = extract_fragment_by_indices(prod_mol, prod_indices)
    raw_prod_info = {
        "smiles": Chem.MolToSmiles(raw_prod_frag) if raw_prod_frag else "", 
        "inchi": Chem.MolToInchi(raw_prod_frag) if raw_prod_frag else ""
    }
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