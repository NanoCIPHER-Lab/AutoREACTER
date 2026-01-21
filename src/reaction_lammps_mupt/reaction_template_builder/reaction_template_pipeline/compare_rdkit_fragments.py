from rdkit import Chem

def extract_fragment_by_indices(mol, atom_indices_to_keep):
    """
    Removes a list of atoms from an RDKit molecule by their indices.
    """
    # 1. Convert to an RWMol object for editing
    rwmol = Chem.RWMol(mol)

    # 2. Sort indices in reverse (descending) order
    all_indices = [a.GetIdx() for a in mol.GetAtoms()]
    atom_indices_to_remove = [idx for idx in all_indices if idx not in atom_indices_to_keep]

    sorted_indices = sorted(atom_indices_to_remove, reverse=True)

    # 3. Iterate and remove atoms
    for idx in sorted_indices:
        rwmol.RemoveAtom(idx) # or rwmol.RemoveAtomIdx(idx)

    # 4. Convert back to a standard Mol object and sanitize
    # The sanitization step is important to fix valence and bond info
    new_mol = rwmol.GetMol()
    Chem.SanitizeMol(new_mol)
    
    return new_mol

def cap_open_valences_with_fr(new_mol_fragment, francium_atomic_num = 87):
    """
    Cap all open valences in an RDKit molecular fragment with placeholder atoms.

    This function scans the given fragment for atoms whose current valence is lower
    than their default valence and adds new atoms with the specified atomic number,
    connected by single bonds, to satisfy the open valences.

    :param new_mol_fragment: RDKit Mol fragment whose open valences will be capped.
    :param fr_atomic_num: Atomic number to use for the capping atoms (default: 87).
    :return: A dict with keys:
        - "object": the capped RDKit Mol object
        - "smiles": SMILES string of the capped molecule
        - "inchi": InChI string of the capped molecule
    """
    # 1. Convert to RWMol for editing
    rw_mol = Chem.RWMol(new_mol_fragment)
    # 2. Iterate over atoms to find open valences
    for atom in rw_mol.GetAtoms():
        atom_idx = atom.GetIdx()
        # 3. Determine open valences
        default_valence = Chem.GetPeriodicTable().GetDefaultValence(atom.GetAtomicNum())
        current_valence = sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds())
        try:
            base_valence = default_valence[0] if isinstance(default_valence, tuple) else default_valence
            open_valences = max(0, int(base_valence) - int(current_valence))
        except (TypeError, ValueError):
            # If we cannot determine a sensible integer valence, skip capping for this atom
            continue
        if open_valences == 0:
            continue

        for _ in range(open_valences):
            # 4. Add a "fr" atom (using atomic number for Francium as a placeholder)
            fr_atom = Chem.Atom(fr_atomic_num)
            fr_idx = rw_mol.AddAtom(fr_atom)
            rw_mol.AddBond(atom_idx, fr_idx, Chem.BondType.SINGLE)
    # 5. Convert back to Mol and sanitize
    capped_mol = rw_mol.GetMol()
    Chem.SanitizeMol(capped_mol)
    capped_mol_info = {
        "object": capped_mol,
        "smiles": Chem.MolToSmiles(capped_mol)    ,
        "inchi": Chem.MolToInchi(capped_mol),
    }  
    return capped_mol_info

def compare_fragments(mol1_info, mol2_info):
    """
    Compares two RDKit molecules and returns True if they are identical,
    otherwise returns False.
    """
    # Check for None inputs
    if not mol1_info or not mol2_info:
        return False
    
    # Compare by SMILES or InChI
    if mol1_info.get("smiles") == mol2_info.get("smiles"):
        return True

    if mol1_info.get("inchi") == mol2_info.get("inchi"):
        return True
    return False

def compare_rdkit_fragments(processed_dict, new_fragment_dict):
    """
    Check whether a new RDKit fragment has already been processed.

    This function compares a candidate fragment (``new_fragment_dict``) against
    all fragments stored in ``processed_dict`` using :func:`compare_fragments`.
    If an identical fragment is found, it returns ``True`` and leaves
    ``processed_dict`` unchanged. If no match is found, the new fragment is
    appended to ``processed_dict`` under a new integer key
    (``len(processed_dict) + 1``), and the function returns ``(False, processed_dict)``.

    :param dict processed_dict: Mapping from fragment identifiers (typically
        integers) to fragment information dictionaries. Each fragment info
        dictionary is expected to contain at least ``"smiles"`` and/or
        ``"inchi"`` entries, as produced by upstream processing functions.
    :param dict new_fragment_dict: Fragment information dictionary for the
        candidate fragment to compare, in the same format as the values in
        ``processed_dict``.

    :returns: If an identical fragment is already present in
        ``processed_dict``, returns ``True``. If the fragment is new, it is
        added to ``processed_dict`` and the function returns a tuple
        ``(False, processed_dict)``.
    :rtype: bool or (bool, dict)
    """
    for proc_id, proc_frag_info in processed_dict.items():
        if compare_fragments(proc_frag_info, new_fragment_dict):
            return True, processed_dict
    processed_dict[len(processed_dict) + 1] = new_fragment_dict
    return False, processed_dict



