from rdkit import Chem

def compare_products(reactions_dict, _prod2):
    """
    Compares a new product molecule against a collection of existing products 
    to identify duplicates, ignoring atom mapping.

    Args:
        reactions_dict (dict): A dictionary where values are dictionaries 
                               potentially containing a "product" RDKit molecule object.
        _prod2 (rdkit.Chem.rdchem.Mol): The RDKit molecule object to check for duplicates.

    Returns:
        bool: Returns False if a duplicate product is found, True otherwise.
    """
    # Create a copy of the input molecule to avoid modifying the original object
    prod2 = Chem.Mol(_prod2)
    
    for key, dict_2 in reactions_dict.items():
        # Retrieve the product molecule from the nested dictionary
        prod = dict_2.get("product")
        
        if prod is not None:
            # Create a copy of the stored product for comparison
            prod1 = Chem.Mol(prod)

            # Remove atom map numbers from both molecules.
            # Atom maps are often used in reactions to track atoms but should be 
            # ignored when checking if two chemical structures are identical.
            for atom in prod1.GetAtoms():
                atom.SetAtomMapNum(0)
            for atom in prod2.GetAtoms():
                atom.SetAtomMapNum(0)

            # Convert both molecules to canonical SMILES strings.
            # Canonicalization ensures that the same molecule always results in the same string.
            smi1 = Chem.MolToSmiles(prod1, canonical=True)
            smi2 = Chem.MolToSmiles(prod2, canonical=True)

            # Perform string comparison to detect duplicates
            if smi1 == smi2:
                # Debug print to indicate a duplicate was found (can be commented out in production)
                # print("DUPLICATE PRODUCT FOUND")
                return False
                
    # No duplicates found after checking all entries
    return True

def compare_rdkit_molecules_canonical(data_smiles_list, mol_smi_2):
    """
    Compares two RDKit molecule SMILES strings to determine if they
    represent the same chemical structure using canonical SMILES.

    Args:
        data_smiles_list (list): List of SMILES strings to compare against.
        mol_smi_2 (str): SMILES string of the second molecule.

    Returns:
        bool: True if molecules are chemically identical, False otherwise.
    """
    if not mol_smi_2:
        return data_smiles_list, False
    mol2 = Chem.MolFromSmiles(mol_smi_2)
    if mol2 is None:
        return data_smiles_list, False
        
    for mol_smi_1 in data_smiles_list:
        mol1 = Chem.MolFromSmiles(mol_smi_1)
        
        # Handle cases where SMILES might be invalid
        if mol1 is None or mol2 is None:
            return data_smiles_list, False  # or raise a ValueError
        # Generate canonical SMILES and compare them
        canonical_smi_1 = Chem.MolToSmiles(mol1, canonical=True)
        canonical_smi_2 = Chem.MolToSmiles(mol2, canonical=True)

        if canonical_smi_1 == canonical_smi_2:
            return data_smiles_list, True
    
    data_smiles_list.append(mol_smi_2)
    return data_smiles_list, False

def format_detected_reactions_dict(detected_reactions, non_monomer_molecules_to_retain = None):
    """
    Aggregates and formats information from a dictionary of detected reactions 
    into a single summary dictionary.

    This function compiles unique reaction names, SMARTS references, and 
    mechanism references into comma-separated strings.

    Args:
        detected_reactions (dict): A dictionary containing reaction data, 
                                   where each value is a dictionary with keys 
                                   like "reaction_name" and "reference".

    Returns:
        dict: A dictionary containing aggregated strings for:
              - "reactions_names"
              - "smart_references"
              - "mechanism_references"
    """
    reactions_names = []
    smart_references = []
    mechanism_references = []

    seen_reactions = set()
    seen_smarts = set()
    seen_mechanisms = set()
    
    formatted_dict = {}
    data_smiles_list = []
    
    for key, reaction in detected_reactions.items():
        # Deduplicate reaction names safely (no substring bugs)
        reactions_name = reaction.get("reaction_name")
        if reactions_name and reactions_name not in seen_reactions:
            seen_reactions.add(reactions_name)
            reactions_names.append(reactions_name)
    
        # Reference is optional; do not skip monomer processing if missing
        reference = reaction.get("reference") or {}
        smarts_ref = reference.get("smarts")
        mech_refs = reference.get("reaction_and_mechanism")
    
        if smarts_ref and smarts_ref not in seen_smarts:
            seen_smarts.add(smarts_ref)
            smart_references.append(smarts_ref)
    
        if mech_refs:
            mech_block = ", ".join(mech_refs)
            if mech_block and mech_block not in seen_mechanisms:
                seen_mechanisms.add(mech_block)
                mechanism_references.append(mech_block)
    
        # Always process monomers (if present)
        m1 = reaction.get("monomer_1")
        m2 = reaction.get("monomer_2")
        m1_smiles = m1.get("smiles") if isinstance(m1, dict) else None
        m2_smiles = m2.get("smiles") if isinstance(m2, dict) else None
    
        if m1_smiles:
            data_smiles_list, _ = compare_rdkit_molecules_canonical(data_smiles_list, m1_smiles)
        if m2_smiles:
            data_smiles_list, _ = compare_rdkit_molecules_canonical(data_smiles_list, m2_smiles)
    
    formatted_dict = {
        "reactions_names": ", ".join(reactions_names),
        "smart_references": ", ".join(smart_references),
        "mechanism_references": ", ".join(mechanism_references),
    }
    return formatted_dict , data_smiles_list

def prep_for_3d_molecule_generation(data_smiles_list, molecule_dict_csv_path_dict):
    molecule_dict = {}
    for i, smiles in enumerate(data_smiles_list):
        key = f"data_{i+1}"
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        molecule = Chem.AddHs(molecule)
        molecule_dict[key] = molecule
    for i, (key, value) in enumerate(molecule_dict_csv_path_dict.items()):
        reactant = value.get("reactant")
        product = value.get("product")
        if reactant and product:
            molecule_dict[f"pre_{i+1}"] = reactant
            molecule_dict[f"post_{i+1}"] = product
    return molecule_dict
