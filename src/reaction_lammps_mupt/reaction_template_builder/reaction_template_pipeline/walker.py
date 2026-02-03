"""
Molecular Environment Extractor
-------------------------------
This script provides utilities for walking through a molecular graph starting from 
specific atom indices. It is designed to identify local chemical environments 
(shells of neighbors) and map these atoms between reactant and product states 
using RDKit.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import json

def get_new_neighbors(molecule_object, atom_index, visited):
    """
    Return neighbor atom indices of a given atom that are not present in any previously visited shells.
    
    Parameters:
        molecule_object (rdkit.Chem.rdchem.Mol): RDKit molecule containing the atom.
        atom_index (int): Index of the atom whose neighbors are queried.
        visited (dict): Dictionary whose values are lists of atom indices already visited.
    
    Returns:
        list: Neighbor atom indices (int) that do not appear in any of the lists in `visited`.
    """
    neighbors = []
    # Retrieve the atom object by its index
    atom = molecule_object.GetAtomWithIdx(atom_index)
    
    # Flatten all indices currently stored in the visited dictionary to check against
    values_list = [x for v in visited.values() for x in v]  
    
    for neighbor in atom.GetNeighbors():
        # Only add the neighbor if its index has not been encountered in any shell yet
        if neighbor.GetIdx() not in values_list:
            neighbors.append(neighbor.GetIdx())
    return neighbors

def reactant_atom_walker(combined_reactant_molecule_object, start_atom_indexes, max_bonds=4):
    """
    Perform a breadth-first walk of the molecule from given seed atom indices up to a specified bond distance.
    
    Parameters:
        combined_reactant_molecule_object (rdkit.Chem.rdchem.Mol): Molecule to traverse.
        start_atom_indexes (int or iterable[int]): Starting atom index or iterable of start indices.
        max_bonds (int): Maximum bond steps (shell depth) to explore.
    
    Returns:
        reactant_template_indexes (list): Flattened list of all atom indices encountered during the walk.
        edge_atoms (list): Atom indices that are exactly at the specified max_bonds distance.
    """
    # Ensure start_atom_indexes is a list even if a single integer is passed
    try:
        start_atom_indexes = list(start_atom_indexes)
    except TypeError:
        start_atom_indexes = [start_atom_indexes]
    
    neighbor_index = 0
    
    # Initialize the shell dictionary. 
    # Key 1 contains starting atoms, keys 2 to max_bonds+1 will contain subsequent shells.
    template_indexes = {i: [] for i in range(1, max_bonds + 1)}
    template_indexes[1] = start_atom_indexes
    processed_indexes = []

    # Expand the search shell by shell
    while neighbor_index < max_bonds - 1: # Loop until reaching max_bonds
        neighbor_index += 1 

        # Iterate through atoms found in the current shell to find their neighbors
        for atom_idx in template_indexes[neighbor_index]:
            # Get neighbors not present in any previous shell
            new_neighbors = get_new_neighbors(combined_reactant_molecule_object, atom_idx, template_indexes)
            
            # Add unique new neighbors to the next shell
            current_next_shell = template_indexes[neighbor_index + 1]
            for n in new_neighbors:
                if n not in current_next_shell:
                    template_indexes[neighbor_index + 1].append(n)
            
            processed_indexes.append(atom_idx)
    print(f"template_indexes: {template_indexes}")
    # Atoms at the furthest distance reached to create the edge atom section in the map file
    edge_atoms = template_indexes[max_bonds]
    print("Edge atoms at max distance shell:", edge_atoms)
    # Flatten all shells into a single list of indices representing the local environment
    reactant_template_indexes = [x for v in template_indexes.values() for x in v]
    return reactant_template_indexes, edge_atoms

def product_atom_walker(template_indexes, MAP_dict):
    """
    Maps a list of reactant atom indices to their corresponding indices in the product 
    molecule using a mapping dictionary.

    Args:
        template_indexes (list): List of reactant atom indices.
        MAP_dict (dict): Dictionary mapping reactant indices (keys) to product indices (values).

    Returns:
        dict: A dictionary of {reactant_index: product_index} for only template atoms.
    """
    template_mapped_dict = {}
    for i in template_indexes:
        print(f"Mapping reactant atom index {i}")
        # Check if the reactant atom index exists in our mapping dictionary
        atom = MAP_dict.get(i)
        if atom is not None:
            template_mapped_dict[i] = atom
    return template_mapped_dict
            
def reaction_atom_walker(combined_reactant_molecule_object, start_atom_indexes, MAP_dict, max_bonds=4):
    """
    Orchestrates a local-environment walk on the combined reactant molecule and maps the discovered reactant atom indices to product atom indices.
    
    Parameters:
        combined_reactant_molecule_object (rdkit.Chem.rdchem.Mol): Combined reactant molecule used for the graph walk.
        start_atom_indexes (list|int): Starting atom index or iterable of starting atom indices that seed the traversal.
        MAP_dict (dict): Mapping from reactant atom indices to product atom indices.
        max_bonds (int): Maximum bond distance (shell depth) to explore.
    
    Returns:
        tuple: (template_mapped_dict, edge_atoms)
            template_mapped_dict (dict): Mapping of reactant atom indices found in the walked template to their corresponding product atom indices.
            edge_atoms (list): Atom indices located exactly at the specified `max_bonds` distance from the start atoms.
    """
    # 1. Identify the local environment in the reactant
    reactant_template_indexes, edge_atoms = reactant_atom_walker(combined_reactant_molecule_object, start_atom_indexes, max_bonds)
    
    # 2. Map those reactant atoms to product atoms
    template_mapped_dict = product_atom_walker(reactant_template_indexes, MAP_dict)
    
    return template_mapped_dict, edge_atoms

if __name__ == "__main__":
    # Define reactants and reaction SMARTS (Esterification example)
    reactant_smiles1 = "C1=CC=C(C(=C1)C(=O)O)O"  # Salicylic acid
    reactant_smiles2 = "OCCC(O)=O"                
    reaction_smarts = "[O;!$(OC=*):1]-[H:3].[CX3:2](=[O:5])[OX2H1:4]>>[OX2:1]-[CX3:2](=[O:5]).[O:4]-[H:3]"
    
    # Example mapping dictionary: Reactant Atom Index -> Product Atom Index.
    # This mapping was manually derived for this specific salicylic acid esterification
    # example by inspecting the atom indices in the combined reactant and product
    # RDKit molecules; it encodes the correspondence used in this demo workflow.
    mapped_atoms = {
        0: 10,
        1: 6,
        2: 4,
        3: 3,
        4: 5,
        5: 8,
        6: 9,
        7: 13,
        8: 14,
        9: 0,
        10: 15,
        11: 11,
        12: 7,
        13: 12,
        14: 16,
        15: 26,
        16: 21,
        17: 18,
        18: 17,
        19: 1,
        20: 25,
        21: 2,
        22: 24,
        23: 22,
        24: 23,
        25: 19,
        26: 20,
        27: 27,
    }
    # Load and add hydrogens to reactants to ensure full graph connectivity
    reactant1 = Chem.MolFromSmiles(reactant_smiles1)
    reactant1 = Chem.AddHs(reactant1)
    reactant2 = Chem.MolFromSmiles(reactant_smiles2)
    reactant2 = Chem.AddHs(reactant2)
    
    # Initialize the chemical reaction from SMARTS
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)
    
    # Run the reaction
    products = rxn.RunReactants((reactant1, reactant2))
    
    # Combine individual reactant molecules into one object for easier indexing
    combined_reactants = Chem.CombineMols(reactant1, reactant2)
    
    # Process the first set of products generated
    if products:
        products_set = products[0]
        combined_products = Chem.CombineMols(*products_set)
        
        # Walk the graph starting from atoms 9 and 19 (indices in the combined reactant)
        template_mapped_dict, edge_atoms = reaction_atom_walker(combined_reactants, [9, 19], mapped_atoms)
        
        # Display results
        print("Reactant to Product Mapping (Local Environment):")
        print(json.dumps(template_mapped_dict, indent=2))
        print("\nEdge Atoms in Reactant (Max Distance Shell):", edge_atoms)