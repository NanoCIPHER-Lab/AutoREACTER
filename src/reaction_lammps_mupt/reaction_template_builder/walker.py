from rdkit import Chem
from rdkit.Chem import AllChem

def get_new_neighbors(molecule_object, atom_index, visited):
        neighbors = []
        atom = molecule_object.GetAtomWithIdx(atom_index)
        values_list = [x for v in visited.values() for x in v]  
        for neighbor in atom.GetNeighbors():
            # Fixed: check if neighbor is already in the visited values
            if neighbor.GetIdx() not in values_list:
                neighbors.append(neighbor.GetIdx())
        return neighbors

def reactant_atom_walker(reactant_molecule_object, start_atom_indexs, max_bonds=4):
    '''
    Walks the molecular graph starting from start_atom_index up to max_bonds away.
    Returns a dictionary with shells of neighbors.
    '''
    try:
        start_atom_indexs = list(start_atom_indexs)
    except TypeError:
        start_atom_indexs = [start_atom_indexs]
        pass
    # if isinstance(start_atom_index, int):
    #     start_atom_index = [start_atom_index]
    
    neighbor_index = 0
    # Initialize shell dictionary: 1=starting atom, 2-4=neighbor shells
    # Fixed: Created dynamic dict based on max_bonds to prevent KeyErrors
    template_indexes = {i: [] for i in range(1, max_bonds + 2)}
    template_indexes[1] = start_atom_indexs
    processed_indexes = []

    # Expand neighbors shell by shell
    while neighbor_index < max_bonds:
        neighbor_index += 1
        if neighbor_index == max_bonds:
            break
        
        # For each atom in current shell, find new neighbors
        # Fixed: iterate over a copy to be safe
        for neighbor in template_indexes[neighbor_index]:
            # Pass the whole dict so we check all history
            new_neighbors = get_new_neighbors(reactant_molecule_object, neighbor, template_indexes)
            if new_neighbors is None:
                return None
            
            # Add new neighbors to next shell
            # Fixed: Check duplicates before extending
            current_next_shell = template_indexes[neighbor_index + 1]
            for n in new_neighbors:
                if n not in current_next_shell:
                    template_indexes[neighbor_index + 1].append(n)
            
            processed_indexes.append(neighbor)
    edge_atoms = template_indexes[max_bonds]
    reactant_template_indexes = [x for v in template_indexes.values() for x in v]
    
    # Flatten all shells into single list of atom indices
    return reactant_template_indexes, edge_atoms

def product_atom_walker(template_indexes, MAP_dict):
    template_mapped_dict = {}
    for i in template_indexes:
        atom = MAP_dict.get(i)
        if atom:
            template_mapped_dict[i] = atom
    return template_mapped_dict
            
def reaction_atom_walker(molecule_object, start_atom_indexs, MAP_dict, max_bonds=4):
    '''
    Walks the molecular graph starting from start_atom_index up to max_bonds away.
    Returns a dictionary with shells of neighbors mapped to product atom indices.
    '''
    reactant_template_indexes, edge_atoms = reactant_atom_walker(molecule_object, start_atom_indexs, max_bonds)
    template_mapped_dict = product_atom_walker(reactant_template_indexes, MAP_dict)
    return template_mapped_dict, edge_atoms

if __name__ == "__main__":
    reactant_smiles1 = "C1=CC=C(C(=C1)C(=O)O)O" 
    reactant_smiles2 = "OCCC(O)=O" 
    reaction_smarts = "[O;!$(OC=*):1]-[H:3].[CX3:2](=[O:5])[OX2H1:4]>>[OX2:1]-[CX3:2](=[O:5]).[O:4]-[H:3]"
    mapped_atoms = {9: 0, 19: 1, 3: 3, 15: 26, 18: 17, 20: 25, 21: 2, 2: 4, 4: 5, 17: 18, 25: 19, 26: 20, 27: 27, 1: 6, 12: 7, 5: 8, 6: 9, 16: 21, 23: 22, 24: 23}
    # Load and add hydrogens to reactants
    reactant1 = Chem.MolFromSmiles(reactant_smiles1)
    reactant1 = Chem.AddHs(reactant1)
    reactant2 = Chem.MolFromSmiles(reactant_smiles2)
    reactant2 = Chem.AddHs(reactant2)
    
    # Create reaction from SMARTS
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)
    products = rxn.RunReactants((reactant1, reactant2))
    delete_atom = True  # Expect byproduct (water)
    combined_reactants = Chem.CombineMols(reactant1, reactant2)
    products = products[0]
    combined_products = Chem.CombineMols(*products)
    template_mapped_dict, edge_atoms = reaction_atom_walker(combined_reactants, [9, 19], mapped_atoms)
    import json
    print("Reactant Template Indexes:", json.dumps(template_mapped_dict, indent=2))
    print("Edge Atoms in Reactant:", edge_atoms)
    