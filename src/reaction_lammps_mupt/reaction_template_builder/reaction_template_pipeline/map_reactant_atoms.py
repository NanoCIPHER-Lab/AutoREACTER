"""
Still, there is no method to verify that all possible reactions can be implemented correctly.
So, for now, we will proceed with only one or two reactions.

TODO:
- Add a verification method for all reactions.

If there are multiple reactions between the same reactants, this is the place to handle that.

Example:
    reaction_smarts = "[O;!$(OC=*):1]-[H:2].[CX3:3](=[O])[OX2H1:4]>>[OX2:1]-[CX3:3](=[O]).[O:4]-[H:2]"
    reactant_smiles1 = "OC(F)c1cc(C(O)Cl)cc(C(O)Br)c1"
    reactant_smiles2 = "O=C(O)c1c(F)c(C(=O)O)c(Br)c(C(=O)O)c1Cl"

To handle such cases, we need to enumerate all possible matches, process them one-by-one,
and assign manual mapping based on each match.

This may also require modifications to the reaction selector logic inside
`run_reaction_template_pipeline.py` to support iterating through multiple match candidates.
"""

from rdkit import Chem
from rdkit.Chem import AllChem


def map_reactant_atoms(reactant1, reactant2, rxn, delete_atom=False):
    """
    Maps atom indices in two reactant molecules to the atom map numbers defined
    in the reaction SMARTS template using substructure matching. Runs the reaction
    to produce products, reveals the template map numbers in the products, combines
    reactants and products into single molecules, and optionally extracts map numbers
    from a byproduct for deletion.

    This ensures proper atom tracking across the reaction for downstream analysis,
    such as identifying changed atoms or handling byproducts.

    Args:
        reactant1 (Chem.Mol): First reactant molecule (with hydrogens added).
        reactant2 (Chem.Mol): Second reactant molecule (with hydrogens added).
        rxn (AllChem.ReactionFromSmarts): The chemical reaction object from SMARTS.
        delete_atom (bool): If True, expects a byproduct in the second position of products
            and collects its atom map numbers for deletion handling. Raises error if only
            one product found.

    Returns:
        tuple:
            - combined_reactants (Chem.Mol): Combined molecule of both reactants.
            - combined_products (Chem.Mol): Combined molecule(s) of the first product set.
            - byproduct_map_numbers (list[int]): List of atom map numbers from byproduct
              (empty if delete_atom=False).

    Raises:
        RuntimeError: If no products produced or (if delete_atom=True) only one product found.

    Notes:
        - Initial map numbers: reactant1 atoms get 1001+, reactant2 get 2001+ (for identification).
        - Assumes single product set; uses first set for combination, reveals maps on all.
    """
    # Assign unique initial map numbers to all atoms in reactants for identification
    # reactant1: 1001 + idx, reactant2: 2001 + idx
    for atom in reactant1.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1001)
    for atom in reactant2.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 2001)

    def smart_mapping(reactant, smarts_template, match_tuple):
        """
        Applies atom map numbers from the reaction's reactant SMARTS template
        to the corresponding matched atoms in the reactant molecule.

        Args:
            reactant (Chem.Mol): The reactant molecule.
            smarts_template (Chem.Mol): The reactant part of the reaction SMARTS.
            match_tuple (tuple): Atom indices match from GetSubstructMatch.
        """
        if not match_tuple:
            return

        # Map from SMARTS atom indices to their map numbers (non-zero)
        SMARTS_index_to_Map_Number = {}
        for smarts_atom in smarts_template.GetAtoms():
            map_num = smarts_atom.GetAtomMapNum()
            if map_num != 0:
                SMARTS_index_to_Map_Number[smarts_atom.GetIdx()] = map_num

        # Assign map numbers to matched atoms in reactant
        for smarts_pos, map_num in SMARTS_index_to_Map_Number.items():
            if smarts_pos < len(match_tuple):
                atom_index_in_mol = match_tuple[smarts_pos]
                atom = reactant.GetAtomWithIdx(atom_index_in_mol)
                atom.SetAtomMapNum(map_num)

    # Prepare mapping info for both reactants using substructure matches
    mapping_dict = {
        "reactant1": {
            "reactant": reactant1,
            "smarts_template": rxn.GetReactants()[0],
            "match": reactant1.GetSubstructMatch(rxn.GetReactants()[0])
        },
        "reactant2": {
            "reactant": reactant2,
            "smarts_template": rxn.GetReactants()[1],
            "match": reactant2.GetSubstructMatch(rxn.GetReactants()[1])
        }
    }

    # Apply smart mapping to both reactants
    smart_mapping(
        mapping_dict["reactant1"]["reactant"],
        mapping_dict["reactant1"]["smarts_template"],
        mapping_dict["reactant1"]["match"]
    )
    smart_mapping(
        mapping_dict["reactant2"]["reactant"],
        mapping_dict["reactant2"]["smarts_template"],
        mapping_dict["reactant2"]["match"]
    )

    reactants_tuple = (reactant1, reactant2)

    def reveal_template_map_numbers(mol: Chem.Mol) -> None:
        """
        Reveals the original template atom map numbers in a post-reaction product molecule.
        RDKit sets 'old_mapno' property on atoms after reaction; this copies it to AtomMapNum
        for visualization and tracking.

        Args:
            mol (Chem.Mol): Product molecule from reaction.
        """
        for atom in mol.GetAtoms():
            if atom.HasProp('old_mapno'):  # RDKit "magic" property post-reaction
                map_num = atom.GetIntProp('old_mapno')
                atom.SetAtomMapNum(map_num)

    # Execute the reaction and get all possible product sets
    products_sets = list(rxn.RunReactants(reactants_tuple))
    if not products_sets:
        raise RuntimeError("Reaction produced no products. SMARTS or mapping failed.")

    # Take the first product set as main products
    products = products_sets[0]

    # Reveal map numbers on all products across all sets (for completeness)
    for products in products_sets:
        for product in products:
            reveal_template_map_numbers(product)

    # Handle byproduct if requested
    byproduct_map_numbers = []
    if delete_atom:
        if len(products) == 1:
            raise RuntimeError("Reaction expected to produce a byproduct to delete, but only one product found.")
        byproduct = products[1]  # Assume byproduct is second in tuple
        for atom in byproduct.GetAtoms():
            byproduct_map_numbers.append(atom.GetAtomMapNum())

    # Combine reactants and products into single molecules for easier handling
    combined_reactants = Chem.CombineMols(reactant1, reactant2)
    combined_products = Chem.CombineMols(*products)
    return combined_reactants, combined_products, byproduct_map_numbers

def map_product_atoms(combined_reactants, combined_products, byproduct_map_numbers, delete_atom):
    MAP_dict = {}
    atom_count_reactants = 0
    atom_count_products = 0
    initiator_atom = []
    byproduct_atom = []
    for r_atom in combined_reactants.GetAtoms():
        atom_count_reactants += 1
        for p_atom in combined_products.GetAtoms():
            if r_atom.GetAtomMapNum() == 1 or r_atom.GetAtomMapNum() == 2:
                if r_atom.GetIdx() not in initiator_atom:
                    initiator_atom.append(r_atom.GetIdx())
            if delete_atom and r_atom.GetAtomMapNum() in byproduct_map_numbers:
                if r_atom.GetIdx() not in byproduct_atom:
                    byproduct_atom.append(r_atom.GetIdx())
            if r_atom.GetAtomMapNum() == p_atom.GetAtomMapNum() and r_atom.GetIdx() not in MAP_dict and p_atom.GetIdx() not in MAP_dict.values():
                print(f"Reactant atom {r_atom.GetIdx()} is mapped to product atom {p_atom.GetIdx()}")
                MAP_dict[r_atom.GetIdx()] = p_atom.GetIdx()
                r_atom.SetAtomMapNum(0)
                p_atom.SetAtomMapNum(0)
                atom_count_products += 1
    if atom_count_reactants != atom_count_products:
        raise ValueError(f"Mismatch in Number of mapped atoms between reactants and products. {atom_count_reactnts} vs {atom_count_products}"
                         "  Contact Developers.")

    for molecule in [combined_reactants, combined_products]:
        for atom in molecule.GetAtoms():
            if atom.GetAtomMapNum() != 0:
                raise ValueError(
                    f"Unmapped atom with map number {atom.GetAtomMapNum()} found after reaction processing. "
                    "Please contact developers."
                )
    return MAP_dict, initiator_atom, byproduct_atom


if __name__ == "__main__":
    # Example usage: Esterification- reaction 
    reactant_smiles1 = "C1=CC=C(C(=C1)C(=O)O)O" 
    reactant_smiles2 = "OCCC(O)=O" 
    reaction_smarts = "[O;!$(OC=*):1]-[H:3].[CX3:2](=[O:5])[OX2H1:4]>>[OX2:1]-[CX3:2](=[O:5]).[O:4]-[H:3]"
    
    # Load and add hydrogens to reactants
    reactant1 = Chem.MolFromSmiles(reactant_smiles1)
    reactant1 = Chem.AddHs(reactant1)
    reactant2 = Chem.MolFromSmiles(reactant_smiles2)
    reactant2 = Chem.AddHs(reactant2)
    
    # Create reaction from SMARTS
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)
    delete_atom = True  # Expect byproduct (water)
    
    # Run mapping and reaction
    combined_reactants, combined_products, byproduct_map_numbers = map_reactant_atoms(
        reactant1, reactant2, rxn, delete_atom
    )
    
    # Output results
    print("Combined Reactants SMILES:", Chem.MolToSmiles(combined_reactants))
    print("Combined Products SMILES:", Chem.MolToSmiles(combined_products))
    print("Byproduct Map Numbers:", byproduct_map_numbers)

    MAP_dict, initiator_atom, byproduct_atom = map_product_atoms(
        combined_reactants, combined_products, byproduct_map_numbers, delete_atom
    )
    print("Atom Mapping Dictionary:", MAP_dict)
    print("Initiator Atoms:", initiator_atom)
    print("Byproduct Atoms:", byproduct_atom)
