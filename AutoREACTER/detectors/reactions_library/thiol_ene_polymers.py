REACTIONS = {
    'Dithiol and Diene Thiol-Ene Click Polymerization': {
        'same_reactants': False,
        'reactant_1': 'dithiol',
        'reactant_2': 'diene',
        'product': 'poly_thioether_chain',
        'delete_atom': False,

        # Thiol-ene addition:
        #
        #     S-H + CH2=C
        #          ↓
        #     S-CH2-C-H
        #
        # AutoREACTER/LAMMPS convention:
        #     maps 1 and 2 are the initiator atoms
        #     the new bond is 1-2
        #
        # Important:
        #     Product map 2 is written as [C:2], not [CH2:2].
        #     When reactants contain explicit hydrogens (Chem.AddHs),
        #     [CH2:2] would add an additional hydrogen specification
        #     while RDKit also preserves the mapped atom's existing H
        #     neighbors, producing an invalid carbon valence.
        #
        # Map 1 = thiol sulfur
        # Map 2 = terminal alkene carbon
        # Map 3 = substituted alkene carbon
        # Map 5 = transferred thiol hydrogen
        'reaction': (
            '[SX2H1:1]-[H:5].'
            '[CH2:2]=[C;!R:3]'
            '>>'
            '[SX2:1]-[C:2]-[C:3]-[H:5]'
        ),

        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },

        'comments': None
    }
}