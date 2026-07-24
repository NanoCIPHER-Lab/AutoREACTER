REACTIONS = {
    'Dithiol and Diene Thiol-Ene Click Polymerization': {
        'same_reactants': False,
        'reactant_1': 'dithiol',
        'reactant_2': 'diene',
        'product': 'poly_thioether_chain',
        'delete_atom': False,
        'reaction': '[SX2H1:1]-[H:5].[CH2:2]=[C;!R:3]>>[SX2:1]-[CH2:2]-[C:3]-[H:5]',
        'reference': {
            'smarts': None, 
            'reaction_and_mechanism': None
        },
        'comments': None
    }
}