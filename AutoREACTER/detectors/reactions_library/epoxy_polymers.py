REACTIONS = {
    'Primary Amine and Epoxide Polyaddition (Epoxy-Amine, First Addition)': {   
        'same_reactants': False,
        'reactant_1': 'primary_amine',
        'reactant_2': 'di_epoxide',
        'product': 'secondary_amine_hydroxyl_product',
        'delete_atom': False,
        # FIXED: N attacks the less hindered CH2 (:2), O stays on the more hindered CH (:3)
        'reaction': '[NX3H2:1]-[H:6].[CH2;X4:2]1[OX2:5][CH1;X4:3]1>>[NX3:1][C:2][C:3][OX2:5][H:6]',
        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },
        'comments': None
    },
    'Secondary Amine and Epoxide Polyaddition (Epoxy-Amine, Second Addition / Crosslink)': {
        'same_reactants': False,
        'reactant_1': 'secondary_amine',
        'reactant_2': 'di_epoxide',
        'product': 'tertiary_amine_crosslink_product',
        'delete_atom': False,
        # FIXED: N attacks the less hindered CH2 (:2), O stays on the more hindered CH (:3)
        'reaction': '[NX3H1:1]-[H:6].[CH2;X4:2]1[OX2:5][CH1;X4:3]1>>[NX3:1][C:2][C:3][OX2:5][H:6]',
        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },
        'comments': None
    }
}