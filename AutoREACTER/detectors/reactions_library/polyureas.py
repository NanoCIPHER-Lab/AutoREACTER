REACTIONS = {
    'Di-Amine and Di-Isocyanate Polyaddition(Polyurea Formation)': {
        'same_reactants': False,
        'reactant_1': 'di_amine',
        'reactant_2': 'di_isocyanate',
        'product': 'polyurea_chain',
        'delete_atom': False,
        'reaction': '[NX3;H2,H1;!$([N][C,S]=*):1]-[H:3].[NX2:4]=[CX2:2]=[OX1:5]>>[NX3:1]-[CX3:2](=[OX1:5])-[NX3:4]-[H:3]',
        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },
        'comments': None
    }
}