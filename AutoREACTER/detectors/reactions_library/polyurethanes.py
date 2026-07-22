REACTIONS = {
    'Diol and Di-Isocyanate Polyaddition(Polyurethane Formation)': {
        'same_reactants': False,
        'reactant_1': 'diol',
        'reactant_2': 'di_isocyanate',
        'product': 'polyurethane_chain',
        'delete_atom': False,
        'reaction': '[OX2H1;!$([O][C,S]=*):1]-[H:3].[NX2:4]=[CX2:2]=[OX1:5]>>[OX2:1]-[CX3:2](=[OX1:5])-[NX3:4]-[H:3]',
        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },
        'comments': None
    },
    
    'Dithiol and Di-Isocyanate Polyaddition(Polythiourethane Formation)': {
        'same_reactants': False,
        'reactant_1': 'dithiol',
        'reactant_2': 'di_isocyanate',
        'product': 'polythiourethane_chain',
        'delete_atom': False,
        'reaction': '[SX2H1;!$([S][C,S]=*):1]-[H:3].[NX2:4]=[CX2:2]=[OX1:5]>>[SX2:1]-[CX3:2](=[OX1:5])-[NX3:4]-[H:3]',
        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },
        'comments': 'UNTESTED: chemically consistent replacement for the commented entry whose label said epoxide/isocyanate but whose SMARTS was isocyanate + O/S-H.'
    }
}