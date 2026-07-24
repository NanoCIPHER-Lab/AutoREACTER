REACTIONS = {
    'Dichlorosilane Hydrolysis to Silanol': {
        'same_reactants': False,
        'reactant_1': 'dichlorosilane',
        'reactant_2': 'water',
        'product': 'silanediol',
        'delete_atom': True,
        'reaction': '[Si:1]-[Cl:2].[OX2H2:3](-[H:4])-[H:5]>>[Si:1]-[OX2:3]-[H:4].[Cl:2]-[H:5]',
        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },
        'comments': None
    },
    'Silanediol Polycondensation(Polysiloxane Formation)': {
        'same_reactants': True,
        'reactant_1': 'silanediol',
        'product': 'polysiloxane_chain',
        'delete_atom': True,
        'reaction': '[Si:1]-[OX2H1:2]-[H:5].[Si:3]-[OX2H1:4]-[H:6]>>[Si:1]-[OX2:2]-[Si:3].[OX2:4](-[H:5])-[H:6]',
        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },
        'comments': None
    },
    'Silanediol and Silanediol Copolycondensation(Polysiloxane Formation)': {
        'same_reactants': False,
        'reactant_1': 'silanediol',
        'reactant_2': 'silanediol',
        'product': 'polysiloxane_chain',
        'delete_atom': True,
        'reaction': '[Si:1]-[OX2H1:2]-[H:5].[Si:3]-[OX2H1:4]-[H:6]>>[Si:1]-[OX2:2]-[Si:3].[OX2:4](-[H:5])-[H:6]',
        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },
        'comments': None
    }
}