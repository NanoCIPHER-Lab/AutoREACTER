
REACTIONS = {

    'Dichlorosilane Hydrolysis to Silanol': {
        'same_reactants': False,
        'reactant_1': 'dichlorosilane',
        'reactant_2': 'water',
        'product': 'silanediol',
        'delete_atom': True,

        'reaction': (
            '[Si:1]-[Cl:3].'
            '[O;X2:2](-[H:4])-[H:5]'
            '>>'
            '[Si:1]-[O:2]-[H:4].'
            '[Cl:3]-[H:5]'
        ),

        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },

        'notes': (
            'Atom maps 1 and 2 are reserved as AutoREACTER/LAMMPS '
            'bond/react initiator atoms.'
        )
    },


    'Silanediol Polycondensation(Polysiloxane Formation)': {
        'same_reactants': True,
        'reactant_1': 'silanediol',
        'product': 'polysiloxane_chain',
        'delete_atom': True,

        'reaction': (
            '[Si:3]-[O;X2:1]-[H:5].'
            '[Si:2]-[O;X2:4]-[H:6]'
            '>>'
            '[Si:3]-[O:1]-[Si:2].'
            '[O:4](-[H:5])-[H:6]'
        ),

        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },

        'notes': None
    },


    'Silanediol and Silanediol Copolycondensation(Polysiloxane Formation)': {
        'same_reactants': False,
        'reactant_1': 'silanediol',
        'reactant_2': 'silanediol',
        'product': 'polysiloxane_chain',
        'delete_atom': True,

        'reaction': (
            '[Si:3]-[O;X2:1]-[H:5].'
            '[Si:2]-[O;X2:4]-[H:6]'
            '>>'
            '[Si:3]-[O:1]-[Si:2].'
            '[O:4](-[H:5])-[H:6]'
        ),

        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },

        'notes': None
    }
}
