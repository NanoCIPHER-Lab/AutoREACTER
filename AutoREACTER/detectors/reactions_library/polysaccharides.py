
REACTIONS = {

    'Cellulose Glycosidic Bond Formation': {
        'same_reactants': True,
        'reactant_1': 'saccharide',
        'product': 'saccharide',
        'delete_atom': True,

        # Simplified glucose condensation:
        #
        #     Glucose + Glucose
        #             ↓
        #     Cellulose chain + H2O
        #
        # Map 1 = Anomeric carbon C1 (unchanged)
        # Map 2 = Hydrogen attached to acceptor C4
        # Map 3 = Leaving anomeric hydroxyl oxygen
        # Map 4 = Acceptor C4 carbon
        # Map 9 = Hydrogen transferred from C4 hydroxyl
        # Map 10 = Hydrogen on leaving hydroxyl
        # Map 11 = Acceptor C4 hydroxyl oxygen
        #
        # Maps 1 and 2 belong to different reactant molecules.
        # Actual new glycosidic bond is 1-11.

        'reaction': (
            '[C;R;H1:3]([O;H1:1]-[H:10])-[O;R:5].'
            '[O;H1:2](-[H:9])-[C;R;H1:4](-[H:11])-'
            '[C;R;H1:6](-[CH2:7]-[O;H1:8])'
            '>>'
            '[C:3](-[O:2]-[C:4](-[H:11])-'
            '[C:6](-[C:7]-[O:8]))-[O:5].'
            '[O:1](-[H:9])-[H:10]'
        ),


        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },

        'comments': (
            'Simplified cellulose formation through '
            '1,4-glycosidic bond formation between glucose units.'
        )
    },

}
