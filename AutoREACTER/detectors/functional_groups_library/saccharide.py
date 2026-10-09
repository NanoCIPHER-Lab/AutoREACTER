
FUNCTIONAL_GROUPS = {

    'saccharide': {
        'functionality_type': 'di_different',

        # Anomeric carbon C1 with a free hydroxyl
        'smarts_1': (
            '[C;R;H1:1]([OX2H1])-[O;R]'
        ),

        # C4 hydroxyl oxygen, identified using neighboring C5-CH2OH
        'smarts_2': (
            '[OX2H1:2]-[C;R;H1]-'
            '[C;R;H1](-[CH2]-[OX2H1])'
        ),

        'group_name': 'saccharide',
        'comments': None
    },

}
