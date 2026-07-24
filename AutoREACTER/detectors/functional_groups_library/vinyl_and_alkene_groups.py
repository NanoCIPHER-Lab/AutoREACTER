"""Functional-group definitions organized by the vinyl and alkene groups motif.

Each entry is defined once in the entire functional-group registry. Reaction libraries reference the entry's group_name value."""

FUNCTIONAL_GROUPS = {
    'vinyl_monomer': 
        {
            'functionality_type': 'vinyl',
            'smarts_1': '[CH2]=[C;!R]',
            'group_name': 'vinyl',
            'comments': None
        },
    'diene_monomer': 
        {
            'functionality_type': 'di_identical',
            'smarts_1': '[CH2]=[C;!R]',
            'group_name': 'diene',
            'comments': None
        },
}
#  'bis_alkene_monomer': {'functionality_type': 'di_identical',
#                         'smarts_1': '[C]=[C]',
#                         'group_name': 'bis_alkene',
#                         'comments': 'One alkene site; a bis-alkene supplies two identical matches.'}}
