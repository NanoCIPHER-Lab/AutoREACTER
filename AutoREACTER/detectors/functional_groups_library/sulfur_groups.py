"""Functional-group definitions organized by the sulfur groups motif.

Each entry is defined once in the entire functional-group registry. Reaction libraries reference the entry's group_name value."""

FUNCTIONAL_GROUPS = {
    'dithiol_monomer': 
        {
            'functionality_type': 'di_identical',
            'smarts_1': '[SX2H1;!$([S][C,S]=*):1]',
            'group_name': 'dithiol',
            'comments': None
        },
#  'sodium_sulfide_monomer': {'functionality_type': 'mono',
#                             'smarts_1': '[S-2]',
#                             'group_name': 'sodium_sulfide',
#                             'comments': 'Sulfide dianion; sodium counterions are included in the reaction SMARTS.'}
}
