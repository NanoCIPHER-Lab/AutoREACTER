"""Functional-group definitions organized by the carboxyl and carbonyl groups motif.

Each entry is defined once in the entire functional-group registry. Reaction libraries reference the entry's group_name value."""

FUNCTIONAL_GROUPS = {
    'di_carboxylic_acid_monomer': {
        'functionality_type': 'di_identical',
        'smarts_1': '[CX3:1](=[O])[OX2H1]',
        'group_name': 'di_carboxylic_acid',
        'comments': None
        },
    'di_carboxylic_acid_halide_monomer': {
        'functionality_type': 'di_identical',
        'smarts_1': '[CX3:1](=[O])[Cl,Br,I]',
        'group_name': 'di_carboxylic_acid_halide',
        'comments': None
    },
    'di_carboxylic_ester_monomer': {
        'functionality_type': 'di_identical',
        'smarts_1': '[CX3:1](=[O])[OX2H0][#6]',
        'group_name': 'di_carboxylic_ester',
        'comments': None
    },
    'phosgene_monomer': {
        'functionality_type': 'di_identical',
        'smarts_1': '[CX3](=[OX1])[Cl]',
        'group_name': 'phosgene',
        'comments': None,
    },
    'diphenyl_carbonate_monomer': {
        'functionality_type': 'di_identical',
        'smarts_1': '[CX3](=[OX1])[OX2][c]',
        'group_name': 'diphenyl_carbonate',
        'comments': None
    }
#  'formaldehyde_monomer': {'functionality_type': 'mono',
#                           'smarts_1': '[CH2]=[OX1]',
#                           'group_name': 'formaldehyde',
#                           'comments': None}}
}