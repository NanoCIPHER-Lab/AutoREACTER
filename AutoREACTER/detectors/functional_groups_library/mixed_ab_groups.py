FUNCTIONAL_GROUPS = {
    'hydroxy_carboxylic_acid_monomer': {
        'functionality_type': 'di_different',
        'smarts_1': '[OX2H1;!$([O][C,S]=*):1]',
        'smarts_2': '[CX3:2](=[O])[OX2H1]',
        'group_name': 'hydroxy_carboxylic_acid',
        'comments': None
    },

    'hydroxy_acid_halides_monomer': {
        'functionality_type': 'di_different',
        'smarts_1': '[OX2H1;!$([O][C,S]=*):1]',
        'smarts_2': '[CX3:2](=[O])[Cl,Br,I]',
        'group_name': 'hydroxy_acid_halide',
        'comments': None
    },

    'amino_acid_monomer': {
        'functionality_type': 'di_different',
        'smarts_1': '[NX3;H2,H1;!$([N][C,S]=*):1]',
        'smarts_2': '[CX3:2](=[O])[OX2H1]',
        'group_name': 'amino_acid',
        'comments': None
    },

    'carboxylic_acid_acid_halide_monomer': {
        'functionality_type': 'di_different',
        'smarts_1': '[CX3:1](=[O])[OX2H1]',
        'smarts_2': '[CX3:2](=[O])[Cl,Br,I]',
        'group_name': 'carboxylic_acid_acid_halide',
        'comments': None
    },

    'hydroxy_thiol_monomer': {
        'functionality_type': 'di_different',
        'smarts_1': '[OX2H1;!$([O][C,S]=*):1]',
        'smarts_2': '[SX2H1;!$([S][C,S]=*):2]',
        'group_name': 'hydroxy_thiol',
        'comments': None
    }
}