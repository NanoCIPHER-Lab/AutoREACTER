FUNCTIONAL_GROUPS = {
    'primary_amine_monomer': {
        'functionality_type': 'mono',
        'smarts_1': '[NX3H2;!$(NC=O);!$(NC=[N,O,S])]',
        'group_name': 'primary_amine',
        'comments': None
    },
    'secondary_amine_monomer': {
        'functionality_type': 'mono',
        'smarts_1': '[NX3H1;!$(NC=O);!$(NC=[N,O,S])]',
        'group_name': 'secondary_amine',
        'comments': None
    },
    'di_amine_monomer': {
        'functionality_type': 'di_identical',
        'smarts_1': '[NX3;H2,H1;!$([N][C,S]=*):1]',
        'group_name': 'di_amine',
        'comments': None
    },
    'di_primary_amine_monomer': {
        'functionality_type': 'di_identical',
        'smarts_1': '[NX3H2;!$([N][C,S]=*)]',
        'group_name': 'di_primary_amine',
        'comments': None
    },
    # 'tetra_amine_monomer': {
    #     'functionality_type': 'di_identical',
    #     'smarts_1': '[c]([NX3H2;!$([N][C,S]=*)])[c][NX3H2;!$([N][C,S]=*)]',
    #     'group_name': 'tetra_amine',
    #     'comments': None
    # }
}