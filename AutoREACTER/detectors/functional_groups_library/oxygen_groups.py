FUNCTIONAL_GROUPS = {
    'diol_monomer': {
        'functionality_type': 'di_identical',
        'smarts_1': '[OX2H1;!$([O][C,S]=*):1]',
        'group_name': 'diol',
        'comments': None
    },
    # 'initiator_monomer': {
    #     'functionality_type': 'mono',
    #     'smarts_1': '[OX2H1;!$([O][C,S]=*)]',
    #     'group_name': 'lactone_initiator',
    #     'comments': None
    # },
    'water_monomer': {
        'functionality_type': 'mono',
        'smarts_1': '[OH2:1]',
        'group_name': 'water',
        'comments': None
    }
}