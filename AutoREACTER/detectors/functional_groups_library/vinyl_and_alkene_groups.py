FUNCTIONAL_GROUPS = {
    'vinyl_monomer': {
        'functionality_type': 'vinyl',
        'smarts_1': '[CH2]=[C;!R]',
        'group_name': 'vinyl',
        'comments': None
    },
    'diene_monomer': {
        'functionality_type': 'di_identical',
        'smarts_1': '[CH2]=[C;!R]',
        'group_name': 'diene',
        'comments': None
    },
    # 'bis_alkene_monomer': {
    #     'functionality_type': 'di_identical',
    #     'smarts_1': '[C]=[C]',
    #     'group_name': 'bis_alkene',
    #     'comments': None
    # }
    'tetrafluoroethylene_monomer': {
        'functionality_type': 'vinyl',
        'smarts_1': '[CX3](-[F])(-[F])=[CX3](-[F])(-[F])',
        'group_name': 'tetrafluoroethylene',
        'comments': None
    }
}