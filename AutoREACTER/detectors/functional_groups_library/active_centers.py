"""Functional-group definitions organized by the active centers motif.

Each entry is defined once in the entire functional-group registry. Reaction libraries reference the entry's group_name value."""

FUNCTIONAL_GROUPS = {
    'vinyl_chain_end_radical': 
        {
            'functionality_type': 'vinyl',
            'smarts_1': '[C;!R;D3;v3]',
            'group_name': 'vinyl_chain_end_radical',
            'comments': None
        },
    # 'romp_alkylidene_motif': 
    #     {
    #         'functionality_type': 'mono',
    #         'smarts_1': '[Ru]=[C]',
    #         'group_name': 'romp_alkylidene',
    #         'comments': 'UNTESTED simplified ruthenium alkylidene motif used for both ROMP initiation '
    #                     'and propagation. Atom maps 1 and 2 are assigned only in the reaction SMARTS.'
    #     }
    }
