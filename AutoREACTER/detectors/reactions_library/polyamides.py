REACTIONS = {
    'Amino Acid Polycondensation (Polyamidation)': {
        'same_reactants': True,
        'reactant_1': 'amino_acid',
        'product': 'polyamide_chain',
        'delete_atom': True,
        'reaction': '[NX3;H2,H1;!$([N][C,S]=*):1]-[H:3].[CX3:2](=[O:4])[OX2H1:5]>>[NX3:1]-[CX3:2](=[O:4]).[O:5]-[H:3]',
        'reference': {
            'smarts': ['https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329'],
            'reaction_and_mechanism': [
                'https://pubs.acs.org/doi/10.1021/ed048pA734.1',
                'https://pubs.acs.org/doi/10.1021/ed073pA312'
            ]
        },
        'comments': None
    },
    'Amino Acid and Amino Acid Polycondensation (Polyamidation)': {
        'same_reactants': False,
        'reactant_1': 'amino_acid',
        'reactant_2': 'amino_acid',
        'product': 'polyamide_chain',
        'delete_atom': True,
        'reaction': '[NX3;H2,H1;!$([N][C,S]=*):1]-[H:3].[CX3:2](=[O:4])[OX2H1:5]>>[NX3:1]-[CX3:2](=[O:4]).[O:5]-[H:3]',
        'reference': {
            'smarts': ['https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329'],
            'reaction_and_mechanism': [
                'https://pubs.acs.org/doi/10.1021/ed048pA734.1',
                'https://pubs.acs.org/doi/10.1021/ed073pA312'
            ]
        },
        'comments': None
    },
    'Di-Amine and Di-Carboxylic Acid Polycondensation (Polyamidation)': {
        'same_reactants': False,
        'reactant_1': 'di_amine',
        'reactant_2': 'di_carboxylic_acid',
        'product': 'polyamide_chain',
        'delete_atom': True,
        'reaction': '[NX3;H2,H1;!$([N][C,S]=*):1]-[H:3].[CX3:2](=[O:4])[OX2H1:5]>>[NX3:1]-[CX3:2](=[O:4]).[O:5]-[H:3]',
        'reference': {
            'smarts': ['https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329'],
            'reaction_and_mechanism': [
                'https://pubs.acs.org/doi/10.1021/ed048pA734.1'
            ]
        },
        'comments': None
    },
    'Di-Amine and Di-Carboxylic Acid Halide Polycondensation (Polyamidation)': {
        'same_reactants': False,
        'reactant_1': 'di_amine',
        'reactant_2': 'di_carboxylic_acid_halide',
        'product': 'polyamide_chain',
        'delete_atom': True,
        'reaction': '[NX3;H2,H1;!$([N][C,S]=*):1]-[H:3].[CX3:2](=[O:4])[Cl,Br,I:5]>>[NX3:1]-[CX3:2](=[O:4]).[Cl,Br,I:5]-[H:3]',
        'reference': {
            'smarts': ['https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329'],
            'reaction_and_mechanism': [
                'https://pubs.acs.org/doi/10.1021/ed048pA734'
            ]
        },
        'comments': None
    },
    'Hydrolytic Initiation of Caprolactam': {
        'same_reactants': False,
        'reactant_1': 'water',
        'reactant_2': 'lactam',
        'product': 'polyamide_chain',
        'delete_atom': False,
        'reaction': '[O:1]-[H:12].[CX3:2]1(=[OX1:3])-[CX4:7]-[CX4:8]-[CX4:9]-[CX4:10]-[CX4:11]-[NX3:4]1-[H:6]>>[O:1]-[CX3:2](=[OX1:3])-[CX4:7]-[CX4:8]-[CX4:9]-[CX4:10]-[CX4:11]-[NX3:4](-[H:6])-[H:12]',
        'reference': {
            'smarts': None,
            'reaction_and_mechanism': ['https://doi.org/10.1002/047147875X']
        },
        'comments': None,
        'Notes': 'Validated on 2026-07-26, Passed'
    },
    # 'Caprolactam Ring-Opening Polyamidation': {
    #     'same_reactants': True,
    #     'reactant_1': 'lactam',
    #     'product': 'polyamide_chain',
    #     'delete_atom': False,
    #     # Reactant 2 explicitly maps the 5 CH2 groups (maps 7 through 11) and uses '1' for the ring closure.
    #     # Product SMARTS removes the '.' and explicitly connects the opened chain.
    #     'reaction': '[NX3:1]-[H:5].[CX3:2]1(=[OX1:3])-[CX4:7]-[CX4:8]-[CX4:9]-[CX4:10]-[CX4:11]-[NX3:4]1-[H:6]>>[NX3:1]-[CX3:2](=[OX1:3])-[CX4:7]-[CX4:8]-[CX4:9]-[CX4:10]-[CX4:11]-[NX3:4](-[H:5])-[H:6]',
    #     'reference': {
    #         'smarts': None,
    #         'reaction_and_mechanism': None
    #     },
    #     'comments': None
    # },
    # 'Lactam and Lactam Ring-Opening Copolyamidation': {
    #     'same_reactants': False,
    #     'reactant_1': 'lactam',
    #     'reactant_2': 'lactam',
    #     'product': 'polyamide_chain',
    #     'delete_atom': False,
    #     'reaction': '[NX3H1;R:1]-[H:2].[CX3;R:3](=[OX1:4])[NX3H1;R:5]-[H:6]>>[NX3:1]-[CX3:3](=[OX1:4]).[NX3:5](-[H:2])-[H:6]',
    #     'reference': {
    #         'smarts': None,
    #         'reaction_and_mechanism': None
    #     },
    #     'comments': None
    # }
}