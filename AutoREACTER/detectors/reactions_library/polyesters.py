REACTIONS = {
    'Hydroxy Carboxylic Acid Polycondensation(Polyesterification)': 
        {
            'same_reactants': True,
            'reactant_1': 'hydroxy_carboxylic_acid',
            'product': 'polyester_chain',
            'delete_atom': True,
            'reaction': '[OX2;!$([O][C,S]=*):1]-[H:3].[CX3:2](=[O:5])[OX2:4]>>[OX2:1]-[CX3:2](=[O:5]).[O:4]-[H:3]',
            'reference': {
                'smarts': 'https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329',
                'reaction_and_mechanism': [
                    'https://pubs.acs.org/doi/10.1021/ed048pA734.1',
                    'https://pubs.acs.org/doi/10.1021/ed073pA312'
                ]
            },
            'comments': 'Fixed SMARTS to prevent RDKit explicit node valence crashes.'
        },
    'Hydroxy Carboxylic and Hydroxy Carboxylic Polycondensation(Polyesterification)': 
        {
            'same_reactants': False,
            'reactant_1': 'hydroxy_carboxylic_acid',
            'reactant_2': 'hydroxy_carboxylic_acid',
            'product': 'polyester_chain',
            'delete_atom': True,
            'reaction': '[OX2;!$([O][C,S]=*):1]-[H:3].[CX3:2](=[O:5])[OX2:4]>>[OX2:1]-[CX3:2](=[O:5]).[O:4]-[H:3]',
            'reference': {
                'smarts': 'https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329',
                'reaction_and_mechanism': [
                    'https://pubs.acs.org/doi/10.1021/ed048pA734.1',
                    'https://pubs.acs.org/doi/10.1021/ed073pA312'
                ]
            },
            'comments': 'Fixed SMARTS to prevent RDKit explicit node valence crashes.'
        },
    'Hydroxy Acid Halides Polycondensation(Polyesterification)': 
        {
            'same_reactants': True,
            'reactant_1': 'hydroxy_acid_halide',
            'product': 'polyester_chain',
            'delete_atom': True,
            'reaction': '[OX2;!$([O][C,S]=*):1]-[H:3].[CX3:2](=[O:5])[Cl,Br,I:4]>>[OX2:1]-[CX3:2](=[O:5]).[Cl,Br,I:4]-[H:3]',
            'reference': {
                'smarts': None,
                'reaction_and_mechanism': ['https://pubs.acs.org/doi/10.1021/ed073pA312']
            },
            'comments': 'Fixed SMARTS to prevent RDKit explicit node valence crashes.'
        },
    'Hydroxy Acid Halides Hydroxy Acid Halides Polycondensation(Polyesterification)': 
        {
            'same_reactants': False,
            'reactant_1': 'hydroxy_acid_halide',
            'reactant_2': 'hydroxy_acid_halide',
            'product': 'polyester_chain',
            'delete_atom': True,
            'reaction': '[OX2;!$([O][C,S]=*):1]-[H:3].[CX3:2](=[O:5])[Cl,Br,I:4]>>[OX2:1]-[CX3:2](=[O:5]).[Cl,Br,I:4]-[H:3]',
            'reference': {
                'smarts': None,
                'reaction_and_mechanism': ['https://pubs.acs.org/doi/10.1021/ed073pA312']
            },
            'comments': 'Fixed SMARTS to prevent RDKit explicit node valence crashes.'
        },
    'Diol and Di-Carboxylic Acid Polycondensation(Polyesterification)': 
        {
            'same_reactants': False,
            'reactant_1': 'diol',
            'reactant_2': 'di_carboxylic_acid',
            'product': 'polyester_chain',
            'delete_atom': True,
            'reaction': '[CX3:1](=[O:3])[OX2:4].[OX2;!$([O][C,S]=*):2]-[H:5]>>[CX3:1](=[O:3])-[OX2:2].[O:4]-[H:5]',
            'reference': {
                'smarts': ['https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329'],
                'reaction_and_mechanism': [
                    'https://pubs.acs.org/doi/10.1021/ed048pA734.1',
                    'https://pubs.acs.org/doi/10.1021/ed073pA312'
                ]
            },
            'comments': 'Fixed SMARTS to prevent RDKit explicit node valence crashes.'
        },
    'Diol and Di-Acid Halide Polycondensation(Polyesterification)': 
        {
            'same_reactants': False,
            'reactant_1': 'diol',
            'reactant_2': 'di_carboxylic_acid_halide',
            'product': 'polyester_chain',
            'delete_atom': True,
            'reaction': '[CX3:1](=[O:3])[Cl,Br,I:4].[OX2;!$([O][C,S]=*):2]-[H:5]>>[CX3:1](=[O:3])-[OX2:2].[Cl,Br,I:4]-[H:5]',
            'reference': {
                'smarts': ['https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329'],
                'reaction_and_mechanism': [
                    'https://pubs.acs.org/doi/10.1021/ed048pA734.1',
                    'https://pubs.acs.org/doi/10.1021/ed073pA312'
                ]
            },
            'comments': 'Fixed SMARTS to prevent RDKit explicit node valence crashes.'
        },
    'Diol and Di-Carboxylic Ester Polycondensation(Transesterification)': 
        {
            'same_reactants': False,
            'reactant_1': 'diol',
            'reactant_2': 'di_carboxylic_ester',
            'product': 'polyester_chain',
            'delete_atom': True,
            'reaction': '[OX2;!$([O][C,S]=*):1]-[H:3].[CX3:2](=[O:5])[OX2:4][#6:6]>>[OX2:1]-[CX3:2](=[O:5]).[OX2:4](-[H:3])-[#6:6]',
            'reference': {
                'smarts': None,
                'reaction_and_mechanism': None
            },
            'comments': None
        },
    # Skipped for later additinos after proper validations
    # 'Lactone Ring-Opening Polyesterification': 
    #     {
    #         'same_reactants': False,
    #         'reactant_1': 'lactone',
    #         'reactant_2': 'lactone_initiator',
    #         'product': 'polyester_chain',
    #         'delete_atom': False,
    #         'reaction': '[OX2:1]-[H:2].[CX3:3]1(=[OX1:4])-[CX4:6]-[CX4:7]-[CX4:8]-[CX4:9]-[CX4:10]-[OX2:5]1>>[OX2:1]-[CX3:3](=[OX1:4])-[CX4:6]-[CX4:7]-[CX4:8]-[CX4:9]-[CX4:10]-[OX2:5]-[H:2]',
    #         'reference': {
    #             'smarts': None, 
    #             'reaction_and_mechanism': None
    #         },
    #         'comments': None
    #     },
    # 'Cyclic Anhydride and Epoxide Polyesterification': 
    #     {
    #         'same_reactants': False,
    #         'reactant_1': 'cyclic_anhydride',
    #         'reactant_2': 'epoxide',
    #         'product': 'polyester_chain',
    #         'delete_atom': False,
    #         'reaction': '[CX3:1]1(=[OX1:2])-[OX2:3]-[CX3:4](=[OX1:5])-[CX4:9]-[CX4:10]1.[CX4:6]2-[OX2:7]-[CX4:8]2>>[CX3:1](=[OX1:2])-[OX2:7]-[CX4:6]-[CX4:8]-[OX2:3]-[CX3:4](=[OX1:5])-[CX4:9]-[CX4:10]',
    #         'reference': {
    #             'smarts': None, 
    #             'reaction_and_mechanism': None
    #         },
    #         'comments': None
    #     }
}