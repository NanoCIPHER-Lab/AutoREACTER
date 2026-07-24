REACTIONS = {
    'Dithiol and Di-Carboxylic Acid Halide Polycondensation(Polythioesterification)': {
        'same_reactants': False,
        'reactant_1': 'dithiol',
        'reactant_2': 'di_carboxylic_acid_halide',
        'product': 'polythioester_chain',
        'delete_atom': True,
        'reaction': '[CX3:1](=[O:3])[Cl,Br,I:4].[SX2H1;!$([S][C,S]=*):2]-[H:5]>>[CX3:1](=[O:3])-[SX2:2].[Cl,Br,I:4]-[H:5]',
        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },
        'comments': None
    },
    'Dithiol and Di-Carboxylic Acid Polycondensation(Polythioesterification)': {
        'same_reactants': False,
        'reactant_1': 'dithiol',
        'reactant_2': 'di_carboxylic_acid',
        'product': 'polythioester_chain',
        'delete_atom': True,
        'reaction': '[CX3:1](=[O:3])[OX2H1:4].[SX2H1;!$([S][C,S]=*):2]-[H:5]>>[CX3:1](=[O:3])-[SX2:2].[O:4]-[H:5]',
        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },
        'comments': None
    },
    'Hydroxy-Thiol and Di-Carboxylic Acid Halide Polycondensation through Hydroxy Group': {
        'same_reactants': False,
        'reactant_1': 'hydroxy_thiol',
        'reactant_2': 'di_carboxylic_acid_halide',
        'product': 'mixed_polyester_polythioester_chain',
        'delete_atom': True,
        'reaction': '[CX3:1](=[O:3])[Cl,Br,I:4].[OX2H1;!$([O][C,S]=*):2]-[H:5]>>[CX3:1](=[O:3])-[OX2:2].[Cl,Br,I:4]-[H:5]',
        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },
        'comments': None
    },
    'Hydroxy-Thiol and Di-Carboxylic Acid Halide Polycondensation through Thiol Group': {
        'same_reactants': False,
        'reactant_1': 'hydroxy_thiol',
        'reactant_2': 'di_carboxylic_acid_halide',
        'product': 'mixed_polyester_polythioester_chain',
        'delete_atom': True,
        'reaction': '[CX3:1](=[O:3])[Cl,Br,I:4].[SX2H1;!$([S][C,S]=*):2]-[H:5]>>[CX3:1](=[O:3])-[SX2:2].[Cl,Br,I:4]-[H:5]',
        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },
        'comments': None
    }
}