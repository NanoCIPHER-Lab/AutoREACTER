REACTIONS = {
    'Carboxylic Acid and Acid Halide Polycondensation (Polyanhydride Formation)': {
        'same_reactants': True,
        'reactant_1': 'carboxylic_acid_acid_halide',
        'product': 'polyanhydride_chain',
        'delete_atom': True,
        'reaction': '[CX3:1](=[O:3])[Cl,Br,I:4].[CX3:2](=[O:5])[OX2:6]-[H:7]>>[CX3:1](=[O:3])-[OX2:6]-[CX3:2](=[O:5]).[Cl,Br,I:4]-[H:7]',
        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },
        'comments': None
    },
    'Carboxylic Acid and Acid Halide Copolycondensation (Polyanhydride Copolymerization)': {
        'same_reactants': False,
        'reactant_1': 'carboxylic_acid_acid_halide',
        'reactant_2': 'carboxylic_acid_acid_halide',
        'product': 'polyanhydride_chain',
        'delete_atom': True,
        'reaction': '[CX3:1](=[O:3])[Cl,Br,I:4].[CX3:2](=[O:5])[OX2:6]-[H:7]>>[CX3:1](=[O:3])-[OX2:6]-[CX3:2](=[O:5]).[Cl,Br,I:4]-[H:7]',
        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },
        'comments': None
    }
}