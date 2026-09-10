REACTIONS = {
    'Carboxylic Acid and Acid Halide Polycondensation (Polyanhydride Formation)': {
        'same_reactants': True,
        'reactant_1': 'carboxylic_acid_acid_halide',
        'product': 'polyanhydride_chain',
        'delete_atom': True,
        'reaction': '[CX3:1](=[O:3])[Cl,Br,I:4].[CX3:6](=[O:5])[OX2:2]-[H:7]>>[CX3:1](=[O:3])-[OX2:2]-[CX3:6](=[O:5]).[Cl,Br,I:4]-[H:7]',
        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },
        'notes': 'Atom maps 1 and 2 are reserved as AutoREACTER/LAMMPS bond/react initiator atoms.'
    },

    'Carboxylic Acid and Acid Halide Copolycondensation (Polyanhydride Copolymerization)': {
        'same_reactants': False,
        'reactant_1': 'carboxylic_acid_acid_halide',
        'reactant_2': 'carboxylic_acid_acid_halide',
        'product': 'polyanhydride_chain',
        'delete_atom': True,
        'reaction': '[CX3:1](=[O:3])[Cl,Br,I:4].[CX3:6](=[O:5])[OX2:2]-[H:7]>>[CX3:1](=[O:3])-[OX2:2]-[CX3:6](=[O:5]).[Cl,Br,I:4]-[H:7]',
        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },
        'notes': 'Atom maps 1 and 2 are reserved as AutoREACTER/LAMMPS bond/react initiator atoms.'
    }
}