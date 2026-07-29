REACTIONS = {
    'Vinyl Addition Polymerization Initiation': {
    'same_reactants': True,
    'reactant_1': 'vinyl',
    'product': 'vinyl_chain_end_radical',
    'delete_atom': False,
    'reaction': '[CH2:1]=[C;!R:3].[CH2:2]=[C;!R:4]>>[CH2:1](-[C:3])-[CH2:2]-[C:4]',
    'reference': {'smarts': None, 'reaction_and_mechanism': None},
    'comments': None
},

'Vinyl Addition Polymerization Propagation': {
    'same_reactants': False,
    'reactant_1': 'vinyl',
    'reactant_2': 'vinyl_chain_end_radical',
    'product': 'vinyl_chain_end_radical',
    'delete_atom': False,
    'reaction': '[CH2:2]=[C;!R:3].[C;!R;D3;v3:1]>>[C:1]-[CH2:2]-[C:3]',
    'reference': {'smarts': None, 'reaction_and_mechanism': None},
    'comments': None
},

'Vinyl Copolymerization': {
    'same_reactants': False,
    'reactant_1': 'vinyl',
    'reactant_2': 'vinyl',
    'product': 'copolyvinyl_chain',
    'delete_atom': False,
    'reaction': '[CH2:1]=[C;!R:2].[CH2:3]=[C;!R:4]>>[CH2:1]-[C:2]-[CH2:3]-[C:4]',
    'reference': {'smarts': None, 'reaction_and_mechanism': None},
    'comments': 'General vinyl copolymerization; supports terminal vinyl and methacrylate-style substituted vinyls'
},
    # 'Vinyl Radical Coupling Termination': {
    #     'same_reactants': True,
    #     'reactant_1': 'vinyl_chain_end_radical',
    #     'product': 'vinyl_terminated_chain',
    #     'delete_atom': False,
    #     'reaction': '[C;!R;D3;v3;+0:1].[C;!R;D3;v3;+0:2]>>[C:1]-[C:2]',
    #     'reference': {'smarts': None, 'reaction_and_mechanism': None},
    #     'comments': None
    # },
    # 'Vinyl Addition Polymerization': {
    #     'same_reactants': True,
    #     'reactant_1': 'vinyl',
    #     'product': 'polyvinyl_chain',
    #     'delete_atom': False,
    #     'reaction': '[CH2:1]=[CH;H1,H0;!R:2].[CH2:3]=[CH;H1,H0;!R:4]>>[CH2:1]-[CH:2]-[CH2:3]-[CH:4]',
    #     'reference': {'smarts': None, 'reaction_and_mechanism': None},
    #     'comments': None
    # },
    # 'Vinyl Copolymerization': {
    #     'same_reactants': False,
    #     'reactant_1': 'vinyl',
    #     'reactant_2': 'vinyl',
    #     'product': 'copolyvinyl_chain',
    #     'delete_atom': False,
    #     'reaction': '[CH2:1]=[CH;H1,H0;!R:2].[CH2:3]=[CH;H1,H0;!R:4]>>[CH2:1]-[CH:2]-[CH2:3]-[CH:4]',
    #     'reference': {'smarts': None, 'reaction_and_mechanism': None},
    #     'comments': None
    # },
    # Later Implementation
    # 'Cyclic Olefin Addition Polymerization': {
    #     'same_reactants': True,
    #     'reactant_1': 'cyclic_olefin',
    #     'product': 'polycyclic_chain',
    #     'delete_atom': False,
    #     'reaction': '[CX3;R:1]=[CX3;R:2].[CX3;R:3]=[CX3;R:4]>>[CX4:1]-[CX4:2]-[CX4:3]-[CX4:4]',
    #     'reference': {'smarts': None, 'reaction_and_mechanism': None},
    #     'comments': None
    # },
    # 'Cyclic Olefin and Vinyl Copolymerization': {
    #     'same_reactants': False,
    #     'reactant_1': 'vinyl',
    #     'reactant_2': 'cyclic_olefin',
    #     'product': 'copolycyclicvinyl_chain',
    #     'delete_atom': False,
    #     'reaction': '[CH2:1]=[C;!R:2].[C;R:3]=[C;R:4]>>[C:1]-[C:2]-[C:3]-[C:4]',
    #     'reference': {'smarts': None, 'reaction_and_mechanism': None},
    #     'comments': None
    # },
    # 'Cyclic Olefin Copolymerization': {
    #     'same_reactants': False,
    #     'reactant_1': 'cyclic_olefin',
    #     'reactant_2': 'cyclic_olefin',
    #     'product': 'copolycyclic_chain',
    #     'delete_atom': False,
    #     'reaction': '[CX3;R:1]=[CX3;R:2].[CX3;R:3]=[CX3;R:4]>>[CX4:1]-[CX4:2]-[CX4:3]-[CX4:4]',
    #     'reference': {'smarts': None, 'reaction_and_mechanism': None},
    #     'comments': None
    # }

    'Tetrafluoroethylene Addition Polymerization Initiation': {
        'same_reactants': True,
        'reactant_1': 'tetrafluoroethylene',
        'product': 'vinyl_chain_end_radical',
        'delete_atom': False,
        'reaction': '[CX3:1](-[F:5])(-[F:6])=[C;!R:3](-[F:7])(-[F:8]).[CX3:2](-[F:9])(-[F:10])=[C;!R:4](-[F:11])(-[F:12])>>[CX3:1](-[F:5])(-[F:6])(-[C:3](-[F:7])(-[F:8]))-[CX3:2](-[F:9])(-[F:10])-[C;!R;D3;v3:4](-[F:11])(-[F:12])',
        'reference': {'smarts': None, 'reaction_and_mechanism': None},
        'comments': 'Self-initiation cheat for TFE polymerization'
    },
    'Tetrafluoroethylene Addition Polymerization Propagation': {
        'same_reactants': False,
        'reactant_1': 'tetrafluoroethylene',
        'reactant_2': 'vinyl_chain_end_radical',
        'product': 'vinyl_chain_end_radical',
        'delete_atom': False,
        'reaction': '[CX3:2](-[F:5])(-[F:6])=[C;!R:3](-[F:7])(-[F:8]).[C;!R;D3;v3:1](-[F:9])(-[F:10])>>[C:1](-[F:9])(-[F:10])-[CX3:2](-[F:5])(-[F:6])-[C;!R;D3;v3:3](-[F:7])(-[F:8])',
        'reference': {'smarts': None, 'reaction_and_mechanism': None},
        'comments': 'Propagation step carrying the radical center forward for TFE'
    }
}