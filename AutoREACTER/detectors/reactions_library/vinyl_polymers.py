REACTIONS = {

    # =========================================================================
    # Vinyl Addition Polymerization
    # =========================================================================

    'Vinyl Addition Polymerization Initiation': {
        'same_reactants': True,
        'reactant_1': 'vinyl',
        'product': 'vinyl_chain_end_radical',
        'delete_atom': False,

        # H-T self-initiation cheat:
        #
        #     tail=head + tail=head
        #          ↓
        #     cap-tail-head-tail-head*
        #
        # AutoREACTER/LAMMPS rule:
        #     maps 1 and 2 are the initiator atoms.
        #     the new bond is 1-2.
        #
        # Map 3 = first vinyl terminal CH2 tail; capped as CH3
        # Map 1 = first vinyl substituted head
        # Map 2 = second vinyl terminal CH2 tail
        # Map 4 = second vinyl substituted head; becomes active radical
        'reaction': (
            '[CH2:3]=[C;!R:1].'
            '[CH2:2]=[C;!R:4]'
            '>>'
            '[CH3:3]-[C:1]-[CH2:2]-[C;!R:4]'
        ),

        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },

        'comments': None,

        'notes': (
            'Self-initiation cheat for styrene/vinyl. H-T topology. '
            'New bond is maps 1-2. Map 3 is capped as CH3 so the '
            'post-initiation template matches propagation. Map 4 is the '
            'new active chain-end radical.'
        ),
    },


    'Vinyl Addition Polymerization Propagation': {
        'same_reactants': False,
        'reactant_1': 'vinyl',
        'reactant_2': 'vinyl_chain_end_radical',
        'product': 'vinyl_chain_end_radical',
        'delete_atom': False,

        # H-T propagation:
        #
        #     chain-head* + tail=head
        #          ↓
        #     chain-head-tail-head*
        #
        # AutoREACTER/LAMMPS rule:
        #     maps 1 and 2 are the initiator atoms.
        #     the new bond is 1-2.
        #
        # Map 1 = existing radical chain end/head.
        #         This is forced by the loop-detected
        #         vinyl_chain_end_radical index.
        #
        # Map 2 = incoming vinyl terminal CH2 tail
        # Map 3 = incoming vinyl substituted head; becomes new radical
        'reaction': (
            '[CH2:2]=[C;!R:3].'
            '[C;!R:1]'
            '>>'
            '[C:1]-[CH2:2]-[C;!R:3]'
        ),

        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },

        'comments': None,

        'notes': (
            'Head-to-tail vinyl propagation. Existing chain-end radical map '
            '1 bonds to incoming vinyl tail map 2. Map 3 becomes the new '
            'active chain-end radical. Map 1 is selected from explicit '
            'radical detection during loop progression.'
        ),
    },


    # =========================================================================
    # Optional Vinyl Radical Coupling Termination
    # =========================================================================
    #
    # Users can uncomment this reaction if they want to include vinyl radical
    # coupling termination in the simulation.
    #
    # It is normally controlled separately from initiation and propagation
    # so that vinyl growth is not terminated too early.
    #
    # 'Vinyl Radical Coupling Termination': {
    #     'same_reactants': True,
    #     'reactant_1': 'vinyl_chain_end_radical',
    #     'product': 'vinyl_terminated_chain',
    #     'delete_atom': False,
    #
    #     # H-H radical coupling termination:
    #     #
    #     #     chain-head* + *head-chain
    #     #          ↓
    #     #     chain-head-head-chain
    #     #
    #     # This is intentionally the only head-to-head vinyl step.
    #     # Initiation and propagation remain head-to-tail.
    #     #
    #     # AutoREACTER/LAMMPS rule:
    #     #     maps 1 and 2 are the initiator atoms.
    #     #     the new termination bond is 1-2.
    #     #
    #     # Map 1 = active radical head carbon on one vinyl chain end.
    #     # Map 2 = active radical head carbon on another vinyl chain end.
    #     #
    #     # The D3/v3/+0 constraints make this target the
    #     # methacrylate-style chain-end radical center instead of
    #     # an ordinary acyclic carbon.
    #     'reaction': (
    #         '[C;!R;D3;v3;+0:1].'
    #         '[C;!R;D3;v3;+0:2]'
    #         '>>'
    #         '[C:1]-[C:2]'
    #     ),
    #
    #     'reference': {
    #         'smarts': None,
    #         'reaction_and_mechanism': None
    #     },
    #
    #     'comments': None,
    #
    #     'notes': (
    #         'Head-to-head radical coupling termination for '
    #         'PMMA/TEGDMA-style vinyl chain ends. Maps 1 and 2 are '
    #         'the two active radical head carbons and form the new '
    #         '1-2 termination bond. This reaction should be controlled '
    #         'separately from initiation and propagation, usually turned '
    #         'on late or pulsed so it does not kill vinyl growth too early.'
    #     ),
    # },


    # =========================================================================
    # Vinyl Copolymerization Initiation
    # =========================================================================

    'Vinyl Copolymerization Initiation': {
        'same_reactants': False,
        'reactant_1': 'vinyl',
        'reactant_2': 'vinyl',
        'product': 'vinyl_chain_end_radical',
        'delete_atom': False,

        # H-T branchable-vinyl copolymerization initiation:
        #
        #     branchable tail=head + vinyl tail=head
        #              ↓
        #     capped-branchable-head-tail-vinyl-head*
        #
        # This intentionally forces the rare branchable vinyl into the seed.
        #
        # AutoREACTER/LAMMPS rule:
        #     maps 1 and 2 are the initiator atoms.
        #     the new bond is 1-2.
        #
        # Map 3 = branchable vinyl terminal CH2 tail; capped as CH3
        # Map 1 = branchable vinyl substituted head
        # Map 2 = incoming normal vinyl terminal CH2 tail
        # Map 4 = incoming normal vinyl substituted head; becomes active radical
        'reaction': (
            '[CH2:3]=[C;!R:1].'
            '[CH2:2]=[C;!R:4]'
            '>>'
            '[CH3:3]-[C:1]-[CH2:2]-[C;!R;D3;v3;+0:4]'
        ),

        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },

        'comments': None,

        'notes': (
            'Head-to-tail vinyl copolymerization initiation. '
            'Maps 1 and 2 form the new bond. Map 3 is capped as CH3 '
            'and map 4 becomes the active chain-end radical.'
        ),
    },


    # =========================================================================
    # Currently Disabled Cyclic / Legacy Vinyl Reactions
    # =========================================================================

    # 'Vinyl Addition Polymerization': {
    #     'same_reactants': True,
    #     'reactant_1': 'vinyl',
    #     'product': 'polyvinyl_chain',
    #     'delete_atom': False,
    #     'reaction': (
    #         '[CH2:1]=[CH;H1,H0;!R:2].'
    #         '[CH2:3]=[CH;H1,H0;!R:4]'
    #         '>>'
    #         '[CH2:1]-[CH:2]-[CH2:3]-[CH:4]'
    #     ),
    #     'reference': {
    #         'smarts': None,
    #         'reaction_and_mechanism': None
    #     },
    #     'comments': None
    # },


    # 'Vinyl Copolymerization': {
    #     'same_reactants': False,
    #     'reactant_1': 'vinyl',
    #     'reactant_2': 'vinyl',
    #     'product': 'copolyvinyl_chain',
    #     'delete_atom': False,
    #     'reaction': (
    #         '[CH2:1]=[CH;H1,H0;!R:2].'
    #         '[CH2:3]=[CH;H1,H0;!R:4]'
    #         '>>'
    #         '[CH2:1]-[CH:2]-[CH2:3]-[CH:4]'
    #     ),
    #     'reference': {
    #         'smarts': None,
    #         'reaction_and_mechanism': None
    #     },
    #     'comments': None
    # },


    # -------------------------------------------------------------------------
    # Cyclic olefin reactions intentionally remain disabled
    # -------------------------------------------------------------------------

    # 'Cyclic Olefin Addition Polymerization': {
    #     'same_reactants': True,
    #     'reactant_1': 'cyclic_olefin',
    #     'product': 'polycyclic_chain',
    #     'delete_atom': False,
    #
    #     'reaction': (
    #         '[CX3;R:1]=[CX3;R:2].'
    #         '[CX3;R:3]=[CX3;R:4]'
    #         '>>'
    #         '[CX4:1]-[CX4:2]-[CX4:3]-[CX4:4]'
    #     ),
    #
    #     'reference': {
    #         'smarts': None,
    #         'reaction_and_mechanism': None
    #     },
    #
    #     'comments': None
    # },


    # 'Cyclic Olefin and Vinyl Copolymerization': {
    #     'same_reactants': False,
    #     'reactant_1': 'vinyl',
    #     'reactant_2': 'cyclic_olefin',
    #     'product': 'copolycyclicvinyl_chain',
    #     'delete_atom': False,
    #
    #     'reaction': (
    #         '[CH2:1]=[C;!R:2].'
    #         '[C;R:3]=[C;R:4]'
    #         '>>'
    #         '[C:1]-[C:2]-[C:3]-[C:4]'
    #     ),
    #
    #     'reference': {
    #         'smarts': None,
    #         'reaction_and_mechanism': None
    #     },
    #
    #     'comments': None
    # },


    # 'Cyclic Olefin Copolymerization': {
    #     'same_reactants': False,
    #     'reactant_1': 'cyclic_olefin',
    #     'reactant_2': 'cyclic_olefin',
    #     'product': 'copolycyclic_chain',
    #     'delete_atom': False,
    #
    #     'reaction': (
    #         '[CX3;R:1]=[CX3;R:2].'
    #         '[CX3;R:3]=[CX3;R:4]'
    #         '>>'
    #         '[CX4:1]-[CX4:2]-[CX4:3]-[CX4:4]'
    #     ),
    #
    #     'reference': {
    #         'smarts': None,
    #         'reaction_and_mechanism': None
    #     },
    #
    #     'comments': None
    # },


    # =========================================================================
    # Tetrafluoroethylene / PTFE
    # =========================================================================

    'Tetrafluoroethylene Initiation': {
        'same_reactants': False,
        'reactant_1': 'tetrafluoroethylene',
        'reactant_2': 'tetrafluoroethylene',
        'product': 'vinyl_chain_end_radical',
        'delete_atom': False,

        # AutoREACTER/LAMMPS rule:
        #
        #     map 1 + map 2
        #          ↓
        #       new 1-2 bond
        #
        # Map 1 = existing initiating/active atom
        # Map 2 = one TFE carbon
        # Map 3 = second TFE carbon
        # Maps 4-7 = fluorines
        'reaction': (
            '[*:1].'
            '[CX3:2](-[F:4])(-[F:5])='
            '[CX3:3](-[F:6])(-[F:7])'
            '>>'
            '[*:1]-'
            '[CX4:2](-[F:4])(-[F:5])-'
            '[CX4:3](-[F:6])(-[F:7])'
        ),

        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },

        'comments': (
            'Initiation of TFE to form the first active center'
        ),

        'notes': (
            'The new TFE initiation bond is between atom maps 1 and 2.'
        ),
    },


    'Tetrafluoroethylene Propagation': {
        'same_reactants': False,
        'reactant_1': 'vinyl_chain_end_radical',
        'reactant_2': 'tetrafluoroethylene',
        'product': 'vinyl_chain_end_radical',
        'delete_atom': False,

        # PTFE propagation:
        #
        #     chain-1-2    3=8
        #
        #          ↓
        #
        #     chain-1-2-3-8
        #
        # IMPORTANT:
        #
        # Map 1-2 is already bonded in the reactant.
        #
        # The new propagation bond is:
        #
        #     2-3
        #
        # Therefore this reaction must override the registry default
        # initiator_atom_maps=(1, 2).
        'initiator_atom_maps': (2, 3),

        'reaction': (
            '[*:1]-'
            '[CX4:2](-[F:4])(-[F:5]).'
            '[CX3:3](-[F:6])(-[F:7])='
            '[CX3:8](-[F:9])(-[F:10])'
            '>>'
            '[*:1]-'
            '[CX4:2](-[F:4])(-[F:5])-'
            '[CX4:3](-[F:6])(-[F:7])-'
            '[CX4:8](-[F:9])(-[F:10])'
        ),

        'reference': {
            'smarts': None,
            'reaction_and_mechanism': None
        },

        'comments': (
            'Propagation step for PTFE chain growth'
        ),

        'notes': (
            'The new propagation bond is atom maps 2-3. '
            'Maps 1-2 are already connected in the incoming chain, '
            'so initiator_atom_maps is explicitly set to (2, 3).'
        ),
    },


    # =========================================================================
    # Optional PTFE Termination
    # =========================================================================

    # 'Tetrafluoroethylene Termination (Recombination)': {
    #     'same_reactants': True,
    #     'reactant_1': 'vinyl_chain_end_radical',
    #     'product': 'ptfe_chain',
    #     'delete_atom': False,
    #
    #     # Two growing chains meet and recombine at their active
    #     # tetrafluoroethylene chain centers.
    #
    #     'reaction': (
    #         '[*:1]-[CX4:2](-[F:4])(-[F:5]).'
    #         '[*:3]-[CX4:6](-[F:7])(-[F:8])'
    #         '>>'
    #         '[*:1]-'
    #         '[CX4:2](-[F:4])(-[F:5])-'
    #         '[CX4:6](-[F:7])(-[F:8])-'
    #         '[*:3]'
    #     ),
    #
    #     'reference': {
    #         'smarts': None,
    #         'reaction_and_mechanism': None
    #     },
    #
    #     'comments': (
    #         'UNTESTED: Radical recombination termination for PTFE'
    #     )
    # },
}