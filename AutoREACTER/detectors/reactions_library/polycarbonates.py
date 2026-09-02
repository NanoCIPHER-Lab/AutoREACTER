REACTIONS = {
    'Diol and Phosgene Polycondensation(Polycarbonate Formation)': 
        {
            'same_reactants': False,
            'reactant_1': 'diol',
            'reactant_2': 'phosgene',
            'product': 'polycarbonate_chain',
            'delete_atom': True,
            'reaction': '[OX2:1]-[H:4].[CX3:2](=[OX1:5])[Cl:3]>>[OX2:1]-[CX3:2](=[OX1:5]).[Cl:3]-[H:4]',
            'reference': 
                {
                    'smarts': None,
                    'reaction_and_mechanism': None
                },
            'comments': None,
        },
    'Diol and Diphenyl Carbonate Polycondensation(Transcarbonation)': 
        {
            'same_reactants': False,
            'reactant_1': 'diol',
            'reactant_2': 'diphenyl_carbonate',
            'product': 'polycarbonate_chain',
            'delete_atom': True,
            'reaction': '[OX2:1]-[H:4].[CX3:2](=[OX1:5])[OX2:3][c:6]>>[OX2:1]-[CX3:2](=[OX1:5]).[OX2:3](-[H:4])-[c:6]',
            'reference': 
                {
                    'smarts': None,
                    'reaction_and_mechanism': None
                },
            'comments': None
        }
    }