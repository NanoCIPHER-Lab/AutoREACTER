# REACTIONS = {
#     'Phenol and Formaldehyde Hydroxymethylation': {
#         'same_reactants': False,
#         'reactant_1': 'phenol',
#         'reactant_2': 'formaldehyde',
#         'product': 'hydroxymethyl_phenol',
#         'delete_atom': False,
#         'reaction': '[c:1]-[H:4].[CH2:2]=[OX1:3]>>[c:1]-[C:2]-[OX2:3]-[H:4]',
#         'reference': {
#             'smarts': None,
#             'reaction_and_mechanism': None
#         },
#         'comments': None
#     },
#     'Hydroxymethyl Phenol and Phenol Condensation (Methylene Bridge Formation)': {
#         'same_reactants': False,
#         'reactant_1': 'hydroxymethyl_phenol',
#         'reactant_2': 'phenol',
#         'product': 'phenol_formaldehyde_chain',
#         'delete_atom': True,
#         'reaction': '[c:1]-[CH2:2]-[OX2H1:3]-[H:6].[c:4]-[H:5]>>[c:1]-[C:2]-[c:4].[OX2:3](-[H:5])-[H:6]',
#         'reference': {
#             'smarts': None,
#             'reaction_and_mechanism': None
#         },
#         'comments': None
#     }
# }