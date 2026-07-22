
# REACTIONS = {
#     'Epoxide Ring-Opening Polyetherification': 
#         {
#             'same_reactants': False,
#             'reactant_1': 'epoxide',
#             'reactant_2': 'initiator',
#             'product': 'polyether_chain',
#             'delete_atom': False,
#             'reaction': '[OX2H1:1]-[H:2].[CX4:3]1[OX2:4][CX4:5]1>>[OX2:1]-[CX4:3]-[CX4:5]-[OX2:4]-[H:2]',
#             'reference': 
#                 {
#                     'smarts': None, 
#                     'reaction_and_mechanism': None
#                 },
#             'comments': None
#         },
#     'Cyclic Anhydride and Epoxide Polyetherification': 
#         {
#             'same_reactants': False,
#             'reactant_1': 'cyclic_anhydride',
#             'reactant_2': 'epoxide',
#             'product': 'polyester_chain',
#             'delete_atom': False,
#             'reaction': '[CX3;R:1](=[OX1:2])[OX2;R:3][CX3;R:4](=[OX1:5]).[CX4:6]1[OX2:7][CX4:8]1>>[CX3:1](=[OX1:2])-[OX2:7]-[CX4:8]-[CX4:6]-[OX2:3]-[CX3:4](=[OX1:5])',
#             'reference': 
#                 {
#                     'smarts': None, 
#                     'reaction_and_mechanism': None
#                 },
#             'comments': None
#         },
#     'Hindered Phenol Polyetherification': 
#         {
#             'same_reactants': True,
#             'reactant_1': 'hindered_phenol',
#             'product': 'polyether_chain',
#             'delete_atom': False,
#             'reaction': '[c:1]-[OX2H1:2]-[H:5].[cH1:3]-[H:4]>>[c:1]-[OX2:2]-[c:3].[H:4]-[H:5]',
#             'reference': 
#                 {
#                     'smarts': None, 
#                     'reaction_and_mechanism': None
#                 },
#             'comments': None
#         },
#     'Hindered Phenol Hindered Phenol Polyetherification':
#         {
#             'same_reactants': False,
#             'reactant_1': 'hindered_phenol',
#             'reactant_2': 'hindered_phenol',
#             'product': 'polyether_chain',
#             'delete_atom': False,
#             'reaction': '[c:1]-[OX2H1:2]-[H:5].[cH1:3]-[H:4]>>[c:1]-[OX2:2]-[c:3].[H:4]-[H:5]',
#             'reference': 
#                 {
#                     'smarts': None, 
#                     'reaction_and_mechanism': None
#                 },
#             'comments': None
#         },
#     'Bis(p-halogenatedaryl)sulfone Diol (without thiol) polycondensation': 
#         {
#             'same_reactants': False,
#             'reactant_1': 'bis(p-halogenatedaryl)sulfone',
#             'reactant_2': 'diol',
#             'product': 'polyether_chain',
#             'delete_atom': True,
#             'reaction': '[c:1]([F,Cl,Br,I:3]).[OX2H1:2]-[H:4]>>[c:1]-[OX2:2].[F,Cl,Br,I:3]-[H:4]',
#             'reference': 
#                 {
#                     'smarts': None,
#                     'reaction_and_mechanism': None
#                 },
#             'comments': None
#         },
#     'Bis(bis(p-fluoroaryl)ketone Diol (without thiol) polycondensation': 
#         {
#             'same_reactants': False,
#             'reactant_1': 'bis(p-fluoroaryl)ketone_monomer',
#             'reactant_2': 'diol',
#             'product': 'polyether_chain',
#             'delete_atom': True,
#             'reaction': '[c:1]([F:3]).[OX2H1:2]-[H:4]>>[c:1]-[OX2:2].[F:3]-[H:4]',
#             'reference': 
#             {
#                 'smarts': None,
#                 'reaction_and_mechanism': None
#             },
#             'comments': None
#         }
# }