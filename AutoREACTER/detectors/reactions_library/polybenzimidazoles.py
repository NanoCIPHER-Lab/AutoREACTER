# """
# This is an EDGE case for polybenzimidazole formation reactions. This needs to be futher studied and realease the reaction only after thorough validation.
# """

# REACTIONS = {'Tetra-Amine and Di-Carboxylic Acid Polycondensation (PBI Formation)': 
#                 {
#                     'same_reactants': False,
#                     'reactant_1': 'tetra_amine',
#                     'reactant_2': 'di_carboxylic_acid',
#                     'product': 'polybenzimidazole_chain',
#                     'delete_atom': True,
#                     'reaction': '[c:7]([NX3H2:1](-[H:6])-[H:9])-[c:8]([NX3H2:2](-[H:10])-[H:11]).[CX3:3](=[OX1:4])[OX2H1:5]-[H:12]>>[c:7]1-[NX3:1](-[H:6])-[CX3:3]=[NX2:2]-[c:8]-1.[OX2:4](-[H:9])-[H:10].[OX2:5](-[H:11])-[H:12]',
#                     'reference': {'smarts': None,
#                                 'reaction_and_mechanism': None},
#                     'comments': 'UNTESTED new chemistry; one benzimidazole ring-forming event with two mapped waters.'
#                 }
#             }
