import json
from typing import Dict, Any

class SecondaryReactionLibrary:
    def __init__(self):
        self.secondary_reactions = {
            "Di-Epoxide and Secondary Amine Polyamination": {
                "same_reactants": False,
                "reactant_1": "di_epoxide",
                "reactant_2": "secondary_amine",
                "product": "polyamine_chain",
                "delete_atom": False,
                "reaction": "[C:4]1[O:3][C:2]1.[N:1]-[H:5]>>[H:5]-[O:3]-[C:4]-[C:2]-[N:1]",
                "reference": {},
            }
        }