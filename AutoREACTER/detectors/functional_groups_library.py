"""
This module defines a library of functional groups relevant to polymer chemistry, 
particularly for the detection of monomers in reaction simulations. Each functional group is characterized by 
its functionality type (e.g., 'vinyl', 'mono', 'di_different', 'di_identical'), 
SMARTS patterns for substructure matching, and group names for identification. This library serves as a reference 
for the FunctionalGroupsDetector to identify and classify monomers based on their chemical structure.
"""

class FunctionalGroupsLibrary:
    def __init__(self):
        self.monomer_types = {
            "hydroxy_carboxylic_acid_monomer": {
                "functionality_type": "di_different",
                "smarts_1": "[OX2H1;!$(OC=*):1]",
                "smarts_2": "[CX3:2](=[O])[OX2H1]",
                "group_name": "hydroxy_carboxylic_acid"
            },
            "hydroxy_acid_halides_monomer": {
                "functionality_type": "di_different",
                "smarts_1": "[OX2H1;!$(OC=*):1]",
                "smarts_2": "[CX3:2](=[O])[Cl,Br,I]",
                "group_name": "hydroxy_acid_halide",
                "comments": "Hydroxy acid halides are highly reactive and less commonly used monomers for polyesterification compared to hydroxy carboxylic acids."
            },  
            "diol_monomer": {
                "functionality_type": "di_identical",
                "smarts_1": "[O,S;X2;H1;!$([O,S]C=*):3]",
                "group_name": "diol",
                "comments": None,
            },

            "amino_acid_monomer": {
                "functionality_type": "di_different",
                "smarts_1": "[NX3;H2,H1;!$(OC=*):1]",
                "smarts_2": "[CX3:2](=[O])[OX2H1]",
                "group_name": "amino_acid",
                "comments": None,
            },
            "di_amine_monomer": { 
                "functionality_type": "di_identical",
                "smarts_1": "[N&X3;H2,H1;!$(NC=*):3]",
                "group_name": "di_amine",
                "comments": None,
            },
            "di_carboxylic_acid_monomer": {
                "functionality_type": "di_identical",
                "smarts_1": "[CX3:2](=[O])[OX2H1:1]",
                "group_name": "di_carboxylic_acid",
                "comments": None,
            },
            "di_carboxylic_acid_halide_monomer": {
                "functionality_type": "di_identical",
                "smarts_1": "[CX3:2](=[O])[Cl,Br,I:1]",
                "group_name": "di_carboxylic_acid_halide",
                "comments": None,
            },
            "di_epoxide_monomer": {
                "functionality_type": "di_identical",
                "smarts_1": "[C:1]1[O:2][C:3]1",
                "group_name": "di_epoxide"
            },
            
            "di_amine_monomer": {
                "functionality_type": "di_identical",
                "smarts_1": "[N;H2,H1:4]",
                "group_name": "di_amine"
            }
            # "di_cyclic_anhydride_monomer": {
            #     "functionality_type": "di_identical",
            #     "smarts_1": "[CX3,c;R:1](=[OX1])[OX2,o;R][CX3,c;R:2](=[OX1])",
            #     "group_name": "di_cyclic_anhydride"
            # },
            # "di_isocyanate_monomer": {
            #     "functionality_type": "di_identical",
            #     "smarts_1": "[NX2:1]=[CX2]=[OX1,SX1:2]",
            #     "group_name": "di_isocyanate"
            # },
            # "di_epoxide_monomer": {
            #     "functionality_type": "di_identical",
            #     "smarts_1": "[CX4;H2,H1,H0;R:1]1[OX2;R:2][CX4;H1,H0;R:3]1",
            #     "group_name": "di_epoxide"
            # }
            # "vinyl_monomer": {
            #     "functionality_type": "vinyl",
            #     "smarts_1": "[C]=[C;D1]",
            #     "group_name": "vinyl"
            # },
            # "cyclic_olefin_monomer": {
            #     "functionality_type": "vinyl",
            #     "smarts_1": "[CX3;R:1]=[CX3;R:2]",
            #     "group_name": "cyclic_olefin"
            # },
            # "lactone_monomer": {
            #     "functionality_type": "mono",
            #     "smarts_1": "[CX3;R:1](=[OX1])[OX2;R:2]",
            #     "group_name": "lactone"
            # },
            # "cyclic_anhydride_monomer": {
            #     "functionality_type": "mono",
            #     "smarts_1": "[C,c;R:1][CX3,c;R](=[OX1])[OX2,o;R][CX3,c;R](=[OX1])[C,c;R:2]",
            #     "group_name": "cyclic_anhydride"
            # },
            # "epoxide_monomer": {
            #     "functionality_type": "mono",
            #     "smarts_1": "[CX4;R:3]1[OX2;R:4][CX4;R:5]1",
            #     "group_name": "epoxide"
            # },
            # "lactam_monomer": {
            #     "functionality_type": "mono",
            #     "smarts_1": "[CX3;R:1](=[OX1])[NX3;R:2]",
            #     "group_name": "lactam"
            # },
            # "di_amine_monomer": {
            #     "functionality_type": "di_identical",
            #     "smarts_1": "[N&X3;H2,H1;!$(NC=*):3]",

            #     "group_name": "di_amine"
            # },
            # "primery_di_amine_monomer": {
            #     "functionality_type": "di_identical",
            #     "smarts_1": "[C,c:6][NX3;H2;!$(N[C,S]=*)]",
            #     "group_name": "di_primery_amine"
            # },
            # "di_cyclic_anhydride_monomer": {
            #     "functionality_type": "di_identical",
            #     "smarts_1": "[CX3,c;R:1](=[OX1])[OX2,o;R][CX3,c;R:2](=[OX1])",
            #     "group_name": "di_cyclic_anhydride"
            # },
            # "di_isocyanate_monomer": {
            #     "functionality_type": "di_identical",
            #     "smarts_1": "[NX2:1]=[CX2]=[OX1,SX1:2]",
            #     "group_name": "di_isocyanate"
            # },
            # "di_epoxide_monomer": {
            #     "functionality_type": "di_identical",
            #     "smarts_1": "[CX4;H2,H1,H0;R:1]1[OX2;R:2][CX4;H1,H0;R:3]1",
            #     "group_name": "di_epoxide"
            # }
            # need to add more functional groups here from "J. Chem. Inf. Model. 2023, 63, 5539−5548"
        }
        # is there monomers with both COCl and COOH groups?
