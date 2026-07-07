"""
Functional group library for epoxy polymerization chemistry.

This library focuses on monomers and curing-agent functional groups relevant to
epoxy polymerization / epoxy curing systems. It includes epoxides and common
epoxy-reactive groups such as primary amines, secondary amines, thiols, alcohols,
carboxylic acids, and cyclic anhydrides.

Each entry defines:
    - functionality_type:
        "mono"          : one reactive functional group
        "di_identical"  : two identical reactive functional groups
        "di_different"  : two different reactive functional groups

    - smarts_1 / smarts_2:
        SMARTS patterns used for substructure matching

    - group_name:
        functional group label used by the detector

    - comments:
        optional chemistry notes
"""



class FunctionalGroupsLibrary:
    def __init__(self):
        self.monomer_types = {

            # ============================================================
            # Hydroxy / Carboxylic Acid AB-Type Monomers
            # ============================================================

            "hydroxy_carboxylic_acid_monomer": {
                "functionality_type": "di_different",
                "smarts_1": "[OX2H1;!$([O][C,S]=*):1]",
                "smarts_2": "[CX3:2](=[O])[OX2H1]",
                "group_name": "hydroxy_carboxylic_acid",
                "comments": None,
            },

            "hydroxy_acid_halides_monomer": {
                "functionality_type": "di_different",
                "smarts_1": "[OX2H1;!$([O][C,S]=*):1]",
                "smarts_2": "[CX3:2](=[O])[Cl,Br,I]",
                "group_name": "hydroxy_acid_halide",
                "comments": "Hydroxy acid halides are highly reactive and less commonly used monomers for polyesterification compared to hydroxy carboxylic acids."
            },

            # ============================================================
            # Alcohol / Thiol Functional Monomers
            # ============================================================

            "diol_monomer": {
                "functionality_type": "di_identical",
                "smarts_1": "[OX2H1;!$([O][C,S]=*):1]",
                "group_name": "diol",
                "comments": None,
            },

            "dithiol_monomer": {
                "functionality_type": "di_identical",
                "smarts_1": "[SX2H1;!$([S][C,S]=*):1]",
                "group_name": "dithiol",
                "comments": None,
            },

            "hydroxy_thiol_monomer": {
                "functionality_type": "di_different",
                "smarts_1": "[OX2H1;!$([O][C,S]=*):1]",
                "smarts_2": "[SX2H1;!$([S][C,S]=*):2]",
                "group_name": "hydroxy_thiol",
                "comments": None,
            },

            # ============================================================
            # Amine / Amino Acid Monomers
            # ============================================================

            "amino_acid_monomer": {
                "functionality_type": "di_different",
                "smarts_1": "[NX3;H2,H1;!$([N][C,S]=*):1]",
                "smarts_2": "[CX3:2](=[O])[OX2H1]",
                "group_name": "amino_acid",
                "comments": None,
            },

            "di_amine_monomer": { 
                "functionality_type": "di_identical",
                "smarts_1": "[NX3;H2,H1;!$([N][C,S]=*):1]",
                "group_name": "di_amine",
                "comments": None,
            },

            # ============================================================
            # Carboxylic Acid / Acid Halide / Ester Monomers
            # ============================================================

            "di_carboxylic_acid_monomer": {
                "functionality_type": "di_identical",
                "smarts_1": "[CX3:1](=[O])[OX2H1]",
                "group_name": "di_carboxylic_acid",
                "comments": None,
            },

            "di_carboxylic_acid_halide_monomer": {
                "functionality_type": "di_identical",
                "smarts_1": "[CX3:1](=[O])[Cl,Br,I]",
                "group_name": "di_carboxylic_acid_halide",
                "comments": None,
            },

            "carboxylic_acid_acid_halide_monomer": {
                "functionality_type": "di_different",
                "smarts_1": "[CX3:1](=[O])[OX2H1]",
                "smarts_2": "[CX3:2](=[O])[Cl,Br,I]",
                "group_name": "carboxylic_acid_acid_halide",
                "comments": "Mixed COOH/acid-halide AB monomer. Edge case; forms polyanhydride-type linkage, not polyester.",
            },

            "di_carboxylic_ester_monomer": {
                "functionality_type": "di_identical",
                "smarts_1": "[CX3:1](=[O])[OX2H0][#6]",
                "group_name": "di_carboxylic_ester",
                "comments": None,
            },

            # ============================================================
            # Isocyanate Monomers
            # ============================================================

            "di_isocyanate_monomer": {
                "functionality_type": "di_identical",
                "smarts_1": "[NX2]=[CX2:1]=[OX1]",
                "group_name": "di_isocyanate",
                "comments": None,
            },

            # ============================================================
            # Epoxy / Amine Functional Monomers
            # Relevant for epoxy-amine polymerization
            #
            # Polymer-forming epoxy monomer:
            #     must contain two epoxide groups, i.e. diepoxy
            #
            # Primary monoamine:
            #     one -NH2 group has two active hydrogens
            #     can react with two epoxide groups in two stages
            #
            # Stage 1:
            #     primary amine + epoxide -> secondary amine
            #
            # Stage 2:
            #     secondary amine + epoxide -> tertiary amine
            # ============================================================

            "di_epoxy_monomer": {
                "functionality_type": "di_identical",
                "smarts_1": "[OX2r3:1]1[#6r3:2][#6r3:3]1",
                "group_name": "di_epoxide",
                "comments": "Difunctional epoxide monomer. Required on the epoxy side for epoxy-amine polymerization.",
            },

            "primary_amine_monomer": {
                "functionality_type": "mono",
                "smarts_1": "[NX3;H2;!$([N][C,S]=*):1]",
                "group_name": "primary_amine",
                "comments": "Mono primary amine. One -NH2 group has two active hydrogens and can react with two epoxide groups.",
            },

            "secondary_amine_monomer": {
                "functionality_type": "mono",
                "smarts_1": "[NX3;H1;!$([N][C,S]=*):1]",
                "group_name": "secondary_amine",
                "comments": "Mono secondary amine. Represents the second-stage reactive amine after primary amine reacts once with epoxide. By itself, a mono secondary amine reacts only once with epoxide and is not a true polymer-forming monomer.",
            },

            # ============================================================
            # Commented functional groups
            # ============================================================

            # ------------------------------------------------------------
            # Cyclic Anhydride / Epoxide Functional Groups
            # ------------------------------------------------------------

            # "di_cyclic_anhydride_monomer": {
            #     "functionality_type": "di_identical",
            #     "smarts_1": "[CX3,c;R:1](=[OX1])[OX2,o;R][CX3,c;R:2](=[OX1])",
            #     "group_name": "di_cyclic_anhydride"
            # },

            # "di_epoxide_monomer": {
            #     "functionality_type": "di_identical",
            #     "smarts_1": "[CX4;H2,H1,H0;R:1]1[OX2;R:2][CX4;H1,H0;R:3]1",
            #     "group_name": "di_epoxide"
            # },

            # ------------------------------------------------------------
            # Vinyl / Olefin Functional Groups
            # ------------------------------------------------------------

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

            # ------------------------------------------------------------
            # Ring-Opening Functional Groups
            # ------------------------------------------------------------

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
        }            #     "smarts_1": "[CX3;R:1](=[OX1])[NX3;R:2]",
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
        # is there monomers with both COCl and COOH groups?
