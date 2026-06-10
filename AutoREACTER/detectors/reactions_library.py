import json
from typing import Dict, Any


"""
TODO: Missing Polymerization Mechanisms

Step-Growth Condensation / Addition
- Polycarbonates: Diols + Phosgene or Diphenyl Carbonate
- Polyureas: Diamines + Diisocyanates
- Aromatic Polyimides: Dianhydrides + Diamines
- Polybenzimidazoles (PBI): Tetraamines + Dicarboxylates
- Phenol-Formaldehyde (Bakelite): Phenol + Formaldehyde

Aromatic / High-Performance Polymers
- Aromatic Polyethers (PEEK/Sulfones): Activated dihalides + Bisphenols
- Spiro Polymers

Sulfur / Silicon-Based Polymers
- Polysiloxanes (Silicones): Hydrolysis/Condensation of Dichlorosilanes
- Polysulfides: Dihalides + Sodium Sulfide

Click / Ring-Opening / Metathesis Polymerizations
- Thiol-Ene Click Polymerizations
- Ring-Opening Metathesis Polymerization (ROMP)
- Cycloaddition (Four-Center) Reactions

Architectural / Supramolecular Polymers
- Dendritic Polymers: Random Hyperbranched and Dendrimers
- Pseudopolyrotaxanes and Polyrotaxanes

Special Polymerization Environments
- Enzymatic Polymerizations: In Vivo / In Vitro biocatalysis
- Polymerization in Supercritical Carbon Dioxide
"""


class ReactionLibrary:
    def __init__(self):
        self.reactions = {

            # ============================================================
            # Polyesterification: Hydroxy Acids / Acid Halides
            # ============================================================

            "Hydroxy Carboxylic Acid Polycondensation(Polyesterification)": {
                "same_reactants": True,
                "reactant_1": "hydroxy_carboxylic_acid",
                "product": "polyester_chain",
                "delete_atom": True,
                "reaction": "[OX2H1;!$([O][C,S]=*):1]-[H:3].[CX3:2](=[O:5])[OX2H1:4]>>[OX2:1]-[CX3:2](=[O:5]).[O:4]-[H:3]",
                "reference": {
                    "smarts": "https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329",
                    "reaction_and_mechanism": [
                        "https://pubs.acs.org/doi/10.1021/ed048pA734.1",
                        "https://pubs.acs.org/doi/10.1021/ed073pA312"
                    ]
                },
                "comments": None
            },

            "Hydroxy Carboxylic and Hydroxy Carboxylic Polycondensation(Polyesterification)": {
                "same_reactants": False,
                "reactant_1": "hydroxy_carboxylic_acid",
                "reactant_2": "hydroxy_carboxylic_acid",
                "product": "polyester_chain",
                "delete_atom": True,
                "reaction": "[OX2H1;!$([O][C,S]=*):1]-[H:3].[CX3:2](=[O:5])[OX2H1:4]>>[OX2:1]-[CX3:2](=[O:5]).[O:4]-[H:3]",
                "reference": {
                    "smarts": "https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329",
                    "reaction_and_mechanism": [
                        "https://pubs.acs.org/doi/10.1021/ed048pA734.1",
                        "https://pubs.acs.org/doi/10.1021/ed073pA312"
                    ]
                },
                "comments": None
            },

            "Hydroxy Acid Halides Polycondensation(Polyesterification)": {
                "same_reactants": True,
                "reactant_1": "hydroxy_acid_halide",
                "product": "polyester_chain",
                "delete_atom": True,
                "reaction": "[OX2H1;!$([O][C,S]=*):1]-[H:3].[CX3:2](=[O:5])[Cl,Br,I:4]>>[OX2:1]-[CX3:2](=[O:5]).[Cl,Br,I:4]-[H:3]",
                "reference": {
                    "smarts": None,
                    "reaction_and_mechanism": [
                        "https://pubs.acs.org/doi/10.1021/ed073pA312"
                    ]
                },
                "comments": None
            },

            "Hydroxy Acid Halides Hydroxy Acid Halides Polycondensation(Polyesterification)": {
                "same_reactants": False,
                "reactant_1": "hydroxy_acid_halide",
                "reactant_2": "hydroxy_acid_halide",
                "product": "polyester_chain",
                "delete_atom": True,
                "reaction": "[OX2H1;!$([O][C,S]=*):1]-[H:3].[CX3:2](=[O:5])[Cl,Br,I:4]>>[OX2:1]-[CX3:2](=[O:5]).[Cl,Br,I:4]-[H:3]",
                "reference": {
                    "smarts": None,
                    "reaction_and_mechanism": [
                        "https://pubs.acs.org/doi/10.1021/ed073pA312"
                    ]
                },
                "comments": None
            },

            # ============================================================
            # Polyesterification: Diols + Diacids / Diacid Halides / Esters
            # ============================================================

            "Diol and Di-Carboxylic Acid Polycondensation(Polyesterification)": {
                "same_reactants": False,
                "reactant_1": "diol",
                "reactant_2": "di_carboxylic_acid",
                "product": "polyester_chain",
                "delete_atom": True,
                "reaction": "[CX3:1](=[O:3])[OX2H1:4].[OX2H1;!$([O][C,S]=*):2]-[H:5]>>[CX3:1](=[O:3])-[OX2:2].[O:4]-[H:5]",
                "reference": {
                    "smarts": [
                        "https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329"
                    ],
                    "reaction_and_mechanism": [
                        "https://pubs.acs.org/doi/10.1021/ed048pA734.1",
                        "https://pubs.acs.org/doi/10.1021/ed073pA312"
                    ]
                },
                "comments": None
            },

            "Diol and Di-Acid Halide Polycondensation(Polyesterification)": {
                "same_reactants": False,
                "reactant_1": "diol",
                "reactant_2": "di_carboxylic_acid_halide",
                "product": "polyester_chain",
                "delete_atom": True,
                "reaction": "[CX3:1](=[O:3])[Cl,Br,I:4].[OX2H1;!$([O][C,S]=*):2]-[H:5]>>[CX3:1](=[O:3])-[OX2:2].[Cl,Br,I:4]-[H:5]",
                "reference": {
                    "smarts": [
                        "https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329"
                    ],
                    "reaction_and_mechanism": [
                        "https://pubs.acs.org/doi/10.1021/ed048pA734.1",
                        "https://pubs.acs.org/doi/10.1021/ed073pA312"
                    ]
                },
                "comments": None
            },

            "Diol and Di-Carboxylic Ester Polycondensation(Transesterification)": {
                "same_reactants": False,
                "reactant_1": "diol",
                "reactant_2": "di_carboxylic_ester",
                "product": "polyester_chain",
                "delete_atom": True,
                "reaction": "[OX2H1;!$([O][C,S]=*):1]-[H:3].[CX3:2](=[O:5])[OX2H0:4][#6:6]>>[OX2:1]-[CX3:2](=[O:5]).[OX2:4](-[H:3])-[#6:6]",
                "reference": {
                    "smarts": None,
                    "reaction_and_mechanism": None
                },
                "comments": None
            },

            # ============================================================
            # Polyanhydride Formation
            # ============================================================

            "Carboxylic Acid and Acid Halide Polycondensation(Polyanhydride Formation)": {
                "same_reactants": True,
                "reactant_1": "carboxylic_acid_acid_halide",
                "product": "polyanhydride_chain",
                "delete_atom": True,
                "reaction": "[CX3:1](=[O:3])[Cl,Br,I:4].[CX3:2](=[O:5])[OX2H1:6]-[H:7]>>[CX3:1](=[O:3])-[OX2:6]-[CX3:2](=[O:5]).[Cl,Br,I:4]-[H:7]",
                "reference": {
                    "smarts": None,
                    "reaction_and_mechanism": None
                },
                "comments": None
            },

            # ============================================================
            # Polythioesterification
            # ============================================================

            "Dithiol and Di-Carboxylic Acid Halide Polycondensation(Polythioesterification)": {
                "same_reactants": False,
                "reactant_1": "dithiol",
                "reactant_2": "di_carboxylic_acid_halide",
                "product": "polythioester_chain",
                "delete_atom": True,
                "reaction": "[CX3:1](=[O:3])[Cl,Br,I:4].[SX2H1;!$([S][C,S]=*):2]-[H:5]>>[CX3:1](=[O:3])-[SX2:2].[Cl,Br,I:4]-[H:5]",
                "reference": {
                    "smarts": None,
                    "reaction_and_mechanism": None
                },
                "comments": None
            },

            "Dithiol and Di-Carboxylic Acid Polycondensation(Polythioesterification)": {
                "same_reactants": False,
                "reactant_1": "dithiol",
                "reactant_2": "di_carboxylic_acid",
                "product": "polythioester_chain",
                "delete_atom": True,
                "reaction": "[CX3:1](=[O:3])[OX2H1:4].[SX2H1;!$([S][C,S]=*):2]-[H:5]>>[CX3:1](=[O:3])-[SX2:2].[O:4]-[H:5]",
                "reference": {
                    "smarts": None,
                    "reaction_and_mechanism": None
                },
                "comments": "Possible thioesterification with water elimination, but generally less straightforward than acid-halide route."
            },

            # ============================================================
            # Polyamidation
            # ============================================================

            "Amino Acid Polycondensation (Polyamidation)": {
                "same_reactants": True,
                "reactant_1": "amino_acid",
                "product": "polyamide_chain",
                "delete_atom": True,
                "reaction": "[NX3;H2,H1;!$([N][C,S]=*):1]-[H:3].[CX3:2](=[O:4])[OX2H1:5]>>[NX3:1]-[CX3:2](=[O:4]).[O:5]-[H:3]",
                "reference": {
                    "smarts": [
                        "https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329"
                    ],
                    "reaction_and_mechanism": [
                        "https://pubs.acs.org/doi/10.1021/ed048pA734.1",
                        "https://pubs.acs.org/doi/10.1021/ed073pA312"
                    ]
                },
                "comments": None
            },

            "Amino Acid and Amino Acid Polycondensation (Polyamidation)": {
                "same_reactants": False,
                "reactant_1": "amino_acid",
                "reactant_2": "amino_acid",
                "product": "polyamide_chain",
                "delete_atom": True,
                "reaction": "[NX3;H2,H1;!$([N][C,S]=*):1]-[H:3].[CX3:2](=[O:4])[OX2H1:5]>>[NX3:1]-[CX3:2](=[O:4]).[O:5]-[H:3]",
                "reference": {
                    "smarts": [
                        "https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329"
                    ],
                    "reaction_and_mechanism": [
                        "https://pubs.acs.org/doi/10.1021/ed048pA734.1",
                        "https://pubs.acs.org/doi/10.1021/ed073pA312"
                    ]
                },
                "comments": None
            },

            "Di-Amine and Di-Carboxylic Acid Polycondensation (Polyamidation)": {
                "same_reactants": False,
                "reactant_1": "di_amine",
                "reactant_2": "di_carboxylic_acid",
                "product": "polyamide_chain",
                "delete_atom": True,
                "reaction": "[NX3;H2,H1;!$([N][C,S]=*):1]-[H:3].[CX3:2](=[O:4])[OX2H1:5]>>[NX3:1]-[CX3:2](=[O:4]).[O:5]-[H:3]",
                "reference": {
                    "smarts": [
                        "https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329"
                    ],
                    "reaction_and_mechanism": [
                        "https://pubs.acs.org/doi/10.1021/ed048pA734"
                    ]
                },
                "comments": None
            },

            "Di-Amine and Di-Carboxylic Acid Halide Polycondensation (Polyamidation)": {
                "same_reactants": False,
                "reactant_1": "di_amine",
                "reactant_2": "di_carboxylic_acid_halide",
                "product": "polyamide_chain",
                "delete_atom": True,
                "reaction": "[NX3;H2,H1;!$([N][C,S]=*):1]-[H:3].[CX3:2](=[O:4])[Cl,Br,I:5]>>[NX3:1]-[CX3:2](=[O:4]).[Cl,Br,I:5]-[H:3]",
                "reference": {
                    "smarts": [
                        "https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329"
                    ],
                    "reaction_and_mechanism": [
                        "https://pubs.acs.org/doi/10.1021/ed048pA734"
                    ]
                },
                "comments": None
            },

            # ============================================================
            # Mixed Polyester / Polythioester Formation
            # ============================================================

            "Hydroxy-Thiol and Di-Carboxylic Acid Halide Polycondensation through Hydroxy Group": {
                "same_reactants": False,
                "reactant_1": "hydroxy_thiol",
                "reactant_2": "di_carboxylic_acid_halide",
                "product": "mixed_polyester_polythioester_chain",
                "delete_atom": True,
                "reaction": "[CX3:1](=[O:3])[Cl,Br,I:4].[OX2H1;!$([O][C,S]=*):2]-[H:5]>>[CX3:1](=[O:3])-[OX2:2].[Cl,Br,I:4]-[H:5]",
                "reference": {
                    "smarts": None,
                    "reaction_and_mechanism": None
                },
                "comments": None
            },

            "Hydroxy-Thiol and Di-Carboxylic Acid Halide Polycondensation through Thiol Group": {
                "same_reactants": False,
                "reactant_1": "hydroxy_thiol",
                "reactant_2": "di_carboxylic_acid_halide",
                "product": "mixed_polyester_polythioester_chain",
                "delete_atom": True,
                "reaction": "[CX3:1](=[O:3])[Cl,Br,I:4].[SX2H1;!$([S][C,S]=*):2]-[H:5]>>[CX3:1](=[O:3])-[SX2:2].[Cl,Br,I:4]-[H:5]",
                "reference": {
                    "smarts": None,
                    "reaction_and_mechanism": None
                },
                "comments": None
            },

            # ============================================================
            # Polyurethane Formation
            # ============================================================

            "Diol and Di-Isocyanate Polyaddition(Polyurethane Formation)": {
                "same_reactants": False,
                "reactant_1": "diol",
                "reactant_2": "di_isocyanate",
                "product": "polyurethane_chain",
                "delete_atom": False,
                "reaction": "[OX2H1;!$([O][C,S]=*):1]-[H:3].[NX2:4]=[CX2:2]=[OX1:5]>>[OX2:1]-[CX3:2](=[OX1:5])-[NX3:4]-[H:3]",
                "reference": {
                    "smarts": None,
                    "reaction_and_mechanism": None
                },
                "comments": None
            },

            # ============================================================
            # Commented reactions
            # ============================================================

            # "Vinyl Addition Polymerization": {
            #     "same_reactants": True,
            #     "reactant_1": "vinyl",
            #     "product": "polyvinyl_chain",
            #     "delete_atom": False,
            #     "reaction": "[CH2:1]=[CH;H1,H0;!R:2].[CH2:3]=[CH;H1,H0;!R:4]>>[CH2:1]-[CH:2]-[CH2:3]-[CH:4]"
            # },

            # "Cyclic Olefin Addition Polymerization": {
            #     "same_reactants": True,
            #     "reactant_1": "cyclic_olefin",
            #     "product": "polycyclic_chain",
            #     "delete_atom": False,
            #     "reaction": "[CX3;R:1]=[CX3;R:2].[CX3;R:3]=[CX3;R:4]>>[CX4:1]-[CX4:2]-[CX4:3]-[CX4:4]"
            # },

            # "Vinyl Copolymerization": {
            #     "same_reactants": False,
            #     "reactant_1": "vinyl",
            #     "reactant_2": "vinyl",
            #     "product": "copolyvinyl_chain",
            #     "delete_atom": False,
            #     "reaction": "[CH2:1]=[CH;H1,H0;!R:2].[CH2:3]=[CH;H1,H0;!R:4]>>[CH2:1]-[CH:2]-[CH2:3]-[CH:4]"
            # },

            # "Cyclic Olefin and Vinyl Copolymerization": {
            #     "same_reactants": False,
            #     "reactant_1": "vinyl",
            #     "reactant_2": "cyclic_olefin",
            #     "product": "copolycyclicvinyl_chain",
            #     "delete_atom": False,
            #     "reaction": "[CH2:1]=[CH;H1,H0;!R:2].[CX3;R:3]=[CX3;R:4]>>[CX4:1]-[CX4:2]-[CH2:3]-[CH:4]"
            # },

            # "Cyclic Olefin Copolymerization": {
            #     "same_reactants": False,
            #     "reactant_1": "cyclic_olefin",
            #     "reactant_2": "cyclic_olefin",
            #     "product": "copolycyclic_chain",
            #     "delete_atom": False,
            #     "reaction": "[CX3;R:1]=[CX3;R:2].[CX3;R:3]=[CX3;R:4]>>[CX4:1]-[CX4:2]-[CX4:3]-[CX4:4]"
            # },

            # "Lactone Ring-Opening Polyesterification": {
            #     "same_reactants": False,
            #     "reactant_1": "lactone",
            #     "reactant_2": "initiator",
            #     "product": "polyester_chain",
            #     "delete_atom": False
            # },

            # "Cyclic Anhydride and Epoxide Polyesterification": {
            #     "same_reactants": False,
            #     "reactant_1": "cyclic_anhydride_monomer",
            #     "reactant_2": "diol_monomer",
            #     "product": "polyester_chain",
            #     "delete_atom": False
            # },

            # "Cyclic Anhydride and Epoxide Polyetherification": {
            #     "same_reactants": False,
            #     "reactant_1": "cyclic_anhydride",
            #     "reactant_2": "epoxide",
            #     "product": "polyether_chain",
            #     "delete_atom": False
            # },

            # "Epoxide Ring-Opening Polyetherification": {
            #     "same_reactants": False,
            #     "reactant_1": "epoxide",
            #     "reactant_2": "initiator",
            #     "product": "polyether_chain",
            #     "delete_atom": False
            # },

            # "Hindered Phenol Polyetherification": {
            #     "same_reactants": True,
            #     "reactant_1": "hindered_phenol",
            #     "product": "polyether_chain",
            #     "delete_atom": False
            # },

            # "Hindered Phenol Hindered Phenol Polyetherification": {
            #     "same_reactants": False,
            #     "reactant_1": "hindered_phenol",
            #     "reactant_2": "hindered_phenol",
            #     "product": "polyether_chain",
            #     "delete_atom": False
            # },

            # "Bis(p-halogenatedaryl)sulfone Diol (without thiol) polycondensation": {
            #     "same_reactants": False,
            #     "reactant_1": "bis(p-halogenatedaryl)sulfone",
            #     "reactant_2": "diol",
            #     "product": "polyether_chain",
            #     "delete_atom": True
            # },

            # "Bis(bis(p-fluoroaryl)ketone Diol (without thiol) polycondensation": {
            #     "same_reactants": False,
            #     "reactant_1": "bis(p-fluoroaryl)ketone_monomer",
            #     "reactant_2": "diol_monomer",
            #     "product": "polyether_chain",
            #     "delete_atom": True
            # },

            # "Lactam Ring-Opening Polyamidation": {
            #     "same_reactants": True,
            #     "reactant_1": "lactam_monomer",
            #     "product": "polyamide_chain",
            #     "delete_atom": False
            # },

            # "Di-cyclic Anhydride and Di-Primary Amine Polycondensation (Polyimidation)": {
            #     "same_reactants": False,
            #     "reactant_1": "di_cyclic_anhydride_monomer",
            #     "reactant_2": "di_amine_monomer",
            #     "product": "polyimide_chain",
            #     "delete_atom": True
            # },

            # "Di-Epoxide and Di-Isocyanate Polyamination": {
            #     "same_reactants": False,
            #     "reactant_1": "di_epoxide_monomer",
            #     "reactant_2": "di_isocyanate_monomer",
            #     "product": "polyamine_chain",
            #     "delete_atom": False,
            #     "reaction": "[NX2:3]=[CX2:4]=[OX1,SX1:5].[OX2,SX2;H1;!$([O,S]C=*):6]>>[NX3:3][CX3:4](=[OX1,SX1:5])[OX2,SX2;!$([O,S]C=*):6]"
            # },
        }