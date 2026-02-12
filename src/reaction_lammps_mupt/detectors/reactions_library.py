import json
from typing import Dict, Any
# Dictionary defining various polymerization reactions.
# Each reaction includes metadata such as whether reactants are the same,
# the types of reactants, the product type, whether to delete atoms,
# and the SMARTS string for the reaction pattern using RDKit syntax.

"""
TODO: Missing Polymerization Mechanisms
- Polycarbonates: Diols + Phosgene or Diphenyl Carbonate
- Polyureas: Diamines + Diisocyanates
- Aromatic Polyethers (PEEK/Sulfones): Activated dihalides + Bisphenols
- Aromatic Polyimides: Dianhydrides + Diamines
- Polybenzimidazoles (PBI): Tetraamines + Dicarboxylates
- Phenol-Formaldehyde (Bakelite): Phenol + Formaldehyde
- Polysiloxanes (Silicones): Hydrolysis/Condensation of Dichlorosilanes
- Polysulfides: Dihalides + Sodium Sulfide
- Thiol-Ene Click Polymerizations
- Ring-Opening Metathesis Polymerization (ROMP)
- Dendritic Polymers: Random Hyperbranched and Dendrimers
- Enzymatic Polymerizations: In Vivo / In Vitro biocatalysis
- Cycloaddition (Four-Center) Reactions
- Spiro Polymers
- Pseudopolyrotaxanes and Polyrotaxanes
- Polymerization in Supercritical Carbon Dioxide
"""

class ReactionLibrary:
    def __init__(self):
        self.reactions = {
            "Hydroxy Carboxylic Acid Polycondensation(Polyesterification)": {
                "same_reactants": True,
                "reactant_1": "hydroxy_carboxylic_acid",
                "product": "polyester_chain",
                "delete_atom": True,
                "reaction": "[O;!$(OC=*):1]-[H:3].[CX3:2](=[O:5])[OX2H1:4]>>[OX2:1]-[CX3:2](=[O:5]).[O:4]-[H:3]",
                "reference": {"smarts": "https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329",
                            "reaction_and_mechanism": ["https://pubs.acs.org/doi/10.1021/ed048pA734.1", 
                                                        "https://pubs.acs.org/doi/10.1021/ed073pA312"]}
            },
            "Hydroxy Carboxylic and Hydroxy Carboxylic Polycondensation(Polyesterification)": {
                "same_reactants": False,
                "reactant_1": "hydroxy_carboxylic_acid",
                "reactant_2": "hydroxy_carboxylic_acid",
                "product": "polyester_chain",
                "delete_atom": True,
                "reaction": "[O;!$(OC=*):1]-[H:3].[CX3:2](=[O:5])[OX2H1:4]>>[OX2:1]-[CX3:2](=[O:5]).[O:4]-[H:3]",
                "reference": {"smarts": "https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329",
                            "reaction_and_mechanism": ["https://pubs.acs.org/doi/10.1021/ed048pA734.1", 
                                                        "https://pubs.acs.org/doi/10.1021/ed073pA312"]}
            },
            "Hydroxy Acid Halides Hydroxy Acid Halides Polycondensation(Polyesterification)": {
                "same_reactants": False,
                "reactant_1": "hydroxy_acid_halides_monomer",
                "reactant_2": "hydroxy_acid_halides_monomer",
                "product": "polyester_chain",
                "delete_atom": True,
                "reaction": "[O;!$(OC=*):1]-[H:3].[CX3:2](=[O:5])[Cl,Br,I:4]>>[OX2:1]-[CX3:2](=[O:5]).[Cl,Br,I:4]-[H:3]",
                "reference": {"smarts": None,
                            "reaction_and_mechanism": ["https://pubs.acs.org/doi/10.1021/ed073pA312"]},
                "comments": "This is not a commonly documented reaction, but it is a possibke polyesterification between hydroxy acids and carboxylic acids. Since hydroxy acid halides are highly reactive,monomers are less known and not widely used in practice.",
            },
            "Hydroxy Acid Halides Polycondensation(Polyesterification)": {
                "same_reactants": True,
                "reactant_1": "hydroxy_acid_halides_monomer",
                "product": "polyester_chain",
                "delete_atom": True,
                "reaction": "[O;!$(OC=*):1]-[H:3].[CX3:2](=[O:5])[Cl,Br,I:4]>>[OX2:1]-[CX3:2](=[O:5]).[Cl,Br,I:4]-[H:3]",
                "reference": {"smarts": None,
                            "reaction_and_mechanism": ["https://pubs.acs.org/doi/10.1021/ed073pA312"]},
                "comments": "This is not a commonly documented reaction, but it is a possibke polyesterification between hydroxy acids and carboxylic acids. Since hydroxy acid halides are highly reactive,monomers are less known and not widely used in practice."
            },
            "Diol and Di-Acid Halide Polycondensation(Polyesterification)": {
                "same_reactants": False,
                "reactant_1": "diol",
                "reactant_2": "di_carboxylic_acid_halide",
                "product": "polyester_chain",   
                "delete_atom": True,
                "reaction": "[CX3:1](=[O:3])[Cl,Br,I:4].[O,S;X2;H1;!$([O,S]C=*):2]-[H:5]>>[CX3:1](=[O:3])-[O,S;X2;!$([O,S]C=*):2].[Cl,Br,I:4]-[H:5]",
                "reference": {"smarts": ["https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329"],
                            "reaction_and_mechanism": ["https://pubs.acs.org/doi/10.1021/ed048pA734.1", 
                                                        "https://pubs.acs.org/doi/10.1021/ed073pA312"]},
                "comments": None
            },
            "Diol and Di-Carboxylic Acid Polycondensation(Polyesterification)": {
                "same_reactants": False,
                "reactant_1": "diol",
                "reactant_2": "di_carboxylic_acid",
                "product": "polyester_chain",   
                "delete_atom": True,
                "reaction": "[CX3:1](=[O:3])[OX2H1:4].[O,S;X2;H1;!$([O,S]C=*):2]-[H:5]>>[CX3:1](=[O:3])-[O,S;X2;!$([O,S]C=*):2].[O:4]-[H:5]",
                "reference": {"smarts": ["https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329"],
                            "reaction_and_mechanism": ["https://pubs.acs.org/doi/10.1021/ed048pA734.1", 
                                                        "https://pubs.acs.org/doi/10.1021/ed073pA312"]},
                "comments": None
            },
            "Amino Acid Polycondensation (Polyamidation)": {
                "same_reactants": True,
                "reactant_1": "amino_acid",
                "product": "polyamide_chain",   
                "delete_atom": True,
                "reaction": "[NX3;H2,H1;!$(OC=*):1]-[H:3].[CX3:2](=[O:4])[OX2H1:5]>>[NX3:1]-[CX3:2](=[O:4]).[O:5]-[H:3]",
                "reference": {"smarts": ["https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329"],
                            "reaction_and_mechanism": ["https://pubs.acs.org/doi/10.1021/ed048pA734.1", 
                                                        "https://pubs.acs.org/doi/10.1021/ed073pA312"]},
                "comments": None
            },
            "Amino Acid and Amino Acid Polycondensation (Polyamidation)": {
                "same_reactants": False,
                "reactant_1": "amino_acid",
                "reactant_2": "amino_acid",
                "product": "polyamide_chain",   
                "delete_atom": True,
                "reaction": "[NX3;H2,H1;!$(OC=*):1]-[H:3].[CX3:2](=[O:4])[OX2H1:5]>>[NX3:1]-[CX3:2](=[O:4]).[O:5]-[H:3]",
                "reference": {"smarts": ["https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329"],
                            "reaction_and_mechanism": ["https://pubs.acs.org/doi/10.1021/ed048pA734.1", 
                                                        "https://pubs.acs.org/doi/10.1021/ed073pA312"]},
                "comments": None
            },
            "Di-Amine and Di-Carboxylic Acid Polycondensation (Polyamidation)": {
                "same_reactants": False,
                "reactant_1": "di_amine",
                "reactant_2": "di_carboxylic_acid",
                "product": "polyamide_chain",
                "delete_atom": True,
                "reaction": "[CX3:1](=[O:3])[OX2:4].[N;H2,H1;!$(NC=*):2][H:5]>>[CX3:1](=[O:3])[NX3;!$(NC=*):2].[OX2:4]([H:5])",
                "reference": {"smarts": ["https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329"],
                            "reaction_and_mechanism": ["https://pubs.acs.org/doi/10.1021/ed048pA734"]},
            },
            "Di-Amine and Di-Carboxylic Acid Halide Polycondensation (Polyamidation)": {
                "same_reactants": False,
                "reactant_1": "di_amine",
                "reactant_2": "di_carboxylic_acid_halide",
                "product": "polyamide_chain",
                "delete_atom": True,
                "reaction": "[CX3:1](=[O:3])[Cl,Br,I:4].[N;H2,H1;!$(NC=*):2][H:5]>>[CX3:1](=[O:3])[NX3;!$(NC=*):2].[Cl,Br,I:4]([H:5])",
                "reference": {"smarts": ["https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329"],
                            "reaction_and_mechanism": ["https://pubs.acs.org/doi/10.1021/ed048pA734"]},
            }
            # "Di-Amine and Di-Carboxylic Acid Polycondensation (Polyamidation)": {
            #     "same_reactants": False,
            #     "reactant_1": "di_amine_monomer",
            #     "reactant_2": "di_carboxylic_acid_monomer",
            #     "product": "polyamide_chain",
            #     "delete_atom": True,
            #     "reaction": "[CX3:2](=[O])[OX2H1,Cl,Br:1].[N&X3;H2,H1;!$(NC=*):3]>>[CX3:2](=[O])-[NX3;!$(NC=*):3].[OX2H1,Cl,Br:1]" # same product issue as above in polyesterification
            # },
            # "Di-cyclic Anhydride and Di-Primery ammine Polycondensation (Polyimidation)": { # Do not have a clear idea about this reaction yet
            #     "same_reactants": False,
            #     "reactant_1": "di_cyclic_anhydride_monomer",
            #     "reactant_2": "di_isocyanate_monomer",
            #     "product": "polyurethane_chain",
            #     "delete_atom": True
            # },
        #     "Di-Isocyanate and Diol Polyurethane Formation": {
        #         "same_reactants": False,
        #         "reactant_1": "di_isocyanate_monomer",
        #         "reactant_2": "diol_monomer",
        #         "product": "polyurethane_chain",
        #         "delete_atom": True,
        #         "reaction": "[NX3;H1,H0;!$(N[C,S]=*):1].[CX4;H2,H1;!$([CX4](=O)):2]>>[NX3:1]-[CX4:2](=O)"
        #     },
        #     "Di-Epoxide and Di-Isocyanate Polyamination": {
        #         "same_reactants": False,
        #         "reactant_1": "di_epoxide_monomer",
        #         "reactant_2": "di_isocyanate_monomer",
        #         "product": "polyamine_chain",
        #         "delete_atom": False,
        #         "reaction": "[NX2:3]=[CX2:4]=[OX1,SX1:5].[OX2,SX2;H1;!$([O,S]C=*):6]>>[NX3:3][CX3:4](=[OX1,SX1:5])[OX2,SX2;!$([O,S]C=*):6]"
        #     }
        # }
        # "Vinyl Addition Polymerization": {
            #     "same_reactants": True,
            #     "reactant_1": "vinyl",
            #     "product": "polyvinyl_chain",
            #     "delete_atom": False,
            #     "reaction" : "[CH2:1]=[CH;H1,H0;!R:2].[CH2:3]=[CH;H1,H0;!R:4]>>[CH2:1]-[CH:2]-[CH2:3]-[CH:4]" 
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
            #     "reaction" : "[CH2:1]=[CH;H1,H0;!R:2].[CH2:3]=[CH;H1,H0;!R:4]>>[CH2:1]-[CH:2]-[CH2:3]-[CH:4]"
            # },
            # "Cyclic Olefin and Vinyl Copolymerization": {
            #     "same_reactants": False,
            #     "reactant_1": "vinyl",
            #     "reactant_2": "cyclic_olefin",
            #     "product": "copolycyclicvinyl_chain",
            #     "delete_atom": False,
            #     "reaction": "[CH2:1]=[CH;H1,H0;!R:2].[CX3;R:3]=[CX3;R:4]>>[CX4:1]-[CX4:2]-[CH2:3]-[CH:4]" # has to give reactants as in the order 
            # },
            # "Cyclic Olefin Copolymerization": {
            #     "same_reactants": False,
            #     "reactant_1": "cyclic_olefin",
            #     "reactant_2": "cyclic_olefin",
            #     "product": "copolycyclic_chain",
            #     "delete_atom": False,
            #     "reaction": "[CX3;R:1]=[CX3;R:2].[CX3;R:3]=[CX3;R:4]>>[CX4:1]-[CX4:2]-[CX4:3]-[CX4:4]"
            # },
            # "Lactone Ring-Opening Polyesterification": { # does not work yet
            #     "reactant_1": "lactone",
            #     "reactant_2": "initiator",
            #     "product": "polyester_chain"
            # },
            # "Cyclic Anhydride and Epoxide Polyesterification": {
            #     "same_reactants": False,
            #     "reactant_1": "cyclic_anhydride_monomer",
            #     "reactant_2": "diol_monomer",
            #     "product": "polyester_chain",
            #     "delete_atom": False # Not sure about this
            # },
            # "Cyclic Anhydride and Epoxide Polyetherification": { # Not sure untill tested
            #     "same_reactants": False,
            #     "reactant_1": "cyclic_anhydride",
            #     "reactant_2": "epoxide",   
            #     "product": "polyether_chain",
            #     "delete_atom": False # Not sure about this
            # },
            # "Epoxide Ring-Opening Polyetherification": { # Do not have a clear idea about this reaction yet
            #     "same_reactants": False,
            #     "reactant_1": "epoxide",
            #     "reactant_2": "initiator",
            #     "product": "polyether_chain",
            #     "delete_atom": False
            # },
            # "Hindered Phenol Polyetherification": { # Same as above
            #     "same_reactants": True,
            #     "reactant_1": "hindered_phenol",
            #     "product": "polyether_chain",
            #     "delete_atom": False
            # },
            # "Hindered Phenol Hindered Phenol Polyetherification": { # Same as above
            #     "same_reactants": False,
            #     "reactant_1": "hindered_phenol",
            #     "reactant_2": "hindered_phenol",
            #     "product": "polyether_chain",
            #     "delete_atom": False
            # },
            # "Bis(p-halogenatedaryl)sulfone Diol (without thiol) polycondensation": { # Same as above
            #     "same_reactants": False,
            #     "reactant_1": "bis(p-halogenatedaryl)sulfone",
            #     "reactant_2": "diol",
            #     "product": "polyether_chain",
            #     "delete_atom": True
            # },
            # "Bis(bis(p-fluoroaryl)ketone Diol (without thiol) polycondensation": { # Same as above
            #     "same_reactants": False,
            #     "reactant_1": "bis(p-fluoroaryl)ketone_monomer",
            #     "reactant_2": "diol_monomer",
            #     "product": "polyether_chain",
            #     "delete_atom": True
            # },
            # "Bis(bis(p-fluoroaryl)ketone Diol (without thiol) polycondensation": { # Same as above
            #     "same_reactants": False,
            #     "reactant_1": "bis(p-fluoroaryl)ketone_monomer",
            #     "reactant_2": "diol_monomer",
            #     "product": "polyether_chain",
            #     "delete_atom": True
            # },
            # "Lactam Ring-Opening Polyamidation": { # does not work yet
            #     "same_reactants": True,
            #     "reactant_1": "lactam_monomer",
            #     "product": "polyamide_chain",
            #     "delete_atom": False
            # },
            # "Di-Amine and Di-Carboxylic Acid Polycondensation (Polyamidation)": {
            #     "same_reactants": False,
            #     "reactant_1": "di_amine_monomer",
            #     "reactant_2": "di_carboxylic_acid_monomer",
            #     "product": "polyamide_chain",
            #     "delete_atom": True,
            #     "reaction": "[CX3:2](=[O])[OX2H1,Cl,Br:1].[N&X3;H2,H1;!$(NC=*):3]>>[CX3:2](=[O])-[NX3;!$(NC=*):3].[OX2H1,Cl,Br:1]" # same product issue as above in polyesterification
            # },
            # "Di-cyclic Anhydride and Di-Primery ammine Polycondensation (Polyimidation)": { # Do not have a clear idea about this reaction yet
            #     "same_reactants": False,
            #     "reactant_1": "di_cyclic_anhydride_monomer",
            #     "reactant_2": "di_isocyanate_monomer",
            #     "product": "polyurethane_chain",
            #     "delete_atom": True
            # },
        #     "Di-Isocyanate and Diol Polyurethane Formation": {
        #         "same_reactants": False,
        #         "reactant_1": "di_isocyanate_monomer",
        #         "reactant_2": "diol_monomer",
        #         "product": "polyurethane_chain",
        #         "delete_atom": True,
        #         "reaction": "[NX3;H1,H0;!$(N[C,S]=*):1].[CX4;H2,H1;!$([CX4](=O)):2]>>[NX3:1]-[CX4:2](=O)"
        #     },
        #     "Di-Epoxide and Di-Isocyanate Polyamination": {
        #         "same_reactants": False,
        #         "reactant_1": "di_epoxide_monomer",
        #         "reactant_2": "di_isocyanate_monomer",
        #         "product": "polyamine_chain",
        #         "delete_atom": False,
        #         "reaction": "[NX2:3]=[CX2:4]=[OX1,SX1:5].[OX2,SX2;H1;!$([O,S]C=*):6]>>[NX3:3][CX3:4](=[OX1,SX1:5])[OX2,SX2;!$([O,S]C=*):6]"
        #     }
        # }
        }
