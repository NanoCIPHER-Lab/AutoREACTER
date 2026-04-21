import json
from typing import Dict, Any
# Dictionary defining various polymerization reactions.
# Each reaction includes metadata such as whether reactants are the same,
# the types of reactants, the product type, whether to delete atoms,
# and the SMARTS string for the reaction pattern using RDKit syntax.
"""
Move to below structure for better organization
    reactions/
        polyester.py
        polyamide.py
        polyurethane.py
        polyurea.py
        epoxy.py
        vinyl.py
        ring_opening.py
        aromatic_condensation.py
This will allow us to modularly manage and expand the reaction library, 
making it easier to maintain and update specific reaction types without affecting the entire library. 
Each file can focus on a specific class of polymerization reactions, providing detailed SMARTS patterns 
and references for those reactions. Same will go with functional groups library, we can have separate files 
for different classes of functional groups,
    functional_groups/
        polyester.py
        polyamide.py
        polyurethane.py
        polyurea.py
        epoxy.py
        vinyl.py
        ring_opening.py
        aromatic_condensation.py

May be a data class for reaction definition can be created to standardize the structure of reaction entries,
@dataclass
class ReactionDef:
    reaction_id: str
    display_name: str
    family: str
    mechanism: str
    reactant_1: str
    reactant_2: Optional[str]
    product: str
    smarts: str
    delete_atom: bool
    symmetric: bool = False
    byproduct: Optional[str] = None
    references: Dict[str, List[str]] = field(default_factory=dict)
    comments: Optional[str] = None
which can be used to create instances for each reaction, 
ensuring consistency and making it easier to manage the reaction 
library rather than having a flat dictionary structure.
"""

"""
TODO: Missing Polymerization Mechanisms
- Polycarbonates: Diols + Phosgene or Diphenyl Carbonate
- Polyurethanes: Diols + Diisocyanates
- Polyureas: Diamines + Diisocyanates
- Aromatic Polyethers (PEEK/Sulfones): Activated dihalides + Bisphenols
- Aromatic Polyimides: Dianhydrides + Diamines
- Polybenzimidazoles (PBI): Tetraamines + Dicarboxylates
- Phenol-Formaldehyde (Bakelite): Phenol + Formaldehyde
- Urea-Formaldehyde Resins
- Melamine-Formaldehyde Resins
- Polyacetals / Polyketals: Diols + Aldehydes/Ketones
- Polyazomethines (Schiff-Base Polymers): Diamines + Dialdehydes
- Polyhydrazones / Polyoximes: Dicarbonyls + Dihydrazides / Dihydroxylamines
- Polythioesters: Dithiols + Diacid Halides
- Polythiocarbonates / Polythiourethanes: Thiols + Carbonyl / Isocyanate Derivatives
- Polyphosphazenes
- Polyphosphoesters / Polyphosphonates
- Sol-Gel / Silsesquioxane Network Formation
- Lactone Ring-Opening Polymerization
- Lactide Ring-Opening Polymerization
- Lactam Ring-Opening Polymerization
- Epoxide / Oxirane Ring-Opening Polymerization
- Cyclic Ether Ring-Opening Polymerization
- Cyclic Carbonate Ring-Opening Polymerization
- Cyclic Phosphate / Phosphonate Ring-Opening Polymerization
- Cyclic Siloxane Ring-Opening Polymerization
- N-Carboxyanhydride (NCA) Polymerization: Polypeptides
- O-Carboxyanhydride (OCA) Polymerization
- Ring-Opening Copolymerization (ROCOP): Epoxide + CO2
- Ring-Opening Copolymerization (ROCOP): Epoxide + Cyclic Anhydride
- Ring-Opening Copolymerization (ROCOP): Epoxide + COS / CS2 / Isothiocyanates
- Ring-Opening Alternating Polymerization (ROAP): Cyclic Anhydride + Cyclic Ether
- Thiol-Ene Click Polymerizations
- Thiol-Yne Polymerizations
- Azide-Alkyne Cycloaddition (CuAAC / SPAAC) Polymerizations
- Thiol-Michael Polymerizations
- Aza-Michael Polymerizations
- SuFEx Polymerizations
- Diels-Alder Polymerizations
- Imine / Hydrazone / Oxime Dynamic Covalent Polymerizations
- Boronic Ester / Boroxine Dynamic Covalent Polymerizations
- Disulfide-Exchange Polymerizations
- Polysiloxanes (Silicones): Hydrolysis/Condensation of Dichlorosilanes
- Polysilanes
- Polycarbosilanes
- Polysulfides: Dihalides + Sodium Sulfide
- Polysulfones / Polysulfoxides
- Polysulfates / Polysulfonates
- Thiol-Ene Click Polymerizations
- Ring-Opening Metathesis Polymerization (ROMP)
- Dendritic Polymers: Random Hyperbranched and Dendrimers
- Enzymatic Polymerizations: In Vivo / In Vitro biocatalysis
- Cycloaddition (Four-Center) Reactions
- Spiro Polymers
- Pseudopolyrotaxanes and Polyrotaxanes
- Polymerization in Supercritical Carbon Dioxide
- CF2=CF2 (Tetrafluoroethylene) Polymerization (different from CH2=CHY mechanism)
"""

class ReactionLibrary:
    def __init__(self):
        self.reactions = {
            
            # ==========================================
            # POLYESTERIFICATION REACTIONS
            # ==========================================
            "Hydroxy Carboxylic Acid Polycondensation(Polyesterification)": {
                "same_reactants": True,
                "reactant_1": "hydroxy_carboxylic_acid",
                "product": "polyester_chain",
                "delete_atom": True,
                "reaction": "[O;!$(OC=*):1]-[H:3].[CX3:2](=[O:4])[O:5]-[H:6]>>[O:1]-[CX3:2](=[O:4]).[H:3]-[O:5]-[H:6]",
                "reference": {
                    "smarts": "https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329",
                    "reaction_and_mechanism": [
                        "https://pubs.acs.org/doi/10.1021/ed048pA734.1", 
                        "https://pubs.acs.org/doi/10.1021/ed073pA312"
                    ]
                }
            },
            "Hydroxy Carboxylic and Hydroxy Carboxylic Polycondensation(Polyesterification)": {
                "same_reactants": False,
                "reactant_1": "hydroxy_carboxylic_acid",
                "reactant_2": "hydroxy_carboxylic_acid",
                "product": "polyester_chain",
                "delete_atom": True,
                "reaction": "[O;!$(OC=*):1]-[H:3].[CX3:2](=[O:4])[O:5]-[H:6]>>[O:1]-[CX3:2](=[O:4]).[H:3]-[O:5]-[H:6]",
                "reference": {
                    "smarts": "https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329",
                    "reaction_and_mechanism": [
                        "https://pubs.acs.org/doi/10.1021/ed048pA734.1", 
                        "https://pubs.acs.org/doi/10.1021/ed073pA312"
                    ]
                }
            },
            "Hydroxy Acid Halides Polycondensation(Polyesterification)": {
                "same_reactants": True,
                "reactant_1": "hydroxy_acid_halide",
                "product": "polyester_chain",
                "delete_atom": True,
                "reaction": "[O;!$(OC=*):1]-[H:3].[CX3:2](=[O:4])[Cl,Br,I:5]>>[O:1]-[CX3:2](=[O:4]).[Cl,Br,I:5]-[H:3]",
                "reference": {
                    "smarts": None,
                    "reaction_and_mechanism": ["https://pubs.acs.org/doi/10.1021/ed073pA312"]
                },
                "comments": "This is not a commonly documented reaction, but it is a possible polyesterification between hydroxy acids and carboxylic acids. Since hydroxy acid halides are highly reactive, monomers are less known and not widely used in practice."
            },
            "Hydroxy Acid Halides Hydroxy Acid Halides Polycondensation(Polyesterification)": {
                "same_reactants": False,
                "reactant_1": "hydroxy_acid_halide",
                "reactant_2": "hydroxy_acid_halide",
                "product": "polyester_chain",
                "delete_atom": True,
                "reaction": "[O;!$(OC=*):1]-[H:3].[CX3:2](=[O:4])[Cl,Br,I:5]>>[O:1]-[CX3:2](=[O:4]).[Cl,Br,I:5]-[H:3]",
                "reference": {
                    "smarts": None,
                    "reaction_and_mechanism": ["https://pubs.acs.org/doi/10.1021/ed073pA312"]
                },
                "comments": "This is not a commonly documented reaction, but it is a possible polyesterification between hydroxy acids and carboxylic acids. Since hydroxy acid halides are highly reactive, monomers are less known and not widely used in practice."
            },
            "Diol and Di-Carboxylic Acid Polycondensation(Polyesterification)": {
                "same_reactants": False,
                "reactant_1": "diol",
                "reactant_2": "di_carboxylic_acid",
                "product": "polyester_chain",   
                "delete_atom": True,
                "reaction": "[O,S;X2;!$([O,S]C=*):1]-[H:3].[CX3:2](=[O:4])[O:5]-[H:6]>>[O,S;X2;!$([O,S]C=*):1]-[CX3:2](=[O:4]).[H:3]-[O:5]-[H:6]",
                "reference": {
                    "smarts": ["https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329"],
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
                "reaction": "[O,S;X2;!$([O,S]C=*):1]-[H:3].[CX3:2](=[O:4])[Cl,Br,I:5]>>[O,S;X2;!$([O,S]C=*):1]-[CX3:2](=[O:4]).[Cl,Br,I:5]-[H:3]",
                "reference": {
                    "smarts": ["https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329"],
                    "reaction_and_mechanism": [
                        "https://pubs.acs.org/doi/10.1021/ed048pA734.1", 
                        "https://pubs.acs.org/doi/10.1021/ed073pA312"
                    ]
                },
                "comments": None
            },

            # ==========================================
            # POLYAMIDATION REACTIONS
            # ==========================================
            "Amino Acid Polycondensation (Polyamidation)": {
                "same_reactants": True,
                "reactant_1": "amino_acid",
                "product": "polyamide_chain",   
                "delete_atom": True,
                "reaction": "[NX3;H2,H1;!$(OC=*):1]-[H:3].[CX3:2](=[O:4])[O:5]-[H:6]>>[NX3:1]-[CX3:2](=[O:4]).[H:3]-[O:5]-[H:6]",
                "reference": {
                    "smarts": ["https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329"],
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
                "reaction": "[NX3;H2,H1;!$(OC=*):1]-[H:3].[CX3:2](=[O:4])[O:5]-[H:6]>>[NX3:1]-[CX3:2](=[O:4]).[H:3]-[O:5]-[H:6]",
                "reference": {
                    "smarts": ["https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329"],
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
                "reaction": "[N;H2,H1;!$(NC=*):1]-[H:3].[CX3:2](=[O:4])[O:5]-[H:6]>>[N;!$(NC=*):1]-[CX3:2](=[O:4]).[H:3]-[O:5]-[H:6]",
                "reference": {
                    "smarts": ["https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329"],
                    "reaction_and_mechanism": ["https://pubs.acs.org/doi/10.1021/ed048pA734"]
                },
            },
            "Di-Amine and Di-Carboxylic Acid Halide Polycondensation (Polyamidation)": {
                "same_reactants": False,
                "reactant_1": "di_amine",
                "reactant_2": "di_carboxylic_acid_halide",
                "product": "polyamide_chain",
                "delete_atom": True,
                "reaction": "[N;H2,H1;!$(NC=*):1]-[H:3].[CX3:2](=[O:4])[Cl,Br,I:5]>>[N;!$(NC=*):1]-[CX3:2](=[O:4]).[Cl,Br,I:5]-[H:3]",
                "reference": {
                    "smarts": ["https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329"],
                    "reaction_and_mechanism": ["https://pubs.acs.org/doi/10.1021/ed048pA734"]
                },
            },

            # ==========================================
            # POLYAMINATION / EPOXY CURING
            # ==========================================
            "Di-Epoxide and Di-Amine Polyamination": {
                "same_reactants": False,
                "reactant_1": "di_epoxide",
                "reactant_2": "di_amine",
                "product": "polyamine_chain",
                "delete_atom": False,
                "reaction": "[C:4]1[O:3][C:2]1.[N:1]-[H:5]>>[H:5]-[O:3]-[C:4]-[C:2]-[N:1]",
                "reference": {},
            }

            # ==========================================
            # COMMENTED OUT REACTIONS
            # ==========================================
            
            # --- POLYIMIDATION ---
            # "Di-cyclic Anhydride and Di-Primary Amine Polycondensation (Polyimidation)": { 
            #     "same_reactants": False,
            #     "reactant_1": "di_cyclic_anhydride_monomer",
            #     "reactant_2": "di_amine_monomer",
            #     "product": "polyimide_chain",
            #     "delete_atom": True,
            #     "reaction": "[N;H2:1](-[H:3])-[H:4].[CX3;R:2](=[O:5])[OX2;R:6][CX3;R:7](=[O:8])>>[CX3:2](=[O:5])-[N:1]-[CX3:7](=[O:8]).[H:3]-[O:6]-[H:4]"
            # },

            # --- POLYURETHANES / UREAS ---
            # "Di-Isocyanate and Diol Polyurethane Formation": {
            #     "same_reactants": False,
            #     "reactant_1": "di_isocyanate_monomer",
            #     "reactant_2": "diol_monomer",
            #     "product": "polyurethane_chain",
            #     "delete_atom": False,
            #     "reaction": "[O,S;!$([O,S]C=*):1]-[H:5].[C:2](=[O,S:3])=[N:4]>>[O,S:1]-[C:2](=[O,S:3])-[N:4]-[H:5]"
            # },
            # "Di-Epoxide and Di-Isocyanate Polyamination": {
            #     "same_reactants": False,
            #     "reactant_1": "di_epoxide_monomer",
            #     "reactant_2": "di_isocyanate_monomer",
            #     "product": "polyamine_chain",
            #     "delete_atom": False,
            #     "reaction": "[C:4]1[O:3][C:2]1.[N:1]=[C:6]=[O,S:5]>>[O:3]1-[C:4]-[C:2]-[N:1](-[C:6](=[O,S:5]))-1" # Simplified placeholder
            # },

            # --- VINYL / ADDITION POLYMERIZATION ---
            # "Vinyl Addition Polymerization": {
            #     "same_reactants": True,
            #     "reactant_1": "vinyl",
            #     "product": "polyvinyl_chain",
            #     "delete_atom": False,
            #     "reaction" : "[CH2:1]=[CH;H1,H0;!R:3].[CH2:4]=[CH;H1,H0;!R:2]>>[CH2:4]-[CH:2]-[CH2:1]-[CH:3]" 
            # },
            # "Vinyl Copolymerization": {
            #     "same_reactants": False,
            #     "reactant_1": "vinyl",
            #     "reactant_2": "vinyl",
            #     "product": "copolyvinyl_chain",
            #     "delete_atom": False,
            #     "reaction" : "[CH2:1]=[CH;H1,H0;!R:3].[CH2:4]=[CH;H1,H0;!R:2]>>[CH2:4]-[CH:2]-[CH2:1]-[CH:3]"
            # },
            # "Cyclic Olefin Addition Polymerization": {
            #     "same_reactants": True,
            #     "reactant_1": "cyclic_olefin",   
            #     "product": "polycyclic_chain",
            #     "delete_atom": False,
            #     "reaction": "[CX3;R:1]=[CX3;R:3].[CX3;R:4]=[CX3;R:2]>>[CX4:4]-[CX4:2]-[CX4:1]-[CX4:3]"
            # },
            # "Cyclic Olefin Copolymerization": {
            #     "same_reactants": False,
            #     "reactant_1": "cyclic_olefin",
            #     "reactant_2": "cyclic_olefin",
            #     "product": "copolycyclic_chain",
            #     "delete_atom": False,
            #     "reaction": "[CX3;R:1]=[CX3;R:3].[CX3;R:4]=[CX3;R:2]>>[CX4:4]-[CX4:2]-[CX4:1]-[CX4:3]"
            # },
            # "Cyclic Olefin and Vinyl Copolymerization": {
            #     "same_reactants": False,
            #     "reactant_1": "vinyl",
            #     "reactant_2": "cyclic_olefin",
            #     "product": "copolycyclicvinyl_chain",
            #     "delete_atom": False,
            #     "reaction": "[CH2:1]=[CH;H1,H0;!R:3].[CX3;R:4]=[CX3;R:2]>>[CX4:4]-[CX4:2]-[CH2:1]-[CH:3]" 
            # },

            # --- RING-OPENING / POLYETHERS ---
            # "Lactone Ring-Opening Polyesterification": {
            #     "same_reactants": False,
            #     "reactant_1": "lactone",
            #     "reactant_2": "initiator",
            #     "product": "polyester_chain",
            #     "delete_atom": False,
            #     "reaction": "[O:1]-[H:3].[CX3;R:2](=[O:4])[OX2;R:5]-[C:6]>>[O:1]-[CX3:2](=[O:4]).[O:5](-[C:6])-[H:3]"
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
            #     "delete_atom": False,
            #     "reaction": "[O:1]-[H:5].[C:4]1[O:3][C:2]1>>[H:5]-[O:3]-[C:4]-[C:2]-[O:1]"
            # },
            # "Hindered Phenol Polyetherification": { 
            #     "same_reactants": True,
            #     "reactant_1": "hindered_phenol",
            #     "product": "polyether_chain",
            #     "delete_atom": True,
            #     "reaction": "[c:1]-[H:3].[c:4]-[O:2]-[H:5]>>[c:1]-[O:2]-[c:4].[H:3]-[H:5]"
            # },
            # "Hindered Phenol Hindered Phenol Polyetherification": { 
            #     "same_reactants": False,
            #     "reactant_1": "hindered_phenol",
            #     "reactant_2": "hindered_phenol",
            #     "product": "polyether_chain",
            #     "delete_atom": True,
            #     "reaction": "[c:1]-[H:3].[c:4]-[O:2]-[H:5]>>[c:1]-[O:2]-[c:4].[H:3]-[H:5]"
            # },
            
            # --- ADVANCED POLYCONDENSATIONS (PEEK/SULFONES) ---
            # "Bis(p-halogenatedaryl)sulfone Diol (without thiol) polycondensation": { 
            #     "same_reactants": False,
            #     "reactant_1": "bis(p-halogenatedaryl)sulfone",
            #     "reactant_2": "diol",
            #     "product": "polyether_chain",
            #     "delete_atom": True,
            #     "reaction": "[O:1]-[H:3].[c:2]-[Cl,F:4]>>[O:1]-[c:2].[Cl,F:4]-[H:3]"
            # },
            # "Bis(p-fluoroaryl)ketone Diol (without thiol) polycondensation": { 
            #     "same_reactants": False,
            #     "reactant_1": "bis(p-fluoroaryl)ketone_monomer",
            #     "reactant_2": "diol_monomer",
            #     "product": "polyether_chain",
            #     "delete_atom": True,
            #     "reaction": "[O:1]-[H:3].[c:2]-[F:4]>>[O:1]-[c:2].[F:4]-[H:3]"
            # },
            # "Lactam Ring-Opening Polyamidation": { 
            #     "same_reactants": True,
            #     "reactant_1": "lactam_monomer",
            #     "product": "polyamide_chain",
            #     "delete_atom": False,
            #     "reaction": "[N:1]-[H:3].[CX3;R:2](=[O:4])[NX3;R:5]-[C:6]>>[N:1]-[CX3:2](=[O:4]).[N:5](-[C:6])-[H:3]"
            # }
        }
            # ADVANCED POLYMERIZATION REACTIONS

            # --- CLICK CHEMISTRY (CuAAC) ---
            # "Azide-Alkyne Cycloaddition (Click Polymerization)": {
            #     "same_reactants": False,
            #     "reactant_1": "di_azide",
            #     "reactant_2": "di_alkyne",
            #     "product": "polytriazole_chain",
            #     "delete_atom": False,
            #     # Initiators: N:1 and C:2. Forms a 1,2,3-triazole ring.
            #     "reaction": "[N:1]=[N+:3]=[N-:4].[C:5]#[C:2]>>[N:1]1-[C:2]=[C:5]-[N:4]=[N:3]-1"
            # },

            # --- SUFEX POLYMERIZATION ---
            # "Sulfur-Fluoride Exchange (SuFEx) Polysulfate Formation": {
            #     "same_reactants": False,
            #     "reactant_1": "bis_sulfonyl_fluoride",
            #     "reactant_2": "bis_silyl_ether",
            #     "product": "polysulfate_chain",
            #     "delete_atom": True,
            #     # Initiators: S:1 and O:2. Byproduct: Fluorosilane gas (F:6-Si:7)
            #     "reaction": "[S:1](=[O:4])(=[O:5])-[F:6].[O:2]-[Si:7](-[C:8])(-[C:9])-[C:10]>>[S:1](=[O:4])(=[O:5])-[O:2].[F:6]-[Si:7](-[C:8])(-[C:9])-[C:10]"
            # },

            # --- SCHIFF BASE (AZOMETHINES) ---
            # "Di-Amine and Di-Aldehyde Polycondensation (Schiff Base)": {
            #     "same_reactants": False,
            #     "reactant_1": "di_amine",
            #     "reactant_2": "di_aldehyde",
            #     "product": "polyazomethine_chain",
            #     "delete_atom": True,
            #     # Initiators: N:1 and C:2. Byproduct: H2O (H:5 + O:3 + H:6)
            #     "reaction": "[N:1](-[H:5])-[H:6].[CX3H1:2]=[O:3]>>[N:1]=[C:2].[H:5]-[O:3]-[H:6]"
            # },

            # --- POLYTHIOESTERS ---
            # "Di-Thiol and Di-Acid Halide Polycondensation": {
            #     "same_reactants": False,
            #     "reactant_1": "di_thiol",
            #     "reactant_2": "di_carboxylic_acid_halide",
            #     "product": "polythioester_chain",
            #     "delete_atom": True,
            #     # Initiators: S:1 and C:2. Byproduct: HX (X:5 + H:3)
            #     "reaction": "[S:1]-[H:3].[CX3:2](=[O:4])[Cl,Br,I:5]>>[S:1]-[CX3:2](=[O:4]).[Cl,Br,I:5]-[H:3]"
            # },

            # --- ADVANCED RING OPENING ---
            # "N-Carboxyanhydride (NCA) Ring-Opening Polyamidation": {
            #     "same_reactants": False,
            #     "reactant_1": "di_amine", # Acts as initiator
            #     "reactant_2": "nca",
            #     "product": "polypeptide_chain",
            #     "delete_atom": True,
            #     # Initiators N:1 and C:2. Ring opens, losing CO2 gas (O:5=C:6=O:7)
            #     "reaction": "[N:1]-[H:3].[NX3:8]1-[CX4:9]-[CX3:2](=[O:4])-[O:5]-[CX3:6](=[O:7])-1>>[N:1]-[CX3:2](=[O:4])-[CX4:9]-[N:8]-[H:3].[O:5]=[C:6]=[O:7]"
            # },

            # --- CONJUGATED POLYMERS ---
            # "Thiophene Oxidative Coupling (PEDOT synthesis equivalent)": {
            #     "same_reactants": True,
            #     "reactant_1": "thiophene",
            #     "product": "polythiophene_chain",
            #     "delete_atom": True,
            #     # Initiators: c:1 and c:2. Byproduct: H2 gas (H:3 + H:4)
            #     "reaction": "[c:1]-[H:3].[c:2]-[H:4]>>[c:1]-[c:2].[H:3]-[H:4]"
            # }




if __name__ == "__main__":
    import json
    from rdkit import Chem
    from rdkit.Chem import AllChem

    # Initialize your library
    reaction_library = ReactionLibrary()
    
    print("Starting SMARTS validation test...\n" + "="*50)

    length = 0

    for reaction_name, reaction_data in reaction_library.reactions.items():

        # Safely get the reaction SMARTS (returns None if the key doesn't exist or is None)
        rxn_smarts = reaction_data.get("reaction")
        
        print(f"Testing: {reaction_name}")
        
        # 1. Check if SMARTS exists
        if rxn_smarts is None:
            raise ValueError(f" No SMARTS string defined for reaction: {reaction_name}")
        else:
            print("SMARTS string found.")
            print("-" * 50)
            
        # 2. Attempt to parse the SMARTS
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
        
        # 3. Validate parsing
        if rxn is None:
            raise ValueError(f"Failed to parse reaction SMARTS!\nSMARTS: {rxn_smarts}")
        else:
            print("Successfully parsed reaction SMARTS.")
            # Optional: Print metadata if you want to see it
            # print(f"Metadata: {json.dumps(reaction_data, indent=2)}")
            
        print("-" * 50)
        length += 1
        
    print(f"All {length} defined SMARTS strings parsed successfully!")