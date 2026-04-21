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
            
            # ==========================================
            # ACIDS & ACID HALIDES
            # ==========================================
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

            # ==========================================
            # ALCOHOLS & PHENOLS
            # ==========================================
            "diol_monomer": {
                "functionality_type": "di_identical",
                "smarts_1": "[O,S;X2;H1;!$([O,S]C=*):1]",
                "group_name": "diol",
                "comments": None,
            },

            # ==========================================
            # AMINES & AMINO ACIDS
            # ==========================================
            "di_amine_monomer": { 
                "functionality_type": "di_identical",
                "smarts_1": "[N;H2,H1:1]",
                "group_name": "di_amine",
                "comments": None,
            },
            "amino_acid_monomer": {
                "functionality_type": "di_different",
                "smarts_1": "[NX3;H2,H1;!$(OC=*):1]",
                "smarts_2": "[CX3:2](=[O])[OX2H1]",
                "group_name": "amino_acid",
                "comments": None,
            },

            # ==========================================
            # EPOXIDES
            # ==========================================
            "di_epoxide_monomer": {
                "functionality_type": "di_identical",
                "smarts_1": "[C:1]1[O:2][C:3]1",
                "group_name": "di_epoxide"
            }

            # COMMENTED OUT FUNCTIONAL GROUPS FOR FUTURE ADDITION (NEED TO VALIDATE SMARTS FIRST)
            
            # --- CYCLICS & ISOCYANATES ---
            # "di_cyclic_anhydride_monomer": {
            #     "functionality_type": "di_identical",
            #     "smarts_1": "[CX3,c;R:1](=[OX1])[OX2,o;R][CX3,c;R:2](=[OX1])",
            #     "group_name": "di_cyclic_anhydride"
            # },
            # "di_isocyanate_monomer": {
            #     "functionality_type": "di_identical",
            #     "smarts_1": "[NX2:1]=[CX2:2]=[OX1,SX1:3]",
            #     "group_name": "di_isocyanate"
            # },

            # --- SPECIFIC AMINES ---
            # "primary_di_amine_monomer": {
            #     "functionality_type": "di_identical",
            #     "smarts_1": "[C,c:2][NX3;H2;!$(N[C,S]=*):1]",
            #     "group_name": "di_primary_amine"
            # },

            # --- VINYLS & OLEFINS ---
            # "vinyl_monomer": {
            #     "functionality_type": "vinyl",
            #     "smarts_1": "[CH2:1]=[CH;!R:2]",
            #     "group_name": "vinyl"
            # },
            # "cyclic_olefin_monomer": {
            #     "functionality_type": "vinyl",
            #     "smarts_1": "[CX3;R:1]=[CX3;R:2]",
            #     "group_name": "cyclic_olefin"
            # },

            # --- RING OPENING MONOMERS ---
            # "lactone_monomer": {
            #     "functionality_type": "mono",
            #     "smarts_1": "[CX3;R:1](=[OX1])[OX2;R:2]",
            #     "group_name": "lactone"
            # },
            # "lactam_monomer": {
            #     "functionality_type": "mono",
            #     "smarts_1": "[CX3;R:1](=[OX1])[NX3;R:2]",
            #     "group_name": "lactam"
            # },
            # "cyclic_anhydride_monomer": {
            #     "functionality_type": "mono",
            #     "smarts_1": "[C,c;R:1][CX3,c;R](=[OX1])[OX2,o;R][CX3,c;R](=[OX1])[C,c;R:2]",
            #     "group_name": "cyclic_anhydride"
            # },
            # "epoxide_monomer": {
            #     "functionality_type": "mono",
            #     "smarts_1": "[CX4;R:1]1[OX2;R:2][CX4;R:3]1",
            #     "group_name": "epoxide"
            # },
            
            # need to add more functional groups here from "J. Chem. Inf. Model. 2023, 63, 5539−5548"
        }
            # ADVANCED POLYMERIZATION MONOMERS
            
            # --- CLICK CHEMISTRY & SUFEX ---
            # "di_azide_monomer": {
            #     "functionality_type": "di_identical",
            #     "smarts_1": "[N:1]=[N+:2]=[N-:3]", 
            #     "group_name": "di_azide"
            # },
            # "di_alkyne_monomer": {
            #     "functionality_type": "di_identical",
            #     "smarts_1": "[CX2:1]#[CX2:2]",
            #     "group_name": "di_alkyne"
            # },
            # "bis_sulfonyl_fluoride_monomer": {
            #     "functionality_type": "di_identical",
            #     "smarts_1": "[SX4](=[OX1])(=[OX1])[F:1]",
            #     "group_name": "bis_sulfonyl_fluoride"
            # },
            # "bis_silyl_ether_monomer": {
            #     "functionality_type": "di_identical",
            #     "smarts_1": "[O:1]-[Si](-[C])(-[C])-[C]",
            #     "group_name": "bis_silyl_ether"
            # },

            # --- SCHIFF BASE & THIOESTERS ---
            # "di_aldehyde_monomer": {
            #     "functionality_type": "di_identical",
            #     "smarts_1": "[CX3H1:1]=[O:2]",
            #     "group_name": "di_aldehyde"
            # },
            # "di_thiol_monomer": {
            #     "functionality_type": "di_identical",
            #     "smarts_1": "[S;X2;H1;!$(SC=*):1]",
            #     "group_name": "di_thiol"
            # },

            # --- ADVANCED RING OPENING & CONJUGATED ---
            # "nca_monomer": {
            #     "functionality_type": "mono",
            #     "smarts_1": "[NX3;H1:1]1-[CX4:2]-[CX3:3](=[O:4])-[O:5]-[CX3:6](=[O:7])-1",
            #     "group_name": "nca",
            #     "comments": "N-Carboxyanhydride for synthetic polypeptide synthesis"
            # },
            # "thiophene_monomer": {
            #     "functionality_type": "di_identical",
            #     "smarts_1": "[c:1]1[s][c][c][c:2]1",
            #     "group_name": "thiophene"
            # },
        # is there monomers with both COCl and COOH groups?

if __name__ == "__main__":
    from rdkit import Chem
    
    # Initialize the functional groups library
    fg_library = FunctionalGroupsLibrary()
    
    print("Starting Functional Group SMARTS validation test...\n" + "="*50)

    length = 0

    for monomer_name, monomer_data in fg_library.monomer_types.items():
        print(f"Testing: {monomer_name}")
        
        # 1. Test smarts_1 (Required for all functional groups)
        smarts_1 = monomer_data.get("smarts_1")
        if smarts_1 is None:
            raise ValueError(f"No 'smarts_1' string defined for monomer: {monomer_name}")
            
        mol_1 = Chem.MolFromSmarts(smarts_1)
        if mol_1 is None:
            raise ValueError(f"Failed to parse 'smarts_1'!\nSMARTS: {smarts_1}")
        
        # 2. Test smarts_2 (Optional: Only required for 'di_different' types)
        smarts_2 = monomer_data.get("smarts_2")
        if smarts_2 is not None:
            mol_2 = Chem.MolFromSmarts(smarts_2)
            if mol_2 is None:
                raise ValueError(f"Failed to parse 'smarts_2'!\nSMARTS: {smarts_2}")
            print("Successfully parsed both smarts_1 and smarts_2.")
        else:
            print("Successfully parsed smarts_1.")
            
        print("-" * 50)
        length += 1
        
    print(f"All {length} defined functional group SMARTS strings parsed successfully!")