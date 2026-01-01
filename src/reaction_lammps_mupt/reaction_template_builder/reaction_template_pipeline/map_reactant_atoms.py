import rdkit
from rdkit import Chem

class MapReactantAtoms:
    """Still there is no method to verify the all possible reactions 
    can be implemented correctly. So, for now, we will just process with one or two reactions.
    TODO: Add verification method for all reactions.
    If there is a multiple reaction between the same reactants the this is the place to handle that.
    example:
reaction_smarts = "[O;!$(OC=*):1]-[H:2].[CX3:3](=[O])[OX2H1:4]>>[OX2:1]-[CX3:3](=[O]).[O:4]-[H:2]"
reactant_smiles1 = "OC(F)c1cc(C(O)Cl)cc(C(O)Br)c1"
reactant_smiles2 = "O=C(O)c1c(F)c(C(=O)O)c(Br)c(C(=O)O)c1Cl"
    to handle such cases we need to check all the possible matches and select go through them one by one.
    and assign manual mapping based on that.
    There may be some modifications reaction_selector in "run_reaction_template_pipeline.py" function as well to handle such cases.
    IF its needed same reactants case should be bring here as well.
    """
    def __init__(self, reactant1, reactant2, reaction, delete_atom=False):
        self.reactant1 = reactant1
        self.reactant2 = reactant2
        self.reaction = reaction
        self.delete_atom = delete_atom
        self.intial_atom_mapping(self.reactant1, self.reactant2)
        self.mapping_dictionary = {
            "reactant1": {
                "reactant": reactant1,
                "smarts_template": self.reaction.GetReactants()[0],
                "match": reactant1.GetSubstructMatch(self.reaction.GetReactants()[0])
            },
            "reactant2": {
                "reactant": reactant2,
                "smarts_template": self.reaction.GetReactants()[1],
                "match": reactant2.GetSubstructMatch(self.reaction.GetReactants()[1])
            }
        }

    def intial_atom_mapping(self, reactant1, reactant2):
        for atom in reactant1.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + 1001)
        for atom in reactant2.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + 2001)

    def smart_mapping(self, reactant, smarts_template, match_tuple):
        if not match_tuple:
            return
        SMARTS_index_to_Map_Number = {}
        for smarts_atom in smarts_template.GetAtoms():
            map_num = smarts_atom.GetAtomMapNum()
            if map_num != 0:
                SMARTS_index_to_Map_Number[smarts_atom.GetIdx()] = map_num

        for smarts_pos, map_num in SMARTS_index_to_Map_Number.items():
            if smarts_pos < len(match_tuple):
                atom_index_in_mol = match_tuple[smarts_pos]
                atom = reactant.GetAtomWithIdx(atom_index_in_mol)
                atom.SetAtomMapNum(map_num)

    def reveal_template_map_numbers(mol : Chem.Mol) -> None:
        '''Make map numbers assigned to products of reaction visible for display'''
        for atom in mol.GetAtoms():
            if atom.HasProp('old_mapno'): # RDKit sets certain "magic" properties post-reaction, including map numbers and reactant atom indices
                map_num = atom.GetIntProp('old_mapno')
                atom.SetAtomMapNum(map_num)