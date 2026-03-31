def smart_mapping(reactant, smarts_template, match_tuple):
        """
        Apply atom mapping from SMARTS template to reactant molecule based on match indices.
        
        This function transfers atom map numbers from a SMARTS template to corresponding
        atoms in the reactant molecule using the provided match tuple that maps SMARTS
        atom indices to reactant atom indices.
        
        Args:
            reactant (Chem.Mol): The reactant molecule to apply mapping to
            smarts_template (Chem.Mol): SMARTS pattern molecule with atom map numbers
            match_tuple (tuple): Tuple mapping SMARTS atom indices to reactant atom indices
            
        Returns:
            None: Modifies the reactant molecule in-place by setting atom map numbers
            
        Note:
            Only atoms with non-zero map numbers in the SMARTS template are processed.
            The function does nothing if match_tuple is empty.
        """
        if not match_tuple:
            return
        
        # Create mapping from SMARTS atom indices to their map numbers
        SMARTS_index_to_Map_Number = {}
        for smarts_atom in smarts_template.GetAtoms():
            map_num = smarts_atom.GetAtomMapNum()
            if map_num != 0:  # Only process atoms with explicit mapping
                SMARTS_index_to_Map_Number[smarts_atom.GetIdx()] = map_num
        
        # Apply mapping to corresponding atoms in reactant
        for smarts_pos, map_num in SMARTS_index_to_Map_Number.items():
            if smarts_pos < len(match_tuple):
                atom_index_in_mol = match_tuple[smarts_pos]
                atom = reactant.GetAtomWithIdx(atom_index_in_mol)
                atom.SetAtomMapNum(map_num)