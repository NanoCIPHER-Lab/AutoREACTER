
from dataclasses import dataclass
from typing import Any, List, Literal, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import rdmolops, Draw

from AutoREACTER.input_parser import SimulationSetup
from AutoREACTER.input_parser import MonomerEntry
from AutoREACTER.detectors.reaction_detector import ReactionInstance
from AutoREACTER.detectors.reaction_detector  import FunctionalGroupInfo
from AutoREACTER.detectors.reaction_detector  import MonomerRole


class NonReactantsDetector:

    def _same_molecule(self, reactants_list: List[str], mol2_smi: str) -> List[str]:
        """
        Checks if a molecule (SMILES) exists in the list, 
        if not, appends it and returns the updated list.
        """
        mol2 = Chem.MolFromSmiles(mol2_smi)
        if mol2 is None:
            return reactants_list # Or handle error: invalid SMILES
            
        can_smi2 = Chem.MolToSmiles(mol2, canonical=True)
        detected = False
        
        for mol_smi in reactants_list:
            mol1 = Chem.MolFromSmiles(mol_smi)
            if mol1 is not None:
                if Chem.MolToSmiles(mol1, canonical=True) == can_smi2:
                    detected = True
                    break
                    
        if not detected:
            reactants_list.append(mol2_smi)
            
        return reactants_list
        
    def _same_molecule_for_initaials(self, reactants_list: List[str], mol2: Chem.Mol) -> bool:
        if mol2 is None:
            return False
        
        # 1. In case - Remove Hs from query molecule and get its canonical SMILES
        mol2_clean = rdmolops.RemoveHs(mol2)
        can_smi2 = Chem.MolToSmiles(mol2_clean, canonical=True)
        
        for mol_smi in reactants_list:
            mol1 = Chem.MolFromSmiles(mol_smi)
            if mol1 is not None:
                # 2. Also remove Hs from input list molecules for fair comparison
                mol1 = rdmolops.RemoveHs(mol1)
                can_smi1 = Chem.MolToSmiles(mol1, canonical=True)
                
                if can_smi1 == can_smi2:
                    return True

    def non_monomer_detector(self, 
                             simulation_setup: SimulationSetup, 
                             reaction_instances: List[ReactionInstance]) -> SimulationSetup: 
        """
        Detects monomers that do not participate in any of the identified reactions.
        Args:
            simulation_setup: The simulation setup containing monomers and reactions.
            reaction_instances: A list of identified and filterd reactions by the user
        Returns:
            SimulationSetup: The updated simulation setup with non-reactant monomers identified.
            This function iterates through the monomers in the simulation setup and checks if they are involved in any of the detected reactions. 
            If a monomer is not found in any reaction, It will be added to a non-reactants list. 
            Then ask the used if they want to retain these non-reactant monomers in the simulation. 
            Inputs will be "N" for no and "A" for all "S" for selected.
            SimulationSetup.monomers: list[MonomerEntry] . status: StatusType = "active" will be updated to "discarded" for non-reactant monomers if the user choose to discard them.
            Finally, the function returns the updated SimulationSetup with non-reactant monomers either retained or discarded based on user input.
            """
        print("""
        Selection Guide for Non-Reactant Monomers:
        - N: Discard all non-reactant monomers from the simulation.
        - A: Retain all non-reactant monomers in the simulation.
        - S: Select specific non-reactant monomers to retain. You will be prompted to enter the IDs of the monomers you wish to keep, separated by commas.""")
        reactants_list = []
        non_reactants_list = []
        for reaction in reaction_instances:
            mol1 = reaction.monomer_1.smiles
            reactants_list = self._same_molecule(reactants_list, mol1)
            if reaction.monomer_2:
                mol2 = reaction.monomer_2.smiles
                reactants_list = self._same_molecule(reactants_list, mol2)
        for monomer in simulation_setup.monomers:
            is_in_reactions = self._same_molecule_for_initaials(reactants_list, monomer.rdkit_mol)
            if not is_in_reactions:
                non_reactants_list.append(monomer)
        
        return non_reactants_list

    def non_reactants_to_visualization(self, non_reactants_list: List[MonomerEntry]) -> Draw.MolsToGridImage:
        for monomer in non_reactants_list:
            print(f"ID: {monomer.id}, Name: {monomer.name}, SMILES: {monomer.smiles}")
        
        img = Draw.MolsToGridImage([monomer.rdkit_mol for monomer in non_reactants_list],
                                   legends=[monomer.name for monomer in non_reactants_list],
                                   molsPerRow=3,
                                   subImgSize=(400,400))
        return img
    
    def non_reactant_selection(self, simulation_setup: SimulationSetup, non_reactants_list: List[MonomerEntry]) -> SimulationSetup:
        """
        Handles user selection of non-reactant monomers to retain or discard in the simulation.
        """
        if non_reactants_list:
            print("The following monomers do not participate in any detected reactions:")
            for monomer in non_reactants_list:
                print(f"ID: {monomer.id}, Name: {monomer.name}, SMILES: {monomer.smiles}")
            
            while True:
                if len(non_reactants_list) == 1:
                    print("Only one non-reactant monomer detected")
                    user_input = input(f"Do you want to keep {non_reactants_list[0].name} in the simulation? (N/A): ").strip().upper()
                    if user_input not in ['N', 'A']:
                        print("Invalid input. Please enter N or A.")
                        continue
                    if user_input == 'N':
                        non_reactants_list[0].status = 'discarded'
                        break
                    elif user_input == 'A':
                        break
                else:        
                    user_input = input("Do you want to retain these non-reactant monomers in the simulation? (N/A/S): ").strip().upper()
                    if user_input not in ['N', 'A', 'S']:
                        print("Invalid input. Please enter N, A, or S.")
                        continue

                    if user_input == 'N':
                        # Discard all non-reactants
                        for monomer in non_reactants_list:
                            monomer.status = 'discarded'
                        break
                    elif user_input == 'A':
                        # Retain all non-reactants
                        break
                    elif user_input == 'S':
                        # Selectively retain
                        selected_ids = input("Enter the IDs of the monomers you want to retain, separated by commas: ")
                        selected_ids = [int(id.strip()) for id in selected_ids.split(',')]
                        for monomer in non_reactants_list:
                            if monomer.id not in selected_ids:
                                monomer.status = 'discarded'
                        break
                    
        return simulation_setup