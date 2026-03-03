from typing import List
from rdkit import Chem
from rdkit.Chem import rdmolops, Draw
from PIL.Image import Image
from dataclasses import replace

# AutoREACTER internal imports for simulation configuration and reaction detection
from AutoREACTER.input_parser import SimulationSetup
from AutoREACTER.input_parser import MonomerEntry
from AutoREACTER.detectors.reaction_detector import ReactionInstance


class NonReactantsDetector:
    """
    A detector class responsible for identifying monomers that do not participate
    in any detected chemical reactions within a simulation.
    
    This class provides methods to:
    - Detect non-reacting monomers by comparing reaction participants with monomer lists
    - Visualize non-reacting monomers using RDKit grid image generation
    - Handle user selection for retaining or discarding non-reacting monomers
    
    Attributes:
        None: This class primarily operates on input parameters rather than storing state.
    
    Example:
        >>> detector = NonReactantsDetector()
        >>> updated_setup = detector.non_monomer_detector(setup, reactions)
        >>> detector.non_reactant_selection(setup, non_reactants)
    """

    def _same_molecule(self, reactants_list: List[str], mol2_smi: str) -> List[str]:
        """
        Checks if a molecule (given as SMILES string) exists in the reactants list.
        
        This method performs canonical SMILES comparison to ensure that molecules
        with the same structure but different string representations are correctly
        identified as identical.
        
        Args:
            reactants_list (List[str]): A list of SMILES strings representing molecules
                                       that participate in detected reactions.
            mol2_smi (str): The SMILES string of the molecule to check for membership.
        
        Returns:
            List[str]: The updated reactants_list with the new molecule appended
                       if it was not already present. If the input SMILES is invalid
                       (cannot be parsed), the original list is returned unchanged.
        
        Raises:
            None: Invalid SMILES are silently handled by returning the original list.
        
        Note:
            Canonical SMILES are used for comparison to ensure structural equivalence
            regardless of string representation differences (e.g., branching order).
        """
        # Parse the input SMILES string into an RDKit Mol object
        mol2 = Chem.MolFromSmiles(mol2_smi)
        
        # Return original list if SMILES is invalid (None indicates parsing failure)
        if mol2 is None:
            return reactants_list  # Or handle error: invalid SMILES
            
        # Convert to canonical SMILES for consistent comparison across different representations
        can_smi2 = Chem.MolToSmiles(mol2, canonical=True)
        detected = False
        
        # Iterate through existing reactants to check for duplicate molecules
        for mol_smi in reactants_list:
            mol1 = Chem.MolFromSmiles(mol_smi)
            if mol1 is not None:
                # Compare canonical SMILES to determine if molecules are identical
                if Chem.MolToSmiles(mol1, canonical=True) == can_smi2:
                    detected = True
                    break
                    
        # Append molecule to list if it wasn't found (not a duplicate)
        if not detected:
            reactants_list.append(mol2_smi)
            
        return reactants_list
        
    def _same_molecule_for_initaials(self, reactants_list: List[str], mol2: Chem.Mol) -> bool:
        """
        Checks if an RDKit Mol object represents a molecule already in the reactants list.
        
        This method differs from _same_molecule in that it accepts an RDKit Mol object
        directly rather than a SMILES string, and it also removes explicit hydrogens
        before comparison to ensure fair structural comparison.
        
        Args:
            reactants_list (List[str]): A list of SMILES strings representing molecules
                                       that participate in detected reactions.
            mol2 (Chem.Mol): The RDKit Mol object to check for membership.
        
        Returns:
            bool: True if the molecule exists in the reactants list (structurally 
                  equivalent), False otherwise. Returns False if the input mol is None.
        
        Note:
            - Explicit hydrogens are removed from both the query molecule and list
              molecules before comparison to ensure structural equivalence.
            - This method handles the typo "initaials" in the function name (preserved
              for backward compatibility).
        """
        # Handle None input gracefully
        if mol2 is None:
            return False
        
        # Step 1: Remove Hs from query molecule and get its canonical SMILES
        # Removing hydrogens ensures fair comparison regardless of how hydrogens were added
        mol2_clean = rdmolops.RemoveHs(mol2)
        can_smi2 = Chem.MolToSmiles(mol2_clean, canonical=True)
        
        # Iterate through the reactants list for comparison
        for mol_smi in reactants_list:
            mol1 = Chem.MolFromSmiles(mol_smi)
            if mol1 is not None:
                # Step 2: Also remove Hs from input list molecules for fair comparison
                mol1 = rdmolops.RemoveHs(mol1)
                can_smi1 = Chem.MolToSmiles(mol1, canonical=True)
                
                # Check for structural equivalence
                if can_smi1 == can_smi2:
                    return True
        
        # Return False if no matching molecule was found in the list
        return False

    def non_monomer_detector(self, 
                             simulation_setup: SimulationSetup, 
                             reaction_instances: List[ReactionInstance]) -> List[MonomerEntry]: 
        """
        Detects monomers that do not participate in any of the identified reactions.
        
        This method iterates through all monomers in the simulation setup and identifies
        those that are not involved in any of the detected reaction instances. It builds
        a list of reactants from the reaction instances and compares each monomer against
        this list to determine non-reactivity.
        
        Args:
            simulation_setup (SimulationSetup): The simulation setup containing monomers
                                                and reaction configuration.
            reaction_instances (List[ReactionInstance]): A list of identified and filtered
                                                          reactions by the user.
        
        Returns:
            List[MonomerEntry]: A list of monomer entries that do not participate in any
                                detected reactions (non-reactant monomers).
        
        Example:
            >>> setup = SimulationSetup(...)
            >>> reactions = [ReactionInstance(...), ...]
            >>> detector = NonReactantsDetector()
            >>> non_reactants = detector.non_monomer_detector(setup, reactions)
            >>> # non_reactants contains monomers not found in any reaction
        """

        # Initialize lists to track reactants and non-reactants
        reactants_list = []  # Stores SMILES of all molecules participating in reactions
        non_reactants_list = []  # Stores monomer entries that don't participate in reactions
        
        # Build the reactants list from all reaction instances
        # Each reaction involves up to two monomers (monomer_1 and monomer_2)
        for reaction in reaction_instances:
            # Add the first monomer from each reaction
            mol1 = reaction.monomer_1.smiles
            reactants_list = self._same_molecule(reactants_list, mol1)
            
            # Add the second monomer if it exists (some reactions may have only one participant)
            if reaction.monomer_2:
                mol2 = reaction.monomer_2.smiles
                reactants_list = self._same_molecule(reactants_list, mol2)
        
        # Check each monomer in the simulation setup against the reactants list
        for monomer in simulation_setup.monomers:
            # Use RDKit Mol object for comparison (more reliable than SMILES strings)
            is_in_reactions = self._same_molecule_for_initaials(reactants_list, monomer.rdkit_mol)
            if not is_in_reactions:
                # Monomer is not found in any reaction - add to non-reactants list
                non_reactants_list.append(monomer)
        
        # Return the list of non-reactant monomers
        return non_reactants_list

    def non_reactants_to_visualization(self, non_reactants_list: List[MonomerEntry]) -> Image | None:
        """
        Generates a grid image visualization of non-reacting monomers.
        
        This method creates a visual representation of all non-reacting monomers
        using RDKit's MolsToGridImage function, allowing users to visually identify
        which molecules are not participating in reactions.
        
        Args:
            non_reactants_list (List[MonomerEntry]): A list of monomer entries that do not
                                                     participate in any reactions.
        
        Returns:
            PIL.Image.Image | None: A grid image object containing the molecular
                                     structures of all non-reacting monomers, or None if no non-reacting monomers are found.
        
        Note:
            The grid image displays molecules with their names as legends and arranges
            them in a grid with 3 molecules per row. Each molecule image is 400x400 pixels.
        """
        # Print details of each non-reacting monomer to console
        if not non_reactants_list:
            print("No non-reacting monomers detected.")
            return None  # Or handle as appropriate (e.g., return an empty image)
        
        for monomer in non_reactants_list:
            print(f"ID: {monomer.id}, Name: {monomer.name}, SMILES: {monomer.smiles}")
        
        # Generate grid image using RDKit
        # - Extract RDKit Mol objects from monomer entries
        # - Use monomer names as legend text
        # - Configure grid layout: 3 molecules per row, 400x400 pixel images
        img = Draw.MolsToGridImage([monomer.rdkit_mol for monomer in non_reactants_list],
                                   legends=[monomer.name for monomer in non_reactants_list],
                                   molsPerRow=3,
                                   subImgSize=(400,400))
        return img
    
    def non_reactant_selection(self, simulation_setup: SimulationSetup, non_reactants_list: List[MonomerEntry]) -> SimulationSetup:
        """
        Handles user selection of non-reactant monomers to retain or discard in the simulation.
        
        This method presents the user with options to either keep or remove non-reacting
        monomers from the simulation. It supports three selection modes:
        - N (No): Discard all non-reactant monomers
        - A (All): Retain all non-reactant monomers
        - S (Select): Allow selective retention based on monomer IDs
        
        The method updates the 'status' attribute of selected monomers to 'discarded'
        for those that should be removed from the simulation.
        
        Args:
            simulation_setup (SimulationSetup): The simulation setup object whose monomers
                                                will be potentially modified.
            non_reactants_list (List[MonomerEntry]): A list of monomer entries identified
                                                     as non-reacting.
        
        Returns:
            SimulationSetup: The updated simulation setup with non-reactant monomers
                              either retained or marked as discarded based on user input.
        
        Note:
            - This method replaces MonomerEntry objects in simulation_setup.monomers with
              updated copies (via dataclasses.replace) to reflect 'discarded' status.
            - For single non-reactant monomers, only N and A options are presented.
            - The status values follow the StatusType convention: 'active' or 'discarded'.
        """
        if not non_reactants_list:
            print("No non-reacting monomers to select from.")
            return simulation_setup  # No changes needed if there are no non-reactants
        
        # Display selection guide for user interaction with non-reactant monomers
        print("""
Selection Guide for Non-Reactant Monomers:
- N: Discard all non-reactant monomers from the simulation.
- A: Retain all non-reactant monomers in the simulation.
- S: Select specific non-reactant monomers to retain. You will be prompted 
to enter the IDs of the monomers you wish to keep, separated by commas.
Example: If you want to keep monomers with IDs 1 and 3, you would enter: 1,3
              """)
        # Proceed only if there are non-reactant monomers to process
        if non_reactants_list:
            # Display information about non-participating monomers
            print("The following monomers do not participate in any detected reactions:")
            for monomer in non_reactants_list:
                print(f"ID: {monomer.id}, Name: {monomer.name}, SMILES: {monomer.smiles}")
            
            # Main selection loop - continues until valid input is received
            while True:
                # Handle special case: only one non-reactant monomer
                if len(non_reactants_list) == 1:
                    print("Only one non-reactant monomer detected")
                    # Simplified prompt for single monomer scenario
                    user_input = input(f"Do you want to keep {non_reactants_list[0].name} in the simulation? (N/A): ").strip().upper()
                    if user_input not in ['N', 'A']:
                        print("Invalid input. Please enter N or A.")
                        continue
                    if user_input == 'N':
                        # Discard the single non-reactant monomer: mark it as 'discarded'
                        target_id = non_reactants_list[0].id
                        simulation_setup.monomers = [
                            replace(m, status='discarded') if m.id == target_id else m
                            for m in simulation_setup.monomers
                        ]
                        break
                    elif user_input == 'A':
                        # Keep the single monomer (retain)
                        break
                else:
                    # Multi-monomer scenario - offer all three options        
                    user_input = input("Do you want to retain these non-reactant monomers in the simulation? (N/A/S): ").strip().upper()
                    if user_input not in ['N', 'A', 'S']:
                        print("Invalid input. Please enter N, A, or S.")
                        continue

                    if user_input == 'N':
                        # Discard all non-reactants: replace each with 'discarded' status
                        discard_ids = {m.id for m in non_reactants_list}
                        simulation_setup.monomers = [
                            replace(m, status='discarded') if m.id in discard_ids else m
                            for m in simulation_setup.monomers
                        ]
                        break
                    elif user_input == 'A':
                        # Retain all non-reactants: do nothing (keep 'active' status)
                        break
                    elif user_input == 'S':
                        # Selective retention: user specifies which IDs to keep
                        # Others will be marked as 'discarded'
                        selected_ids = input("Enter the IDs of the monomers you want to retain, separated by commas: ")
                        # Parse input into list of integers ( monomer IDs)
                        selected_ids = [int(id.strip()) for id in selected_ids.split(',')]
                        discard_ids = {m.id for m in non_reactants_list if m.id not in selected_ids}
                        simulation_setup.monomers = [
                            replace(m, status='discarded') if m.id in discard_ids else m
                            for m in simulation_setup.monomers
                        ]
                        break
                    
        # Return the modified simulation setup
        return simulation_setup
