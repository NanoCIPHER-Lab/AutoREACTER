"""
Module for detecting reactions and handling non-reactant monomers.

This module provides functionality to analyze a set of input monomers, detect
potential chemical reactions based on functional groups, and identify monomers
that do not participate in any detected reactions. It also includes an interactive
workflow to allow users to decide whether to retain non-reactant molecules in
the simulation.
"""

import json

# Attempt to import detector modules. Handles different import paths depending on
# whether the script is run as a module, part of a package, or standalone.
try:
    from functional_groups_detector import functional_groups_detector
    from reaction_detector import reaction_selector
except (ImportError, ModuleNotFoundError):
    
    from .functional_groups_detector import functional_groups_detector
    from .reaction_detector import reaction_selector



def find_non_reactant_monomers(reactions: dict, input_dict: dict, interactive: bool = True) -> list:
    """
    Identifies monomers from the input that are not participating in any detected reactions
    and prompts the user to decide whether to retain them in the simulation.

    Args:
        reactions (dict): A dictionary containing detected reaction data. Expected keys
            include 'monomer_1', 'monomer_2', and 'smiles' for each reaction entry.
        input_dict (dict): The original input dictionary containing a 'monomers' key
            mapping monomer IDs to SMILES strings.
        interactive (bool): If False, automatically retain all non-reactants without prompting.

    Returns:
        list: A list of SMILES strings representing the non-reactant monomers selected
              by the user to be retained. Returns an empty list if the user chooses
              to proceed with reactants only.
    """
    # 1) Collect all unique SMILES strings that participate in detected reactions
    reactant_smiles = set()

    for reaction_data in reactions.values():
        # Extract SMILES for monomer 1, monomer 2, and the reaction itself
        m1 = reaction_data.get("monomer_1", {}).get("smiles")
        m2 = reaction_data.get("monomer_2", {}).get("smiles")
        rxn_smiles = reaction_data.get("smiles")

        # Add valid strings to the set of reactants
        if isinstance(m1, str):
            reactant_smiles.add(m1)
        if isinstance(m2, str):
            reactant_smiles.add(m2)
        if isinstance(rxn_smiles, str):
            reactant_smiles.add(rxn_smiles)

    # 2) Build a dictionary of monomers that are not found in the reactant set
    monomers = input_dict.get("monomers", {})
    non_reactants = {}

    # Re-index non-reactants starting from ID 1
    new_id = 1
    for _, smi in monomers.items():
        if smi not in reactant_smiles:
            non_reactants[new_id] = smi
            new_id += 1

    # 3) Handle user interaction if non-reactant monomers are found
    if non_reactants:
        print(
            "\nThere are non-reactant monomers/molecules in the input.\n"
            "There may be reactions possible with the given monomers in the user inputs,\n"
            "but some of them are not detected for a reaction.\n"
            "You can choose to retain some or all of these non-reactant molecules in the simulation.\n"
            "Non-reactant molecules:"
        )
        # Display the list of non-reactant molecules to the user
        for mid, molecule in non_reactants.items():
            print(f"{mid}. {molecule}")

        # Ask user if they want to exclude all non-reactants
        select = input(
            "Do you want to proceed with monomers only (no non-monomer molecules)? (y/n): "
        ).strip().lower()

        if select == "y":
            print("Proceeding with monomers only. No non-monomer molecules will be retained.")
            non_reactants_list = []  # Return empty list if user declines non-reactants
        else:
            # Ask user to specify which non-reactants to keep by ID
            selected_non_reactants = input(
                "Please specify the monomer IDs (comma-separated) you wish to retain as non-monomer molecules: "
            ).strip()

            # Parse the comma-separated input into a set of IDs
            selected_ids = {s.strip() for s in selected_non_reactants.split(",") if s.strip()}

            # Filter the non-reactants dictionary based on user selection
            # strings only (SMILES)
            non_reactants_list = [
                smi
                for mid, smi in non_reactants.items()
                if str(mid) in selected_ids
            ]

        return non_reactants_list
    
    # Return empty list if no non-reactants were found
    return []


def detect_reactions(input_dict, interactive=True) -> tuple:
    """
    Detects chemical reactions and identifies non-reactant monomers based on the provided input.

    This function serves as the main workflow controller. It validates the input,
    detects functional groups, selects appropriate reactions, and identifies any
    monomers that did not participate in the detected reactions.

    Args:
        input_dict (dict): A dictionary containing the 'monomers' key, which maps
                           monomer IDs (int/str) to SMILES strings (str).

    Returns:
        tuple: A tuple containing:
            - dict: The detected reactions.
            - list: A list of SMILES strings for non-reactant monomers to be retained.

    Raises:
        ValueError: If the input dictionary is missing the 'monomers' key or if it is empty.
    """
    # Extract monomers from input and validate presence
    monomer_dict = input_dict.get("monomers", {})
    if not monomer_dict:
        raise ValueError("Input dictionary must contain a 'monomers' key with monomer data.")
    
    # Step 1: Detect functional groups within the monomers
    fg_results = functional_groups_detector(monomer_dict)
    
    # Step 2: Select reactions based on the detected functional groups
    reactions = reaction_selector(fg_results)
    
    # Step 3: Identify and handle monomers that are not part of any detected reaction
    non_reactants_list = find_non_reactant_monomers(reactions, input_dict, interactive=interactive)
    
    return reactions, non_reactants_list


if __name__ == "__main__":
    # Example usage of the module with sample monomer data
    sample_inputs = {
        "monomers": {
            1: "O=C(O)c1cc(O)cc(C(=O)O)c1",  # Example monomer 1
            2: "O=C(O)CCCC(O)CCCO",          # Example monomer 2
            3: "CCO",                        # Example monomer 3 (Ethanol)
        }
    }
    
    # Run the detection workflow
    detected_reactions, non_monomers = detect_reactions(sample_inputs)
    
    # Output results
    print("Detected Reactions:", json.dumps(detected_reactions, indent=2))
    if non_monomers:
        print("Non-monomer molecules to retain:", non_monomers)
