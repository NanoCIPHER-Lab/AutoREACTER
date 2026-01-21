"""
TODO:

If no reaction-capable molecules are detected among the provided monomers, but the reaction database indicates that reactions are theoretically possible, prompt the user to decide whether non-monomer molecules should be retained.

- If the user chooses to keep non-monomer molecules:
  - Generate the simulation scripts without automatic reaction handling.
  - Require the user to manually provide reaction templates and the corresponding LUNAR map files and input scripts.

Additionally, handle the case where reactions are theoretically possible, but none are detected for the specific set of provided monomers.

In this situation:
- First, ask the user whether to proceed with:
  1) monomers only, or  
  2) a mixture of monomers and non-monomer molecules.
- If the user selects monomers only:
  - Generate scripts without automatic reaction handling.
  - Request manual reaction templates and map files from the user.
"""
import json

try:
    from functional_groups_detector import functional_groups_detector
    from reaction_detector import reaction_selector
except ImportError or ModuleNotFoundError:
    try:
        from detectors.functional_groups_detector import functional_groups_detector
        from detectors.reaction_detector import reaction_selector
    except ImportError or ModuleNotFoundError:
        from reaction_lammps_mupt.detectors.functional_groups_detector import functional_groups_detector
        from reaction_lammps_mupt.detectors.reaction_detector import reaction_selector


def find_non_reactant_monomers(reactions: dict, input_dict: dict) -> dict:
    # 1) Collect all reactant SMILES
    reactant_smiles = set()

    for reaction_data in reactions.values():
        m1 = reaction_data.get("monomer_1", {}).get("smiles")
        m2 = reaction_data.get("monomer_2", {}).get("smiles")
        rxn_smiles = reaction_data.get("smiles")

        if isinstance(m1, str):
            reactant_smiles.add(m1)
        if isinstance(m2, str):
            reactant_smiles.add(m2)
        if isinstance(rxn_smiles, str):
            reactant_smiles.add(rxn_smiles)

    # 2) Build re-indexed non-reactant dict (start from 1)
    monomers = input_dict.get("monomers", {})
    non_reactants = {}

    new_id = 1
    for _, smi in monomers.items():
        if smi not in reactant_smiles:
            non_reactants[new_id] = smi
            new_id += 1

    if non_reactants:
        print(
            "\nThere are non-reactant monomers/molecules in the input.\n" \
            "There are may be reactions possible with the given monomers in the user inputs,\n" \
            "but some of them are not detected for a reaction.\n" \
            "You can choose to retain some or all of these non-reactant molecules in the simulation.\n"
            "Non-reactant molecules:"
        )
        for mid, molecule in non_reactants.items():
            print(f"{mid}. {molecule}")

        select = input(
            "Do you want to proceed with monomers only (no non-monomer molecules)? (y/n): "
        ).strip().lower()

        if select == "y":
            print("Proceeding with monomers only. No non-monomer molecules will be retained.")
            non_reactants_list = []   # strings only
        else:
            selected_non_reactants = input(
                "Please specify the monomer IDs (comma-separated) you wish to retain as non-monomer molecules: "
            ).strip()

            selected_ids = {s.strip() for s in selected_non_reactants.split(",") if s.strip()}

            # strings only (SMILES)
            non_reactants_list = [
                smi
                for mid, smi in non_reactants.items()
                if str(mid) in selected_ids
            ]

            
        return non_reactants_list

def detect_reactions(input_dict) -> dict[str, list]:
    monomer_dict = input_dict.get("monomers", {})
    if not monomer_dict:
        raise ValueError("Input dictionary must contain a 'monomers' key with monomer data.")
    fg_results = functional_groups_detector(monomer_dict)
    reactions = reaction_selector(fg_results)
    non_reactants_list = find_non_reactant_monomers(reactions, input_dict)
    return reactions, non_reactants_list

if __name__ == "__main__":
    sample_inputs = {
        "monomers": {
            1: "O=C(O)c1cc(O)cc(C(=O)O)c1",
            2: "O=C(O)CCCC(O)CCCO",
            3: "CCO",
        }
    }
    detected_reactions, non_monomers = detect_reactions(sample_inputs)
    print("Detected Reactions:", json.dumps(detected_reactions, indent=2))
    if non_monomers:
        print("Non-monomer molecules to retain:", non_monomers)
