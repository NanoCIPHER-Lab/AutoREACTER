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

def detect_reactions(input_dict) -> dict[str, list]:
    monomer_dict = input_dict.get("monomers", {})
    if not monomer_dict:
        raise ValueError("Input dictionary must contain a 'monomers' key with monomer data.")
    fg_results = functional_groups_detector(monomer_dict)
    reactions = reaction_selector(fg_results)
    return reactions

if __name__ == "__main__":
    sample_inputs = {
        "monomers": {
            "1": "C=CC(=O)O",
            "2": "C=CC(=O)N",
        }
    }
    detected_reactions = detect_reactions(sample_inputs)
    print("Detected Reactions:", detected_reactions)