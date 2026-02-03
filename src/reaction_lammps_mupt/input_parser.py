"""
TODO (Input Parsing Layer)
- [ ] Parse user inputs (CLI / JSON / GUI later) and validate required keys.
- [ ] Enforce: either "Number of monomers" OR "stoichiometric_ratio" must be provided (not both unless you define a rule).
- [ ] If "stoichiometric_ratio" is given, compute "Number of monomers" based on:
      - "Number of total atoms" target(s) (likely multiple) and
      - monomer atom counts / molar masses from RDKit
- [ ] Derive box size / initial density estimates from monomer counts + density.
- [ ] Output a clean, consistently-formatted dictionary to feed into main.py.
"""

from rdkit import Chem




class InputParser:
    """
    Purpose:
    - Validate and preprocess the user 'inputs' dictionary
    - Return a clean dict ready for main.py
    """

    def __init__(self, inputs_dict: dict) -> None:
        self.inputs = inputs_dict
        self.validated_inputs = self.validate_inputs(inputs_dict)

    def component_check(self, inputs: dict) -> None:
        required_keys = ["simulation_name", "temperature", "density", "monomers", "number_of_monomers"]
        for key in required_keys:
            if key not in inputs:
                raise KeyError(f"Missing required key: {key!r} in inputs dictionary.")

        if not isinstance(inputs["monomers"], dict) or len(inputs["monomers"]) == 0:
            raise ValueError("The 'monomers' key must map to a non-empty dictionary.")

        if not isinstance(inputs["number_of_monomers"], dict) or len(inputs["number_of_monomers"]) == 0:
            raise ValueError("The 'number_of_monomers' key must map to a non-empty dictionary.")

        if len(inputs["monomers"]) != len(inputs["number_of_monomers"]):
            raise ValueError(
                f"Length mismatch: 'monomers' has {len(inputs['monomers'])} entries "
                f"but 'number_of_monomers' has {len(inputs['number_of_monomers'])} entries. They must match."
            )

    def _validate_smiles(self, smiles: str) -> None:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles!r}. Terminating program.")

    def validate_smiles_rdkit(self, inputs: dict) -> None:
        if "monomers" not in inputs:
            raise KeyError('Missing required key: "monomers"')

        for monomer_number, smiles_string in inputs["monomers"].items():
            if not isinstance(smiles_string, str) or not smiles_string.strip():
                raise ValueError(f"Monomer {monomer_number}: SMILES must be a non-empty string.")
            self._validate_smiles(smiles_string)

    def validate_no_duplicate_smiles(self, inputs: dict) -> None:
        smiles_map = inputs.get("monomers", {})
        seen = {}
        for monomer_number, smi in smiles_map.items():
            if smi in seen:
                raise ValueError(
                    f"Duplicate SMILES detected for monomers {seen[smi]} and {monomer_number}: {smi!r}"
                )
            seen[smi] = monomer_number

    def validate_inputs(self, inputs: dict) -> dict:
        self.component_check(inputs)
        self.validate_smiles_rdkit(inputs)
        self.validate_no_duplicate_smiles(inputs)
        return inputs

    def to_dict(self) -> dict:
        """Return the validated inputs dict (ready for main.py)."""
        print("Validation successful. Returning cleaned inputs dictionary.")
        return self.validated_inputs


if __name__ == "__main__":
    inputs = {
        "simulation_name": "MySimulation",
        "temperature": [300, 400, 500],
        "density": 0.8,
        "monomers": {
            1: "CCO",
            2: "CCN",
        },
        "number_of_monomers": {
            1: 1000,
            2: 1000,
        },  # or
        "stoichiometric_ratio": {
            1: 1,
            2: 1,
        },
        "number_of_total_atoms": [10000, 100000],
    }
    parser = InputParser(inputs)
    print(parser.inputs)   
