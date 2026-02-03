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


from rdkit import Chem


class InputParser:
    """
    Purpose:
    - Validate and preprocess the user 'inputs' dictionary
    - Return a clean dict ready for main.py
    """

    def __init__(self, inputs_dict: dict) -> None:
        """
        Initialize the InputParser with a raw inputs dictionary and perform validation.
        
        Parameters:
            inputs_dict (dict): User-provided inputs to validate and normalize; stored as `self.inputs` and the validated result as `self.validated_inputs`.
        """
        self.inputs = inputs_dict
        self.validated_inputs = self.validate_inputs(inputs_dict)

    def component_check(self, inputs: dict) -> None:
        """
        Validate presence and basic structure of required top-level input fields.
        
        Parameters:
            inputs (dict): User-provided inputs dictionary expected to contain keys
                "simulation_name", "temperature", "density", "monomers", and "number_of_monomers".
                "monomers" and "number_of_monomers" must each be non-empty dictionaries with matching entry counts.
        
        Raises:
            KeyError: If any required key is missing from `inputs`.
            ValueError: If "monomers" or "number_of_monomers" is not a non-empty dictionary,
                        or if their lengths differ.
        """
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
        """
        Validate that a SMILES string can be parsed by RDKit.
        
        Parameters:
            smiles (str): SMILES string to validate.
        
        Raises:
            ValueError: If the SMILES string cannot be converted to an RDKit molecule.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles!r}. Terminating program.")

    def validate_smiles_rdkit(self, inputs: dict) -> None:
        """
        Validate that `inputs` contains a "monomers" mapping and that each monomer's SMILES is a non-empty string that RDKit can parse.
        
        Parameters:
            inputs (dict): Input dictionary expected to contain a "monomers" mapping of monomer identifiers to SMILES strings.
        
        Raises:
            KeyError: If the "monomers" key is missing from `inputs`.
            ValueError: If a monomer's SMILES is not a non-empty string.
            ValueError: If RDKit fails to parse any SMILES (propagated from `_validate_smiles`).
        """
        if "monomers" not in inputs:
            raise KeyError('Missing required key: "monomers"')

        for monomer_number, smiles_string in inputs["monomers"].items():
            if not isinstance(smiles_string, str) or not smiles_string.strip():
                raise ValueError(f"Monomer {monomer_number}: SMILES must be a non-empty string.")
            self._validate_smiles(smiles_string)

    def validate_no_duplicate_smiles(self, inputs: dict) -> None:
        """
        Ensure no two monomer entries share the same SMILES string.
        
        Parameters:
            inputs (dict): Input dictionary containing a "monomers" key whose value is a mapping of monomer identifiers to SMILES strings.
        
        Raises:
            ValueError: If the same SMILES string appears for more than one monomer, indicating which monomer identifiers conflict.
        """
        smiles_map = inputs.get("monomers", {})
        seen = {}
        for monomer_number, smi in smiles_map.items():
            if smi in seen:
                raise ValueError(
                    f"Duplicate SMILES detected for monomers {seen[smi]} and {monomer_number}: {smi!r}"
                )
            seen[smi] = monomer_number

    def validate_inputs(self, inputs: dict) -> dict:
        """
        Validate the provided inputs dictionary and return it once all checks pass.
        
        Parameters:
            inputs (dict): Raw user inputs to validate (must contain required keys and valid monomer SMILES).
        
        Returns:
            dict: The same inputs dictionary after successful validation, suitable for downstream processing.
        """
        self.component_check(inputs)
        self.validate_smiles_rdkit(inputs)
        self.validate_no_duplicate_smiles(inputs)
        return inputs

    def to_dict(self) -> dict:
        """
        Provide the validated inputs dictionary prepared for downstream processing.
        
        Returns:
            dict: Cleaned inputs dictionary containing validated and ready-to-use values for main.py.
        """
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