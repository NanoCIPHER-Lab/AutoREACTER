"""
TODO (Input Parsing Layer)
- [ ] Parse user inputs (CLI / JSON / GUI later) and validate required keys.
- [ ] Enforce: either "count" (per monomer) OR "stoichiometric_ratio" must be provided (not both unless you define a rule).
- [ ] If "stoichiometric_ratio" is given, compute "count" based on:
      - "Number of total atoms" target(s) (likely multiple) and
      - monomer atom counts / molar masses from RDKit
- [ ] Derive box size / initial density estimates from monomer counts + density.
- [ ] Output a clean, consistently-formatted dictionary to feed into main.py.
"""

from rdkit import Chem
import logging
logger = logging.getLogger(__name__)


class InputParser:
    """
    Purpose:
    - Validate and preprocess the user 'inputs' dictionary
    - Return a clean dict ready for main.py

    Expected monomer schema (unified format)::

        "monomers": {
            "adipic_acid": {
                "smiles": "O=C(O)CCCCC(=O)O",
                "count": 500
            },
            ...
        }
    """

    def __init__(self, inputs_dict: dict) -> None:
        self.inputs = inputs_dict
        self.validated_inputs = self.validate_inputs(inputs_dict)

    def component_check(self, inputs: dict) -> None:
        required_keys = ["simulation_name", "temperature", "density", "monomers"]
        for key in required_keys:
            if key not in inputs:
                raise KeyError(f"Missing required key: {key!r} in inputs dictionary.")

        if not isinstance(inputs["monomers"], dict) or len(inputs["monomers"]) == 0:
            raise ValueError("The 'monomers' key must map to a non-empty dictionary.")

        for name, entry in inputs["monomers"].items():
            if not isinstance(entry, dict):
                raise ValueError(
                    f"Monomer {name!r} must be a dictionary with 'smiles' and 'count' keys. "
                    f"Got: {entry!r}"
                )
            if "smiles" not in entry:
                raise KeyError(f"Monomer {name!r} is missing required key 'smiles'.")
            if "count" not in entry:
                raise KeyError(f"Monomer {name!r} is missing required key 'count'.")

    def validate_inputs(self, inputs: dict) -> dict:
        self.component_check(inputs)
        self._validate_numeric_fields(inputs)
        self.validate_smiles_rdkit(inputs)
        self.validate_no_duplicate_smiles(inputs)
        return inputs

    def _validate_smiles(self, smiles: str) -> str:
        """
        Validate a SMILES string with RDKit and return a canonical SMILES.
        Raises ValueError if invalid.
        """
        s = (smiles or "").strip()
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles!r}")
        # Canonicalize to remove whitespace/ordering differences
        return Chem.MolToSmiles(mol, canonical=True)

    def _validate_numeric_fields(self, inputs: dict) -> None:
        # Density
        density = inputs.get("density", None)
        if not isinstance(density, (int, float)) or density <= 0:
            raise ValueError(f"'density' must be a positive number. Got: {density!r}")

        # Temperature: allow scalar or list
        temps = inputs.get("temperature", None)
        if temps is None:
            raise ValueError("'temperature' is required (number or list of numbers).")

        temps_list = temps if isinstance(temps, list) else [temps]
        for t in temps_list:
            if not isinstance(t, (int, float)) or t <= 0:
                raise ValueError(f"Temperature values must be positive numbers. Got: {t!r}")

        # count: validated per-monomer entry
        for name, entry in inputs.get("monomers", {}).items():
            if not isinstance(entry, dict):
                continue  # structural error caught in component_check
            count = entry.get("count")
            # bool is subclass of int, so guard it explicitly
            if isinstance(count, bool) or not isinstance(count, int) or count <= 0:
                raise ValueError(
                    f"Monomer count for {name!r} must be a positive integer. Got: {count!r}"
                )

    def validate_smiles_rdkit(self, inputs: dict) -> None:
        if "monomers" not in inputs:
            raise KeyError('Missing required key: "monomers"')

        for name, entry in inputs["monomers"].items():
            smiles_string = entry.get("smiles", "")
            if not isinstance(smiles_string, str) or not smiles_string.strip():
                raise ValueError(f"Monomer {name!r}: SMILES must be a non-empty string.")
            canonical = self._validate_smiles(smiles_string)
            # store cleaned version so everything downstream is consistent
            inputs["monomers"][name]["smiles"] = canonical

    def validate_no_duplicate_smiles(self, inputs: dict) -> None:
        smiles_map = inputs.get("monomers", {})
        smiles_to_name: dict[str, str] = {}
        for name, entry in smiles_map.items():
            smi = entry.get("smiles", "") if isinstance(entry, dict) else ""
            normalized = (smi or "").strip()
            if normalized in smiles_to_name:
                raise ValueError(
                    f"Duplicate SMILES detected for monomers {smiles_to_name[normalized]!r} and {name!r}: {normalized!r}"
                )
            smiles_to_name[normalized] = name

    def to_dict(self) -> dict:
        """Return the validated inputs dict (ready for main.py)."""
        logger.info("Validation successful. Returning cleaned inputs dictionary.")
        return self.validated_inputs


if __name__ == "__main__":
    inputs = {
        "simulation_name": "MySimulation",
        "temperature": [300, 400, 500],
        "density": 0.8,
        "monomers": {
            "ethanol": {
                "smiles": "CCO",
                "count": 1000,
            },
            "ethylamine": {
                "smiles": "CCN",
                "count": 1000,
            },
        },
    }
    parser = InputParser(inputs)
    print(parser.inputs)
