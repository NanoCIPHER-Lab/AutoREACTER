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

"""
TODO (Input Parsing Layer)
from dataclasses import dataclass
from typing import Literal

@dataclass(slots=True, frozen=True)
class MonomerRole:
    spec: MoleculeSpec
    functionality_type: str
    fg_name: str
    fg_smarts_1: str
    fg_count_1: int
    fg_smarts_2: str | None = None
    fg_count_2: int | None = None

@dataclass(slots=True, frozen=True)
class ReactionDefinition:
    reaction_name: str
    reaction_smarts: str
    same_reactants: bool
    delete_atom: bool
    references: dict
    monomer_1: MonomerRole
    monomer_2: MonomerRole | None = None
"""

from rdkit import Chem
import logging
logger = logging.getLogger(__name__)


class InputParser:
    """
    Purpose:
    - Validate and preprocess the user 'inputs' dictionary
    - Return a clean dict ready for main.py
    """
    def validate_inputs(self, inputs: dict) -> dict:
        """
        Validate and preprocess the user 'inputs' dictionary.
        Fails with informative exceptions if any required keys are missing or if values are invalid.
        :param self: The instance of the InputParser class.
        :param inputs: The dictionary of user inputs to validate.
        :type inputs: dict
        :return: The validated and preprocessed inputs dictionary.
        :rtype: dict
        """
        self.component_check(inputs)
        self._validate_numeric_fields(inputs)
        self.validate_smiles_rdkit(inputs)
        self.validate_no_duplicate_smiles(inputs)
        return inputs
    
    
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
      
        # number_of_monomers: dict of id -> count
        num_monomers = inputs.get("number_of_monomers", None)
        if not isinstance(num_monomers, dict) or not num_monomers:
            raise ValueError("'number_of_monomers' must be a non-empty dict of monomer_id -> positive int.")
     
        for monomer_id, count in num_monomers.items():
            # bool is subclass of int, so guard it explicitly
            if isinstance(count, bool) or not isinstance(count, int) or count <= 0:
                raise ValueError(
                    f"Monomer count for {monomer_id!r} must be a positive integer. Got: {count!r}"
               )


    def validate_smiles_rdkit(self, inputs: dict) -> None:
        if "monomers" not in inputs:
            raise KeyError('Missing required key: "monomers"')

        for monomer_number, smiles_string in inputs["monomers"].items():
            if not isinstance(smiles_string, str) or not smiles_string.strip():
                raise ValueError(f"Monomer {monomer_number}: SMILES must be a non-empty string.")
            canonical = self._validate_smiles(smiles_string)
            # store cleaned version so everything downstream is consistent
            inputs["monomers"][monomer_number] = canonical


    def validate_no_duplicate_smiles(self, inputs: dict) -> None:
        smiles_map = inputs.get("monomers", {})
        seen: dict[str, int] = {}
        for monomer_number, smi in smiles_map.items():
            normalized = (smi or "").strip()
            if normalized in seen:
                raise ValueError(
                    f"Duplicate SMILES detected for monomers {seen[normalized]} and {monomer_number}: {normalized!r}"
                )
            seen[normalized] = monomer_number



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
