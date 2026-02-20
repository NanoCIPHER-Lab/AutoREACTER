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

from dataclasses import dataclass
from rdkit import Chem
import logging

logger = logging.getLogger(__name__)  # Module-level logger for future diagnostics.

class InputError(Exception):
    """Base class for all input-related errors."""

class InputSchemaError(InputError):
    """Missing keys or wrong top-level structure."""

class InputConflictError(InputError):
    """Mutually exclusive/invalid combinations of input fields."""

class NumericFieldError(InputError):
    """Invalid numeric values (density, temperature, counts)."""

class SmilesValidationError(InputError):
    """Invalid SMILES or RDKit parsing/canonicalization failure."""

class DuplicateMonomerError(InputError):
    """Duplicate monomer definitions detected."""

@dataclass(slots=True)
class Inputs:
    """Container that stores normalized and validated simulation inputs."""
    simulation_name: str
    temperature: list[float]  # Always normalized to a list internally for consistency.
    density: float
    monomers: dict[int, str]  # Maps monomer_id -> canonical SMILES string.
    number_of_monomers: dict[int, int]  # Maps monomer_id -> count.

class InputParser:
    """
    Validate, normalize, and prepare raw user inputs for the main simulation pipeline.
    This class is responsible for enforcing required keys, value contract checks, and
    any canonicalization steps (e.g., SMILES normalization).
    """

    def validate_inputs(self, inputs: dict) -> Inputs:
        """
        Validate the provided inputs dictionary and return a structured Inputs container.
        
        :param inputs: Raw dictionary from the user.
        :return: Inputs dataclass with normalized fields for downstream use.
        :raises KeyError: If a required field is missing.
        :raises ValueError: If any numeric value is invalid or required structures are malformed.
        """
        self.component_check(inputs)  # Ensure that the required structural keys exist.
        self._validate_numeric_fields(inputs)  # Guard numeric constraints prior to normalization.
        self.validate_smiles_rdkit(inputs)  # Normalize SMILES strings via RDKit.
        self.validate_no_duplicate_smiles(inputs)  # Prevent duplicate monomer definitions.
        logger.debug("All input validations passed successfully.")
        print("All input validations passed successfully.")  # Feedback for CLI usage.
        inputs = Inputs(
            simulation_name=inputs["simulation_name"],
            temperature=inputs["temperature"] if isinstance(inputs["temperature"], list) else [inputs["temperature"]],
            density=inputs["density"],
            monomers=inputs["monomers"],
            number_of_monomers=inputs["number_of_monomers"],
        )
        return inputs

    def component_check(self, inputs: dict) -> None:
        """
        Ensure that required high-level keys exist and that monomer dictionaries match in length.
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

    def _validate_smiles(self, smiles: str) -> str:
        """
        Validate a SMILES string via RDKit and return the canonicalized representation.

        :param smiles: Raw SMILES string input.
        :return: Canonical SMILES string consistent across the pipeline.
        :raises ValueError: If RDKit cannot parse the string.
        """
        s = (smiles or "").strip()
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            raise SmilesValidationError(f"Invalid SMILES string: {smiles!r}. RDKit failed to parse it.")
        # Canonicalize to remove whitespace/ordering differences for downstream comparison.
        return Chem.MolToSmiles(mol, canonical=True)

    def _validate_numeric_fields(self, inputs: dict) -> None:
        """
        Ensure density, temperature, and monomer counts meet the expected numeric requirements.
        """
        density = inputs.get("density", None)
        if not isinstance(density, (int, float)) or density <= 0:
            raise NumericFieldError(f"'density' must be a positive number. Got: {density!r}")

        temps = inputs.get("temperature", None)
        if temps is None:
            raise NumericFieldError("'temperature' is required (number or list of numbers).")

        temps_list = temps if isinstance(temps, list) else [temps]
        for t in temps_list:
            if not isinstance(t, (int, float)) or t <= 0:
                raise NumericFieldError(f"Temperature values must be positive numbers. Got: {t!r}")

        num_monomers = inputs.get("number_of_monomers", None)
        if not isinstance(num_monomers, dict) or not num_monomers:
            raise NumericFieldError("'number_of_monomers' must be a non-empty dict of monomer_id -> positive int.")
        for monomer_id, count in num_monomers.items():
            # Guard against bool values because bool is a subclass of int.
            if isinstance(count, bool) or not isinstance(count, int) or count <= 0:
                raise NumericFieldError(
                    f"Monomer count for {monomer_id!r} must be a positive integer. Got: {count!r}"
                )

    def validate_smiles_rdkit(self, inputs: dict) -> None:
        """
        Normalize and replace SMILES strings with their RDKit canonical equivalents.
        """
        if "monomers" not in inputs:
            raise KeyError('Missing required key: "monomers"')

        for monomer_number, smiles_string in inputs["monomers"].items():
            if not isinstance(smiles_string, str) or not smiles_string.strip():
                raise SmilesValidationError(f"Monomer {monomer_number}: SMILES must be a non-empty string.")
            canonical = self._validate_smiles(smiles_string)
            inputs["monomers"][monomer_number] = canonical  # Update in-place for downstream consistency.

    def validate_no_duplicate_smiles(self, inputs: dict) -> None:
        """
        Confirm that there are no duplicate SMILES strings between different monomer entries.
        """
        smiles_map = inputs.get("monomers", {})
        seen: dict[str, int] = {}
        for monomer_number, smi in smiles_map.items():
            normalized = (smi or "").strip()
            if normalized in seen:
                raise DuplicateMonomerError(
                    f"Duplicate SMILES detected for monomers {seen[normalized]} and {monomer_number}: {normalized!r}"
                )
            seen[normalized] = monomer_number

if __name__ == "__main__":
    # Sample input data for quick manual verification of the parser.
    inputs = {
        "simulation_name": "MySimulation",
        "temperature": [300, 400, 500],
        "density": 0.8,
        "monomers": {
            1: "ClC(=O)c1cc(cc(c1)C(Cl)=O)C(Cl)=O",  # Example monomer 1 - Trimesoyl chloride (TMC)
            2: "C1=CC(=CC(=C1)N)N",          # Example monomer 2 - m-Phenylenediamine (MPD)
            3: "CCO",   # Example monomer 3 - Ethanol
        },  
        "number_of_monomers": {
            1: 1000,
            2: 1000,
            3: 500,
        },  # or
        "stoichiometric_ratio": {
            1: 1,
            2: 1,
        },
        "number_of_total_atoms": [10000, 100000],
    }
    parser = InputParser()
    print(parser.validate_inputs(inputs))
