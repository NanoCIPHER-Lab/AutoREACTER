from __future__ import annotations
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


"""
NOTE: Use this method to store monomers in a consistent way. 
This is just an example of how to create multiple instances of a class in a loop and store them in a list. 
You can adapt this pattern to your specific use case.

class MyClass:
    def __init__(self, name):
        self.name = name
    
    def __repr__(self):
        return f"<{self.name}>"

# Create 5 instances in a loop and store them in a list
instance_list = []
for i in range(1, 6):
    instance_list.append(MyClass(f"name_{i}"))

# Accessing the objects
print(instance_list)
print(instance_list[0].name)
"""



from dataclasses import dataclass, asdict
from typing import Any
from rdkit import Chem
from rdkit.Chem import Descriptors
import logging
from typing import Literal

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

StatusType = Literal["active", "filtered", "consumed"]
CompositionMethodType = Literal["counts", "stoichiometric_ratio"]  # Placeholder for future support.

@dataclass(slots=True, frozen=True)
class MonomerEntry:
    """Canonical internal monomer representation."""
    id: int # Unique identifier for the monomer (e.g., 1, 2, 3)
    data_id: str # Original identifier from input (e.g., data_1, data_2, etc.)
    name: str | None # Optional human-readable name (not required, but can be helpful for debugging or future features) or data_1, data_2, etc. 
    smiles: str  # RDKit-canonical SMILES
    count: dict | None  # None only if stoichiometric mode 
    ratio: float | None  # None only if counts mode 
    atom_count: int # Derived from RDKit 
    molar_mass: float # Derived from RDKit
    status: StatusType = "active" # Default to "active" for all monomers at initialization


@dataclass(slots=True)
class SimulationSetup:
    """Container that stores normalized and validated simulation inputs."""
    simulation_name: str # A simple name for the simulation, used for output organization and logging.
    temperature: list[float]  # Always normalized to a list internally for consistency.
    density: float # Overall target density for the system (g/cm^3) will be calculated later based on monomer counts and molar masses. 
    monomers: list[MonomerEntry]  # Canonical internal representation (coupled smiles + count)
    composition_method: CompositionMethodType | None = None  # "counts" or "stoichiometric_ratio"
    composition : dict[str, Any] | None = None # Placeholder for future support: raw composition dict from inputs for flexible downstream use (e.g., {"method": "counts", "targets": [{"tag": "10k"}, {"tag": "100k"}]})
    stoichiometric_ratio: dict[int, float] | None = None # Placeholder for future support: monomer_id -> ratio (e.g., {1: 1, 2: 1} for a 1:1 ratio of monomer 1 to monomer 2)
    number_of_total_atoms: list[int] | None = None  # Placeholder for future support: total atom counts per monomer
    box_estimates: float | None = None  # Placeholder for future box size estimation based on monomer counts + density. will be 1/4 of the target box density length.


class InputParser:
    """
    Validate, normalize, and prepare raw user inputs for the main simulation pipeline.
    This class is responsible for enforcing required keys, value contract checks, and
    any canonicalization steps (e.g., SMILES normalization).
    """
    def validate_inputs(self, inputs: dict) -> SimulationSetup:
        self.validate_basic_format(inputs)  # Check top-level structure and determine composition method.
        simulation_name = inputs["simulation_name"]
        temperature = self._validate_temperature(inputs["temperature"])
        density = self._validate_density(inputs["density"])
        composition_method = self._get_inputs_mode(inputs["composition"])
        composition_dict = self._validate_composition(inputs["composition"], composition_method)
        monomers = self._validate_monomer_entry(inputs, composition_method, composition_dict)


        return SimulationSetup(
            simulation_name=simulation_name,
            temperature=temperature,
            density=density,
            monomers=monomers,
            composition_method=composition_method,
            composition=composition_dict
        )

    def validate_basic_format(self, inputs: dict) -> None:
        if not isinstance(inputs, dict):
            raise InputSchemaError(f"Expected input to be a dictionary. Got {type(inputs).__name__} instead.")
        required_keys = ["simulation_name", "temperature", "density", "composition", "monomers"]
        for key in required_keys:
            if key not in inputs:
                raise InputSchemaError(f"Missing required key: {key!r} in inputs dictionary.")
        
        return None

    def _get_inputs_mode(self, composition_dict: dict) -> CompositionMethodType:
        if not isinstance(composition_dict, dict):
            raise InputSchemaError(
                f"'composition' must be a dictionary. Got {type(composition_dict).__name__} instead."
            )
        method = composition_dict.get("method")

        allowed: set[CompositionMethodType] = {"counts", "stoichiometric_ratio"}

        if method not in allowed:
            raise InputSchemaError(f"Unsupported composition method: {method!r}")

        return method
        
    def _validate_temperature(self, temp: Any) -> list[float]:
        if isinstance(temp, (int, float)):
            temp_list = [temp]
        elif isinstance(temp, list) and all(isinstance(t, (int, float)) for t in temp):
            temp_list = temp
        else:
            raise NumericFieldError(f"'temperature' must be a number or a list of numbers. Got: {temp!r}")
        for t in temp_list:
            if t <= 0:
                raise NumericFieldError(f"Temperature values must be positive. Got: {t!r}")
        return temp_list

    def _validate_density(self, density: Any) -> float:
        if not isinstance(density, (int, float)) or density <= 0:
            raise NumericFieldError(f"'density' must be a positive number. Got: {density!r}")
        return density
    
    def _validate_composition(
        self,
        composition_dict: dict,
        method: CompositionMethodType
    ) -> dict:

        targets = composition_dict.get("targets")

        if not isinstance(targets, list) or len(targets) == 0:
            raise InputSchemaError(
                "'composition' must include a non-empty 'targets' list."
            )

        seen_tags: set[str] = set()

        for target in targets:
            if not isinstance(target, dict):
                raise InputSchemaError(
                    f"Each target must be a dictionary. Got: {target!r}"
                )

            tag = target.get("tag")
            if not isinstance(tag, str) or not tag.strip():
                raise InputSchemaError(
                    f"Each target must include a non-empty string 'tag'. Got: {target!r}"
                )

            if tag in seen_tags:
                raise InputSchemaError(
                    f"Duplicate target tag detected: {tag!r}"
                )

            seen_tags.add(tag)

            # Stoichiometric mode requires total_atoms
            if method == "stoichiometric_ratio":
                total_atoms = target.get("total_atoms")
                if not isinstance(total_atoms, int) or total_atoms <= 0:
                    raise NumericFieldError(
                        f"'total_atoms' must be a positive integer for stoichiometric mode. Got: {total_atoms!r}"
                    )

            # Counts mode must NOT include total_atoms (optional strictness)
            if method == "counts":
                if "total_atoms" in target:
                    raise InputSchemaError(
                        "'total_atoms' should not be provided in 'counts' mode."
                    )

        return composition_dict
    
    
    def _validate_monomer_entry(self, 
                                inputs: dict, 
                                method: CompositionMethodType,
                                composition_dict: dict
                                ) -> list[MonomerEntry]:
        """
        Validate that monomer entries are well-formed and consistent with expected structure.
        """
        validated_monomers: list[MonomerEntry] = []
        seen_monomer_list: list = []

        monomers = inputs.get("monomers")
        if not isinstance(monomers, list):
            if isinstance(monomers, dict):
                monomers = [monomers]  # Allow single monomer dict for convenience, but normalize to list internally.
            else:
                raise ValueError(f"'monomers' must be a list of monomer definitions. Got: {type(monomers).__name__} instead.")
        
        for monomer_id, monomer_dict in enumerate(monomers, start=1):
            id = monomer_id
            if "name" in monomer_dict:
                name = monomer_dict["name"]
                if not isinstance(name, str) or not name.strip():
                    name = str("data_" + str(id))
            else:
                name = str("data_" + str(id))
            smiles = monomer_dict.get("smiles", None)
            if not isinstance(smiles, str) or not smiles.strip():
                raise SmilesValidationError(f"Monomer {id}: 'smiles' must be a non-empty string. Got: {smiles!r}")
            else:
                smiles, mol = self._validate_smiles(smiles)  # Validate and canonicalize SMILES immediately.
                seen_monomer_list = self.validate_no_duplicate_smiles(smiles, seen_monomer_list=seen_monomer_list)  # Check for duplicates against previously validated monomers.

            if method == "counts":
                count_info = monomer_dict.get("count")

                if not isinstance(count_info, dict):
                    raise InputSchemaError(
                        f"Monomer {monomer_id}: 'count' must be a dictionary mapping target tags to counts."
                    )

                # Extract valid tags from composition
                valid_tags = {t["tag"] for t in composition_dict["targets"]}

                count_tags = set(count_info.keys())

                # Check missing tags
                missing = valid_tags - count_tags
                if missing:
                    raise InputSchemaError(
                        f"Monomer {monomer_id}: missing count entries for targets: {missing}"
                    )

                # Check extra tags
                extra = count_tags - valid_tags
                if extra:
                    raise InputSchemaError(
                        f"Monomer {monomer_id}: unknown target tags in count: {extra}"
                    )

                # Validate count values
                for tag, value in count_info.items():
                    if not isinstance(value, int) or value <= 0:
                        raise NumericFieldError(
                            f"Monomer {monomer_id}, target '{tag}': count must be a positive integer. Got: {value!r}"
                        )

                count = count_info
                ratio = None
            
            elif method == "stoichiometric_ratio":
                ratio = monomer_dict.get("ratio", None)
                if not isinstance(ratio, (int, float)) or ratio <= 0:
                    raise NumericFieldError(
                        f"Monomer {monomer_id}: 'ratio' must be a positive number in stoichiometric mode. Got: {ratio!r}"
                    )
                count = None
            
            else:
                raise InputSchemaError(f"Unsupported composition method: {method!r}")
            
            # Derive atom count and molar mass from RDKit
            mw = Descriptors.MolWt(mol)
            atom_count = mol.GetNumAtoms()
            validated_monomers.append(MonomerEntry(
                id=id,
                data_id=str("data_" + str(id)),
                name=name,
                smiles=smiles,
                molar_mass=mw,
                atom_count=atom_count,
                count=count,
                ratio=ratio
            ))
        return validated_monomers
    
    def _int_to_dict(self, integer: int) -> dict:
        """
        Convert an integer to a dictionary format for counts mode.
        This is a placeholder for future support of stoichiometric mode where counts may be derived from ratios.
        """
        return {"_": integer}

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
        smiles = Chem.MolToSmiles(mol, canonical=True)
        return smiles, mol

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

    def validate_no_duplicate_smiles(self, current_monomer: str, seen_monomer_list: list) -> list:
        """
        Confirm that there are no duplicate SMILES strings between different monomer entries.
        """
        for monomer in seen_monomer_list:
            if current_monomer == monomer:
                raise DuplicateMonomerError(
                    f"Duplicate monomer detected with SMILES: {current_monomer!r}. Each monomer must have a unique SMILES string."
                )
        seen_monomer_list.append(current_monomer)
        return seen_monomer_list

if __name__ == "__main__":
    # Sample input data for quick manual verification of the parser.
    inputs = {
        "simulation_name": "Example_Count_Mode",
        "temperature": [300, 400, 500],
        "density": 0.8,

        "composition": {
            "method": "counts",
            "targets": [
            {"tag": "10k"},
            {"tag": "100k"}
            ]
        },
        
        "monomers": [ 
            {
            "name": "tmc",
            "smiles": "ClC(=O)c1cc(cc(c1)C(Cl)=O)C(Cl)=O",
            "count": {
                "10k": 220,
                "100k": 2200
            }
            },
            {
            "name": "mpd",
            "smiles": "C1=CC(=CC(=C1)N)N",
            "count": {
                "10k": 220,
                "100k": 2200
            }
            },
            {
            "name": "ethanol",
            "smiles": "CCO",
            "count": {
                "10k": 110,
                "100k": 1100
            }
            }
        ]
    }
    inputs_stoichiometric =   {
        "simulation_name": "Example_Stoichiometric_Mode",
        "temperature": [300, 400, 500],
        "density": 0.8,

        "composition": {
            "method": "stoichiometric_ratio",
            "targets": [
            {"tag": "10k", "total_atoms": 10000},
            {"tag": "100k", "total_atoms": 100000}
            ]
        },

        "monomers": [
            {
            "name": "tmc",
            "smiles": "ClC(=O)c1cc(cc(c1)C(Cl)=O)C(Cl)=O",
            "ratio": 1.0
            },
            {
            "name": "mpd",
            "smiles": "C1=CC(=CC(=C1)N)N",
            "ratio": 1.0
            },
            {
            "name": "ethanol",
            "smiles": "CCO",
            "ratio": 0.5
            }
        ]
    }
    parser = InputParser()
    print(parser.validate_inputs(inputs))
    print(parser.validate_inputs(inputs_stoichiometric))
