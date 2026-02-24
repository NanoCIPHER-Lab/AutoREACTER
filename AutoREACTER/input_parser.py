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

from dataclasses import dataclass, asdict
from typing import Any
from rdkit import Chem
from rdkit.Chem import Descriptors
import logging
from typing import Literal

# Module-level logger for future diagnostics.
logger = logging.getLogger(__name__)


class InputError(Exception):
    """Base class for all input-related errors in the parsing pipeline."""
    pass


class InputSchemaError(InputError):
    """Raised when there are missing keys or wrong top-level structure in the input."""
    pass


class InputConflictError(InputError):
    """Raised when there are mutually exclusive or invalid combinations of input fields."""
    pass


class NumericFieldError(InputError):
    """Raised when numeric values (density, temperature, counts) are invalid."""
    pass


class SmilesValidationError(InputError):
    """Raised when an invalid SMILES string is provided or RDKit parsing fails."""
    pass


class DuplicateMonomerError(InputError):
    """Raised when duplicate monomer definitions are detected."""
    pass


# Type aliases for clarity and validation
StatusType = Literal["active", "discarded"]
CompositionMethodType = Literal["counts", "stoichiometric_ratio"]  # Placeholder for future support.


@dataclass(slots=True, frozen=True)
class MonomerEntry:
    """
    Canonical internal representation of a monomer.
    
    This dataclass stores validated and processed information about a single monomer
    species, including its chemical identity (SMILES), physical properties derived 
    from RDKit, and its quantity in the simulation.
    
    Attributes:
        id: Unique integer identifier for the monomer (e.g., 1, 2, 3).
        data_id: Original identifier string from input (e.g., "data_1").
        name: Optional human-readable name for LAMMPS labeling (defaults to "data_{id}" if not provided).
        smiles: RDKit-canonical SMILES string representing the molecular structure.
        count: Dictionary mapping target tags to integer counts (used in 'counts' mode).
        ratio: Float representing the stoichiometric ratio (used in 'stoichiometric_ratio' mode).
        atom_count: Total number of atoms in the monomer, derived from RDKit.
        molar_mass: Molecular weight of the monomer, derived from RDKit.
        status: Current status of the monomer (default is "active").
    """
    id: int
    data_id: str 
    name: str | None 
    smiles: str  
    count: dict | None  # None only if stoichiometric mode 
    ratio: float | None  # None only if counts mode 
    atom_count: int 
    molar_mass: float 
    status: StatusType = "active" 
    rdkit_mol: Chem.Mol | None = None 


@dataclass(slots=True)
class SimulationSetup:
    """
    Container that stores normalized and validated simulation inputs.
    
    This class acts as the primary data transfer object (DTO) passed to the 
    simulation engine after input parsing is complete.
    Most of these data will consumed in the LAMMPS input generation step, but having them structured here allows for
    easier data management and potential future extensions.
    monomers section will be used heavyly at first stages.
    
    Attributes:
        simulation_name: A simple name for the simulation, used for output organization.
        temperature: List of temperature values (Kelvin). Normalized to list internally.
        density: Overall target density for the system (g/cm^3).
        monomers: List of MonomerEntry objects representing the system composition.
        composition_method: The method used to define composition ("counts" or "stoichiometric_ratio").
        composition: Raw composition dictionary from inputs for flexible downstream use.
        stoichiometric_ratio: Mapping of monomer IDs to ratios.
        number_of_total_atoms: List of total atom counts per target (place holder to be calculated).
        box_estimates: Estimated box size based on density (placeholder).
    """
    simulation_name: str
    temperature: list[float]  
    density: float 
    monomers: list[MonomerEntry] 
    composition_method: CompositionMethodType | None = None 
    composition : dict[str, Any] | None = None 
    stoichiometric_ratio: dict[int, float] | None = None 
    number_of_total_atoms: list[int] | None = None 
    box_estimates: float | None = None 


class InputParser:
    """
    Validates, normalizes, and prepares raw user inputs for the main simulation pipeline.
    
    This class is responsible for enforcing required keys, checking value contracts, 
    and performing canonicalization steps (e.g., SMILES normalization via RDKit).
    It transforms a raw dictionary input into a structured SimulationSetup object.
    """

    def validate_inputs(self, inputs: dict) -> SimulationSetup:
        """
        Main entry point for input validation.
        
        Orchestrates the validation of all input fields and constructs the final
        SimulationSetup object.

        Args:
            inputs: Raw dictionary input containing simulation parameters.

        Returns:
            A validated SimulationSetup object ready for the simulation engine.
        """
        # Check top-level structure and determine composition method.
        self.validate_basic_format(inputs)  
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

    def molecule_representation_of_initial_molecules(self, simulation_setup: SimulationSetup) -> list[Chem.Mol]:
        """
        Extracts RDKit Mol objects for all monomers in the simulation setup.
        
        This method can be used for downstream processing that requires direct access
        to the molecular structures, such as reaction definition or advanced analysis.

        Args:
            simulation_setup: The validated SimulationSetup object.

        Returns:
            A list of RDKit Mol objects corresponding to each monomer entry.
        """
        initial_molecules = [monomer.rdkit_mol for monomer in simulation_setup.monomers if monomer.rdkit_mol is not None]
        initial_molecules_legends = [monomer.name for monomer in simulation_setup.monomers]
        return initial_molecules, initial_molecules_legends

    def validate_basic_format(self, inputs: dict) -> None:
        """
        Validates the top-level structure of the input dictionary.
        
        Ensures the input is a dictionary and contains all required keys.

        Args:
            inputs: The raw input dictionary.

        Raises:
            InputSchemaError: If the input is not a dict or is missing required keys.
        """
        if not isinstance(inputs, dict):
            raise InputSchemaError(f"Expected input to be a dictionary. Got {type(inputs).__name__} instead.")
        
        required_keys = ["simulation_name", "temperature", "density", "composition", "monomers"]
        for key in required_keys:
            if key not in inputs:
                raise InputSchemaError(f"Missing required key: {key!r} in inputs dictionary.")
        
        return None

    def _get_inputs_mode(self, composition_dict: dict) -> CompositionMethodType:
        """
        Determines the composition method from the composition dictionary.
        
        Args:
            composition_dict: The 'composition' sub-dictionary from the inputs.

        Returns:
            The composition method string ("counts" or "stoichiometric_ratio").

        Raises:
            InputSchemaError: If the composition dict is invalid or method is unsupported.
        """
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
        """
        Validates and normalizes temperature input.
        
        Accepts a single number or a list of numbers. Ensures all values are positive.

        Args:
            temp: Temperature input (int, float, or list).

        Returns:
            A list of float temperature values.

        Raises:
            NumericFieldError: If the format is invalid or values are non-positive.
        """
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
        """
        Validates the density input.

        Args:
            density: Density value (g/cm^3).

        Returns:
            Validated float density value.

        Raises:
            NumericFieldError: If density is not a positive number.
        """
        if not isinstance(density, (int, float)) or density <= 0:
            raise NumericFieldError(f"'density' must be a positive number. Got: {density!r}")
        return density
    
    def _validate_composition(
        self,
        composition_dict: dict,
        method: CompositionMethodType
    ) -> dict:
        """
        Validates the composition dictionary structure and targets.
        
        Checks that targets exist, have valid tags, and conform to the requirements
        of the specific composition method (e.g., presence of 'total_atoms').

        Args:
            composition_dict: The composition dictionary to validate.
            method: The composition method being used.

        Returns:
            The validated composition dictionary.

        Raises:
            InputSchemaError: If structure is invalid or tags are missing/duplicated.
            NumericFieldError: If numeric fields like 'total_atoms' are invalid.
        """
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
        Validates monomer entries and constructs MonomerEntry objects.
        
        Iterates through the monomers list, validates SMILES strings, checks for 
        duplicates, validates counts/ratios based on the composition method, and 
        calculates physical properties using RDKit.

        Args:
            inputs: The raw input dictionary.
            method: The composition method ("counts" or "stoichiometric_ratio").
            composition_dict: Validated composition dictionary containing target tags.

        Returns:
            A list of validated MonomerEntry objects.

        Raises:
            ValueError: If monomers field is not a list or dict.
            SmilesValidationError: If SMILES are invalid.
            DuplicateMonomerError: If duplicate SMILES are found.
            InputSchemaError: If count/ratio structure is invalid.
            NumericFieldError: If count/ratio values are invalid.
        """
        validated_monomers: list[MonomerEntry] = []
        seen_monomer_list: list = []

        monomers = inputs.get("monomers")
        if not isinstance(monomers, list):
            if isinstance(monomers, dict):
                # Allow single monomer dict for convenience, but normalize to list internally.
                monomers = [monomers] 
            else:
                raise ValueError(f"'monomers' must be a list of monomer definitions. Got: {type(monomers).__name__} instead.")
        
        for monomer_id, monomer_dict in enumerate(monomers, start=1):
            id = monomer_id
            
            # Handle naming
            if "name" in monomer_dict:
                name = monomer_dict["name"]
                if not isinstance(name, str) or not name.strip():
                    name = str("data_" + str(id))
            else:
                name = str("data_" + str(id))
            
            # Validate SMILES
            smiles = monomer_dict.get("smiles", None)
            if not isinstance(smiles, str) or not smiles.strip():
                raise SmilesValidationError(f"Monomer {id}: 'smiles' must be a non-empty string. Got: {smiles!r}")
            else:
                # Validate and canonicalize SMILES immediately.
                smiles, mol = self._validate_smiles(smiles)  
                # Check for duplicates against previously validated monomers.
                seen_monomer_list = self.validate_no_duplicate_smiles(smiles, seen_monomer_list=seen_monomer_list)  

            # Handle Counts Mode
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
            
            # Handle Stoichiometric Ratio Mode
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
                ratio=ratio,
                rdkit_mol=mol
            ))
            
        return validated_monomers
    
    def _int_to_dict(self, integer: int) -> dict:
        """
        Convert an integer to a dictionary format for counts mode.
        
        Note: This is a placeholder for future support of stoichiometric mode 
        where counts may be derived from ratios.

        Args:
            integer: The integer value.

        Returns:
            A dictionary with a placeholder key "_".
        """
        return {"_": integer}

    def _validate_smiles(self, smiles: str) -> tuple[str, Chem.Mol]:
        """
        Validates a SMILES string via RDKit and returns the canonicalized representation.

        Args:
            smiles: Raw SMILES string input.

        Returns:
            A tuple containing:
                - Canonical SMILES string.
                - RDKit Mol object.
        
        Raises:
            SmilesValidationError: If RDKit cannot parse the string.
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
        Ensures density, temperature, and monomer counts meet expected numeric requirements.
        
        Note: This method appears to be a legacy or alternative validation routine 
        not currently hooked into the main `validate_inputs` flow.

        Args:
            inputs: Input dictionary.

        Raises:
            NumericFieldError: If any numeric constraints are violated.
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
        Confirms that there are no duplicate SMILES strings between different monomer entries.

        Args:
            current_monomer: The SMILES string of the current monomer being validated.
            seen_monomer_list: List of previously seen SMILES strings.

        Returns:
            The updated list of seen SMILES strings.
        
        Raises:
            DuplicateMonomerError: If the current SMILES is already in the seen list.
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
    
    # Example 1: Counts Mode
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
    
    # Example 2: Stoichiometric Mode
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
    print("Validating Counts Mode Input:")
    print(parser.validate_inputs(inputs))
    
    print("\nValidating Stoichiometric Mode Input:")
    print(parser.validate_inputs(inputs_stoichiometric))
