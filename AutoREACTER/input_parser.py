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

from dataclasses import dataclass, asdict
from typing import Any
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
import logging
from typing import Literal
from PIL.Image import Image

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
CompositionMethodType = Literal["counts", "stoichiometric_ratio"]
ForceFieldType = Literal[
    "PCFF-IFF", "PCFF", "Compass",
    "CVFF-IFF", "CVFF", "DRIEDING", 
]
# Note: Since Clay-FF and OPLS-AA are not fully supported within LUNAR, they are consequently unsupported in AutoREACTER.

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
        rdkit_mol: The RDKit Mol object corresponding to the SMILES string.
        count: Dictionary mapping target tags to integer counts (used in 'counts' mode).
        ratio: Float representing the stoichiometric ratio (used in 'stoichiometric_ratio' mode).
        status: Current status of the monomer (default is "active").
    """
    id: int
    data_id: str 
    name: str | None 
    smiles: str  
    count: dict | None  # None only if stoichiometric mode 
    ratio: float | None  # None only if counts mode 
    rdkit_mol: Chem.Mol | None = None # Store the RDKit Mol object
    status: StatusType = "active"

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
        density: List of density values in g/cm^3.
        force_field: Optional string specifying the force field to use (e.g., "PCFF", "OPLS-AA"). Default is "PCFF".
        monomers: List of MonomerEntry objects representing the system composition.
        composition_method: The method used to define composition ("counts" or "ratio").
        composition: Raw composition dictionary from inputs for flexible downstream use.
        stoichiometric_ratio: Mapping of monomer IDs to ratios.
        number_of_total_atoms: List of total atom counts per target (place holder to be calculated).
        box_estimates: Estimated box size based on density (placeholder).
    """
    simulation_name: str
    temperature: list[float]
    density: list[float]
    force_field: str | None
    monomers: list[MonomerEntry]
    composition_method: CompositionMethodType | None = None
    composition: dict[str, Any] | None = None
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

        replicas_dict = inputs["replicas"]
        composition_method = self._get_inputs_mode(replicas_dict)

        validated_replicas = self._validate_replicas(
            replicas_dict, composition_method
        )

        monomers = self._validate_monomer_entry(
            inputs,
            composition_method,
            validated_replicas["systems"],
        )

        force_field = self._validate_force_field(
            inputs.get("force_field", None)
        )

        return SimulationSetup(
            simulation_name=simulation_name,
            temperature=validated_replicas["temperatures"],
            density=validated_replicas["density"],
            monomers=monomers,
            composition_method=composition_method,
            composition=validated_replicas,
            force_field=force_field
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
    
    def initial_molecules_image_grid(self, simulation_setup: SimulationSetup) -> Image:
        """
        Creates a grid image of the initial molecules for visualization purposes.
        Args:
            simulation_setup: The validated SimulationSetup object.

        Returns:
            A Image object containing the grid of molecule images.
        """
        initial_molecules = [monomer.rdkit_mol for monomer in simulation_setup.monomers if monomer.rdkit_mol is not None]
        initial_molecules_legends = [monomer.name for monomer in simulation_setup.monomers]
        return Draw.MolsToGridImage(initial_molecules, molsPerRow=3, subImgSize=(400, 400), legends=initial_molecules_legends)

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
            raise InputSchemaError(
            f"Expected input to be a dictionary. Got {type(inputs).__name__} instead."
            )
        
        required_keys = ["simulation_name", "replicas", "monomers"]
        for key in required_keys:
            if key not in inputs:
                raise InputSchemaError(
                    f"Missing required key: {key!r} in inputs dictionary."
                )
        
        return None

    def _get_inputs_mode(self, replicas_dict: dict) -> CompositionMethodType:
        """
        Determines the composition method from the composition dictionary.
        
        Args:
            composition_dict: The 'composition' sub-dictionary from the inputs.

        Returns:
            The composition method string ("counts" or "ratio").

        Raises:
            InputSchemaError: If the composition dict is invalid or method is unsupported.
        """
        if not isinstance(replicas_dict, dict):
            raise InputSchemaError(
                f"'replicas' must be a dictionary. Got {type(replicas_dict).__name__} instead."
            )

        method = replicas_dict.get("method")
        allowed = {"counts", "ratio"}

        if method not in allowed:
            raise InputSchemaError(
                f"Unsupported replicas method: {method!r}"
            )

        return "counts" if method == "counts" else "stoichiometric_ratio"
        
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
        if isinstance(density, (int, float)):
            density_list = [float(density)]

        elif isinstance(density, list) and all(isinstance(d, (int, float)) for d in density):
            density_list = [float(d) for d in density]

        else:
            raise NumericFieldError(
                f"'density' must be a number or list of numbers. Got: {density!r}"
            )

        for d in density_list:
            if d <= 0:
                raise NumericFieldError(
                    f"Density values must be positive. Got: {d!r}"
                )

        return density_list
    
    def _validate_force_field(self, force_field: Any) -> str:
        """
        Validates the force field input.

        Args:
            force_field: Force field name as a string.

        Returns:
            Validated force field name.

        Raises:
            InputSchemaError: If force field is not a string or is empty.
        """
        allowed: set[ForceFieldType] = {
            "PCFF-IFF", "PCFF", "Compass", "CVFF-IFF", "CVFF", "Clay-FF", "DRIEDING", "OPLS-AA"
        }
        
        if force_field is None:
            return "PCFF"  # Default value

        if force_field not in allowed:
            raise InputSchemaError(f"Unsupported force field: {force_field!r}")
        return force_field
    
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
                                systems: list[dict]
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
            systems: List of system dictionaries from the replicas section, used to validate counts/ratios.

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
        seen_smiles: list[str] = []

        monomers = inputs.get("monomers")
        if not isinstance(monomers, list):
            raise InputSchemaError("'monomers' must be a list.")

        # Collect all tags
        system_tags = [system["tag"] for system in systems]

        for monomer_id, monomer_dict in enumerate(monomers, start=1):

            name = monomer_dict.get("name")

            if not isinstance(name, str) or not name.strip():
                name = f"data_{monomer_id}" # Default name if not provided or invalid
                
            smiles_raw = monomer_dict.get("smiles")

            smiles, mol = self._validate_smiles(smiles_raw)
            seen_smiles = self.validate_no_duplicate_smiles(
                smiles, seen_smiles
            )

            if method == "counts":

                count_map = {}

                for system in systems:
                    tag = system["tag"]
                    monomer_counts = system["monomer_counts"]

                    if name not in monomer_counts:
                        raise InputSchemaError(
                            f"Monomer {name!r} missing in system '{tag}'."
                        )

                    value = monomer_counts[name]
                    if not isinstance(value, int) or value < 0:
                        raise NumericFieldError(
                            f"Invalid count for monomer {name!r} in system '{tag}'."
                        )

                    count_map[tag] = value

                count = count_map
                ratio = None

            else:
                # ratio defined per system (must be consistent)
                first_system = systems[0]
                ratios = first_system["monomer_ratios"]

                if name not in ratios:
                    raise InputSchemaError(
                        f"Monomer {name!r} missing in ratio definition."
                    )

                ratio_value = ratios[name]

                if not isinstance(ratio_value, (int, float)) or ratio_value < 0:
                    raise NumericFieldError(
                        f"Invalid ratio for monomer {name!r}."
                    )

                ratio = float(ratio_value)
                count = None

            validated_monomers.append(
                MonomerEntry(
                    id=monomer_id,
                    data_id=f"data_{monomer_id}",
                    name=name,
                    smiles=smiles,
                    count=count,
                    ratio=ratio,
                    rdkit_mol=mol,
                )
            )

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
    
    def _validate_replicas(self, replicas_dict: dict, method: CompositionMethodType) -> dict:
        """
        Validates the 'replicas' section of the input.

        Ensures structure, numeric fields, and system-level composition fields
        are correct based on the selected composition method.

        Args:
            replicas_dict: The 'replicas' dictionary from user input.
            method: Composition method ("counts" or "stoichiometric_ratio").

        Returns:
            Normalized replicas dictionary.
        """

        temperatures = self._validate_temperature(
            replicas_dict.get("temperatures")
        )

        density = self._validate_density(
            replicas_dict.get("density")
        )

        systems = replicas_dict.get("systems")

        if not isinstance(systems, list) or not systems:
            raise InputSchemaError(
                "'replicas.systems' must be a non-empty list."
            )

        seen_tags: set[str] = set()

        # used for ratio consistency check
        reference_ratios: dict | None = None

        for system in systems:

            if not isinstance(system, dict):
                raise InputSchemaError(
                    "Each system must be a dictionary."
                )

            tag = system.get("tag")

            if not isinstance(tag, str) or not tag.strip():
                raise InputSchemaError(
                    "Each system must include a non-empty 'tag'."
                )

            if tag in seen_tags:
                raise InputSchemaError(
                    f"Duplicate system tag detected: {tag!r}"
                )

            seen_tags.add(tag)

            # -------- ratio mode --------
            if method == "stoichiometric_ratio":

                total_atoms = system.get("total_atoms")

                if not isinstance(total_atoms, int) or total_atoms <= 0:
                    raise NumericFieldError(
                        "'total_atoms' must be a positive integer in ratio mode."
                    )

                ratios = system.get("monomer_ratios")

                if not isinstance(ratios, dict) or not ratios:
                    raise InputSchemaError(
                        "'monomer_ratios' must be provided in ratio mode."
                    )

                # validate ratio values
                for monomer, value in ratios.items():
                    if not isinstance(value, (int, float)) or value < 0:
                        raise NumericFieldError(
                            f"Invalid ratio value for monomer {monomer!r}: {value!r}"
                        )

                # enforce ratio consistency across systems
                if reference_ratios is None:
                    reference_ratios = ratios
                else:
                    if ratios != reference_ratios:
                        raise InputSchemaError(
                            "All systems must use identical 'monomer_ratios'."
                        )

            # -------- counts mode --------
            if method == "counts":

                counts = system.get("monomer_counts")

                if not isinstance(counts, dict) or not counts:
                    raise InputSchemaError(
                        "'monomer_counts' must be provided in counts mode."
                    )

                for monomer, value in counts.items():

                    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                        raise NumericFieldError(
                            f"Invalid count for monomer {monomer!r}: {value!r}"
                        )

                if "total_atoms" in system:
                    raise InputSchemaError(
                        "'total_atoms' should not appear in counts mode."
                    )

        return {
            "method": method,
            "temperatures": temperatures,
            "density": density,
            "systems": systems
        }

if __name__ == "__main__":
    # Sample input data for quick manual verification of the parser.
        
    # Example 1: Counts Mode
    inputs = {
        "simulation_name": "Example_Count_Mode",

        "replicas": {
            "temperatures": [300, 400, 500],
            "density": 0.8,
            "method": "counts",
            "systems": [
                {
                    "tag": "10k",
                    "monomer_counts": {
                        "tmc": 220,
                        "mpd": 220,
                        "ethanol": 110
                    }
                },
                {
                    "tag": "100k",
                    "monomer_counts": {
                        "tmc": 2200,
                        "mpd": 2200,
                        "ethanol": 1100
                    }
                }
            ]
        },

        "monomers": [
            {
                "name": "tmc",
                "smiles": "ClC(=O)c1cc(cc(c1)C(Cl)=O)C(Cl)=O"
            },
            {
                "name": "mpd",
                "smiles": "C1=CC(=CC(=C1)N)N"
            },
            {
                "name": "ethanol",
                "smiles": "CCO"
            }
        ]
    }
        
    # Example 2: Stoichiometric Mode
    inputs_stoichiometric = {

        "simulation_name": "Example_Ratio_Mode",

        "replicas": {
            "temperatures": [300, 400],
            "density": [0.8],
            "method": "ratio",

            "systems": [

                {
                    "tag": "10k_base",
                    "total_atoms": 10000,

                    "monomer_ratios": {
                        "tmc": 1.0,
                        "mpd": 1.0,
                        "ethanol": 0.5
                    }
                },

                {
                    "tag": "100k_base",
                    "total_atoms": 100000,

                    "monomer_ratios": {
                        "tmc": 1.0,
                        "mpd": 1.0,
                        "ethanol": 0.5
                    }
                }

            ]
        },

        "monomers": [

            {
                "name": "tmc",
                "smiles": "ClC(=O)c1cc(cc(c1)C(Cl)=O)C(Cl)=O"
            },

            {
                "name": "mpd",
                "smiles": "C1=CC(=CC(=C1)N)N"
            },

            {
                "name": "ethanol",
                "smiles": "CCO"
            }

        ]
    }
    input_ff = {
        "simulation_name": "Example_Count_Mode",

        "replicas": {
            "temperatures": [300, 400, 500],
            "density": 0.8,
            "method": "counts",

            "systems": [
            {
                "tag": "10k",
                "monomer_counts": {
                "tmc": 220,
                "mpd": 220,
                "ethanol": 110
                }
            },
            {
                "tag": "100k",
                "monomer_counts": {
                "tmc": 2200,
                "mpd": 2200,
                "ethanol": 1100
                }
            }
            ]
        },

        "force_field": "PCFF",

        "monomers": [
            {
            "name": "tmc",
            "smiles": "ClC(=O)c1cc(cc(c1)C(Cl)=O)C(Cl)=O"
            },
            {
            "name": "mpd",
            "smiles": "C1=CC(=CC(=C1)N)N"
            },
            {
            "name": "ethanol",
            "smiles": "CCO"
            }
        ]
    }
    parser = InputParser()
    print("Validating Counts Mode Input:")
    print(parser.validate_inputs(inputs))
    
    print("\nValidating Ratio Mode Input:")
    print(parser.validate_inputs(inputs_stoichiometric))

    print("\nValidating Force Field Input:")
    print(parser.validate_inputs(input_ff))


