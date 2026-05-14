from __future__ import annotations
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, Optional

from PIL.Image import Image
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw

# Module-level logger for future diagnostics.
logger = logging.getLogger(__name__)


class InputError(Exception):
    """Base class for all input-related errors in the parsing pipeline."""


class InputSchemaError(InputError):
    """Raised when there are missing keys or wrong top-level structure in the input."""


class InputConflictError(InputError):
    """Raised when there are mutually exclusive or invalid combinations of input fields."""


class NumericFieldError(InputError):
    """Raised when numeric values such as density, temperature, or counts are invalid."""


class SmilesValidationError(InputError):
    """Raised when an invalid SMILES string is provided or RDKit parsing fails."""


class DuplicateMonomerError(InputError):
    """Raised when duplicate monomer definitions are detected."""


# Type aliases for clarity and validation.
CompositionMethodType = Literal["counts", "ratio"]
ForceFieldType = Literal[
    "PCFF-IFF",
    "PCFF",
    "Compass",
    "CVFF-IFF",
    "CVFF",
    "Clay-FF",
    "DRIEDING",
    "OPLS-AA",
]


@dataclass(slots=True)
class MonomerEntry:
    """
    Canonical internal representation of a monomer.

    This dataclass stores validated and processed information about a single
    monomer species, including chemical identity, RDKit-derived properties,
    and its quantity in each simulation replica.

    Attributes:
        id: Unique integer identifier for the monomer, e.g., 1, 2, 3.
        data_id: Original internal identifier string, e.g., "data_1".
        name: Human-readable name used for LAMMPS labeling.
        smiles: RDKit-canonical SMILES string representing the molecular structure.
        count: Dictionary mapping replica tags to integer counts in counts mode.
        ratio: Molar ratio of the monomer in ratio mode.
        rdkit_mol: RDKit Mol object corresponding to the monomer SMILES.
        molecule_3Dmol_path: Optional file path to the 3D .mol representation.
        num_atoms: Number of atoms in the monomer with hydrogens included.
        molecular_weight: Molecular weight of the monomer from RDKit.
        status: Boolean indicating whether the monomer should be included.
    """

    id: int
    data_id: str
    name: str | None
    smiles: str
    count: dict | None  # None only if ratio mode.
    ratio: float | None  # None only if counts mode.
    rdkit_mol: Chem.Mol | None = None
    molecule_3Dmol_path: Optional[Path] = None
    num_atoms: int | None = None
    molecular_weight: float | None = None
    status: bool = True


@dataclass(slots=True)
class Replica:
    """
    Structured representation of one simulation case.

    This class was originally introduced when the input schema used the
    top-level JSON key "replicas". The public input schema now uses
    "simulations", but this internal class is still kept as Replica because it
    represents one validated simulation condition/case used throughout the
    workflow.

    Attributes:
        tag: Unique identifier for the simulation case.
        temperature: Target temperature in Kelvin.
        density: Target density in g/cm^3.
        monomer_counts: Mapping of monomer names to counts in counts mode.
        monomer_ratios: Mapping of monomer names to ratios in ratio mode.
        total_atoms: Total target atom count, required for ratio mode.
        initial_box_volume: Estimated initial box volume in cubic Angstroms.
        initial_box_length: Estimated cubic box length in Angstroms.
    """

    tag: str
    temperature: float
    density: float
    monomer_counts: dict | None = None
    monomer_ratios: dict | None = None
    total_atoms: int | None = None
    initial_box_volume: float | None = None
    initial_box_length: float | None = None


@dataclass(slots=True)
class SimulationSetup:
    """
    Container that stores normalized and validated simulation inputs.

    This class acts as the primary data transfer object passed into later
    AutoREACTER workflow stages. Most of these data are consumed during LAMMPS
    input generation, but keeping them structured here makes downstream stages
    cleaner and easier to extend.

    Attributes:
        simulation_name: Name of the simulation, used for output organization.
        temperature: List of temperature values in Kelvin.
        density: List of density values in g/cm^3.
        force_field: Force field name. Defaults to "PCFF" if not provided.
        monomers: List of MonomerEntry objects representing the system composition.
        replicas: List of validated Replica objects.
        composition_method: Composition method, either "counts" or "ratio".
        composition: Normalized composition dictionary from the input simulations.
        ratio: Optional mapping of monomer IDs to ratios.
        number_of_total_atoms: Optional list of total atom targets.
        box_estimates: Optional estimated box size placeholder.
    """

    simulation_name: str
    temperature: list[float]
    density: list[float]
    force_field: str | None
    monomers: list[MonomerEntry]
    replicas: list[Replica] | None = None
    composition_method: CompositionMethodType | None = None
    composition: dict[str, Any] | None = None
    ratio: dict[int, float] | None = None
    number_of_total_atoms: list[int] | None = None
    box_estimates: float | None = None


class InputParser:
    """
    Validates, normalizes, and prepares raw user inputs for AutoREACTER.

    This class enforces required keys, checks value contracts, canonicalizes
    SMILES strings through RDKit, and transforms a raw dictionary into a
    structured SimulationSetup object.
    """

    def validate_inputs(self, inputs: dict) -> SimulationSetup:
        """
        Main entry point for input validation.

        Args:
            inputs: Raw dictionary input containing simulation parameters.

        Returns:
            A validated SimulationSetup object ready for downstream workflow steps.
        """
        self.validate_basic_format(inputs)

        simulation_name = inputs["simulation_name"]
        simulations_list = inputs["simulations"]
        composition_method = self._get_inputs_mode(simulations_list)

        validated_replicas = self._validate_replicas(
            simulations_list,
            composition_method,
        )

        self._validate_system_monomer_keys(
            inputs,
            validated_replicas["systems"],
            composition_method,
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
            replicas=validated_replicas["replicas"],
            composition_method=composition_method,
            composition=validated_replicas,
            force_field=force_field,
        )

    def molecule_representation_of_initial_molecules(
        self,
        simulation_setup: SimulationSetup,
    ) -> tuple[list[Chem.Mol], list[str]]:
        """
        Extracts RDKit Mol objects and names for all monomers.

        Args:
            simulation_setup: Validated SimulationSetup object.

        Returns:
            A tuple containing RDKit Mol objects and molecule legends.
        """
        initial_molecules = [
            monomer.rdkit_mol
            for monomer in simulation_setup.monomers
            if monomer.rdkit_mol is not None
        ]
        initial_molecules_legends = [
            monomer.name or monomer.data_id
            for monomer in simulation_setup.monomers
        ]
        return initial_molecules, initial_molecules_legends

    def initial_molecules_image_grid(self, simulation_setup: SimulationSetup) -> Image:
        """
        Creates a grid image of the initial molecules for visualization.

        Args:
            simulation_setup: Validated SimulationSetup object.

        Returns:
            PIL Image object containing the molecule grid.
        """
        initial_molecules, initial_molecules_legends = (
            self.molecule_representation_of_initial_molecules(simulation_setup)
        )

        return Draw.MolsToGridImage(
            initial_molecules,
            molsPerRow=3,
            subImgSize=(400, 400),
            legends=initial_molecules_legends,
        )

    def validate_basic_format(self, inputs: dict) -> None:
        """
        Validates the top-level structure of the input dictionary.

        Args:
            inputs: Raw input dictionary.

        Raises:
            InputSchemaError: If input is not a dict or is missing required keys.
        """
        if not isinstance(inputs, dict):
            raise InputSchemaError(
                f"Expected input to be a dictionary. Got {type(inputs).__name__} instead."
            )

        required_keys = ["simulation_name", "simulations", "monomers"]
        for key in required_keys:
            if key not in inputs:
                raise InputSchemaError(
                    f"Missing required key: {key!r} in inputs dictionary."
                )

    def _get_inputs_mode(self, simulations_list: list) -> CompositionMethodType:
        """
        Determines the composition method from the simulations list.

        Args:
            simulations_list: The 'simulations' list from the input dictionary.

        Returns:
            Composition method, either "counts" or "ratio".

        Raises:
            InputSchemaError: If the list is invalid or mode cannot be inferred.
            InputConflictError: If a replica contains both counts and ratios.
        """
        if not isinstance(simulations_list, list) or len(simulations_list) == 0:
            raise InputSchemaError(
                f"'simulations' must be a non-empty list. Got {type(simulations_list).__name__} instead."
            )

        detected_modes: set[CompositionMethodType] = set()

        for replica in simulations_list:
            if not isinstance(replica, dict):
                raise InputSchemaError(
                    f"Each replica must be a dictionary. Got {type(replica).__name__} instead."
                )

            has_counts = "monomer_counts" in replica
            has_ratios = "monomer_ratios" in replica

            if has_counts and has_ratios:
                raise InputConflictError(
                    "Each replica must contain either 'monomer_counts' or 'monomer_ratios', not both."
                )

            if has_counts:
                detected_modes.add("counts")
            elif has_ratios:
                detected_modes.add("ratio")
            else:
                raise InputSchemaError(
                    "Could not infer composition method. Replicas must contain either "
                    "'monomer_counts' or 'monomer_ratios'."
                )

        if len(detected_modes) != 1:
            raise InputConflictError(
                "All simulations must use the same composition method. Do not mix counts and ratio modes."
            )

        return detected_modes.pop()

    def _validate_temperature(self, temp: Any) -> float:
        """
        Validates and normalizes temperature input.

        Args:
            temp: Temperature input.

        Returns:
            Validated float temperature value.

        Raises:
            NumericFieldError: If the value is not a positive number.
        """
        if isinstance(temp, bool) or not isinstance(temp, (int, float)):
            raise NumericFieldError(f"'temperature' must be a number. Got: {temp!r}")

        if temp <= 0:
            raise NumericFieldError(
                f"Temperature values must be positive. Got: {temp!r}"
            )

        return float(temp)

    def _validate_density(self, density: Any) -> float:
        """
        Validates the density input.

        Args:
            density: Density value in g/cm^3.

        Returns:
            Validated float density value.

        Raises:
            NumericFieldError: If density is not a positive number.
        """
        if isinstance(density, bool) or not isinstance(density, (int, float)):
            raise NumericFieldError(
                f"'density' must be a number. Got: {density!r}"
            )

        density_value = float(density)

        if density_value <= 0:
            raise NumericFieldError(
                f"Density values must be positive. Got: {density!r}"
            )

        return density_value

    def _validate_force_field(self, force_field: Any) -> str:
        """
        Validates the force field input.

        Args:
            force_field: Force field name as a string.

        Returns:
            Validated force field name.

        Raises:
            InputSchemaError: If force field is unsupported.
        """
        allowed: set[ForceFieldType] = {
            "PCFF-IFF",
            "PCFF",
            "Compass",
            "CVFF-IFF",
            "CVFF",
            "Clay-FF",
            "DRIEDING",
            "OPLS-AA",
        }

        if force_field is None:
            return "PCFF"

        if not isinstance(force_field, str) or not force_field.strip():
            raise InputSchemaError(
                f"'force_field' must be a non-empty string. Got: {force_field!r}"
            )

        if force_field not in allowed:
            raise InputSchemaError(f"Unsupported force field: {force_field!r}")

        return force_field

    def _validate_composition(
        self,
        composition_dict: dict,
        method: CompositionMethodType,
    ) -> dict:
        """
        Validates a legacy composition dictionary structure.

        Note:
            This method is kept for compatibility with older input schemas. The
            current parser path validates composition through the 'simulations' section.

        Args:
            composition_dict: Composition dictionary to validate.
            method: Composition method, either "counts" or "ratio".

        Returns:
            The validated composition dictionary.
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

            if method == "ratio":
                total_atoms = target.get("total_atoms")
                if isinstance(total_atoms, bool) or not isinstance(total_atoms, int) or total_atoms <= 0:
                    raise NumericFieldError(
                        f"'total_atoms' must be a positive integer for ratio mode. Got: {total_atoms!r}"
                    )

            if method == "counts" and "total_atoms" in target:
                raise InputSchemaError(
                    "'total_atoms' should not be provided in 'counts' mode."
                )

        return composition_dict

    def _validate_system_monomer_keys(
        self,
        inputs: dict,
        systems: list[dict],
        method: CompositionMethodType,
    ) -> None:
        """
        Ensures that monomer keys in each system match the defined monomers.

        Args:
            inputs: Raw input dictionary containing the top-level 'monomers' list.
            systems: List of system dictionaries from the simulations section.
            method: Composition method, either "counts" or "ratio".

        Raises:
            InputSchemaError: If systems contain missing or unknown monomer keys.
        """
        monomers = inputs.get("monomers", [])
        allowed_names: set[str] = set()

        for monomer_id, monomer in enumerate(monomers, start=1):
            name = monomer.get("name")

            if not isinstance(name, str) or not name.strip():
                name = f"data_{monomer_id}"
                monomer["name"] = name

            allowed_names.add(name)

        field_name = "monomer_counts" if method == "counts" else "monomer_ratios"

        for system in systems:
            tag = system["tag"]
            provided_names = set(system.get(field_name, {}).keys())

            extra = provided_names - allowed_names
            missing = allowed_names - provided_names

            if extra:
                raise InputSchemaError(
                    f"System '{tag}' has unknown monomer(s) in '{field_name}': {sorted(extra)}"
                )

            if missing:
                raise InputSchemaError(
                    f"System '{tag}' is missing monomer(s) in '{field_name}': {sorted(missing)}"
                )

    def _validate_monomer_entry(
        self,
        inputs: dict,
        method: CompositionMethodType,
        systems: list[dict],
    ) -> list[MonomerEntry]:
        """
        Validates monomer entries and constructs MonomerEntry objects.

        Args:
            inputs: Raw input dictionary.
            method: Composition method, either "counts" or "ratio".
            systems: Validated system dictionaries from the simulations section.

        Returns:
            A list of validated MonomerEntry objects.
        """
        validated_monomers: list[MonomerEntry] = []
        seen_smiles: list[str] = []

        monomers = inputs.get("monomers")
        if not isinstance(monomers, list):
            raise InputSchemaError("'monomers' must be a list.")

        for monomer_id, monomer_dict in enumerate(monomers, start=1):
            if not isinstance(monomer_dict, dict):
                raise InputSchemaError(
                    f"Each monomer entry must be a dictionary. Got: {monomer_dict!r}"
                )

            name = monomer_dict.get("name")
            if not isinstance(name, str) or not name.strip():
                name = f"data_{monomer_id}"

            smiles_raw = monomer_dict.get("smiles")
            smiles, mol = self._validate_smiles(smiles_raw)
            seen_smiles = self.validate_no_duplicate_smiles(smiles, seen_smiles)

            num_atoms, molecular_weight = self._derive_molecule_properties(mol)

            if method == "counts":
                count_map: dict[str, int] = {}

                for system in systems:
                    tag = system["tag"]
                    monomer_counts = system["monomer_counts"]

                    if name not in monomer_counts:
                        raise InputSchemaError(
                            f"Monomer {name!r} missing in system '{tag}'."
                        )

                    value = monomer_counts[name]
                    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                        raise NumericFieldError(
                            f"Invalid count for monomer {name!r} in system '{tag}'."
                        )

                    count_map[tag] = value

                count = count_map
                ratio = None

            else:
                first_system = systems[0]
                ratios = first_system["monomer_ratios"]

                if name not in ratios:
                    raise InputSchemaError(
                        f"Monomer {name!r} missing in ratio definition."
                    )

                ratio_value = ratios[name]

                if isinstance(ratio_value, bool) or not isinstance(ratio_value, (int, float)) or ratio_value < 0:
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
                    num_atoms=num_atoms,
                    molecular_weight=molecular_weight,
                )
            )

        return validated_monomers

    def _derive_molecule_properties(self, mol: Chem.Mol) -> tuple[int, float]:
        """
        Derives simple molecule properties from an RDKit Mol object.

        Args:
            mol: RDKit molecule.

        Returns:
            Tuple of number of atoms with hydrogens included and molecular weight.
        """
        mol_with_h = Chem.AddHs(mol)
        num_atoms = int(mol_with_h.GetNumAtoms())
        molecular_weight = float(Descriptors.MolWt(mol))
        return num_atoms, molecular_weight

    def _int_to_dict(self, integer: int) -> dict:
        """
        Converts an integer to a dictionary format for counts mode.

        Note:
            This is kept as a compatibility placeholder for future ratio-mode
            conversion where counts may be derived from ratios and total atoms.
        """
        return {"_": integer}

    def _validate_smiles(self, smiles: Any) -> tuple[str, Chem.Mol]:
        """
        Validates a SMILES string through RDKit and canonicalizes it.

        Args:
            smiles: Raw SMILES string input.

        Returns:
            Canonical SMILES string and RDKit Mol object.

        Raises:
            SmilesValidationError: If RDKit cannot parse the string.
        """
        if not isinstance(smiles, str) or not smiles.strip():
            raise SmilesValidationError(
                f"SMILES must be a non-empty string. Got: {smiles!r}"
            )

        smiles_clean = smiles.strip()
        mol = Chem.MolFromSmiles(smiles_clean)

        if mol is None:
            raise SmilesValidationError(
                f"Invalid SMILES string: {smiles!r}. RDKit failed to parse it."
            )

        canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        return canonical_smiles, mol

    def _validate_numeric_fields(self, inputs: dict) -> None:
        """
        Legacy numeric validation helper.

        Note:
            This method is kept for compatibility with older input schemas and is
            not currently used in the main validate_inputs flow.
        """
        density = inputs.get("density", None)
        if isinstance(density, bool) or not isinstance(density, (int, float)) or density <= 0:
            raise NumericFieldError(
                f"'density' must be a positive number. Got: {density!r}"
            )

        temps = inputs.get("temperature", None)
        if temps is None:
            raise NumericFieldError(
                "'temperature' is required as a number or list of numbers."
            )

        temps_list = temps if isinstance(temps, list) else [temps]
        for temp in temps_list:
            if isinstance(temp, bool) or not isinstance(temp, (int, float)) or temp <= 0:
                raise NumericFieldError(
                    f"Temperature values must be positive numbers. Got: {temp!r}"
                )

        num_monomers = inputs.get("number_of_monomers", None)
        if not isinstance(num_monomers, dict) or not num_monomers:
            raise NumericFieldError(
                "'number_of_monomers' must be a non-empty dict of monomer_id -> positive int."
            )

        for monomer_id, count in num_monomers.items():
            if isinstance(count, bool) or not isinstance(count, int) or count <= 0:
                raise NumericFieldError(
                    f"Monomer count for {monomer_id!r} must be a positive integer. Got: {count!r}"
                )

    def validate_no_duplicate_smiles(
        self,
        current_monomer: str,
        seen_monomer_list: list[str],
    ) -> list[str]:
        """
        Confirms that there are no duplicate canonical SMILES strings.

        Args:
            current_monomer: Canonical SMILES string for the current monomer.
            seen_monomer_list: List of previously seen canonical SMILES strings.

        Returns:
            Updated list of seen SMILES strings.
        """
        if current_monomer in seen_monomer_list:
            raise DuplicateMonomerError(
                f"Duplicate monomer detected with SMILES: {current_monomer!r}. "
                "Each monomer must have a unique SMILES string."
            )

        seen_monomer_list.append(current_monomer)
        return seen_monomer_list

    def _validate_single_replica(
        self,
        replica: Replica,
        method: CompositionMethodType,
    ) -> None:
        """
        Validates one Replica object after initial construction.
        """
        if not isinstance(replica.tag, str) or not replica.tag.strip():
            raise InputSchemaError(
                "Each system must include a non-empty 'tag'."
            )

        if isinstance(replica.temperature, bool) or not isinstance(replica.temperature, (int, float)) or replica.temperature <= 0:
            raise NumericFieldError(
                f"Replica '{replica.tag}' has invalid temperature."
            )

        if isinstance(replica.density, bool) or not isinstance(replica.density, (int, float)) or replica.density <= 0:
            raise NumericFieldError(
                f"Replica '{replica.tag}' has invalid density."
            )

        if method == "counts":
            if not isinstance(replica.monomer_counts, dict) or not replica.monomer_counts:
                raise InputSchemaError(
                    "'monomer_counts' must be provided in counts mode."
                )

            for monomer, value in replica.monomer_counts.items():
                if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                    raise NumericFieldError(
                        f"Invalid count for monomer {monomer!r}: {value!r}"
                    )

        elif method == "ratio":
            if isinstance(replica.total_atoms, bool) or not isinstance(replica.total_atoms, int) or replica.total_atoms <= 0:
                raise NumericFieldError(
                    "'total_atoms' must be a positive integer in ratio mode."
                )

            if not isinstance(replica.monomer_ratios, dict) or not replica.monomer_ratios:
                raise InputSchemaError(
                    "'monomer_ratios' must be provided in ratio mode."
                )

            for monomer, value in replica.monomer_ratios.items():
                if isinstance(value, bool) or not isinstance(value, (int, float)) or value < 0:
                    raise NumericFieldError(
                        f"Invalid ratio value for monomer {monomer!r}: {value!r}"
                    )

    def _validate_replicas(
        self,
        systems: list[dict],
        method: CompositionMethodType,
    ) -> dict:
        """
        Validates the 'simulations' section of the input.

        Args:
            systems: List of replica dictionaries.
            method: Composition method, either "counts" or "ratio".

        Returns:
            Normalized replicas dictionary.
        """
        if not isinstance(systems, list) or not systems:
            raise InputSchemaError(
                "'simulations' must be a non-empty list."
            )

        seen_tags: set[str] = set()
        temperatures: list[float] = []
        density: list[float] = []
        replicas: list[Replica] = []
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

            system["temperature"] = self._validate_temperature(
                system.get("temperature")
            )
            system["density"] = self._validate_density(
                system.get("density")
            )

            temperatures.append(system["temperature"])
            density.append(system["density"])

            if method == "ratio":
                total_atoms = system.get("total_atoms")

                if isinstance(total_atoms, bool) or not isinstance(total_atoms, int) or total_atoms <= 0:
                    raise NumericFieldError(
                        "'total_atoms' must be a positive integer in ratio mode."
                    )

                ratios = system.get("monomer_ratios")

                if not isinstance(ratios, dict) or not ratios:
                    raise InputSchemaError(
                        "'monomer_ratios' must be provided in ratio mode."
                    )

                for monomer, value in ratios.items():
                    if isinstance(value, bool) or not isinstance(value, (int, float)) or value < 0:
                        raise NumericFieldError(
                            f"Invalid ratio value for monomer {monomer!r}: {value!r}"
                        )

                if reference_ratios is None:
                    reference_ratios = ratios
                elif ratios != reference_ratios:
                    raise InputSchemaError(
                        "All systems must use identical 'monomer_ratios'."
                    )

                replica = Replica(
                    tag=system["tag"],
                    temperature=system["temperature"],
                    density=system["density"],
                    monomer_counts=None,
                    monomer_ratios=ratios,
                    total_atoms=total_atoms,
                )

                self._validate_single_replica(replica, method)
                replicas.append(replica)

            elif method == "counts":
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

                replica = Replica(
                    tag=system["tag"],
                    temperature=system["temperature"],
                    density=system["density"],
                    monomer_counts=counts,
                    monomer_ratios=None,
                    total_atoms=None,
                )

                self._validate_single_replica(replica, method)
                replicas.append(replica)

        return {
            "method": method,
            "temperatures": temperatures,
            "density": density,
            "systems": systems,
            "replicas": replicas,
        }


if __name__ == "__main__":
    # Sample input data for quick manual verification of the parser.

    # Example 1: Counts Mode
    inputs = {
        "simulation_name": "Example_Count_Mode",
        "simulations": [
            {
                "tag": "10k",
                "temperature": 300,
                "density": 0.8,
                "monomer_counts": {
                    "tmc": 220,
                    "mpd": 220,
                    "ethanol": 110,
                },
            },
            {
                "tag": "100k",
                "temperature": 400,
                "density": 0.8,
                "monomer_counts": {
                    "tmc": 2200,
                    "mpd": 2200,
                    "ethanol": 1100,
                },
            },
        ],
        "monomers": [
            {
                "name": "tmc",
                "smiles": "ClC(=O)c1cc(cc(c1)C(Cl)=O)C(Cl)=O",
            },
            {
                "name": "mpd",
                "smiles": "C1=CC(=CC(=C1)N)N",
            },
            {
                "name": "ethanol",
                "smiles": "CCO",
            },
        ],
    }

    # Example 2: Ratio Mode
    inputs_ratio = {
        "simulation_name": "Example_Ratio_Mode",
        "simulations": [
            {
                "tag": "10k_base",
                "temperature": 300,
                "density": 0.8,
                "total_atoms": 10000,
                "monomer_ratios": {
                    "tmc": 1.0,
                    "mpd": 1.0,
                    "ethanol": 0.5,
                },
            },
            {
                "tag": "100k_base",
                "temperature": 400,
                "density": 0.8,
                "total_atoms": 100000,
                "monomer_ratios": {
                    "tmc": 1.0,
                    "mpd": 1.0,
                    "ethanol": 0.5,
                },
            },
        ],
        "monomers": [
            {
                "name": "tmc",
                "smiles": "ClC(=O)c1cc(cc(c1)C(Cl)=O)C(Cl)=O",
            },
            {
                "name": "mpd",
                "smiles": "C1=CC(=CC(=C1)N)N",
            },
            {
                "name": "ethanol",
                "smiles": "CCO",
            },
        ],
    }

    input_ff = {
        "simulation_name": "Example_Count_Mode_FF",
        "force_field": "PCFF",
        "simulations": [
            {
                "tag": "10k",
                "temperature": 300,
                "density": 0.8,
                "monomer_counts": {
                    "tmc": 220,
                    "mpd": 220,
                    "ethanol": 110,
                },
            },
            {
                "tag": "100k",
                "temperature": 400,
                "density": 0.8,
                "monomer_counts": {
                    "tmc": 2200,
                    "mpd": 2200,
                    "ethanol": 1100,
                },
            },
        ],
        "monomers": [
            {
                "name": "tmc",
                "smiles": "ClC(=O)c1cc(cc(c1)C(Cl)=O)C(Cl)=O",
            },
            {
                "name": "mpd",
                "smiles": "C1=CC(=CC(=C1)N)N",
            },
            {
                "name": "ethanol",
                "smiles": "CCO",
            },
        ],
    }

    parser = InputParser()

    print("Validating Counts Mode Input:")
    print(parser.validate_inputs(inputs))

    print("\nValidating Ratio Mode Input:")
    print(parser.validate_inputs(inputs_ratio))

    print("\nValidating Force Field Input:")
    print(parser.validate_inputs(input_ff))
