"""
Comprehensive unit tests for :mod:`AutoREACTER.input_parser`.

The suite exercises the current counts-mode and ratio-mode schemas, numeric
validation, force-field normalization, SMILES handling, monomer construction,
loop configuration, and the public ``validate_inputs`` workflow.

Run from the repository root with either:

    python -m unittest tests.test_input_parser -v

or:

    pytest -q tests/test_input_parser.py
"""

from __future__ import annotations

import copy
import unittest
from unittest.mock import patch

from rdkit import Chem

from AutoREACTER.input_parser import (
    CompatibilityError,
    DuplicateMonomerError,
    InputConflictError,
    InputParser,
    InputSchemaError,
    MonomerEntry,
    NumericFieldError,
    Simulation,
    SimulationSetup,
    SmilesValidationError,
)


class InputParserTestCase(unittest.TestCase):
    """Shared fixtures and helpers for input-parser tests."""

    def setUp(self) -> None:
        """Create a fresh parser before every test."""
        self.parser = InputParser()

    @staticmethod
    def counts_input() -> dict:
        """Return a fresh, valid counts-mode input dictionary."""
        return {
            "simulation_name": "counts_example",
            "force_field": "pcff",
            "loop": False,
            "simulations": [
                {
                    "tag": "small",
                    "temperature": 300,
                    "density": 0.80,
                    "monomer_counts": {
                        "ethanol": 2,
                        "water": 1,
                    },
                },
                {
                    "tag": "large",
                    "temperature": 350.0,
                    "density": 0.95,
                    "monomer_counts": {
                        "ethanol": 20,
                        "water": 10,
                    },
                },
            ],
            "monomers": [
                {
                    "name": "ethanol",
                    "smiles": "CCO",
                },
                {
                    "name": "water",
                    "smiles": "O",
                },
            ],
        }

    @staticmethod
    def ratio_input() -> dict:
        """Return a fresh, valid ratio-mode input dictionary."""
        return {
            "simulation_name": "ratio_example",
            "force_field": "PCFF-IFF",
            "simulations": [
                {
                    "tag": "small",
                    "temperature": 300,
                    "density": 0.80,
                    "total_atoms": 1000,
                    "monomer_ratios": {
                        "ethanol": 2.0,
                        "water": 1.0,
                    },
                },
                {
                    "tag": "large",
                    "temperature": 350,
                    "density": 0.95,
                    "total_atoms": 10000,
                    "monomer_ratios": {
                        "ethanol": 2.0,
                        "water": 1.0,
                    },
                },
            ],
            "monomers": [
                {
                    "name": "ethanol",
                    "smiles": "CCO",
                },
                {
                    "name": "water",
                    "smiles": "O",
                },
            ],
        }


class TestBasicFormat(InputParserTestCase):
    """Tests for top-level input structure validation."""

    def test_validate_basic_format_accepts_required_keys(self) -> None:
        inputs = {
            "simulation_name": "demo",
            "simulations": [],
            "monomers": [],
        }

        self.assertIsNone(self.parser.validate_basic_format(inputs))

    def test_validate_basic_format_rejects_non_dictionary(self) -> None:
        for value in (None, [], "input", 12):
            with self.subTest(value=value):
                with self.assertRaises(InputSchemaError):
                    self.parser.validate_basic_format(value)

    def test_validate_basic_format_rejects_each_missing_key(self) -> None:
        valid = {
            "simulation_name": "demo",
            "simulations": [],
            "monomers": [],
        }

        for missing_key in tuple(valid):
            with self.subTest(missing_key=missing_key):
                inputs = valid.copy()
                inputs.pop(missing_key)

                with self.assertRaises(InputSchemaError):
                    self.parser.validate_basic_format(inputs)


class TestCompositionModeDetection(InputParserTestCase):
    """Tests for counts/ratio mode inference."""

    def test_get_inputs_mode_detects_counts(self) -> None:
        simulations = [
            {
                "tag": "a",
                "monomer_counts": {"ethanol": 1},
            }
        ]

        self.assertEqual(self.parser._get_inputs_mode(simulations), "counts")

    def test_get_inputs_mode_detects_ratio(self) -> None:
        simulations = [
            {
                "tag": "a",
                "monomer_ratios": {"ethanol": 1.0},
            }
        ]

        self.assertEqual(self.parser._get_inputs_mode(simulations), "ratio")

    def test_get_inputs_mode_rejects_empty_or_non_list_value(self) -> None:
        for value in (None, {}, (), [], "simulations"):
            with self.subTest(value=value):
                with self.assertRaises(InputSchemaError):
                    self.parser._get_inputs_mode(value)

    def test_get_inputs_mode_rejects_non_dictionary_simulation(self) -> None:
        with self.assertRaises(InputSchemaError):
            self.parser._get_inputs_mode(["not-a-dictionary"])

    def test_get_inputs_mode_rejects_missing_composition_field(self) -> None:
        with self.assertRaises(InputSchemaError):
            self.parser._get_inputs_mode(
                [
                    {
                        "tag": "a",
                        "temperature": 300,
                        "density": 0.8,
                    }
                ]
            )

    def test_get_inputs_mode_rejects_both_fields_in_one_simulation(self) -> None:
        with self.assertRaises(InputConflictError):
            self.parser._get_inputs_mode(
                [
                    {
                        "tag": "a",
                        "monomer_counts": {"ethanol": 1},
                        "monomer_ratios": {"ethanol": 1.0},
                    }
                ]
            )

    def test_get_inputs_mode_rejects_mixed_modes(self) -> None:
        with self.assertRaises(InputConflictError):
            self.parser._get_inputs_mode(
                [
                    {
                        "tag": "counts",
                        "monomer_counts": {"ethanol": 1},
                    },
                    {
                        "tag": "ratio",
                        "monomer_ratios": {"ethanol": 1.0},
                    },
                ]
            )


class TestNumericValidation(InputParserTestCase):
    """Tests for temperature and density helpers."""

    def test_validate_temperature_returns_float(self) -> None:
        self.assertEqual(self.parser._validate_temperature(300), 300.0)
        self.assertEqual(self.parser._validate_temperature(300.5), 300.5)

    def test_validate_temperature_rejects_non_numeric_and_boolean(self) -> None:
        for value in (None, "300", True, False, [], {}):
            with self.subTest(value=value):
                with self.assertRaises(NumericFieldError):
                    self.parser._validate_temperature(value)

    def test_validate_temperature_rejects_zero_and_negative_values(self) -> None:
        for value in (0, -1, -273.15):
            with self.subTest(value=value):
                with self.assertRaises(NumericFieldError):
                    self.parser._validate_temperature(value)

    def test_validate_density_returns_float(self) -> None:
        self.assertEqual(self.parser._validate_density(1), 1.0)
        self.assertEqual(self.parser._validate_density(0.85), 0.85)

    def test_validate_density_rejects_non_numeric_and_boolean(self) -> None:
        for value in (None, "0.85", True, False, [], {}):
            with self.subTest(value=value):
                with self.assertRaises(NumericFieldError):
                    self.parser._validate_density(value)

    def test_validate_density_rejects_zero_and_negative_values(self) -> None:
        for value in (0, -0.1, -1):
            with self.subTest(value=value):
                with self.assertRaises(NumericFieldError):
                    self.parser._validate_density(value)


class TestForceFieldValidation(InputParserTestCase):
    """Tests for force-field defaults, aliases, and compatibility."""

    def test_validate_force_field_defaults_to_pcff(self) -> None:
        self.assertEqual(self.parser._validate_force_field(None), "PCFF")

    def test_validate_force_field_normalizes_supported_aliases(self) -> None:
        aliases = {
            "pcff": "PCFF",
            " PCFF-IFF ": "PCFF-IFF",
            "cvff": "CVFF",
            "cvff-iff": "CVFF-IFF",
            "clayff": "Clay-FF",
            "clay-ff": "Clay-FF",
            "dreiding": "DREIDING",
            "drieding": "DREIDING",
        }

        for raw_value, expected in aliases.items():
            with self.subTest(raw_value=raw_value):
                self.assertEqual(
                    self.parser._validate_force_field(raw_value),
                    expected,
                )

    def test_validate_force_field_rejects_empty_or_non_string_value(self) -> None:
        for value in ("", "   ", 12, True, [], {}):
            with self.subTest(value=value):
                with self.assertRaises(InputSchemaError):
                    self.parser._validate_force_field(value)

    def test_validate_force_field_rejects_unknown_name(self) -> None:
        with self.assertRaises(InputSchemaError):
            self.parser._validate_force_field("not-a-force-field")

    def test_validate_force_field_rejects_currently_incompatible_fields(self) -> None:
        for value in ("oplsaa", "opls", "opls-aa", "gaff"):
            with self.subTest(value=value):
                with self.assertRaises(CompatibilityError):
                    self.parser._validate_force_field(value)


class TestSmilesValidation(InputParserTestCase):
    """Tests for SMILES parsing, canonicalization, and duplicates."""

    def test_validate_smiles_returns_canonical_smiles_and_molecule(self) -> None:
        smiles, mol = self.parser._validate_smiles(" C(C)O ")

        self.assertEqual(smiles, "CCO")
        self.assertIsInstance(mol, Chem.Mol)

    def test_validate_smiles_rejects_empty_or_non_string_values(self) -> None:
        for value in (None, "", "   ", 10, True, [], {}):
            with self.subTest(value=value):
                with self.assertRaises(SmilesValidationError):
                    self.parser._validate_smiles(value)

    def test_validate_smiles_rejects_invalid_smiles(self) -> None:
        with self.assertRaises(SmilesValidationError):
            self.parser._validate_smiles("not_a_smiles")

    def test_validate_no_duplicate_smiles_appends_unique_value(self) -> None:
        seen = ["O"]

        result = self.parser.validate_no_duplicate_smiles("CCO", seen)

        self.assertIs(result, seen)
        self.assertEqual(result, ["O", "CCO"])

    def test_validate_no_duplicate_smiles_rejects_duplicate_value(self) -> None:
        with self.assertRaises(DuplicateMonomerError):
            self.parser.validate_no_duplicate_smiles("CCO", ["CCO"])


class TestMoleculeProperties(InputParserTestCase):
    """Tests for RDKit-derived monomer properties."""

    def test_derive_molecule_properties_includes_hydrogen_atoms(self) -> None:
        mol = Chem.MolFromSmiles("CCO")
        self.assertIsNotNone(mol)

        num_atoms, molecular_weight = self.parser._derive_molecule_properties(
            mol
        )

        self.assertEqual(num_atoms, 9)
        self.assertAlmostEqual(molecular_weight, 46.069, places=2)

    def test_int_to_dict_retains_legacy_shape(self) -> None:
        self.assertEqual(self.parser._int_to_dict(7), {"_": 7})


class TestLegacyCompositionValidation(InputParserTestCase):
    """Tests for the legacy composition validation helper."""

    def test_validate_composition_accepts_counts_targets(self) -> None:
        composition = {
            "targets": [
                {"tag": "a"},
                {"tag": "b"},
            ]
        }

        self.assertIs(
            self.parser._validate_composition(composition, "counts"),
            composition,
        )

    def test_validate_composition_accepts_ratio_total_atoms(self) -> None:
        composition = {
            "targets": [
                {"tag": "a", "total_atoms": 1000},
                {"tag": "b", "total_atoms": 2000},
            ]
        }

        self.assertIs(
            self.parser._validate_composition(composition, "ratio"),
            composition,
        )

    def test_validate_composition_rejects_missing_or_empty_targets(self) -> None:
        for composition in ({}, {"targets": None}, {"targets": []}):
            with self.subTest(composition=composition):
                with self.assertRaises(InputSchemaError):
                    self.parser._validate_composition(
                        composition,
                        "counts",
                    )

    def test_validate_composition_rejects_duplicate_tags(self) -> None:
        with self.assertRaises(InputSchemaError):
            self.parser._validate_composition(
                {
                    "targets": [
                        {"tag": "same"},
                        {"tag": "same"},
                    ]
                },
                "counts",
            )

    def test_validate_composition_rejects_invalid_ratio_total_atoms(self) -> None:
        for value in (None, 0, -1, 1.5, True):
            with self.subTest(value=value):
                with self.assertRaises(NumericFieldError):
                    self.parser._validate_composition(
                        {
                            "targets": [
                                {
                                    "tag": "a",
                                    "total_atoms": value,
                                }
                            ]
                        },
                        "ratio",
                    )

    def test_validate_composition_rejects_total_atoms_in_counts_mode(self) -> None:
        with self.assertRaises(InputSchemaError):
            self.parser._validate_composition(
                {
                    "targets": [
                        {
                            "tag": "a",
                            "total_atoms": 1000,
                        }
                    ]
                },
                "counts",
            )


class TestSimulationValidation(InputParserTestCase):
    """Tests for individual and grouped simulation validation."""

    def test_validate_single_simulation_accepts_counts_mode(self) -> None:
        simulation = Simulation(
            tag="a",
            temperature=300.0,
            density=0.8,
            monomer_counts={"ethanol": 0},
        )

        self.assertIsNone(
            self.parser._validate_single_simulation(simulation, "counts")
        )

    def test_validate_single_simulation_accepts_ratio_mode(self) -> None:
        simulation = Simulation(
            tag="a",
            temperature=300.0,
            density=0.8,
            monomer_ratios={"ethanol": 0.0},
            total_atoms=1000,
        )

        self.assertIsNone(
            self.parser._validate_single_simulation(simulation, "ratio")
        )

    def test_validate_single_simulation_rejects_invalid_common_fields(self) -> None:
        cases = (
            Simulation(
                tag="",
                temperature=300.0,
                density=0.8,
                monomer_counts={"ethanol": 1},
            ),
            Simulation(
                tag="a",
                temperature=0,
                density=0.8,
                monomer_counts={"ethanol": 1},
            ),
            Simulation(
                tag="a",
                temperature=300.0,
                density=0,
                monomer_counts={"ethanol": 1},
            ),
        )

        expected_errors = (
            InputSchemaError,
            NumericFieldError,
            NumericFieldError,
        )

        for simulation, expected_error in zip(cases, expected_errors):
            with self.subTest(simulation=simulation):
                with self.assertRaises(expected_error):
                    self.parser._validate_single_simulation(
                        simulation,
                        "counts",
                    )

    def test_validate_single_simulation_rejects_invalid_count(self) -> None:
        simulation = Simulation(
            tag="a",
            temperature=300.0,
            density=0.8,
            monomer_counts={"ethanol": -1},
        )

        with self.assertRaises(NumericFieldError):
            self.parser._validate_single_simulation(simulation, "counts")

    def test_validate_single_simulation_rejects_invalid_ratio(self) -> None:
        simulation = Simulation(
            tag="a",
            temperature=300.0,
            density=0.8,
            monomer_ratios={"ethanol": -0.1},
            total_atoms=1000,
        )

        with self.assertRaises(NumericFieldError):
            self.parser._validate_single_simulation(simulation, "ratio")

    def test_validate_simulations_normalizes_counts_mode(self) -> None:
        systems = copy.deepcopy(self.counts_input()["simulations"])

        result = self.parser._validate_simulations(systems, "counts")

        self.assertEqual(result["method"], "counts")
        self.assertEqual(result["temperatures"], [300.0, 350.0])
        self.assertEqual(result["density"], [0.8, 0.95])
        self.assertEqual(len(result["simulations"]), 2)
        self.assertTrue(
            all(
                isinstance(simulation, Simulation)
                for simulation in result["simulations"]
            )
        )

    def test_validate_simulations_normalizes_ratio_mode(self) -> None:
        systems = copy.deepcopy(self.ratio_input()["simulations"])

        result = self.parser._validate_simulations(systems, "ratio")

        self.assertEqual(result["method"], "ratio")
        self.assertEqual(result["temperatures"], [300.0, 350.0])
        self.assertEqual(result["density"], [0.8, 0.95])
        self.assertEqual(
            [simulation.total_atoms for simulation in result["simulations"]],
            [1000, 10000],
        )

    def test_validate_simulations_rejects_duplicate_tags(self) -> None:
        systems = copy.deepcopy(self.counts_input()["simulations"])
        systems[1]["tag"] = systems[0]["tag"]

        with self.assertRaises(InputSchemaError):
            self.parser._validate_simulations(systems, "counts")

    def test_validate_simulations_rejects_total_atoms_in_counts_mode(self) -> None:
        systems = copy.deepcopy(self.counts_input()["simulations"])
        systems[0]["total_atoms"] = 1000

        with self.assertRaises(InputSchemaError):
            self.parser._validate_simulations(systems, "counts")

    def test_validate_simulations_rejects_different_ratios_between_systems(
        self,
    ) -> None:
        systems = copy.deepcopy(self.ratio_input()["simulations"])
        systems[1]["monomer_ratios"]["water"] = 2.0

        with self.assertRaises(InputSchemaError):
            self.parser._validate_simulations(systems, "ratio")


class TestSystemMonomerKeys(InputParserTestCase):
    """Tests for agreement between monomer definitions and systems."""

    def test_validate_system_monomer_keys_accepts_exact_names(self) -> None:
        inputs = self.counts_input()
        systems = inputs["simulations"]

        self.assertIsNone(
            self.parser._validate_system_monomer_keys(
                inputs,
                systems,
                "counts",
            )
        )

    def test_validate_system_monomer_keys_assigns_default_name(self) -> None:
        inputs = {
            "monomers": [
                {
                    "smiles": "CCO",
                }
            ]
        }
        systems = [
            {
                "tag": "a",
                "monomer_counts": {
                    "data_1": 1,
                },
            }
        ]

        self.parser._validate_system_monomer_keys(
            inputs,
            systems,
            "counts",
        )

        self.assertEqual(inputs["monomers"][0]["name"], "data_1")

    def test_validate_system_monomer_keys_rejects_unknown_name(self) -> None:
        inputs = self.counts_input()
        systems = copy.deepcopy(inputs["simulations"])
        systems[0]["monomer_counts"]["unknown"] = 1

        with self.assertRaises(InputSchemaError):
            self.parser._validate_system_monomer_keys(
                inputs,
                systems,
                "counts",
            )

    def test_validate_system_monomer_keys_rejects_missing_name(self) -> None:
        inputs = self.counts_input()
        systems = copy.deepcopy(inputs["simulations"])
        systems[0]["monomer_counts"].pop("water")

        with self.assertRaises(InputSchemaError):
            self.parser._validate_system_monomer_keys(
                inputs,
                systems,
                "counts",
            )


class TestMonomerEntryValidation(InputParserTestCase):
    """Tests for construction of MonomerEntry objects."""

    def test_validate_monomer_entry_builds_counts_map(self) -> None:
        inputs = self.counts_input()
        systems = copy.deepcopy(inputs["simulations"])

        monomers = self.parser._validate_monomer_entry(
            inputs,
            "counts",
            systems,
        )

        self.assertEqual(len(monomers), 2)
        self.assertTrue(
            all(isinstance(monomer, MonomerEntry) for monomer in monomers)
        )
        self.assertEqual(
            monomers[0].count,
            {
                "small": 2,
                "large": 20,
            },
        )
        self.assertIsNone(monomers[0].ratio)
        self.assertEqual(monomers[0].smiles, "CCO")
        self.assertIsInstance(monomers[0].rdkit_mol, Chem.Mol)
        self.assertGreater(monomers[0].num_atoms, 0)
        self.assertGreater(monomers[0].molecular_weight, 0.0)

    def test_validate_monomer_entry_builds_ratio_value(self) -> None:
        inputs = self.ratio_input()
        systems = copy.deepcopy(inputs["simulations"])

        monomers = self.parser._validate_monomer_entry(
            inputs,
            "ratio",
            systems,
        )

        self.assertEqual(monomers[0].ratio, 2.0)
        self.assertIsNone(monomers[0].count)

    def test_validate_monomer_entry_rejects_non_list_monomers(self) -> None:
        inputs = self.counts_input()
        inputs["monomers"] = {}

        with self.assertRaises(InputSchemaError):
            self.parser._validate_monomer_entry(
                inputs,
                "counts",
                inputs["simulations"],
            )

    def test_validate_monomer_entry_rejects_non_dictionary_entry(self) -> None:
        inputs = self.counts_input()
        inputs["monomers"] = ["CCO"]

        with self.assertRaises(InputSchemaError):
            self.parser._validate_monomer_entry(
                inputs,
                "counts",
                inputs["simulations"],
            )

    def test_validate_monomer_entry_rejects_negative_count(self) -> None:
        inputs = self.counts_input()
        inputs["simulations"][0]["monomer_counts"]["ethanol"] = -1

        with self.assertRaises(NumericFieldError):
            self.parser._validate_monomer_entry(
                inputs,
                "counts",
                inputs["simulations"],
            )

    def test_validate_monomer_entry_rejects_negative_ratio(self) -> None:
        inputs = self.ratio_input()
        inputs["simulations"][0]["monomer_ratios"]["ethanol"] = -1.0

        with self.assertRaises(NumericFieldError):
            self.parser._validate_monomer_entry(
                inputs,
                "ratio",
                inputs["simulations"],
            )


class TestLegacyNumericFields(InputParserTestCase):
    """Tests for the retained legacy numeric validator."""

    def test_validate_numeric_fields_accepts_valid_values(self) -> None:
        inputs = {
            "density": 0.8,
            "temperature": [300, 350.0],
            "number_of_monomers": {
                "data_1": 1,
                "data_2": 2,
            },
        }

        self.assertIsNone(self.parser._validate_numeric_fields(inputs))

    def test_validate_numeric_fields_rejects_invalid_density(self) -> None:
        inputs = {
            "density": 0,
            "temperature": 300,
            "number_of_monomers": {"data_1": 1},
        }

        with self.assertRaises(NumericFieldError):
            self.parser._validate_numeric_fields(inputs)

    def test_validate_numeric_fields_rejects_invalid_temperature(self) -> None:
        inputs = {
            "density": 0.8,
            "temperature": [300, False],
            "number_of_monomers": {"data_1": 1},
        }

        with self.assertRaises(NumericFieldError):
            self.parser._validate_numeric_fields(inputs)

    def test_validate_numeric_fields_rejects_invalid_monomer_count(self) -> None:
        inputs = {
            "density": 0.8,
            "temperature": 300,
            "number_of_monomers": {"data_1": 0},
        }

        with self.assertRaises(NumericFieldError):
            self.parser._validate_numeric_fields(inputs)


class TestLoopValidation(InputParserTestCase):
    """Tests for loop configuration."""

    def test_validate_loop_defaults_to_enabled(self) -> None:
        self.assertEqual(self.parser._validate_loop({}), (True, None))

    def test_validate_loop_accepts_booleans(self) -> None:
        self.assertEqual(
            self.parser._validate_loop({"loop": True}),
            (True, None),
        )
        self.assertEqual(
            self.parser._validate_loop({"loop": False}),
            (False, None),
        )

    def test_validate_loop_accepts_positive_integer(
        self,
    ) -> None:
        self.assertEqual(
            self.parser._validate_loop({"loop": 4}),
            (True, 4),
        )

    def test_validate_loop_accepts_supported_keywords(self) -> None:
        for keyword in ("loop", "repeat", "iterations", "do_loop"):
            with self.subTest(keyword=keyword):
                self.assertEqual(
                    self.parser._validate_loop({"loop": keyword}),
                    (True, None),
                )

    def test_validate_loop_rejects_invalid_values(self) -> None:
        for value in (0, -1, 1.5, None, "yes", [], {}):
            with self.subTest(value=value):
                with self.assertRaises(InputSchemaError):
                    self.parser._validate_loop({"loop": value})


class TestValidateInputsWorkflow(InputParserTestCase):
    """End-to-end tests for the public validation entry point."""

    def test_validate_inputs_builds_counts_setup(self) -> None:
        inputs = self.counts_input()

        result = self.parser.validate_inputs(inputs)

        self.assertIsInstance(result, SimulationSetup)
        self.assertEqual(result.simulation_name, "counts_example")
        self.assertEqual(result.composition_method, "counts")
        self.assertEqual(result.force_field, "PCFF")
        self.assertFalse(result.loop)
        self.assertIsNone(result.max_loop_count)
        self.assertEqual(result.temperature, [300.0, 350.0])
        self.assertEqual(result.density, [0.8, 0.95])
        self.assertEqual(len(result.monomers), 2)
        self.assertEqual(len(result.simulations), 2)
        self.assertEqual(
            result.monomers[0].count,
            {
                "small": 2,
                "large": 20,
            },
        )

    def test_validate_inputs_builds_ratio_setup(self) -> None:
        inputs = self.ratio_input()

        result = self.parser.validate_inputs(inputs)

        self.assertIsInstance(result, SimulationSetup)
        self.assertEqual(result.simulation_name, "ratio_example")
        self.assertEqual(result.composition_method, "ratio")
        self.assertEqual(result.force_field, "PCFF-IFF")
        self.assertTrue(result.loop)
        self.assertEqual(
            [simulation.total_atoms for simulation in result.simulations],
            [1000, 10000],
        )
        self.assertEqual(result.monomers[0].ratio, 2.0)
        self.assertIsNone(result.monomers[0].count)

    @patch("AutoREACTER.input_parser.time.sleep")
    def test_validate_inputs_stores_integer_loop_limit(
        self,
        sleep_mock,
    ) -> None:
        inputs = self.counts_input()
        inputs["loop"] = 3

        result = self.parser.validate_inputs(inputs)

        self.assertTrue(result.loop)
        self.assertEqual(result.max_loop_count, 3)
        sleep_mock.assert_called_once_with(5)

    def test_validate_inputs_assigns_default_monomer_name(self) -> None:
        inputs = {
            "simulation_name": "default_name",
            "simulations": [
                {
                    "tag": "a",
                    "temperature": 300,
                    "density": 0.8,
                    "monomer_counts": {
                        "data_1": 1,
                    },
                }
            ],
            "monomers": [
                {
                    "smiles": "CCO",
                }
            ],
        }

        result = self.parser.validate_inputs(inputs)

        self.assertEqual(result.monomers[0].name, "data_1")
        self.assertEqual(result.monomers[0].data_id, "data_1")

    def test_validate_inputs_rejects_unknown_system_monomer(self) -> None:
        inputs = self.counts_input()
        inputs["simulations"][0]["monomer_counts"]["unknown"] = 1

        with self.assertRaises(InputSchemaError):
            self.parser.validate_inputs(inputs)

    def test_validate_inputs_rejects_missing_system_monomer(self) -> None:
        inputs = self.counts_input()
        inputs["simulations"][0]["monomer_counts"].pop("water")

        with self.assertRaises(InputSchemaError):
            self.parser.validate_inputs(inputs)

    def test_validate_inputs_rejects_canonically_duplicate_smiles(self) -> None:
        inputs = {
            "simulation_name": "duplicates",
            "simulations": [
                {
                    "tag": "a",
                    "temperature": 300,
                    "density": 0.8,
                    "monomer_counts": {
                        "ethanol_a": 1,
                        "ethanol_b": 1,
                    },
                }
            ],
            "monomers": [
                {
                    "name": "ethanol_a",
                    "smiles": "CCO",
                },
                {
                    "name": "ethanol_b",
                    "smiles": "C(C)O",
                },
            ],
        }

        with self.assertRaises(DuplicateMonomerError):
            self.parser.validate_inputs(inputs)


class TestMoleculeRepresentation(InputParserTestCase):
    """Tests for initial-molecule extraction and image-grid integration."""

    def test_molecule_representation_returns_molecules_and_legends(self) -> None:
        setup = self.parser.validate_inputs(self.counts_input())

        molecules, legends = (
            self.parser.molecule_representation_of_initial_molecules(setup)
        )

        self.assertEqual(len(molecules), 2)
        self.assertTrue(all(isinstance(mol, Chem.Mol) for mol in molecules))
        self.assertEqual(legends, ["ethanol", "water"])

    @patch("AutoREACTER.input_parser.Draw.MolsToGridImage")
    def test_initial_molecules_image_grid_delegates_to_rdkit(
        self,
        grid_mock,
    ) -> None:
        setup = self.parser.validate_inputs(self.counts_input())

        class DummySession:
            inputs = setup

        expected_image = object()
        grid_mock.return_value = expected_image

        result = self.parser.initial_molecules_image_grid(DummySession())

        self.assertIs(result, expected_image)
        grid_mock.assert_called_once()

        args, kwargs = grid_mock.call_args
        self.assertEqual(len(args[0]), 2)
        self.assertEqual(kwargs["molsPerRow"], 3)
        self.assertEqual(kwargs["subImgSize"], (400, 400))
        self.assertEqual(kwargs["legends"], ["ethanol", "water"])


if __name__ == "__main__":
    unittest.main(verbosity=2)