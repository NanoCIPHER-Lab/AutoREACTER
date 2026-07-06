import unittest
#this should import the main class from the code
from AutoREACTER.input_parser import (
    InputParser,
    NumericFieldError,
    InputSchemaError,
    CompatibilityError,
    DuplicateMonomerError,
    InputConflictError,
    SmilesValidationError
)

class TestInputParser(unittest.TestCase) :
    
    def setUp(self):
        """Runs automatically before every test to give us a fresh parser."""
        self.parser = InputParser()
    #===============
    # Tests for:  _validate_temperature
    #===============

    def test_validate_temperature_valid(self):
        """Good Ending: A valid positive temp should return as a float."""
        result = self.parser._validate_temperature(300)
        self.assertEqual(result, 300.0)
    
    def test_validate_temperature_negative(self):
        """Bad Ending: A negative temp has to raise a NumericFieldError."""
        with self.assertRaises(NumericFieldError):
            self.parser._validate_temperature(-150)

    def test_validate_temperature_zero(self):
        """Bad Ending: Absolute zero has to raise a NumericFieldError."""
        with self.assertRaises(NumericFieldError):
            self.parser._validate_temperature(0)

    def test_validate_temperature_wrong_type(self):
        """Bad Ending: If a string or boolean is given instead of a number it has to fail."""
        with self.assertRaises(NumericFieldError):
            self.parser._validate_temperature("room_temperature")

    #===============
    # Tests for:  _validate_density
    #===============   
    def test_validate_density_valid(self):
        """Good Ending: Valid positive density returns as float."""
        result = self.parser._validate_density(0.85)
        self.assertEqual(result, 0.85)
    
    def test_validate_density_negative_or_zero(self):
        """Bad Ending: Rejects zero or negative density values."""
        with self.assertRaises(NumericFieldError):
            self.parser._validate_density(0)
        with self.assertRaises(NumericFieldError):
            self.parser._validate_density(-0.5)
    
    #================
    # Tests for: _validate_force_field
    #================

    def test_validate_force_field_default(self):
        """Good Ending: If None then defaults to 'PCFF'."""
        result = self.parser._validate_force_field(None)
        self.assertEqual(result, "PCFF")
    
    def test_validate_force_field_canonical(self):
        """Good path: Normalizes force-field capitalization and aliases."""
        result = self.parser._validate_force_field("pcff-iff")
        self.assertEqual(result, "PCFF-IFF")

    def test_validate_force_field_unsupported (self):
        """Bad Ending: Passing a completely random name should raise an InputSchemaError."""
        with self.assertRaises(InputSchemaError):
            self.parser._validate_force_field("NotAForceField")

    def test_validate_force_field_incompatible(self):
        """Bad Ending: 'OPLSAA' is recognized but incompatible, so it should raise a CompatibilityError."""
        with self.assertRaises(CompatibilityError):
            self.parser._validate_force_field("oplsaa")
    
    # ==========================================
    # Tests for: validate_no_duplicate_smiles
    # ==========================================

    def test_validate_no_duplicate_smiles_good_path(self):
        """Good Ending: Adding a unique SMILES appends it to the tracker list."""
        tracker_list = ["CCO", "C1=CC(=CC(=C1)N)N"]
        new_smiles = "ClC(=O)c1cc(cc(c1)C(Cl)=O)C(Cl)=O"
        
        result = self.parser.validate_no_duplicate_smiles(new_smiles, tracker_list)
        
        # Verify the list grew to 3 items and includes this new molecule
        self.assertEqual(len(result), 3)
        self.assertIn(new_smiles, result)

    def test_validate_no_duplicate_smiles_raises_error(self):
        """Bad Ending: Adding an existing SMILES triggers DuplicateMonomerError."""
        tracker_list = ["CCO", "C1=CC(=CC(=C1)N)N"]
        duplicate_smiles = "CCO"  # Already present!
        
        with self.assertRaises(DuplicateMonomerError):
            self.parser.validate_no_duplicate_smiles(duplicate_smiles, tracker_list)

    # ==========================================
    # Tests for: _validate_smiles
    # ==========================================

    def test_validate_smiles_valid(self):
        """Good Ending: A valid SMILES should return canonical SMILES and an RDKit Mol."""
        smiles, mol = self.parser._validate_smiles("CCO")

        self.assertEqual(smiles, "CCO")
        self.assertIsNotNone(mol)

    def test_validate_smiles_empty_raises_error(self):
        """Bad Ending: Empty SMILES should raise SmilesValidationError."""
        with self.assertRaises(SmilesValidationError):
            self.parser._validate_smiles("")

    def test_validate_smiles_invalid_raises_error(self):
        """Bad Ending: Invalid SMILES should raise SmilesValidationError."""
        with self.assertRaises(SmilesValidationError):
            self.parser._validate_smiles("not_a_smiles")

    # ==========================================
    # Tests for: validate_basic_format
    # ==========================================

    def test_validate_basic_format_valid(self):
        """Good Ending: A minimal valid top-level input should pass basic format validation."""
        inputs = {
            "simulation_name": "test_sim",
            "simulations": [],
            "monomers": [],
        }

        self.assertIsNone(self.parser.validate_basic_format(inputs))

    def test_validate_basic_format_missing_key(self):
        """Bad Ending: Missing required top-level keys should raise InputSchemaError."""
        inputs = {
            "simulation_name": "test_sim",
            "simulations": [],
        }

        with self.assertRaises(InputSchemaError):
            self.parser.validate_basic_format(inputs)

    def test_validate_basic_format_wrong_type(self):
        """Bad Ending: Non-dictionary input should raise InputSchemaError."""
        with self.assertRaises(InputSchemaError):
            self.parser.validate_basic_format(["not", "a", "dict"])

    # ==========================================
    # Tests for: _get_inputs_mode
    # ==========================================

    def test_get_inputs_mode_counts(self):
        """Good Ending: Simulations using monomer_counts should be detected as counts mode."""
        simulations = [
            {
                "tag": "10k",
                "temperature": 300,
                "density": 0.8,
                "monomer_counts": {"tmc": 1},
            }
        ]

        result = self.parser._get_inputs_mode(simulations)

        self.assertEqual(result, "counts")

    def test_get_inputs_mode_ratio(self):
        """Good Ending: Simulations using monomer_ratios should be detected as ratio mode."""
        simulations = [
            {
                "tag": "10k",
                "temperature": 300,
                "density": 0.8,
                "total_atoms": 10000,
                "monomer_ratios": {"tmc": 1.0},
            }
        ]

        result = self.parser._get_inputs_mode(simulations)

        self.assertEqual(result, "ratio")

    def test_get_inputs_mode_mixed_modes_raises_error(self):
        """Bad Ending: Mixing counts and ratio modes should raise InputConflictError."""
        simulations = [
            {
                "tag": "counts_system",
                "temperature": 300,
                "density": 0.8,
                "monomer_counts": {"tmc": 1},
            },
            {
                "tag": "ratio_system",
                "temperature": 300,
                "density": 0.8,
                "total_atoms": 10000,
                "monomer_ratios": {"tmc": 1.0},
            },
        ]

        with self.assertRaises(InputConflictError):
            self.parser._get_inputs_mode(simulations)

    def test_get_inputs_mode_both_counts_and_ratios_raises_error(self):
        """Bad Ending: One simulation cannot contain both counts and ratios."""
        simulations = [
            {
                "tag": "bad_system",
                "temperature": 300,
                "density": 0.8,
                "monomer_counts": {"tmc": 1},
                "monomer_ratios": {"tmc": 1.0},
            }
        ]

        with self.assertRaises(InputConflictError):
            self.parser._get_inputs_mode(simulations)

    # ==========================================
    # Tests for: validate_inputs
    # ==========================================

    def test_validate_inputs_counts_mode_valid(self):
        """Good Ending: Full valid counts-mode input should produce a SimulationSetup object."""
        inputs = {
            "simulation_name": "test_counts",
            "simulations": [
                {
                    "tag": "10k",
                    "temperature": 300,
                    "density": 0.8,
                    "monomer_counts": {
                        "tmc": 1,
                        "mpd": 1,
                    },
                }
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
            ],
        }

        result = self.parser.validate_inputs(inputs)

        self.assertEqual(result.simulation_name, "test_counts")
        self.assertEqual(result.composition_method, "counts")
        self.assertEqual(result.force_field, "PCFF")
        self.assertEqual(len(result.monomers), 2)
        self.assertEqual(len(result.simulations), 1)

    def test_validate_inputs_ratio_mode_valid(self):
        """Good Ending: Full valid ratio-mode input should produce a SimulationSetup object."""
        inputs = {
            "simulation_name": "test_ratio",
            "simulations": [
                {
                    "tag": "10k",
                    "temperature": 300,
                    "density": 0.8,
                    "total_atoms": 10000,
                    "monomer_ratios": {
                        "tmc": 1.0,
                        "mpd": 1.0,
                    },
                }
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
            ],
        }

        result = self.parser.validate_inputs(inputs)

        self.assertEqual(result.simulation_name, "test_ratio")
        self.assertEqual(result.composition_method, "ratio")
        self.assertEqual(result.force_field, "PCFF")
        self.assertEqual(len(result.monomers), 2)
        self.assertEqual(len(result.simulations), 1)

    def test_validate_inputs_unknown_monomer_in_system_raises_error(self):
        """Bad Ending: System composition cannot reference undefined monomers."""
        inputs = {
            "simulation_name": "bad_unknown_monomer",
            "simulations": [
                {
                    "tag": "10k",
                    "temperature": 300,
                    "density": 0.8,
                    "monomer_counts": {
                        "tmc": 1,
                        "unknown": 1,
                    },
                }
            ],
            "monomers": [
                {
                    "name": "tmc",
                    "smiles": "ClC(=O)c1cc(cc(c1)C(Cl)=O)C(Cl)=O",
                },
            ],
        }

        with self.assertRaises(InputSchemaError):
            self.parser.validate_inputs(inputs)

    def test_validate_inputs_missing_monomer_in_system_raises_error(self):
        """Bad Ending: Every defined monomer must appear in each system composition."""
        inputs = {
            "simulation_name": "bad_missing_monomer",
            "simulations": [
                {
                    "tag": "10k",
                    "temperature": 300,
                    "density": 0.8,
                    "monomer_counts": {
                        "tmc": 1,
                    },
                }
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
            ],
        }

        with self.assertRaises(InputSchemaError):
            self.parser.validate_inputs(inputs)

    def test_validate_inputs_duplicate_monomer_smiles_raises_error(self):
        """Bad Ending: Duplicate monomer SMILES should raise DuplicateMonomerError."""
        inputs = {
            "simulation_name": "bad_duplicate_smiles",
            "simulations": [
                {
                    "tag": "10k",
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
                    "smiles": "CCO",
                },
            ],
        }

        with self.assertRaises(DuplicateMonomerError):
            self.parser.validate_inputs(inputs)


if __name__ == "__main__":
    unittest.main() 
