"""
Comprehensive tests for the input_parser.py module.

Tests cover:
- Input validation
- SMILES validation with RDKit
- Duplicate SMILES detection
- Component checking
- Error handling and edge cases
"""
import pytest
from unittest.mock import patch, MagicMock

from reaction_lammps_mupt.input_parser import InputParser


class TestInputParserInitialization:
    """Tests for InputParser initialization."""

    def test_init_with_valid_inputs(self):
        """Test initialization with valid inputs."""
        inputs = {
            "simulation_name": "TestSim",
            "temperature": [300],
            "density": 0.8,
            "monomers": {1: "CCO", 2: "CCN"},
            "number_of_monomers": {1: 100, 2: 100},
        }
        parser = InputParser(inputs)
        assert parser.inputs == inputs
        assert parser.validated_inputs == inputs

    def test_init_calls_validate_inputs(self):
        """Test that initialization calls validate_inputs."""
        inputs = {
            "simulation_name": "TestSim",
            "temperature": [300],
            "density": 0.8,
            "monomers": {1: "CCO", 2: "CCN"},
            "number_of_monomers": {1: 100, 2: 100},
        }
        with patch.object(InputParser, "validate_inputs") as mock_validate:
            mock_validate.return_value = inputs
            parser = InputParser(inputs)
            mock_validate.assert_called_once_with(inputs)


class TestComponentCheck:
    """Tests for component_check method."""

    def test_component_check_all_required_keys_present(self):
        """Test component_check with all required keys."""
        inputs = {
            "simulation_name": "TestSim",
            "temperature": [300, 400],
            "density": 0.8,
            "monomers": {1: "CCO", 2: "CCN"},
            "number_of_monomers": {1: 100, 2: 100},
        }
        parser = InputParser.__new__(InputParser)
        # Should not raise any exception
        parser.component_check(inputs)

    def test_component_check_missing_simulation_name(self):
        """Test component_check with missing simulation_name."""
        inputs = {
            "temperature": [300],
            "density": 0.8,
            "monomers": {1: "CCO"},
            "number_of_monomers": {1: 100},
        }
        parser = InputParser.__new__(InputParser)
        with pytest.raises(KeyError, match="Missing required key: 'simulation_name'"):
            parser.component_check(inputs)

    def test_component_check_missing_temperature(self):
        """Test component_check with missing temperature."""
        inputs = {
            "simulation_name": "Test",
            "density": 0.8,
            "monomers": {1: "CCO"},
            "number_of_monomers": {1: 100},
        }
        parser = InputParser.__new__(InputParser)
        with pytest.raises(KeyError, match="Missing required key: 'temperature'"):
            parser.component_check(inputs)

    def test_component_check_missing_density(self):
        """Test component_check with missing density."""
        inputs = {
            "simulation_name": "Test",
            "temperature": [300],
            "monomers": {1: "CCO"},
            "number_of_monomers": {1: 100},
        }
        parser = InputParser.__new__(InputParser)
        with pytest.raises(KeyError, match="Missing required key: 'density'"):
            parser.component_check(inputs)

    def test_component_check_missing_monomers(self):
        """Test component_check with missing monomers."""
        inputs = {
            "simulation_name": "Test",
            "temperature": [300],
            "density": 0.8,
            "number_of_monomers": {1: 100},
        }
        parser = InputParser.__new__(InputParser)
        with pytest.raises(KeyError, match="Missing required key: 'monomers'"):
            parser.component_check(inputs)

    def test_component_check_missing_number_of_monomers(self):
        """Test component_check with missing number_of_monomers."""
        inputs = {
            "simulation_name": "Test",
            "temperature": [300],
            "density": 0.8,
            "monomers": {1: "CCO"},
        }
        parser = InputParser.__new__(InputParser)
        with pytest.raises(KeyError, match="Missing required key: 'number_of_monomers'"):
            parser.component_check(inputs)

    def test_component_check_monomers_not_dict(self):
        """Test component_check when monomers is not a dict."""
        inputs = {
            "simulation_name": "Test",
            "temperature": [300],
            "density": 0.8,
            "monomers": ["CCO", "CCN"],  # Should be dict
            "number_of_monomers": {1: 100},
        }
        parser = InputParser.__new__(InputParser)
        with pytest.raises(ValueError, match="'monomers' key must map to a non-empty dictionary"):
            parser.component_check(inputs)

    def test_component_check_monomers_empty_dict(self):
        """Test component_check when monomers is empty dict."""
        inputs = {
            "simulation_name": "Test",
            "temperature": [300],
            "density": 0.8,
            "monomers": {},
            "number_of_monomers": {},
        }
        parser = InputParser.__new__(InputParser)
        with pytest.raises(ValueError, match="'monomers' key must map to a non-empty dictionary"):
            parser.component_check(inputs)

    def test_component_check_number_of_monomers_not_dict(self):
        """Test component_check when number_of_monomers is not a dict."""
        inputs = {
            "simulation_name": "Test",
            "temperature": [300],
            "density": 0.8,
            "monomers": {1: "CCO"},
            "number_of_monomers": [100],  # Should be dict
        }
        parser = InputParser.__new__(InputParser)
        with pytest.raises(ValueError, match="'number_of_monomers' key must map to a non-empty dictionary"):
            parser.component_check(inputs)

    def test_component_check_length_mismatch(self):
        """Test component_check when monomers and number_of_monomers have different lengths."""
        inputs = {
            "simulation_name": "Test",
            "temperature": [300],
            "density": 0.8,
            "monomers": {1: "CCO", 2: "CCN", 3: "CCC"},
            "number_of_monomers": {1: 100, 2: 100},  # Length mismatch
        }
        parser = InputParser.__new__(InputParser)
        with pytest.raises(ValueError, match="Length mismatch"):
            parser.component_check(inputs)


class TestSMILESValidation:
    """Tests for SMILES validation methods."""

    def test_validate_smiles_valid(self):
        """Test _validate_smiles with valid SMILES."""
        parser = InputParser.__new__(InputParser)
        # Should not raise exception
        parser._validate_smiles("CCO")  # Ethanol
        parser._validate_smiles("c1ccccc1")  # Benzene
        parser._validate_smiles("CC(=O)O")  # Acetic acid

    def test_validate_smiles_invalid(self):
        """Test _validate_smiles with invalid SMILES."""
        parser = InputParser.__new__(InputParser)
        with pytest.raises(ValueError, match="Invalid SMILES string"):
            parser._validate_smiles("INVALID_SMILES")

    def test_validate_smiles_rdkit_all_valid(self):
        """Test validate_smiles_rdkit with all valid SMILES."""
        inputs = {
            "monomers": {
                1: "CCO",
                2: "CCN",
                3: "c1ccccc1",
            }
        }
        parser = InputParser.__new__(InputParser)
        # Should not raise exception
        parser.validate_smiles_rdkit(inputs)

    def test_validate_smiles_rdkit_missing_monomers_key(self):
        """Test validate_smiles_rdkit with missing monomers key."""
        inputs = {}
        parser = InputParser.__new__(InputParser)
        with pytest.raises(KeyError, match='Missing required key: "monomers"'):
            parser.validate_smiles_rdkit(inputs)

    def test_validate_smiles_rdkit_invalid_smiles(self):
        """Test validate_smiles_rdkit with invalid SMILES string."""
        inputs = {
            "monomers": {
                1: "CCO",
                2: "INVALID_SMILES",
            }
        }
        parser = InputParser.__new__(InputParser)
        with pytest.raises(ValueError, match="Invalid SMILES string: 'INVALID_SMILES'"):
            parser.validate_smiles_rdkit(inputs)

    def test_validate_smiles_rdkit_empty_smiles(self):
        """Test validate_smiles_rdkit with empty SMILES string."""
        inputs = {
            "monomers": {
                1: "CCO",
                2: "",
            }
        }
        parser = InputParser.__new__(InputParser)
        with pytest.raises(ValueError, match="SMILES must be a non-empty string"):
            parser.validate_smiles_rdkit(inputs)

    def test_validate_smiles_rdkit_whitespace_only(self):
        """Test validate_smiles_rdkit with whitespace-only SMILES."""
        inputs = {
            "monomers": {
                1: "CCO",
                2: "   ",
            }
        }
        parser = InputParser.__new__(InputParser)
        with pytest.raises(ValueError, match="SMILES must be a non-empty string"):
            parser.validate_smiles_rdkit(inputs)

    def test_validate_smiles_rdkit_not_string(self):
        """Test validate_smiles_rdkit when SMILES is not a string."""
        inputs = {
            "monomers": {
                1: "CCO",
                2: 12345,  # Not a string
            }
        }
        parser = InputParser.__new__(InputParser)
        with pytest.raises(ValueError, match="SMILES must be a non-empty string"):
            parser.validate_smiles_rdkit(inputs)


class TestDuplicateDetection:
    """Tests for duplicate SMILES detection."""

    def test_validate_no_duplicate_smiles_all_unique(self):
        """Test validate_no_duplicate_smiles with all unique SMILES."""
        inputs = {
            "monomers": {
                1: "CCO",
                2: "CCN",
                3: "CCC",
            }
        }
        parser = InputParser.__new__(InputParser)
        # Should not raise exception
        parser.validate_no_duplicate_smiles(inputs)

    def test_validate_no_duplicate_smiles_duplicates(self):
        """Test validate_no_duplicate_smiles with duplicate SMILES."""
        inputs = {
            "monomers": {
                1: "CCO",
                2: "CCN",
                3: "CCO",  # Duplicate
            }
        }
        parser = InputParser.__new__(InputParser)
        with pytest.raises(ValueError, match="Duplicate SMILES detected"):
            parser.validate_no_duplicate_smiles(inputs)

    def test_validate_no_duplicate_smiles_case_sensitive(self):
        """Test that duplicate detection is case-sensitive."""
        inputs = {
            "monomers": {
                1: "CCO",
                2: "cco",  # Different case
            }
        }
        parser = InputParser.__new__(InputParser)
        # Should not raise exception (case-sensitive comparison)
        parser.validate_no_duplicate_smiles(inputs)

    def test_validate_no_duplicate_smiles_empty_monomers(self):
        """Test validate_no_duplicate_smiles with empty monomers dict."""
        inputs = {"monomers": {}}
        parser = InputParser.__new__(InputParser)
        # Should not raise exception
        parser.validate_no_duplicate_smiles(inputs)

    def test_validate_no_duplicate_smiles_missing_monomers_key(self):
        """Test validate_no_duplicate_smiles with missing monomers key."""
        inputs = {}
        parser = InputParser.__new__(InputParser)
        # Should not raise exception (gracefully handles missing key)
        parser.validate_no_duplicate_smiles(inputs)

    def test_validate_no_duplicate_smiles_multiple_duplicates(self):
        """Test validate_no_duplicate_smiles with multiple duplicates."""
        inputs = {
            "monomers": {
                1: "CCO",
                2: "CCN",
                3: "CCO",  # Duplicate of 1
            }
        }
        parser = InputParser.__new__(InputParser)
        with pytest.raises(ValueError, match="Duplicate SMILES detected for monomers 1 and 3: 'CCO'"):
            parser.validate_no_duplicate_smiles(inputs)


class TestValidateInputs:
    """Tests for the main validate_inputs method."""

    def test_validate_inputs_all_checks_pass(self):
        """Test validate_inputs when all checks pass."""
        inputs = {
            "simulation_name": "TestSim",
            "temperature": [300, 400],
            "density": 0.8,
            "monomers": {1: "CCO", 2: "CCN"},
            "number_of_monomers": {1: 100, 2: 100},
        }
        parser = InputParser(inputs)
        assert parser.validated_inputs == inputs

    def test_validate_inputs_calls_all_validators(self):
        """Test that validate_inputs calls all validation methods."""
        inputs = {
            "simulation_name": "TestSim",
            "temperature": [300],
            "density": 0.8,
            "monomers": {1: "CCO"},
            "number_of_monomers": {1: 100},
        }
        parser = InputParser.__new__(InputParser)
        parser.inputs = inputs

        with patch.object(parser, "component_check") as mock_component, \
             patch.object(parser, "validate_smiles_rdkit") as mock_smiles, \
             patch.object(parser, "validate_no_duplicate_smiles") as mock_duplicates:

            result = parser.validate_inputs(inputs)

            mock_component.assert_called_once_with(inputs)
            mock_smiles.assert_called_once_with(inputs)
            mock_duplicates.assert_called_once_with(inputs)
            assert result == inputs

    def test_validate_inputs_fails_on_component_check(self):
        """Test validate_inputs fails when component_check fails."""
        inputs = {
            "simulation_name": "TestSim",
            # Missing other required keys
        }
        with pytest.raises(KeyError):
            InputParser(inputs)

    def test_validate_inputs_fails_on_smiles_validation(self):
        """Test validate_inputs fails when SMILES validation fails."""
        inputs = {
            "simulation_name": "TestSim",
            "temperature": [300],
            "density": 0.8,
            "monomers": {1: "INVALID"},
            "number_of_monomers": {1: 100},
        }
        with pytest.raises(ValueError, match="Invalid SMILES string"):
            InputParser(inputs)

    def test_validate_inputs_fails_on_duplicate_detection(self):
        """Test validate_inputs fails when duplicate SMILES detected."""
        inputs = {
            "simulation_name": "TestSim",
            "temperature": [300],
            "density": 0.8,
            "monomers": {1: "CCO", 2: "CCO"},
            "number_of_monomers": {1: 100, 2: 100},
        }
        with pytest.raises(ValueError, match="Duplicate SMILES detected"):
            InputParser(inputs)


class TestToDict:
    """Tests for the to_dict method."""

    def test_to_dict_returns_validated_inputs(self):
        """Test that to_dict returns validated_inputs."""
        inputs = {
            "simulation_name": "TestSim",
            "temperature": [300],
            "density": 0.8,
            "monomers": {1: "CCO"},
            "number_of_monomers": {1: 100},
        }
        parser = InputParser(inputs)
        result = parser.to_dict()
        assert result == parser.validated_inputs

    def test_to_dict_prints_success_message(self, capsys):
        """Test that to_dict prints success message."""
        inputs = {
            "simulation_name": "TestSim",
            "temperature": [300],
            "density": 0.8,
            "monomers": {1: "CCO"},
            "number_of_monomers": {1: 100},
        }
        parser = InputParser(inputs)
        parser.to_dict()

        captured = capsys.readouterr()
        assert "Validation successful" in captured.out


class TestComplexScenarios:
    """Tests for complex and edge case scenarios."""

    def test_multiple_monomers_valid(self):
        """Test validation with multiple monomers."""
        inputs = {
            "simulation_name": "MultiMonomer",
            "temperature": [300, 350, 400],
            "density": 1.0,
            "monomers": {
                1: "CCO",
                2: "CCN",
                3: "CCC",
                4: "c1ccccc1",
                5: "CC(=O)O",
            },
            "number_of_monomers": {
                1: 100,
                2: 200,
                3: 300,
                4: 400,
                5: 500,
            },
        }
        parser = InputParser(inputs)
        assert len(parser.validated_inputs["monomers"]) == 5

    def test_complex_smiles_structures(self):
        """Test validation with complex SMILES structures."""
        inputs = {
            "simulation_name": "Complex",
            "temperature": [300],
            "density": 0.8,
            "monomers": {
                1: "O=C(O)c1cc(O)cc(C(=O)O)c1",  # Aromatic with multiple functional groups
                2: "C1CCC(=O)OCC1",  # Cyclic structure
                3: "CC(C)(C)c1ccc(O)cc1",  # Branched aromatic
            },
            "number_of_monomers": {1: 100, 2: 100, 3: 100},
        }
        parser = InputParser(inputs)
        assert len(parser.validated_inputs["monomers"]) == 3

    def test_string_keys_for_monomers(self):
        """Test that string keys work for monomers dict."""
        inputs = {
            "simulation_name": "StringKeys",
            "temperature": [300],
            "density": 0.8,
            "monomers": {"1": "CCO", "2": "CCN"},
            "number_of_monomers": {"1": 100, "2": 100},
        }
        parser = InputParser(inputs)
        assert parser.validated_inputs is not None

    def test_additional_optional_keys(self):
        """Test that additional optional keys are preserved."""
        inputs = {
            "simulation_name": "TestSim",
            "temperature": [300],
            "density": 0.8,
            "monomers": {1: "CCO"},
            "number_of_monomers": {1: 100},
            "stoichiometric_ratio": {1: 1},  # Optional key
            "number_of_total_atoms": [10000],  # Optional key
        }
        parser = InputParser(inputs)
        assert "stoichiometric_ratio" in parser.validated_inputs
        assert "number_of_total_atoms" in parser.validated_inputs

    def test_very_large_monomer_dictionary(self):
        """Test validation with a large number of monomers."""
        num_monomers = 100
        inputs = {
            "simulation_name": "LargeScale",
            "temperature": [300],
            "density": 0.8,
            "monomers": {i: f"C{'C' * i}O" for i in range(1, num_monomers + 1)},
            "number_of_monomers": {i: 10 for i in range(1, num_monomers + 1)},
        }
        parser = InputParser(inputs)
        assert len(parser.validated_inputs["monomers"]) == num_monomers

    def test_unicode_in_simulation_name(self):
        """Test that unicode characters in simulation_name are handled."""
        inputs = {
            "simulation_name": "测试模拟",
            "temperature": [300],
            "density": 0.8,
            "monomers": {1: "CCO"},
            "number_of_monomers": {1: 100},
        }
        parser = InputParser(inputs)
        assert parser.validated_inputs["simulation_name"] == "测试模拟"

    def test_negative_density(self):
        """Test that negative density is accepted (validation not enforced)."""
        inputs = {
            "simulation_name": "NegativeDensity",
            "temperature": [300],
            "density": -0.8,  # Physically invalid but not checked
            "monomers": {1: "CCO"},
            "number_of_monomers": {1: 100},
        }
        # Should not raise exception (density validation not implemented)
        parser = InputParser(inputs)
        assert parser.validated_inputs["density"] == -0.8

    def test_empty_temperature_list(self):
        """Test with empty temperature list."""
        inputs = {
            "simulation_name": "EmptyTemp",
            "temperature": [],  # Empty but present
            "density": 0.8,
            "monomers": {1: "CCO"},
            "number_of_monomers": {1: 100},
        }
        # Should not raise exception (temperature list validation not enforced)
        parser = InputParser(inputs)
        assert parser.validated_inputs["temperature"] == []


class TestRegressionScenarios:
    """Tests for regression scenarios and boundary conditions."""

    def test_canonical_smiles_not_normalized(self):
        """Test that different representations of same molecule are treated as different."""
        inputs = {
            "simulation_name": "SMILES",
            "temperature": [300],
            "density": 0.8,
            "monomers": {
                1: "CCO",
                2: "OCC",  # Same molecule, different SMILES
            },
            "number_of_monomers": {1: 100, 2: 100},
        }
        # Should not raise exception (no canonical normalization)
        parser = InputParser(inputs)
        assert len(parser.validated_inputs["monomers"]) == 2

    def test_monomer_zero_count(self):
        """Test with zero count for monomers."""
        inputs = {
            "simulation_name": "ZeroCount",
            "temperature": [300],
            "density": 0.8,
            "monomers": {1: "CCO", 2: "CCN"},
            "number_of_monomers": {1: 0, 2: 100},  # Zero count
        }
        # Should not raise exception (count validation not enforced)
        parser = InputParser(inputs)
        assert parser.validated_inputs["number_of_monomers"][1] == 0

    def test_special_characters_in_smiles(self):
        """Test SMILES with special characters."""
        inputs = {
            "simulation_name": "SpecialChars",
            "temperature": [300],
            "density": 0.8,
            "monomers": {
                1: "C[C@@H](N)C(=O)O",  # Stereochemistry
                2: "[Na+].[Cl-]",  # Ionic
            },
            "number_of_monomers": {1: 100, 2: 100},
        }
        parser = InputParser(inputs)
        assert len(parser.validated_inputs["monomers"]) == 2

    def test_immutability_of_original_inputs(self):
        """Test that original inputs dict is not modified."""
        inputs = {
            "simulation_name": "Immutable",
            "temperature": [300],
            "density": 0.8,
            "monomers": {1: "CCO"},
            "number_of_monomers": {1: 100},
        }
        original_inputs = inputs.copy()
        parser = InputParser(inputs)

        # Original should remain unchanged
        assert inputs == original_inputs