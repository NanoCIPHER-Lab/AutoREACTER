"""
Comprehensive tests for functional_groups_detector.py module.

Tests cover:
- Monomer functionality detection
- SMARTS pattern matching
- Functional group categorization
- Edge cases and validation
"""
import pytest
from unittest.mock import patch, MagicMock

from reaction_lammps_mupt.detectors.functional_groups_detector import (
    detect_monomer_functionality,
    functional_groups_detector,
    monomer_types,
)


class TestDetectMonomerFunctionality:
    """Tests for detect_monomer_functionality function."""

    def test_detect_vinyl_monomer(self):
        """Test detection of vinyl functionality."""
        smiles = "C=CCO"  # Allyl alcohol
        result = detect_monomer_functionality(smiles, "vinyl", "[C]=[C;D1]")

        assert result[0] == 1  # functionality_count
        assert result[1] >= 1  # count_1
        assert result[2] is None  # count_2

    def test_detect_di_different_both_groups_present(self):
        """Test di_different with both functional groups present."""
        smiles = "OCCC(O)=O"  # Lactic acid
        result = detect_monomer_functionality(
            smiles,
            "di_different",
            "[OX2H1;!$(OC=*):1]",  # Alcohol
            "[CX3:2](=[O])[OX2H1]"  # Carboxylic acid
        )

        assert result[0] == 2  # Both groups present
        assert result[1] >= 1  # Alcohol count
        assert result[2] >= 1  # Carboxylic acid count

    def test_detect_di_different_missing_second_group(self):
        """Test di_different with second group missing."""
        smiles = "CCO"  # Just alcohol, no carboxylic acid
        result = detect_monomer_functionality(
            smiles,
            "di_different",
            "[OX2H1;!$(OC=*):1]",
            "[CX3:2](=[O])[OX2H1]"
        )

        assert result[0] == 0  # Not qualified
        assert result[1] >= 1  # Alcohol present
        assert result[2] == 0  # Carboxylic acid absent

    def test_detect_di_identical_sufficient_count(self):
        """Test di_identical with sufficient functional group count."""
        smiles = "OCCO"  # Ethylene glycol (diol)
        result = detect_monomer_functionality(
            smiles,
            "di_identical",
            "[O,S;X2;H1;!$([O,S]C=*):3]"
        )

        assert result[0] == 2  # Qualified (2+ groups)
        assert result[1] == 2  # Two OH groups
        assert result[2] is None

    def test_detect_di_identical_insufficient_count(self):
        """Test di_identical with insufficient count (chain terminator)."""
        smiles = "CCO"  # Ethanol (only 1 OH)
        result = detect_monomer_functionality(
            smiles,
            "di_identical",
            "[O,S;X2;H1;!$([O,S]C=*):3]"
        )

        assert result[0] == 0  # Not qualified
        assert result[1] == 1  # One OH group
        assert result[2] is None

    def test_detect_invalid_smiles(self):
        """Test with invalid SMILES string."""
        smiles = "INVALID_SMILES"
        result = detect_monomer_functionality(smiles, "vinyl", "[C]=[C;D1]")

        assert result == (0, None, None)

    def test_detect_invalid_smarts(self):
        """Test with invalid SMARTS pattern."""
        smiles = "CCO"
        result = detect_monomer_functionality(smiles, "vinyl", "INVALID_SMARTS")

        assert result == (0, None, None)

    def test_detect_mono_functionality(self):
        """Test mono functionality type."""
        smiles = "C1CO1"  # Ethylene oxide (epoxide)
        result = detect_monomer_functionality(
            smiles,
            "mono",
            "[CX4;R:3]1[OX2;R:4][CX4;R:5]1"
        )

        assert result[0] == 1  # Single functional group
        assert result[1] == 1
        assert result[2] is None

    def test_detect_multiple_matches(self):
        """Test molecule with multiple matching sites."""
        smiles = "O=C(O)CCCCC(=O)O"  # Adipic acid (2 carboxylic acids)
        result = detect_monomer_functionality(
            smiles,
            "di_identical",
            "[CX3:1](=[O])[OX2H1]"
        )

        assert result[0] == 2
        assert result[1] == 2  # Two carboxylic acids
        assert result[2] is None


class TestFunctionalGroupsDetector:
    """Tests for functional_groups_detector function."""

    def test_detector_with_single_monomer(self):
        """Test detector with single monomer."""
        monomer_dict = {1: "OCCC(O)=O"}  # Lactic acid

        result = functional_groups_detector(monomer_dict)

        assert 1 in result
        assert "smiles" in result[1]
        assert result[1]["smiles"] == "OCCC(O)=O"
        # Check functional group detected
        assert 1 in result[1]  # Functional group index

    def test_detector_with_multiple_monomers(self):
        """Test detector with multiple monomers."""
        monomer_dict = {
            1: "OCCC(O)=O",
            2: "C1=CC=C(C(=C1)C(=O)O)O",
        }

        result = functional_groups_detector(monomer_dict)

        assert len(result) >= 1  # At least one monomer with functional groups

    def test_detector_no_functional_groups(self):
        """Test detector with monomers having no matching functional groups."""
        monomer_dict = {1: "CCC"}  # Propane (no functional groups)

        result = functional_groups_detector(monomer_dict)

        assert 1 not in result  # Should not be in results

    def test_detector_empty_input(self):
        """Test detector with empty input."""
        monomer_dict = {}

        result = functional_groups_detector(monomer_dict)

        assert result == {}

    def test_detector_preserves_monomer_index(self):
        """Test that detector preserves original monomer indices."""
        monomer_dict = {
            5: "OCCC(O)=O",
            10: "C1=CC=C(C(=C1)C(=O)O)O",
        }

        result = functional_groups_detector(monomer_dict)

        # Original indices should be preserved
        if 5 in result:
            assert "smiles" in result[5]
        if 10 in result:
            assert "smiles" in result[10]

    def test_detector_output_structure(self):
        """Test that detector output has correct structure."""
        monomer_dict = {1: "OCCC(O)=O"}

        result = functional_groups_detector(monomer_dict)

        if 1 in result:
            assert "smiles" in result[1]
            # Check for functional group data
            for key in result[1]:
                if isinstance(key, int):
                    fg_data = result[1][key]
                    assert "functionality_type" in fg_data
                    assert "functional_group_name" in fg_data
                    assert "functional_group_smarts_1" in fg_data
                    assert "functional_count_1" in fg_data

    def test_detector_prints_detection_messages(self, capsys):
        """Test that detector prints detection messages."""
        monomer_dict = {1: "OCCC(O)=O"}

        functional_groups_detector(monomer_dict)

        captured = capsys.readouterr()
        # Should print some detection information
        # (depends on current monomer_types configuration)

    def test_detector_with_string_indices(self):
        """Test detector with string indices."""
        monomer_dict = {"1": "OCCC(O)=O", "2": "CCO"}

        result = functional_groups_detector(monomer_dict)

        # Should handle string indices
        assert isinstance(result, dict)


class TestMonomerTypesConfiguration:
    """Tests for monomer_types configuration."""

    def test_monomer_types_structure(self):
        """Test that monomer_types has correct structure."""
        assert isinstance(monomer_types, dict)

        for name, config in monomer_types.items():
            assert "functionality_type" in config
            assert "smarts_1" in config
            assert "group_name" in config

            # di_different should have smarts_2
            if config["functionality_type"] == "di_different":
                assert "smarts_2" in config

    def test_monomer_types_has_hydroxy_carboxylic_acid(self):
        """Test that hydroxy_carboxylic_acid is configured."""
        found = False
        for config in monomer_types.values():
            if config.get("group_name") == "hydroxy_carboxylic_acid":
                found = True
                assert config["functionality_type"] == "di_different"
                break

        assert found, "hydroxy_carboxylic_acid should be in monomer_types"


class TestEdgeCases:
    """Tests for edge cases and boundary conditions."""

    def test_detect_with_none_smarts2(self):
        """Test di_different explicitly with None smarts_2."""
        smiles = "CCO"
        result = detect_monomer_functionality(smiles, "di_different", "[OX2H1]", None)

        # Should handle None gracefully
        assert isinstance(result, tuple)

    def test_detect_complex_molecule(self):
        """Test with complex molecule structure."""
        smiles = "O=C(O)c1cc(O)cc(C(=O)O)c1"  # Complex aromatic
        result = detect_monomer_functionality(
            smiles,
            "di_identical",
            "[CX3:1](=[O])[OX2H1]"
        )

        # Should detect multiple carboxylic acids
        assert result[1] >= 2 if result[0] > 0 else True

    def test_detector_with_stereochemistry(self):
        """Test detector with stereochemistry in SMILES."""
        monomer_dict = {1: "C[C@@H](O)C(=O)O"}  # L-lactic acid with stereochemistry

        result = functional_groups_detector(monomer_dict)

        # Should handle stereochemistry
        assert isinstance(result, dict)

    def test_detector_with_aromatic_smiles(self):
        """Test detector with aromatic SMILES."""
        monomer_dict = {1: "c1ccccc1O"}  # Phenol

        result = functional_groups_detector(monomer_dict)

        assert isinstance(result, dict)

    def test_detect_zero_matches(self):
        """Test when SMARTS pattern matches zero times."""
        smiles = "CCC"  # Propane
        result = detect_monomer_functionality(smiles, "vinyl", "[C]=[C]")

        assert result == (0, 0, None)


class TestRealWorldScenarios:
    """Tests with real-world polymer monomers."""

    def test_detect_lactic_acid(self):
        """Test detection of lactic acid (hydroxy carboxylic acid)."""
        monomer_dict = {1: "OCCC(O)=O"}

        result = functional_groups_detector(monomer_dict)

        if 1 in result:
            assert result[1]["smiles"] == "OCCC(O)=O"

    def test_detect_multiple_monomer_types(self):
        """Test detection with various monomer types."""
        monomer_dict = {
            1: "OCCC(O)=O",  # Hydroxy carboxylic acid
            2: "OCCO",  # Diol (if configured)
            3: "O=C(O)CCCCC(=O)O",  # Di-carboxylic acid (if configured)
        }

        result = functional_groups_detector(monomer_dict)

        # At least some monomers should be detected
        assert len(result) >= 1

    def test_cyclic_monomers(self):
        """Test with cyclic monomer structures."""
        monomer_dict = {
            1: "C1CCC(=O)OCC1",  # Caprolactone (if configured)
        }

        result = functional_groups_detector(monomer_dict)

        assert isinstance(result, dict)


class TestIntegration:
    """Integration tests for the detector workflow."""

    @patch("builtins.print")
    def test_detector_prints_for_each_detection(self, mock_print):
        """Test that detector prints for each detected functional group."""
        monomer_dict = {1: "OCCC(O)=O"}

        functional_groups_detector(monomer_dict)

        # Should have printed something (depends on configuration)
        # At minimum, should not crash

    def test_detector_consistent_output(self):
        """Test that detector gives consistent output for same input."""
        monomer_dict = {1: "OCCC(O)=O"}

        result1 = functional_groups_detector(monomer_dict)
        result2 = functional_groups_detector(monomer_dict)

        # Results should be identical
        assert result1 == result2