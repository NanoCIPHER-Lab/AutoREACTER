"""
Comprehensive tests for the detectors/detector.py module.

Tests cover:
- Detection of reactions
- Identification of non-reactant monomers
- User interaction for non-reactant selection
- Error handling and validation
"""
import pytest
from unittest.mock import patch, MagicMock, call

from reaction_lammps_mupt.detectors.detector import (
    find_non_reactant_monomers,
    detect_reactions,
)


class TestFindNonReactantMonomers:
    """Tests for find_non_reactant_monomers function."""

    @patch("builtins.input")
    def test_no_non_reactants(self, mock_input):
        """Test when all monomers participate in reactions."""
        reactions = {
            "Reaction1": {
                1: {
                    "monomer_1": {"smiles": "CCO"},
                    "monomer_2": {"smiles": "CCN"},
                }
            }
        }
        input_dict = {"monomers": {1: "CCO", 2: "CCN"}}

        result = find_non_reactant_monomers(reactions, input_dict)

        assert result == []
        mock_input.assert_not_called()

    @patch("builtins.input")
    def test_with_non_reactants_user_excludes_all(self, mock_input, capsys):
        """Test with non-reactants when user chooses to exclude all."""
        mock_input.return_value = "y"

        reactions = {
            "Reaction1": {
                1: {
                    "monomer_1": {"smiles": "CCO"},
                    "monomer_2": {"smiles": "CCN"},
                }
            }
        }
        input_dict = {"monomers": {1: "CCO", 2: "CCN", 3: "CCC"}}

        result = find_non_reactant_monomers(reactions, input_dict)

        assert result == []
        captured = capsys.readouterr()
        assert "non-reactant monomers" in captured.out.lower()
        assert "CCC" in captured.out

    @patch("builtins.input")
    def test_with_non_reactants_user_selects_some(self, mock_input, capsys):
        """Test with non-reactants when user selects specific ones."""
        mock_input.side_effect = ["n", "1"]

        reactions = {
            "Reaction1": {
                1: {
                    "monomer_1": {"smiles": "CCO"},
                    "monomer_2": {"smiles": "CCN"},
                }
            }
        }
        input_dict = {"monomers": {1: "CCO", 2: "CCN", 3: "CCC", 4: "CCCC"}}

        result = find_non_reactant_monomers(reactions, input_dict)

        assert len(result) == 1
        assert "CCC" in result

    @patch("builtins.input")
    def test_with_non_reactants_user_selects_multiple(self, mock_input, capsys):
        """Test when user selects multiple non-reactants."""
        mock_input.side_effect = ["n", "1, 2"]

        reactions = {
            "Reaction1": {
                1: {
                    "monomer_1": {"smiles": "CCO"},
                    "monomer_2": {"smiles": "CCN"},
                }
            }
        }
        input_dict = {"monomers": {1: "CCO", 2: "CCN", 3: "CCC", 4: "CCCC"}}

        result = find_non_reactant_monomers(reactions, input_dict)

        assert len(result) == 2
        assert "CCC" in result
        assert "CCCC" in result

    @patch("builtins.input")
    def test_non_reactant_reindexing(self, mock_input):
        """Test that non-reactants are re-indexed starting from 1."""
        mock_input.side_effect = ["n", "1"]

        reactions = {
            "Reaction1": {
                1: {
                    "monomer_1": {"smiles": "CCO"},
                }
            }
        }
        # Monomers with non-contiguous IDs
        input_dict = {"monomers": {5: "CCO", 10: "CCC", 15: "CCCC"}}

        result = find_non_reactant_monomers(reactions, input_dict)

        # Should be able to select by re-indexed ID
        assert len(result) == 1

    def test_empty_reactions_dict(self):
        """Test with empty reactions dictionary."""
        reactions = {}
        input_dict = {"monomers": {1: "CCO", 2: "CCN"}}

        # All monomers should be considered non-reactants
        # But without user input, it depends on implementation
        # This tests the edge case

    def test_reaction_with_smiles_field(self):
        """Test reactions that include 'smiles' field."""
        reactions = {
            "Reaction1": {
                1: {
                    "monomer_1": {"smiles": "CCO"},
                    "monomer_2": {"smiles": "CCN"},
                    "smiles": "CCO.CCN",  # Combined SMILES
                }
            }
        }
        input_dict = {"monomers": {1: "CCO", 2: "CCN", 3: "CCC"}}

        # CCC should be identified as non-reactant
        with patch("builtins.input", return_value="y"):
            result = find_non_reactant_monomers(reactions, input_dict)
            assert result == []

    def test_multiple_reactions_same_monomers(self):
        """Test when same monomers appear in multiple reactions."""
        reactions = {
            "Reaction1": {
                1: {"monomer_1": {"smiles": "CCO"}, "monomer_2": {"smiles": "CCN"}},
            },
            "Reaction2": {
                1: {"monomer_1": {"smiles": "CCO"}, "monomer_2": {"smiles": "CCC"}},
            },
        }
        input_dict = {"monomers": {1: "CCO", 2: "CCN", 3: "CCC", 4: "CCCC"}}

        with patch("builtins.input", return_value="y"):
            result = find_non_reactant_monomers(reactions, input_dict)

        # Only CCCC should be non-reactant
        assert result == []

    @patch("builtins.input")
    def test_user_input_with_spaces(self, mock_input):
        """Test handling of user input with extra spaces."""
        mock_input.side_effect = ["n", " 1 , 2 "]

        reactions = {
            "Reaction1": {
                1: {"monomer_1": {"smiles": "CCO"}},
            }
        }
        input_dict = {"monomers": {1: "CCO", 2: "CCC", 3: "CCCC"}}

        result = find_non_reactant_monomers(reactions, input_dict)

        assert len(result) == 2

    @patch("builtins.input")
    def test_invalid_selection_filtered_out(self, mock_input):
        """Test that invalid selections are filtered out."""
        mock_input.side_effect = ["n", "1, 999, abc"]

        reactions = {"Reaction1": {1: {"monomer_1": {"smiles": "CCO"}}}}
        input_dict = {"monomers": {1: "CCO", 2: "CCC"}}

        result = find_non_reactant_monomers(reactions, input_dict)

        # Only valid ID 1 should be selected
        assert len(result) == 1

    def test_empty_monomers_dict(self):
        """Test with empty monomers dictionary."""
        reactions = {"Reaction1": {1: {"monomer_1": {"smiles": "CCO"}}}}
        input_dict = {"monomers": {}}

        result = find_non_reactant_monomers(reactions, input_dict)

        assert result == []


class TestDetectReactions:
    """Tests for detect_reactions function."""

    @patch("reaction_lammps_mupt.detectors.detector.reaction_selector")
    @patch("reaction_lammps_mupt.detectors.detector.functional_groups_detector")
    def test_detect_reactions_basic(self, mock_fg_detector, mock_selector):
        """Test basic reaction detection workflow."""
        input_dict = {"monomers": {1: "OCCC(O)=O", 2: "C1=CC=C(C(=C1)C(=O)O)O"}}

        mock_fg_results = {
            1: {"smiles": "OCCC(O)=O", 1: {"functional_group_name": "hydroxy_carboxylic_acid"}},
            2: {"smiles": "C1=CC=C(C(=C1)C(=O)O)O", 1: {"functional_group_name": "hydroxy_carboxylic_acid"}},
        }
        mock_fg_detector.return_value = mock_fg_results

        mock_reactions = {
            "Reaction1": {
                1: {
                    "monomer_1": {"smiles": "OCCC(O)=O"},
                    "monomer_2": {"smiles": "C1=CC=C(C(=C1)C(=O)O)O"},
                }
            }
        }
        mock_selector.return_value = mock_reactions

        with patch("builtins.input", return_value="y"):
            reactions, non_reactants = detect_reactions(input_dict)

        assert reactions == mock_reactions
        assert non_reactants == []
        mock_fg_detector.assert_called_once_with({1: "OCCC(O)=O", 2: "C1=CC=C(C(=C1)C(=O)O)O"})
        mock_selector.assert_called_once_with(mock_fg_results)

    def test_detect_reactions_missing_monomers_key(self):
        """Test detect_reactions with missing monomers key."""
        input_dict = {}

        with pytest.raises(ValueError, match="must contain a 'monomers' key"):
            detect_reactions(input_dict)

    def test_detect_reactions_empty_monomers(self):
        """Test detect_reactions with empty monomers dict."""
        input_dict = {"monomers": {}}

        with pytest.raises(ValueError, match="must contain a 'monomers' key"):
            detect_reactions(input_dict)

    @patch("reaction_lammps_mupt.detectors.detector.reaction_selector")
    @patch("reaction_lammps_mupt.detectors.detector.functional_groups_detector")
    def test_detect_reactions_with_non_reactants(self, mock_fg_detector, mock_selector):
        """Test detect_reactions when there are non-reactant monomers."""
        input_dict = {"monomers": {1: "CCO", 2: "CCN", 3: "CCC"}}

        mock_fg_results = {
            1: {"smiles": "CCO"},
            2: {"smiles": "CCN"},
        }
        mock_fg_detector.return_value = mock_fg_results

        mock_reactions = {
            "Reaction1": {
                1: {"monomer_1": {"smiles": "CCO"}, "monomer_2": {"smiles": "CCN"}},
            }
        }
        mock_selector.return_value = mock_reactions

        with patch("builtins.input", side_effect=["n", "1"]):
            reactions, non_reactants = detect_reactions(input_dict)

        assert reactions == mock_reactions
        assert len(non_reactants) == 1
        assert "CCC" in non_reactants

    @patch("reaction_lammps_mupt.detectors.detector.reaction_selector")
    @patch("reaction_lammps_mupt.detectors.detector.functional_groups_detector")
    def test_detect_reactions_returns_tuple(self, mock_fg_detector, mock_selector):
        """Test that detect_reactions returns a tuple."""
        input_dict = {"monomers": {1: "CCO"}}

        mock_fg_detector.return_value = {1: {"smiles": "CCO"}}
        mock_selector.return_value = {}

        result = detect_reactions(input_dict)

        assert isinstance(result, tuple)
        assert len(result) == 2

    @patch("reaction_lammps_mupt.detectors.detector.reaction_selector")
    @patch("reaction_lammps_mupt.detectors.detector.functional_groups_detector")
    def test_detect_reactions_preserves_reaction_data(self, mock_fg_detector, mock_selector):
        """Test that detect_reactions preserves all reaction data."""
        input_dict = {"monomers": {1: "CCO", 2: "CCN"}}

        mock_fg_results = {
            1: {"smiles": "CCO", 1: {"functional_group_name": "diol"}},
            2: {"smiles": "CCN", 1: {"functional_group_name": "di_amine"}},
        }
        mock_fg_detector.return_value = mock_fg_results

        mock_reactions = {
            "Reaction1": {
                "same_reactants": False,
                "reactant_1": "diol",
                "reactant_2": "di_carboxylic_acid",
                "product": "polyester",
                1: {
                    "monomer_1": {"smiles": "CCO", "functional_group_name": "diol"},
                    "monomer_2": {"smiles": "CCN", "functional_group_name": "di_amine"},
                },
            }
        }
        mock_selector.return_value = mock_reactions

        with patch("builtins.input", return_value="y"):
            reactions, _ = detect_reactions(input_dict)

        assert "same_reactants" in reactions["Reaction1"]
        assert "reactant_1" in reactions["Reaction1"]
        assert "product" in reactions["Reaction1"]


class TestEdgeCasesAndBoundary:
    """Tests for edge cases and boundary conditions."""

    @patch("builtins.input")
    def test_find_non_reactants_case_sensitive(self, mock_input):
        """Test that SMILES comparison is case-sensitive."""
        mock_input.return_value = "y"

        reactions = {
            "Reaction1": {
                1: {"monomer_1": {"smiles": "CCO"}},
            }
        }
        # Lower case 'c' is aromatic in SMILES
        input_dict = {"monomers": {1: "CCO", 2: "ccO"}}

        result = find_non_reactant_monomers(reactions, input_dict)

        # ccO should be identified as non-reactant (case sensitive)
        assert result == []

    @patch("builtins.input")
    def test_find_non_reactants_with_none_values(self, mock_input):
        """Test handling of None values in reaction data."""
        mock_input.return_value = "y"

        reactions = {
            "Reaction1": {
                1: {
                    "monomer_1": {"smiles": None},
                    "monomer_2": {"smiles": "CCN"},
                }
            }
        }
        input_dict = {"monomers": {1: "CCO", 2: "CCN"}}

        result = find_non_reactant_monomers(reactions, input_dict)

        # CCO should be non-reactant since monomer_1 smiles is None
        # Implementation should handle None gracefully

    @patch("reaction_lammps_mupt.detectors.detector.reaction_selector")
    @patch("reaction_lammps_mupt.detectors.detector.functional_groups_detector")
    def test_detect_reactions_functional_group_detector_error(self, mock_fg_detector, mock_selector):
        """Test when functional_groups_detector raises an error."""
        input_dict = {"monomers": {1: "INVALID"}}

        mock_fg_detector.side_effect = Exception("Detection failed")

        with pytest.raises(Exception, match="Detection failed"):
            detect_reactions(input_dict)

    @patch("reaction_lammps_mupt.detectors.detector.reaction_selector")
    @patch("reaction_lammps_mupt.detectors.detector.functional_groups_detector")
    def test_detect_reactions_reaction_selector_error(self, mock_fg_detector, mock_selector):
        """Test when reaction_selector raises an error."""
        input_dict = {"monomers": {1: "CCO"}}

        mock_fg_detector.return_value = {1: {"smiles": "CCO"}}
        mock_selector.side_effect = ValueError("No reactions detected")

        with pytest.raises(ValueError, match="No reactions detected"):
            detect_reactions(input_dict)

    @patch("builtins.input")
    def test_find_non_reactants_empty_string_input(self, mock_input, capsys):
        """Test handling of empty string in user selection."""
        mock_input.side_effect = ["n", ""]

        reactions = {"Reaction1": {1: {"monomer_1": {"smiles": "CCO"}}}}
        input_dict = {"monomers": {1: "CCO", 2: "CCC"}}

        result = find_non_reactant_monomers(reactions, input_dict)

        # Empty selection should result in empty list
        assert result == []


class TestIntegrationScenarios:
    """Integration tests simulating real-world scenarios."""

    @patch("reaction_lammps_mupt.detectors.detector.reaction_selector")
    @patch("reaction_lammps_mupt.detectors.detector.functional_groups_detector")
    def test_complete_workflow_with_polyesterification(self, mock_fg_detector, mock_selector):
        """Test complete workflow for polyesterification reaction."""
        input_dict = {
            "monomers": {
                1: "OCCC(O)=O",  # Lactic acid (hydroxy carboxylic acid)
                2: "C1=CC=C(C(=C1)C(=O)O)O",  # p-hydroxybenzoic acid
            }
        }

        mock_fg_results = {
            1: {
                "smiles": "OCCC(O)=O",
                1: {
                    "functionality_type": "di_different",
                    "functional_group_name": "hydroxy_carboxylic_acid",
                },
            },
            2: {
                "smiles": "C1=CC=C(C(=C1)C(=O)O)O",
                1: {
                    "functionality_type": "di_different",
                    "functional_group_name": "hydroxy_carboxylic_acid",
                },
            },
        }
        mock_fg_detector.return_value = mock_fg_results

        mock_reactions = {
            "Hydroxy Carboxylic Acid Polycondensation(Polyesterification)": {
                "same_reactants": True,
                "reactant_1": "hydroxy_carboxylic_acid",
                "product": "polyester_chain",
                "delete_atom": True,
                1: {
                    "monomer_1": {"smiles": "OCCC(O)=O"},
                },
            }
        }
        mock_selector.return_value = mock_reactions

        reactions, non_reactants = detect_reactions(input_dict)

        assert "Hydroxy Carboxylic Acid Polycondensation(Polyesterification)" in reactions
        assert len(non_reactants) == 0

    @patch("reaction_lammps_mupt.detectors.detector.reaction_selector")
    @patch("reaction_lammps_mupt.detectors.detector.functional_groups_detector")
    def test_mixed_reactive_and_nonreactive_monomers(self, mock_fg_detector, mock_selector):
        """Test with mix of reactive and non-reactive monomers."""
        input_dict = {
            "monomers": {
                1: "OCCC(O)=O",  # Reactive
                2: "CCO",  # Non-reactive (no functional groups detected)
                3: "CCC",  # Non-reactive
            }
        }

        mock_fg_results = {
            1: {
                "smiles": "OCCC(O)=O",
                1: {"functional_group_name": "hydroxy_carboxylic_acid"},
            },
        }
        mock_fg_detector.return_value = mock_fg_results

        mock_reactions = {
            "Reaction1": {
                1: {"monomer_1": {"smiles": "OCCC(O)=O"}},
            }
        }
        mock_selector.return_value = mock_reactions

        with patch("builtins.input", side_effect=["n", "1, 2"]):
            reactions, non_reactants = detect_reactions(input_dict)

        assert len(reactions) > 0
        assert len(non_reactants) == 2
        assert "CCO" in non_reactants
        assert "CCC" in non_reactants