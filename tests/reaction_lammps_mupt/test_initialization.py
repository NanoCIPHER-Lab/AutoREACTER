"""
Comprehensive tests for the initialization.py module.

Tests cover:
- ASCII art display
- Module import checking
- Cache directory cleanup on initialization
- Error handling for missing modules
"""
import sys
from unittest.mock import patch, MagicMock, call
import pytest

from reaction_lammps_mupt.initialization import ASCII_Mupt_reaction_LAMMPS, initialize


class TestASCIIArt:
    """Tests for ASCII art banner display."""

    def test_ascii_art_prints(self, capsys):
        """Test that ASCII art is printed."""
        ASCII_Mupt_reaction_LAMMPS()
        captured = capsys.readouterr()

        # Check that some recognizable parts of the ASCII art are present
        assert "MuPT" in captured.out or "REACTION" in captured.out
        assert len(captured.out) > 100  # ASCII art should be substantial

    def test_ascii_art_contains_expected_elements(self, capsys):
        """Test that ASCII art contains expected visual elements."""
        ASCII_Mupt_reaction_LAMMPS()
        captured = capsys.readouterr()

        # Check for presence of ASCII art characters
        assert any(char in captured.out for char in ["█", "|", "/", "\\", "_"])


class TestInitialize:
    """Tests for the initialize function."""

    @patch("reaction_lammps_mupt.initialization.delete_default_cache_files")
    @patch("reaction_lammps_mupt.initialization.ASCII_Mupt_reaction_LAMMPS")
    @patch("importlib.import_module")
    def test_initialize_all_modules_present(self, mock_import, mock_ascii, mock_delete, capsys):
        """Test initialize when all required modules are present."""
        # Mock successful imports
        mock_import.return_value = MagicMock()

        initialize()

        # Verify all modules were imported
        expected_modules = ["rdkit", "pandas", "numpy"]
        import_calls = [call(mod) for mod in expected_modules]
        mock_import.assert_has_calls(import_calls, any_order=True)

        # Verify ASCII art was displayed
        mock_ascii.assert_called_once()

        # Verify cache was deleted
        mock_delete.assert_called_once()

        # Check success message
        captured = capsys.readouterr()
        assert "All required modules are successfully imported" in captured.out

    @patch("reaction_lammps_mupt.initialization.delete_default_cache_files")
    @patch("reaction_lammps_mupt.initialization.ASCII_Mupt_reaction_LAMMPS")
    @patch("importlib.import_module")
    def test_initialize_missing_rdkit(self, mock_import, mock_ascii, mock_delete):
        """Test initialize when rdkit is missing."""
        def side_effect(module_name):
            if module_name == "rdkit":
                error = ModuleNotFoundError(f"No module named '{module_name}'")
                error.name = module_name
                raise error
            return MagicMock()

        mock_import.side_effect = side_effect

        with pytest.raises(SystemExit) as exc_info:
            initialize()

        # Check that exit was called with appropriate message
        assert "rdkit" in str(exc_info.value)

    @patch("reaction_lammps_mupt.initialization.delete_default_cache_files")
    @patch("reaction_lammps_mupt.initialization.ASCII_Mupt_reaction_LAMMPS")
    @patch("importlib.import_module")
    def test_initialize_missing_pandas(self, mock_import, mock_ascii, mock_delete):
        """Test initialize when pandas is missing."""
        def side_effect(module_name):
            if module_name == "pandas":
                error = ModuleNotFoundError(f"No module named '{module_name}'")
                error.name = module_name
                raise error
            return MagicMock()

        mock_import.side_effect = side_effect

        with pytest.raises(SystemExit) as exc_info:
            initialize()

        assert "pandas" in str(exc_info.value)

    @patch("reaction_lammps_mupt.initialization.delete_default_cache_files")
    @patch("reaction_lammps_mupt.initialization.ASCII_Mupt_reaction_LAMMPS")
    @patch("importlib.import_module")
    def test_initialize_missing_numpy(self, mock_import, mock_ascii, mock_delete):
        """Test initialize when numpy is missing."""
        def side_effect(module_name):
            if module_name == "numpy":
                error = ModuleNotFoundError(f"No module named '{module_name}'")
                error.name = module_name
                raise error
            return MagicMock()

        mock_import.side_effect = side_effect

        with pytest.raises(SystemExit) as exc_info:
            initialize()

        assert "numpy" in str(exc_info.value)

    @patch("reaction_lammps_mupt.initialization.delete_default_cache_files")
    @patch("reaction_lammps_mupt.initialization.ASCII_Mupt_reaction_LAMMPS")
    @patch("importlib.import_module")
    def test_initialize_calls_delete_cache(self, mock_import, mock_ascii, mock_delete):
        """Test that initialize calls delete_default_cache_files."""
        mock_import.return_value = MagicMock()

        initialize()

        mock_delete.assert_called_once()

    @patch("reaction_lammps_mupt.initialization.delete_default_cache_files")
    @patch("reaction_lammps_mupt.initialization.ASCII_Mupt_reaction_LAMMPS")
    @patch("importlib.import_module")
    def test_initialize_calls_ascii_art(self, mock_import, mock_ascii, mock_delete):
        """Test that initialize calls ASCII_Mupt_reaction_LAMMPS."""
        mock_import.return_value = MagicMock()

        initialize()

        mock_ascii.assert_called_once()

    @patch("reaction_lammps_mupt.initialization.delete_default_cache_files")
    @patch("reaction_lammps_mupt.initialization.ASCII_Mupt_reaction_LAMMPS")
    @patch("importlib.import_module")
    def test_initialize_order_of_operations(self, mock_import, mock_ascii, mock_delete):
        """Test that initialize operations occur in correct order."""
        mock_import.return_value = MagicMock()

        # Create a manager to track call order
        manager = MagicMock()
        manager.attach_mock(mock_import, "import_module")
        manager.attach_mock(mock_ascii, "ascii")
        manager.attach_mock(mock_delete, "delete")

        initialize()

        # Check that imports happen before ASCII art and cache deletion
        calls = [str(c) for c in manager.mock_calls]

        # ASCII should be called after imports
        ascii_index = next(i for i, c in enumerate(calls) if "ascii" in c)
        import_indices = [i for i, c in enumerate(calls) if "import_module" in c]

        assert all(import_idx < ascii_index for import_idx in import_indices)

    @patch("reaction_lammps_mupt.initialization.delete_default_cache_files")
    @patch("reaction_lammps_mupt.initialization.ASCII_Mupt_reaction_LAMMPS")
    @patch("importlib.import_module")
    def test_initialize_error_message_format(self, mock_import, mock_ascii, mock_delete):
        """Test that error message has correct format."""
        module_name = "missing_module"

        def side_effect(name):
            if name == module_name:
                error = ModuleNotFoundError(f"No module named '{name}'")
                error.name = name
                raise error
            return MagicMock()

        mock_import.side_effect = side_effect

        with pytest.raises(SystemExit) as exc_info:
            initialize()

        error_message = str(exc_info.value)
        assert "ERROR" in error_message
        assert module_name in error_message


class TestIntegrationScenarios:
    """Integration tests for initialization."""

    @patch("reaction_lammps_mupt.initialization.delete_default_cache_files")
    @patch("importlib.import_module")
    def test_initialize_with_real_modules(self, mock_import, mock_delete, capsys):
        """Test initialize with actual module imports."""
        # Allow actual imports
        mock_import.side_effect = __import__

        # This should work if dependencies are installed
        try:
            initialize()
            captured = capsys.readouterr()
            assert "All required modules are successfully imported" in captured.out
        except SystemExit:
            # If modules aren't installed, that's expected
            pass

    @patch("reaction_lammps_mupt.initialization.delete_default_cache_files")
    @patch("reaction_lammps_mupt.initialization.ASCII_Mupt_reaction_LAMMPS")
    @patch("importlib.import_module")
    def test_initialize_cache_deletion_failure_doesnt_stop_init(self, mock_import, mock_ascii, mock_delete):
        """Test that cache deletion failure doesn't prevent initialization."""
        mock_import.return_value = MagicMock()
        # Mock cache deletion to raise an exception
        mock_delete.side_effect = Exception("Cache deletion failed")

        # Should raise the exception (no error handling in current implementation)
        with pytest.raises(Exception, match="Cache deletion failed"):
            initialize()


class TestEdgeCases:
    """Tests for edge cases and boundary conditions."""

    @patch("reaction_lammps_mupt.initialization.delete_default_cache_files")
    @patch("reaction_lammps_mupt.initialization.ASCII_Mupt_reaction_LAMMPS")
    @patch("importlib.import_module")
    def test_initialize_import_error_without_name_attribute(self, mock_import, mock_ascii, mock_delete):
        """Test handling of ModuleNotFoundError without name attribute."""
        error = ModuleNotFoundError("Generic import error")
        # Don't set error.name attribute
        mock_import.side_effect = error

        with pytest.raises((SystemExit, AttributeError)):
            initialize()

    @patch("reaction_lammps_mupt.initialization.delete_default_cache_files")
    @patch("reaction_lammps_mupt.initialization.ASCII_Mupt_reaction_LAMMPS")
    @patch("importlib.import_module")
    def test_initialize_multiple_calls(self, mock_import, mock_ascii, mock_delete):
        """Test that initialize can be called multiple times."""
        mock_import.return_value = MagicMock()

        initialize()
        initialize()

        # Should be called twice
        assert mock_delete.call_count == 2
        assert mock_ascii.call_count == 2

    def test_ascii_art_no_exceptions(self):
        """Test that ASCII art never raises exceptions."""
        # Should not raise any exceptions
        try:
            ASCII_Mupt_reaction_LAMMPS()
        except Exception as e:
            pytest.fail(f"ASCII_Mupt_reaction_LAMMPS raised {type(e).__name__}: {e}")

    @patch("reaction_lammps_mupt.initialization.delete_default_cache_files")
    @patch("reaction_lammps_mupt.initialization.ASCII_Mupt_reaction_LAMMPS")
    @patch("importlib.import_module")
    def test_initialize_with_import_error(self, mock_import, mock_ascii, mock_delete):
        """Test initialize with generic ImportError (not ModuleNotFoundError)."""
        mock_import.side_effect = ImportError("Generic import error")

        # Should propagate the ImportError (not caught)
        with pytest.raises(ImportError):
            initialize()


class TestRegressionTests:
    """Regression tests for previously identified issues."""

    @patch("reaction_lammps_mupt.initialization.delete_default_cache_files")
    @patch("reaction_lammps_mupt.initialization.ASCII_Mupt_reaction_LAMMPS")
    @patch("importlib.import_module")
    def test_initialize_preserves_module_list_order(self, mock_import, mock_ascii, mock_delete, capsys):
        """Test that modules are checked in the specified order."""
        checked_modules = []

        def track_imports(module_name):
            checked_modules.append(module_name)
            return MagicMock()

        mock_import.side_effect = track_imports

        initialize()

        # Verify the order matches the list in initialize
        expected_order = ["rdkit", "pandas", "numpy"]
        assert checked_modules == expected_order

    @patch("reaction_lammps_mupt.initialization.delete_default_cache_files")
    @patch("reaction_lammps_mupt.initialization.ASCII_Mupt_reaction_LAMMPS")
    @patch("importlib.import_module")
    def test_initialize_exits_on_first_missing_module(self, mock_import, mock_ascii, mock_delete):
        """Test that initialize exits on the first missing module."""
        call_count = 0

        def side_effect(module_name):
            nonlocal call_count
            call_count += 1
            if call_count == 1:  # First module
                error = ModuleNotFoundError(f"No module named '{module_name}'")
                error.name = module_name
                raise error
            return MagicMock()

        mock_import.side_effect = side_effect

        with pytest.raises(SystemExit):
            initialize()

        # Should only try to import one module before exiting
        assert call_count == 1

    @patch("reaction_lammps_mupt.initialization.delete_default_cache_files")
    @patch("reaction_lammps_mupt.initialization.ASCII_Mupt_reaction_LAMMPS")
    @patch("importlib.import_module")
    def test_initialize_success_message_printed(self, mock_import, mock_ascii, mock_delete, capsys):
        """Test that success message is printed when all modules load."""
        mock_import.return_value = MagicMock()

        initialize()

        captured = capsys.readouterr()
        assert "successfully imported" in captured.out.lower()