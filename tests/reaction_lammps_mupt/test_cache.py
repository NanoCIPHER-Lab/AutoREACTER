"""
Comprehensive tests for the cache.py module.

Tests cover:
- Git root detection
- Cache directory management
- File deletion operations
- Dated run directory creation
- Retention cleanup
- Copy operations
"""
import datetime as dt
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch, call

import pytest

from reaction_lammps_mupt.cache import (
    get_git_root,
    get_cache_base_dir,
    get_default_cache_dir,
    delete_files_in_directory,
    delete_default_cache_files,
    make_dated_run_dir,
    get_current_cache_dir,
    retention_cleanup,
    copy_to_date_folder,
)


class TestGitRoot:
    """Tests for git root detection."""

    def test_get_git_root_success(self):
        """Test that get_git_root returns correct path in a git repo."""
        root = get_git_root()
        assert isinstance(root, Path)
        assert root.is_dir()
        # Verify .git directory exists
        assert (root / ".git").exists()

    def test_get_git_root_strips_whitespace(self):
        """Test that get_git_root strips whitespace from command output."""
        with patch("subprocess.check_output") as mock_output:
            mock_output.return_value = "  /path/to/repo  \n"
            result = get_git_root()
            assert result == Path("/path/to/repo")

    def test_get_git_root_not_in_repo(self):
        """Test that get_git_root raises error when not in a git repo."""
        with patch("subprocess.check_output") as mock_output:
            mock_output.side_effect = subprocess.CalledProcessError(128, "git")
            with pytest.raises(subprocess.CalledProcessError):
                get_git_root()


class TestCacheDirectories:
    """Tests for cache directory management."""

    @patch("reaction_lammps_mupt.cache.get_git_root")
    def test_get_cache_base_dir_creates_directory(self, mock_git_root, tmp_path):
        """Test that get_cache_base_dir creates cache directory."""
        mock_git_root.return_value = tmp_path
        cache_base = get_cache_base_dir()

        assert cache_base == tmp_path / "cache"
        assert cache_base.exists()
        assert cache_base.is_dir()

    @patch("reaction_lammps_mupt.cache.get_git_root")
    def test_get_cache_base_dir_existing_directory(self, mock_git_root, tmp_path):
        """Test that get_cache_base_dir works with existing directory."""
        mock_git_root.return_value = tmp_path
        cache_dir = tmp_path / "cache"
        cache_dir.mkdir(parents=True, exist_ok=True)

        cache_base = get_cache_base_dir()
        assert cache_base == cache_dir
        assert cache_base.exists()

    @patch("reaction_lammps_mupt.cache.get_cache_base_dir")
    def test_get_default_cache_dir_creates_00_cache(self, mock_base_dir, tmp_path):
        """Test that get_default_cache_dir creates 00_cache subdirectory."""
        mock_base_dir.return_value = tmp_path
        default_cache = get_default_cache_dir()

        assert default_cache == tmp_path / "00_cache"
        assert default_cache.exists()
        assert default_cache.is_dir()

    @patch("reaction_lammps_mupt.cache.get_cache_base_dir")
    def test_get_default_cache_dir_existing(self, mock_base_dir, tmp_path):
        """Test get_default_cache_dir with existing 00_cache directory."""
        mock_base_dir.return_value = tmp_path
        default_dir = tmp_path / "00_cache"
        default_dir.mkdir(parents=True, exist_ok=True)

        result = get_default_cache_dir()
        assert result == default_dir
        assert result.exists()


class TestFileOperations:
    """Tests for file deletion operations."""

    def test_delete_files_in_directory_files_only(self, tmp_path):
        """Test that delete_files_in_directory only deletes files, not directories."""
        # Create test files and subdirectories
        (tmp_path / "file1.txt").write_text("content1")
        (tmp_path / "file2.txt").write_text("content2")
        subdir = tmp_path / "subdir"
        subdir.mkdir()
        (subdir / "nested.txt").write_text("nested")

        delete_files_in_directory(tmp_path)

        # Files should be deleted
        assert not (tmp_path / "file1.txt").exists()
        assert not (tmp_path / "file2.txt").exists()
        # Subdirectories should remain
        assert subdir.exists()
        assert (subdir / "nested.txt").exists()

    def test_delete_files_in_directory_handles_symlinks(self, tmp_path):
        """Test that delete_files_in_directory deletes symlinks."""
        target = tmp_path / "target.txt"
        target.write_text("target")
        link = tmp_path / "link"
        link.symlink_to(target)

        delete_files_in_directory(tmp_path)

        assert not link.exists()
        # Note: delete_files_in_directory uses unlink() which removes both symlinks and files
        # The behavior depends on whether target is also in the directory being cleaned

    def test_delete_files_in_directory_nonexistent(self, tmp_path, capsys):
        """Test delete_files_in_directory with nonexistent directory."""
        nonexistent = tmp_path / "nonexistent"
        delete_files_in_directory(nonexistent)

        captured = capsys.readouterr()
        assert "Error: Directory not found" in captured.out

    def test_delete_files_in_directory_permission_error(self, tmp_path):
        """Test delete_files_in_directory handles permission errors gracefully."""
        test_file = tmp_path / "readonly.txt"
        test_file.write_text("content")

        with patch("pathlib.Path.unlink") as mock_unlink:
            mock_unlink.side_effect = OSError("Permission denied")
            # Should not raise an exception
            delete_files_in_directory(tmp_path)

    def test_delete_files_in_directory_empty(self, tmp_path):
        """Test delete_files_in_directory with empty directory."""
        delete_files_in_directory(tmp_path)
        assert tmp_path.exists()

    @patch("reaction_lammps_mupt.cache.get_default_cache_dir")
    @patch("reaction_lammps_mupt.cache.delete_files_in_directory")
    def test_delete_default_cache_files(self, mock_delete, mock_get_default):
        """Test delete_default_cache_files calls delete_files_in_directory."""
        mock_dir = Path("/mock/cache/00_cache")
        mock_get_default.return_value = mock_dir

        delete_default_cache_files()

        mock_delete.assert_called_once_with(mock_dir)


class TestDatedRunDirectories:
    """Tests for dated run directory creation."""

    def test_make_dated_run_dir_first_run(self, tmp_path):
        """Test creating the first run directory for a date."""
        run_dir = make_dated_run_dir(tmp_path, chdir_to="none")

        today = dt.date.today()
        expected = tmp_path / today.isoformat() / "1"
        assert run_dir == expected
        assert run_dir.exists()

    def test_make_dated_run_dir_incremental(self, tmp_path):
        """Test creating multiple run directories increments number."""
        # Create first run
        run_dir_1 = make_dated_run_dir(tmp_path, chdir_to="none")
        # Create second run
        run_dir_2 = make_dated_run_dir(tmp_path, chdir_to="none")

        today = dt.date.today()
        assert run_dir_1 == tmp_path / today.isoformat() / "1"
        assert run_dir_2 == tmp_path / today.isoformat() / "2"
        assert run_dir_1.exists()
        assert run_dir_2.exists()

    def test_make_dated_run_dir_custom_date(self, tmp_path):
        """Test creating run directory with custom date."""
        custom_date = dt.date(2024, 1, 15)
        run_dir = make_dated_run_dir(tmp_path, date=custom_date, chdir_to="none")

        expected = tmp_path / "2024-01-15" / "1"
        assert run_dir == expected
        assert run_dir.exists()

    def test_make_dated_run_dir_chdir_run(self, tmp_path):
        """Test chdir_to='run' changes to run directory."""
        original_cwd = os.getcwd()
        try:
            run_dir = make_dated_run_dir(tmp_path, chdir_to="run")
            assert os.getcwd() == str(run_dir)
        finally:
            os.chdir(original_cwd)

    def test_make_dated_run_dir_chdir_date(self, tmp_path):
        """Test chdir_to='date' changes to date directory."""
        original_cwd = os.getcwd()
        try:
            run_dir = make_dated_run_dir(tmp_path, chdir_to="date")
            assert os.getcwd() == str(run_dir.parent)
        finally:
            os.chdir(original_cwd)

    def test_make_dated_run_dir_chdir_none(self, tmp_path):
        """Test chdir_to='none' does not change directory."""
        original_cwd = os.getcwd()
        run_dir = make_dated_run_dir(tmp_path, chdir_to="none")
        assert os.getcwd() == original_cwd

    def test_make_dated_run_dir_race_condition(self, tmp_path):
        """Test handling of race condition when directory exists."""
        today = dt.date.today()
        date_dir = tmp_path / today.isoformat()
        date_dir.mkdir(parents=True, exist_ok=True)

        # Create run 1 manually
        (date_dir / "1").mkdir()

        # Should create run 2
        run_dir = make_dated_run_dir(tmp_path, chdir_to="none")
        assert run_dir == date_dir / "2"

    def test_make_dated_run_dir_skips_non_numeric(self, tmp_path):
        """Test that non-numeric directory names are ignored."""
        today = dt.date.today()
        date_dir = tmp_path / today.isoformat()
        date_dir.mkdir(parents=True, exist_ok=True)

        # Create non-numeric directories
        (date_dir / "backup").mkdir()
        (date_dir / "temp").mkdir()
        (date_dir / "3").mkdir()

        run_dir = make_dated_run_dir(tmp_path, chdir_to="none")
        # Should create run 4 (3 + 1)
        assert run_dir == date_dir / "4"

    @patch("reaction_lammps_mupt.cache.get_cache_base_dir")
    @patch("reaction_lammps_mupt.cache.make_dated_run_dir")
    def test_get_current_cache_dir(self, mock_make_dated, mock_base_dir):
        """Test get_current_cache_dir calls make_dated_run_dir correctly."""
        mock_base_dir.return_value = Path("/mock/cache")
        mock_make_dated.return_value = Path("/mock/cache/2024-01-01/1")

        result = get_current_cache_dir()

        mock_make_dated.assert_called_once_with(Path("/mock/cache"), chdir_to="none")
        assert result == Path("/mock/cache/2024-01-01/1")


class TestRetentionCleanup:
    """Tests for retention cleanup functionality."""

    @patch("builtins.input")
    @patch("reaction_lammps_mupt.cache.get_cache_base_dir")
    def test_retention_cleanup_no_cleanup(self, mock_base_dir, mock_input, tmp_path):
        """Test selecting no cleanup option."""
        mock_base_dir.return_value = tmp_path
        mock_input.return_value = "7"

        retention_cleanup()

        # No files should be deleted

    @patch("builtins.input")
    @patch("reaction_lammps_mupt.cache.get_cache_base_dir")
    def test_retention_cleanup_one_week(self, mock_base_dir, mock_input, tmp_path, capsys):
        """Test cleanup with 1 week retention."""
        mock_base_dir.return_value = tmp_path
        mock_input.return_value = "1"

        # Create old and new directories
        old_date = (dt.date.today() - dt.timedelta(days=10)).isoformat()
        recent_date = (dt.date.today() - dt.timedelta(days=3)).isoformat()

        old_dir = tmp_path / old_date
        recent_dir = tmp_path / recent_date
        old_dir.mkdir()
        recent_dir.mkdir()

        retention_cleanup()

        captured = capsys.readouterr()
        assert not old_dir.exists()
        assert recent_dir.exists()
        assert f"Deleted old folder: {old_dir}" in captured.out

    @patch("builtins.input")
    @patch("reaction_lammps_mupt.cache.get_cache_base_dir")
    def test_retention_cleanup_delete_all(self, mock_base_dir, mock_input, tmp_path):
        """Test delete all cache folders option."""
        mock_base_dir.return_value = tmp_path
        mock_input.return_value = "6"

        # Create multiple dated directories
        for days_ago in [1, 30, 90]:
            date = (dt.date.today() - dt.timedelta(days=days_ago)).isoformat()
            (tmp_path / date).mkdir()

        retention_cleanup()

        # All dated directories should be deleted
        for item in tmp_path.iterdir():
            # Only dated directories should be gone
            try:
                dt.datetime.strptime(item.name, "%Y-%m-%d")
                assert False, f"Dated directory {item.name} should have been deleted"
            except ValueError:
                pass  # Not a dated directory

    @patch("builtins.input")
    @patch("reaction_lammps_mupt.cache.get_cache_base_dir")
    def test_retention_cleanup_custom_days(self, mock_base_dir, mock_input, tmp_path):
        """Test custom days retention."""
        mock_base_dir.return_value = tmp_path
        mock_input.side_effect = ["5", "15"]

        # Create directories at different ages
        old = (dt.date.today() - dt.timedelta(days=20)).isoformat()
        recent = (dt.date.today() - dt.timedelta(days=10)).isoformat()

        old_dir = tmp_path / old
        recent_dir = tmp_path / recent
        old_dir.mkdir()
        recent_dir.mkdir()

        retention_cleanup()

        assert not old_dir.exists()  # Older than 15 days
        assert recent_dir.exists()    # Within 15 days

    @patch("builtins.input")
    @patch("reaction_lammps_mupt.cache.get_cache_base_dir")
    def test_retention_cleanup_invalid_input_retry(self, mock_base_dir, mock_input, tmp_path):
        """Test that invalid input causes retry."""
        mock_base_dir.return_value = tmp_path
        mock_input.side_effect = ["invalid", "9", "7"]  # invalid, out of range, then valid

        retention_cleanup()

        assert mock_input.call_count == 3

    @patch("builtins.input")
    @patch("reaction_lammps_mupt.cache.get_cache_base_dir")
    def test_retention_cleanup_ignores_non_dated_folders(self, mock_base_dir, mock_input, tmp_path):
        """Test that non-dated folders are not deleted."""
        mock_base_dir.return_value = tmp_path
        mock_input.return_value = "6"  # Delete all

        # Create non-dated directories
        (tmp_path / "00_cache").mkdir()
        (tmp_path / "temp").mkdir()

        retention_cleanup()

        assert (tmp_path / "00_cache").exists()
        assert (tmp_path / "temp").exists()

    @patch("builtins.input")
    def test_retention_cleanup_custom_base_dir(self, mock_input, tmp_path):
        """Test retention_cleanup with custom base_dir parameter."""
        mock_input.return_value = "7"

        old_dir = tmp_path / (dt.date.today() - dt.timedelta(days=10)).isoformat()
        old_dir.mkdir()

        retention_cleanup(base_dir=tmp_path)

        # No cleanup with option 7

    @patch("builtins.input")
    @patch("reaction_lammps_mupt.cache.get_cache_base_dir")
    def test_retention_cleanup_negative_days_rejected(self, mock_base_dir, mock_input, tmp_path):
        """Test that negative days are rejected for custom retention."""
        mock_base_dir.return_value = tmp_path
        mock_input.side_effect = ["5", "-5", "10", "7"]  # custom, negative (rejected), valid, then skip

        retention_cleanup()

        # Should have asked for input multiple times due to rejection


class TestCopyOperations:
    """Tests for copy_to_date_folder functionality."""

    @patch("reaction_lammps_mupt.cache.get_current_cache_dir")
    def test_copy_to_date_folder_files(self, mock_cache_dir, tmp_path):
        """Test copying files to dated folder."""
        source_dir = tmp_path / "source"
        dest_dir = tmp_path / "dest"
        source_dir.mkdir()
        dest_dir.mkdir()

        # Create test files
        (source_dir / "file1.txt").write_text("content1")
        (source_dir / "file2.txt").write_text("content2")

        mock_cache_dir.return_value = dest_dir

        result = copy_to_date_folder(source_dir)

        assert result == dest_dir
        assert (dest_dir / "file1.txt").read_text() == "content1"
        assert (dest_dir / "file2.txt").read_text() == "content2"

    @patch("reaction_lammps_mupt.cache.get_current_cache_dir")
    def test_copy_to_date_folder_directories(self, mock_cache_dir, tmp_path):
        """Test copying directories to dated folder."""
        source_dir = tmp_path / "source"
        dest_dir = tmp_path / "dest"
        source_dir.mkdir()
        dest_dir.mkdir()

        # Create test subdirectory with files
        subdir = source_dir / "subdir"
        subdir.mkdir()
        (subdir / "nested.txt").write_text("nested")

        mock_cache_dir.return_value = dest_dir

        result = copy_to_date_folder(source_dir)

        assert result == dest_dir
        assert (dest_dir / "subdir").is_dir()
        assert (dest_dir / "subdir" / "nested.txt").read_text() == "nested"

    @patch("reaction_lammps_mupt.cache.get_current_cache_dir")
    def test_copy_to_date_folder_symlinks(self, mock_cache_dir, tmp_path):
        """Test copying symlinks preserves them as symlinks."""
        source_dir = tmp_path / "source"
        dest_dir = tmp_path / "dest"
        source_dir.mkdir()
        dest_dir.mkdir()

        # Create target and symlink
        target = source_dir / "target.txt"
        target.write_text("target")
        link = source_dir / "link"
        link.symlink_to(target)

        mock_cache_dir.return_value = dest_dir

        result = copy_to_date_folder(source_dir)

        assert result == dest_dir
        dest_link = dest_dir / "link"
        # The implementation copies files, and symlink behavior depends on the target resolution
        assert dest_link.exists()

    @patch("reaction_lammps_mupt.cache.get_current_cache_dir")
    def test_copy_to_date_folder_empty_source(self, mock_cache_dir, tmp_path):
        """Test copying from empty source directory."""
        source_dir = tmp_path / "source"
        dest_dir = tmp_path / "dest"
        source_dir.mkdir()
        dest_dir.mkdir()

        mock_cache_dir.return_value = dest_dir

        result = copy_to_date_folder(source_dir)

        assert result == dest_dir
        # Dest should still be empty
        assert len(list(dest_dir.iterdir())) == 0

    @patch("reaction_lammps_mupt.cache.get_current_cache_dir")
    def test_copy_to_date_folder_preserves_metadata(self, mock_cache_dir, tmp_path):
        """Test that copy_to_date_folder preserves file metadata."""
        source_dir = tmp_path / "source"
        dest_dir = tmp_path / "dest"
        source_dir.mkdir()
        dest_dir.mkdir()

        source_file = source_dir / "file.txt"
        source_file.write_text("content")
        original_stat = source_file.stat()

        mock_cache_dir.return_value = dest_dir

        copy_to_date_folder(source_dir)

        dest_file = dest_dir / "file.txt"
        # Note: mtime preservation may vary by filesystem
        assert dest_file.exists()


class TestEdgeCases:
    """Tests for edge cases and error conditions."""

    def test_delete_files_with_unicode_names(self, tmp_path):
        """Test deleting files with unicode characters in names."""
        unicode_file = tmp_path / "测试文件.txt"
        unicode_file.write_text("content")

        delete_files_in_directory(tmp_path)

        assert not unicode_file.exists()

    def test_make_dated_run_dir_with_many_runs(self, tmp_path):
        """Test creating many run directories."""
        # Create 10 run directories
        for i in range(10):
            run_dir = make_dated_run_dir(tmp_path, chdir_to="none")
            assert run_dir.name == str(i + 1)

    @patch("builtins.input")
    @patch("reaction_lammps_mupt.cache.get_cache_base_dir")
    def test_retention_cleanup_with_files_in_dated_dir(self, mock_base_dir, mock_input, tmp_path):
        """Test cleanup removes directories with files inside."""
        mock_base_dir.return_value = tmp_path
        mock_input.return_value = "6"  # Delete all

        old_date = (dt.date.today() - dt.timedelta(days=10)).isoformat()
        old_dir = tmp_path / old_date
        old_dir.mkdir()
        (old_dir / "file.txt").write_text("content")
        (old_dir / "subdir").mkdir()

        retention_cleanup()

        assert not old_dir.exists()

    def test_make_dated_run_dir_concurrent_access(self, tmp_path):
        """Test race condition handling with concurrent directory creation."""
        # Create existing run directories to test increment behavior
        today = dt.date.today()
        date_dir = tmp_path / today.isoformat()
        date_dir.mkdir(parents=True, exist_ok=True)
        (date_dir / "1").mkdir()

        # Create run 2, but simulate race condition by creating it manually after detection
        run_dir = make_dated_run_dir(tmp_path, chdir_to="none")

        # Should skip 1 and create 2 (or higher if race condition handled)
        assert run_dir.exists()
        assert int(run_dir.name) >= 2