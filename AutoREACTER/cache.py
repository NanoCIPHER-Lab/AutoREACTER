from pathlib import Path
import sys
import re
import shutil
import subprocess
import datetime as dt
import tempfile
from typing import Optional, Union
from AutoREACTER.reaction_preparation.lunar_client.REACTER_files_builder import REACTERFiles


class GetCacheDir:
    """
    Manages the base cache directory for the AutoREACTER workflow.
    
    This class determines the root of the git repository (or falls back to a
    path relative to the execution location) and creates a standardized cache
    directory structure with a staging area.
    """

    def __init__(
        self,
        clear_staging: bool = False,
    ):
        """
        Initialize cache directory structure.

        Args:
            clear_staging:
                If True, clear all contents inside the staging directory.
                The staging directory itself is preserved.
            base_dir:
                Optional base directory to use instead of auto-detecting.
        """
        # Staging directory for temporary files before moving to dated run folders
        self.staging_dir = Path(tempfile.gettempdir()) / "AutoREACTER_staging"
        self.staging_dir.mkdir(parents=True, exist_ok=True)

        if clear_staging:
            self.clear_staging_dir()

    def clear_staging_dir(self) -> None:
        """
        Clear all files and folders inside the staging cache directory.

        This only clears the temporary staging directory contents.
        It does not delete dated run folders.
        """
        self.staging_dir.mkdir(parents=True, exist_ok=True)

        failed = []
        for item in self.staging_dir.iterdir():
            try:
                if item.is_symlink() or item.is_file():
                    item.unlink()
                elif item.is_dir():
                    shutil.rmtree(item)
            except Exception as e:
                print(f"[WARN] Failed to remove staging cache item {item}: {e}")
                failed.append(item)

        if failed:
            print(f"[WARN] Staging cache partially cleared ({len(failed)} item(s) could not be removed): {self.staging_dir}")
        else:
            print(f"[OK] Cleared staging cache: {self.staging_dir}")


class RunDirectoryManager:
    """
    Handles creation and management of dated run directories for simulations.
    
    Organizes output into a structure like: cache/YYYY-MM-DD/1/, cache/YYYY-MM-DD/2/, etc.
    Provides utilities for moving files and updating REACTER file references.
    """

    def __init__(self, base_dir: Path):
        """
        Initialize with a base directory for storing run folders.
        
        Args:
            base_dir: Base path where dated run directories will be created
        """
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(parents=True, exist_ok=True)

    def _is_empty(self, path: Path) -> bool:
        """Check if a directory contains no files or subdirectories."""
        return not any(path.iterdir())

    def make_dated_run_dir(self) -> Path:
        """
        Create and return a dated run directory following the pattern: {base}/{today}/{run_number}.
        
        Reuses the latest run directory if it is empty. Otherwise, increments the run number.
        
        Returns:
            Path: Path to the final run directory
        """
        today = dt.date.today().isoformat()
        date_dir = self.base_dir / today
        date_dir.mkdir(parents=True, exist_ok=True)

        # Find existing run numbers (directories named with integers)
        existing_runs = [
            int(p.name) for p in date_dir.iterdir()
            if p.is_dir() and p.name.isdigit()
        ]

        if existing_runs:
            latest_run = date_dir / str(max(existing_runs))
            if self._is_empty(latest_run):
                print(f"[INFO] Final directory: {latest_run}")
                return latest_run

        # Create new run directory with incremented number
        run_number = max(existing_runs, default=0) + 1
        run_dir = date_dir / str(run_number)
        run_dir.mkdir()

        print(f"[INFO] Final directory: {run_dir}")
        return run_dir

    def remove_path(self, path: Path) -> None:
        """
        Remove a file, symlink, or directory recursively.
        
        Args:
            path: Path to remove
        """
        if path.is_symlink() or path.is_file():
            path.unlink()
        elif path.is_dir():
            shutil.rmtree(path)

    def move_into_run(self, source_dir: Path, dest_dir: Path) -> Path:
        """
        Move all contents from source_dir into dest_dir, overwriting if necessary.
        
        Args:
            source_dir: Directory containing files to move
            dest_dir: Destination run directory
            
        Returns:
            Path: The destination directory
        """
        for item in source_dir.iterdir():
            target = dest_dir / item.name

            if target.exists():
                self.remove_path(target)

            shutil.move(str(item), str(target))

        print(f"[OK] Moved files → {dest_dir}")
        return dest_dir

    def move_reacter_files(
        self, 
        reacter_files: REACTERFiles, 
        staging_dir: Path, 
        final_dir: Path
    ) -> REACTERFiles:
        """
        Move REACTER-related files from staging area into the final run directory
        and update all internal Path references in the REACTERFiles object.
        
        Args:
            reacter_files: REACTERFiles object containing various file paths
            staging_dir: Current staging directory (usually .../00_cache)
            final_dir: Target run directory
            
        Returns:
            REACTERFiles: Updated object with remapped paths
            
        Raises:
            FileNotFoundError: If the expected REACTER files directory doesn't exist
        """
        old_base = staging_dir / "lunar" / "REACTER_files"

        if not old_base.exists():
            raise FileNotFoundError(f"REACTER files not found: {old_base}")

        new_base = self.move_into_run(old_base, final_dir)

        def remap(p: Path) -> Path:
            """Remap a path from old_base to new_base while preserving relative structure."""
            try:
                return new_base / p.relative_to(old_base)
            except ValueError:
                return p

        # Update all file references to point to the new location
        reacter_files.force_field_data = remap(reacter_files.force_field_data)
        reacter_files.in_file = remap(reacter_files.in_file)

        for mol in reacter_files.molecule_files:
            mol.molecule_files.lmp_molecule_file = remap(
                mol.molecule_files.lmp_molecule_file
            )

        for tpl in reacter_files.template_files:
            tpl.map_file = remap(tpl.map_file)
            tpl.pre_reaction_file.lmp_molecule_file = remap(
                tpl.pre_reaction_file.lmp_molecule_file
            )
            tpl.post_reaction_file.lmp_molecule_file = remap(
                tpl.post_reaction_file.lmp_molecule_file
            )

        print(f"[OK] REACTER files moved → {new_base}")
        return reacter_files

