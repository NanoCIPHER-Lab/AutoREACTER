from pathlib import Path
import os
import re
import shutil
import subprocess
import datetime as dt
from AutoREACTER.reaction_preparation.lunar_client.REACTER_files_builder import REACTERFiles

class GetCacheDir:
    """
    Manages the base cache directory for the AutoREACTER workflow.
    
    This class determines the root of the git repository (or falls back to a
    path relative to the script) and creates a standardized cache directory
    structure with a staging area.
    """

    def __init__(self):
        """Initialize cache directory structure."""
        self.git_root = self.get_git_root()
        self.cache_base_dir = self.git_root / "cache"
        self.cache_base_dir.mkdir(parents=True, exist_ok=True)

        # Staging directory for temporary files before moving to dated run folders
        self.staging_dir = self.cache_base_dir / "00_cache"
        os.makedirs(self.staging_dir, exist_ok=True)

    def get_git_root(self) -> Path:
        """
        Return the root directory of the current git repository.
        
        Falls back to a path derived from this script's location if git
        command fails (e.g., when running outside a git repository).
        
        Returns:
            Path: Root directory path
        """
        try:
            out = subprocess.check_output(
                ["git", "rev-parse", "--show-toplevel"],
                text=True
            ).strip()
            return Path(out)
        except Exception:
            script_dir = Path(__file__).resolve().parent
        
            for parent in [script_dir] + list(script_dir.parents):
                if (parent / ".git").exists() or (parent / "pyproject.toml").exists():
                    return parent
        
            return script_dir.parent


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


class RetentionCleanup:
    """
    Provides interactive cleanup functionality for old cache directories.
    
    Allows users to delete simulation run folders older than a chosen retention
    period (1 week, 1 month, 3 months) or delete everything.
    """

    def __init__(self, base_dir: Path):
        """
        Initialize cleanup manager with target base directory.
        
        Args:
            base_dir: Cache base directory to clean up
        """
        self.base_dir = Path(base_dir)

    def run(self, mode="skip"):
        """
        Run the cleanup process based on the specified mode.
        Modes are:
            - "skip": Do not perform any cleanup
            - "all": Delete all dated run directories
            - int (number of days): Delete directories older than this many days
        Args:
            mode: Cleanup mode - "skip", int(number of days), or "all"
        Returns:
            None
        """
        print(f"\n[INFO] Cache directory: {self.base_dir}")

        if mode == "skip":
            print("[INFO] Cache cleanup skipped.")
            return

        pattern = re.compile(r"\d{4}-\d{2}-\d{2}")
        protected = {"00_cache"}
        
        if mode == "all":
            for folder in self.base_dir.iterdir():
                if (
                    folder.is_dir()
                    and folder.name not in protected
                    and pattern.match(folder.name)
                ):
                    shutil.rmtree(folder)

        try:
            days = int(mode)
        except Exception:
            print(f"[WARN] Invalid cleanup mode: {mode}")
            return

        cutoff = dt.date.today() - dt.timedelta(days=days)

        for folder in self.base_dir.iterdir():
            if not folder.is_dir() or folder.name == "00_cache":
                continue

            try:
                folder_date = dt.datetime.strptime(folder.name, "%Y-%m-%d").date()
            except ValueError:
                continue

            if folder_date < cutoff:
                shutil.rmtree(folder)
                print(f"[OK] Deleted: {folder}")
