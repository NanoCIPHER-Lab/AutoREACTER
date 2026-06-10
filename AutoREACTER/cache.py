from pathlib import Path
import shutil
import tempfile
import uuid
from AutoREACTER.reaction_preparation.ff_wrapper.REACTER_files_builder import REACTERFiles


class GetCacheDir:
    def __init__(self, clear_staging: bool = True):
        """
        Initializes the cache manager.
        
        Args:
            clear_staging: 
                If True, clear all contents inside the staging directory.
                The staging directory itself is preserved.
        """
        # Generate a unique staging directory for this run to prevent concurrent conflicts
        run_id = uuid.uuid4().hex[:8]
        # self.staging_dir = Path(tempfile.gettempdir()) / f"AutoREACTER_staging_{run_id}"
        self.staging_dir = Path(tempfile.gettempdir()) / f"AutoREACTER_staging"
        self.staging_dir.mkdir(parents=True, exist_ok=True)
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
    
    def remove_path(self, path: Path) -> None:
        """
        Remove a file, symlink, or directory recursively.
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

