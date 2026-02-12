from __future__ import annotations

from dataclasses import dataclass
import datetime as dt
import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import Literal
from venv import logger

logger = logging.getLogger(__name__)

@dataclass(slots=True)
class GetCacheDir:
    """Repocontext provides methods to get paths relative to the git repository root, such as cache directories."""

    # Note: For simplicity, this implementation does not actually change the working directory.
    # Instead, it provides methods to get paths relative to the repo root. This avoids issues
    # with changing the global working directory in a multi-threaded or complex application.

    @staticmethod
    def get_git_root() -> Path:
        """
        Return the root directory of the current git repo.
        Raises subprocess.CalledProcessError if not in a git repo.
        """
        out = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip()
        if not out:
            script_dir = Path(__file__).resolve().parent
            out = script_dir.parent.parent.parent  # Fallback: <git_root>/src/reaction_lammps_mupt/cache.py
            logger.warning(f"Warning: 'git rev-parse' returned empty. Falling back to {out} as git root based on script location.")
        return Path(out)
        
    @property
    def cache_base_dir(self) -> Path:
        """
        <git_root>/cache (created if missing)
        """
        base = self.get_git_root() / "cache"
        base.mkdir(parents=True, exist_ok=True)
        return base
    
    @property
    def staging_dir(self) -> Path:
        """
        <git_root>/cache/00_cache (created if missing)

        NOTE: This does NOT auto-delete files. Call clear_staging_dir() explicitly.
        """
        d = self.cache_base_dir / "00_cache"
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d, exist_ok=True)
        return d

@dataclass(slots=True, frozen=True)
class RunDirectoryManager:
    """
    Creates dated run folders and copies artifacts into an existing run folder.

    Key behavior:
    - You create a run folder ONCE (make_dated_run_dir).
    - You can copy into that same run folder unlimited times (copy_into_run).
    - Copy semantics are REPLACE:
        * files overwrite
        * directories are deleted then copied fresh
        * symlinks are replaced
    """

    @staticmethod
    def make_dated_run_dir(
        base_dir: Path,
        *,
        date: dt.date | None = None,
        chdir_to: Literal["run", "date", "none"] = "none",
    ) -> Path:
        """
        Create:
            base_dir/YYYY-MM-DD/<N>/

        Where <N> is the next integer run folder (1, 2, 3, ...).

        chdir_to:
          - "run":  chdir into the new run folder
          - "date": chdir into the date folder
          - "none": do not change cwd
        """
        base_dir = Path(base_dir)
        base_dir.mkdir(parents=True, exist_ok=True)

        target_date = date or dt.date.today()
        date_dir = base_dir / target_date.isoformat()
        date_dir.mkdir(parents=True, exist_ok=True)

        max_run_number = 0
        for p in date_dir.iterdir():
            if p.is_dir() and p.name.isdigit():
                max_run_number = max(max_run_number, int(p.name))

        next_run_number = max_run_number + 1
        while True:
            run_dir = date_dir / str(next_run_number)
            try:
                run_dir.mkdir()
                break
            except FileExistsError:
                next_run_number += 1

        if chdir_to == "run":
            os.chdir(run_dir)
        elif chdir_to == "date":
            os.chdir(date_dir)

        return run_dir

    @staticmethod
    def _remove_path(path: Path) -> None:
        """Remove file/dir/symlink at path if it exists."""
        if not path.exists() and not path.is_symlink():
            return

        # symlink must be checked before is_dir() in some cases
        if path.is_symlink():
            path.unlink()
            return

        if path.is_file():
            path.unlink()
            return

        if path.is_dir():
            shutil.rmtree(path)
            return

        # Fallback for odd filesystem types
        try:
            path.unlink()
        except Exception:
            shutil.rmtree(path, ignore_errors=True)

    @staticmethod
    def copy_into_run(
        source_dir: Path,
        dest_run_dir: Path,
        *,
        overwrite: bool = True,
        preserve_symlinks: bool = True,
        preserve_metadata: bool = True,
    ) -> Path:
        """
        Copy everything from source_dir into an EXISTING dest_run_dir.

        overwrite=True means:
          - if target exists, it is deleted first (dirs removed entirely),
            then copied fresh.

        Returns dest_run_dir for convenience.
        """
        source_dir = Path(source_dir)
        dest_run_dir = Path(dest_run_dir)

        if not source_dir.is_dir():
            raise NotADirectoryError(f"source_dir is not a directory: {source_dir}")

        dest_run_dir.mkdir(parents=True, exist_ok=True)

        for item in source_dir.iterdir():
            target = dest_run_dir / item.name

            if overwrite and (target.exists() or target.is_symlink()):
                RunDirectoryManager._remove_path(target)

            if item.is_symlink():
                if not preserve_symlinks:
                    # If you don't want symlinks, copy the linked contents instead
                    resolved = item.resolve()
                    if resolved.is_dir():
                        shutil.copytree(resolved, target)
                    else:
                        shutil.copy2(resolved, target) if preserve_metadata else shutil.copy(resolved, target)
                else:
                    target.symlink_to(item.readlink())

            elif item.is_file():
                if preserve_metadata:
                    shutil.copy2(item, target)
                else:
                    shutil.copy(item, target)

            elif item.is_dir():
                # copytree requires target not exist (we already removed it if overwrite)
                shutil.copytree(
                    item,
                    target,
                    symlinks=preserve_symlinks,
                    copy_function=shutil.copy2 if preserve_metadata else shutil.copy,
                )

            else:
                # Ignore special files (sockets, devices, etc.)
                # If you want, raise here instead.
                continue

        return dest_run_dir

@dataclass(slots=True)
class RetentionCleanup:
    @staticmethod
    def run(base_dir: Path | None = None) -> None:
        """
        Interactive cleanup of dated cache folders under base_dir (default: <git_root>/cache).

        Deletes folders whose names match YYYY-MM-DD and are older than chosen retention.
        """
        base_dir = Path(base_dir) if base_dir is not None else ()

        print(
            f"\nNOTE: Cached files are stored in {base_dir}, organized by date. "
            "You can clear old folders here to save space."
        )

        while True:
            choice = input(
                "\nDo you want to clear old files?\n"
                "1. 1 week\n"
                "2. 1 month\n"
                "3. 3 months\n"
                "4. 6 months\n"
                "5. Custom (enter number of days)\n"
                "6. Delete ALL cache folders\n"
                "7. No cleanup\n"
                "Enter choice (1-7): "
            ).strip()
            if choice in {"1", "2", "3", "4", "5", "6", "7"}:
                break
            print("Invalid choice. Please enter a number between 1 and 7.")

        retention_map: dict[str, dt.timedelta | str | None] = {
            "1": dt.timedelta(weeks=1),
            "2": dt.timedelta(days=30),
            "3": dt.timedelta(days=90),
            "4": dt.timedelta(days=180),
            "5": "custom",
            "6": "all",
            "7": None,
        }
        retention = retention_map[choice]

        if retention is None:
            print("Skipping cleanup.")
            return

        if retention == "custom":
            while True:
                try:
                    days = int(input("Enter number of days to retain: ").strip())
                    if days < 0:
                        print("Days must be >= 0.")
                        continue
                    retention = dt.timedelta(days=days)
                    break
                except ValueError:
                    print("Invalid number. Please enter an integer.")

        today = dt.date.today()

        for folder in base_dir.iterdir():
            if not folder.is_dir():
                continue

            # Only operate on folders named like YYYY-MM-DD
            try:
                folder_date = dt.datetime.strptime(folder.name, "%Y-%m-%d").date()
            except ValueError:
                continue

            should_delete = (retention == "all") or (today - folder_date > retention)
            if should_delete:
                shutil.rmtree(folder)
                print(f"Deleted old folder: {folder}")



if __name__ == "__main__":
    # Staging cache directory where your code dumps temporary outputs
    get_cache_dir = GetCacheDir()
    cache_dir = get_cache_dir.staging_dir
    print(f"Default cache directory: {cache_dir}")
    dated_cache_dir = RunDirectoryManager.make_dated_run_dir(get_cache_dir.cache_base_dir, chdir_to="none")
    print(f"Created dated run directory: {dated_cache_dir}")

    # # Example: create some test files in staging
    for i in range(5):
        (cache_dir / f"test_file_{i}.txt").write_text(f"This is test file {i}\n")
    run_dir = RunDirectoryManager.copy_into_run(cache_dir, dated_cache_dir)
    print(f"Created run directory: {run_dir}")

    RetentionCleanup.run(get_cache_dir.cache_base_dir)


