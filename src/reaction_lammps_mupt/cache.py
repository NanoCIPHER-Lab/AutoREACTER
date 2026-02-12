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

# class RunDirectoryManager:




def delete_files_in_directory(directory_path: Path) -> None:
    """
    Delete only files/symlinks directly inside directory_path.
    (Does NOT delete subfolders.)
    """
    directory_path = Path(directory_path)

    if not directory_path.is_dir():
        print(f"Error: Directory not found at {directory_path}")
        return

    for entry in directory_path.iterdir(): 
        if entry.is_file() or entry.is_symlink():
            try:
                entry.unlink()
            except OSError:
                pass


def delete_default_cache_files() -> None:
    """
    Convenience: clear out files from <git_root>/cache/0_cache.
    """
    repo_ctx = RepoContext()
    staging_dir = repo_ctx.staging_dir
    delete_files_in_directory(staging_dir)


def make_dated_run_dir(
    base_dir: Path,
    *,
    date: dt.date | None = None,
    chdir_to: Literal["run", "date", "none"] = "run",
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

    # Find max existing run folder number under this date folder
    max_run_number = 0
    for p in date_dir.iterdir():
        if p.is_dir() and p.name.isdigit():
            max_run_number = max(max_run_number, int(p.name))

    # Attempt to create the next run folder; handle rare race condition
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


# def get_current_cache_dir() -> Path:
#     """
#     Create (or pick next) dated run folder under <git_root>/cache and return it.
#     Does NOT change cwd.
#     """
#     return make_dated_run_dir(get_cache_base_dir(), chdir_to="none")


# def retention_cleanup(base_dir: Path | None = None) -> None:
#     """
#     Interactive cleanup of dated cache folders under base_dir (default: <git_root>/cache).

#     Deletes folders whose names match YYYY-MM-DD and are older than chosen retention.
#     """
#     base_dir = Path(base_dir) if base_dir is not None else get_cache_base_dir()

#     print(
#         f"\nNOTE: Cached files are stored in {base_dir}, organized by date. "
#         "You can clear old folders here to save space."
#     )

#     while True:
#         choice = input(
#             "\nDo you want to clear old files?\n"
#             "1. 1 week\n"
#             "2. 1 month\n"
#             "3. 3 months\n"
#             "4. 6 months\n"
#             "5. Custom (enter number of days)\n"
#             "6. Delete ALL cache folders\n"
#             "7. No cleanup\n"
#             "Enter choice (1-7): "
#         ).strip()
#         if choice in {"1", "2", "3", "4", "5", "6", "7"}:
#             break
#         print("Invalid choice. Please enter a number between 1 and 7.")

#     retention_map: dict[str, dt.timedelta | str | None] = {
#         "1": dt.timedelta(weeks=1),
#         "2": dt.timedelta(days=30),
#         "3": dt.timedelta(days=90),
#         "4": dt.timedelta(days=180),
#         "5": "custom",
#         "6": "all",
#         "7": None,
#     }
#     retention = retention_map[choice]

#     if retention is None:
#         print("Skipping cleanup.")
#         return

#     if retention == "custom":
#         while True:
#             try:
#                 days = int(input("Enter number of days to retain: ").strip())
#                 if days < 0:
#                     print("Days must be >= 0.")
#                     continue
#                 retention = dt.timedelta(days=days)
#                 break
#             except ValueError:
#                 print("Invalid number. Please enter an integer.")

#     today = dt.date.today()

#     for folder in base_dir.iterdir():
#         if not folder.is_dir():
#             continue

#         # Only operate on folders named like YYYY-MM-DD
#         try:
#             folder_date = dt.datetime.strptime(folder.name, "%Y-%m-%d").date()
#         except ValueError:
#             continue

#         should_delete = (retention == "all") or (today - folder_date > retention)
#         if should_delete:
#             shutil.rmtree(folder)
#             print(f"Deleted old folder: {folder}")

# def copy_to_date_folder(source_dir: Path) -> Path:
#     """
#     Copy everything from source_dir (e.g., cache/0_cache) into a newly created
#     dated run folder under cache/YYYY-MM-DD/N.

#     Returns the destination run folder path.
#     """
#     source_dir = Path(source_dir)
#     dest_dir = get_current_cache_dir()

#     for item in source_dir.iterdir():
#         target = dest_dir / item.name

#         if item.is_file():
#             # copy2 preserves metadata (mtime, etc.) which is often desirable
#             shutil.copy2(item, target)

#         elif item.is_dir():
#             # shutil.copytree is the correct way to copy a directory
#             shutil.copytree(item, target)

#         elif item.is_symlink():
#             # Preserve symlinks as symlinks (not the file contents)
#             target.symlink_to(item.readlink())

#     return dest_dir

if __name__ == "__main__":
    # Optional: clean up old dated cache folders interactively
    # retention_cleanup()

    # Staging cache directory where your code dumps temporary outputs
    get_cache_dir = GetCacheDir()
    cache_dir = get_cache_dir.staging_dir
    print(f"Default cache directory: {cache_dir}")

    # # Example: create some test files in staging
    # for i in range(5):
    #     (cache_dir / f"test_file_{i}.txt").write_text(f"This is test file {i}\n")

    # # Copy staging contents into a fresh dated run folder
    # current_cache_dir = copy_to_date_folder(default_cache_dir)
    # print(f"Copied from [{default_cache_dir}] to [{current_cache_dir}]")

