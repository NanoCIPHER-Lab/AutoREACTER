from __future__ import annotations

import datetime as dt
import os
import shutil
import subprocess
from pathlib import Path
from typing import Literal


def get_git_root() -> Path:
    """
    Locate the root directory of the current Git repository.
    
    Returns:
        Path: Path to the repository root directory.
    
    Raises:
        subprocess.CalledProcessError: If the current working directory is not inside a Git repository.
    """
    out = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip()
    return Path(out)


def get_cache_base_dir() -> Path:
    """
    Get the repository cache base directory at <git_root>/cache.
    
    Creates the directory and any missing parent directories if it does not exist.
    
    Returns:
        Path: Path to the cache base directory under the repository root.
    """
    base_dir = get_git_root() / "cache"
    base_dir.mkdir(parents=True, exist_ok=True)
    return base_dir


def get_default_cache_dir() -> Path:
    """
    Get the staging cache directory under the repository cache, creating it if missing.
    
    Returns:
        Path: Path to the staging cache directory located at <git_root>/cache/00_cache.
    """
    default_dir = get_cache_base_dir() / "00_cache"
    default_dir.mkdir(parents=True, exist_ok=True)
    return default_dir

def delete_files_in_directory(directory_path: Path) -> None:
    """
    Remove regular files and symlinks that are direct children of the given directory while leaving subdirectories intact.
    
    If the provided path does not refer to an existing directory, an error message is printed and no removal is attempted. Individual unlink failures are ignored so the function continues processing other entries.
    
    Parameters:
        directory_path (Path | str): Path to the target directory whose immediate files and symlinks will be removed.
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
    Remove files and symlinks directly inside the staging cache directory.
    
    The staging cache directory is located at <git_root>/cache/00_cache. This function does not remove subdirectories or their contents.
    """
    delete_files_in_directory(get_default_cache_dir())


def make_dated_run_dir(
    base_dir: Path,
    *,
    date: dt.date | None = None,
    chdir_to: Literal["run", "date", "none"] = "run",
) -> Path:
    """
    Create a new dated run directory under base_dir with an incremented numeric run number.
    
    Parameters:
        base_dir (Path): Base directory where the dated folder (YYYY-MM-DD) will be created.
        date (datetime.date | None): Date to use for the dated folder; uses today's date when None.
        chdir_to (Literal["run", "date", "none"]): If "run", change the current working directory to the new run folder;
            if "date", change to the date folder; if "none", do not change the working directory.
    
    Returns:
        Path: Path to the newly created run directory (base_dir / YYYY-MM-DD / N).
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


def get_current_cache_dir() -> Path:
    """
    Selects or creates the next dated run directory under the repository cache without changing the current working directory.
    
    Returns:
        Path: The path to the created or selected dated run directory (format: <git_root>/cache/YYYY-MM-DD/<N>).
    """
    return make_dated_run_dir(get_cache_base_dir(), chdir_to="none")


def retention_cleanup(base_dir: Path | None = None) -> None:
    """
    Interactively remove dated cache directories under the cache base or a provided base directory.
    
    Prompts the user to choose a retention policy (1 week, 1 month, 3 months, 6 months, custom number of days, delete all, or skip).
    Only directories whose names match the YYYY-MM-DD pattern are considered; directories older than the chosen retention (or all directories when "delete all" is selected) are removed with shutil.rmtree and a deletion message is printed. If no base_dir is provided, the cache base directory returned by get_cache_base_dir() is used.
    
    Parameters:
        base_dir (Path | None): Base directory containing dated cache folders. If None, uses the repository cache base.
    """
    base_dir = Path(base_dir) if base_dir is not None else get_cache_base_dir()

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

def copy_to_date_folder(source_dir: Path) -> Path:
    """
    Copy the contents of source_dir into a newly created dated run folder under the cache and return the destination path.
    
    Files are copied with metadata preserved, directories are copied recursively, and symbolic links are recreated as symlinks at the destination.
    
    Parameters:
        source_dir (Path | str): Source directory whose immediate contents will be copied into the new dated run folder.
    
    Returns:
        Path: Path to the created dated run directory (cache/YYYY-MM-DD/<N>).
    """
    source_dir = Path(source_dir)
    dest_dir = get_current_cache_dir()

    for item in source_dir.iterdir():
        target = dest_dir / item.name

        if item.is_file():
            # copy2 preserves metadata (mtime, etc.) which is often desirable
            shutil.copy2(item, target)

        elif item.is_dir():
            # shutil.copytree is the correct way to copy a directory
            shutil.copytree(item, target)

        elif item.is_symlink():
            # Preserve symlinks as symlinks (not the file contents)
            target.symlink_to(item.readlink())

    return dest_dir

if __name__ == "__main__":
    # Optional: clean up old dated cache folders interactively
    retention_cleanup()

    # Staging cache directory where your code dumps temporary outputs
    default_cache_dir = get_default_cache_dir()
    print(f"Default cache directory: {default_cache_dir}")

    # Example: create some test files in staging
    for i in range(5):
        (default_cache_dir / f"test_file_{i}.txt").write_text(f"This is test file {i}\n")

    # Copy staging contents into a fresh dated run folder
    current_cache_dir = copy_to_date_folder(default_cache_dir)
    print(f"Copied from [{default_cache_dir}] to [{current_cache_dir}]")