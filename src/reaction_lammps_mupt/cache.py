from __future__ import annotations

import subprocess
import datetime
import os
import shutil
from pathlib import Path
from typing import Literal

import os

def delete_files_in_directory(directory_path):
    # Ensure the directory exists
    if not os.path.isdir(directory_path):
        print(f"Error: Directory not found at {directory_path}")
        return

    # List all entries in the directory
    for filename in os.listdir(directory_path):
        file_path = os.path.join(directory_path, filename)

        # Check if the path points to a file and not a subdirectory (optional but safe)
        if os.path.isfile(file_path) or os.path.islink(file_path):
            try:
                os.remove(file_path)
                pass  # Successfully deleted the file
            except OSError as e:
                pass  # Handle errors (e.g., permission issues) if needed
        else:
            pass  # If it's a directory, we can choose to ignore it or handle it as needed

def get_git_root() -> Path:
    out = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip()
    return Path(out)

def make_dated_run_dir(
    base_dir: Path,
    *,
    date: datetime.date | None = None,
    chdir_to: Literal["run", "date", "none"] = "run",
) -> Path:
    """
    Create a date-based folder under base_dir, then create a numeric run folder inside it.

    Structure:
        base_dir/
            YYYY-MM-DD/
                1/
                2/
                3/
                ...

    Args:
        base_dir: Root cache dir (e.g., Path("lunar/cache"))
        date: Specific date to use (defaults to today)
        chdir_to:
            - "run": chdir into the newly created run folder (default)
            - "date": chdir into the YYYY-MM-DD folder
            - "none": don't change cwd

    Returns:
        Path to the newly created run folder (base_dir/YYYY-MM-DD/N).
    """
    base_dir = Path(base_dir)
    base_dir.mkdir(parents=True, exist_ok=True)

    target_date = date or datetime.date.today()
    date_dir = base_dir / target_date.isoformat()
    date_dir.mkdir(parents=True, exist_ok=True)

    # Find next integer folder name: 1..infinite
    max_run_number = 0
    for p in date_dir.iterdir():
        if p.is_dir() and p.name.isdigit():
            run_number = int(p.name)
            if run_number > max_run_number:
                max_run_number = run_number

    next_run_number = max_run_number + 1

    # Create run_dir safely (handles rare race condition)
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
    Get the current cache directory based on today's date and existing run folders.

    Returns:
        Path to the current run folder (base_dir/YYYY-MM-DD/N).
    """
    base_dir = get_git_root() 
    base_dir = Path(base_dir) / "cache"
    base_dir.mkdir(parents=True, exist_ok=True)
    current_cache_dir = make_dated_run_dir(base_dir, chdir_to="none")
    return current_cache_dir

def get_default_cache_dir() -> Path:
    """
    Get the current cache directory and change cwd to it.

    Returns:
        Path to the current run folder (base_dir/YYYY-MM-DD/N).
    """
    base_dir = get_git_root() 
    base_dir = Path(base_dir) / "cache" / "0_cache"
    base_dir.mkdir(parents=True, exist_ok=True)
    delete_files_in_directory(base_dir)
    return base_dir

def retention_cleanup(base_dir: Path | None = None) -> None:
    """
    Ask user whether to clear old dated folders in a cache directory.

    Deletes date folders named YYYY-MM-DD (and everything inside them).
    """
    if base_dir is None:
        base_dir = get_git_root() / "cache"
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

    retention_map = {
        "1": datetime.timedelta(weeks=1),
        "2": datetime.timedelta(days=30),
        "3": datetime.timedelta(days=90),
        "4": datetime.timedelta(days=180),
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
                retention = datetime.timedelta(days=days)
                break
            except ValueError:
                print("Invalid number. Please enter an integer.")

    now = datetime.date.today()
    base_dir = Path(base_dir)

    for folder in base_dir.iterdir():
        if not folder.is_dir():
            continue
        try:
            folder_date = datetime.datetime.strptime(folder.name, "%Y-%m-%d").date()
        except ValueError:
            continue  # ignore non-date folders

        if retention == "all" or (now - folder_date > retention):
            shutil.rmtree(folder)
            print(f"Deleted old folder: {folder}")
            
def copy_to_date_folder(default_cache_dir) -> None:
    """
    Copy all files from base_dir to a date-based folder (base_dir/YYYY-MM-DD/).

    Args:
        base_dir: Directory containing files to move.
        date: Date to use for the target folder name.
    """
    current_cache_dir = get_current_cache_dir()
    for item in default_cache_dir.iterdir():
        target = current_cache_dir / item.name
        if item.is_file():
            shutil.copy(str(item), str(target))
        elif item.is_dir():
            shutil.copy(str(item), str(target))
    return current_cache_dir


if __name__ == "__main__":
    # Example usage
    retention_cleanup()
    default_cache_dir = get_default_cache_dir()
    print(f"Default cache directory (cleaned): {default_cache_dir}")
    # Example inputs for testing
    for i in range(5):
        with open(default_cache_dir / f"test_file_{i}.txt", "w") as f:
            f.write(f"This is test file {i}\n")
    current_cache_dir = copy_to_date_folder(default_cache_dir)
    print(f"Current cache directory (with files moved):from [{default_cache_dir}] to [{current_cache_dir}]")

