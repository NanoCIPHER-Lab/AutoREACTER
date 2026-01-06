from platformdirs import user_cache_dir
import os, datetime, shutil
from pathlib import Path


app_name = "ReactionLammpsMuPT"
app_author = "" # Optional, can be omitted on some platforms

cache_dir = user_cache_dir(f"{app_name}", app_author)

# Ensure the directory exists
os.makedirs(cache_dir, exist_ok=True)
lunar_cache_dir = os.path.join(cache_dir, "lunar")
os.makedirs(lunar_cache_dir, exist_ok=True)

# Use the cache_dir path for storing files (e.g., SQLite DB, pickle files)
print(f"Cache location: {cache_dir}")
"""
TODO:
# from here needs to be added to the functions
# in later steps we will handle the lammps scripts also
# becasue that part havent fininshed yet these functions will no be called anywhere yet
# if user did not specifcally mention the their desired location we will use this location to give them inputs 
# cache_dir = user_cache_dir(f"{app_name}", app_author)
# if this location is the working location we will clean the folder before run
# not the cache folder the folder before that
"""
def clear_cache_folder(cache_folder: Path) -> None:
        """
        Clear all files in script_inputs/000_cache without removing the folder itself.
        """
        if not cache_folder.exists():
            return

        for item in cache_folder.iterdir():
            try:
                if item.is_file():
                    item.unlink()  # delete file
                elif item.is_dir():
                    shutil.rmtree(item)  # delete directory recursively
            except Exception as e:
                print(f"Failed to remove {item}: {e}")

def retention_cleanup(base_dir: Path) -> None:
        """
        Ask user whether to clear old dated folders in a cache directory.

        Args:
            base_dir (Path): Path to cache folder (e.g., lunar/cache).
        """
        print(
            "\nNOTE: Cached files are stored in lunar/cache, organized by date. "
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
        retention = retention_map.get(choice)

        if retention is None:
            print("Skipping cleanup.")
            return

        if retention == "custom":
            while True:
                try:
                    days = int(input("Enter number of days to retain: ").strip())
                    retention = datetime.timedelta(days=days)
                    break
                except ValueError:
                    print("Invalid number. Please enter an integer.")

        now = datetime.date.today()

        for folder in Path(base_dir).iterdir():
            if folder.is_dir():
                try:
                    folder_date = datetime.datetime.strptime(folder.name, "%Y-%m-%d").date()
                except ValueError:
                    continue

                if retention == "all" or (now - folder_date > retention):
                    shutil.rmtree(folder)
                    print(f"Deleted old folder: {folder}")
