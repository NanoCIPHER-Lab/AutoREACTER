from platformdirs import user_cache_dir
import os, datetime, shutil
from pathlib import Path

current_dir = Path(__file__).parent.parent.resolve()
# should implement this for now we will use the project directory
# cache_dir = user_cache_dir("ReactionLammpsMuPT", "MUPT") 
os.makedirs(current_dir / "lunar" / "cache", exist_ok=True)
cache_dir = current_dir / "lunar" / "cache"
print(f"[cache] current_dir   = {current_dir}")
print(f"[cache] cache_dir     = {cache_dir}")

APP_NAME = "ReactionLammpsMuPT"
APP_AUTHOR = "MUPT"  # fine
cache_dir = Path(user_cache_dir(APP_NAME, APP_AUTHOR))
lunar_cache_dir = cache_dir / "lunar"

lunar_cache_dir.mkdir(parents=True, exist_ok=True)

print(f"[cache] cache_dir      = {cache_dir}")
print(f"[cache] lunar_cache_dir= {lunar_cache_dir}")


"""
TODO:
# from here needs to be added to the functions
# in later steps we will handle the lammps scripts also
# becasue that part havent finished yet these functions will no be called anywhere yet
# if user did not specifically mention the their desired location we will use this location to give them inputs 
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
