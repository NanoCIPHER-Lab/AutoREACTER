# LUNAR_client/api_wrapper.py
"""Main API wrapper for interacting with LUNAR."""

import subprocess
import os, sys
import subprocess
from pathlib import Path
import shutil
import datetime
import time

from .....input_parser import InputParser
print("Imported InputParser from input_parser.py")
from .....cache import cache_dir
print(f"Using cache directory: {cache_dir}")
if __package__ is None or __package__ == "":
    sys.path.append(os.path.dirname(os.path.dirname(__file__)).parent.parent)
    from LUNAR_client.locate_LUNAR import get_LUNAR_loc
    print("Imported locate_LUNAR directly.")
else:
    from .locate_LUNAR import get_LUNAR_loc
    print("Imported locate_LUNAR as a package.")


cwd = os.getcwd()
original_cwd = os.getcwd()
work_dir = Path(__file__).parent
os.chdir(work_dir)
LUNAR_LOCATION = get_LUNAR_loc(use_gui=True)

class LunarHandler:
    """
    Handles operations related to the 'lunar' subpackage, including:

    - Validating SMILES strings using RDKit
    - Managing cache directory cleanup based on retention policy
    - Moving/copying cache files into organized subfolders
    - Executing the 'api.py' script from the lunar subpackage
    """

    def __init__(self, number_of_monomers: int):
        """
        Initialize LunarHandler.

        Args:
            number_of_monomers (int): Number of SMILES monomers to collect and process.
        """

        self.clear_cache_folder()
        # Show loading animation for user feedback
        self.loading_screen()

        # Base path of this script (script_forge) and the lunar subpackage
        self.script_forge_path = Path(__file__).parent
        self.subpackage_path = self.script_forge_path / "lunar"

        # Clean up retention folders and reorganize cache
        self.retention_cleanup(self.subpackage_path / "cache")
        self.move_files(self.subpackage_path / "cache")

        # Path to api.py inside lunar
        subpackage_path = self.subpackage_path / "api.py"

        self.monomer_smiles = []
        for i in range(number_of_monomers):
            # Prompt and validate SMILES
            smiles_string = self.validate_smiles_rdkit(i + 1)
            self.monomer_smiles.append(smiles_string)

            # Call api.py with the validated SMILES string
            subprocess.run([sys.executable, subpackage_path, smiles_string])

            # Copy cache outputs into script_vault
            self.copy_to_cache_folder_script_vault(f"{smiles_string}_")

            # Add spacing after each iteration for readability
            if i >= 1:
                print("\n\n")
        os.chdir(cwd)

    def validate_smiles_rdkit(self, monomer_number: int) -> str:
        """
        Prompt user for a SMILES string and validate with RDKit.

        Args:
            monomer_number (int): The index of the current monomer being requested.

        Returns:
            str: A validated SMILES string.
        """
        try:
            from rdkit import Chem
        except ModuleNotFoundError:
            message = (
                "ERROR: RDKit is not installed in this environment.\n\n"
                "👉 To install RDKit, run one of the following:\n"
                "   conda install -c conda-forge rdkit   (recommended)\n"
                "   pip install rdkit-pypi              (alternative)\n\n"
                "Exiting program..."
            )
            sys.exit(message)

        while True:
            smiles_string = input(
                f"Please enter the SMILES string for monomer {monomer_number}: "
            )
            mol = Chem.MolFromSmiles(smiles_string)
            if mol is None:
                print("Invalid SMILES string. Please try again.")
            else:
                return smiles_string

    def retention_cleanup(self, base_dir: Path) -> None:
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

    def loading_screen(self, name="LUNAR") -> None:
        """Show ASCII art banner with loading spinner."""
        banner = r"""
    █████       █████  █████ ██████   █████   █████████   ███████████  
    ░░███       ░░███  ░░███ ░░██████ ░░███   ███░░░░░███ ░░███░░░░░███ 
    ░███        ░███   ░███  ░███░███ ░███  ░███    ░███  ░███    ░███ 
    ░███        ░███   ░███  ░███░░███░███  ░███████████  ░██████████  
    ░███        ░███   ░███  ░███ ░░██████  ░███░░░░░███  ░███░░░░░███ 
    ░███      █ ░███   ░███  ░███  ░░█████  ░███    ░███  ░███    ░███ 
    ███████████ ░░████████   █████  ░░█████ █████   █████ █████   █████
    ░░░░░░░░░░░   ░░░░░░░░   ░░░░░    ░░░░░ ░░░░░   ░░░░░ ░░░░░   ░░░░░ 
        """
        print(banner)
        print(f"Loading {name} ", end="", flush=True)

        animation = ["-", "\\", "|", "/"]
        for i in range(10):
            time.sleep(0.1)
            sys.stdout.write("\b" + animation[i % len(animation)])
            sys.stdout.flush()
        print("\nReady!")

    def move_files(self, base_dir: Path) -> None:
        """
        Move all files from a cache directory into today's dated subfolder.
        If duplicates exist, create subfolders 1, 2, 3... to prevent overwrites.
        """
        today = datetime.date.today().strftime("%Y-%m-%d")
        target_root = Path(base_dir) / today
        target_root.mkdir(parents=True, exist_ok=True)

        for item in Path(base_dir).iterdir():
            if item.is_file():
                subfolder = target_root / "1"
                while True:
                    subfolder.mkdir(parents=True, exist_ok=True)
                    dest_file = subfolder / item.name
                    if dest_file.exists():
                        subfolder = target_root / str(int(subfolder.name) + 1)
                        continue
                    shutil.move(str(item), str(dest_file))
                    break

    def copy_to_cache_folder_script_vault(self, monomer_string: str) -> None:
        """
        Copy cache files related to a given monomer into the script_vault/000_cache folder.
        """
        destination_folder = Path(__file__).parent.parent / "script_inputs"
        source_folder = self.script_forge_path / "lunar" / "cache"
        all_files = os.listdir(source_folder)

        folder_name = "000_cache"
        dest_path = destination_folder / folder_name
        dest_path.mkdir(parents=True, exist_ok=True)

        for filename in all_files:
            if filename.startswith(monomer_string):
                source_file = os.path.join(source_folder, filename)
                dest_file = os.path.join(dest_path, filename)

                # Copy the file first
                shutil.copy(source_file, dest_file)

                # If it matches the target suffix, rename it twice
                if filename.endswith("IFF_merged.lmpmol"):
                    base_name = monomer_string  # prefix for renamed files

                    # Build new filenames
                    typed_file = os.path.join(dest_path, f"{base_name}typed_IFF_merged.lmpmol")
                    molecule_file = os.path.join(dest_path, f"{base_name}.molecule")

                    # Rename (move) the copied file to the first name
                    os.rename(dest_file, typed_file)

                    # Also create a copy with .molecule extension
                    shutil.copy(typed_file, molecule_file)

                    # --- Only make .molecule copy in parent folder ---
                    parent_path = os.path.abspath(os.path.join(dest_path, ".."))
                    molecule_parent = os.path.join(parent_path, f"{base_name}.molecule")
                    shutil.copy(typed_file, molecule_parent)


                if filename.endswith("IFF_merged_coff.data"):
                    base_name = monomer_string  # prefix for renamed files

                    # Build new filenames
                    typed_file = os.path.join(dest_path, f"{base_name}IFF_merged_coff.data")
                    molecule_file = os.path.join(dest_path, f"{base_name}.data")

                    # Rename (move) the copied file to the first name
                    os.rename(dest_file, typed_file)

                    # Also create a copy with .molecule extension
                    shutil.copy(typed_file, molecule_file)

                    # --- Only make .data copy in parent folder ---
                    parent_path = os.path.abspath(os.path.join(dest_path, ".."))
                    molecule_parent = os.path.join(parent_path, f"{base_name}.data")
                    shutil.copy(typed_file, molecule_parent)

    def clear_cache_folder(self) -> None:
        """
        Clear all files in script_inputs/000_cache without removing the folder itself.
        """
        cache_folder = Path(__file__).parent.parent / "script_inputs" / "000_cache"

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

os.chdir(original_cwd)

# if __name__ == "__main__":
#     # Example usage: process 2 monomers
#     handler = LunarHandler(number_of_monomers=2)