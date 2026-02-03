import importlib
import sys
try: # for running as a package
    from cache import delete_default_cache_files as delete_cache_dir
except ImportError: # for running as a script
    from reaction_lammps_mupt.cache import delete_default_cache_files as delete_cache_dir


def ASCII_Mupt_reaction_LAMMPS():
    """
    Display the ASCII art banner for Mupt reaction LAMMPS.
    
    Prints a multi-line ASCII banner to standard output.
    """
    ascii_art = r"""
.___  ___.  __    __  .______   .___________.              .______       _______     ___       ______ .___________. __    ______   .__   __. 
|   \/   | |  |  |  | |   _  \  |           |              |   _  \     |   ____|   /   \     /      ||           ||  |  /  __  \  |  \ |  | 
|  \  /  | |  |  |  | |  |_)  | `---|  |----`    ______    |  |_)  |    |  |__     /  ^  \   |  ,----'`---|  |----`|  | |  |  |  | |   \|  | 
|  |\/|  | |  |  |  | |   ___/      |  |        |______|   |      /     |   __|   /  /_\  \  |  |         |  |     |  | |  |  |  | |  . `  | 
|  |  |  | |  `--'  | |  |          |  |                   |  |\  \----.|  |____ /  _____  \ |  `----.    |  |     |  | |  `--'  | |  |\   | 
|__|  |__|  \______/  | _|          |__|                   | _| `._____||_______/__/     \__\ \______|    |__|     |__|  \______/  |__| \__| 
"""
    print(ascii_art)


def initialize():
    """
    Verify required dependencies are present, display the startup banner, and remove default cache files.
    
    Checks that the `rdkit`, `pandas`, and `numpy` modules can be imported. If all are available, prints a confirmation message and displays the ASCII banner via ASCII_Mupt_reaction_LAMMPS(). If any required module is missing, exits the process with an error message naming the missing module. After successful checks, removes default cache files by calling delete_cache_dir().
    """
    try:
        for module in ['rdkit', 'pandas', 'numpy']:
            importlib.import_module(module)
        print("All required modules are successfully imported.")
        ASCII_Mupt_reaction_LAMMPS()
    except ModuleNotFoundError as e:
        message = (
            f"ERROR: Required module not found: {e.name}\n\n"
            "👉 Please install the missing module and try again.\n\n"
            "Exiting program..."
        )
        sys.exit(message)  
    
    delete_cache_dir()

if __name__ == "__main__":
    initialize()