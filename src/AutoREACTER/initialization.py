import importlib
import sys
try:  # for running as a package
    from .cache import delete_default_cache_files as delete_cache_dir
except ImportError:  # for running as a script
    from cache import delete_default_cache_files as delete_cache_dir


def ASCII_Mupt_reaction_LAMMPS():
    ascii_art = r"""                                                                                                                                                                                                                                                       
      .o.                       .                  ooooooooo.   oooooooooooo       .o.         .oooooo.   ooooooooooooo oooooooooooo ooooooooo.   
     .888.                    .o8                  `888   `Y88. `888'     `8      .888.       d8P'  `Y8b  8'   888   `8 `888'     `8 `888   `Y88. 
    .8"888.     oooo  oooo  .o888oo  .ooooo.        888   .d88'  888             .8"888.     888               888       888          888   .d88' 
   .8' `888.    `888  `888    888   d88' `88b       888ooo88P'   888oooo8       .8' `888.    888               888       888oooo8     888ooo88P'  
  .88ooo8888.    888   888    888   888   888       888`88b.     888    "      .88ooo8888.   888               888       888    "     888`88b.    
 .8'     `888.   888   888    888 . 888   888       888  `88b.   888       o  .8'     `888.  `88b    ooo       888       888       o  888  `88b.  
o88o     o8888o  `V88V"V8P'   "888" `Y8bod8P'      o888o  o888o o888ooooood8 o88o     o8888o  `Y8bood8P'      o888o     o888ooooood8 o888o  o888o 
"""
    print(ascii_art)


def initialize():
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
