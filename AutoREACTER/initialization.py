"""
Author: "Janitha Mahanthe"
Version: v0.1-beta.0
"""
version = "v0.2-beta.0"
import importlib


class Initialization:
    def __init__(self):
        self.moldule_imports()
        self.ASCII_Mupt_reaction_LAMMPS()
        self.print_version()
        

    @classmethod
    def moldule_imports(cls):
        try:
            for module in ['rdkit', 'pandas', 'numpy', 'networkx']:
                importlib.import_module(module)
            print("All required modules are successfully imported.")
        except ModuleNotFoundError as e:
            message = (
                f"ERROR: Required module not found: {e.name}\n\n"
                "👉 Please install the missing module and try again.\n\n"
                "Exiting program..."
            )
            raise RuntimeError(message)
        
    @classmethod
    def ASCII_Mupt_reaction_LAMMPS(cls):
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
    
    @classmethod
    def print_version(cls):
        print(f"AutoREACTER version: {version}")

if __name__ == "__main__":
    Initialization()