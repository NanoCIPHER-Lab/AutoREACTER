import importlib
import sys
try:  # for running as a package
    from .cache import delete_default_cache_files as delete_cache_dir
except ImportError:  # for running as a script
    from cache import delete_default_cache_files as delete_cache_dir


def ASCII_Mupt_reaction_LAMMPS():
    ascii_art = r"""                                                                                                                                                                                                                                                       
                                                                                                                                                                                                                                                       
               AAA                                     tttt                                RRRRRRRRRRRRRRRRR   EEEEEEEEEEEEEEEEEEEEEE               AAA                  CCCCCCCCCCCCCTTTTTTTTTTTTTTTTTTTTTTTEEEEEEEEEEEEEEEEEEEEEERRRRRRRRRRRRRRRRR   
              A:::A                                 ttt:::t                                R::::::::::::::::R  E::::::::::::::::::::E              A:::A              CCC::::::::::::CT:::::::::::::::::::::TE::::::::::::::::::::ER::::::::::::::::R  
             A:::::A                                t:::::t                                R::::::RRRRRR:::::R E::::::::::::::::::::E             A:::::A           CC:::::::::::::::CT:::::::::::::::::::::TE::::::::::::::::::::ER::::::RRRRRR:::::R 
            A:::::::A                               t:::::t                                RR:::::R     R:::::REE::::::EEEEEEEEE::::E            A:::::::A         C:::::CCCCCCCC::::CT:::::TT:::::::TT:::::TEE::::::EEEEEEEEE::::ERR:::::R     R:::::R
           A:::::::::A        uuuuuu    uuuuuuttttttt:::::ttttttt       ooooooooooo          R::::R     R:::::R  E:::::E       EEEEEE           A:::::::::A       C:::::C       CCCCCCTTTTTT  T:::::T  TTTTTT  E:::::E       EEEEEE  R::::R     R:::::R
          A:::::A:::::A       u::::u    u::::ut:::::::::::::::::t     oo:::::::::::oo        R::::R     R:::::R  E:::::E                       A:::::A:::::A     C:::::C                      T:::::T          E:::::E               R::::R     R:::::R
         A:::::A A:::::A      u::::u    u::::ut:::::::::::::::::t    o:::::::::::::::o       R::::RRRRRR:::::R   E::::::EEEEEEEEEE            A:::::A A:::::A    C:::::C                      T:::::T          E::::::EEEEEEEEEE     R::::RRRRRR:::::R 
        A:::::A   A:::::A     u::::u    u::::utttttt:::::::tttttt    o:::::ooooo:::::o       R:::::::::::::RR    E:::::::::::::::E           A:::::A   A:::::A   C:::::C                      T:::::T          E:::::::::::::::E     R:::::::::::::RR  
       A:::::A     A:::::A    u::::u    u::::u      t:::::t          o::::o     o::::o       R::::RRRRRR:::::R   E:::::::::::::::E          A:::::A     A:::::A  C:::::C                      T:::::T          E:::::::::::::::E     R::::RRRRRR:::::R 
      A:::::AAAAAAAAA:::::A   u::::u    u::::u      t:::::t          o::::o     o::::o       R::::R     R:::::R  E::::::EEEEEEEEEE         A:::::AAAAAAAAA:::::A C:::::C                      T:::::T          E::::::EEEEEEEEEE     R::::R     R:::::R
     A:::::::::::::::::::::A  u::::u    u::::u      t:::::t          o::::o     o::::o       R::::R     R:::::R  E:::::E                  A:::::::::::::::::::::AC:::::C                      T:::::T          E:::::E               R::::R     R:::::R
    A:::::AAAAAAAAAAAAA:::::A u:::::uuuu:::::u      t:::::t    tttttto::::o     o::::o       R::::R     R:::::R  E:::::E       EEEEEE    A:::::AAAAAAAAAAAAA:::::AC:::::C       CCCCCC        T:::::T          E:::::E       EEEEEE  R::::R     R:::::R
   A:::::A             A:::::Au:::::::::::::::uu    t::::::tttt:::::to:::::ooooo:::::o     RR:::::R     R:::::REE::::::EEEEEEEE:::::E   A:::::A             A:::::AC:::::CCCCCCCC::::C      TT:::::::TT      EE::::::EEEEEEEE:::::ERR:::::R     R:::::R
  A:::::A               A:::::Au:::::::::::::::u    tt::::::::::::::to:::::::::::::::o     R::::::R     R:::::RE::::::::::::::::::::E  A:::::A               A:::::ACC:::::::::::::::C      T:::::::::T      E::::::::::::::::::::ER::::::R     R:::::R
 A:::::A                 A:::::Auu::::::::uu:::u      tt:::::::::::tt oo:::::::::::oo      R::::::R     R:::::RE::::::::::::::::::::E A:::::A                 A:::::A CCC::::::::::::C      T:::::::::T      E::::::::::::::::::::ER::::::R     R:::::R
AAAAAAA                   AAAAAAA uuuuuuuu  uuuu        ttttttttttt     ooooooooooo        RRRRRRRR     RRRRRRREEEEEEEEEEEEEEEEEEEEEEAAAAAAA                   AAAAAAA   CCCCCCCCCCCCC      TTTTTTTTTTT      EEEEEEEEEEEEEEEEEEEEEERRRRRRRR     RRRRRRR
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
