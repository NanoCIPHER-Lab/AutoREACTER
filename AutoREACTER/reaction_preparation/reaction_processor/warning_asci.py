
import time


def ascii_art(message: str) -> None:
    message = message.upper()
    print(f"WARNING: {message}")
    print(
"""    
 ____      ____  _       _______     ____  _____  _____  ____  _____   ______    _   _   _   
|_  _|    |_  _|/ \     |_   __ \   |_   \|_   _||_   _||_   \|_   _|.' ___  |  | | | | | |  
  \ \  /\  / / / _ \      | |__) |    |   \ | |    | |    |   \ | | / .'   \_|  | | | | | |  
   \ \/  \/ / / ___ \     |  __ /     | |\ \| |    | |    | |\ \| | | |   ____  | | | | | |  
    \  /\  /_/ /   \ \_  _| |  \ \_  _| |_\   |_  _| |_  _| |_\   |_\ `.___]  | |_| |_| |_|  
     \/  \/|____| |____||____| |___||_____|\____||_____||_____|\____|`._____.'  (_) (_) (_)  
                                                                                                                                                                                                                      
"""
    )


def print_warning() -> None:
    message = (
        "Entering the reaction progression loop is still in the beta phase. "
        "Caution: results can be chemically inaccurate."
    )
    ascii_art(message)
    time.sleep(5)
    
