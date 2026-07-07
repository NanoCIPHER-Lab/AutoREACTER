
def ascii_art(message: str) -> None:
    message = message.upper()
    print(f"WARNING: {message}")
    """                  _             _ _ _ 
| |  | |                (_)           | | | |
| |  | | __ _ _ __ _ __  _ _ __   __ _| | | |
| |/\| |/ _` | '__| '_ \| | '_ \ / _` | | | |
\  /\  / (_| | |  | | | | | | | | (_| |_|_|_|
 \/  \/ \__,_|_|  |_| |_|_|_| |_|\__, (_|_|_)
                                  __/ |      
                                 |___/       
"""


def print_warning() -> None:
    message = "Warning " \
              "Entering to the reaction progression Loop still in the Beta phase" \
              "Caution: Can be chemically inaccurate"
    ascii_art(message)
    
