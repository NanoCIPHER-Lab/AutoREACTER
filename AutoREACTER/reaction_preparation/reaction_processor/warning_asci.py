
def ascii_art(message: str) -> None:
    message = message.upper()
    print(f"WARNING: {message}")
    print(
        """                  _             _ _ _
| |  | |                (_)           | | | |
| |  | | __ _ _ __ _ __  _ _ __   __ _| | | |
| |/\| |/ _` | '__| '_ \| | '_ \ / _` | | | |
\  /\  / (_| | |  | | | | | | | | (_| |_|_|_|
 \/  \/ \__,_|_|  |_| |_|_|_| |_|\__, (_|_|_)
                                  __/ |
                                 |___/
"""
    )


def print_warning() -> None:
    message = (
        "Entering the reaction progression loop is still in the beta phase. "
        "Caution: results can be chemically inaccurate."
    )
    ascii_art(message)
    
