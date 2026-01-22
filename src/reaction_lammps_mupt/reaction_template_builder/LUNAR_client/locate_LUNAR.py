r"""
Normal usage (prompts once, then saved in config.py)
from locate_LUNAR import get_LUNAR_loc

LUNAR_DIR = get_LUNAR_loc()

Force user to re-pick (reset)
from locate_LUNAR import reset_LUNAR_loc, get_LUNAR_loc

reset_LUNAR_loc()
LUNAR_DIR = get_LUNAR_loc(force_prompt=True)


Set programmatically (no prompt)
from locate_LUNAR import set_LUNAR_loc

set_LUNAR_loc(r"C:\Users\Janitha\lunar2")
"""

import tkinter as tk
from tkinter import filedialog, messagebox
from pathlib import Path
import os

try:
    import config
except ImportError:
    from lunar_client import config
# Global flag to determine whether to use GUI for input prompts. Defaults to False (CLI mode).
USE_GUI = False


def _normalize_path(p):
    """
    Normalize a path string by stripping whitespace, removing surrounding quotes,
    expanding user and environment variables, and converting Windows paths to WSL format
    if running in a WSL environment (detected by /mnt/c existence).

    Args:
        p (str): The input path string.

    Returns:
        str or None: The normalized and resolved path as a string, or None if input is empty.
    """
    if not p:
        return None

    p = p.strip()

    # Remove surrounding quotes if present.
    if (p.startswith('"') and p.endswith('"')) or (p.startswith("'") and p.endswith("'")):
        p = p[1:-1].strip()

    # Expand user home directory (~) and environment variables.
    p = os.path.expanduser(p)
    p = os.path.expandvars(p)

    # Convert Windows drive paths to WSL format if in a WSL environment.
    if p.startswith(("C:\\", "D:\\", "E:\\", "F:\\")) and os.path.exists("/mnt/c"):
        drive = p[0].lower()
        rest = p[2:].replace("\\", "/")
        p = f"/mnt/{drive}{rest}"

    # Resolve the path to an absolute path.
    return str(Path(p).resolve())


def _ask_gui():
    """
    Prompt the user to select a directory using a GUI file dialog (tkinter).

    Returns:
        str or None: The selected and normalized folder path, or None if canceled.
    """
    root = tk.Tk()
    root.withdraw()  # Hide the root window.
    folder = filedialog.askdirectory(title="Select LUNAR root folder")
    root.destroy()
    if not folder:
        return None
    return _normalize_path(folder)


def _ask_cli():
    """
    Prompt the user to enter a directory path via command-line input.

    Returns:
        str or None: The entered and normalized folder path, or None if empty (canceled).
    """
    folder = input("Enter LUNAR root folder path (Enter to cancel): ").strip()
    if not folder:
        return None
    return _normalize_path(folder)


def _is_valid_dir(p):
    """
    Check if the given path exists and is a valid directory.

    Args:
        p (str): The path to check.

    Returns:
        bool: True if the path is a valid directory, False otherwise.
    """
    if not p:
        return False
    return Path(p).is_dir()


def _write_config_py(lunar_path):
    """
    Write the LUNAR root directory path to the config.py file.

    Args:
        lunar_path (str or None): The path to write, or None to unset.
    """
    cfg_path = Path(config.__file__).resolve()
    # Use repr() to properly quote the string in the Python assignment.
    text = f"LUNAR_ROOT_DIR = {lunar_path!r}\n"
    cfg_path.write_text(text, encoding="utf-8")


def set_LUNAR_loc(folder_path):
    """
    Set the LUNAR root directory path programmatically without prompting the user.
    Validates the path and updates the config.py file.

    Args:
        folder_path (str): The path to the LUNAR root directory.

    Returns:
        str: The normalized and validated path.

    Raises:
        ValueError: If the provided path is not a valid directory.
    """
    folder_path = _normalize_path(folder_path)
    if not _is_valid_dir(folder_path):
        raise ValueError(f"Not a valid directory: {folder_path}")
    config.LUNAR_ROOT_DIR = folder_path
    _write_config_py(folder_path)
    return folder_path


def reset_LUNAR_loc():
    """
    Reset the LUNAR root directory to None, clearing the saved path in config.py.

    Returns:
        None.
    """
    config.LUNAR_ROOT_DIR = None
    _write_config_py(None)
    return None


def get_LUNAR_loc(force_prompt=False, use_gui=None):
    """
    Get the LUNAR root directory path. If a valid path is already saved in config,
    return it. Otherwise, prompt the user to select or enter a path (looping until valid).

    Args:
        force_prompt (bool): If True, ignore saved path and always prompt. Defaults to False.
        use_gui (bool or None): If True, use GUI for prompting; False for CLI. If None,
                                defaults to the global USE_GUI value.

    Returns:
        str or None: The valid LUNAR root path, or None if the user cancels.
    """
    if use_gui is None:
        use_gui = USE_GUI

    # Use saved path if available, valid, and not forcing a prompt.
    if not force_prompt and config.LUNAR_ROOT_DIR and _is_valid_dir(config.LUNAR_ROOT_DIR):
        print("Using saved LUNAR root directory:", config.LUNAR_ROOT_DIR)
        return config.LUNAR_ROOT_DIR

    # Loop until a valid path is provided or canceled.
    while True:
        new_path = _ask_gui() if use_gui else _ask_cli()
        if not new_path:
            return None

        if _is_valid_dir(new_path):
            config.LUNAR_ROOT_DIR = new_path
            _write_config_py(new_path)
            return new_path

        # Handle invalid path feedback.
        if use_gui:
            root_msg = tk.Tk()
            root_msg.withdraw()
            messagebox.showerror("Invalid Folder", f"Not a valid directory:\n\n{new_path}")
            root_msg.destroy()
        else:
            print("Invalid directory. Try again.")


if __name__ == "__main__":
    # Entry point: Prompt for LUNAR location if run as a script.
    loc = get_LUNAR_loc(force_prompt=False, use_gui=True)
