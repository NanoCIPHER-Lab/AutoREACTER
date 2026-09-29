"""
LUNAR Utility Module

Provides cross-platform path resolution, file movement, and CLI feedback 
helpers for the LUNAR preparation workflow.
"""

import ntpath
import posixpath
import os
import re
import sys
import time
import shutil
import platform
from pathlib import Path

def is_wsl() -> bool:
    """Check if running inside Windows Subsystem for Linux."""
    return ("microsoft" in platform.release().lower()) or ("WSL_INTEROP" in os.environ)

def normalize_path(path):
    """Normalize a user-supplied path for the current environment.

    POSIX/WSL behavior is implemented with posixpath so it is identical on
    every host OS; only a genuine native Windows run uses pathlib.
    """
    raw = str(path).strip().strip("\"'").strip() or "."

    # Native Windows (not WSL): use pathlib.
    if not is_wsl() and platform.system() == "Windows":
        # WSL-style path -> drive path: /mnt/c/x or mnt/c/x -> C:\\x
        m = re.match(r"^/?mnt/([A-Za-z])(?:/+(.*))?$", raw.replace("\\", "/"))
        if m:
            rest = (m.group(2) or "").replace("/", "\\")
            return ntpath.normpath(f"{m.group(1).upper()}:\\{rest}")
        return str(Path(raw).resolve())

    posix = raw.replace("\\", "/")

    # WSL: map "D:/x" -> "/mnt/d/x".
    m = re.match(r"^([A-Za-z]):/*(.*)$", posix)
    if m and is_wsl():
        posix = "/mnt/" + m.group(1).lower() + "/" + m.group(2)

    if not posixpath.isabs(posix):
        posix = posixpath.join(os.getcwd().replace("\\", "/"), posix)

    posix = posixpath.normpath(posix)

    # Resolve symlinks only on a real POSIX host and only if it exists.
    if os.name == "posix" and os.path.exists(posix):
        posix = os.path.realpath(posix)

    return posix


def get_ending_integer(s: str) -> int | None:
    """Extract trailing digits from a string (e.g., 'pre12' -> 12)."""
    match = re.search(r'\d+$', s)
    return int(match.group()) if match else None

def move_merge_outputs(src_dir: Path, dst_dir: Path):
    """
    Move generated files from all2lmp output to bond_react_merge directory.
    
    Transfers merged data files, molecule files, force field definitions,
    and log files between processing stages.
    """
    src_dir = Path(src_dir)
    dst_dir = Path(dst_dir)
    if not dst_dir.is_dir():
        dst_dir.mkdir(parents=True, exist_ok=True)

    # Move all relevant output files
    for pattern in ("*_merged.data", "*_merged.lmpmol", "force_field.data", "log.lammps", "*.log"):
        for f in src_dir.glob(pattern):
            shutil.move(str(f), str(dst_dir / f.name))

def loading_screen(name: str = "LUNAR") -> None:
    """
    Display an ASCII art banner with a loading spinner animation.
    
    Provides visual feedback during initialization of long-running 
    LUNAR operations.
    
    Args:
        name: The operation name to display in the loading message.
    """
    # ASCII art banner for LUNAR
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

    # Simple spinner animation
    animation = ["-", "\\", "|", "/"]
    for i in range(10):
        time.sleep(0.1)
        sys.stdout.write("\b" + animation[i % len(animation)])
        sys.stdout.flush()
    
    print("\nReady!")