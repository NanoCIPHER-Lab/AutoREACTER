"""
LUNAR Utility Module

Provides cross-platform path resolution, file movement, and CLI feedback 
helpers for the LUNAR preparation workflow.
"""

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

def normalize_path(p: str | Path) -> str:
    """
    Convert paths between Windows and WSL formats as needed.
    
    Handles the conversion of Windows paths (C:/path) to WSL paths 
    (/mnt/c/path) and vice versa, ensuring compatibility with LUNAR 
    scripts regardless of the execution environment.
    
    Args:
        p: Input path string or Path object
        
    Returns:
        Normalized path string appropriate for the current platform
    """
    p = str(p).strip().strip('"').strip("'")
    p = p.replace("\\", "/")

    # Convert Windows paths to WSL format
    if is_wsl():
        m = re.match(r"^([A-Za-z]):/(.*)$", p)
        if m:
            drive = m.group(1).lower()
            rest = m.group(2)
            return f"/mnt/{drive}/{rest}"
        return os.path.normpath(p)

    host = platform.system().lower()

    # Convert WSL paths to Windows format
    if host == "windows":
        m = re.match(r"^/?mnt/([A-Za-z])/(.*)$", p)
        if m:
            drive = m.group(1).upper()
            rest = m.group(2).replace("/", "\\")
            return f"{drive}:\\{rest}"
        return os.path.normpath(p.replace("/", "\\"))

    return os.path.normpath(p)

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