"""
Pytest configuration and fixtures for reaction_lammps_mupt tests.
"""
import sys
from pathlib import Path

# Add src directory to Python path
src_path = Path(__file__).parent.parent / "src"
sys.path.insert(0, str(src_path))