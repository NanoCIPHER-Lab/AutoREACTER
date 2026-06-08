"""
Force Field Locator Module

Top-level utility for finding and validating force field parameter files 
(e.g., .frc for LUNAR, or .xml for Foyer). 
"""

from pathlib import Path
import importlib.resources as pkg_resources
from typing import Optional

def get_force_field_file(force_field: str, lunar_location: Optional[Path] = None) -> Path:
    """
    Map a force field name to its corresponding parameter file.
    
    Resolves paths for internal AutoREACTER force fields (in FF_files) 
    as well as external engine-specific files (like LUNAR's install directory).
    
    Args:
        force_field: Name of the force field (e.g., "PCFF-IFF", "CVFF", "OPLSAA")
        lunar_location: Path to the LUNAR installation directory (required for some LUNAR FFs)
        
    Returns:
        Absolute Path object to the verified parameter file.
        
    Raises:
        ValueError: If the force field file is not found or unsupported.
        RuntimeError: If there is an issue with pkg_resources.
        NotImplementedError: If the force field is planned but not yet implemented.
    """

    # 1. LUNAR Force Fields (.frc files)
    if force_field in {"PCFF", "PCFF-IFF", "Compass", "CVFF", "CVFF-IFF", "DRIEDING"}:
        paths = {}
        
        # Load internal PCFF via pkg_resources
        try:
            internal_ff_dir = Path(
                pkg_resources.files("AutoREACTER").joinpath("reaction_preparation", "ff_wrapper", "FF_files")
            )
            pcff_frc = internal_ff_dir / "pcff.frc"
            paths["PCFF-IFF"] = pcff_frc
            paths["PCFF"] = pcff_frc
        except Exception as e:
            raise RuntimeError(f"Error locating FF_files directory using pkg_resources: {e}")

        # Add external LUNAR-specific files
        if lunar_location:
            lunar_base = Path(lunar_location) / "frc_files"
            paths.update({
                "Compass": lunar_base / "compass_published.frc",
                "CVFF-IFF": lunar_base / "cvff_aug.frc",
                "CVFF": lunar_base / "cvff.frc",
                "DRIEDING": lunar_base / "all2lmp_dreiding.frc",
            })

        if force_field not in paths:
            raise ValueError(f"Required LUNAR installation not found to resolve: {force_field}")

        ff_path = paths[force_field]
        
        if ff_path.is_file():
            return ff_path
            
        raise ValueError(f"Force field file not found at {ff_path}. Check LUNAR installation.")

    # 2. Foyer Force Fields (.xml files) - PLACEHOLDER
    elif force_field in {"OPLSAA", "GAFF"}:
        # Keep this as a working placeholder for future implementation
        raise NotImplementedError(f"Foyer force field '{force_field}' support is currently in development.")

    else:
        raise ValueError(f"Unsupported force field requested: '{force_field}'.")