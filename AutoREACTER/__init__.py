"""
AutoREACTER

Automated generation of LAMMPS/REACTER-ready reaction templates
for polymerization workflows.
"""

__version__ = "0.2.2b0"

from .input_parser import (
    InputParser,
    SimulationSetup,
    MonomerEntry,
    Simulation,
    InputError,
    InputSchemaError,
    InputConflictError,
    NumericFieldError,
    SmilesValidationError,
    DuplicateMonomerError,
)

from .detectors.functional_groups_detector import (
    FunctionalGroupsDetector,
    FunctionalGroupInfo,
    MonomerRole,
    FunctionalGroupVisualization,
)

__all__ = [
    "__version__",
    "InputParser",
    "SimulationSetup",
    "MonomerEntry",
    "Replica",
    "InputError",
    "InputSchemaError",
    "InputConflictError",
    "NumericFieldError",
    "SmilesValidationError",
    "DuplicateMonomerError",
    "FunctionalGroupsDetector",
    "FunctionalGroupInfo",
    "MonomerRole",
    "FunctionalGroupVisualization",
]
