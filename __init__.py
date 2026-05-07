"""
AutoREACTER

Automated generation of LAMMPS/REACTER-ready reaction templates
for polymerization workflows.
"""

__version__ = "0.3.0b1"

from AutoREACTER.input_parser import (
    InputParser,
    SimulationSetup,
    MonomerEntry,
    Replica,
    InputError,
    InputSchemaError,
    InputConflictError,
    NumericFieldError,
    SmilesValidationError,
    DuplicateMonomerError,
)

from AutoREACTER.detectors.functional_groups_detector import (
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
