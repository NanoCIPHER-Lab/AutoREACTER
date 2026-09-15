def __init__(self):
    """Load the RMG kinetics database."""

    try:
        from rmgpy import settings
        from rmgpy.data.kinetics.database import KineticsDatabase
        from rmgpy.molecule.molecule import Molecule

    except ImportError as exc:
        raise ImportError(
            "RMG is required for kinetics estimation. "
            "Please install RMG-Py and configure RMG-database "
            "before enabling kinetics."
        ) from exc

    self._rmg_molecule_class = Molecule
    self.database = KineticsDatabase()

    kinetics_path = os.path.join(
        settings["database.directory"],
        "kinetics",
    )

    if not os.path.isdir(kinetics_path):
        raise KineticsEstimationError(
            f"RMG kinetics database not found: {kinetics_path}"
        )

    self.database.load(
        kinetics_path,
        families="all",
        libraries=[],
    )