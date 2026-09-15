from dataclasses import dataclass


@dataclass(slots=True)
class LammpsArrheniusConstraint:
    """LAMMPS-ready Arrhenius reaction constraint.

    These values must already be converted into the form expected by
    LAMMPS fix bond/react. They are intentionally kept separate from
    the raw RMG Arrhenius parameters.
    """

    A: float
    n: float
    Ea: float


@dataclass(slots=True)
class KineticsEstimate:
    """Kinetic information estimated by RMG for one simulation."""

    temperature: float
    pressure: float
    rate_coefficient: float

    kinetics_type: str
    rate_units: str | None = None

    # Raw RMG Arrhenius parameters in SI units.
    A: float | None = None
    n: float | None = None
    Ea: float | None = None
    T0: float | None = None

    # Populated only after a valid RMG -> LAMMPS conversion.
    lammps_arrhenius: LammpsArrheniusConstraint | None = None
