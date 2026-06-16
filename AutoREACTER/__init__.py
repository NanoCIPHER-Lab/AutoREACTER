"""
Public API for the AutoREACTER (ARX) pipeline.

This module exposes a **global-session** interface: call :func:`load` once to
bootstrap a workflow from an input file, then use the remaining functions to
inspect molecules, select reactions and non-reactants, and finally run the
full simulation pipeline.  Each function delegates to the currently active
:class:`ARXCLI` instance behind the scenes.

Usage (e.g. in a Jupyter notebook)::

    import arx

    arx.load("input.json")
    arx.show_molecules()
    arx.show_functional_groups()
    arx.show_reactions()
    arx.select_reactions()
    arx.show_non_reactants()
    arx.select_non_reactants()
    arx.process()

All public symbols are listed in :data:`__all__`.

:data:`__all__` includes:
    - ``load``
    - ``show_molecules``
    - ``show_functional_groups``
    - ``show_reactions``
    - ``select_reactions``
    - ``show_non_reactants``
    - ``select_non_reactants``
    - ``process``
"""

from .arx_cli import ARXCLI

# ---------------------------------------------------------------------------
# Global session handle – holds the *single* active workflow instance.
# The `load()` function is the only way to create (or replace) it.
# ---------------------------------------------------------------------------
_active_workflow = None  # type: ARXCLI | None


def _ensure_workflow() -> ARXCLI:
    """
    Return the active :class:`ARXCLI` instance, raising if none exists.

    This is the single choke-point used by every public function before
    delegating.  It guarantees that the user called :func:`load` first.

    Returns
    -------
    ARXCLI
        The currently active workflow session.

    Raises
    ------
    RuntimeError
        If :func:`load` has not been called yet (i.e. no session is active).
    """
    if _active_workflow is None:
        raise RuntimeError(
            "No active session. Please run `arx.load('your_file.json')` first."
        )
    return _active_workflow


# ===================================================================
# PUBLIC API
# ===================================================================


def load(input_file: str) -> None:
    """
    Create (or replace) the global ARX workflow session.

    Reads the JSON/YAML input file, initialises internal state, performs
    functional-group detection, and sets up the working directories.  Once
    this returns successfully the other ``arx.*`` functions are ready to use.

    Parameters
    ----------
    input_file : str
        Path to the input configuration file (JSON or YAML).

    Returns
    -------
    None
        Sets the module-level ``_active_workflow`` as a side-effect.
    """
    global _active_workflow
    _active_workflow = ARXCLI(input_file)


def show_molecules() -> None:
    """
    Display the 2D structures of all molecules in the current session.

    This triggers molecule visualisation (RDKit images) and is safe to call
    multiple times – the underlying detection runs only once.

    Returns
    -------
    None
        Images are saved to the session's output directory.
    """
    return _ensure_workflow().show_molecules()


def show_functional_groups() -> None:
    """
    Display the detected functional groups for each molecule.

    Functional-group detection must have completed; if it hasn't, this method
    triggers it automatically (idempotent).  The visualisation highlights the
    matched SMARTS patterns on each molecule.

    Returns
    -------
    None
        Images are saved to the session's output directory.
    """
    return _ensure_workflow().show_functional_groups()


def show_reactions() -> None:
    """
    Display the chemical reactions identified by the pipeline.

    If reaction detection has not run yet, it is triggered first.  The output
    includes reaction templates and their mappings onto the input molecules.

    Returns
    -------
    None
        Reaction diagrams are saved to the session's output directory.
    """
    return _ensure_workflow().show_reactions()


def select_reactions() -> None:
    """
    Interactively select which reaction(s) to proceed with.

    If multiple candidate reactions were identified, the user is prompted to
    choose.  When only one reaction is found it is auto-selected.  This step
    must complete before non-reactants can be selected or :func:`process` can
    be called.

    Returns
    -------
    None
        Updates the internal :class:`ErrorHandler` state machine.
    """
    _ensure_workflow().select_reactions()


def show_non_reactants() -> None:
    """
    Display the non-reactant species detected in the simulation box.

    Non-reactant detection depends on reaction selection having already
    occurred; if it hasn't, an error will be raised.

    Returns
    -------
    None
        Visualisations are saved to the session's output directory.
    """
    return _ensure_workflow().show_non_reactants()


def select_non_reactants() -> None:
    """
    Interactively select which non-reactant species to include.

    This may prompt the user to pick from a list of detected species.  Must
    be called after :func:`select_reactions` and before :func:`process`.

    Returns
    -------
    None
        Updates the internal :class:`ErrorHandler` state machine.
    """
    _ensure_workflow().select_non_reactants()


def process() -> None:
    """
    Execute the full simulation setup pipeline.

    This is the final step.  It performs, in order:

    1. Reaction-template preparation
    2. 3D geometry generation
    3. Force-field assignment (via :class:`FFWrapper`)
    4. REACTER file creation
    5. LAMMPS input generation

    All preceding steps (reaction selection, non-reactant selection) must have
    completed first; otherwise a :class:`RuntimeError` is raised.

    Returns
    -------
    None
        Writes all simulation artifacts to the session's output directory.
    """
    _ensure_workflow().process()


# Expose only these commands to the user
__all__ = [
    "load",
    "show_molecules",
    "show_functional_groups",
    "show_reactions",
    "select_reactions",
    "show_non_reactants",
    "select_non_reactants",
    "process",
]
