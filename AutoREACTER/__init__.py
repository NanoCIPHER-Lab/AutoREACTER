"""
AutoREACTER

AutoREACTER is a tool for automated reaction-based molecular system generation.
"""
__version__ = "0.2.3"

__title__ = "AutoREACTER"
__author__ = "Janitha Mahanthe, Jacob Gissinger"

__package__ = "AutoREACTER"
__release__ = __version__
__license__ = "MIT"

__authors__ = [
    "Janitha Mahanthe",
    "Jacob Gissinger",
]

__author__ = ", ".join(__authors__)

# Import all public API symbols
from importlib.metadata import version, PackageNotFoundError
from pathlib import Path
from typing import Optional




"""
Public API for the AutoREACTER (ARX) pipeline.

This module exposes a **global-session** interface: call :func:`run` once to
bootstrap a workflow from an input file, then use the remaining functions to
inspect molecules, select reactions and non-reactants, and finally run the
full simulation pipeline.  Each function delegates to the currently active
:class:`ARXCLI` instance behind the scenes.

Usage (e.g. in a Jupyter notebook)::

    import arx

    arx.run("input.json")
    arx.show_molecules()
    arx.show_functional_groups()
    arx.show_reactions()
    arx.select_reactions()
    arx.show_non_reactants()
    arx.select_non_reactants()
    arx.process()

All public symbols are listed in :data:`__all__`.

:data:`__all__` includes:
    - ``run``
    - ``show_molecules``
    - ``show_functional_groups``
    - ``show_reactions``
    - ``select_reactions``
    - ``show_non_reactants``
    - ``select_non_reactants``
    - ``process``
"""

from .arx_cli import ARXCLI
from PIL import Image

# ---------------------------------------------------------------------------
# Global session handle – holds the *single* active workflow instance.
# The `run()` function is the only way to create (or replace) it.
# ---------------------------------------------------------------------------
_active_workflow = None  # type: ARXCLI | None


def _ensure_workflow() -> ARXCLI:
    """
    Return the active :class:`ARXCLI` instance, raising if none exists.

    This is the single choke-point used by every public function before
    delegating.  It guarantees that the user called :func:`run` first.

    Returns
    -------
    ARXCLI
        The currently active workflow session.

    Raises
    ------
    RuntimeError
        If :func:`run` has not been called yet (i.e. no session is active).
    """
    if _active_workflow is None:
        raise RuntimeError(
            "No active session. Please run `arx.run('your_file.json')` first."
        )
    return _active_workflow


# ===================================================================
# PUBLIC API
# ===================================================================

def run(input_file):
    """
    Create or replace the active AutoREACTER workflow session.

    Parameters
    ----------
    input_file : str or pathlib.Path
        Path to the AutoREACTER input JSON file.

    Returns
    -------
    ARXCLI
        Active AutoREACTER workflow object.
    """
    global _active_workflow

    input_file = Path(input_file).expanduser().resolve()

    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    _active_workflow = ARXCLI(input_file)
    return _active_workflow


def show_molecules() -> Image:
    """
    Display the 2D structures of all molecules in the current session.

    This triggers molecule visualisation (RDKit images) and is safe to call
    multiple times – the underlying detection runs only once.

    Returns
    -------
    Image
        Visualisation of all molecules, saved to the session's output directory.
    """
    return _ensure_workflow().show_molecules()


def show_functional_groups() -> Image:
    """
    Display the detected functional groups for each molecule.

    Functional-group detection must have completed; if it hasn't, this method
    triggers it automatically (idempotent).  The visualisation highlights the
    matched SMARTS patterns on each molecule.

    Returns
    -------
    Image
        Visualisation of functional groups, saved to the session's output directory.
    """
    return _ensure_workflow().show_functional_groups()


def show_reactions() -> Image:
    """
    Display the chemical reactions identified by the pipeline.

    If reaction detection has not run yet, it is triggered first.  The output
    includes reaction templates and their mappings onto the input molecules.

    Returns
    -------
    Image
        Visualisation of reactions, saved to the session's output directory.
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


def show_non_reactants() -> Image:
    """
    Display the non-reactant species detected in the simulation box.

    Non-reactant detection depends on reaction selection having already
    occurred; if it hasn't, an error will be raised.

    Returns
    -------
    Image
        Visualisation of non-reactant species, saved to the session's output directory.
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

def prepare_reactions() -> None:
    """
    Prepare the reaction templates for simulation.

    This is an intermediate step that performs reaction-template preparation
    without running the full pipeline.  It is not intended for end-users but
    may be useful for debugging or development.

    Returns
    -------
    None
        Writes prepared reaction templates to the session's output directory.
        Updates the internal :class:`ErrorHandler` state machine.
    """
    _ensure_workflow().prepare_reactions()


def show_reaction_templates(highlight_type: Optional[str] = "template") -> Image:
    """
    Visualise the prepared reaction templates.

    This is an intermediate step that visualises the reaction templates after
    preparation.  It is not intended for end-users but may be useful for
    debugging or development.

    Returns
    -------
    Image
        Visualisation of prepared reaction templates, saved to the session's output directory.
    """
    highlight_type = highlight_type.lower() if highlight_type else "template"
    highlight_types = [
        "template",    # Highlight the reaction template itself (i.e. the changed bonds).
        "edge",        # Highlight the edge atoms (dihedral distance away from the reaction center).
        "initiators",  # Highlight the initiator atoms/bonds that trigger the reaction.
        "delete"       # Highlight the atoms/bonds that are deleted in the reaction.
        ]
    if highlight_type not in highlight_types:
        raise ValueError(f"Invalid highlight_type: {highlight_type}. Must be one of {highlight_types}.")
    return _ensure_workflow().show_reaction_templates(highlight_type=highlight_type)


def process() -> None:
    """
    Execute the full simulation setup pipeline.

    This is the final step.  It performs, in order:

    1. 3D geometry generation
    2. Force-field assignment (via :class:`FFWrapper`)
    3. REACTER file creation
    4. LAMMPS input generation

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
    "__title__",
    "__version__",
    "__release__",
    "__authors__",
    "__license__",
    "run",
    "show_molecules",
    "show_functional_groups",
    "show_reactions",
    "select_reactions",
    "show_non_reactants",
    "select_non_reactants",
    "prepare_reactions",
    "show_reaction_templates",
    "process",
]
