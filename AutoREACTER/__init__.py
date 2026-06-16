from .arx_cli import ARXCLI

def _ensure_workflow():
    """Internal check to ensure the user loaded a file first."""
    if _active_workflow is None:
        raise RuntimeError("No active session. Please run `arx.load('your_file.json')` first.")
    return _active_workflow

# --- PUBLIC API ---

def load(input_file: str):
    global _active_workflow
    _active_workflow = ARXCLI(input_file)

def show_molecules():
    return _ensure_workflow().show_molecules()

def show_functional_groups():
    return _ensure_workflow().show_functional_groups()

def show_reactions():
    return _ensure_workflow().show_reactions()

def select_reactions():
    _ensure_workflow().select_reactions()

def show_non_reactants():
    return _ensure_workflow().show_non_reactants()

def select_non_reactants():
    _ensure_workflow().select_non_reactants()

def process():
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
    "process"
]