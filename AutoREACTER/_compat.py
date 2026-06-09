import sys
import collections
import collections.abc
import numpy as np

def apply_legacy_patches():
    """
    Injects removed aliases back into numpy and collections at runtime, 
    and bridges the OpenMM simtk namespace for older versions of Foyer/mBuild.
    """
    # 1. Restore removed Collections aliases
    _missing_classes = [
        "MutableSet", "MutableMapping", "Mapping", "MutableSequence",
        "Sequence", "Set", "Iterable", "Iterator", "Callable",
        "Container", "Hashable", "ItemsView", "KeysView", "ValuesView"
    ]
    for _name in _missing_classes:
        if not hasattr(collections, _name) and hasattr(collections.abc, _name):
            setattr(collections, _name, getattr(collections.abc, _name))

    # 2. Restore removed NumPy aliases
    if not hasattr(np, "float"): np.float = float
    if not hasattr(np, "int"): np.int = int
    if not hasattr(np, "complex"): np.complex = complex
    if not hasattr(np, "bool"): np.bool = np.bool_
    if not hasattr(np, "object"): np.object = np.object_
    if not hasattr(np, "str"): np.str = np.str_

    # 3. Robust OpenMM 'simtk' shim for Python 3.12+
    try:
        import sys
        import types
        import openmm
        import openmm.app
        import openmm.app.element
        import openmm.unit

        # Create pure, fake modules to bypass strict filesystem import checks
        simtk = types.ModuleType("simtk")
        simtk_openmm = types.ModuleType("simtk.openmm")
        simtk_openmm_app = types.ModuleType("simtk.openmm.app")

        # Copy the contents of the real modules into our fake ones
        simtk_openmm.__dict__.update(openmm.__dict__)
        simtk_openmm_app.__dict__.update(openmm.app.__dict__)
        
        # Manually wire the internal tree together
        simtk.openmm = simtk_openmm
        simtk.unit = openmm.unit
        simtk_openmm.app = simtk_openmm_app
        simtk_openmm_app.element = openmm.app.element

        # Register them in sys.modules so the 'import' statements find them instantly
        sys.modules["simtk"] = simtk
        sys.modules["simtk.openmm"] = simtk_openmm
        sys.modules["simtk.openmm.app"] = simtk_openmm_app
        sys.modules["simtk.openmm.app.element"] = openmm.app.element
        sys.modules["simtk.unit"] = openmm.unit
    except ImportError:
        pass