"""Aggregate and validate all motif-based functional-group modules."""

from .oxygen_groups import FUNCTIONAL_GROUPS as OXYGEN_GROUPS
from .sulfur_groups import FUNCTIONAL_GROUPS as SULFUR_GROUPS
from .nitrogen_groups import FUNCTIONAL_GROUPS as NITROGEN_GROUPS
from .carboxyl_and_carbonyl_groups import FUNCTIONAL_GROUPS as CARBOXYL_AND_CARBONYL_GROUPS
from .mixed_ab_groups import FUNCTIONAL_GROUPS as MIXED_AB_GROUPS
from .ring_groups import FUNCTIONAL_GROUPS as RING_GROUPS
from .vinyl_and_alkene_groups import FUNCTIONAL_GROUPS as VINYL_AND_ALKENE_GROUPS
from .aromatic_groups import FUNCTIONAL_GROUPS as AROMATIC_GROUPS
from .silicon_groups import FUNCTIONAL_GROUPS as SILICON_GROUPS
# from .halide_groups import FUNCTIONAL_GROUPS as HALIDE_GROUPS
from .heterocumulene_groups import FUNCTIONAL_GROUPS as HETEROCUMULENE_GROUPS
from .active_centers import FUNCTIONAL_GROUPS as ACTIVE_CENTERS

_FUNCTIONAL_GROUP_MODULES = [
    OXYGEN_GROUPS,
    SULFUR_GROUPS,
    NITROGEN_GROUPS,
    CARBOXYL_AND_CARBONYL_GROUPS,
    MIXED_AB_GROUPS,
    RING_GROUPS,
    VINYL_AND_ALKENE_GROUPS,
    AROMATIC_GROUPS,
    SILICON_GROUPS,
    # HALIDE_GROUPS,
    HETEROCUMULENE_GROUPS,
    ACTIVE_CENTERS,
]


def load_functional_groups() -> dict:
    """Return one flat functional-group dictionary with duplicate protection."""
    merged = {}
    group_name_to_key = {}

    for module in _FUNCTIONAL_GROUP_MODULES:
        for entry_key, entry in module.items():
            if entry_key in merged:
                raise ValueError(f"Duplicate functional-group key: {entry_key}")

            group_name = entry["group_name"]
            if group_name in group_name_to_key:
                previous = group_name_to_key[group_name]
                raise ValueError(
                    f"Duplicate group_name {group_name!r} in {previous!r} and {entry_key!r}"
                )

            merged[entry_key] = entry
            group_name_to_key[group_name] = entry_key

    return merged


FUNCTIONAL_GROUPS = load_functional_groups()


class FunctionalGroupsLibrary:
    """Backward-compatible class exposing ``self.monomer_types``."""

    def __init__(self):
        self.monomer_types = load_functional_groups()
