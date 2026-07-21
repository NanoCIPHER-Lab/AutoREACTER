"""Aggregate and validate all polymer-family reaction modules."""

from .polyesters import REACTIONS as POLYESTERS
from .polyethers import REACTIONS as POLYETHERS
from .polyamides import REACTIONS as POLYAMIDES
from .polyanhydrides import REACTIONS as POLYANHYDRIDES
from .polythioesters import REACTIONS as POLYTHIOESTERS
from .polyurethanes import REACTIONS as POLYURETHANES
from .polyureas import REACTIONS as POLYUREAS
from .epoxy_polymers import REACTIONS as EPOXY_POLYMERS
from .vinyl_polymers import REACTIONS as VINYL_POLYMERS
from .polycarbonates import REACTIONS as POLYCARBONATES
from .polyimides import REACTIONS as POLYIMIDES
# from .polybenzimidazoles import REACTIONS as POLYBENZIMIDAZOLES
# from .phenolic_resins import REACTIONS as PHENOLIC_RESINS
from .polysiloxanes import REACTIONS as POLYSILOXANES
from .polysulfides import REACTIONS as POLYSULFIDES
from .thiol_ene_polymers import REACTIONS as THIOL_ENE_POLYMERS
# from .metathesis_polymers import REACTIONS as METATHESIS_POLYMERS
# from .cycloaddition_polymers import REACTIONS as CYCLOADDITION_POLYMERS

_REACTION_MODULES = [
    POLYESTERS,
    POLYETHERS,
    POLYAMIDES,
    POLYANHYDRIDES,
    POLYTHIOESTERS,
    POLYURETHANES,
    POLYUREAS,
    EPOXY_POLYMERS,
    VINYL_POLYMERS,
    POLYCARBONATES,
    POLYIMIDES,
    #POLYBENZIMIDAZOLES,
    # PHENOLIC_RESINS,
    POLYSILOXANES,
    POLYSULFIDES,
    THIOL_ENE_POLYMERS,
    # METATHESIS_POLYMERS,
    # CYCLOADDITION_POLYMERS,
]


def load_reactions() -> dict:
    """Return one flat reaction dictionary with duplicate-name protection."""
    merged = {}

    for module in _REACTION_MODULES:
        for reaction_name, reaction in module.items():
            if reaction_name in merged:
                raise ValueError(f"Duplicate reaction name: {reaction_name}")
            merged[reaction_name] = reaction

    return merged


REACTIONS = load_reactions()


class ReactionLibrary:
    """Backward-compatible class exposing ``self.reactions``."""

    def __init__(self):
        self.reactions = load_reactions()
