from pathlib import Path
from dataclasses import dataclass
from typing import Optional, TYPE_CHECKING
from AutoREACTER.input_parser import SimulationSetup
from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import ReactionMetadata

if TYPE_CHECKING:
    from AutoREACTER.session import ARXSession

# =============================================================================
# Shared Data Structures (Imported by both LUNAR and Foyer wrappers)
# =============================================================================

@dataclass(slots=True)
class DataFiles:
    data_file: Path
    lmp_molecule_file: Path

@dataclass(slots=True)
class MoleculeFile:
    id: str
    molecule_files: Optional[DataFiles]

@dataclass(slots=True)
class TemplateFile:
    reaction_id: Optional[int]
    pre_reaction_file: Optional[DataFiles]
    post_reaction_file: Optional[DataFiles]

@dataclass(slots=True)
class FFFiles:
    """The final compiled output from the chosen Force Field engine."""
    force_field_data: Path
    in_file: Optional[Path]  # LUNAR's generated input script (not applicable for Foyer)
    molecule_files: list[MoleculeFile]
    template_files: list[TemplateFile]


# =============================================================================
# The Central Wrapper Manager
# =============================================================================

class FFWrapper:
    """
    Central router that delegates molecule preparation to the correct 
    force field engine (LUNAR or Foyer) based on user input.
    """
    def __init__(self, ARX: "ARXSession"):
        self.session = ARX
        self.inputs = ARX.inputs

    def generate_force_field_files(
        self,
        updated_inputs: SimulationSetup,
        prepared_reactions: list[ReactionMetadata]
    ) -> FFFiles:
        
        force_field = updated_inputs.force_field.lower() if updated_inputs.force_field else ""

        # Route to FOYER
        if force_field in ["oplsaa", "opls", "opls-aa", "gaff"]: # this should be normalized in the input parser, not in this wrapper
            from AutoREACTER.reaction_preparation.ff_wrapper.foyer_client.foyer_api_wrapper import FoyerAPIWrapper
            print(f"Routing to Foyer for force field: {updated_inputs.force_field}")
            foyer_wrapper = FoyerAPIWrapper(ARX=self.session)
            return foyer_wrapper.run_workflow(updated_inputs, prepared_reactions)

        # Route to LUNAR
        elif force_field in ["pcff", "pcff-iff", "compass", "cvff", "cvff-iff", "drieding"]:
            from AutoREACTER.reaction_preparation.ff_wrapper.lunar_client.lunar_api_wrapper import LunarAPIWrapper
            print(f"Routing to LUNAR for force field: {updated_inputs.force_field}")
            lunar_wrapper = LunarAPIWrapper(ARX=self.session)
            return lunar_wrapper.lunar_workflow(updated_inputs, prepared_reactions)

        else:
            raise ValueError(f"Unsupported force field requested: {updated_inputs.force_field}")