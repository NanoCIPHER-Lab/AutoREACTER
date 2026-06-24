from pathlib import Path
from dataclasses import dataclass
from typing import Optional, TYPE_CHECKING
from AutoREACTER.input_parser import SimulationSetup
from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import ReactionMetadata

if TYPE_CHECKING:
    from AutoREACTER.session import Session


# Shared data structures for molecule and template files across both LUNAR and FOYER outputs and future extensions. 
# This allows the FFWrapper to return a consistent output format regardless of the underlying engine used.

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
    
    molecule_files: list[MoleculeFile]
    template_files: list[TemplateFile]
    force_field_data: Path
    in_file: Optional[Path] = None


class FFWrapper:
    """
    Central router that delegates molecule preparation to the correct 
    force field engine (LUNAR or Foyer) based on user input.
    """
    def __init__(self, ARX: "Session"):
        self.session = ARX
        self.inputs = ARX.inputs

    def generate_force_field_files(
        self,
        session: "Session") -> None:
        """Generates force field files by routing to the appropriate engine based on user input.
        Args:
            session: The current Session containing validated inputs and reaction metadata.
        Returns:
            None
        """
        updated_inputs = session.inputs
        prepared_reactions = session.reaction_metadata
        if updated_inputs is None:
            updated_inputs = self.inputs
            
        if updated_inputs is None:
            raise ValueError("No inputs provided to FFWrapper and session inputs are None.")

        # Since input_parser already normalized this, we just grab the canonical name.
        # It should never be None at this stage if the parser did its job, but we fallback safely.
        ff_name = updated_inputs.force_field or "PCFF"

        # Define engine routing based on strict canonical names from the parser
        foyer_supported = {"OPLSAA", "GAFF"}
        lunar_supported = {
            "PCFF", "PCFF-IFF", "compass", 
            "CVFF", "CVFF-IFF", "DREIDING", "Clay-FF"
        }

        # Route to FOYER
        if ff_name in foyer_supported:
            from AutoREACTER.reaction_preparation.ff_wrapper.foyer_client.foyer_api_wrapper import FoyerAPIWrapper
            print(f"Routing to Foyer for force field: {ff_name}")
            foyer_wrapper = FoyerAPIWrapper(ARX=self.session, prepared_reactions_with_3d_mols=prepared_reactions)
            session.ff_files = foyer_wrapper.final_foyer_files
            return None

        # Route to LUNAR
        elif ff_name in lunar_supported:
            from AutoREACTER.reaction_preparation.ff_wrapper.lunar_client.lunar_api_wrapper import LunarAPIWrapper
            print(f"Routing to LUNAR for force field: {ff_name}")
            lunar_wrapper = LunarAPIWrapper(ARX=self.session)
            session.ff_files = lunar_wrapper.lunar_workflow(updated_inputs, prepared_reactions)
            return None
        else:
            raise ValueError(f"Unsupported or unrecognized force field requested: {ff_name}")