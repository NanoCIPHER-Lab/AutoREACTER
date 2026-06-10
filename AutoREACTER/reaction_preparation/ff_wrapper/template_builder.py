from pathlib import Path
from typing import List

from AutoREACTER.input_parser import SimulationSetup
from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import ReactionMetadata
from AutoREACTER.reaction_preparation.lunar_client.molecule_3d_preparation import Molecule3DPreparation
from AutoREACTER.reaction_preparation.lunar_client.lunar_api_wrapper import LunarAPIWrapper


class TemplateBuilder:
    """
    The TemplateBuilder class is responsible for constructing reaction templates based on the prepared reactions and monomer information. It takes the processed reaction data and generates the necessary files and structures for use in simulations.

    Attributes:
        cache_dir: Path to the directory where intermediate files and templates will be stored.
        simulation_setup: An instance of SimulationSetup containing the parsed input data and configuration for the simulation.
    """

    def __init__(self, cache_dir: Path, ):
        self.cache_dir = cache_dir


    def build_templates(self, simulation_setup: SimulationSetup, prepared_reactions: List[ReactionMetadata]) -> None:
        """
        Build reaction templates from the prepared reactions.

        Args:
            simulation_setup: An instance of SimulationSetup containing the parsed input data and configuration for the simulation.
            prepared_reactions: A list of ReactionMetadata objects containing information about each prepared reaction.

        Returns:
            None. The function generates template files in the specified cache directory.
        """
        molecule3dpreparation = Molecule3DPreparation(self.cache_dir)
        updated_inputs_with_3d_mols, prepared_reactions_with_3d_mols = molecule3dpreparation.prepare_molecule_3d_geometry(simulation_setup, prepared_reactions)
        lunar_api_wrapper = LunarAPIWrapper(self.cache_dir)
        lunar_api_wrapper.lunar_workflow(updated_inputs_with_3d_mols, prepared_reactions_with_3d_mols)
