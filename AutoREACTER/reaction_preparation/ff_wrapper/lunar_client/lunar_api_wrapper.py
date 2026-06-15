"""
LUNAR API Wrapper Module

This module provides the high-level orchestration for the LUNAR 
preparation workflow. It delegates all heavy lifting to specialized 
submodules (locator, executor, builder).
"""

import os
from pathlib import Path
from typing import TYPE_CHECKING

from AutoREACTER.input_parser import SimulationSetup
from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import ReactionMetadata

# Import Central Wrapper Components
from AutoREACTER.reaction_preparation.ff_wrapper.ff_wrapper import FFFiles
from AutoREACTER.reaction_preparation.ff_wrapper.ff_locator import get_force_field_file

# Import Local LUNAR Components
from AutoREACTER.reaction_preparation.ff_wrapper.lunar_client.locate_lunar import get_LUNAR_loc
from AutoREACTER.reaction_preparation.ff_wrapper.lunar_client.lunar_executor import LunarExecutor
from AutoREACTER.reaction_preparation.ff_wrapper.lunar_client.merge_builder import write_bond_react_merge_input
from AutoREACTER.reaction_preparation.ff_wrapper.lunar_client.lunar_utils import loading_screen

if TYPE_CHECKING:
    from AutoREACTER.session import Session

class LunarAPIWrapper:
    """Orchestrates the multi-step LUNAR preparation workflow."""

    def __init__(self, ARX: "Session"):
        self.session = ARX
        self.inputs = ARX.inputs

        # Setup standard AutoREACTER cache directories
        self.cache_dir = Path(ARX.staging_dir) / "lunar"
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        # Locate LUNAR Installation
        self.LUNAR_LOCATION = Path(get_LUNAR_loc(use_gui=False))

        # Initialize the Executor (The Muscle)
        self.executor = LunarExecutor(
            lunar_location=self.LUNAR_LOCATION, 
            cache_dir=self.cache_dir
        )

    def lunar_workflow(
        self,
        updated_inputs_with_3d_mols: SimulationSetup,
        prepared_reactions_with_3d_mols: list[ReactionMetadata],
    ) -> FFFiles:
        """Executes the complete four-stage LUNAR preparation pipeline."""
        
        loading_screen("LUNAR Workflow")
        
        force_field_name = updated_inputs_with_3d_mols.force_field

        # Stage 0: Locate Force Field Parameters
        ff_file_path = get_force_field_file(
            force_field=force_field_name, 
            lunar_location=self.LUNAR_LOCATION
        )

        # Stage 1: Assign atom types 
        atom_typing_results = self.executor.run_atom_typing(
            updated_inputs=updated_inputs_with_3d_mols,
            prepared_reactions=prepared_reactions_with_3d_mols,
            force_field=force_field_name
        )

        # Stage 2: Convert to LAMMPS format
        all2lmp_results = self.executor.run_all2lmp(
            atom_typing_results=atom_typing_results,
            frc_file=ff_file_path
        )

        # Stage 3: Generate merge instruction file
        merge_input_file = write_bond_react_merge_input(
            cache_bond_react_merge=self.executor.cache_bond_react_merge,
            cache_all2lmp=self.executor.cache_all2lmp,
            all2lmp_results=all2lmp_results
        )

        # Stage 4: Merge molecules and templates into final files 
        final_lunar_files = self.executor.run_bond_react_merge(
            merge_input_file_path=merge_input_file,
            all2lmp_results=all2lmp_results
        )

        return final_lunar_files