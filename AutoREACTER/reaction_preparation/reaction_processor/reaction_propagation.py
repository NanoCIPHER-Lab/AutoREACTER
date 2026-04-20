
from dataclasses import dataclass
from typing import List, Optional
from rdkit import Chem

from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import ReactionMetadata
from AutoREACTER.detectors.functional_groups_detector import FunctionalGroupInfo

@dataclass(slots=True)
class TemplateIndexedMolecule:
    mol: Chem.Mol
    indexes: List[int]
    name : Optional[str] = None
class ReactionPropagation:

    def run_propagation_loop(self, reactions_metadata: list[ReactionMetadata]) -> list[ReactionMetadata]:
        template_indexed_molecules = self._prepare_for_second_loop(reactions_metadata)


    def _prepare_for_second_loop(self, reaction_metadata_list: list[ReactionMetadata]) -> list[TemplateIndexedMolecule]:
        """
        JUST A PLACEHOLDER FOR NOW - THIS FUNCTION CAN BE EXPANDED TO PERFORM ANY NECESSARY PREPARATION STEPS BEFORE THE SECOND LOOP OF PROCESSING.
        Prepares reaction metadata for a second loop of processing by ensuring all necessary fields are populated.
        
        Args:
            reaction_metadata_list: List of ReactionMetadata objects to prepare

        Returns:
            List of TemplateIndexedMolecule objects containing molecules and their corresponding template indices
        """
        template_indexed_molecules = []
        
        for reaction in reaction_metadata_list:
            if not reaction.activity_stats:
                continue  # Skip reactions marked as duplicates
            
            if reaction.product_combined_RDmol:
                template_product_to_reactant_to_product_mapping = {v: k for k, v in reaction.template_reactant_to_product_mapping.items()}
                template_indexed_molecules.append(
                    TemplateIndexedMolecule(
                        mol=reaction.product_combined_RDmol,
                        indexes=list(template_product_to_reactant_to_product_mapping.keys()),
                        name=f"reaction_{reaction.reaction_id}_product"
                    )
                )
        
        return template_indexed_molecules
