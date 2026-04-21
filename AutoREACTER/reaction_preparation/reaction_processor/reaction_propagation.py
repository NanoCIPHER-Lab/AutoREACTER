from dataclasses import dataclass
from typing import List, Optional
from rdkit import Chem

from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import ReactionMetadata
from AutoREACTER.detectors.functional_groups_detector import (
    FunctionalGroupsDetector,
    MonomerRole,
)
from AutoREACTER.detectors.reaction_detector import (
    ReactionDetector,
    ReactionInstance,
)

MAX_ITERATIONS = 6  # Maximum iterations for reaction pooling to prevent infinite loops
@dataclass(slots=True)
class TemplateIndexedMolecule:
    mol: Chem.Mol
    indexes: List[int]
    name: Optional[str] = None

class ReactantPool:
    """
    A pool of reactants for possible reactions for multiple reaction loops.
    """
    loop_no : int
    pool : List[MonomerRole] 

class InfiniteReactionLoopError(Exception):
    """Raised when the reaction pooling process exceeds a reasonable number of iterations, indicating a potential infinite loop."""


class ReactionPropagation:
    def __init__(self):
        self.fg_detector = FunctionalGroupsDetector()
        self.rxn_detector = ReactionDetector()

    def run_propagation_loop(
        self,
        reactions_metadata: list[ReactionMetadata],
        monomer_roles: list[MonomerRole],
    ) -> list[ReactionInstance]:
        """
        Propagation loop:
        1. Convert prior reaction products into indexed molecules
        2. Detect functional groups only at relevant mapped indices
        3. Detect which polymerization reactions are now possible
        """
        for iteration in range(MAX_ITERATIONS + 1):
            if iteration >= MAX_ITERATIONS:
                raise InfiniteReactionLoopError(
                    f"Exceeded maximum iterations ({MAX_ITERATIONS}) in reaction propagation loop. "
                    "Potential infinite loop detected. "
                    "Try increasing MAX_ITERATIONS if this is a false positive. "
                    "If the issue persists, please raise an issue with the problematic reaction instances at: "
                    "https://github.com/NanoCIPHER-Lab/AutoREACTER/issues"
                )
            template_indexed_molecules = self._prepare_for_second_loop(reactions_metadata)

            if monomer_roles is None:
                raise ValueError("Monomer roles must be provided for the second loop of processing.")

            # Second-pass / index-based functional group detection
            detected_roles = self.fg_detector.index_based_functional_groups_detector(
                template_indexed_molecules
            )

            # If nothing reactive remains, propagation stops here
            if not detected_roles:
                return reactions_metadata

            # Detect reactions from the newly detected indexed roles
            reaction_instances = self.rxn_detector.index_based_reaction_detector(
                detected_roles
            )

            return reaction_instances

    def _prepare_for_second_loop(
        self,
        reaction_metadata_list: list[ReactionMetadata]
    ) -> list[TemplateIndexedMolecule]:
        """
        Prepare reaction products for propagation by extracting the post-reaction
        molecule and the mapped atom indices that should be used in the next loop.
        """
        template_indexed_molecules = []

        for reaction in reaction_metadata_list:
            if not reaction.activity_stats:
                continue  # Skip duplicates / inactive reactions

            if reaction.product_combined_RDmol:
                template_product_to_reactant_to_product_mapping = {
                    v: k for k, v in reaction.template_reactant_to_product_mapping.items()
                }

                template_indexed_molecules.append(
                    TemplateIndexedMolecule(
                        mol=reaction.product_combined_RDmol,
                        indexes=list(template_product_to_reactant_to_product_mapping.keys()),
                        name=f"reaction_{reaction.reaction_id}_product",
                    )
                )

        return template_indexed_molecules