from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from rdkit import Chem
from rdkit.Chem.rdchem import MolSanitizeException

from AutoREACTER.detectors.functional_groups_detector import (
    FunctionalGroupsDetector,
)
from AutoREACTER.detectors.reaction_detector import ReactionDetector
from AutoREACTER.reaction_preparation.deduplication_detector import (
    DeduplicationDetector,
)

if TYPE_CHECKING:
    from AutoREACTER.detectors.functional_groups_detector import MonomerRole
    from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import (
        ReactionInstance,
        ReactionMetadata,
    )
    from AutoREACTER.session import Session


# Prevent the progression process from continuing indefinitely when reactions
# keep producing additional detectable functional groups.
MAX_LOOP = 5


@dataclass(slots=True)
class MonomerRoleforIndexBasedFGDetection:
    """Describe a molecule prepared for index-based functional-group detection.

    The stored atom indexes refer to positions in the reaction template and
    allow functional groups to be associated with their original reactants.
    """

    smiles: str
    name: str
    indexes_in_template: list[int]
    is_monomer: bool = False
    is_looped: bool = False
    rdkit_mol: Chem.Mol | None = None


@dataclass(slots=True)
class ReactionProgressionSession:
    """Track state that is shared across reaction-progression iterations."""

    monomer_roles: list["MonomerRole"] = field(default_factory=list)
    iteration: int = 0


class ReactionProgression:
    """Coordinate iterative functional-group detection and reaction generation.

    Each iteration uses products from previously prepared reactions as potential
    monomers. Newly detected functional groups are converted into reaction
    instances, prepared into reaction metadata, and deduplicated before the
    next iteration begins.
    """

    def __init__(self, session: "Session"):
        """Initialize detectors and attach progression state to a session.

        Args:
            session: Session containing monomer roles and reaction metadata.
        """
        self.session = session
        self.session.reaction_progression_session = (
            ReactionProgressionSession()
        )

        self.fg_detector = FunctionalGroupsDetector()
        self.rxn_detector = ReactionDetector()
        self.deduplication_detector = DeduplicationDetector()

    def reaction_progression(
        self,
        max_loop: int = MAX_LOOP,
    ) -> list["ReactionMetadata"]:
        """Run the reaction-progression loop until no progress is possible.

        Each iteration detects functional groups in generated products, finds
        compatible reactions, prepares the reactions, and removes duplicates.
        The loop stops when no new functional groups or reactions are found,
        when the reaction pool does not grow, or when ``max_loop`` is reached.

        Args:
            max_loop: Maximum number of progression iterations to execute.

        Returns:
            The prepared and deduplicated reaction metadata.
        """
        iteration = 0
        monomer_roles_in_loop = list["MonomerRole"](
            self.session.monomer_roles
        )
        all_prepared_reactions: list["ReactionMetadata"] = list(
            self.session.reaction_metadata
        )

        while iteration < max_loop:
            iteration += 1
            self.session.reaction_progression_session.iteration = iteration

            if iteration == 1:
                # Convert the initial monomer SMILES strings into RDKit
                # molecules before the first functional-group search.
                self._populate_monomer_roles()
            else:
                print(
                    f"Starting iteration {iteration} "
                    "of the reaction progression loop."
                )

            # Roles seen in earlier iterations are marked so detectors can
            # distinguish already-processed molecules from newly added ones.
            self._set_is_looped_flag(monomer_roles_in_loop)

            initial_reaction_pool_size = self._count_active_reactions(
                self.session.reaction_metadata
            )

            print(
                f"Initial reaction pool size at iteration {iteration}: "
                f"{initial_reaction_pool_size}"
            )

            roles_for_fg_detection = (
                self._prepare_products_for_idx_based_fg_detection()
            )

            fg_detection_results = (
                self.fg_detector.index_based_functional_groups_detector(
                    roles_for_fg_detection
                )
            )

            if not fg_detection_results:
                print(
                    f"No new functional groups detected in iteration "
                    f"{iteration}. Ending the reaction progression loop."
                )
                break

            monomer_roles_in_loop.extend(fg_detection_results)
            self.session.monomer_roles = monomer_roles_in_loop

            reaction_instances = (
                self.rxn_detector.index_based_reaction_detector(
                    monomer_roles_in_loop
                )
            )

            if not reaction_instances:
                print(
                    f"No new reactions detected in iteration {iteration}. "
                    "Ending the reaction progression loop."
                )
                break

            prepared_reactions = self._index_based_reaction_preparation(
                reaction_instances=reaction_instances
            )

            all_prepared_reactions.extend(prepared_reactions)
            self.session.reaction_metadata = all_prepared_reactions

            # Deduplication occurs after preparation because equivalent
            # products may be generated through different reaction paths.
            all_prepared_reactions = (
                self.deduplication_detector.compare_graphs_mol(
                    all_prepared_reactions
                )
            )
            self.session.reaction_metadata = all_prepared_reactions

            deduplicated_reaction_count = self._count_active_reactions(
                all_prepared_reactions
            )

            should_break = self._loop_break_condition(
                size_before=initial_reaction_pool_size,
                size_after=deduplicated_reaction_count,
            )

            if should_break:
                return self._store_reactions(all_prepared_reactions)

        return all_prepared_reactions

    def _index_based_reaction_preparation(
        self,
        reaction_instances: list["ReactionInstance"],
    ) -> list["ReactionMetadata"]:
        """Convert detected reaction instances into prepared reaction metadata.

        The import is local to avoid importing the reaction-preparation module
        during module initialization, which also helps prevent circular imports.
        """
        from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import (
            PrepareReactions,
        )

        reaction_preparer = PrepareReactions(self.session)

        return reaction_preparer._prepare_reactions_stage(
            reaction_instances
        )

    def _prepare_products_for_idx_based_fg_detection(
        self,
    ) -> list[MonomerRoleforIndexBasedFGDetection]:
        """Prepare active reaction products for index-based FG detection.

        Products are converted into cleaned SMILES strings and RDKit molecules.
        Their template atom indexes are retained so newly detected functional
        groups can be traced back to the reaction that produced them.
        """
        prepared_monomer_roles: list[
            MonomerRoleforIndexBasedFGDetection
        ] = []

        for reaction in self.session.reaction_metadata:
            if not reaction.activity_stats:
                continue

            product_mol = reaction.product_combined_RDmol

            indexes_in_template, product_mol = self._get_product_idxs(
                reaction.template_reactant_to_product_mapping,
                product_mol,
            )

            sanitized_mol = self._sanitize_molecule(product_mol)

            if sanitized_mol is None:
                print(
                    f"Skipping reaction product {reaction.reaction_id}: "
                    "RDKit molecule sanitization failed."
                )
                continue

            prepared_monomer_roles.append(
                MonomerRoleforIndexBasedFGDetection(
                    smiles=self._get_product_smiles(sanitized_mol),
                    name=f"new_{reaction.reaction_id}",
                    indexes_in_template=indexes_in_template,
                    rdkit_mol=sanitized_mol,
                )
            )

        return prepared_monomer_roles

    def _store_reactions(
        self,
        reactions: list["ReactionMetadata"],
    ) -> list["ReactionMetadata"]:
        """Save reaction metadata to the session and return it."""
        self.session.reaction_metadata = reactions
        return reactions

    def _sanitize_molecule(
        self,
        mol: Chem.Mol,
    ) -> Chem.Mol | None:
        """Clean and sanitize an RDKit molecule.

        Atom-map numbers and isotope labels are removed before sanitization.
        Products that RDKit cannot sanitize are discarded and represented by
        ``None``.

        Args:
            mol: Molecule to clean and sanitize.

        Returns:
            The sanitized molecule, or ``None`` if sanitization fails.
        """
        cleaned_mol = self._clean_product(mol)

        try:
            Chem.SanitizeMol(cleaned_mol)
        except Exception as error:
            print(f"Could not sanitize product molecule: {error}")
            return None

        return cleaned_mol

    def _clean_product(self, mol: Chem.Mol) -> Chem.Mol:
        """Return a copy of ``mol`` without atom maps or isotope labels.

        The input molecule is copied and therefore remains unchanged.

        Args:
            mol: Molecule whose atom annotations should be removed.

        Returns:
            A cleaned copy of the input molecule.
        """
        cleaned_mol = Chem.Mol(mol)

        for atom in cleaned_mol.GetAtoms():
            atom.SetAtomMapNum(0)
            atom.SetIsotope(0)

        return cleaned_mol

    def _get_product_smiles(self, mol: Chem.Mol) -> str:
        """Convert a cleaned product molecule to canonical SMILES.

        Args:
            mol: Product molecule to serialize.

        Returns:
            The product SMILES, or an empty string if conversion fails.
        """
        cleaned_mol = self._clean_product(mol)

        try:
            return Chem.MolToSmiles(cleaned_mol)
        except Exception:
            return ""

    def _get_product_idxs(
        self,
        template_reactant_to_product_mapping: dict[int, int],
        mol: Chem.Mol,
    ) -> tuple[list[int], Chem.Mol]:
        """Return product indexes and the molecule containing those indexes.

        If the product contains multiple disconnected fragments, only the
        fragment with the greatest number of heavy atoms is retained. Indexes
        are remapped to match the retained fragment.

        Args:
            template_reactant_to_product_mapping:
                Mapping from template reactant atom indexes to product indexes.
            mol: Combined product molecule.

        Returns:
            A tuple containing remapped product indexes and the selected
            product molecule.
        """
        product = Chem.Mol(mol)
        product_idxs = list(
            template_reactant_to_product_mapping.values()
        )

        if len(Chem.GetMolFrags(product)) > 1:
            product, product_idxs = self._keep_largest_fragment(
                product,
                product_idxs,
            )

        return product_idxs, product

    def _keep_largest_fragment(
        self,
        mol: Chem.Mol,
        product_idxs: list[int],
    ) -> tuple[Chem.Mol, list[int]]:
        """Keep the largest disconnected fragment and remap atom indexes.

        Args:
            mol: Molecule containing one or more disconnected fragments.
            product_idxs: Atom indexes referring to the original molecule.

        Returns:
            The largest fragment and the indexes remapped to that fragment.

        Raises:
            ValueError: If the molecule contains no fragments.
        """
        fragment_atom_mappings: list[tuple[int, ...]] = []

        fragments = Chem.GetMolFrags(
            mol,
            asMols=True,
            sanitizeFrags=True,
            fragsMolAtomMapping=fragment_atom_mappings,
        )

        if not fragments:
            raise ValueError(
                "No fragments found in the product molecule."
            )

        largest_fragment_position = max(
            range(len(fragments)),
            key=lambda position: fragments[position].GetNumHeavyAtoms(),
        )

        largest_fragment = fragments[largest_fragment_position]

        # RDKit provides each retained fragment's original atom indexes.
        # Build the inverse mapping to translate indexes into fragment space.
        original_atom_idxs = fragment_atom_mappings[
            largest_fragment_position
        ]
        original_to_new_idx = {
            original_idx: new_idx
            for new_idx, original_idx in enumerate(original_atom_idxs)
        }

        # Ignore mapped indexes that belong to discarded fragments.
        remapped_product_idxs = [
            original_to_new_idx[product_idx]
            for product_idx in product_idxs
            if product_idx in original_to_new_idx
        ]

        return largest_fragment, remapped_product_idxs

    def _set_is_looped_flag(
        self,
        monomer_roles: list["MonomerRole"],
    ) -> None:
        """Mark supplied monomer roles as processed by the current loop."""
        for monomer_role in monomer_roles:
            monomer_role.is_looped = True

    def _populate_monomer_roles(self) -> None:
        """Create RDKit molecules for all roles identified as monomers."""
        for monomer in self.session.monomer_roles:
            if monomer.is_monomer:
                monomer.rdkit_mol = self._smiles_to_rdkit_mol(
                    monomer.smiles
                )

    def _smiles_to_rdkit_mol(
        self,
        smiles: str,
    ) -> Chem.Mol | None:
        """Parse a SMILES string into an RDKit molecule.

        Args:
            smiles: SMILES representation of the molecule.

        Returns:
            The parsed molecule, or ``None`` if RDKit cannot parse the string.
        """
        return Chem.MolFromSmiles(smiles)

    def _loop_break_condition(
        self,
        size_before: int,
        size_after: int,
    ) -> bool:
        """Return whether the active reaction pool failed to grow.

        A non-growing pool indicates that the latest iteration did not add
        useful reaction products and further progression is unlikely to help.

        Args:
            size_before: Active reaction count before the iteration.
            size_after: Active reaction count after deduplication.

        Returns:
            ``True`` when the pool stayed the same size or became smaller.
        """
        if size_after <= size_before:
            print(
                "Breaking the loop as the pool did not grow "
                f"(before={size_before}, after={size_after})."
            )
            return True

        return False

    def _count_active_reactions(
        self,
        reactions: list["ReactionMetadata"],
    ) -> int:
        """Count reactions that contain activity statistics.

        Args:
            reactions: Reaction metadata objects to inspect.

        Returns:
            The number of active reactions.
        """
        return sum(
            bool(reaction.activity_stats)
            for reaction in reactions
        )
