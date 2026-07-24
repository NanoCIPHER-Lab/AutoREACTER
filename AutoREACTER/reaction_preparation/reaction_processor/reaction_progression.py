"""
Iterative reaction-progression engine for AutoREACTER.

This module drives the discovery of follow-up reactions by repeatedly
re-detecting functional groups in products produced during earlier
reaction-generation rounds. Detected functional groups are turned into
new reaction instances, prepared into fully described reaction metadata,
deduplicated, and then fed back into the next loop iteration. A hard
iteration cap and a pool-growth guard prevent unbounded execution.
"""

from dataclasses import dataclass, field
import time
from typing import TYPE_CHECKING

from rdkit import Chem

from AutoREACTER.detectors.functional_groups_detector import (
    FunctionalGroupsDetector,
)
from AutoREACTER.detectors.reaction_detector import ReactionDetector
from AutoREACTER.reaction_preparation.deduplication_detector import (
    DeduplicationDetector,
)
from AutoREACTER.reaction_preparation.reaction_processor.warning_asci import (
    print_warning,
)

if TYPE_CHECKING:
    from AutoREACTER.detectors.functional_groups_detector import MonomerRole
    from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import (
        ReactionInstance,
        ReactionMetadata,
    )
    from AutoREACTER.session import Session


# Prevent progression from continuing indefinitely when newly generated
# products keep exposing additional detectable functional groups.
MAX_LOOP = 5


@dataclass(slots=True)
class MonomerRoleforIndexBasedFGDetection:
    """Describe a molecule prepared for index-based functional-group detection.

    These lightweight containers bridge the gap between prepared reaction
    products and the index-based functional-group detector. The stored atom
    indexes refer to positions in the parent reaction template, which makes
    it possible to trace any newly detected functional group back to the
    product (and ultimately the reactant) that produced it.

    Attributes:
        smiles: Canonical SMILES of the prepared product.
        name: Generated identifier for this potential monomer.
        indexes_in_template: Atom indexes within the reaction template.
        is_monomer: Whether this molecule should be treated as a monomer.
        is_looped: Whether this molecule has already been processed by a
            progression iteration.
        rdkit_mol: Optional RDKit molecule object for downstream use.
    """

    smiles: str
    name: str
    indexes_in_template: list[int]
    is_monomer: bool = False
    is_looped: bool = False
    rdkit_mol: Chem.Mol | None = None


@dataclass(slots=True)
class ReactionProgressionSession:
    """Track state shared across reaction-progression iterations.

    This object is attached to the main ``Session`` so that other stages of
    the pipeline can inspect or update progression-related bookkeeping.

    Attributes:
        monomer_roles: Roles that have been produced by progression and
            are eligible for subsequent functional-group detection.
        iteration: Current progression loop iteration (1-indexed).
    """

    monomer_roles: list["MonomerRole"] = field(default_factory=list)
    iteration: int = 0


class ReactionProgression:
    """Coordinate iterative functional-group detection and reaction generation.

    The progression workflow repeatedly broadens the set of considered
    monomers by using products from the previous round. In each iteration:

    1. Active products are sanitized and converted into monomer-like roles.
    2. Functional groups are re-detected in those products by index.
    3. Compatible reactions are detected from the expanded monomer pool.
    4. Reaction instances are prepared into full ``ReactionMetadata``.
    5. Radical centers are annotated before deduplication.
    6. Duplicate reaction products are merged.
    7. The loop continues if the active reaction pool has grown.

    The class also owns the cleanup and sanitization of RDKit product
    molecules, including special handling for radical carbons that are
    deliberately under-valent in the underlying reaction SMARTS.
    """

    def __init__(self, session: "Session", preparer=None):
        """Initialize detectors and attach progression state to a session.

        Args:
            session: The active AutoREACTER ``Session`` that contains the
                current monomer roles and reaction metadata.
            preparer: An existing ``PrepareReactions`` instance used to
                convert ``ReactionInstance`` objects into ``ReactionMetadata``.
        """
        self.session = session
        self.preparer = preparer

        # Attach a dedicated progression sub-session to the main session so
        # that loop state is visible elsewhere in the pipeline.
        self.session.reaction_progression_session = (
            ReactionProgressionSession()
        )

        self.fg_detector = FunctionalGroupsDetector()
        self.rxn_detector = ReactionDetector()
        self.deduplication_detector = DeduplicationDetector()

        # Display the runtime warning banner (e.g. about experimental
        # progression behavior or licensing).
        print_warning()

    def reaction_progression(
        self,
        max_loop: int = MAX_LOOP,
    ) -> list["ReactionMetadata"]:
        """Run progression until no further useful reactions are generated.

        Each iteration detects functional groups in generated products, finds
        compatible reactions, prepares them, and removes duplicates. The loop
        stops when no new functional groups or reactions are found, when the
        active reaction pool does not grow, or when ``max_loop`` is reached.

        Args:
            max_loop: Maximum number of progression iterations.

        Returns:
            Prepared and deduplicated reaction metadata accumulated across
            all iterations.
        """
        # Respect a user-supplied global loop limit if one was configured.
        if self.session.inputs.max_loop_count is not None:
            max_loop = self.session.inputs.max_loop_count
            print(
                f"Overriding default max_loop of {MAX_LOOP} with "
                f"user-specified max_loop_count of {max_loop}."
            )
            # Avoid sleeping in library code; callers control pacing.

        iteration = 0
        # Start from the monomer roles already present in the session.
        monomer_roles_in_loop = list(self.session.monomer_roles)
        # Accumulate reaction metadata across iterations for deduplication.
        all_prepared_reactions = list(self.session.reaction_metadata)

        while iteration < max_loop:
            iteration += 1
            self.session.reaction_progression_session.iteration = iteration

            if iteration == 1:
                # On the first pass, convert input monomer SMILES into
                # explicit RDKit molecules so detectors can work on them.
                self._populate_monomer_roles()
            else:
                print(
                    f"Starting iteration {iteration} "
                    "of the reaction progression loop."
                )

            # Mark roles from earlier iterations so detectors can distinguish
            # already processed molecules from newly added molecules.
            self._set_is_looped_flag(monomer_roles_in_loop)

            # Snapshot the pool size so we can decide whether this iteration
            # produced any genuinely new chemistry.
            initial_reaction_pool_size = self._count_active_reactions(
                self.session.reaction_metadata
            )
            print(
                f"Initial reaction pool size at iteration {iteration}: "
                f"{initial_reaction_pool_size}"
            )

            # Convert active products into the form required by the
            # index-based functional-group detector.
            roles_for_fg_detection = (
                self._prepare_products_for_idx_based_fg_detection()
            )
            fg_detection_results = (
                self.fg_detector.index_based_functional_groups_detector(
                    roles_for_fg_detection
                )
            )

            # No new functional groups means no new chemistry is possible.
            if not fg_detection_results:
                print(
                    f"No new functional groups detected in iteration "
                    f"{iteration}. Ending the reaction progression loop."
                )
                break

            # Expand the monomer pool with functional groups found in this
            # iteration's products.
            monomer_roles_in_loop.extend(fg_detection_results)
            self.session.monomer_roles = monomer_roles_in_loop

            # Search for reactions that involve the expanded monomer set.
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

            if not isinstance(reaction_instances, list):
                reaction_instances = list(reaction_instances)

            # Convert raw reaction instances into fully prepared metadata.
            prepared_reactions = self._index_based_reaction_preparation(
                reaction_instances=reaction_instances
            )

            # Radical identity must be available before NetworkX
            # deduplication, because equivalent products may differ only in
            # how radical centers are represented.
            self._annotate_radicals_before_deduplication(
                prepared_reactions
            )

            all_prepared_reactions.extend(prepared_reactions)
            self.session.reaction_metadata = all_prepared_reactions

            # Equivalent products can be generated through different reaction
            # paths, so deduplication occurs after preparation.
            all_prepared_reactions = (
                self.deduplication_detector.compare_graphs_mol(
                    all_prepared_reactions
                )
            )
            self.session.reaction_metadata = all_prepared_reactions

            deduplicated_reaction_count = self._count_active_reactions(
                all_prepared_reactions
            )

            # If the pool did not grow, further iterations are unlikely to
            # yield new chemistry.
            if self._loop_break_condition(
                size_before=initial_reaction_pool_size,
                size_after=deduplicated_reaction_count,
            ):
                return self._store_reactions(all_prepared_reactions)

        return all_prepared_reactions

    def _index_based_reaction_preparation(
        self,
        reaction_instances: list["ReactionInstance"],
    ) -> list["ReactionMetadata"]:
        """Convert detected reaction instances into reaction metadata.

        This is a thin wrapper around the preparer's loop-aware preparation
        stage, which builds full ``ReactionMetadata`` records (products,
        mappings, activity statistics, etc.) from the raw instances.

        Args:
            reaction_instances: Reaction instances produced by the detector.

        Returns:
            Fully prepared reaction metadata.
        """
        return self.preparer._prepare_reactions_stage(
            reaction_instances,
            loop=True,
        )

    def _prepare_products_for_idx_based_fg_detection(
        self,
    ) -> list[MonomerRoleforIndexBasedFGDetection]:
        """Prepare active products for index-based functional-group detection.

        Active reaction products are converted into cleaned SMILES strings
        and sanitized RDKit molecules. Their template atom indexes are
        retained so that functional groups detected in the product can be
        traced back to the reaction that produced them.

        Returns:
            A list of roles ready for index-based functional-group detection.
        """
        prepared_monomer_roles: list[
            MonomerRoleforIndexBasedFGDetection
        ] = []

        for reaction in self.session.reaction_metadata:
            # Skip reactions that were deactivated during preparation.
            if not reaction.activity_stats:
                continue

            product_mol = reaction.product_combined_RDmol
            indexes_in_template, product_mol = self._get_product_idxs(
                reaction.template_reactant_to_product_mapping,
                product_mol,
            )

            sanitized_mol, success = self._sanitize_molecule(product_mol)

            # Record radical metadata when sanitization succeeds; otherwise
            # mark the product as non-radical.
            if success and sanitized_mol is not None:
                self._set_reaction_radical_metadata(
                    reaction,
                    sanitized_mol,
                )
            else:
                reaction.is_radical = False
                reaction.radical_atom_idxs = ()

            if not success:
                print(
                    f"Skipping reaction product {reaction.reaction_id}: "
                    "RDKit molecule sanitization failed."
                )

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
        """Persist reaction metadata in the session and return it.

        Args:
            reactions: Final deduplicated reaction metadata.

        Returns:
            The same list, now stored on the session.
        """
        self.session.reaction_metadata = reactions
        return reactions

    def _sanitize_molecule(
        self,
        mol: Chem.Mol,
    ) -> tuple[Chem.Mol | None, bool]:
        """Clean and sanitize a product while preserving radical centers.

        RDKit ``RunReactants`` output is often not directly sanitizable,
        especially when the reaction SMARTS intentionally leaves a carbon
        under-valent (radical). This method strips tracking properties,
        recomputes ring information, and adjusts explicit hydrogens and
        radical electrons until the molecule either sanitizes cleanly or
        is returned in the best possible state.

        Args:
            mol: Raw product molecule from reaction preparation.

        Returns:
            A tuple of (sanitized molecule or best-effort molecule,
            success flag indicating whether RDKit sanitization succeeded).
        """
        cleaned_mol = self._clean_product(mol)
        patched_mol = Chem.RWMol(cleaned_mol)

        # Rebuild valence and ring state without raising on the first error.
        patched_mol.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(patched_mol)

        for atom in patched_mol.GetAtoms():
            # Focus on carbons, where radical centers are expected.
            if atom.GetAtomicNum() != 6:
                continue

            # Disable automatic implicit-H addition so we can manage valence
            # explicitly for the radical center.
            atom.SetNoImplicit(True)

            heavy_valence = int(
                sum(
                    bond.GetValenceContrib(atom)
                    for bond in atom.GetBonds()
                )
            )
            explicit_hs = atom.GetNumExplicitHs()
            radical_electrons = atom.GetNumRadicalElectrons()

            # A carbon that has reached valence four is no longer radical.
            if (
                heavy_valence + explicit_hs >= 4
                and radical_electrons > 0
            ):
                atom.SetNumRadicalElectrons(0)
                radical_electrons = 0

            # Reduce explicit hydrogens when reaction output is over-valent.
            if (
                heavy_valence
                + explicit_hs
                + radical_electrons
                > 4
            ):
                explicit_hs = max(
                    0,
                    4 - heavy_valence - radical_electrons,
                )
                atom.SetNumExplicitHs(explicit_hs)

            # Preserve a neutral, trivalent carbon as the new radical center.
            if (
                heavy_valence + explicit_hs == 3
                and atom.GetFormalCharge() == 0
            ):
                atom.SetNumRadicalElectrons(1)

        patched_mol = patched_mol.GetMol()
        patched_mol.ClearComputedProps()

        # First attempt: sanitize the valence-patched molecule.
        try:
            Chem.SanitizeMol(patched_mol)
            return patched_mol, True
        except Exception:
            pass

        # Second attempt: explicitly mark under-valent carbons as radicals
        # and try sanitization again.
        radical_fixed_mol = self._fix_radical_and_sanitize(
            patched_mol
        )

        try:
            Chem.SanitizeMol(radical_fixed_mol)
            return radical_fixed_mol, True
        except Exception:
            # If sanitization still fails, return the best-effort molecule
            # with updated caches so downstream code can still inspect it.
            radical_fixed_mol.UpdatePropertyCache(strict=False)
            Chem.FastFindRings(radical_fixed_mol)
            return radical_fixed_mol, False

    def _fix_radical_and_sanitize(
        self,
        raw_mol: Chem.Mol,
        query: str = "[CH;X3;v3]",
    ) -> Chem.Mol:
        """Represent deliberately under-valent carbons as radicals.

        ``RunReactants`` output can be unsanitized. The radical carbon is
        deliberately under-valent in the reaction SMARTS, so this method adds
        the radical electron needed for RDKit sanitization and SMILES
        round-tripping.

        Args:
            raw_mol: Molecule that may contain an under-valent radical carbon.
            query: SMARTS used to identify the radical carbon. Defaults to a
                neutral carbon with one hydrogen, three explicit connections,
                and total valence three.

        Returns:
            Molecule with radical valence represented explicitly.
        """
        mol = Chem.RWMol(raw_mol)
        mol.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(mol)

        query_mol = Chem.MolFromSmarts(query)
        hits = mol.GetSubstructMatches(query_mol)

        for match in hits:
            atom = mol.GetAtomWithIdx(match[0])
            atom.SetNoImplicit(True)

            # Ensure the matched carbon has exactly one explicit hydrogen.
            if atom.GetTotalNumHs() != 1:
                atom.SetNumExplicitHs(1)

            # Add the single radical electron that completes the valence
            # representation.
            atom.SetNumRadicalElectrons(1)

        return mol.GetMol()

    def _clean_product(self, mol: Chem.Mol) -> Chem.Mol:
        """Return a copy of a molecule stripped of preparation artifacts.

        Atom maps, isotope labels, and internal RDKit tracking properties are
        removed so that downstream SMILES are canonical and do not carry
        state from the reaction engine.

        Args:
            mol: Molecule to clean.

        Returns:
            A new molecule with atom maps, isotopes, and tracking properties
            cleared.
        """
        cleaned_mol = Chem.Mol(mol)

        for atom in cleaned_mol.GetAtoms():
            atom.SetAtomMapNum(0)
            atom.SetIsotope(0)

            if atom.HasProp("old_mapno"):
                atom.ClearProp("old_mapno")
            if atom.HasProp("react_atom_idx"):
                atom.ClearProp("react_atom_idx")

        return cleaned_mol

    def _get_product_smiles(self, mol: Chem.Mol) -> str:
        """Convert a cleaned product molecule to canonical SMILES.

        Args:
            mol: Product molecule (may be ``None`` if sanitization failed).

        Returns:
            Canonical SMILES string, or an empty string if conversion fails.
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

        When a product contains disconnected fragments, only the fragment with
        the greatest number of heavy atoms is retained and the indexes are
        remapped into that fragment.

        Args:
            template_reactant_to_product_mapping: Mapping from reactant atom
                indexes to product atom indexes in the original template.
            mol: Product molecule, possibly multi-fragment.

        Returns:
            A tuple of (list of remapped product indexes, retained fragment).
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
        """Keep the largest heavy-atom fragment and remap its atom indexes.

        Args:
            mol: Multi-fragment product molecule.
            product_idxs: Product-side atom indexes to retain.

        Returns:
            A tuple of (largest fragment molecule, product indexes remapped
            into that fragment).

        Raises:
            ValueError: If no fragments could be extracted from the molecule.
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

        # Select the fragment with the most heavy atoms.
        largest_fragment_position = max(
            range(len(fragments)),
            key=lambda position: (
                fragments[position].GetNumHeavyAtoms()
            ),
        )
        largest_fragment = fragments[largest_fragment_position]
        original_atom_idxs = fragment_atom_mappings[
            largest_fragment_position
        ]

        # Build a map from original atom indexes to their positions in the
        # largest fragment.
        original_to_new_idx = {
            original_idx: new_idx
            for new_idx, original_idx in enumerate(original_atom_idxs)
        }
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
        """Mark supplied monomer roles as processed by the current loop.

        Args:
            monomer_roles: Monomer roles to flag as looped.
        """
        for monomer_role in monomer_roles:
            monomer_role.is_looped = True

    def _populate_monomer_roles(self) -> None:
        """Create RDKit molecules for roles identified as monomers.

        Input monomers are typically supplied as SMILES; this method ensures
        that each one has an associated RDKit molecule before detection runs.
        """
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
            smiles: SMILES string to parse.

        Returns:
            The parsed RDKit molecule, or ``None`` if parsing fails.
        """
        return Chem.MolFromSmiles(smiles)

    def _loop_break_condition(
        self,
        size_before: int,
        size_after: int,
    ) -> bool:
        """Return whether the active reaction pool failed to grow.

        Args:
            size_before: Number of active reactions before deduplication.
            size_after: Number of active reactions after deduplication.

        Returns:
            ``True`` if the pool did not grow and the loop should stop.
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
        """Count reactions included in activity statistics.

        Args:
            reactions: Reaction metadata to inspect.

        Returns:
            Number of reactions with truthy ``activity_stats``.
        """
        return sum(
            bool(reaction.activity_stats)
            for reaction in reactions
        )

    def _set_reaction_radical_metadata(
        self,
        reaction: "ReactionMetadata",
        sanitized_product: Chem.Mol,
    ) -> None:
        """Store product radical atoms in reactant-index space.

        Deduplication relabels product atoms into reactant-index space, so
        radical indexes are converted through
        ``product_to_reactant_mapping``.

        Args:
            reaction: Reaction metadata to annotate.
            sanitized_product: Sanitized product molecule in which radical
                electrons have already been assigned.
        """
        product_radical_idxs = {
            atom.GetIdx()
            for atom in sanitized_product.GetAtoms()
            if atom.GetNumRadicalElectrons() > 0
        }
        radical_reactant_idxs = {
            reaction.product_to_reactant_mapping[product_idx]
            for product_idx in product_radical_idxs
            if product_idx in reaction.product_to_reactant_mapping
        }

        reaction.is_radical = bool(radical_reactant_idxs)
        reaction.radical_atom_idxs = tuple(
            sorted(radical_reactant_idxs)
        )

    def _annotate_radicals_before_deduplication(
        self,
        reactions: list["ReactionMetadata"],
    ) -> None:
        """Sanitize products and record radical atoms before deduplication.

        This must run before ``compare_graphs_mol`` because the deduplication
        step relies on consistent radical annotation to distinguish otherwise
        isomorphic products.

        Args:
            reactions: Newly prepared reaction metadata to annotate.
        """
        for reaction in reactions:
            if not reaction.activity_stats:
                continue

            product_mol = reaction.product_combined_RDmol

            if product_mol is None:
                reaction.is_radical = False
                reaction.radical_atom_idxs = ()
                continue

            sanitized_mol, success = self._sanitize_molecule(
                product_mol
            )

            if not success or sanitized_mol is None:
                reaction.is_radical = False
                reaction.radical_atom_idxs = ()
                continue

            self._set_reaction_radical_metadata(
                reaction,
                sanitized_mol,
            )
