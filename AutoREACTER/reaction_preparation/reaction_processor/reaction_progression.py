MAX_LOOP = 5  # Maximum number of reaction progression iterations.

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from rdkit import Chem

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


@dataclass(slots=True)
class MonomerRoleforIndexBasedFGDetection:
    """
    Represents a monomer role used for index-based functional-group
    detection.
    """

    smiles: str
    name: str
    indexes_in_template: list[int]
    is_monomer: bool = False
    is_looped: bool = False
    rdkit_mol: Chem.Mol | None = None


@dataclass(slots=True)
class ReactionProgressionSession:
    """
    Stores state associated with reaction progression.
    """

    monomer_roles: list["MonomerRole"] = field(default_factory=list)
    iteration: int = 0


class ReactionProgression:
    def __init__(self, session: "Session"):
        self.session = session
        self.session.reaction_progression_session = (
            ReactionProgressionSession()
        )

        self.fg_detector = FunctionalGroupsDetector()
        self.rxn_detector = ReactionDetector()

    def reaction_progression(
        self,
        max_loop: int = MAX_LOOP,
    ) -> list["ReactionMetadata"]:
        """
        Progress a reaction by repeatedly detecting functional groups,
        detecting reactions, preparing reactions, and removing duplicates.

        Args:
            max_loop:
                Maximum number of progression iterations.

        Returns:
            Prepared reaction metadata accumulated across the progression
            loop.
        """
        iteration = 0
        monomer_roles_in_loop = self.session.monomer_roles.copy()
        all_prepared_reactions: list["ReactionMetadata"] = []

        while iteration < max_loop:
            iteration += 1
            self.session.reaction_progression_session.iteration = iteration

            size_of_initial_reaction_pool = (
                self._length_of_active_reactions()
            )

            if iteration == 1:
                self._populate_monomer_roles()

            if iteration > 1:
                print(
                    f"Starting iteration {iteration} "
                    "of the reaction progression loop."
                )

            self._set_is_monomer_flag()

            initial_reaction_pool_size = (
                self._length_of_active_reactions()
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

            print(len(reaction_instances))

            prepared_reactions = (
                self._index_based_reaction_preparation(
                    reaction_instances=reaction_instances
                )
            )

            print(prepared_reactions)

            all_prepared_reactions.extend(prepared_reactions)

            print(
                f"Total prepared reactions after iteration {iteration}: "
                f"{len(all_prepared_reactions)}"
            )
            print(size_of_initial_reaction_pool)

            deduplication_detector = DeduplicationDetector()

            all_prepared_reactions = (
                deduplication_detector.compare_graphs_mol(
                    all_prepared_reactions
                )
            )

            should_break = self._loop_break_condition(
                size_before=size_of_initial_reaction_pool,
                size_after=len(all_prepared_reactions),
            )

            if should_break:
                self.session.reaction_metadata = all_prepared_reactions
                return all_prepared_reactions

        self.session.reaction_metadata = all_prepared_reactions
        return all_prepared_reactions

    def _index_based_reaction_preparation(
        self,
        reaction_instances: list["ReactionInstance"],
    ) -> list["ReactionMetadata"]:
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
        """
        Prepare generated products for index-based functional-group
        detection.
        """
        prepared_monomer_roles = []
        reaction_metadata = self.session.reaction_metadata

        for reaction in reaction_metadata:
            if not reaction.activity_stats:
                continue

            product_mol = reaction.product_combined_RDmol

            indexes_in_template, product_mol = self._get_product_idxs(
                reaction.template_reactant_to_product_mapping,
                product_mol,
            )

            prepared_monomer_roles.append(
                MonomerRoleforIndexBasedFGDetection(
                    smiles=self._get_product_smiles(product_mol),
                    name=f"new_{reaction.reaction_id}",
                    indexes_in_template=indexes_in_template,
                    is_monomer=False,
                    is_looped=False,
                    rdkit_mol=self._sanitize_molecule(product_mol),
                )
            )

        return prepared_monomer_roles

    def _sanitize_molecule(
        self,
        mol: Chem.Mol,
    ) -> Chem.Mol | None:
        """
        Sanitize an RDKit molecule.

        Args:
            mol:
                RDKit molecule to sanitize.

        Returns:
            The sanitized molecule, or ``None`` when sanitization fails.
        """
        self._clean_product(mol)

        try:
            Chem.SanitizeMol(mol)
            return mol
        except Exception:
            return None

    def _clean_product(self, mol: Chem.Mol) -> Chem.Mol:
        """
        Return a copy of a molecule with atom-map numbers and isotope
        labels removed.

        The input molecule is not modified.
        """
        cleaned_mol = Chem.Mol(mol)

        for atom in cleaned_mol.GetAtoms():
            atom.SetAtomMapNum(0)
            atom.SetIsotope(0)

        return cleaned_mol

    def _get_product_smiles(self, mol: Chem.Mol) -> str:
        """
        Convert a product molecule to SMILES without modifying the
        original molecule.
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
        """
        Retrieve mapped product atom indexes.

        When the product contains disconnected fragments, retain only
        the fragment with the largest number of heavy atoms and remap
        the product indexes to that fragment.

        Returns:
            A tuple containing the product atom indexes and the retained
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

        print(f"Product idxs: {product_idxs}")

        return product_idxs, product

    def _keep_largest_fragment(
        self,
        mol: Chem.Mol,
        product_idxs: list[int],
    ) -> tuple[Chem.Mol, list[int]]:
        """
        Retain the fragment with the largest number of heavy atoms and
        remap product atom indexes to the retained fragment.

        Args:
            mol:
                Molecule containing one or more disconnected fragments.
            product_idxs:
                Product atom indexes referring to the original molecule.

        Returns:
            A tuple containing the largest fragment and the remapped
            product atom indexes.
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
            key=lambda position: (
                fragments[position].GetNumHeavyAtoms()
            ),
        )

        largest_fragment = fragments[largest_fragment_position]

        # Mapping direction:
        # fragment atom index -> original molecule atom index
        original_atom_idxs = fragment_atom_mappings[
            largest_fragment_position
        ]

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

    def _set_is_monomer_flag(self) -> int:
        """
        Mark every session monomer role as looped.

        Returns:
            Number of monomer roles in the session.
        """
        for monomer_role in self.session.monomer_roles:
            monomer_role.is_looped = True

        return len(self.session.monomer_roles)

    def _populate_monomer_roles(self) -> None:
        """
        Populate RDKit molecules for roles marked as monomers.
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
        """
        Convert a SMILES string to an RDKit molecule.
        """
        return Chem.MolFromSmiles(smiles)

    def _loop_break_condition(
        self,
        size_before: int,
        size_after: int,
    ) -> bool:
        """
        Determine whether progression should stop based on pool growth.
        """
        if size_after <= size_before:
            print(
                "Breaking the loop as the pool did not grow "
                f"(before={size_before}, after={size_after})."
            )
            return True

        return False

    def _length_of_active_reactions(self) -> int:
        """
        Return the number of session reactions with activity statistics.
        """
        return sum(
            1
            for reaction in self.session.reaction_metadata
            if reaction.activity_stats
        )