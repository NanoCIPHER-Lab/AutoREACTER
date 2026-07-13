MAX_LOOP = 5  # Maximum number of iterations for the reaction progression loop.

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from rdkit import Chem

from AutoREACTER.detectors.functional_groups_detector import FunctionalGroupsDetector
from AutoREACTER.detectors.reaction_detector import ReactionDetector
from AutoREACTER.reaction_preparation.deduplication_detector import DeduplicationDetector

if TYPE_CHECKING:
    from AutoREACTER.session import Session
    from AutoREACTER.detectors.functional_groups_detector import MonomerRole
    from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import ReactionInstance, ReactionMetadata


@dataclass(slots=True)
class MonomerRoleforIndexBasedFGDetection:
    """
    Represents a monomer role with its associated properties for
    index-based functional group detection.
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
    Placeholder for future reaction progression state.
    """
    monomer_roles: list["MonomerRole"] = field(default_factory=list)
    iteration: int = 0


class ReactionProgression:
    def __init__(self, session: "Session"):
        self.session = session
        self.session.reaction_progression_session = ReactionProgressionSession()
        self.fg_detector = FunctionalGroupsDetector()
        self.rxn_detector = ReactionDetector()

    def reaction_progression(self, max_loop: int = MAX_LOOP) -> list["ReactionMetadata"]:
        """
        Progresses the reaction by iteratively applying detected chemistries
        to the session's molecules.

        Args:
            max_loop (int): Maximum number of iterations for the reaction progression loop.
        """
        iteration = 0
        monomer_roles_in_loop = self.session.monomer_roles.copy()

        # Accumulate prepared reactions across iterations to avoid reprocessing in each loop.
        all_prepared_reactions: list["ReactionMetadata"] = []

        while iteration < max_loop:
            iteration += 1
            self.session.reaction_progression_session.iteration = iteration
            size_of_initial_reaction_pool = self._length_of_active_reactions()
            if iteration == 1:
                self._populate_monomer_roles()

            if iteration > 1:
                print(
                    f"Starting iteration {iteration} "
                    f"of the reaction progression loop."
                )

            size_of_pool = self._set_is_monomer_flag()

            print(self.session.monomer_roles)  # Debug print

            monomer_roles_for_idx_based_fg_detection = (
                self._prepare_products_for_idx_based_fg_detection()
            )

            fg_detection_results = (
                self.fg_detector.index_based_functional_groups_detector(
                    monomer_roles_for_idx_based_fg_detection
                )
            )

            if not fg_detection_results:
                print(
                    f"No new functional groups detected in iteration {iteration}. "
                    f"Ending the reaction progression loop."
                )
                break

            monomer_roles_in_loop.extend(fg_detection_results)

            rxns = self.rxn_detector.index_based_reaction_detector(
                monomer_roles_in_loop
            )

            # Debug prints to trace the reaction progression loop.
            if not rxns:
                print(
                    f"No new reactions detected in iteration {iteration}. "
                    f"Ending the reaction progression loop."
                )
                break

            print(len(rxns))  # Debug print (was len(reaction_instances); now just the new batch)
            print(monomer_roles_in_loop)  # Debug print

            # 
            prepared_reactions = self._index_based_reaction_preparation(
                reaction_instances=rxns
            )
            print(prepared_reactions)  # Debug print

            # Accumulate prepared reactions across iterations.
            all_prepared_reactions.extend(prepared_reactions)
            print(size_of_initial_reaction_pool)
            deduplication_detector = DeduplicationDetector()
            all_prepared_reactions = deduplication_detector.compare_graphs_mol(
                all_prepared_reactions
            )
            
            if self._loop_break_condition(
                size_before=size_of_initial_reaction_pool, size_after=len(all_prepared_reactions)
            ):
                self.session.reaction_metadata = all_prepared_reactions
                return all_prepared_reactions
        self.session.reaction_metadata = all_prepared_reactions
        return all_prepared_reactions

    def _index_based_reaction_preparation(
        self, reaction_instances: list["ReactionInstance"]
    ) -> list["ReactionMetadata"]:
        from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import PrepareReactions

        prepare_reactions = PrepareReactions(self.session)
        prepared_reactions = prepare_reactions._prepare_reactions_stage(reaction_instances)
        return prepared_reactions
    
    def _prepare_products_for_idx_based_fg_detection(
        self,
    ) -> list[MonomerRoleforIndexBasedFGDetection]:
        """
        Prepares generated reaction products for index-based functional group detection.
        """
        monomer_roles_for_idx_based_fg_detection = []
        reaction_metadata = self.session.reaction_metadata

        for reaction in reaction_metadata:
            if not reaction.activity_stats:
                continue
            product_mol = reaction.product_combined_RDmol
            indexes_in_template, product_mol = self._get_product_idxs(
                        reaction.template_reactant_to_product_mapping,
                        product_mol
                    )
            monomer_roles_for_idx_based_fg_detection.append(
                MonomerRoleforIndexBasedFGDetection(
                    smiles=self._get_product_smiles(product_mol),
                    name=f"new_{reaction.reaction_id}",
                    indexes_in_template=indexes_in_template,
                    is_monomer=False,
                    is_looped=False,
                    rdkit_mol=self._sanitize_molecule(product_mol),
                )
            )

        return monomer_roles_for_idx_based_fg_detection

    def _sanitize_molecule(self, mol: Chem.Mol) -> Chem.Mol | None:
        """
        Sanitizes an RDKit molecule object.

        Args:
            mol (Chem.Mol): RDKit molecule object.

        Returns:
            Chem.Mol | None: Sanitized molecule, or None if sanitization fails.
        """
        self._clean_product(mol)

        try:
            Chem.SanitizeMol(mol)
            return mol
        except Exception:
            return None

    def _clean_product(self, mol: Chem.Mol) -> Chem.Mol:
        """
        Return a copy of the molecule with atom-map numbers and isotope
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
        original RDKit molecule.
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
        Retrieve mapped product atom idxs and keep only the largest
        molecular fragment when the product contains multiple fragments.

        Returns:
            A tuple containing:
                - Product atom idxs relative to the returned molecule.
                - The complete product or its largest fragment.
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
        Keep the fragment with the largest number of heavy atoms and
        remap product atom idxs to the retained fragment.

        Args:
            mol:
                Molecule containing one or more disconnected fragments.
            product_idxs:
                Product atom idxs referring to the original molecule.

        Returns:
            A tuple containing:
                - The largest fragment.
                - Product atom idxs remapped to the largest fragment.
        """
        fragment_atom_mappings: list[tuple[int, ...]] = []

        fragments = Chem.GetMolFrags(
            mol,
            asMols=True,
            sanitizeFrags=True,
            fragsMolAtomMapping=fragment_atom_mappings,
        )

        if not fragments:
            raise ValueError("No fragments found in the product molecule.")

        largest_fragment_position = max(
            range(len(fragments)),
            key=lambda position: fragments[position].GetNumHeavyAtoms(),
        )

        largest_fragment = fragments[largest_fragment_position]

        # The atom mapping stores:
        # new fragment idx -> original molecule idx.
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
        Sets the is_looped flag for each monomer role in the session.

        Returns:
            int: Number of monomer roles before functional group detection.
        """
        for monomer_role in self.session.monomer_roles:
            monomer_role.is_looped = True

        return len(self.session.monomer_roles)

    def _populate_monomer_roles(self) -> None:
        """
        Populates RDKit molecule objects for monomers marked as monomers.
        """
        for monomer in self.session.monomer_roles:
            if monomer.is_monomer:
                monomer.rdkit_mol = self._smiles_to_rdkit_mol(monomer.smiles)

    def _smiles_to_rdkit_mol(self, smiles: str) -> Chem.Mol | None:
        """
        Converts a SMILES string to an RDKit molecule object.

        Args:
            smiles (str): SMILES string.

        Returns:
            Chem.Mol | None: RDKit molecule object.
        """
        return Chem.MolFromSmiles(smiles)

    def _loop_break_condition(self, size_before: int, size_after: int) -> bool:
        if size_after <= size_before:
            print(
                f"Breaking the loop as the pool did not grow "
                f"(before={size_before}, after={size_after})."
            )
            return True
        return False
    
    def _length_of_active_reactions(self) -> int:
        """
        Returns the number of reactions in the session that have activity stats.
        """
        return sum(
            1 for reaction in self.session.reaction_metadata
            if reaction.activity_stats
        )