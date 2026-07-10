MAX_LOOP = 5  # Maximum number of iterations for the reaction progression loop.

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from rdkit import Chem

from AutoREACTER.detectors.functional_groups_detector import FunctionalGroupsDetector
from AutoREACTER.detectors.reaction_detector import ReactionDetector

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
        reaction_instances = self.session.reaction_instances.copy()

        while iteration < max_loop:
            iteration += 1
            self.session.reaction_progression_session.iteration = iteration

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
            reaction_instances.extend(rxns)

            print(len(reaction_instances))  # Debug print
            print(monomer_roles_in_loop)  # Debug print
            prepared_reactions = self._index_based_reaction_preparation(
                reaction_instances=reaction_instances
            )
            print(prepared_reactions)  # Debug print
            if self._loop_break_condition(
                size_before=size_of_pool, size_after=len(monomer_roles_in_loop)
            ):
                return prepared_reactions
        return prepared_reactions

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
            product_mol = reaction.product_combined_RDmol

            monomer_roles_for_idx_based_fg_detection.append(
                MonomerRoleforIndexBasedFGDetection(
                    smiles=self._get_product_smiles(product_mol),
                    name=f"new_{reaction.reaction_id}",
                    indexes_in_template=self._get_product_index(
                        reaction.template_reactant_to_product_mapping
                    ),
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
        Removes atom map numbers and isotope labels from a product molecule.
        """
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
            atom.SetIsotope(0)

        return mol

    def _get_product_smiles(self, mol: Chem.Mol) -> str:
        """
        Converts an RDKit molecule object to its corresponding SMILES string.
        """
        self._clean_product(mol)

        try:
            return Chem.MolToSmiles(mol)
        except Exception:
            return ""

    def _get_product_index(
        self,
        template_reactant_to_product_mapping: dict,
    ) -> list[int]:
        """
        Retrieves product indexes from the reactant-to-product mapping.

        Args:
            template_reactant_to_product_mapping (dict): Mapping from reactant
                indices to product indices.

        Returns:
            list[int]: Corresponding product indexes.
        """
        product_indexes = []

        for product_idx in template_reactant_to_product_mapping.values():
            product_indexes.append(product_idx)

        print(f"Product indexes: {product_indexes}")  # Debug print

        return product_indexes

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