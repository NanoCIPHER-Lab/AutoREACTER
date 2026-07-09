MAX_LOOP = 5  # Maximum number of iterations for the reaction progression loop.

from typing import TYPE_CHECKING
from dataclasses import dataclass, field

from rdkit import Chem

from AutoREACTER.detectors.functional_groups_detector import FunctionalGroupsDetector

if TYPE_CHECKING:
    from AutoREACTER.session import Session
    from AutoREACTER.detectors.functional_groups_detector import MonomerRole


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
        self.reaction_progression()

    def reaction_progression(self, max_loop: int = MAX_LOOP) -> None:
        """
        Progresses the reaction by iteratively applying detected chemistries
        to the session's molecules.

        Args:
            max_loop (int): Maximum number of iterations for the reaction progression loop.
        """
        iteration = 0
        monomer_roles_in_loop = self.session.monomer_roles.copy()  # Create a copy to avoid modifying the original list during iteration

        while iteration < max_loop:
            iteration += 1
            
            if iteration == 1:
                self._populate_monomer_roles()

            if iteration > 1:
                print(
                    f"Starting iteration {iteration} "
                    f"of the reaction progression loop."
                )

            size_of_pool = self._set_is_monomer_flag()

            print(self.session.monomer_roles)  # Debug print
            monomer_roles_for_idx_based_fg_detection = self._prepare_products_for_idx_based_fg_detection()
            fg_detection_results = self.fg_detector.index_based_functional_groups_detector(
                monomer_roles_for_idx_based_fg_detection
                )
            monomer_roles_in_loop.extend(fg_detection_results)
            print(monomer_roles_in_loop)  # Debug print
            if self._loop_break_condition(size_of_pool):
                break

    def _prepare_products_for_idx_based_fg_detection(
        self,
    ) -> list[MonomerRoleforIndexBasedFGDetection]:
        """
        Prepares generated reaction products for index-based functional group detection.
        """
        monomer_roles_for_idx_based_fg_detection = []
        reaction_metadata = self.session.reaction_metadata

        for reaction in reaction_metadata:
            monomer_roles_for_idx_based_fg_detection.append(
                MonomerRoleforIndexBasedFGDetection(
                    smiles=self._get_product_smiles(reaction.product_combined_RDmol),
                    name=f"new_{reaction.reaction_id}",
                    indexes_in_template=self._get_product_index(
                        reaction.template_reactant_to_product_mapping
                        ),
                    is_monomer=False,
                    is_looped=False,
                    rdkit_mol=self._sanitize_molecule(reaction.product_combined_RDmol),

                )
            )

        return monomer_roles_for_idx_based_fg_detection
    
    def _sanitize_molecule(self, mol: Chem.Mol) -> Chem.Mol | None:
        """
        Sanitizes an RDKit molecule object to ensure it is chemically valid.

        Args:
            mol (Chem.Mol): The RDKit molecule object to sanitize.

        Returns:
            Chem.Mol | None: The sanitized RDKit molecule object, or None if sanitization fails.
        """
        self._clean_product(mol)  # Clean the product before sanitization
        try:
            Chem.SanitizeMol(mol)
            return mol
        except Exception:
            return None
        
    def _clean_product(self, products: Chem.Mol) -> Chem.Mol | None:
        for atom in products.GetAtoms():
            atom.SetAtomMapNum(0)
            atom.SetIsotope(0)

    def _get_product_smiles(self, mol: Chem.Mol) -> str:
        """
        Converts an RDKit molecule object to its corresponding SMILES string.
        """
        self._clean_product(mol)  # Clean the product before converting to SMILES
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
            smiles (str): The SMILES string to convert.

        Returns:
            Chem.Mol | None: The corresponding RDKit molecule object.
        """
        return Chem.MolFromSmiles(smiles)

    def _loop_break_condition(self, size_of_pool: int) -> bool:
        if size_of_pool <= len(self.session.monomer_roles):
            print(
                f"Breaking the loop as the size of the pool ({size_of_pool}) "
                f"is less than or equal to the number of monomer roles "
                f"({len(self.session.monomer_roles)})."
            )
            return True

        return False