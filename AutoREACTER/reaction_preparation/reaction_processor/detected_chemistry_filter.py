from __future__ import annotations  # 1. Must be the first line

from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from AutoREACTER.session import Session


@dataclass(slots=True)
class DetectedChemistries:
    functional_groups: list[str]  # Unique functional group SMARTS detected in available reactions
    reactions: dict[str, str]  # reaction_name -> reaction_smarts


class DetectedChemistryFilter:
    """
    Collects available detected chemistry from reaction instances.

    This class extracts:
    - unique functional group SMARTS used by detected reactions
    - unique reaction names and their reaction SMARTS
    """

    def __init__(self, session: Session):
        self.reaction_instances = session.reaction_instances or []

    def _add_to_the_list(self, list_to_add: list[str], item: Optional[str]) -> None:
        """
        Adds an item to the list if it is not None and not already present.
        """
        if item is not None and item not in list_to_add:
            list_to_add.append(item)

    def _add_to_the_dict(self, dict_to_add: dict[str, str], key: str, value: str) -> None:
        """
        Adds a key-value pair to the dictionary if the key is not already present.
        """
        if key not in dict_to_add:
            dict_to_add[key] = value

    def filter(self) -> DetectedChemistries:
        """
        Extracts available functional group SMARTS and reaction SMARTS
        from detected reaction instances.

        Returns:
            DetectedChemistries: Available functional groups and reactions.
        """
        available_functional_groups: list[str] = []
        available_reactions: dict[str, str] = {}

        for reaction_instance in self.reaction_instances:
            reaction_name = reaction_instance.reaction_name
            reaction_smarts = reaction_instance.reaction_smarts

            self._add_to_the_dict(
                available_reactions,
                reaction_name,
                reaction_smarts,
            )

            functional_group_1 = reaction_instance.functional_group_1
            self._add_to_the_list(
                available_functional_groups,
                functional_group_1.fg_smarts_1,
            )
            self._add_to_the_list(
                available_functional_groups,
                functional_group_1.fg_smarts_2,
            )

            functional_group_2 = reaction_instance.functional_group_2
            if functional_group_2 is not None:
                self._add_to_the_list(
                    available_functional_groups,
                    functional_group_2.fg_smarts_1,
                )
                self._add_to_the_list(
                    available_functional_groups,
                    functional_group_2.fg_smarts_2,
                )

        return DetectedChemistries(
            functional_groups=available_functional_groups,
            reactions=available_reactions,
        )