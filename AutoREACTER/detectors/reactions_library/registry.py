"""Aggregate and validate all polymer-family reaction modules."""

from __future__ import annotations

try:
    from rdkit.Chem import rdChemReactions
except ImportError:  # pragma: no cover
    rdChemReactions = None


try:
    from .polyesters import REACTIONS as POLYESTERS
    # from .polyethers import REACTIONS as POLYETHERS
    from .polyamides import REACTIONS as POLYAMIDES
    from .polyanhydrides import REACTIONS as POLYANHYDRIDES
    from .polythioesters import REACTIONS as POLYTHIOESTERS
    from .polyurethanes import REACTIONS as POLYURETHANES
    from .polyureas import REACTIONS as POLYUREAS
    from .epoxy_polymers import REACTIONS as EPOXY_POLYMERS
    from .vinyl_polymers import REACTIONS as VINYL_POLYMERS
    from .polycarbonates import REACTIONS as POLYCARBONATES
    # from .polyimides import REACTIONS as POLYIMIDES
    # from .polybenzimidazoles import REACTIONS as POLYBENZIMIDAZOLES
    # from .phenolic_resins import REACTIONS as PHENOLIC_RESINS
    from .polysiloxanes import REACTIONS as POLYSILOXANES
    # from .polysulfides import REACTIONS as POLYSULFIDES
    from .thiol_ene_polymers import REACTIONS as THIOL_ENE_POLYMERS
    # from .metathesis_polymers import REACTIONS as METATHESIS_POLYMERS
    # from .cycloaddition_polymers import REACTIONS as CYCLOADDITION_POLYMERS

except ImportError as e:
    from polyesters import REACTIONS as POLYESTERS
    from polyamides import REACTIONS as POLYAMIDES
    from polyanhydrides import REACTIONS as POLYANHYDRIDES
    from polythioesters import REACTIONS as POLYTHIOESTERS
    from polyurethanes import REACTIONS as POLYURETHANES
    from polyureas import REACTIONS as POLYUREAS
    from epoxy_polymers import REACTIONS as EPOXY_POLYMERS
    from vinyl_polymers import REACTIONS as VINYL_POLYMERS
    from polycarbonates import REACTIONS as POLYCARBONATES
    from polysiloxanes import REACTIONS as POLYSILOXANES
    from thiol_ene_polymers import REACTIONS as THIOL_ENE_POLYMERS


_REACTION_MODULES = [
    POLYESTERS,
    # POLYETHERS,
    POLYAMIDES,
    POLYANHYDRIDES,
    POLYTHIOESTERS,
    POLYURETHANES,
    POLYUREAS,
    EPOXY_POLYMERS,
    VINYL_POLYMERS,
    POLYCARBONATES,
    # POLYIMIDES,
    # POLYBENZIMIDAZOLES,
    # PHENOLIC_RESINS,
    POLYSILOXANES,
    # POLYSULFIDES,
    THIOL_ENE_POLYMERS,
    # METATHESIS_POLYMERS,
    # CYCLOADDITION_POLYMERS,
]


class ReactionLibraryValidationError(ValueError):
    """Raised when a reaction-library SMARTS violates AutoREACTER rules."""


def _atom_maps_in_templates(templates) -> set[int]:
    """Return all nonzero atom-map numbers present in RDKit templates."""
    atom_maps: set[int] = set()

    for template in templates:
        for atom in template.GetAtoms():
            atom_map = atom.GetAtomMapNum()
            if atom_map:
                atom_maps.add(atom_map)

    return atom_maps


def _has_bond_between_atom_maps(
    templates,
    atom_map_1: int,
    atom_map_2: int,
) -> bool:
    """Return True if any template contains a bond between two atom maps."""
    target = {atom_map_1, atom_map_2}

    for template in templates:
        for bond in template.GetBonds():
            begin_map = bond.GetBeginAtom().GetAtomMapNum()
            end_map = bond.GetEndAtom().GetAtomMapNum()

            if {begin_map, end_map} == target:
                return True

    return False


def _validate_reaction_smarts(
    reaction_name: str,
    reaction: dict,
) -> list[str]:
    """
    Validate one reaction-library entry.

    AutoREACTER convention:
        atom maps :1 and :2 are reserved as LAMMPS bond/react initiator atoms.

    By default this validator requires:
        - reaction["reaction"] exists
        - maps 1 and 2 exist in reactants
        - maps 1 and 2 exist in products
        - products contain a bond between map 1 and map 2

    A reaction can override this with:
        "initiator_atom_maps": (a, b)

    A special reaction can skip this check with:
        "validate_initiator_bond": False
    """
    errors: list[str] = []

    smarts = reaction.get("reaction")

    if not smarts:
        return [f"{reaction_name}: missing required key 'reaction'"]

    if reaction.get("validate_initiator_bond", True) is False:
        return errors

    if rdChemReactions is None:
        return [f"{reaction_name}: RDKit is required to validate reaction SMARTS"]

    initiator_atom_maps = reaction.get("initiator_atom_maps", (1, 2))

    if len(initiator_atom_maps) != 2:
        return [
            f"{reaction_name}: initiator_atom_maps must contain exactly two atom maps"
        ]

    initiator_1, initiator_2 = map(int, initiator_atom_maps)

    try:
        rdkit_reaction = rdChemReactions.ReactionFromSmarts(smarts)
    except Exception as error:
        return [f"{reaction_name}: invalid reaction SMARTS: {error}"]

    if rdkit_reaction is None:
        return [f"{reaction_name}: RDKit could not parse reaction SMARTS"]

    reactant_templates = [
        rdkit_reaction.GetReactantTemplate(i)
        for i in range(rdkit_reaction.GetNumReactantTemplates())
    ]

    product_templates = [
        rdkit_reaction.GetProductTemplate(i)
        for i in range(rdkit_reaction.GetNumProductTemplates())
    ]

    reactant_maps = _atom_maps_in_templates(reactant_templates)
    product_maps = _atom_maps_in_templates(product_templates)

    required_maps = {initiator_1, initiator_2}

    missing_reactant_maps = required_maps - reactant_maps
    missing_product_maps = required_maps - product_maps

    if missing_reactant_maps:
        errors.append(
            f"{reaction_name}: initiator atom maps missing from reactants: "
            f"{sorted(missing_reactant_maps)}"
        )

    if missing_product_maps:
        errors.append(
            f"{reaction_name}: initiator atom maps missing from products: "
            f"{sorted(missing_product_maps)}"
        )

    product_has_initiator_bond = _has_bond_between_atom_maps(
        product_templates,
        initiator_1,
        initiator_2,
    )

    if not product_has_initiator_bond:
        errors.append(
            f"{reaction_name}: product does not contain required "
            f"AutoREACTER initiator bond between atom maps "
            f"{initiator_1} and {initiator_2}"
        )

    return errors


def validate_reactions(reactions: dict) -> None:
    """Validate the merged AutoREACTER reaction library."""
    errors: list[str] = []

    for reaction_name, reaction in reactions.items():
        if not isinstance(reaction, dict):
            errors.append(f"{reaction_name}: reaction entry must be a dictionary")
            continue

        errors.extend(_validate_reaction_smarts(reaction_name, reaction))

    if errors:
        message = "\n".join(f"  - {error}" for error in errors)
        raise ReactionLibraryValidationError(
            "Reaction library validation failed:\n" + message
        )


def load_reactions() -> dict:
    """Return one flat reaction dictionary with duplicate-name protection."""
    merged = {}

    for module in _REACTION_MODULES:
        for reaction_name, reaction in module.items():
            if reaction_name in merged:
                raise ValueError(f"Duplicate reaction name: {reaction_name}")
            merged[reaction_name] = reaction

    validate_reactions(merged)

    return merged


REACTIONS = load_reactions()


class ReactionLibrary:
    """Backward-compatible class exposing ``self.reactions``."""

    def __init__(self):
        self.reactions = load_reactions()


if __name__ == "__main__":
    REACTIONS = load_reactions()
    num = 0

    with open("reactions.txt", "w") as f:
        for reaction in REACTIONS.items():
            f.write(str(reaction) + "\n")

    reaction_len = len(REACTIONS)

    import os

    file_abs_path = os.path.abspath("reactions.txt")
    print(
        f"reactions.txt has been written to {file_abs_path}, "
        f"num reactions: {reaction_len}"
    )