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
    from .polysaccharides import REACTIONS as POLYSACCHARIDES

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
    POLYSACCHARIDES,
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
        - Atom maps 1 and 2 are the default REACTER initiators.
        - Initiator atoms must belong to different reactant molecules.
        - Initiator atoms must be preserved in the product.
        - Initiator atoms do not need to form the new bond.

    A reaction can override the default initiator maps using:
        "initiator_atom_maps": (a, b)
    """
    errors: list[str] = []

    smarts = reaction.get("reaction")

    if not smarts:
        return [f"{reaction_name}: missing required key 'reaction'"]

    if rdChemReactions is None:
        return [
            f"{reaction_name}: RDKit is required to validate reaction SMARTS"
        ]

    initiator_atom_maps = reaction.get("initiator_atom_maps", (1, 2))

    if not isinstance(initiator_atom_maps, (tuple, list)):
        return [
            f"{reaction_name}: initiator_atom_maps must contain two atom maps"
        ]

    if len(initiator_atom_maps) != 2:
        return [
            f"{reaction_name}: initiator_atom_maps must contain exactly two atom maps"
        ]

    try:
        initiator_1, initiator_2 = map(int, initiator_atom_maps)
    except (TypeError, ValueError):
        return [
            f"{reaction_name}: initiator atom maps must be integers"
        ]

    if initiator_1 <= 0 or initiator_2 <= 0:
        errors.append(
            f"{reaction_name}: initiator atom maps must be positive"
        )

    if initiator_1 == initiator_2:
        errors.append(
            f"{reaction_name}: initiator atom maps must be different"
        )

    if errors:
        return errors

    try:
        rdkit_reaction = rdChemReactions.ReactionFromSmarts(smarts)
    except Exception as error:
        return [
            f"{reaction_name}: invalid reaction SMARTS: {error}"
        ]

    if rdkit_reaction is None:
        return [
            f"{reaction_name}: RDKit could not parse reaction SMARTS"
        ]

    reactant_templates = [
        rdkit_reaction.GetReactantTemplate(i)
        for i in range(rdkit_reaction.GetNumReactantTemplates())
    ]

    product_templates = [
        rdkit_reaction.GetProductTemplate(i)
        for i in range(rdkit_reaction.GetNumProductTemplates())
    ]

    # Find the reactant molecule containing each initiator atom.
    initiator_locations = {
        initiator_1: [],
        initiator_2: [],
    }

    for molecule_index, template in enumerate(reactant_templates):
        for atom in template.GetAtoms():
            atom_map = atom.GetAtomMapNum()

            if atom_map in initiator_locations:
                initiator_locations[atom_map].append(molecule_index)

    # Each initiator must appear exactly once in the reactants.
    for atom_map, locations in initiator_locations.items():
        if not locations:
            errors.append(
                f"{reaction_name}: initiator atom map {atom_map} "
                "is missing from reactants"
            )
        elif len(locations) != 1:
            errors.append(
                f"{reaction_name}: initiator atom map {atom_map} "
                "must appear exactly once in reactants"
            )

    # Initiators must belong to different reactant molecules.
    locations_1 = initiator_locations[initiator_1]
    locations_2 = initiator_locations[initiator_2]

    if len(locations_1) == 1 and len(locations_2) == 1:
        if locations_1[0] == locations_2[0]:
            errors.append(
                f"{reaction_name}: initiator atom maps "
                f"{initiator_1} and {initiator_2} "
                "must belong to different reactant molecules"
            )

    # Initiator atoms must still exist in the products.
    product_maps = _atom_maps_in_templates(product_templates)

    for atom_map in (initiator_1, initiator_2):
        if atom_map not in product_maps:
            errors.append(
                f"{reaction_name}: initiator atom map {atom_map} "
                "is missing from products"
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