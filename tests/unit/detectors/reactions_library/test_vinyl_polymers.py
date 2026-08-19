import pytest
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from AutoREACTER.detectors.reactions_library import registry
from AutoREACTER.detectors.reactions_library.vinyl_polymers import (
    REACTIONS,
)


INITIATION = "Vinyl Addition Polymerization Initiation"
PROPAGATION = "Vinyl Addition Polymerization Propagation"
CO_INITIATION = "Vinyl Copolymerization Initiation"
TFE_INITIATION = "Tetrafluoroethylene Initiation"
TFE_PROPAGATION = "Tetrafluoroethylene Propagation"


EXPECTED_REACTIONS = {
    INITIATION,
    PROPAGATION,
    CO_INITIATION,
    TFE_INITIATION,
    TFE_PROPAGATION,
}


def reaction_templates(name):
    rxn = rdChemReactions.ReactionFromSmarts(
        REACTIONS[name]["reaction"]
    )

    assert rxn is not None

    reactants = [
        Chem.Mol(
            rxn.GetReactantTemplate(i)
        )
        for i in range(
            rxn.GetNumReactantTemplates()
        )
    ]

    products = [
        Chem.Mol(
            rxn.GetProductTemplate(i)
        )
        for i in range(
            rxn.GetNumProductTemplates()
        )
    ]

    return reactants, products


# =============================================================================
# Structure
# =============================================================================

def test_expected_active_reactions():
    """
    This deliberately catches accidental dictionary nesting.

    All five currently enabled vinyl/TFE reactions are intended to be
    top-level reaction-library entries.
    """
    assert set(REACTIONS) == EXPECTED_REACTIONS


@pytest.mark.parametrize(
    "name",
    sorted(EXPECTED_REACTIONS),
)
def test_each_reaction_is_top_level(name):
    assert name in REACTIONS
    assert isinstance(REACTIONS[name], dict)


@pytest.mark.parametrize(
    "name",
    sorted(EXPECTED_REACTIONS),
)
def test_smarts_parses(name):
    assert rdChemReactions.ReactionFromSmarts(
        REACTIONS[name]["reaction"]
    ) is not None


@pytest.mark.parametrize(
    "name",
    sorted(EXPECTED_REACTIONS),
)
def test_registry_validation_passes(name):
    assert registry._validate_reaction_smarts(
        name,
        REACTIONS[name],
    ) == []


# =============================================================================
# Metadata
# =============================================================================

def test_vinyl_initiation_metadata():
    reaction = REACTIONS[INITIATION]

    assert reaction["same_reactants"] is True
    assert reaction["reactant_1"] == "vinyl"
    assert reaction.get("reactant_2") is None
    assert reaction["product"] == "vinyl_chain_end_radical"
    assert reaction["delete_atom"] is False


def test_vinyl_propagation_metadata():
    reaction = REACTIONS[PROPAGATION]

    assert reaction["same_reactants"] is False
    assert reaction["reactant_1"] == "vinyl"

    assert (
        reaction["reactant_2"]
        == "vinyl_chain_end_radical"
    )

    assert reaction["product"] == "vinyl_chain_end_radical"


def test_copolymerization_initiation_metadata():
    reaction = REACTIONS[CO_INITIATION]

    assert reaction["same_reactants"] is False
    assert reaction["reactant_1"] == "vinyl"
    assert reaction["reactant_2"] == "vinyl"
    assert reaction["product"] == "vinyl_chain_end_radical"


def test_tfe_initiation_metadata():
    reaction = REACTIONS[TFE_INITIATION]

    assert reaction["reactant_1"] == "tetrafluoroethylene"
    assert reaction["reactant_2"] == "tetrafluoroethylene"
    assert reaction["product"] == "vinyl_chain_end_radical"


def test_tfe_propagation_metadata():
    reaction = REACTIONS[TFE_PROPAGATION]

    assert reaction["reactant_1"] == "vinyl_chain_end_radical"
    assert reaction["reactant_2"] == "tetrafluoroethylene"
    assert reaction["product"] == "vinyl_chain_end_radical"


# =============================================================================
# H-T vinyl initiation topology
# =============================================================================

def test_vinyl_initiation_forms_new_bond_between_maps_1_and_2():
    reactants, products = reaction_templates(
        INITIATION
    )

    assert registry._has_bond_between_atom_maps(
        reactants,
        1,
        2,
    ) is False

    assert registry._has_bond_between_atom_maps(
        products,
        1,
        2,
    ) is True


def test_vinyl_initiation_caps_map_3_as_ch3():
    _, products = reaction_templates(
        INITIATION
    )

    map_3_atoms = [
        atom
        for mol in products
        for atom in mol.GetAtoms()
        if atom.GetAtomMapNum() == 3
    ]

    assert len(map_3_atoms) == 1

    atom = map_3_atoms[0]

    assert atom.GetSymbol() == "C"


# =============================================================================
# Vinyl propagation topology
# =============================================================================

def test_vinyl_propagation_forms_new_1_2_bond():
    reactants, products = reaction_templates(
        PROPAGATION
    )

    assert registry._has_bond_between_atom_maps(
        reactants,
        1,
        2,
    ) is False

    assert registry._has_bond_between_atom_maps(
        products,
        1,
        2,
    ) is True


def test_vinyl_propagation_retains_new_active_map_3():
    _, products = reaction_templates(
        PROPAGATION
    )

    product_maps = registry._atom_maps_in_templates(
        products
    )

    assert 3 in product_maps


# =============================================================================
# Vinyl copolymerization initiation
# =============================================================================

def test_copolymerization_initiation_forms_new_1_2_bond():
    reactants, products = reaction_templates(
        CO_INITIATION
    )

    assert registry._has_bond_between_atom_maps(
        reactants,
        1,
        2,
    ) is False

    assert registry._has_bond_between_atom_maps(
        products,
        1,
        2,
    ) is True


def test_copolymerization_product_retains_radical_map_4():
    _, products = reaction_templates(
        CO_INITIATION
    )

    assert 4 in registry._atom_maps_in_templates(
        products
    )


# =============================================================================
# TFE initiation
# =============================================================================

def test_tfe_initiation_forms_new_1_2_bond():
    reactants, products = reaction_templates(
        TFE_INITIATION
    )

    assert registry._has_bond_between_atom_maps(
        reactants,
        1,
        2,
    ) is False

    assert registry._has_bond_between_atom_maps(
        products,
        1,
        2,
    ) is True


def test_tfe_initiation_preserves_all_four_fluorines():
    _, products = reaction_templates(
        TFE_INITIATION
    )

    maps = registry._atom_maps_in_templates(
        products
    )

    assert {
        4,
        5,
        6,
        7,
    }.issubset(maps)


# =============================================================================
# TFE propagation topology
# =============================================================================

def test_tfe_propagation_forms_new_bond_between_maps_2_and_3():
    """
    In the actual SMARTS, 1-2 already exists in the growing chain.

    The propagation step forms the new connection between maps 2 and 3.
    """
    reactants, products = reaction_templates(
        TFE_PROPAGATION
    )

    assert registry._has_bond_between_atom_maps(
        reactants,
        2,
        3,
    ) is False

    assert registry._has_bond_between_atom_maps(
        products,
        2,
        3,
    ) is True


def test_tfe_propagation_maps_1_and_2_are_already_bonded_before_reaction():
    reactants, _ = reaction_templates(
        TFE_PROPAGATION
    )

    assert registry._has_bond_between_atom_maps(
        reactants,
        1,
        2,
    ) is True


def test_tfe_propagation_should_declare_new_initiator_maps():
    """
    The registry defaults initiator_atom_maps to (1, 2).

    But TFE propagation's actual NEW bond is 2-3, while 1-2 already
    exists before the reaction. Therefore this reaction should explicitly
    override the default initiator mapping.
    """
    reaction = REACTIONS[TFE_PROPAGATION]

    assert reaction.get(
        "initiator_atom_maps"
    ) == (2, 3)