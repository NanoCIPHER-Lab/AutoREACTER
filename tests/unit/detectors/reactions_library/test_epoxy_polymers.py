import pytest
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from AutoREACTER.detectors.reactions_library import registry
from AutoREACTER.detectors.reactions_library.epoxy_polymers import (
    REACTIONS,
)


PRIMARY = (
    "Primary Amine and Epoxide Polyaddition "
    "(Epoxy-Amine, First Addition)"
)

SECONDARY = (
    "Secondary Amine and Epoxide Polyaddition "
    "(Epoxy-Amine, Second Addition / Crosslink)"
)


EXPECTED_REACTIONS = {
    PRIMARY,
    SECONDARY,
}


def run_reaction(reaction_name, smiles_1, smiles_2):
    reaction_info = REACTIONS[reaction_name]

    rxn = rdChemReactions.ReactionFromSmarts(
        reaction_info["reaction"]
    )

    assert rxn is not None

    mol_1 = Chem.AddHs(
        Chem.MolFromSmiles(smiles_1)
    )

    mol_2 = Chem.AddHs(
        Chem.MolFromSmiles(smiles_2)
    )

    return rxn.RunReactants(
        (
            mol_1,
            mol_2,
        )
    )


# =============================================================================
# Structure
# =============================================================================

def test_expected_active_reactions():
    assert set(REACTIONS) == EXPECTED_REACTIONS


@pytest.mark.parametrize(
    "reaction_name",
    sorted(EXPECTED_REACTIONS),
)
def test_required_fields_exist(reaction_name):
    reaction = REACTIONS[reaction_name]

    required = {
        "same_reactants",
        "reactant_1",
        "reactant_2",
        "product",
        "delete_atom",
        "reaction",
        "reference",
        "comments",
    }

    assert required.issubset(
        reaction.keys()
    )


@pytest.mark.parametrize(
    "reaction_name",
    sorted(EXPECTED_REACTIONS),
)
def test_reaction_smarts_parses(reaction_name):
    reaction = REACTIONS[reaction_name]

    rxn = rdChemReactions.ReactionFromSmarts(
        reaction["reaction"]
    )

    assert rxn is not None


@pytest.mark.parametrize(
    "reaction_name",
    sorted(EXPECTED_REACTIONS),
)
def test_reaction_passes_registry_validation(reaction_name):
    errors = registry._validate_reaction_smarts(
        reaction_name,
        REACTIONS[reaction_name],
    )

    assert errors == []


# =============================================================================
# Metadata
# =============================================================================

def test_primary_epoxy_amine_metadata():
    reaction = REACTIONS[PRIMARY]

    assert reaction["same_reactants"] is False
    assert reaction["reactant_1"] == "primary_amine"
    assert reaction["reactant_2"] == "di_epoxide"

    assert (
        reaction["product"]
        == "secondary_amine_hydroxyl_product"
    )

    assert reaction["delete_atom"] is False
    assert reaction["comments"] is None


def test_secondary_epoxy_amine_metadata():
    reaction = REACTIONS[SECONDARY]

    assert reaction["same_reactants"] is False
    assert reaction["reactant_1"] == "secondary_amine"
    assert reaction["reactant_2"] == "di_epoxide"

    assert (
        reaction["product"]
        == "tertiary_amine_crosslink_product"
    )

    assert reaction["delete_atom"] is False
    assert reaction["comments"] is None


def test_primary_and_secondary_reactions_use_different_amine_queries():
    primary = REACTIONS[PRIMARY]["reaction"]
    secondary = REACTIONS[SECONDARY]["reaction"]

    assert "[NX3H2:1]" in primary
    assert "[NX3H1:1]" in secondary


# =============================================================================
# Initiator-map convention
# =============================================================================

@pytest.mark.parametrize(
    "reaction_name",
    sorted(EXPECTED_REACTIONS),
)
def test_reaction_contains_initiator_maps_1_and_2(
    reaction_name,
):
    reaction = rdChemReactions.ReactionFromSmarts(
        REACTIONS[reaction_name]["reaction"]
    )

    reactants = [
        Chem.Mol(
            reaction.GetReactantTemplate(i)
        )
        for i in range(
            reaction.GetNumReactantTemplates()
        )
    ]

    products = [
        Chem.Mol(
            reaction.GetProductTemplate(i)
        )
        for i in range(
            reaction.GetNumProductTemplates()
        )
    ]

    assert {
        1,
        2,
    }.issubset(
        registry._atom_maps_in_templates(
            reactants
        )
    )

    assert {
        1,
        2,
    }.issubset(
        registry._atom_maps_in_templates(
            products
        )
    )


@pytest.mark.parametrize(
    "reaction_name",
    sorted(EXPECTED_REACTIONS),
)
def test_product_contains_new_initiator_bond(
    reaction_name,
):
    reaction = rdChemReactions.ReactionFromSmarts(
        REACTIONS[reaction_name]["reaction"]
    )

    products = [
        Chem.Mol(
            reaction.GetProductTemplate(i)
        )
        for i in range(
            reaction.GetNumProductTemplates()
        )
    ]

    assert registry._has_bond_between_atom_maps(
        products,
        1,
        2,
    ) is True


# =============================================================================
# Primary amine + epoxide chemistry
# =============================================================================

def test_primary_amine_reaction_executes():
    """
    Ethylamine contains a primary amine.

    Propylene oxide provides the CH2/CH epoxide environment required by
    this reaction SMARTS.
    """
    products = run_reaction(
        PRIMARY,
        "CCN",
        "CC1CO1",
    )

    assert len(products) > 0


def test_primary_reaction_products_are_sanitizable():
    products = run_reaction(
        PRIMARY,
        "CCN",
        "CC1CO1",
    )

    first_product_set = products[0]

    assert first_product_set

    for product in first_product_set:
        Chem.SanitizeMol(product)


def test_primary_reaction_does_not_accept_secondary_amine():
    products = run_reaction(
        PRIMARY,
        "CNC",
        "CC1CO1",
    )

    assert len(products) == 0


def test_primary_reaction_does_not_accept_nonamine():
    products = run_reaction(
        PRIMARY,
        "CCO",
        "CC1CO1",
    )

    assert len(products) == 0


# =============================================================================
# Secondary amine + epoxide chemistry
# =============================================================================

def test_secondary_amine_reaction_executes():
    """
    Dimethylamine contains one N-H bond and therefore matches NX3H1.
    """
    products = run_reaction(
        SECONDARY,
        "CNC",
        "CC1CO1",
    )

    assert len(products) > 0


def test_secondary_reaction_products_are_sanitizable():
    products = run_reaction(
        SECONDARY,
        "CNC",
        "CC1CO1",
    )

    first_product_set = products[0]

    assert first_product_set

    for product in first_product_set:
        Chem.SanitizeMol(product)


def test_secondary_reaction_does_not_accept_primary_amine():
    products = run_reaction(
        SECONDARY,
        "CCN",
        "CC1CO1",
    )

    assert len(products) == 0


def test_secondary_reaction_does_not_accept_tertiary_amine():
    products = run_reaction(
        SECONDARY,
        "CN(C)C",
        "CC1CO1",
    )

    assert len(products) == 0