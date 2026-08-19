import pytest
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from AutoREACTER.detectors.reactions_library import registry
from AutoREACTER.detectors.reactions_library.polyanhydrides import (
    REACTIONS,
)


HOMO = (
    "Carboxylic Acid and Acid Halide Polycondensation "
    "(Polyanhydride Formation)"
)

CO = (
    "Carboxylic Acid and Acid Halide Copolycondensation "
    "(Polyanhydride Copolymerization)"
)


EXPECTED_REACTIONS = {
    HOMO,
    CO,
}


def run_reaction(reaction_name, smiles_1, smiles_2):
    rxn = rdChemReactions.ReactionFromSmarts(
        REACTIONS[reaction_name]["reaction"]
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
        "product",
        "delete_atom",
        "reaction",
        "reference",
    }

    assert required.issubset(
        reaction.keys()
    )


@pytest.mark.parametrize(
    "reaction_name",
    sorted(EXPECTED_REACTIONS),
)
def test_smarts_parses(reaction_name):
    rxn = rdChemReactions.ReactionFromSmarts(
        REACTIONS[reaction_name]["reaction"]
    )

    assert rxn is not None


@pytest.mark.parametrize(
    "reaction_name",
    sorted(EXPECTED_REACTIONS),
)
def test_registry_validation_passes(reaction_name):
    errors = registry._validate_reaction_smarts(
        reaction_name,
        REACTIONS[reaction_name],
    )

    assert errors == []


def test_homopolymerization_metadata():
    reaction = REACTIONS[HOMO]

    assert reaction["same_reactants"] is True
    assert (
        reaction["reactant_1"]
        == "carboxylic_acid_acid_halide"
    )
    assert reaction.get("reactant_2") is None
    assert reaction["product"] == "polyanhydride_chain"
    assert reaction["delete_atom"] is True


def test_copolymerization_metadata():
    reaction = REACTIONS[CO]

    assert reaction["same_reactants"] is False

    assert (
        reaction["reactant_1"]
        == "carboxylic_acid_acid_halide"
    )

    assert (
        reaction["reactant_2"]
        == "carboxylic_acid_acid_halide"
    )

    assert reaction["product"] == "polyanhydride_chain"
    assert reaction["delete_atom"] is True


@pytest.mark.parametrize(
    "reaction_name",
    sorted(EXPECTED_REACTIONS),
)
def test_initiator_maps_are_present(reaction_name):
    rxn = rdChemReactions.ReactionFromSmarts(
        REACTIONS[reaction_name]["reaction"]
    )

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

    assert {1, 2}.issubset(
        registry._atom_maps_in_templates(
            reactants
        )
    )

    assert {1, 2}.issubset(
        registry._atom_maps_in_templates(
            products
        )
    )


@pytest.mark.parametrize(
    "reaction_name",
    sorted(EXPECTED_REACTIONS),
)
def test_product_contains_initiator_bond(
    reaction_name,
):
    rxn = rdChemReactions.ReactionFromSmarts(
        REACTIONS[reaction_name]["reaction"]
    )

    products = [
        Chem.Mol(
            rxn.GetProductTemplate(i)
        )
        for i in range(
            rxn.GetNumProductTemplates()
        )
    ]

    assert registry._has_bond_between_atom_maps(
        products,
        1,
        2,
    ) is True


def test_homopolymerization_smarts_executes():
    """
    HOOC-(CH2)3-COCl contains both groups required by the
    heterobifunctional polyanhydride definition.
    """
    products = run_reaction(
        HOMO,
        "O=C(O)CCC(=O)Cl",
        "O=C(O)CCC(=O)Cl",
    )

    assert len(products) > 0


def test_homopolymerization_produces_byproduct():
    products = run_reaction(
        HOMO,
        "O=C(O)CCC(=O)Cl",
        "O=C(O)CCC(=O)Cl",
    )

    assert len(products[0]) == 2


def test_products_are_sanitizable():
    products = run_reaction(
        HOMO,
        "O=C(O)CCC(=O)Cl",
        "O=C(O)CCC(=O)Cl",
    )

    for product in products[0]:
        Chem.SanitizeMol(product)


def test_reaction_fails_without_acid_halide():
    products = run_reaction(
        HOMO,
        "O=C(O)CCC(=O)O",
        "O=C(O)CCC(=O)O",
    )

    assert len(products) == 0