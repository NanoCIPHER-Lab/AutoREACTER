import pytest
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from AutoREACTER.detectors.reactions_library import registry
from AutoREACTER.detectors.reactions_library.polythioesters import REACTIONS


THIOL_HALIDE = (
    "Dithiol and Di-Carboxylic Acid Halide "
    "Polycondensation(Polythioesterification)"
)

THIOL_ACID = (
    "Dithiol and Di-Carboxylic Acid "
    "Polycondensation(Polythioesterification)"
)

HYDROXY_PATH = (
    "Hydroxy-Thiol and Di-Carboxylic Acid Halide "
    "Polycondensation through Hydroxy Group"
)

THIOL_PATH = (
    "Hydroxy-Thiol and Di-Carboxylic Acid Halide "
    "Polycondensation through Thiol Group"
)

EXPECTED_REACTIONS = {
    THIOL_HALIDE,
    THIOL_ACID,
    HYDROXY_PATH,
    THIOL_PATH,
}


def run_reaction(name, smiles_1, smiles_2):
    rxn = rdChemReactions.ReactionFromSmarts(
        REACTIONS[name]["reaction"]
    )
    assert rxn is not None

    return rxn.RunReactants(
        (
            Chem.AddHs(Chem.MolFromSmiles(smiles_1)),
            Chem.AddHs(Chem.MolFromSmiles(smiles_2)),
        )
    )


def test_expected_active_reactions():
    assert set(REACTIONS) == EXPECTED_REACTIONS


@pytest.mark.parametrize("name", sorted(EXPECTED_REACTIONS))
def test_smarts_parses(name):
    assert rdChemReactions.ReactionFromSmarts(
        REACTIONS[name]["reaction"]
    ) is not None


@pytest.mark.parametrize("name", sorted(EXPECTED_REACTIONS))
def test_registry_validation_passes(name):
    assert registry._validate_reaction_smarts(
        name,
        REACTIONS[name],
    ) == []


@pytest.mark.parametrize("name", sorted(EXPECTED_REACTIONS))
def test_reactions_delete_leaving_group(name):
    assert REACTIONS[name]["delete_atom"] is True


def test_dithiol_acid_halide_metadata():
    reaction = REACTIONS[THIOL_HALIDE]

    assert reaction["reactant_1"] == "dithiol"
    assert reaction["reactant_2"] == "di_carboxylic_acid_halide"
    assert reaction["product"] == "polythioester_chain"


def test_dithiol_acid_metadata():
    reaction = REACTIONS[THIOL_ACID]

    assert reaction["reactant_1"] == "dithiol"
    assert reaction["reactant_2"] == "di_carboxylic_acid"
    assert reaction["product"] == "polythioester_chain"


def test_hydroxy_thiol_has_two_distinct_reaction_routes():
    assert "OX2H1" in REACTIONS[HYDROXY_PATH]["reaction"]
    assert "SX2H1" in REACTIONS[THIOL_PATH]["reaction"]


def test_dithiol_acid_halide_reaction_executes():
    # SMARTS order = acid halide first, dithiol second
    products = run_reaction(
        THIOL_HALIDE,
        "O=C(Cl)CCC(=O)Cl",
        "SCCS",
    )

    assert len(products) > 0


def test_dithiol_acid_reaction_executes():
    products = run_reaction(
        THIOL_ACID,
        "O=C(O)CCC(=O)O",
        "SCCS",
    )

    assert len(products) > 0


def test_hydroxy_thiol_hydroxy_route_executes():
    products = run_reaction(
        HYDROXY_PATH,
        "O=C(Cl)CCC(=O)Cl",
        "OCCS",
    )

    assert len(products) > 0


def test_hydroxy_thiol_thiol_route_executes():
    products = run_reaction(
        THIOL_PATH,
        "O=C(Cl)CCC(=O)Cl",
        "OCCS",
    )

    assert len(products) > 0


def test_dithiol_acid_halide_produces_byproduct():
    products = run_reaction(
        THIOL_HALIDE,
        "O=C(Cl)CCC(=O)Cl",
        "SCCS",
    )

    assert len(products[0]) == 2