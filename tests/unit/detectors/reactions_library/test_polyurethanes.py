import pytest
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from AutoREACTER.detectors.reactions_library import registry
from AutoREACTER.detectors.reactions_library.polyurethanes import REACTIONS


URETHANE = "Diol and Di-Isocyanate Polyaddition(Polyurethane Formation)"
THIOURETHANE = (
    "Dithiol and Di-Isocyanate "
    "Polyaddition(Polythiourethane Formation)"
)

EXPECTED_REACTIONS = {
    URETHANE,
    THIOURETHANE,
}


def run_reaction(name, smiles_1, smiles_2):
    rxn = rdChemReactions.ReactionFromSmarts(
        REACTIONS[name]["reaction"]
    )

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


def test_polyurethane_metadata():
    reaction = REACTIONS[URETHANE]

    assert reaction["reactant_1"] == "diol"
    assert reaction["reactant_2"] == "di_isocyanate"
    assert reaction["product"] == "polyurethane_chain"
    assert reaction["delete_atom"] is False


def test_polythiourethane_metadata():
    reaction = REACTIONS[THIOURETHANE]

    assert reaction["reactant_1"] == "dithiol"
    assert reaction["reactant_2"] == "di_isocyanate"
    assert reaction["product"] == "polythiourethane_chain"
    assert reaction["delete_atom"] is False


def test_diol_diisocyanate_reaction_executes():
    products = run_reaction(
        URETHANE,
        "OCCO",
        "O=C=NCCCCCCN=C=O",
    )

    assert len(products) > 0


def test_dithiol_diisocyanate_reaction_executes():
    products = run_reaction(
        THIOURETHANE,
        "SCCS",
        "O=C=NCCCCCCN=C=O",
    )

    assert len(products) > 0


def test_polyurethane_rejects_dithiol():
    products = run_reaction(
        URETHANE,
        "SCCS",
        "O=C=NCCCCCCN=C=O",
    )

    assert len(products) == 0


def test_polythiourethane_rejects_diol():
    products = run_reaction(
        THIOURETHANE,
        "OCCO",
        "O=C=NCCCCCCN=C=O",
    )

    assert len(products) == 0