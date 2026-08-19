from rdkit import Chem
from rdkit.Chem import rdChemReactions

from AutoREACTER.detectors.reactions_library import registry
from AutoREACTER.detectors.reactions_library.polyureas import REACTIONS


NAME = "Di-Amine and Di-Isocyanate Polyaddition(Polyurea Formation)"


def run_reaction(smiles_1, smiles_2):
    rxn = rdChemReactions.ReactionFromSmarts(
        REACTIONS[NAME]["reaction"]
    )

    return rxn.RunReactants(
        (
            Chem.AddHs(Chem.MolFromSmiles(smiles_1)),
            Chem.AddHs(Chem.MolFromSmiles(smiles_2)),
        )
    )


def test_expected_active_reaction():
    assert set(REACTIONS) == {NAME}


def test_metadata():
    reaction = REACTIONS[NAME]

    assert reaction["same_reactants"] is False
    assert reaction["reactant_1"] == "di_amine"
    assert reaction["reactant_2"] == "di_isocyanate"
    assert reaction["product"] == "polyurea_chain"
    assert reaction["delete_atom"] is False


def test_smarts_parses():
    assert rdChemReactions.ReactionFromSmarts(
        REACTIONS[NAME]["reaction"]
    ) is not None


def test_registry_validation_passes():
    assert registry._validate_reaction_smarts(
        NAME,
        REACTIONS[NAME],
    ) == []


def test_diamine_diisocyanate_reaction_executes():
    products = run_reaction(
        "NCCN",
        "O=C=NCCCCCCN=C=O",
    )

    assert len(products) > 0


def test_products_are_sanitizable():
    products = run_reaction(
        "NCCN",
        "O=C=NCCCCCCN=C=O",
    )

    for product in products[0]:
        Chem.SanitizeMol(product)


def test_nonamine_does_not_react():
    products = run_reaction(
        "OCCO",
        "O=C=NCCCCCCN=C=O",
    )

    assert len(products) == 0