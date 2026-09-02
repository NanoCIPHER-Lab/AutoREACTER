from rdkit import Chem
from rdkit.Chem import rdChemReactions

from AutoREACTER.detectors.reactions_library import registry
from AutoREACTER.detectors.reactions_library.thiol_ene_polymers import (
    REACTIONS,
)


NAME = "Dithiol and Diene Thiol-Ene Click Polymerization"


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
    assert reaction["reactant_1"] == "dithiol"
    assert reaction["reactant_2"] == "diene"
    assert reaction["product"] == "poly_thioether_chain"
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


def test_dithiol_diene_reaction_executes():
    products = run_reaction(
        "SCCS",
        "C=CC=C",
    )

    assert len(products) > 0


def test_products_are_sanitizable():
    products = run_reaction(
        "SCCS",
        "C=CC=C",
    )

    for product in products[0]:
        Chem.SanitizeMol(product)


def test_dithiol_diene_reaction_rejects_saturated_hydrocarbon():
    products = run_reaction(
        "SCCS",
        "CCCC",
    )

    assert len(products) == 0


def test_reaction_rejects_thioether_without_sh():
    products = run_reaction(
        "CSC",
        "C=CC=C",
    )

    assert len(products) == 0