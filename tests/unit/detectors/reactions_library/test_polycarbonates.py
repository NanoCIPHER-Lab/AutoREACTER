import pytest
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from AutoREACTER.detectors.reactions_library import registry
from AutoREACTER.detectors.reactions_library.polycarbonates import (
    REACTIONS,
)


PHOSGENE = (
    "Diol and Phosgene "
    "Polycondensation(Polycarbonate Formation)"
)

DIPHENYL = (
    "Diol and Diphenyl Carbonate "
    "Polycondensation(Transcarbonation)"
)


EXPECTED_REACTIONS = {
    PHOSGENE,
    DIPHENYL,
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
def test_smarts_parses(reaction_name):
    assert rdChemReactions.ReactionFromSmarts(
        REACTIONS[reaction_name]["reaction"]
    ) is not None


@pytest.mark.parametrize(
    "reaction_name",
    sorted(EXPECTED_REACTIONS),
)
def test_registry_validation_passes(reaction_name):
    assert registry._validate_reaction_smarts(
        reaction_name,
        REACTIONS[reaction_name],
    ) == []


def test_phosgene_reaction_metadata():
    reaction = REACTIONS[PHOSGENE]

    assert reaction["same_reactants"] is False
    assert reaction["reactant_1"] == "diol"
    assert reaction["reactant_2"] == "phosgene"
    assert reaction["product"] == "polycarbonate_chain"
    assert reaction["delete_atom"] is True


def test_diphenyl_carbonate_metadata():
    reaction = REACTIONS[DIPHENYL]

    assert reaction["same_reactants"] is False
    assert reaction["reactant_1"] == "diol"

    assert (
        reaction["reactant_2"]
        == "diphenyl_carbonate"
    )

    assert reaction["product"] == "polycarbonate_chain"
    assert reaction["delete_atom"] is True


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


def test_diol_phosgene_reaction_executes():
    products = run_reaction(
        PHOSGENE,
        "OCCO",
        "O=C(Cl)Cl",
    )

    assert len(products) > 0


def test_diol_phosgene_produces_hcl_byproduct():
    products = run_reaction(
        PHOSGENE,
        "OCCO",
        "O=C(Cl)Cl",
    )

    assert len(products[0]) == 2


def test_diol_phosgene_products_are_sanitizable():
    products = run_reaction(
        PHOSGENE,
        "OCCO",
        "O=C(Cl)Cl",
    )

    for product in products[0]:
        Chem.SanitizeMol(product)


def test_diol_diphenyl_carbonate_reaction_executes():
    products = run_reaction(
        DIPHENYL,
        "OCCO",
        "O=C(Oc1ccccc1)Oc1ccccc1",
    )

    assert len(products) > 0


def test_transcarbonation_produces_two_product_fragments():
    products = run_reaction(
        DIPHENYL,
        "OCCO",
        "O=C(Oc1ccccc1)Oc1ccccc1",
    )

    assert len(products[0]) == 2


def test_transcarbonation_products_are_sanitizable():
    products = run_reaction(
        DIPHENYL,
        "OCCO",
        "O=C(Oc1ccccc1)Oc1ccccc1",
    )

    for product in products[0]:
        Chem.SanitizeMol(product)


def test_phosgene_reaction_rejects_non_diol_alcohol_free_molecule():
    products = run_reaction(
        PHOSGENE,
        "CCC",
        "O=C(Cl)Cl",
    )

    assert len(products) == 0