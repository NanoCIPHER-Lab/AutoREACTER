import pytest
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from AutoREACTER.detectors.reactions_library import registry
from AutoREACTER.detectors.reactions_library.polyesters import (
    REACTIONS,
)


HYDROXY_ACID_HOMO = (
    "Hydroxy Carboxylic Acid "
    "Polycondensation(Polyesterification)"
)

HYDROXY_ACID_CO = (
    "Hydroxy Carboxylic and Hydroxy Carboxylic "
    "Polycondensation(Polyesterification)"
)

HYDROXY_HALIDE_HOMO = (
    "Hydroxy Acid Halides "
    "Polycondensation(Polyesterification)"
)

HYDROXY_HALIDE_CO = (
    "Hydroxy Acid Halides Hydroxy Acid Halides "
    "Polycondensation(Polyesterification)"
)

DIOL_DIACID = (
    "Diol and Di-Carboxylic Acid "
    "Polycondensation(Polyesterification)"
)

DIOL_DIHALIDE = (
    "Diol and Di-Acid Halide "
    "Polycondensation(Polyesterification)"
)


EXPECTED_REACTIONS = {
    HYDROXY_ACID_HOMO,
    HYDROXY_ACID_CO,
    HYDROXY_HALIDE_HOMO,
    HYDROXY_HALIDE_CO,
    DIOL_DIACID,
    DIOL_DIHALIDE,
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


# =============================================================================
# Metadata
# =============================================================================

def test_hydroxy_acid_homo_metadata():
    reaction = REACTIONS[
        HYDROXY_ACID_HOMO
    ]

    assert reaction["same_reactants"] is True

    assert (
        reaction["reactant_1"]
        == "hydroxy_carboxylic_acid"
    )

    assert reaction.get("reactant_2") is None
    assert reaction["product"] == "polyester_chain"
    assert reaction["delete_atom"] is True


def test_hydroxy_acid_co_metadata():
    reaction = REACTIONS[
        HYDROXY_ACID_CO
    ]

    assert reaction["same_reactants"] is False

    assert (
        reaction["reactant_1"]
        == "hydroxy_carboxylic_acid"
    )

    assert (
        reaction["reactant_2"]
        == "hydroxy_carboxylic_acid"
    )


def test_hydroxy_halide_homo_metadata():
    reaction = REACTIONS[
        HYDROXY_HALIDE_HOMO
    ]

    assert reaction["same_reactants"] is True

    assert (
        reaction["reactant_1"]
        == "hydroxy_acid_halide"
    )


def test_hydroxy_halide_co_metadata():
    reaction = REACTIONS[
        HYDROXY_HALIDE_CO
    ]

    assert reaction["same_reactants"] is False

    assert (
        reaction["reactant_1"]
        == "hydroxy_acid_halide"
    )

    assert (
        reaction["reactant_2"]
        == "hydroxy_acid_halide"
    )


def test_diol_diacid_metadata():
    reaction = REACTIONS[DIOL_DIACID]

    assert reaction["same_reactants"] is False
    assert reaction["reactant_1"] == "diol"

    assert (
        reaction["reactant_2"]
        == "di_carboxylic_acid"
    )

    assert reaction["product"] == "polyester_chain"


def test_diol_diacid_halide_metadata():
    reaction = REACTIONS[DIOL_DIHALIDE]

    assert reaction["same_reactants"] is False
    assert reaction["reactant_1"] == "diol"

    assert (
        reaction["reactant_2"]
        == "di_carboxylic_acid_halide"
    )


# =============================================================================
# Initiator bond
# =============================================================================

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


# =============================================================================
# Hydroxy carboxylic acid
# =============================================================================

def test_hydroxy_acid_polycondensation_executes():
    products = run_reaction(
        HYDROXY_ACID_HOMO,
        "CC(O)C(=O)O",
        "CC(O)C(=O)O",
    )

    assert len(products) > 0


def test_hydroxy_acid_polycondensation_has_byproduct():
    products = run_reaction(
        HYDROXY_ACID_HOMO,
        "CC(O)C(=O)O",
        "CC(O)C(=O)O",
    )

    assert len(products[0]) == 2


def test_hydroxy_acid_products_are_sanitizable():
    products = run_reaction(
        HYDROXY_ACID_HOMO,
        "CC(O)C(=O)O",
        "CC(O)C(=O)O",
    )

    for product in products[0]:
        Chem.SanitizeMol(product)


# =============================================================================
# Hydroxy acid halide
# =============================================================================

def test_hydroxy_acid_halide_polycondensation_executes():
    products = run_reaction(
        HYDROXY_HALIDE_HOMO,
        "OCCC(=O)Cl",
        "OCCC(=O)Cl",
    )

    assert len(products) > 0


def test_hydroxy_acid_halide_has_byproduct():
    products = run_reaction(
        HYDROXY_HALIDE_HOMO,
        "OCCC(=O)Cl",
        "OCCC(=O)Cl",
    )

    assert len(products[0]) == 2


# =============================================================================
# Diol + diacid
# =============================================================================

def test_diol_diacid_polycondensation_executes():
    products = run_reaction(
        DIOL_DIACID,
        "O=C(O)CCC(=O)O",
        "OCCO",
    )

    assert len(products) > 0


def test_diol_diacid_products_are_sanitizable():
    products = run_reaction(
        DIOL_DIACID,
        "O=C(O)CCC(=O)O",
        "OCCO",
    )

    for product in products[0]:
        Chem.SanitizeMol(product)


def test_diol_diacid_has_byproduct():
    products = run_reaction(
        DIOL_DIACID,
        "O=C(O)CCC(=O)O",
        "OCCO",
    )

    assert len(products[0]) == 2


# =============================================================================
# Diol + diacid halide
# =============================================================================

def test_diol_diacid_halide_polycondensation_executes():
    products = run_reaction(
        DIOL_DIHALIDE,
        "O=C(Cl)CCC(=O)Cl",
        "OCCO",
    )

    assert len(products) > 0


def test_diol_diacid_halide_products_are_sanitizable():
    products = run_reaction(
        DIOL_DIHALIDE,
        "O=C(Cl)CCC(=O)Cl",
        "OCCO",
    )

    for product in products[0]:
        Chem.SanitizeMol(product)


def test_diol_diacid_halide_has_byproduct():
    products = run_reaction(
        DIOL_DIHALIDE,
        "O=C(Cl)CCC(=O)Cl",
        "OCCO",
    )

    assert len(products[0]) == 2


# =============================================================================
# Negative tests
# =============================================================================

def test_hydroxy_acid_reaction_rejects_molecule_without_hydroxyl():
    products = run_reaction(
        HYDROXY_ACID_HOMO,
        "CC(=O)O",
        "CC(=O)O",
    )

    assert len(products) == 0


def test_diol_diacid_reaction_rejects_hydrocarbon():
    products = run_reaction(
        DIOL_DIACID,
        "O=C(O)CCC(=O)O",
        "CCCC",
    )

    assert len(products) == 0