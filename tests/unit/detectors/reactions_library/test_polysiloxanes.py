import pytest
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from AutoREACTER.detectors.reactions_library import registry
from AutoREACTER.detectors.reactions_library.polysiloxanes import REACTIONS


HYDROLYSIS = "Dichlorosilane Hydrolysis to Silanol"
HOMO = "Silanediol Polycondensation(Polysiloxane Formation)"
CO = "Silanediol and Silanediol Copolycondensation(Polysiloxane Formation)"

EXPECTED_REACTIONS = {
    HYDROLYSIS,
    HOMO,
    CO,
}


def run_reaction(name, smiles_1, smiles_2):
    rxn = rdChemReactions.ReactionFromSmarts(
        REACTIONS[name]["reaction"]
    )
    assert rxn is not None

    mol_1 = Chem.AddHs(Chem.MolFromSmiles(smiles_1))
    mol_2 = Chem.AddHs(Chem.MolFromSmiles(smiles_2))

    return rxn.RunReactants((mol_1, mol_2))


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


def test_hydrolysis_metadata():
    reaction = REACTIONS[HYDROLYSIS]

    assert reaction["same_reactants"] is False
    assert reaction["reactant_1"] == "dichlorosilane"
    assert reaction["reactant_2"] == "water"
    assert reaction["product"] == "silanediol"
    assert reaction["delete_atom"] is True


def test_silanediol_homopolycondensation_metadata():
    reaction = REACTIONS[HOMO]

    assert reaction["same_reactants"] is True
    assert reaction["reactant_1"] == "silanediol"
    assert reaction.get("reactant_2") is None
    assert reaction["product"] == "polysiloxane_chain"
    assert reaction["delete_atom"] is True


def test_silanediol_copolycondensation_metadata():
    reaction = REACTIONS[CO]

    assert reaction["same_reactants"] is False
    assert reaction["reactant_1"] == "silanediol"
    assert reaction["reactant_2"] == "silanediol"


def test_dichlorosilane_hydrolysis_executes():
    products = run_reaction(
        HYDROLYSIS,
        "Cl[Si](C)(C)Cl",
        "O",
    )

    assert len(products) > 0


def test_dichlorosilane_hydrolysis_has_two_products():
    products = run_reaction(
        HYDROLYSIS,
        "Cl[Si](C)(C)Cl",
        "O",
    )

    assert len(products[0]) == 2


def test_silanediol_condensation_executes():
    products = run_reaction(
        HOMO,
        "O[Si](C)(C)O",
        "O[Si](C)(C)O",
    )

    assert len(products) > 0


def test_silanediol_condensation_has_water_fragment():
    products = run_reaction(
        HOMO,
        "O[Si](C)(C)O",
        "O[Si](C)(C)O",
    )

    assert len(products[0]) == 2


def test_hydrolysis_rejects_non_water_second_reactant():
    products = run_reaction(
        HYDROLYSIS,
        "Cl[Si](C)(C)Cl",
        "CCO",
    )

    assert len(products) == 0