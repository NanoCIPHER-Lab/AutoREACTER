import pytest
from rdkit import Chem

from AutoREACTER.detectors.functional_groups_library.oxygen_groups import (
    FUNCTIONAL_GROUPS,
)


# ---------------------------------------------------------------------------
# Basic library structure
# ---------------------------------------------------------------------------

def test_functional_groups_is_dict():
    """oxygen_groups should expose FUNCTIONAL_GROUPS as a dictionary."""
    assert isinstance(FUNCTIONAL_GROUPS, dict)


def test_expected_active_groups_are_present():
    """The currently active oxygen functional groups should be registered."""
    assert set(FUNCTIONAL_GROUPS) == {
        "diol_monomer",
        "water_monomer",
    }


@pytest.mark.parametrize(
    "group_key",
    [
        "diol_monomer",
        "water_monomer",
    ],
)
def test_each_group_has_required_fields(group_key):
    """Every oxygen functional-group definition should contain required fields."""
    entry = FUNCTIONAL_GROUPS[group_key]

    assert "functionality_type" in entry
    assert "smarts_1" in entry
    assert "group_name" in entry
    assert "comments" in entry


def test_diol_metadata():
    """The diol definition should expose the intended metadata."""
    entry = FUNCTIONAL_GROUPS["diol_monomer"]

    assert entry["functionality_type"] == "di_identical"
    assert entry["group_name"] == "diol"
    assert entry["comments"] is None


def test_water_metadata():
    """The water definition should expose the intended metadata."""
    entry = FUNCTIONAL_GROUPS["water_monomer"]

    assert entry["functionality_type"] == "mono"
    assert entry["group_name"] == "water"
    assert entry["comments"] is None


# ---------------------------------------------------------------------------
# SMARTS validity
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "group_key",
    [
        "diol_monomer",
        "water_monomer",
    ],
)
def test_smarts_compiles(group_key):
    """Every active oxygen SMARTS definition must compile with RDKit."""
    smarts = FUNCTIONAL_GROUPS[group_key]["smarts_1"]

    pattern = Chem.MolFromSmarts(smarts)

    assert pattern is not None


# ---------------------------------------------------------------------------
# Diol hydroxyl SMARTS
# ---------------------------------------------------------------------------

def test_diol_pattern_matches_two_hydroxyl_sites_in_ethylene_glycol():
    """
    The diol SMARTS describes an eligible hydroxyl site.

    Ethylene glycol contains two such hydroxyl groups, so two matches
    should be found.
    """
    molecule = Chem.MolFromSmiles("OCCO")
    pattern = Chem.MolFromSmarts(
        FUNCTIONAL_GROUPS["diol_monomer"]["smarts_1"]
    )

    matches = molecule.GetSubstructMatches(pattern)

    assert len(matches) == 2


def test_diol_pattern_matches_one_site_in_ethanol():
    """
    A monoalcohol contains one matching hydroxyl site.

    Whether a molecule qualifies as a diol based on the number of detected
    sites belongs to detector logic, not the SMARTS library itself.
    """
    molecule = Chem.MolFromSmiles("CCO")
    pattern = Chem.MolFromSmarts(
        FUNCTIONAL_GROUPS["diol_monomer"]["smarts_1"]
    )

    matches = molecule.GetSubstructMatches(pattern)

    assert len(matches) == 1


def test_diol_pattern_does_not_match_ether_oxygen():
    """An ether oxygen should not match the hydroxyl SMARTS."""
    molecule = Chem.MolFromSmiles("COC")
    pattern = Chem.MolFromSmarts(
        FUNCTIONAL_GROUPS["diol_monomer"]["smarts_1"]
    )

    assert not molecule.HasSubstructMatch(pattern)


def test_diol_pattern_excludes_carboxylic_acid_hydroxyl():
    """
    The SMARTS explicitly excludes oxygen attached to a carbonyl center.

    Therefore the hydroxyl oxygen of acetic acid should not be detected.
    """
    molecule = Chem.MolFromSmiles("CC(=O)O")
    pattern = Chem.MolFromSmarts(
        FUNCTIONAL_GROUPS["diol_monomer"]["smarts_1"]
    )

    assert not molecule.HasSubstructMatch(pattern)


# ---------------------------------------------------------------------------
# Water SMARTS
# ---------------------------------------------------------------------------

def test_water_pattern_matches_water():
    """The water SMARTS should match a water molecule."""
    molecule = Chem.MolFromSmiles("O")
    pattern = Chem.MolFromSmarts(
        FUNCTIONAL_GROUPS["water_monomer"]["smarts_1"]
    )

    assert molecule.HasSubstructMatch(pattern)


@pytest.mark.parametrize(
    "smiles",
    [
        "CO",       # methanol
        "CCO",      # ethanol
        "OCCO",     # ethylene glycol
        "COC",      # dimethyl ether
        "CC(=O)O",  # acetic acid
    ],
)
def test_water_pattern_does_not_match_non_water_oxygen_compounds(smiles):
    """The water SMARTS should not identify ordinary oxygen compounds as water."""
    molecule = Chem.MolFromSmiles(smiles)
    pattern = Chem.MolFromSmarts(
        FUNCTIONAL_GROUPS["water_monomer"]["smarts_1"]
    )

    assert not molecule.HasSubstructMatch(pattern)
