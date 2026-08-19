import pytest
from rdkit import Chem

from AutoREACTER.detectors.functional_groups_detector import (
    FunctionalGroupsDetector,
)
from AutoREACTER.detectors.functional_groups_library.oxygen_groups import (
    FUNCTIONAL_GROUPS,
)


# ============================================================================
# Library structure
# ============================================================================

def test_functional_groups_is_dictionary():
    assert isinstance(FUNCTIONAL_GROUPS, dict)


def test_expected_active_oxygen_groups():
    assert set(FUNCTIONAL_GROUPS.keys()) == {
        "diol_monomer",
        "water_monomer",
    }


@pytest.mark.parametrize(
    "group_name",
    [
        "diol_monomer",
        "water_monomer",
    ],
)
def test_required_fields_exist(group_name):
    group = FUNCTIONAL_GROUPS[group_name]

    assert "functionality_type" in group
    assert "smarts_1" in group
    assert "group_name" in group
    assert "comments" in group


# ============================================================================
# Metadata
# ============================================================================

def test_diol_metadata():
    group = FUNCTIONAL_GROUPS["diol_monomer"]

    assert group["functionality_type"] == "di_identical"
    assert group["group_name"] == "diol"
    assert group["comments"] is None


def test_water_metadata():
    group = FUNCTIONAL_GROUPS["water_monomer"]

    assert group["functionality_type"] == "mono"
    assert group["group_name"] == "water"
    assert group["comments"] is None


# ============================================================================
# SMARTS validity
# ============================================================================

@pytest.mark.parametrize(
    "group_name",
    [
        "diol_monomer",
        "water_monomer",
    ],
)
def test_smarts_is_valid_rdkit_pattern(group_name):
    smarts = FUNCTIONAL_GROUPS[group_name]["smarts_1"]

    pattern = Chem.MolFromSmarts(smarts)

    assert pattern is not None


# ============================================================================
# Diol SMARTS behavior
# ============================================================================

def test_diol_smarts_finds_two_sites_in_ethylene_glycol():
    group = FUNCTIONAL_GROUPS["diol_monomer"]

    mol = Chem.MolFromSmiles("OCCO")
    pattern = Chem.MolFromSmarts(group["smarts_1"])

    matches = mol.GetSubstructMatches(pattern)

    assert len(matches) == 2


def test_diol_smarts_finds_three_sites_in_glycerol():
    group = FUNCTIONAL_GROUPS["diol_monomer"]

    mol = Chem.MolFromSmiles("OCC(O)CO")
    pattern = Chem.MolFromSmarts(group["smarts_1"])

    matches = mol.GetSubstructMatches(pattern)

    assert len(matches) == 3


def test_diol_smarts_finds_only_one_site_in_ethanol():
    """
    The SMARTS identifies hydroxyl sites.

    Ethanol contains one valid hydroxyl site, but the detector should later
    reject it as a di_identical monomer because at least two sites are required.
    """
    group = FUNCTIONAL_GROUPS["diol_monomer"]

    mol = Chem.MolFromSmiles("CCO")
    pattern = Chem.MolFromSmarts(group["smarts_1"])

    matches = mol.GetSubstructMatches(pattern)

    assert len(matches) == 1


def test_diol_smarts_does_not_match_ether():
    group = FUNCTIONAL_GROUPS["diol_monomer"]

    mol = Chem.MolFromSmiles("COC")
    pattern = Chem.MolFromSmarts(group["smarts_1"])

    matches = mol.GetSubstructMatches(pattern)

    assert len(matches) == 0


def test_diol_smarts_excludes_carboxylic_acid_oh():
    """
    The diol SMARTS explicitly excludes oxygen bonded to a carbonyl carbon.
    """
    group = FUNCTIONAL_GROUPS["diol_monomer"]

    mol = Chem.MolFromSmiles("CC(=O)O")
    pattern = Chem.MolFromSmarts(group["smarts_1"])

    matches = mol.GetSubstructMatches(pattern)

    assert len(matches) == 0


# ============================================================================
# Diol behavior through FunctionalGroupsDetector
# ============================================================================

def test_ethylene_glycol_qualifies_as_diol():
    detector = FunctionalGroupsDetector()
    group = FUNCTIONAL_GROUPS["diol_monomer"]

    mol = Chem.MolFromSmiles("OCCO")

    functionality_count, count_1, count_2, matches = (
        detector.detect_monomer_functionality(
            mol,
            group["functionality_type"],
            group["smarts_1"],
        )
    )

    assert functionality_count == 2
    assert count_1 == 2
    assert count_2 is None
    assert len(matches) == 2


def test_ethanol_does_not_qualify_as_diol():
    detector = FunctionalGroupsDetector()
    group = FUNCTIONAL_GROUPS["diol_monomer"]

    mol = Chem.MolFromSmiles("CCO")

    functionality_count, count_1, count_2, matches = (
        detector.detect_monomer_functionality(
            mol,
            group["functionality_type"],
            group["smarts_1"],
        )
    )

    assert functionality_count == 0
    assert count_1 == 1
    assert count_2 is None
    assert len(matches) == 1


def test_glycerol_qualifies_and_preserves_three_site_count():
    detector = FunctionalGroupsDetector()
    group = FUNCTIONAL_GROUPS["diol_monomer"]

    mol = Chem.MolFromSmiles("OCC(O)CO")

    functionality_count, count_1, count_2, matches = (
        detector.detect_monomer_functionality(
            mol,
            group["functionality_type"],
            group["smarts_1"],
        )
    )

    assert functionality_count == 2
    assert count_1 == 3
    assert count_2 is None
    assert len(matches) == 3


# ============================================================================
# Water SMARTS behavior
# ============================================================================

def test_water_smarts_matches_water():
    group = FUNCTIONAL_GROUPS["water_monomer"]

    mol = Chem.MolFromSmiles("O")
    pattern = Chem.MolFromSmarts(group["smarts_1"])

    matches = mol.GetSubstructMatches(pattern)

    assert len(matches) == 1


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
def test_water_smarts_does_not_match_non_water_oxygen_compounds(smiles):
    group = FUNCTIONAL_GROUPS["water_monomer"]

    mol = Chem.MolFromSmiles(smiles)
    pattern = Chem.MolFromSmarts(group["smarts_1"])

    matches = mol.GetSubstructMatches(pattern)

    assert len(matches) == 0


# ============================================================================
# Water behavior through FunctionalGroupsDetector
# ============================================================================

def test_water_qualifies_as_mono_functional():
    detector = FunctionalGroupsDetector()
    group = FUNCTIONAL_GROUPS["water_monomer"]

    mol = Chem.MolFromSmiles("O")

    functionality_count, count_1, count_2, matches = (
        detector.detect_monomer_functionality(
            mol,
            group["functionality_type"],
            group["smarts_1"],
        )
    )

    assert functionality_count == 1
    assert count_1 == 1
    assert count_2 is None
    assert len(matches) == 1


def test_ethanol_does_not_match_water_functionality():
    detector = FunctionalGroupsDetector()
    group = FUNCTIONAL_GROUPS["water_monomer"]

    mol = Chem.MolFromSmiles("CCO")

    functionality_count, count_1, count_2, matches = (
        detector.detect_monomer_functionality(
            mol,
            group["functionality_type"],
            group["smarts_1"],
        )
    )

    assert functionality_count == 0
    assert count_1 == 0
    assert count_2 is None
    assert matches == ()