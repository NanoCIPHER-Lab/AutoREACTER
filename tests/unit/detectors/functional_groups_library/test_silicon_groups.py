import pytest
from rdkit import Chem

from AutoREACTER.detectors.functional_groups_detector import (
    FunctionalGroupsDetector,
)
from AutoREACTER.detectors.functional_groups_library.silicon_groups import (
    FUNCTIONAL_GROUPS,
)


EXPECTED_GROUPS = {
    "dichlorosilane_monomer",
    "silanediol_monomer",
}


def detect(group_key, smiles):
    detector = FunctionalGroupsDetector()
    group = FUNCTIONAL_GROUPS[group_key]

    return detector.detect_monomer_functionality(
        Chem.MolFromSmiles(smiles),
        group["functionality_type"],
        group["smarts_1"],
    )


def test_expected_active_groups():
    assert set(FUNCTIONAL_GROUPS) == EXPECTED_GROUPS


@pytest.mark.parametrize("group_key", sorted(EXPECTED_GROUPS))
def test_required_fields_exist(group_key):
    group = FUNCTIONAL_GROUPS[group_key]

    assert "functionality_type" in group
    assert "smarts_1" in group
    assert "group_name" in group
    assert "comments" in group


@pytest.mark.parametrize("group_key", sorted(EXPECTED_GROUPS))
def test_all_smarts_compile(group_key):
    assert Chem.MolFromSmarts(
        FUNCTIONAL_GROUPS[group_key]["smarts_1"]
    ) is not None


def test_dichlorosilane_metadata():
    group = FUNCTIONAL_GROUPS["dichlorosilane_monomer"]

    assert group["functionality_type"] == "di_identical"
    assert group["group_name"] == "dichlorosilane"
    assert group["comments"] is None


def test_silanediol_metadata():
    group = FUNCTIONAL_GROUPS["silanediol_monomer"]

    assert group["functionality_type"] == "di_identical"
    assert group["group_name"] == "silanediol"
    assert group["comments"] is None


def test_dichlorosilane_has_two_si_cl_sites():
    functionality_count, count_1, count_2, matches = detect(
        "dichlorosilane_monomer",
        "Cl[Si](C)(C)Cl",
    )

    assert functionality_count == 2
    assert count_1 == 2
    assert count_2 is None
    assert len(matches) == 2


def test_monochlorosilane_is_rejected():
    functionality_count, count_1, count_2, matches = detect(
        "dichlorosilane_monomer",
        "C[Si](C)(C)Cl",
    )

    assert functionality_count == 0
    assert count_1 == 1
    assert count_2 is None
    assert len(matches) == 1


def test_nonchlorinated_silane_does_not_match():
    result = detect(
        "dichlorosilane_monomer",
        "C[Si](C)(C)C",
    )

    assert result == (0, 0, None, ())


def test_silanediol_has_two_si_oh_sites():
    functionality_count, count_1, count_2, matches = detect(
        "silanediol_monomer",
        "O[Si](C)(C)O",
    )

    assert functionality_count == 2
    assert count_1 == 2
    assert count_2 is None
    assert len(matches) == 2


def test_single_silanol_is_rejected_as_silanediol():
    functionality_count, count_1, count_2, matches = detect(
        "silanediol_monomer",
        "C[Si](C)(C)O",
    )

    assert functionality_count == 0
    assert count_1 == 1
    assert count_2 is None
    assert len(matches) == 1


def test_regular_alcohol_does_not_match_silanediol():
    result = detect(
        "silanediol_monomer",
        "CCO",
    )

    assert result == (0, 0, None, ())