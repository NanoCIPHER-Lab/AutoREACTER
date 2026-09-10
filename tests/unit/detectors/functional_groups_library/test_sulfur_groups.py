from rdkit import Chem

from AutoREACTER.detectors.functional_groups_detector import (
    FunctionalGroupsDetector,
)
from AutoREACTER.detectors.functional_groups_library.sulfur_groups import (
    FUNCTIONAL_GROUPS,
)


GROUP_KEY = "dithiol_monomer"


def detect(smiles):
    detector = FunctionalGroupsDetector()
    group = FUNCTIONAL_GROUPS[GROUP_KEY]

    return detector.detect_monomer_functionality(
        Chem.MolFromSmiles(smiles),
        group["functionality_type"],
        group["smarts_1"],
    )


def test_expected_active_groups():
    assert set(FUNCTIONAL_GROUPS) == {
        "dithiol_monomer",
    }


def test_dithiol_metadata():
    group = FUNCTIONAL_GROUPS[GROUP_KEY]

    assert group["functionality_type"] == "di_identical"
    assert group["group_name"] == "dithiol"
    assert group["comments"] is None


def test_dithiol_smarts_compiles():
    assert Chem.MolFromSmarts(
        FUNCTIONAL_GROUPS[GROUP_KEY]["smarts_1"]
    ) is not None


def test_dithiol_has_two_thiol_sites():
    functionality_count, count_1, count_2, matches = detect(
        "SCCS"
    )

    assert functionality_count == 2
    assert count_1 == 2
    assert count_2 is None
    assert len(matches) == 2


def test_monothiol_is_rejected_as_dithiol():
    functionality_count, count_1, count_2, matches = detect(
        "CCS"
    )

    assert functionality_count == 0
    assert count_1 == 1
    assert count_2 is None
    assert len(matches) == 1


def test_thioether_does_not_match_thiol():
    result = detect(
        "CSC"
    )

    assert result == (0, 0, None, ())


def test_disulfide_does_not_match_thiol():
    result = detect(
        "CSSC"
    )

    assert result == (0, 0, None, ())


def test_carbonyl_attached_sulfur_is_excluded():
    result = detect(
        "CC(=O)S"
    )

    assert result == (0, 0, None, ())