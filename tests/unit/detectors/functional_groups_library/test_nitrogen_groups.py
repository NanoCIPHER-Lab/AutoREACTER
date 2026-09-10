import pytest
from rdkit import Chem

from AutoREACTER.detectors.functional_groups_detector import (
    FunctionalGroupsDetector,
)
from AutoREACTER.detectors.functional_groups_library.nitrogen_groups import (
    FUNCTIONAL_GROUPS,
)


EXPECTED_GROUPS = {
    "primary_amine_monomer",
    "secondary_amine_monomer",
    "di_amine_monomer",
    "di_primary_amine_monomer",
}


def detect(group_key, smiles):
    detector = FunctionalGroupsDetector()
    group = FUNCTIONAL_GROUPS[group_key]

    return detector.detect_monomer_functionality(
        Chem.MolFromSmiles(smiles),
        group["functionality_type"],
        group["smarts_1"],
    )


# ============================================================================
# Structure / metadata
# ============================================================================

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
    pattern = Chem.MolFromSmarts(
        FUNCTIONAL_GROUPS[group_key]["smarts_1"]
    )

    assert pattern is not None


def test_primary_amine_metadata():
    group = FUNCTIONAL_GROUPS["primary_amine_monomer"]

    assert group["functionality_type"] == "mono"
    assert group["group_name"] == "primary_amine"
    assert group["comments"] is None


def test_secondary_amine_metadata():
    group = FUNCTIONAL_GROUPS["secondary_amine_monomer"]

    assert group["functionality_type"] == "mono"
    assert group["group_name"] == "secondary_amine"
    assert group["comments"] is None


def test_di_amine_metadata():
    group = FUNCTIONAL_GROUPS["di_amine_monomer"]

    assert group["functionality_type"] == "di_identical"
    assert group["group_name"] == "di_amine"


def test_di_primary_amine_metadata():
    group = FUNCTIONAL_GROUPS["di_primary_amine_monomer"]

    assert group["functionality_type"] == "di_identical"
    assert group["group_name"] == "di_primary_amine"


# ============================================================================
# Primary amine
# ============================================================================

def test_primary_amine_matches_ethylamine():
    functionality_count, count_1, count_2, matches = detect(
        "primary_amine_monomer",
        "CCN",
    )

    assert functionality_count == 1
    assert count_1 == 1
    assert count_2 is None
    assert len(matches) == 1


def test_primary_amine_matches_glycine_nitrogen():
    functionality_count, count_1, _, matches = detect(
        "primary_amine_monomer",
        "NCC(=O)O",
    )

    assert functionality_count == 1
    assert count_1 == 1
    assert len(matches) == 1


def test_primary_amine_does_not_match_secondary_amine():
    result = detect(
        "primary_amine_monomer",
        "CNC",
    )

    assert result == (0, 0, None, ())


def test_primary_amine_does_not_match_amide():
    """
    The SMARTS explicitly excludes N-C(=O) environments.
    """
    result = detect(
        "primary_amine_monomer",
        "CC(=O)N",
    )

    assert result == (0, 0, None, ())


# ============================================================================
# Secondary amine
# ============================================================================

def test_secondary_amine_matches_dimethylamine():
    functionality_count, count_1, count_2, matches = detect(
        "secondary_amine_monomer",
        "CNC",
    )

    assert functionality_count == 1
    assert count_1 == 1
    assert count_2 is None
    assert len(matches) == 1


def test_secondary_amine_does_not_match_primary_amine():
    result = detect(
        "secondary_amine_monomer",
        "CCN",
    )

    assert result == (0, 0, None, ())


def test_secondary_amine_does_not_match_tertiary_amine():
    result = detect(
        "secondary_amine_monomer",
        "CN(C)C",
    )

    assert result == (0, 0, None, ())


def test_secondary_amine_does_not_match_secondary_amide():
    result = detect(
        "secondary_amine_monomer",
        "CC(=O)NC",
    )

    assert result == (0, 0, None, ())


# ============================================================================
# Di-amine
# ============================================================================

def test_ethylenediamine_qualifies_as_diamine():
    functionality_count, count_1, count_2, matches = detect(
        "di_amine_monomer",
        "NCCN",
    )

    assert functionality_count == 2
    assert count_1 == 2
    assert count_2 is None
    assert len(matches) == 2


def test_hexamethylenediamine_qualifies_as_diamine():
    functionality_count, count_1, _, matches = detect(
        "di_amine_monomer",
        "NCCCCCCN",
    )

    assert functionality_count == 2
    assert count_1 == 2
    assert len(matches) == 2


def test_monoamine_is_rejected_as_diamine():
    functionality_count, count_1, count_2, matches = detect(
        "di_amine_monomer",
        "CCN",
    )

    assert functionality_count == 0
    assert count_1 == 1
    assert count_2 is None
    assert len(matches) == 1


def test_mixed_primary_secondary_diamine_qualifies_as_diamine():
    """
    di_amine permits both H2 and H1 nitrogens.

    N-methylethylenediamine has one primary and one secondary amine,
    so both sites should be counted.
    """
    functionality_count, count_1, _, matches = detect(
        "di_amine_monomer",
        "NCCNC",
    )

    assert functionality_count == 2
    assert count_1 == 2
    assert len(matches) == 2


# ============================================================================
# Di-primary amine
# ============================================================================

def test_ethylenediamine_qualifies_as_di_primary_amine():
    functionality_count, count_1, count_2, matches = detect(
        "di_primary_amine_monomer",
        "NCCN",
    )

    assert functionality_count == 2
    assert count_1 == 2
    assert count_2 is None
    assert len(matches) == 2


def test_hexamethylenediamine_qualifies_as_di_primary_amine():
    functionality_count, count_1, _, matches = detect(
        "di_primary_amine_monomer",
        "NCCCCCCN",
    )

    assert functionality_count == 2
    assert count_1 == 2
    assert len(matches) == 2


def test_single_primary_amine_is_rejected_as_di_primary():
    functionality_count, count_1, _, matches = detect(
        "di_primary_amine_monomer",
        "CCN",
    )

    assert functionality_count == 0
    assert count_1 == 1
    assert len(matches) == 1


def test_primary_secondary_pair_is_not_di_primary():
    """
    N-methylethylenediamine contains only one primary amine.

    The secondary nitrogen must not count toward di_primary_amine.
    """
    functionality_count, count_1, _, matches = detect(
        "di_primary_amine_monomer",
        "NCCNC",
    )

    assert functionality_count == 0
    assert count_1 == 1
    assert len(matches) == 1


def test_amide_nitrogen_does_not_count_toward_di_primary():
    functionality_count, count_1, _, matches = detect(
        "di_primary_amine_monomer",
        "NCCC(=O)N",
    )

    assert functionality_count == 0
    assert count_1 == 1
    assert len(matches) == 1