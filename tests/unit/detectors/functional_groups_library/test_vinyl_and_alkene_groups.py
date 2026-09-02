import pytest
from rdkit import Chem

from AutoREACTER.detectors.functional_groups_detector import (
    FunctionalGroupsDetector,
)
from AutoREACTER.detectors.functional_groups_library.vinyl_and_alkene_groups import (
    FUNCTIONAL_GROUPS,
)


EXPECTED_GROUPS = {
    "vinyl_monomer",
    "diene_monomer",
    "tetrafluoroethylene_monomer",
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
    assert Chem.MolFromSmarts(
        FUNCTIONAL_GROUPS[group_key]["smarts_1"]
    ) is not None


def test_vinyl_metadata():
    group = FUNCTIONAL_GROUPS["vinyl_monomer"]

    assert group["functionality_type"] == "vinyl"
    assert group["group_name"] == "vinyl"


def test_diene_metadata():
    group = FUNCTIONAL_GROUPS["diene_monomer"]

    assert group["functionality_type"] == "di_identical"
    assert group["group_name"] == "diene"


def test_tetrafluoroethylene_metadata():
    group = FUNCTIONAL_GROUPS["tetrafluoroethylene_monomer"]

    assert group["functionality_type"] == "vinyl"
    assert group["group_name"] == "tetrafluoroethylene"


# ============================================================================
# Vinyl
# ============================================================================

@pytest.mark.parametrize(
    "smiles",
    [
        "C=C",          # ethylene
        "C=CC",         # propene
        "C=Cc1ccccc1",  # styrene
        "C=C(C)C",      # substituted terminal alkene
    ],
)
def test_terminal_vinyl_groups_qualify(smiles):
    functionality_count, count_1, count_2, matches = detect(
        "vinyl_monomer",
        smiles,
    )

    assert functionality_count == 1
    assert count_1 >= 1
    assert count_2 is None
    assert len(matches) == count_1


def test_internal_alkene_does_not_match_terminal_vinyl():
    result = detect(
        "vinyl_monomer",
        "CC=CC",
    )

    assert result == (0, 0, None, ())


def test_saturated_molecule_does_not_match_vinyl():
    result = detect(
        "vinyl_monomer",
        "CCCC",
    )

    assert result == (0, 0, None, ())


def test_ring_double_bond_does_not_match_vinyl():
    result = detect(
        "vinyl_monomer",
        "C1=CCCCC1",
    )

    assert result == (0, 0, None, ())


# ============================================================================
# Diene
# ============================================================================

def test_butadiene_has_two_terminal_vinyl_sites():
    functionality_count, count_1, count_2, matches = detect(
        "diene_monomer",
        "C=CC=C",
    )

    assert functionality_count == 2
    assert count_1 == 2
    assert count_2 is None
    assert len(matches) == 2


def test_single_terminal_alkene_is_rejected_as_diene():
    functionality_count, count_1, count_2, matches = detect(
        "diene_monomer",
        "C=CCC",
    )

    assert functionality_count == 0
    assert count_1 == 1
    assert count_2 is None
    assert len(matches) == 1


def test_internal_diene_without_terminal_ch2_is_rejected():
    result = detect(
        "diene_monomer",
        "CC=CC=CC",
    )

    assert result == (0, 0, None, ())


def test_molecule_with_more_than_two_terminal_vinyl_sites_preserves_count():
    functionality_count, count_1, count_2, matches = detect(
        "diene_monomer",
        "C=C(C=C)C=C",
    )

    assert functionality_count == 2
    assert count_1 >= 2
    assert count_2 is None
    assert len(matches) == count_1


# ============================================================================
# Tetrafluoroethylene
# ============================================================================

def test_tetrafluoroethylene_qualifies():
    functionality_count, count_1, count_2, matches = detect(
        "tetrafluoroethylene_monomer",
        "FC(F)=C(F)F",
    )

    assert functionality_count == 1
    assert count_1 == 1
    assert count_2 is None
    assert len(matches) == 1


@pytest.mark.parametrize(
    "smiles",
    [
        "C=C",
        "C=CC",
        "C=Cc1ccccc1",
        "FC(F)=CC",
    ],
)
def test_non_tetrafluoroethylene_alkenes_do_not_match(smiles):
    result = detect(
        "tetrafluoroethylene_monomer",
        smiles,
    )

    assert result == (0, 0, None, ())