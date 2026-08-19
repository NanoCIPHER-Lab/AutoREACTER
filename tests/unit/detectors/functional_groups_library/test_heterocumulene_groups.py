from rdkit import Chem

from AutoREACTER.detectors.functional_groups_detector import (
    FunctionalGroupsDetector,
)
from AutoREACTER.detectors.functional_groups_library.heterocumulene_groups import (
    FUNCTIONAL_GROUPS,
)


def detect(smiles):
    detector = FunctionalGroupsDetector()
    group = FUNCTIONAL_GROUPS["di_isocyanate_monomer"]

    return detector.detect_monomer_functionality(
        Chem.MolFromSmiles(smiles),
        group["functionality_type"],
        group["smarts_1"],
    )


# ============================================================================
# Structure / metadata
# ============================================================================

def test_expected_active_groups():
    assert set(FUNCTIONAL_GROUPS) == {
        "di_isocyanate_monomer",
    }


def test_di_isocyanate_metadata():
    group = FUNCTIONAL_GROUPS["di_isocyanate_monomer"]

    assert group["functionality_type"] == "di_identical"
    assert group["group_name"] == "di_isocyanate"
    assert group["comments"] is None


def test_di_isocyanate_smarts_compiles():
    group = FUNCTIONAL_GROUPS["di_isocyanate_monomer"]

    pattern = Chem.MolFromSmarts(group["smarts_1"])

    assert pattern is not None


# ============================================================================
# SMARTS behavior
# ============================================================================

def test_diisocyanate_contains_two_isocyanate_sites():
    group = FUNCTIONAL_GROUPS["di_isocyanate_monomer"]

    mol = Chem.MolFromSmiles("O=C=NCCCCCCN=C=O")
    pattern = Chem.MolFromSmarts(group["smarts_1"])

    matches = mol.GetSubstructMatches(pattern)

    assert len(matches) == 2


def test_monoisocyanate_contains_one_isocyanate_site():
    group = FUNCTIONAL_GROUPS["di_isocyanate_monomer"]

    mol = Chem.MolFromSmiles("CN=C=O")
    pattern = Chem.MolFromSmarts(group["smarts_1"])

    matches = mol.GetSubstructMatches(pattern)

    assert len(matches) == 1


def test_amine_does_not_match_isocyanate():
    group = FUNCTIONAL_GROUPS["di_isocyanate_monomer"]

    mol = Chem.MolFromSmiles("NCCCCCCN")
    pattern = Chem.MolFromSmarts(group["smarts_1"])

    assert mol.GetSubstructMatches(pattern) == ()


def test_carbon_dioxide_does_not_match_isocyanate():
    group = FUNCTIONAL_GROUPS["di_isocyanate_monomer"]

    mol = Chem.MolFromSmiles("O=C=O")
    pattern = Chem.MolFromSmarts(group["smarts_1"])

    assert mol.GetSubstructMatches(pattern) == ()


# ============================================================================
# Detector behavior
# ============================================================================

def test_diisocyanate_qualifies():
    functionality_count, count_1, count_2, matches = detect(
        "O=C=NCCCCCCN=C=O"
    )

    assert functionality_count == 2
    assert count_1 == 2
    assert count_2 is None
    assert len(matches) == 2


def test_monoisocyanate_is_rejected_as_diisocyanate():
    functionality_count, count_1, count_2, matches = detect(
        "CN=C=O"
    )

    assert functionality_count == 0
    assert count_1 == 1
    assert count_2 is None
    assert len(matches) == 1


def test_non_isocyanate_is_rejected():
    result = detect("NCCCCCCN")

    assert result == (0, 0, None, ())