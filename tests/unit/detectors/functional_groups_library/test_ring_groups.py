from rdkit import Chem

from AutoREACTER.detectors.functional_groups_detector import (
    FunctionalGroupsDetector,
)
from AutoREACTER.detectors.functional_groups_library.ring_groups import (
    FUNCTIONAL_GROUPS,
)


GROUP_KEY = "di_epoxy_monomer"


def detect(smiles):
    detector = FunctionalGroupsDetector()
    group = FUNCTIONAL_GROUPS[GROUP_KEY]

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
        "di_epoxy_monomer",
    }


def test_di_epoxy_metadata():
    group = FUNCTIONAL_GROUPS[GROUP_KEY]

    assert group["functionality_type"] == "di_identical"
    assert group["group_name"] == "di_epoxide"
    assert group["comments"] is None


def test_di_epoxy_smarts_compiles():
    group = FUNCTIONAL_GROUPS[GROUP_KEY]

    pattern = Chem.MolFromSmarts(group["smarts_1"])

    assert pattern is not None


# ============================================================================
# Epoxide SMARTS
# ============================================================================

def test_single_epoxide_ring_produces_one_match():
    """
    Ethylene oxide contains one three-membered C-O-C epoxide ring.
    """
    group = FUNCTIONAL_GROUPS[GROUP_KEY]

    mol = Chem.MolFromSmiles("C1CO1")
    pattern = Chem.MolFromSmarts(group["smarts_1"])

    matches = mol.GetSubstructMatches(pattern)

    assert len(matches) == 1


def test_two_separate_epoxide_rings_produce_two_matches():
    """
    Molecule containing two independent oxirane rings should expose
    two epoxide reactive sites.
    """
    group = FUNCTIONAL_GROUPS[GROUP_KEY]

    mol = Chem.MolFromSmiles("C1OC1CC2CO2")
    pattern = Chem.MolFromSmarts(group["smarts_1"])

    matches = mol.GetSubstructMatches(pattern)

    assert len(matches) == 2


def test_acyclic_ether_does_not_match_epoxide():
    group = FUNCTIONAL_GROUPS[GROUP_KEY]

    mol = Chem.MolFromSmiles("COC")
    pattern = Chem.MolFromSmarts(group["smarts_1"])

    assert mol.GetSubstructMatches(pattern) == ()


def test_tetrahydrofuran_does_not_match_epoxide():
    """
    A five-membered cyclic ether must not be mistaken for an epoxide.
    """
    group = FUNCTIONAL_GROUPS[GROUP_KEY]

    mol = Chem.MolFromSmiles("C1CCOC1")
    pattern = Chem.MolFromSmarts(group["smarts_1"])

    assert mol.GetSubstructMatches(pattern) == ()


def test_cyclopropane_does_not_match_epoxide():
    """
    Ring size alone is insufficient; the epoxide oxygen must be present.
    """
    group = FUNCTIONAL_GROUPS[GROUP_KEY]

    mol = Chem.MolFromSmiles("C1CC1")
    pattern = Chem.MolFromSmarts(group["smarts_1"])

    assert mol.GetSubstructMatches(pattern) == ()


# ============================================================================
# Detector behavior
# ============================================================================

def test_single_epoxide_is_rejected_as_di_epoxide():
    functionality_count, count_1, count_2, matches = detect(
        "C1CO1"
    )

    assert functionality_count == 0
    assert count_1 == 1
    assert count_2 is None
    assert len(matches) == 1


def test_two_epoxide_sites_qualify_as_di_epoxide():
    functionality_count, count_1, count_2, matches = detect(
        "C1OC1CC2CO2"
    )

    assert functionality_count == 2
    assert count_1 == 2
    assert count_2 is None
    assert len(matches) == 2


def test_non_epoxide_is_rejected():
    result = detect("COC")

    assert result == (0, 0, None, ())