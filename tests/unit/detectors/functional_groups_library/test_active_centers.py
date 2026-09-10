from rdkit import Chem

from AutoREACTER.detectors.functional_groups_detector import (
    FunctionalGroupsDetector,
)
from AutoREACTER.detectors.functional_groups_library.active_centers import (
    FUNCTIONAL_GROUPS,
)


def test_expected_active_groups():
    assert set(FUNCTIONAL_GROUPS) == {
        "vinyl_chain_end_radical",
    }


def test_vinyl_chain_end_radical_metadata():
    group = FUNCTIONAL_GROUPS["vinyl_chain_end_radical"]

    assert group["functionality_type"] == "vinyl"
    assert group["group_name"] == "vinyl_chain_end_radical"
    assert group["comments"] is None


def test_vinyl_chain_end_radical_smarts_compiles():
    group = FUNCTIONAL_GROUPS["vinyl_chain_end_radical"]

    pattern = Chem.MolFromSmarts(group["smarts_1"])

    assert pattern is not None


def test_vinyl_chain_end_radical_smarts_matches_trivalent_non_ring_carbon():
    """
    A neutral trivalent carbon radical with three carbon neighbors
    should match the active-center SMARTS.
    """
    group = FUNCTIONAL_GROUPS["vinyl_chain_end_radical"]

    mol = Chem.MolFromSmiles("C[C](C)C")
    pattern = Chem.MolFromSmarts(group["smarts_1"])

    matches = mol.GetSubstructMatches(pattern)

    assert len(matches) == 1


def test_vinyl_chain_end_radical_smarts_does_not_match_saturated_carbon():
    """
    A normal tetravalent saturated carbon should not be identified
    as a chain-end radical.
    """
    group = FUNCTIONAL_GROUPS["vinyl_chain_end_radical"]

    mol = Chem.MolFromSmiles("CC(C)C")
    pattern = Chem.MolFromSmarts(group["smarts_1"])

    assert mol.GetSubstructMatches(pattern) == ()


def test_vinyl_chain_end_radical_smarts_does_not_match_alkene_carbon():
    group = FUNCTIONAL_GROUPS["vinyl_chain_end_radical"]

    mol = Chem.MolFromSmiles("C=C(C)C")
    pattern = Chem.MolFromSmarts(group["smarts_1"])

    assert mol.GetSubstructMatches(pattern) == ()


def test_vinyl_chain_end_radical_smarts_excludes_ring_carbon():
    """
    The !R condition should prevent ring atoms from being classified
    as vinyl chain-end radicals.
    """
    group = FUNCTIONAL_GROUPS["vinyl_chain_end_radical"]

    mol = Chem.MolFromSmiles("[C]1(C)CCC1")
    pattern = Chem.MolFromSmarts(group["smarts_1"])

    assert mol.GetSubstructMatches(pattern) == ()


def test_vinyl_chain_end_radical_qualifies_through_detector():
    detector = FunctionalGroupsDetector()
    group = FUNCTIONAL_GROUPS["vinyl_chain_end_radical"]

    mol = Chem.MolFromSmiles("C[C](C)C")

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


def test_saturated_carbon_does_not_qualify_through_detector():
    detector = FunctionalGroupsDetector()
    group = FUNCTIONAL_GROUPS["vinyl_chain_end_radical"]

    mol = Chem.MolFromSmiles("CC(C)C")

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