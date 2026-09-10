import pytest
from rdkit import Chem

from AutoREACTER.detectors.functional_groups_detector import (
    FunctionalGroupsDetector,
)
from AutoREACTER.detectors.functional_groups_library.carboxyl_and_carbonyl_groups import (
    FUNCTIONAL_GROUPS,
)


EXPECTED_GROUPS = {
    "di_carboxylic_acid_monomer",
    "di_carboxylic_acid_halide_monomer",
    "di_carboxylic_ester_monomer",
    "phosgene_monomer",
    "diphenyl_carbonate_monomer",
}


def detect(group_key, smiles):
    detector = FunctionalGroupsDetector()
    group = FUNCTIONAL_GROUPS[group_key]
    mol = Chem.MolFromSmiles(smiles)

    return detector.detect_monomer_functionality(
        mol,
        group["functionality_type"],
        group["smarts_1"],
    )


# ============================================================================
# Structure
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


# ============================================================================
# Dicarboxylic acid
# ============================================================================

def test_dicarboxylic_acid_metadata():
    group = FUNCTIONAL_GROUPS["di_carboxylic_acid_monomer"]

    assert group["functionality_type"] == "di_identical"
    assert group["group_name"] == "di_carboxylic_acid"


def test_adipic_acid_has_two_carboxylic_acid_sites():
    result = detect(
        "di_carboxylic_acid_monomer",
        "O=C(O)CCCCC(=O)O",
    )

    functionality_count, count_1, count_2, matches = result

    assert functionality_count == 2
    assert count_1 == 2
    assert count_2 is None
    assert len(matches) == 2


def test_monocarboxylic_acid_is_rejected_as_dicarboxylic():
    result = detect(
        "di_carboxylic_acid_monomer",
        "CC(=O)O",
    )

    functionality_count, count_1, _, matches = result

    assert functionality_count == 0
    assert count_1 == 1
    assert len(matches) == 1


def test_ester_does_not_match_carboxylic_acid():
    result = detect(
        "di_carboxylic_acid_monomer",
        "CC(=O)OC",
    )

    assert result == (0, 0, None, ())


# ============================================================================
# Dicarboxylic acid halide
# ============================================================================

def test_dicarboxylic_acid_halide_metadata():
    group = FUNCTIONAL_GROUPS[
        "di_carboxylic_acid_halide_monomer"
    ]

    assert group["functionality_type"] == "di_identical"
    assert group["group_name"] == "di_carboxylic_acid_halide"


def test_terephthaloyl_chloride_has_two_acid_halide_sites():
    result = detect(
        "di_carboxylic_acid_halide_monomer",
        "O=C(Cl)c1ccc(C(=O)Cl)cc1",
    )

    functionality_count, count_1, count_2, matches = result

    assert functionality_count == 2
    assert count_1 == 2
    assert count_2 is None
    assert len(matches) == 2


def test_single_acid_chloride_is_rejected_as_di_acid_halide():
    result = detect(
        "di_carboxylic_acid_halide_monomer",
        "CC(=O)Cl",
    )

    functionality_count, count_1, _, matches = result

    assert functionality_count == 0
    assert count_1 == 1
    assert len(matches) == 1


# ============================================================================
# Dicarboxylic ester
# ============================================================================

def test_dicarboxylic_ester_metadata():
    group = FUNCTIONAL_GROUPS["di_carboxylic_ester_monomer"]

    assert group["functionality_type"] == "di_identical"
    assert group["group_name"] == "di_carboxylic_ester"


def test_diester_has_two_ester_sites():
    result = detect(
        "di_carboxylic_ester_monomer",
        "COC(=O)CCCCC(=O)OC",
    )

    functionality_count, count_1, count_2, matches = result

    assert functionality_count == 2
    assert count_1 == 2
    assert count_2 is None
    assert len(matches) == 2


def test_single_ester_is_rejected_as_diester():
    result = detect(
        "di_carboxylic_ester_monomer",
        "CC(=O)OC",
    )

    functionality_count, count_1, _, matches = result

    assert functionality_count == 0
    assert count_1 == 1
    assert len(matches) == 1


def test_carboxylic_acid_does_not_match_ester_pattern():
    result = detect(
        "di_carboxylic_ester_monomer",
        "CC(=O)O",
    )

    assert result == (0, 0, None, ())


# ============================================================================
# Phosgene
# ============================================================================

def test_phosgene_metadata():
    group = FUNCTIONAL_GROUPS["phosgene_monomer"]

    assert group["functionality_type"] == "mono"
    assert group["group_name"] == "phosgene"


def test_phosgene_pattern_matches_phosgene_once():
    result = detect(
        "phosgene_monomer",
        "O=C(Cl)Cl",
    )

    functionality_count, count_1, count_2, matches = result

    assert functionality_count == 1
    assert count_1 == 1
    assert count_2 is None
    assert len(matches) == 1


@pytest.mark.parametrize(
    "smiles",
    [
        "CC(=O)Cl",
        "O=C(O)O",
        "O=C(OC)OC",
    ],
)
def test_phosgene_pattern_rejects_non_phosgene_carbonyls(smiles):
    result = detect(
        "phosgene_monomer",
        smiles,
    )

    assert result == (0, 0, None, ())


# ============================================================================
# Diphenyl carbonate
# ============================================================================

def test_diphenyl_carbonate_metadata():
    group = FUNCTIONAL_GROUPS["diphenyl_carbonate_monomer"]

    assert group["functionality_type"] == "di_identical"
    assert group["group_name"] == "diphenyl_carbonate"


def test_diphenyl_carbonate_has_two_aryl_carbonate_matches():
    result = detect(
        "diphenyl_carbonate_monomer",
        "O=C(Oc1ccccc1)Oc1ccccc1",
    )

    functionality_count, count_1, count_2, matches = result

    assert functionality_count == 2
    assert count_1 == 2
    assert count_2 is None
    assert len(matches) == 2


def test_single_aryl_carbonate_side_is_rejected():
    """
    Methyl phenyl carbonate contains only one O-aryl carbonate match,
    so it cannot satisfy the di_identical >= 2 rule.
    """
    result = detect(
        "diphenyl_carbonate_monomer",
        "O=C(Oc1ccccc1)OC",
    )

    functionality_count, count_1, _, matches = result

    assert functionality_count == 0
    assert count_1 == 1
    assert len(matches) == 1


def test_dimethyl_carbonate_does_not_match_diphenyl_carbonate():
    result = detect(
        "diphenyl_carbonate_monomer",
        "O=C(OC)OC",
    )

    assert result == (0, 0, None, ())