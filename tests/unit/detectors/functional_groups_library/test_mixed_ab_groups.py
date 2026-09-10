import pytest
from rdkit import Chem

from AutoREACTER.detectors.functional_groups_detector import (
    FunctionalGroupsDetector,
)
from AutoREACTER.detectors.functional_groups_library.mixed_ab_groups import (
    FUNCTIONAL_GROUPS,
)


EXPECTED_GROUPS = {
    "hydroxy_carboxylic_acid_monomer",
    "hydroxy_acid_halides_monomer",
    "amino_acid_monomer",
    "carboxylic_acid_acid_halide_monomer",
    "hydroxy_thiol_monomer",
}


def detect(group_key, smiles):
    detector = FunctionalGroupsDetector()
    group = FUNCTIONAL_GROUPS[group_key]

    return detector.detect_monomer_functionality(
        Chem.MolFromSmiles(smiles),
        group["functionality_type"],
        group["smarts_1"],
        group["smarts_2"],
    )


# ============================================================================
# Structure
# ============================================================================

def test_expected_active_groups():
    assert set(FUNCTIONAL_GROUPS) == EXPECTED_GROUPS


@pytest.mark.parametrize("group_key", sorted(EXPECTED_GROUPS))
def test_all_groups_are_di_different(group_key):
    assert (
        FUNCTIONAL_GROUPS[group_key]["functionality_type"]
        == "di_different"
    )


@pytest.mark.parametrize("group_key", sorted(EXPECTED_GROUPS))
def test_required_fields_exist(group_key):
    group = FUNCTIONAL_GROUPS[group_key]

    assert "functionality_type" in group
    assert "smarts_1" in group
    assert "smarts_2" in group
    assert "group_name" in group
    assert "comments" in group


@pytest.mark.parametrize("group_key", sorted(EXPECTED_GROUPS))
def test_both_smarts_patterns_compile(group_key):
    group = FUNCTIONAL_GROUPS[group_key]

    assert Chem.MolFromSmarts(group["smarts_1"]) is not None
    assert Chem.MolFromSmarts(group["smarts_2"]) is not None


# ============================================================================
# Hydroxy carboxylic acid
# ============================================================================

def test_hydroxy_carboxylic_acid_metadata():
    group = FUNCTIONAL_GROUPS[
        "hydroxy_carboxylic_acid_monomer"
    ]

    assert group["group_name"] == "hydroxy_carboxylic_acid"


def test_lactic_acid_qualifies_as_hydroxy_carboxylic_acid():
    functionality_count, count_1, count_2, matches = detect(
        "hydroxy_carboxylic_acid_monomer",
        "CC(O)C(=O)O",
    )

    assert functionality_count == 2
    assert count_1 == 1
    assert count_2 == 1
    assert len(matches) == 2


def test_carboxylic_acid_without_alcohol_is_rejected():
    functionality_count, count_1, count_2, matches = detect(
        "hydroxy_carboxylic_acid_monomer",
        "CC(=O)O",
    )

    assert functionality_count == 0
    assert count_1 == 0
    assert count_2 == 1
    assert len(matches) == 1


def test_diol_without_carboxylic_acid_is_rejected():
    functionality_count, count_1, count_2, matches = detect(
        "hydroxy_carboxylic_acid_monomer",
        "OCCO",
    )

    assert functionality_count == 0
    assert count_1 == 2
    assert count_2 == 0
    assert len(matches) == 2


# ============================================================================
# Hydroxy acid halide
# ============================================================================

def test_hydroxy_acid_halide_metadata():
    group = FUNCTIONAL_GROUPS[
        "hydroxy_acid_halides_monomer"
    ]

    assert group["group_name"] == "hydroxy_acid_halide"


def test_hydroxy_acid_chloride_qualifies():
    functionality_count, count_1, count_2, matches = detect(
        "hydroxy_acid_halides_monomer",
        "OCCC(=O)Cl",
    )

    assert functionality_count == 2
    assert count_1 == 1
    assert count_2 == 1
    assert len(matches) == 2


def test_acid_chloride_without_hydroxyl_is_rejected():
    functionality_count, count_1, count_2, matches = detect(
        "hydroxy_acid_halides_monomer",
        "CCC(=O)Cl",
    )

    assert functionality_count == 0
    assert count_1 == 0
    assert count_2 == 1
    assert len(matches) == 1


def test_alcohol_without_acid_halide_is_rejected():
    functionality_count, count_1, count_2, matches = detect(
        "hydroxy_acid_halides_monomer",
        "CCCO",
    )

    assert functionality_count == 0
    assert count_1 == 1
    assert count_2 == 0
    assert len(matches) == 1


# ============================================================================
# Amino acid
# ============================================================================

def test_amino_acid_metadata():
    group = FUNCTIONAL_GROUPS["amino_acid_monomer"]

    assert group["group_name"] == "amino_acid"


def test_glycine_qualifies_as_amino_acid():
    functionality_count, count_1, count_2, matches = detect(
        "amino_acid_monomer",
        "NCC(=O)O",
    )

    assert functionality_count == 2
    assert count_1 == 1
    assert count_2 == 1
    assert len(matches) == 2


def test_amine_without_carboxylic_acid_is_rejected():
    functionality_count, count_1, count_2, matches = detect(
        "amino_acid_monomer",
        "CCN",
    )

    assert functionality_count == 0
    assert count_1 == 1
    assert count_2 == 0
    assert len(matches) == 1


def test_carboxylic_acid_without_amine_is_rejected():
    functionality_count, count_1, count_2, matches = detect(
        "amino_acid_monomer",
        "CC(=O)O",
    )

    assert functionality_count == 0
    assert count_1 == 0
    assert count_2 == 1
    assert len(matches) == 1


# ============================================================================
# Carboxylic acid + acid halide
# ============================================================================

def test_carboxylic_acid_acid_halide_metadata():
    group = FUNCTIONAL_GROUPS[
        "carboxylic_acid_acid_halide_monomer"
    ]

    assert group["group_name"] == "carboxylic_acid_acid_halide"


def test_acid_acid_chloride_qualifies():
    functionality_count, count_1, count_2, matches = detect(
        "carboxylic_acid_acid_halide_monomer",
        "O=C(O)CCC(=O)Cl",
    )

    assert functionality_count == 2
    assert count_1 == 1
    assert count_2 == 1
    assert len(matches) == 2


def test_dicarboxylic_acid_without_halide_is_rejected():
    functionality_count, count_1, count_2, matches = detect(
        "carboxylic_acid_acid_halide_monomer",
        "O=C(O)CCC(=O)O",
    )

    assert functionality_count == 0
    assert count_1 == 2
    assert count_2 == 0
    assert len(matches) == 2


def test_diacid_chloride_without_acid_is_rejected():
    functionality_count, count_1, count_2, matches = detect(
        "carboxylic_acid_acid_halide_monomer",
        "O=C(Cl)CCC(=O)Cl",
    )

    assert functionality_count == 0
    assert count_1 == 0
    assert count_2 == 2
    assert len(matches) == 2


# ============================================================================
# Hydroxy thiol
# ============================================================================

def test_hydroxy_thiol_metadata():
    group = FUNCTIONAL_GROUPS["hydroxy_thiol_monomer"]

    assert group["group_name"] == "hydroxy_thiol"


def test_hydroxy_thiol_qualifies():
    functionality_count, count_1, count_2, matches = detect(
        "hydroxy_thiol_monomer",
        "OCCS",
    )

    assert functionality_count == 2
    assert count_1 == 1
    assert count_2 == 1
    assert len(matches) == 2


def test_alcohol_without_thiol_is_rejected():
    functionality_count, count_1, count_2, matches = detect(
        "hydroxy_thiol_monomer",
        "CCO",
    )

    assert functionality_count == 0
    assert count_1 == 1
    assert count_2 == 0
    assert len(matches) == 1


def test_thiol_without_hydroxyl_is_rejected():
    functionality_count, count_1, count_2, matches = detect(
        "hydroxy_thiol_monomer",
        "CCS",
    )

    assert functionality_count == 0
    assert count_1 == 0
    assert count_2 == 1
    assert len(matches) == 1


# ============================================================================
# Multiplicity preservation
# ============================================================================

def test_di_different_preserves_unequal_site_counts():
    """
    A di_different monomer only requires >=1 of each type.

    Additional sites must remain visible in the returned counts because
    downstream logic may use them for branching / reaction enumeration.
    """
    functionality_count, count_1, count_2, matches = detect(
        "hydroxy_carboxylic_acid_monomer",
        "O=C(O)CC(O)(CC(=O)O)C(=O)O",
    )

    assert functionality_count == 2
    assert count_1 == 1
    assert count_2 == 3
    assert len(matches) == 4