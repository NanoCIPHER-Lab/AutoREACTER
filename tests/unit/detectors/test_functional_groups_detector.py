from types import SimpleNamespace

import pytest
from rdkit import Chem

from AutoREACTER.detectors.functional_groups_detector import (
    FunctionalGroupInfo,
    FunctionalGroupsDetector,
    FunctionalGroupVisualization,
    MonomerRole,
)


# ============================================================================
# Helpers
# ============================================================================

def make_monomer(name: str, smiles: str):
    """Create a minimal monomer-like object for detector unit tests."""
    return SimpleNamespace(
        name=name,
        smiles=smiles,
        rdkit_mol=Chem.MolFromSmiles(smiles),
    )


def make_session(monomers):
    """Create the minimal Session interface required by the detector."""
    return SimpleNamespace(
        inputs=SimpleNamespace(monomers=monomers),
        monomer_roles=None,
    )


def make_index_role(
    name,
    smiles,
    indexes,
    *,
    is_looped=False,
):
    """Create a minimal role-like object for index-based detection."""
    return SimpleNamespace(
        name=name,
        smiles=smiles,
        rdkit_mol=Chem.MolFromSmiles(smiles),
        indexes_in_template=indexes,
        is_looped=is_looped,
    )


# ============================================================================
# detect_monomer_functionality()
# ============================================================================

def test_mono_single_match_qualifies():
    detector = FunctionalGroupsDetector()
    mol = Chem.MolFromSmiles("CCO")

    result = detector.detect_monomer_functionality(
        mol,
        "mono",
        "[OX2H1]",
    )

    functionality_count, count_1, count_2, matches = result

    assert functionality_count == 1
    assert count_1 == 1
    assert count_2 is None
    assert len(matches) == 1


def test_mono_no_match_is_rejected():
    detector = FunctionalGroupsDetector()
    mol = Chem.MolFromSmiles("CCC")

    result = detector.detect_monomer_functionality(
        mol,
        "mono",
        "[OX2H1]",
    )

    assert result == (0, 0, None, ())


def test_vinyl_single_match_qualifies():
    detector = FunctionalGroupsDetector()
    mol = Chem.MolFromSmiles("C=C")

    functionality_count, count_1, count_2, matches = (
        detector.detect_monomer_functionality(
            mol,
            "vinyl",
            "[CH2]=[CH2]",
        )
    )

    assert functionality_count == 1
    assert count_1 == 1
    assert count_2 is None
    assert len(matches) == 1


def test_vinyl_multiple_sites_still_returns_presence_based_functionality():
    detector = FunctionalGroupsDetector()
    mol = Chem.MolFromSmiles("C=CC=C")

    functionality_count, count_1, count_2, matches = (
        detector.detect_monomer_functionality(
            mol,
            "vinyl",
            "[CH2]=[C]",
        )
    )

    assert functionality_count == 1
    assert count_1 >= 1
    assert count_2 is None
    assert len(matches) == count_1


def test_di_identical_two_sites_qualifies():
    detector = FunctionalGroupsDetector()
    mol = Chem.MolFromSmiles("OCCO")

    functionality_count, count_1, count_2, matches = (
        detector.detect_monomer_functionality(
            mol,
            "di_identical",
            "[OX2H1]",
        )
    )

    assert functionality_count == 2
    assert count_1 == 2
    assert count_2 is None
    assert len(matches) == 2


def test_di_identical_one_site_is_rejected():
    detector = FunctionalGroupsDetector()
    mol = Chem.MolFromSmiles("CCO")

    functionality_count, count_1, count_2, matches = (
        detector.detect_monomer_functionality(
            mol,
            "di_identical",
            "[OX2H1]",
        )
    )

    assert functionality_count == 0
    assert count_1 == 1
    assert count_2 is None
    assert len(matches) == 1


def test_di_identical_more_than_two_sites_qualifies_and_preserves_count():
    detector = FunctionalGroupsDetector()
    mol = Chem.MolFromSmiles("OCC(O)CO")

    functionality_count, count_1, count_2, matches = (
        detector.detect_monomer_functionality(
            mol,
            "di_identical",
            "[OX2H1]",
        )
    )

    assert functionality_count == 2
    assert count_1 == 3
    assert count_2 is None
    assert len(matches) == 3


def test_di_different_requires_both_patterns():
    detector = FunctionalGroupsDetector()
    mol = Chem.MolFromSmiles("NCC(=O)O")

    functionality_count, count_1, count_2, matches = (
        detector.detect_monomer_functionality(
            mol,
            "di_different",
            "[NX3H2]",
            "[CX3](=[O])[OX2H1]",
        )
    )

    assert functionality_count == 2
    assert count_1 == 1
    assert count_2 == 1
    assert len(matches) == 2


def test_di_different_missing_second_pattern_is_rejected():
    detector = FunctionalGroupsDetector()
    mol = Chem.MolFromSmiles("CCN")

    functionality_count, count_1, count_2, matches = (
        detector.detect_monomer_functionality(
            mol,
            "di_different",
            "[NX3H2]",
            "[CX3](=[O])[OX2H1]",
        )
    )

    assert functionality_count == 0
    assert count_1 == 1
    assert count_2 == 0
    assert len(matches) == 1


def test_invalid_primary_smarts_is_rejected():
    detector = FunctionalGroupsDetector()
    mol = Chem.MolFromSmiles("CCO")

    result = detector.detect_monomer_functionality(
        mol,
        "mono",
        "[THIS_IS_INVALID",
    )

    assert result == (0, None, None, None)


def test_invalid_secondary_smarts_is_rejected():
    detector = FunctionalGroupsDetector()
    mol = Chem.MolFromSmiles("NCC(=O)O")

    result = detector.detect_monomer_functionality(
        mol,
        "di_different",
        "[NX3H2]",
        "[THIS_IS_INVALID",
    )

    assert result == (0, None, None, None)


def test_none_molecule_is_rejected():
    detector = FunctionalGroupsDetector()

    result = detector.detect_monomer_functionality(
        None,
        "mono",
        "[OX2H1]",
    )

    assert result == (0, None, None, None)


# ============================================================================
# functional_groups_detector()
# ============================================================================

def test_detector_stores_detected_roles_in_session():
    detector = FunctionalGroupsDetector()

    detector.monomer_types = {
        "test_diol": {
            "functionality_type": "di_identical",
            "smarts_1": "[OX2H1]",
            "group_name": "diol",
            "comments": None,
        }
    }

    monomer = make_monomer("ethylene_glycol", "OCCO")
    session = make_session([monomer])

    result = detector.functional_groups_detector(session)

    assert result is None
    assert len(session.monomer_roles) == 1

    role = session.monomer_roles[0]

    assert isinstance(role, MonomerRole)
    assert role.name == "ethylene_glycol"
    assert role.smiles == "OCCO"
    assert role.is_monomer is True
    assert len(role.functionalities) == 1


def test_detector_records_functional_group_information():
    detector = FunctionalGroupsDetector()

    detector.monomer_types = {
        "test_diol": {
            "functionality_type": "di_identical",
            "smarts_1": "[OX2H1]",
            "group_name": "diol",
            "comments": None,
        }
    }

    session = make_session(
        [make_monomer("ethylene_glycol", "OCCO")]
    )

    detector.functional_groups_detector(session)

    fg = session.monomer_roles[0].functionalities[0]

    assert isinstance(fg, FunctionalGroupInfo)
    assert fg.functionality_type == "di_identical"
    assert fg.fg_name == "diol"
    assert fg.fg_smarts_1 == "[OX2H1]"
    assert fg.fg_count_1 == 2
    assert fg.fg_smarts_2 is None
    assert fg.fg_count_2 is None


def test_detector_can_record_multiple_functionalities_for_same_monomer():
    detector = FunctionalGroupsDetector()

    detector.monomer_types = {
        "amine": {
            "functionality_type": "mono",
            "smarts_1": "[NX3H2]",
            "group_name": "primary_amine",
            "comments": None,
        },
        "acid": {
            "functionality_type": "mono",
            "smarts_1": "[CX3](=[O])[OX2H1]",
            "group_name": "carboxylic_acid",
            "comments": None,
        },
    }

    session = make_session(
        [make_monomer("glycine", "NCC(=O)O")]
    )

    detector.functional_groups_detector(session)

    functionalities = session.monomer_roles[0].functionalities

    assert len(functionalities) == 2

    names = {fg.fg_name for fg in functionalities}

    assert names == {
        "primary_amine",
        "carboxylic_acid",
    }


def test_detector_rejects_di_identical_monomer_with_only_one_site():
    detector = FunctionalGroupsDetector()

    detector.monomer_types = {
        "test_diol": {
            "functionality_type": "di_identical",
            "smarts_1": "[OX2H1]",
            "group_name": "diol",
            "comments": None,
        }
    }

    session = make_session(
        [make_monomer("ethanol", "CCO")]
    )

    with pytest.raises(
        RuntimeError,
        match="No functional groups were detected",
    ):
        detector.functional_groups_detector(session)


def test_detector_raises_when_no_monomer_has_supported_functionality():
    detector = FunctionalGroupsDetector()

    detector.monomer_types = {
        "test_amine": {
            "functionality_type": "mono",
            "smarts_1": "[NX3H2]",
            "group_name": "primary_amine",
            "comments": None,
        }
    }

    session = make_session(
        [
            make_monomer("propane", "CCC"),
            make_monomer("butane", "CCCC"),
        ]
    )

    with pytest.raises(
        RuntimeError,
        match="No functional groups were detected",
    ):
        detector.functional_groups_detector(session)


# ============================================================================
# _functional_groups_detector_for_visualization()
# ============================================================================

def test_visualization_detector_returns_detected_monomer():
    detector = FunctionalGroupsDetector()

    detector.monomer_types = {
        "alcohol": {
            "functionality_type": "mono",
            "smarts_1": "[OX2H1]",
            "group_name": "alcohol",
            "comments": None,
        }
    }

    session = make_session(
        [make_monomer("ethanol", "CCO")]
    )

    result = detector._functional_groups_detector_for_visualization(
        session
    )

    assert len(result) == 1
    assert isinstance(result[0], FunctionalGroupVisualization)
    assert result[0].name == "ethanol"
    assert len(result[0].indexes_to_highlight) == 1


def test_visualization_detector_ignores_nonmatching_monomer():
    detector = FunctionalGroupsDetector()

    detector.monomer_types = {
        "amine": {
            "functionality_type": "mono",
            "smarts_1": "[NX3H2]",
            "group_name": "amine",
            "comments": None,
        }
    }

    session = make_session(
        [make_monomer("propane", "CCC")]
    )

    result = detector._functional_groups_detector_for_visualization(
        session
    )

    assert result == []


# ============================================================================
# _detect_functional_groups_by_index()
# ============================================================================

def test_detect_by_index_returns_true_for_overlapping_atom():
    detector = FunctionalGroupsDetector()
    mol = Chem.MolFromSmiles("CCO")

    oxygen_index = 2

    assert detector._detect_functional_groups_by_index(
        mol,
        "[OX2H1]",
        [oxygen_index],
    ) is True


def test_detect_by_index_returns_false_when_indices_do_not_overlap():
    detector = FunctionalGroupsDetector()
    mol = Chem.MolFromSmiles("CCO")

    assert detector._detect_functional_groups_by_index(
        mol,
        "[OX2H1]",
        [0],
    ) is False


def test_detect_by_index_returns_false_for_empty_indices():
    detector = FunctionalGroupsDetector()
    mol = Chem.MolFromSmiles("CCO")

    assert detector._detect_functional_groups_by_index(
        mol,
        "[OX2H1]",
        [],
    ) is False


def test_detect_by_index_invalid_smarts_returns_false():
    detector = FunctionalGroupsDetector()
    mol = Chem.MolFromSmiles("CCO")

    assert detector._detect_functional_groups_by_index(
        mol,
        "[INVALID",
        [2],
    ) is False


# ============================================================================
# index_based_functional_groups_detector()
# ============================================================================

def test_index_detector_allows_single_di_identical_site():
    """
    Index-based detection intentionally ignores the normal >=2 rule for
    di_identical groups.

    One matching functional group touching the target index is sufficient.
    """
    detector = FunctionalGroupsDetector()

    detector.monomer_types = {
        "diol": {
            "functionality_type": "di_identical",
            "smarts_1": "[OX2H1]",
            "group_name": "diol",
            "comments": None,
        }
    }

    role = make_index_role(
        "ethanol_product",
        "CCO",
        [2],
    )

    result = detector.index_based_functional_groups_detector(
        [role]
    )

    assert result is not False
    assert len(result) == 1

    output_role = result[0]

    assert output_role.is_monomer is False
    assert len(output_role.functionalities) == 1

    fg = output_role.functionalities[0]

    assert fg.functionality_type == "di_identical"
    assert fg.fg_count_1 == 1
    assert fg.fg_1_indexes is not None


def test_index_detector_requires_overlap():
    detector = FunctionalGroupsDetector()

    detector.monomer_types = {
        "alcohol": {
            "functionality_type": "mono",
            "smarts_1": "[OX2H1]",
            "group_name": "alcohol",
            "comments": None,
        }
    }

    role = make_index_role(
        "ethanol_product",
        "CCO",
        [0],
    )

    result = detector.index_based_functional_groups_detector(
        [role]
    )

    assert result is False


def test_index_detector_di_different_requires_hit_for_both_patterns():
    detector = FunctionalGroupsDetector()

    detector.monomer_types = {
        "amino_acid": {
            "functionality_type": "di_different",
            "smarts_1": "[N]",
            "smarts_2": "[OX2H1]",
            "group_name": "amino_acid",
            "comments": None,
        }
    }

    mol = Chem.MolFromSmiles("NCC(=O)O")

    nitrogen_index = next(
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetSymbol() == "N"
    )

    hydroxyl_oxygen_index = next(
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetSymbol() == "O"
        and atom.GetTotalNumHs() == 1
    )

    role = SimpleNamespace(
        name="glycine_product",
        smiles="NCC(=O)O",
        rdkit_mol=mol,
        indexes_in_template=[
            nitrogen_index,
            hydroxyl_oxygen_index,
        ],
        is_looped=False,
    )

    result = detector.index_based_functional_groups_detector(
        [role]
    )

    assert result is not False
    assert len(result) == 1

    fg = result[0].functionalities[0]

    assert fg.fg_count_1 == 1
    assert fg.fg_count_2 == 1
    assert fg.fg_1_indexes is not None
    assert fg.fg_2_indexes is not None


def test_index_detector_di_different_fails_when_only_one_pattern_overlaps():
    detector = FunctionalGroupsDetector()

    detector.monomer_types = {
        "amino_acid": {
            "functionality_type": "di_different",
            "smarts_1": "[N]",
            "smarts_2": "[OX2H1]",
            "group_name": "amino_acid",
            "comments": None,
        }
    }

    mol = Chem.MolFromSmiles("NCC(=O)O")

    nitrogen_index = next(
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetSymbol() == "N"
    )

    role = SimpleNamespace(
        name="glycine_product",
        smiles="NCC(=O)O",
        rdkit_mol=mol,
        indexes_in_template=[nitrogen_index],
        is_looped=False,
    )

    result = detector.index_based_functional_groups_detector(
        [role]
    )

    assert result is False


def test_index_detector_skips_looped_monomers():
    detector = FunctionalGroupsDetector()

    detector.monomer_types = {
        "alcohol": {
            "functionality_type": "mono",
            "smarts_1": "[OX2H1]",
            "group_name": "alcohol",
            "comments": None,
        }
    }

    role = make_index_role(
        "ethanol_product",
        "CCO",
        [2],
        is_looped=True,
    )

    result = detector.index_based_functional_groups_detector(
        [role]
    )

    assert result is False


def test_index_detector_preserves_template_indices():
    detector = FunctionalGroupsDetector()

    detector.monomer_types = {
        "alcohol": {
            "functionality_type": "mono",
            "smarts_1": "[OX2H1]",
            "group_name": "alcohol",
            "comments": None,
        }
    }

    role = make_index_role(
        "ethanol_product",
        "CCO",
        [2],
    )

    result = detector.index_based_functional_groups_detector(
        [role]
    )

    assert result[0].indexes_in_template == [2]


# ============================================================================
# functional_group_highlighted_molecules_image_grid()
# ============================================================================

def test_image_grid_flattens_highlighted_atom_indices(monkeypatch):
    detector = FunctionalGroupsDetector()

    mol = Chem.MolFromSmiles("CCO")

    visualization = FunctionalGroupVisualization(
        monomer=mol,
        name="ethanol",
        indexes_to_highlight=((1, 2), (2,)),
    )

    monkeypatch.setattr(
        detector,
        "_functional_groups_detector_for_visualization",
        lambda session: [visualization],
    )

    captured = {}

    def fake_grid(
        molecules,
        molsPerRow,
        legends,
        subImgSize,
        highlightAtomLists,
    ):
        captured["molecules"] = molecules
        captured["molsPerRow"] = molsPerRow
        captured["legends"] = legends
        captured["subImgSize"] = subImgSize
        captured["highlightAtomLists"] = highlightAtomLists

        return "fake-image"

    monkeypatch.setattr(
        "AutoREACTER.detectors.functional_groups_detector.Draw.MolsToGridImage",
        fake_grid,
    )

    result = detector.functional_group_highlighted_molecules_image_grid(
        SimpleNamespace()
    )

    assert result == "fake-image"
    assert captured["molecules"] == [mol]
    assert captured["legends"] == ["ethanol"]
    assert captured["highlightAtomLists"] == [[1, 2]]
    assert captured["molsPerRow"] == 3
    assert captured["subImgSize"] == (500, 500)
