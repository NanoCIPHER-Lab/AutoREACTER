from dataclasses import dataclass
from types import SimpleNamespace

import pytest
from rdkit import Chem

from AutoREACTER.detectors.non_monomer_detector import (
    NonReactantsDetector,
)


# =============================================================================
# Test helpers
# =============================================================================

@dataclass
class DummyMonomer:
    """
    Minimal dataclass matching what NonReactantsDetector needs.

    A dataclass is intentionally used because non_reactant_selection()
    calls dataclasses.replace().
    """
    id: int
    name: str
    smiles: str
    rdkit_mol: object
    status: bool = True


def make_monomer(
    monomer_id: int,
    name: str,
    smiles: str,
    *,
    status=True,
):
    return DummyMonomer(
        id=monomer_id,
        name=name,
        smiles=smiles,
        rdkit_mol=Chem.MolFromSmiles(smiles),
        status=status,
    )


def make_reaction(smiles_1, smiles_2=None):
    monomer_1 = SimpleNamespace(smiles=smiles_1)

    monomer_2 = (
        SimpleNamespace(smiles=smiles_2)
        if smiles_2 is not None
        else None
    )

    return SimpleNamespace(
        monomer_1=monomer_1,
        monomer_2=monomer_2,
    )


def make_session(
    monomers=None,
    reaction_instances=None,
    non_reactants=None,
):
    return SimpleNamespace(
        inputs=SimpleNamespace(monomers=monomers or []),
        reaction_instances=reaction_instances or [],
        non_reactants=non_reactants or [],
    )


# =============================================================================
# _same_molecule()
# =============================================================================

def test_same_molecule_adds_new_molecule():
    detector = NonReactantsDetector()

    reactants = []

    result = detector._same_molecule(
        reactants,
        "CCO",
    )

    assert result == ["CCO"]


def test_same_molecule_does_not_duplicate_identical_smiles():
    detector = NonReactantsDetector()

    reactants = ["CCO"]

    result = detector._same_molecule(
        reactants,
        "CCO",
    )

    assert result == ["CCO"]


def test_same_molecule_uses_canonical_smiles():
    """
    OCC and CCO represent the same ethanol molecule.

    They should therefore not be stored twice.
    """
    detector = NonReactantsDetector()

    reactants = ["CCO"]

    result = detector._same_molecule(
        reactants,
        "OCC",
    )

    assert result == ["CCO"]


def test_same_molecule_adds_structurally_different_molecule():
    detector = NonReactantsDetector()

    reactants = ["CCO"]

    result = detector._same_molecule(
        reactants,
        "CCN",
    )

    assert result == [
        "CCO",
        "CCN",
    ]


def test_same_molecule_invalid_smiles_is_not_added():
    detector = NonReactantsDetector()

    reactants = ["CCO"]

    result = detector._same_molecule(
        reactants,
        "this_is_not_smiles",
    )

    assert result == ["CCO"]


def test_same_molecule_modifies_original_list():
    """
    Current implementation appends directly to reactants_list.

    The returned object should therefore be the same list object.
    """
    detector = NonReactantsDetector()

    reactants = []

    result = detector._same_molecule(
        reactants,
        "CCO",
    )

    assert result is reactants
    assert reactants == ["CCO"]


# =============================================================================
# _same_molecule_for_initaials()
# =============================================================================

def test_same_molecule_for_initials_finds_identical_molecule():
    detector = NonReactantsDetector()

    reactants = ["CCO"]
    mol = Chem.MolFromSmiles("CCO")

    result = detector._same_molecule_for_initaials(
        reactants,
        mol,
    )

    assert result is True


def test_same_molecule_for_initials_uses_canonical_structure():
    detector = NonReactantsDetector()

    reactants = ["CCO"]
    mol = Chem.MolFromSmiles("OCC")

    result = detector._same_molecule_for_initaials(
        reactants,
        mol,
    )

    assert result is True


def test_same_molecule_for_initials_ignores_explicit_hydrogens():
    """
    Explicit hydrogens should not prevent structural equality.
    """
    detector = NonReactantsDetector()

    reactants = ["CCO"]

    mol = Chem.AddHs(
        Chem.MolFromSmiles("CCO")
    )

    result = detector._same_molecule_for_initaials(
        reactants,
        mol,
    )

    assert result is True


def test_same_molecule_for_initials_returns_false_for_different_molecule():
    detector = NonReactantsDetector()

    reactants = ["CCO"]
    mol = Chem.MolFromSmiles("CCN")

    result = detector._same_molecule_for_initaials(
        reactants,
        mol,
    )

    assert result is False


def test_same_molecule_for_initials_returns_false_for_none():
    detector = NonReactantsDetector()

    result = detector._same_molecule_for_initaials(
        ["CCO"],
        None,
    )

    assert result is False


def test_same_molecule_for_initials_ignores_invalid_smiles_in_reactant_list():
    detector = NonReactantsDetector()

    reactants = [
        "invalid_smiles",
        "CCO",
    ]

    mol = Chem.MolFromSmiles("CCO")

    result = detector._same_molecule_for_initaials(
        reactants,
        mol,
    )

    assert result is True


def test_same_molecule_for_initials_empty_reactant_list_returns_false():
    detector = NonReactantsDetector()

    mol = Chem.MolFromSmiles("CCO")

    result = detector._same_molecule_for_initaials(
        [],
        mol,
    )

    assert result is False


# =============================================================================
# non_monomer_detector()
# =============================================================================

def test_non_monomer_detector_identifies_unused_monomer():
    detector = NonReactantsDetector()

    ethanol = make_monomer(
        1,
        "ethanol",
        "CCO",
    )

    amine = make_monomer(
        2,
        "ethylamine",
        "CCN",
    )

    propane = make_monomer(
        3,
        "propane",
        "CCC",
    )

    reaction = make_reaction(
        "CCO",
        "CCN",
    )

    session = make_session(
        monomers=[
            ethanol,
            amine,
            propane,
        ],
        reaction_instances=[reaction],
    )

    result = detector.non_monomer_detector(session)

    assert result is None
    assert session.non_reactants == [
        propane,
    ]


def test_non_monomer_detector_all_monomers_react():
    detector = NonReactantsDetector()

    ethanol = make_monomer(
        1,
        "ethanol",
        "CCO",
    )

    amine = make_monomer(
        2,
        "ethylamine",
        "CCN",
    )

    session = make_session(
        monomers=[
            ethanol,
            amine,
        ],
        reaction_instances=[
            make_reaction(
                "CCO",
                "CCN",
            )
        ],
    )

    detector.non_monomer_detector(session)

    assert session.non_reactants == []


def test_non_monomer_detector_no_reactions_marks_every_monomer_nonreactant():
    detector = NonReactantsDetector()

    monomers = [
        make_monomer(
            1,
            "ethanol",
            "CCO",
        ),
        make_monomer(
            2,
            "ethylamine",
            "CCN",
        ),
    ]

    session = make_session(
        monomers=monomers,
        reaction_instances=[],
    )

    detector.non_monomer_detector(session)

    assert session.non_reactants == monomers


def test_non_monomer_detector_handles_homopolymerization():
    """
    A ReactionInstance may contain only monomer_1.

    That molecule should still be recognized as a reactant.
    """
    detector = NonReactantsDetector()

    styrene = make_monomer(
        1,
        "styrene",
        "C=Cc1ccccc1",
    )

    solvent = make_monomer(
        2,
        "ethanol",
        "CCO",
    )

    session = make_session(
        monomers=[
            styrene,
            solvent,
        ],
        reaction_instances=[
            make_reaction(
                "C=Cc1ccccc1",
            )
        ],
    )

    detector.non_monomer_detector(session)

    assert session.non_reactants == [
        solvent,
    ]


def test_non_monomer_detector_uses_structural_equivalence_not_raw_string():
    detector = NonReactantsDetector()

    ethanol = make_monomer(
        1,
        "ethanol",
        "OCC",
    )

    session = make_session(
        monomers=[ethanol],
        reaction_instances=[
            make_reaction("CCO")
        ],
    )

    detector.non_monomer_detector(session)

    assert session.non_reactants == []


def test_non_monomer_detector_deduplicates_reactants_across_reactions():
    detector = NonReactantsDetector()

    ethanol = make_monomer(
        1,
        "ethanol",
        "CCO",
    )

    propane = make_monomer(
        2,
        "propane",
        "CCC",
    )

    session = make_session(
        monomers=[
            ethanol,
            propane,
        ],
        reaction_instances=[
            make_reaction("CCO"),
            make_reaction("OCC"),
            make_reaction("CCO"),
        ],
    )

    detector.non_monomer_detector(session)

    assert session.non_reactants == [
        propane,
    ]


# =============================================================================
# non_reactants_to_visualization()
# =============================================================================

def test_non_reactants_visualization_returns_none_for_empty_list():
    detector = NonReactantsDetector()

    session = make_session(
        non_reactants=[],
    )

    result = detector.non_reactants_to_visualization(
        session
    )

    assert result is None


def test_non_reactants_visualization_passes_correct_data_to_rdkit(
    monkeypatch,
):
    detector = NonReactantsDetector()

    monomer_1 = make_monomer(
        1,
        "ethanol",
        "CCO",
    )

    monomer_2 = make_monomer(
        2,
        "ethylamine",
        "CCN",
    )

    session = make_session(
        non_reactants=[
            monomer_1,
            monomer_2,
        ]
    )

    captured = {}

    fake_image = object()

    def fake_grid(
        molecules,
        legends,
        molsPerRow,
        subImgSize,
    ):
        captured["molecules"] = molecules
        captured["legends"] = legends
        captured["molsPerRow"] = molsPerRow
        captured["subImgSize"] = subImgSize

        return fake_image

    monkeypatch.setattr(
        "AutoREACTER.detectors.non_monomer_detector."
        "Draw.MolsToGridImage",
        fake_grid,
    )

    result = detector.non_reactants_to_visualization(
        session
    )

    assert result is fake_image

    assert captured["molecules"] == [
        monomer_1.rdkit_mol,
        monomer_2.rdkit_mol,
    ]

    assert captured["legends"] == [
        "ethanol",
        "ethylamine",
    ]

    assert captured["molsPerRow"] == 3
    assert captured["subImgSize"] == (
        400,
        400,
    )


# =============================================================================
# non_reactant_selection() - no non-reactants
# =============================================================================

def test_selection_does_nothing_when_no_nonreactants():
    detector = NonReactantsDetector()

    monomer = make_monomer(
        1,
        "ethanol",
        "CCO",
    )

    session = make_session(
        monomers=[monomer],
        non_reactants=[],
    )

    result = detector.non_reactant_selection(
        session
    )

    assert result is None
    assert session.inputs.monomers == [
        monomer,
    ]


# =============================================================================
# non_reactant_selection() - one non-reactant
# =============================================================================

def test_single_nonreactant_N_discards_monomer(
    monkeypatch,
):
    detector = NonReactantsDetector()

    monomer = make_monomer(
        1,
        "solvent",
        "CCO",
    )

    session = make_session(
        monomers=[monomer],
        non_reactants=[monomer],
    )

    monkeypatch.setattr(
        "builtins.input",
        lambda prompt: "N",
    )

    detector.non_reactant_selection(session)

    assert session.inputs.monomers[0].status is False


def test_single_nonreactant_A_keeps_monomer(
    monkeypatch,
):
    detector = NonReactantsDetector()

    monomer = make_monomer(
        1,
        "solvent",
        "CCO",
    )

    session = make_session(
        monomers=[monomer],
        non_reactants=[monomer],
    )

    monkeypatch.setattr(
        "builtins.input",
        lambda prompt: "A",
    )

    detector.non_reactant_selection(session)

    assert session.inputs.monomers[0].status is True


def test_single_nonreactant_reprompts_invalid_input(
    monkeypatch,
):
    detector = NonReactantsDetector()

    monomer = make_monomer(
        1,
        "solvent",
        "CCO",
    )

    session = make_session(
        monomers=[monomer],
        non_reactants=[monomer],
    )

    responses = iter(
        [
            "wrong",
            "N",
        ]
    )

    monkeypatch.setattr(
        "builtins.input",
        lambda prompt: next(responses),
    )

    detector.non_reactant_selection(session)

    assert session.inputs.monomers[0].status is False


# =============================================================================
# non_reactant_selection() - multiple non-reactants
# =============================================================================

def test_multiple_nonreactants_N_discards_all(
    monkeypatch,
):
    detector = NonReactantsDetector()

    m1 = make_monomer(
        1,
        "M1",
        "CCO",
    )

    m2 = make_monomer(
        2,
        "M2",
        "CCN",
    )

    reactive = make_monomer(
        3,
        "reactive",
        "CCC",
    )

    session = make_session(
        monomers=[
            m1,
            m2,
            reactive,
        ],
        non_reactants=[
            m1,
            m2,
        ],
    )

    monkeypatch.setattr(
        "builtins.input",
        lambda prompt: "N",
    )

    detector.non_reactant_selection(session)

    status_by_id = {
        m.id: m.status
        for m in session.inputs.monomers
    }

    assert status_by_id == {
        1: False,
        2: False,
        3: True,
    }


def test_multiple_nonreactants_A_keeps_all(
    monkeypatch,
):
    detector = NonReactantsDetector()

    m1 = make_monomer(
        1,
        "M1",
        "CCO",
    )

    m2 = make_monomer(
        2,
        "M2",
        "CCN",
    )

    session = make_session(
        monomers=[
            m1,
            m2,
        ],
        non_reactants=[
            m1,
            m2,
        ],
    )

    monkeypatch.setattr(
        "builtins.input",
        lambda prompt: "A",
    )

    detector.non_reactant_selection(session)

    assert all(
        m.status is True
        for m in session.inputs.monomers
    )


def test_multiple_nonreactants_S_retains_selected_and_discards_others(
    monkeypatch,
):
    detector = NonReactantsDetector()

    m1 = make_monomer(
        1,
        "M1",
        "CCO",
    )

    m2 = make_monomer(
        2,
        "M2",
        "CCN",
    )

    m3 = make_monomer(
        3,
        "M3",
        "CCC",
    )

    session = make_session(
        monomers=[
            m1,
            m2,
            m3,
        ],
        non_reactants=[
            m1,
            m2,
            m3,
        ],
    )

    responses = iter(
        [
            "S",
            "1,3",
        ]
    )

    monkeypatch.setattr(
        "builtins.input",
        lambda prompt: next(responses),
    )

    detector.non_reactant_selection(session)

    status_by_id = {
        m.id: m.status
        for m in session.inputs.monomers
    }

    assert status_by_id == {
        1: True,
        2: False,
        3: True,
    }


def test_selective_mode_accepts_whitespace_and_duplicate_ids(
    monkeypatch,
):
    detector = NonReactantsDetector()

    m1 = make_monomer(
        1,
        "M1",
        "CCO",
    )

    m2 = make_monomer(
        2,
        "M2",
        "CCN",
    )

    session = make_session(
        monomers=[
            m1,
            m2,
        ],
        non_reactants=[
            m1,
            m2,
        ],
    )

    responses = iter(
        [
            "S",
            " 1, 1 ",
        ]
    )

    monkeypatch.setattr(
        "builtins.input",
        lambda prompt: next(responses),
    )

    detector.non_reactant_selection(session)

    status_by_id = {
        m.id: m.status
        for m in session.inputs.monomers
    }

    assert status_by_id == {
        1: True,
        2: False,
    }


def test_selective_mode_reprompts_for_noninteger_input(
    monkeypatch,
):
    detector = NonReactantsDetector()

    m1 = make_monomer(
        1,
        "M1",
        "CCO",
    )

    m2 = make_monomer(
        2,
        "M2",
        "CCN",
    )

    session = make_session(
        monomers=[
            m1,
            m2,
        ],
        non_reactants=[
            m1,
            m2,
        ],
    )

    responses = iter(
        [
            "S",
            "abc",
            "S",
            "1",
        ]
    )

    monkeypatch.setattr(
        "builtins.input",
        lambda prompt: next(responses),
    )

    detector.non_reactant_selection(session)

    status_by_id = {
        m.id: m.status
        for m in session.inputs.monomers
    }

    assert status_by_id == {
        1: True,
        2: False,
    }


def test_selective_mode_reprompts_for_unknown_id(
    monkeypatch,
):
    detector = NonReactantsDetector()

    m1 = make_monomer(
        1,
        "M1",
        "CCO",
    )

    m2 = make_monomer(
        2,
        "M2",
        "CCN",
    )

    session = make_session(
        monomers=[
            m1,
            m2,
        ],
        non_reactants=[
            m1,
            m2,
        ],
    )

    responses = iter(
        [
            "S",
            "99",
            "S",
            "2",
        ]
    )

    monkeypatch.setattr(
        "builtins.input",
        lambda prompt: next(responses),
    )

    detector.non_reactant_selection(session)

    status_by_id = {
        m.id: m.status
        for m in session.inputs.monomers
    }

    assert status_by_id == {
        1: False,
        2: True,
    }


def test_selective_mode_reprompts_for_empty_id_list(
    monkeypatch,
):
    detector = NonReactantsDetector()

    m1 = make_monomer(
        1,
        "M1",
        "CCO",
    )

    m2 = make_monomer(
        2,
        "M2",
        "CCN",
    )

    session = make_session(
        monomers=[
            m1,
            m2,
        ],
        non_reactants=[
            m1,
            m2,
        ],
    )

    responses = iter(
        [
            "S",
            "",
            "S",
            "1",
        ]
    )

    monkeypatch.setattr(
        "builtins.input",
        lambda prompt: next(responses),
    )

    detector.non_reactant_selection(session)

    status_by_id = {
        m.id: m.status
        for m in session.inputs.monomers
    }

    assert status_by_id == {
        1: True,
        2: False,
    }


def test_multiple_nonreactants_reprompts_invalid_main_option(
    monkeypatch,
):
    detector = NonReactantsDetector()

    m1 = make_monomer(
        1,
        "M1",
        "CCO",
    )

    m2 = make_monomer(
        2,
        "M2",
        "CCN",
    )

    session = make_session(
        monomers=[
            m1,
            m2,
        ],
        non_reactants=[
            m1,
            m2,
        ],
    )

    responses = iter(
        [
            "wrong",
            "N",
        ]
    )

    monkeypatch.setattr(
        "builtins.input",
        lambda prompt: next(responses),
    )

    detector.non_reactant_selection(session)

    assert all(
        m.status is False
        for m in session.inputs.monomers
    )


def test_selection_preserves_reactive_monomers(
    monkeypatch,
):
    """
    Discarding non-reactants must not accidentally change a monomer that
    participates in a reaction.
    """
    detector = NonReactantsDetector()

    nonreactant = make_monomer(
        1,
        "solvent",
        "CCO",
    )

    reactive = make_monomer(
        2,
        "reactive",
        "C=C",
    )

    session = make_session(
        monomers=[
            nonreactant,
            reactive,
        ],
        non_reactants=[
            nonreactant,
        ],
    )

    monkeypatch.setattr(
        "builtins.input",
        lambda prompt: "N",
    )

    detector.non_reactant_selection(session)

    status_by_id = {
        m.id: m.status
        for m in session.inputs.monomers
    }

    assert status_by_id == {
        1: False,
        2: True,
    }