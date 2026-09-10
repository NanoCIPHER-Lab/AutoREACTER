from types import SimpleNamespace

import pytest
from PIL import Image
from rdkit import Chem

from AutoREACTER.detectors.functional_groups_detector import (
    FunctionalGroupInfo,
    MonomerRole,
)
from AutoREACTER.detectors.reaction_detector import (
    EmptyReactionListError,
    ReactionDetector,
    ReactionInstance,
    SMARTSerror,
)


# =============================================================================
# Helpers
# =============================================================================

def make_fg(name: str) -> FunctionalGroupInfo:
    return FunctionalGroupInfo(
        functionality_type="test",
        fg_name=name,
        fg_smarts_1="[*]",
        fg_count_1=1,
    )


def make_role(
    name: str,
    smiles: str,
    fg_names,
    *,
    is_looped=False,
) -> MonomerRole:
    return MonomerRole(
        name=name,
        smiles=smiles,
        functionalities=tuple(make_fg(fg) for fg in fg_names),
        is_looped=is_looped,
    )


def make_session(monomer_roles=None, reaction_instances=None):
    return SimpleNamespace(
        monomer_roles=monomer_roles or [],
        reaction_instances=reaction_instances or [],
    )


def homo_reaction(
    reactant="vinyl",
    *,
    reaction_name="vinyl_polymerization",
):
    return {
        reaction_name: {
            "reactant_1": reactant,
            "reactant_2": None,
            "same_reactants": True,
            "reaction": "[*:1]>>[*:1]",
            "delete_atom": False,
            "reference": {"test": "reference"},
        }
    }


def co_reaction(
    reactant_1="A",
    reactant_2="B",
    *,
    reaction_name="A_B_polymerization",
):
    return {
        reaction_name: {
            "reactant_1": reactant_1,
            "reactant_2": reactant_2,
            "same_reactants": False,
            "reaction": "[*:1].[*:2]>>[*:1]-[*:2]",
            "delete_atom": True,
            "reference": {"test": "reference"},
        }
    }


# =============================================================================
# ReactionInstance
# =============================================================================

def test_reaction_instance_stores_required_data():
    monomer = make_role("styrene", "C=Cc1ccccc1", ["vinyl"])
    fg = monomer.functionalities[0]

    instance = ReactionInstance(
        reaction_name="vinyl_polymerization",
        reaction_smarts="[*:1]>>[*:1]",
        delete_atom=False,
        references={"paper": "test"},
        same_reactants=True,
        monomer_1=monomer,
        functional_group_1=fg,
    )

    assert instance.reaction_name == "vinyl_polymerization"
    assert instance.reaction_smarts == "[*:1]>>[*:1]"
    assert instance.delete_atom is False
    assert instance.references == {"paper": "test"}
    assert instance.same_reactants is True
    assert instance.monomer_1 is monomer
    assert instance.functional_group_1 is fg
    assert instance.monomer_2 is None
    assert instance.functional_group_2 is None


# =============================================================================
# _matching_fgs()
# =============================================================================

def test_matching_fgs_returns_matching_functionality():
    detector = ReactionDetector()

    monomer = make_role(
        "test",
        "CC",
        ["vinyl", "diol", "primary_amine"],
    )

    result = detector._matching_fgs(monomer, "diol")

    assert len(result) == 1
    assert result[0].fg_name == "diol"


def test_matching_fgs_returns_all_matching_entries():
    detector = ReactionDetector()

    fg1 = make_fg("vinyl")
    fg2 = make_fg("vinyl")

    monomer = MonomerRole(
        name="test",
        smiles="C=C",
        functionalities=(fg1, fg2),
    )

    result = detector._matching_fgs(monomer, "vinyl")

    assert result == [fg1, fg2]


def test_matching_fgs_returns_empty_list_when_no_match():
    detector = ReactionDetector()

    monomer = make_role("test", "CCO", ["diol"])

    assert detector._matching_fgs(monomer, "vinyl") == []


# =============================================================================
# _seen_pair_key()
# =============================================================================

def test_seen_pair_key_for_homopolymerization():
    detector = ReactionDetector()

    monomer = make_role("styrene", "C=Cc1ccccc1", ["vinyl"])
    fg = monomer.functionalities[0]

    key = detector._seen_pair_key(
        "vinyl_polymerization",
        monomer,
        fg,
    )

    assert key == (
        "vinyl_polymerization",
        "C=Cc1ccccc1",
        "vinyl",
    )


def test_seen_pair_key_for_copolymerization_is_order_independent():
    detector = ReactionDetector()

    monomer_a = make_role("A", "CCO", ["A"])
    monomer_b = make_role("B", "CCN", ["B"])

    key_ab = detector._seen_pair_key(
        "reaction",
        monomer_a,
        monomer_a.functionalities[0],
        monomer_b,
        monomer_b.functionalities[0],
    )

    key_ba = detector._seen_pair_key(
        "reaction",
        monomer_b,
        monomer_b.functionalities[0],
        monomer_a,
        monomer_a.functionalities[0],
    )

    assert key_ab == key_ba


def test_seen_pair_key_distinguishes_reaction_names():
    detector = ReactionDetector()

    monomer_a = make_role("A", "CCO", ["A"])
    monomer_b = make_role("B", "CCN", ["B"])

    key_1 = detector._seen_pair_key(
        "reaction_1",
        monomer_a,
        monomer_a.functionalities[0],
        monomer_b,
        monomer_b.functionalities[0],
    )

    key_2 = detector._seen_pair_key(
        "reaction_2",
        monomer_a,
        monomer_a.functionalities[0],
        monomer_b,
        monomer_b.functionalities[0],
    )

    assert key_1 != key_2


# =============================================================================
# reaction_detector() - homopolymerization
# =============================================================================

def test_homopolymerization_is_detected():
    detector = ReactionDetector()
    detector.reactions = homo_reaction()

    monomer = make_role(
        "styrene",
        "C=Cc1ccccc1",
        ["vinyl"],
    )

    session = make_session([monomer])

    result = detector.reaction_detector(session)

    assert result is None
    assert len(session.reaction_instances) == 1

    rxn = session.reaction_instances[0]

    assert rxn.reaction_name == "vinyl_polymerization"
    assert rxn.same_reactants is True
    assert rxn.monomer_1 is monomer
    assert rxn.functional_group_1.fg_name == "vinyl"
    assert rxn.monomer_2 is None
    assert rxn.functional_group_2 is None


def test_homopolymerization_preserves_reaction_metadata():
    detector = ReactionDetector()

    detector.reactions = {
        "test_reaction": {
            "reactant_1": "vinyl",
            "reactant_2": None,
            "same_reactants": True,
            "reaction": "TEST_SMARTS",
            "delete_atom": True,
            "reference": {"doi": "123"},
        }
    }

    session = make_session(
        [make_role("M", "C=C", ["vinyl"])]
    )

    detector.reaction_detector(session)

    rxn = session.reaction_instances[0]

    assert rxn.reaction_smarts == "TEST_SMARTS"
    assert rxn.delete_atom is True
    assert rxn.references == {"doi": "123"}


def test_duplicate_homo_functionalities_do_not_create_duplicate_reactions():
    detector = ReactionDetector()
    detector.reactions = homo_reaction()

    monomer = make_role(
        "test",
        "C=C",
        ["vinyl", "vinyl"],
    )

    session = make_session([monomer])

    detector.reaction_detector(session)

    assert len(session.reaction_instances) == 1


def test_multiple_different_homo_monomers_create_multiple_instances():
    detector = ReactionDetector()
    detector.reactions = homo_reaction()

    monomer_1 = make_role("M1", "C=C", ["vinyl"])
    monomer_2 = make_role("M2", "C=CC", ["vinyl"])

    session = make_session([monomer_1, monomer_2])

    detector.reaction_detector(session)

    assert len(session.reaction_instances) == 2


# =============================================================================
# reaction_detector() - copolymerization
# =============================================================================

def test_copolymerization_detects_A_plus_B():
    detector = ReactionDetector()
    detector.reactions = co_reaction()

    monomer_a = make_role("A", "CCO", ["A"])
    monomer_b = make_role("B", "CCN", ["B"])

    session = make_session([monomer_a, monomer_b])

    detector.reaction_detector(session)

    assert len(session.reaction_instances) == 1

    rxn = session.reaction_instances[0]

    assert rxn.monomer_1 is monomer_a
    assert rxn.monomer_2 is monomer_b
    assert rxn.functional_group_1.fg_name == "A"
    assert rxn.functional_group_2.fg_name == "B"
    assert rxn.same_reactants is False


def test_copolymerization_reverse_scan_does_not_duplicate_pair():
    detector = ReactionDetector()
    detector.reactions = co_reaction()

    monomer_a = make_role("A", "CCO", ["A"])
    monomer_b = make_role("B", "CCN", ["B"])

    session = make_session(
        [monomer_b, monomer_a]
    )

    detector.reaction_detector(session)

    assert len(session.reaction_instances) == 1


def test_unrelated_monomer_is_ignored():
    detector = ReactionDetector()
    detector.reactions = co_reaction()

    monomer_a = make_role("A", "CCO", ["A"])
    monomer_b = make_role("B", "CCN", ["B"])
    unrelated = make_role("X", "CCC", ["X"])

    session = make_session(
        [monomer_a, monomer_b, unrelated]
    )

    detector.reaction_detector(session)

    assert len(session.reaction_instances) == 1


# =============================================================================
# Same molecule carrying A + B
# =============================================================================

def test_single_monomer_with_both_functional_groups_can_react():
    """
    A heterobifunctional monomer containing both A and B should be detected
    by the special same-monomer A+B branch.
    """
    detector = ReactionDetector()
    detector.reactions = co_reaction()

    monomer = make_role(
        "AB",
        "NCC(=O)O",
        ["A", "B"],
    )

    session = make_session([monomer])

    detector.reaction_detector(session)

    assert len(session.reaction_instances) == 1

    rxn = session.reaction_instances[0]

    assert rxn.monomer_1 is monomer
    assert rxn.monomer_2 is monomer
    assert rxn.functional_group_1.fg_name == "A"
    assert rxn.functional_group_2.fg_name == "B"


def test_same_monomer_AB_reaction_not_duplicated():
    detector = ReactionDetector()
    detector.reactions = co_reaction()

    monomer = make_role(
        "AB",
        "NCC(=O)O",
        ["A", "B"],
    )

    session = make_session([monomer])

    detector.reaction_detector(session)

    assert len(session.reaction_instances) == 1


# =============================================================================
# Empty detection
# =============================================================================

def test_no_matching_reaction_raises_empty_reaction_error():
    detector = ReactionDetector()
    detector.reactions = co_reaction()

    session = make_session(
        [make_role("X", "CCC", ["X"])]
    )

    with pytest.raises(
        EmptyReactionListError,
        match="No reaction instances found",
    ):
        detector.reaction_detector(session)


def test_empty_monomer_role_list_raises_empty_reaction_error():
    detector = ReactionDetector()
    detector.reactions = homo_reaction()

    session = make_session([])

    with pytest.raises(EmptyReactionListError):
        detector.reaction_detector(session)


# =============================================================================
# index_based_reaction_detector() - homo
# =============================================================================

def test_index_homo_fresh_role_is_detected():
    detector = ReactionDetector()
    detector.reactions = homo_reaction()

    monomer = make_role(
        "M",
        "C=C",
        ["vinyl"],
        is_looped=False,
    )

    result = detector.index_based_reaction_detector([monomer])

    assert len(result) == 1


def test_index_homo_looped_role_is_skipped():
    detector = ReactionDetector()
    detector.reactions = homo_reaction()

    monomer = make_role(
        "M",
        "C=C",
        ["vinyl"],
        is_looped=True,
    )

    result = detector.index_based_reaction_detector([monomer])

    assert result == []


# =============================================================================
# index_based_reaction_detector() - co
# =============================================================================

def test_index_co_two_fresh_roles_are_detected():
    detector = ReactionDetector()
    detector.reactions = co_reaction()

    a = make_role("A", "CCO", ["A"])
    b = make_role("B", "CCN", ["B"])

    result = detector.index_based_reaction_detector([a, b])

    assert len(result) == 1


@pytest.mark.parametrize(
    "a_looped,b_looped",
    [
        (True, False),
        (False, True),
    ],
)
def test_index_co_pair_is_kept_when_only_one_role_is_looped(
    a_looped,
    b_looped,
):
    detector = ReactionDetector()
    detector.reactions = co_reaction()

    a = make_role(
        "A",
        "CCO",
        ["A"],
        is_looped=a_looped,
    )

    b = make_role(
        "B",
        "CCN",
        ["B"],
        is_looped=b_looped,
    )

    result = detector.index_based_reaction_detector([a, b])

    assert len(result) == 1


def test_index_co_pair_is_skipped_when_both_roles_are_looped():
    detector = ReactionDetector()
    detector.reactions = co_reaction()

    a = make_role(
        "A",
        "CCO",
        ["A"],
        is_looped=True,
    )

    b = make_role(
        "B",
        "CCN",
        ["B"],
        is_looped=True,
    )

    result = detector.index_based_reaction_detector([a, b])

    assert result == []


def test_index_same_monomer_AB_fresh_role_is_detected():
    detector = ReactionDetector()
    detector.reactions = co_reaction()

    monomer = make_role(
        "AB",
        "NCC(=O)O",
        ["A", "B"],
        is_looped=False,
    )

    result = detector.index_based_reaction_detector([monomer])

    assert len(result) == 1

    assert result[0].monomer_1 is monomer
    assert result[0].monomer_2 is monomer


def test_index_same_monomer_AB_looped_role_is_skipped():
    detector = ReactionDetector()
    detector.reactions = co_reaction()

    monomer = make_role(
        "AB",
        "NCC(=O)O",
        ["A", "B"],
        is_looped=True,
    )

    result = detector.index_based_reaction_detector([monomer])

    assert result == []


def test_index_detector_deduplicates_reaction_pairs():
    detector = ReactionDetector()
    detector.reactions = co_reaction()

    a = make_role("A", "CCO", ["A", "A"])
    b = make_role("B", "CCN", ["B", "B"])

    result = detector.index_based_reaction_detector([a, b])

    assert len(result) == 1


# =============================================================================
# create_reaction_image()
# =============================================================================

def test_create_reaction_image_returns_generated_image(monkeypatch):
    detector = ReactionDetector()

    product = Chem.MolFromSmiles("CO")

    class FakeEngine:
        def RunReactants(self, reactants):
            return ((product,),)

    monkeypatch.setattr(
        "AutoREACTER.detectors.reaction_detector."
        "rdChemReactions.ReactionFromSmarts",
        lambda smarts: FakeEngine(),
    )

    fake_image = object()

    monkeypatch.setattr(
        "AutoREACTER.detectors.reaction_detector."
        "Draw.ReactionToImage",
        lambda reaction, subImgSize: fake_image,
    )

    result = detector.create_reaction_image(
        "C",
        "O",
        "fake_smarts",
        "test_reaction",
    )

    assert result is fake_image


def test_create_reaction_image_tries_reverse_reactant_order(monkeypatch):
    detector = ReactionDetector()

    product = Chem.MolFromSmiles("CO")

    class FakeEngine:
        def __init__(self):
            self.calls = 0

        def RunReactants(self, reactants):
            self.calls += 1

            if self.calls == 1:
                return ()

            return ((product,),)

    engine = FakeEngine()

    monkeypatch.setattr(
        "AutoREACTER.detectors.reaction_detector."
        "rdChemReactions.ReactionFromSmarts",
        lambda smarts: engine,
    )

    monkeypatch.setattr(
        "AutoREACTER.detectors.reaction_detector."
        "Draw.ReactionToImage",
        lambda reaction, subImgSize: "image",
    )

    result = detector.create_reaction_image(
        "C",
        "O",
        "fake_smarts",
        "test",
    )

    assert result == "image"
    assert engine.calls == 2


def test_create_reaction_image_raises_when_both_orders_fail(
    monkeypatch,
):
    detector = ReactionDetector()

    class FakeEngine:
        def RunReactants(self, reactants):
            return ()

    monkeypatch.setattr(
        "AutoREACTER.detectors.reaction_detector."
        "rdChemReactions.ReactionFromSmarts",
        lambda smarts: FakeEngine(),
    )

    with pytest.raises(
        SMARTSerror,
        match="failed",
    ):
        detector.create_reaction_image(
            "C",
            "O",
            "bad_reaction",
            "test",
        )


# =============================================================================
# available_reaction_image_grid()
# =============================================================================

def test_available_reaction_image_grid_stacks_images_vertically(
    monkeypatch,
):
    detector = ReactionDetector()

    m1 = make_role("A", "CCO", ["A"])
    m2 = make_role("B", "CCN", ["B"])

    fg1 = m1.functionalities[0]
    fg2 = m2.functionalities[0]

    reaction_1 = ReactionInstance(
        reaction_name="R1",
        reaction_smarts="x",
        delete_atom=False,
        references={},
        same_reactants=False,
        monomer_1=m1,
        functional_group_1=fg1,
        monomer_2=m2,
        functional_group_2=fg2,
    )

    reaction_2 = ReactionInstance(
        reaction_name="R2",
        reaction_smarts="x",
        delete_atom=False,
        references={},
        same_reactants=False,
        monomer_1=m1,
        functional_group_1=fg1,
        monomer_2=m2,
        functional_group_2=fg2,
    )

    monkeypatch.setattr(
        detector,
        "create_reaction_image",
        lambda *args, **kwargs: Image.new(
            "RGB",
            (100, 50),
            "white",
        ),
    )

    session = make_session(
        reaction_instances=[reaction_1, reaction_2]
    )

    result = detector.available_reaction_image_grid(session)

    assert isinstance(result, Image.Image)
    assert result.size == (180, 100)


def test_available_reaction_image_grid_returns_none_when_empty():
    detector = ReactionDetector()

    session = make_session(reaction_instances=[])

    result = detector.available_reaction_image_grid(session)

    assert result is None


def test_available_reaction_image_grid_skips_failed_visualizations(
    monkeypatch,
):
    detector = ReactionDetector()

    monomer = make_role("M", "C=C", ["vinyl"])

    reaction = ReactionInstance(
        reaction_name="R",
        reaction_smarts="bad",
        delete_atom=False,
        references={},
        same_reactants=True,
        monomer_1=monomer,
        functional_group_1=monomer.functionalities[0],
    )

    def fail(*args, **kwargs):
        raise SMARTSerror("bad SMARTS")

    monkeypatch.setattr(
        detector,
        "create_reaction_image",
        fail,
    )

    session = make_session(reaction_instances=[reaction])

    result = detector.available_reaction_image_grid(session)

    assert result is None


# =============================================================================
# reaction_selection()
# =============================================================================

def test_reaction_selection_raises_for_empty_list():
    detector = ReactionDetector()

    session = make_session(reaction_instances=[])

    with pytest.raises(
        EmptyReactionListError,
        match="Cannot proceed",
    ):
        detector.reaction_selection(session)


def test_reaction_selection_automatically_keeps_single_reaction():
    detector = ReactionDetector()

    monomer = make_role("M", "C=C", ["vinyl"])

    reaction = ReactionInstance(
        reaction_name="R",
        reaction_smarts="x",
        delete_atom=False,
        references={},
        same_reactants=True,
        monomer_1=monomer,
        functional_group_1=monomer.functionalities[0],
    )

    session = make_session(reaction_instances=[reaction])

    detector.reaction_selection(session)

    assert session.reaction_instances == [reaction]


def test_reaction_selection_selects_requested_reactions(monkeypatch):
    detector = ReactionDetector()

    monomer = make_role("M", "C=C", ["vinyl"])
    fg = monomer.functionalities[0]

    reactions = [
        ReactionInstance(
            reaction_name=f"R{i}",
            reaction_smarts="x",
            delete_atom=False,
            references={},
            same_reactants=True,
            monomer_1=monomer,
            functional_group_1=fg,
        )
        for i in range(1, 4)
    ]

    session = make_session(reaction_instances=reactions)

    monkeypatch.setattr(
        "builtins.input",
        lambda prompt: "3,1,3",
    )

    detector.reaction_selection(session)

    assert session.reaction_instances == [
        reactions[0],
        reactions[2],
    ]


def test_reaction_selection_reprompts_after_invalid_input(
    monkeypatch,
):
    detector = ReactionDetector()

    monomer = make_role("M", "C=C", ["vinyl"])
    fg = monomer.functionalities[0]

    reactions = [
        ReactionInstance(
            reaction_name=f"R{i}",
            reaction_smarts="x",
            delete_atom=False,
            references={},
            same_reactants=True,
            monomer_1=monomer,
            functional_group_1=fg,
        )
        for i in range(1, 3)
    ]

    session = make_session(reaction_instances=reactions)

    responses = iter(
        [
            "",
            "abc",
            "9",
            "2",
        ]
    )

    monkeypatch.setattr(
        "builtins.input",
        lambda prompt: next(responses),
    )

    detector.reaction_selection(session)

    assert session.reaction_instances == [
        reactions[1]
    ]