from types import SimpleNamespace

import pytest
from rdkit import Chem

import AutoREACTER.reaction_preparation.reaction_processor.reaction_progression as progression_module
from AutoREACTER.reaction_preparation.reaction_processor.reaction_progression import (
    MAX_LOOP,
    MonomerRoleforIndexBasedFGDetection,
    ReactionProgression,
    ReactionProgressionSession,
)


# =============================================================================
# Helpers
# =============================================================================


def mol(smiles: str) -> Chem.Mol:
    molecule = Chem.MolFromSmiles(smiles)

    assert molecule is not None

    return molecule


def make_role(
    smiles="CC",
    name="monomer",
    *,
    is_monomer=True,
    is_looped=False,
    rdkit_mol=None,
):
    return SimpleNamespace(
        smiles=smiles,
        name=name,
        is_monomer=is_monomer,
        is_looped=is_looped,
        rdkit_mol=rdkit_mol,
    )


def make_reaction(
    reaction_id=1,
    *,
    reactant_smiles="CC",
    product_smiles="CC",
    activity_stats=True,
    delete_atom=False,
    template_mapping=None,
    product_to_reactant_mapping=None,
    is_radical=False,
):
    reactant = mol(reactant_smiles)
    product = mol(product_smiles)

    if template_mapping is None:
        template_mapping = {
            idx: idx
            for idx in range(
                min(
                    reactant.GetNumAtoms(),
                    product.GetNumAtoms(),
                )
            )
        }

    if product_to_reactant_mapping is None:
        product_to_reactant_mapping = {
            product_idx: reactant_idx
            for reactant_idx, product_idx
            in template_mapping.items()
        }

    return SimpleNamespace(
        reaction_id=reaction_id,
        reactant_combined_RDmol=reactant,
        product_combined_RDmol=product,
        template_reactant_to_product_mapping=template_mapping,
        product_to_reactant_mapping=product_to_reactant_mapping,
        activity_stats=activity_stats,
        delete_atom=delete_atom,
        is_radical=is_radical,
        radical_atom_idxs=(),
    )


def make_session(
    *,
    monomer_roles=None,
    reaction_metadata=None,
    deep_search=True,
    reaction_iteration_depth=5,
):
    return SimpleNamespace(
        monomer_roles=list(
            monomer_roles or []
        ),
        reaction_metadata=list(
            reaction_metadata or []
        ),
        reaction_progression_session=None,
        inputs=SimpleNamespace(
            deep_search=deep_search,
            reaction_iteration_depth=reaction_iteration_depth,
        ),
    )


def make_progression(
    session,
    preparer=None,
):
    """
    Build a ReactionProgression without invoking its constructor.

    Most method-level tests do not need real detector construction or the
    user-facing beta warning.
    """
    progression = ReactionProgression.__new__(
        ReactionProgression
    )

    progression.session = session
    progression.preparer = preparer

    if (
        session.reaction_progression_session
        is None
    ):
        session.reaction_progression_session = (
            ReactionProgressionSession()
        )

    progression.fg_detector = (
        SimpleNamespace()
    )

    progression.rxn_detector = (
        SimpleNamespace()
    )

    progression.deduplication_detector = (
        SimpleNamespace()
    )

    return progression


# =============================================================================
# Dataclasses
# =============================================================================


def test_index_based_role_required_fields():
    role = (
        MonomerRoleforIndexBasedFGDetection(
            smiles="CC",
            name="new_1",
            indexes_in_template=[0, 1],
        )
    )

    assert role.smiles == "CC"
    assert role.name == "new_1"

    assert role.indexes_in_template == [
        0,
        1,
    ]


def test_index_based_role_defaults():
    role = (
        MonomerRoleforIndexBasedFGDetection(
            smiles="C",
            name="new_1",
            indexes_in_template=[],
        )
    )

    assert role.is_monomer is False
    assert role.is_looped is False
    assert role.rdkit_mol is None


def test_index_based_role_uses_slots():
    role = (
        MonomerRoleforIndexBasedFGDetection(
            smiles="C",
            name="new_1",
            indexes_in_template=[],
        )
    )

    with pytest.raises(
        AttributeError
    ):
        role.unexpected = True


def test_progression_session_defaults():
    session = ReactionProgressionSession()

    assert session.monomer_roles == []
    assert session.iteration == 0


def test_progression_session_lists_are_independent():
    first = ReactionProgressionSession()
    second = ReactionProgressionSession()

    first.monomer_roles.append(
        object()
    )

    assert second.monomer_roles == []


def test_progression_session_uses_slots():
    session = ReactionProgressionSession()

    with pytest.raises(
        AttributeError
    ):
        session.unexpected = True


def test_max_loop_constant():
    assert MAX_LOOP == 5


# =============================================================================
# Constructor
# =============================================================================


def test_constructor_attaches_progression_session(
    monkeypatch,
):
    session = make_session()

    monkeypatch.setattr(
        progression_module,
        "FunctionalGroupsDetector",
        lambda: "fg-detector",
    )

    monkeypatch.setattr(
        progression_module,
        "ReactionDetector",
        lambda: "rxn-detector",
    )

    monkeypatch.setattr(
        progression_module,
        "DeduplicationDetector",
        lambda: "dedup-detector",
    )

    warnings = []

    monkeypatch.setattr(
        progression_module,
        "print_warning",
        lambda: warnings.append(True),
    )

    preparer = object()

    progression = ReactionProgression(
        session,
        preparer=preparer,
    )

    assert progression.session is session
    assert progression.preparer is preparer

    assert isinstance(
        session.reaction_progression_session,
        ReactionProgressionSession,
    )

    assert progression.fg_detector == (
        "fg-detector"
    )

    assert progression.rxn_detector == (
        "rxn-detector"
    )

    assert (
        progression.deduplication_detector
        == "dedup-detector"
    )

    assert warnings == [
        True
    ]


def test_constructor_replaces_existing_progression_session(
    monkeypatch,
):
    old = ReactionProgressionSession(
        iteration=99
    )

    session = make_session()
    session.reaction_progression_session = old

    monkeypatch.setattr(
        progression_module,
        "FunctionalGroupsDetector",
        lambda: object(),
    )

    monkeypatch.setattr(
        progression_module,
        "ReactionDetector",
        lambda: object(),
    )

    monkeypatch.setattr(
        progression_module,
        "DeduplicationDetector",
        lambda: object(),
    )

    monkeypatch.setattr(
        progression_module,
        "print_warning",
        lambda: None,
    )

    ReactionProgression(session)

    assert (
        session.reaction_progression_session
        is not old
    )

    assert (
        session
        .reaction_progression_session
        .iteration
        == 0
    )


# =============================================================================
# Thin preparation wrapper
# =============================================================================


def test_index_based_reaction_preparation_delegates():
    session = make_session()

    calls = []

    class FakePreparer:
        def _prepare_reactions_stage(
            self,
            reaction_instances,
            loop=False,
        ):
            calls.append(
                (
                    reaction_instances,
                    loop,
                )
            )

            return [
                "prepared"
            ]

    progression = make_progression(
        session,
        preparer=FakePreparer(),
    )

    instances = [
        object(),
        object(),
    ]

    result = (
        progression
        ._index_based_reaction_preparation(
            instances
        )
    )

    assert result == [
        "prepared"
    ]

    assert calls == [
        (
            instances,
            True,
        )
    ]


# =============================================================================
# Basic bookkeeping helpers
# =============================================================================


def test_store_reactions_updates_session():
    session = make_session()

    progression = make_progression(
        session
    )

    reactions = [
        object(),
        object(),
    ]

    result = progression._store_reactions(
        reactions
    )

    assert result is reactions

    assert (
        session.reaction_metadata
        is reactions
    )


def test_count_active_reactions():
    progression = make_progression(
        make_session()
    )

    reactions = [
        SimpleNamespace(
            activity_stats=True
        ),
        SimpleNamespace(
            activity_stats=False
        ),
        SimpleNamespace(
            activity_stats=1
        ),
        SimpleNamespace(
            activity_stats=None
        ),
    ]

    assert (
        progression
        ._count_active_reactions(
            reactions
        )
        == 2
    )


def test_count_active_reactions_empty_list():
    progression = make_progression(
        make_session()
    )

    assert (
        progression
        ._count_active_reactions([])
        == 0
    )


def test_set_is_looped_flag_marks_all_roles():
    roles = [
        make_role(
            is_looped=False
        ),
        make_role(
            is_looped=False
        ),
    ]

    progression = make_progression(
        make_session()
    )

    progression._set_is_looped_flag(
        roles
    )

    assert all(
        role.is_looped
        for role in roles
    )


def test_set_is_looped_flag_empty_list():
    progression = make_progression(
        make_session()
    )

    progression._set_is_looped_flag(
        []
    )


# =============================================================================
# Loop break logic
# =============================================================================


@pytest.mark.parametrize(
    "before, after",
    [
        (1, 1),
        (2, 1),
        (5, 0),
    ],
)
def test_loop_break_condition_true_when_pool_does_not_grow(
    before,
    after,
):
    progression = make_progression(
        make_session()
    )

    assert (
        progression
        ._loop_break_condition(
            before,
            after,
        )
        is True
    )


@pytest.mark.parametrize(
    "before, after",
    [
        (0, 1),
        (1, 2),
        (10, 11),
    ],
)
def test_loop_break_condition_false_when_pool_grows(
    before,
    after,
):
    progression = make_progression(
        make_session()
    )

    assert (
        progression
        ._loop_break_condition(
            before,
            after,
        )
        is False
    )


def test_loop_break_condition_prints_reason(
    capsys,
):
    progression = make_progression(
        make_session()
    )

    progression._loop_break_condition(
        3,
        3,
    )

    output = capsys.readouterr().out

    assert (
        "pool did not grow"
        in output
    )

    assert "before=3" in output
    assert "after=3" in output


# =============================================================================
# Monomer population
# =============================================================================


def test_smiles_to_rdkit_mol_valid():
    progression = make_progression(
        make_session()
    )

    result = (
        progression
        ._smiles_to_rdkit_mol(
            "CCO"
        )
    )

    assert isinstance(
        result,
        Chem.Mol,
    )


def test_smiles_to_rdkit_mol_invalid_returns_none():
    progression = make_progression(
        make_session()
    )

    result = (
        progression
        ._smiles_to_rdkit_mol(
            "not-a-smiles"
        )
    )

    assert result is None


def test_populate_monomer_roles_only_populates_monomers():
    monomer = make_role(
        smiles="CC",
        is_monomer=True,
    )

    generated = make_role(
        smiles="CO",
        is_monomer=False,
    )

    session = make_session(
        monomer_roles=[
            monomer,
            generated,
        ]
    )

    progression = make_progression(
        session
    )

    progression._populate_monomer_roles()

    assert isinstance(
        monomer.rdkit_mol,
        Chem.Mol,
    )

    assert generated.rdkit_mol is None


def test_populate_monomer_roles_invalid_smiles_sets_none():
    role = make_role(
        smiles="not-a-smiles",
        is_monomer=True,
        rdkit_mol=object(),
    )

    session = make_session(
        monomer_roles=[
            role
        ]
    )

    progression = make_progression(
        session
    )

    progression._populate_monomer_roles()

    assert role.rdkit_mol is None


# =============================================================================
# Product cleanup
# =============================================================================


def test_clean_product_returns_copy():
    progression = make_progression(
        make_session()
    )

    molecule = mol("CC")

    cleaned = (
        progression
        ._clean_product(
            molecule
        )
    )

    assert cleaned is not molecule

    assert (
        Chem.MolToSmiles(cleaned)
        == Chem.MolToSmiles(molecule)
    )


def test_clean_product_removes_atom_map_numbers():
    progression = make_progression(
        make_session()
    )

    molecule = mol("CC")

    molecule.GetAtomWithIdx(
        0
    ).SetAtomMapNum(123)

    cleaned = (
        progression
        ._clean_product(
            molecule
        )
    )

    assert (
        cleaned
        .GetAtomWithIdx(0)
        .GetAtomMapNum()
        == 0
    )

    # Original must remain unchanged.
    assert (
        molecule
        .GetAtomWithIdx(0)
        .GetAtomMapNum()
        == 123
    )


def test_clean_product_removes_isotopes():
    progression = make_progression(
        make_session()
    )

    molecule = mol("CC")

    molecule.GetAtomWithIdx(
        0
    ).SetIsotope(13)

    cleaned = (
        progression
        ._clean_product(
            molecule
        )
    )

    assert (
        cleaned
        .GetAtomWithIdx(0)
        .GetIsotope()
        == 0
    )

    assert (
        molecule
        .GetAtomWithIdx(0)
        .GetIsotope()
        == 13
    )


def test_clean_product_removes_tracking_properties():
    progression = make_progression(
        make_session()
    )

    molecule = mol("C")

    atom = molecule.GetAtomWithIdx(
        0
    )

    atom.SetIntProp(
        "old_mapno",
        1,
    )

    atom.SetIntProp(
        "react_atom_idx",
        2,
    )

    cleaned = (
        progression
        ._clean_product(
            molecule
        )
    )

    cleaned_atom = (
        cleaned.GetAtomWithIdx(0)
    )

    assert not cleaned_atom.HasProp(
        "old_mapno"
    )

    assert not cleaned_atom.HasProp(
        "react_atom_idx"
    )


# =============================================================================
# Product SMILES
# =============================================================================


def test_get_product_smiles_returns_canonical_smiles():
    progression = make_progression(
        make_session()
    )

    result = (
        progression
        ._get_product_smiles(
            mol("OCC")
        )
    )

    assert result == Chem.MolToSmiles(
        mol("CCO")
    )


def test_get_product_smiles_removes_mapping_artifacts():
    progression = make_progression(
        make_session()
    )

    molecule = mol("CC")

    molecule.GetAtomWithIdx(
        0
    ).SetAtomMapNum(999)

    result = (
        progression
        ._get_product_smiles(
            molecule
        )
    )

    assert ":" not in result


# =============================================================================
# Product index handling
# =============================================================================


def test_get_product_idxs_single_fragment_preserves_indexes():
    progression = make_progression(
        make_session()
    )

    molecule = mol("CCC")

    mapping = {
        10: 0,
        11: 2,
    }

    idxs, result_mol = (
        progression
        ._get_product_idxs(
            mapping,
            molecule,
        )
    )

    assert idxs == [
        0,
        2,
    ]

    assert (
        Chem.MolToSmiles(
            result_mol
        )
        == Chem.MolToSmiles(
            molecule
        )
    )


def test_get_product_idxs_returns_molecule_copy():
    progression = make_progression(
        make_session()
    )

    molecule = mol("CC")

    _, result = (
        progression
        ._get_product_idxs(
            {
                0: 0,
            },
            molecule,
        )
    )

    assert result is not molecule


def test_keep_largest_fragment_uses_heavy_atom_count():
    progression = make_progression(
        make_session()
    )

    molecule = mol(
        "CCCO.CC"
    )

    largest, _ = (
        progression
        ._keep_largest_fragment(
            molecule,
            [],
        )
    )

    assert (
        largest.GetNumHeavyAtoms()
        == 4
    )


def test_keep_largest_fragment_remaps_product_indices():
    progression = make_progression(
        make_session()
    )

    # First fragment atoms: 0,1
    # Second/larger fragment atoms: 2,3,4
    molecule = mol(
        "CC.CCC"
    )

    largest, remapped = (
        progression
        ._keep_largest_fragment(
            molecule,
            [
                2,
                4,
            ],
        )
    )

    assert (
        largest.GetNumHeavyAtoms()
        == 3
    )

    assert remapped == [
        0,
        2,
    ]


def test_keep_largest_fragment_drops_indices_from_removed_fragment():
    progression = make_progression(
        make_session()
    )

    molecule = mol(
        "CC.CCC"
    )

    _, remapped = (
        progression
        ._keep_largest_fragment(
            molecule,
            [
                0,
                2,
                3,
            ],
        )
    )

    assert remapped == [
        0,
        1,
    ]


def test_get_product_idxs_multifragment_keeps_largest_fragment():
    progression = make_progression(
        make_session()
    )

    molecule = mol(
        "CC.CCC"
    )

    idxs, result = (
        progression
        ._get_product_idxs(
            {
                10: 2,
                11: 4,
            },
            molecule,
        )
    )

    assert (
        result.GetNumHeavyAtoms()
        == 3
    )

    assert idxs == [
        0,
        2,
    ]


# =============================================================================
# Radical metadata
# =============================================================================


def test_set_reaction_radical_metadata_no_radicals():
    progression = make_progression(
        make_session()
    )

    reaction = make_reaction(
        product_smiles="CC",
    )

    progression._set_reaction_radical_metadata(
        reaction,
        mol("CC"),
    )

    assert reaction.is_radical is False

    assert (
        reaction.radical_atom_idxs
        == ()
    )


def test_set_reaction_radical_metadata_maps_product_idx_to_reactant_idx():
    progression = make_progression(
        make_session()
    )

    reaction = make_reaction(
        product_smiles="[CH3]",
        template_mapping={
            5: 0,
        },
        product_to_reactant_mapping={
            0: 5,
        },
    )

    radical_product = mol(
        "[CH3]"
    )

    progression._set_reaction_radical_metadata(
        reaction,
        radical_product,
    )

    assert reaction.is_radical is True

    assert (
        reaction.radical_atom_idxs
        == (5,)
    )


def test_set_reaction_radical_metadata_ignores_unmapped_product_radical():
    progression = make_progression(
        make_session()
    )

    reaction = make_reaction(
        product_smiles="[CH3]",
        template_mapping={
            0: 0,
        },
        product_to_reactant_mapping={},
    )

    progression._set_reaction_radical_metadata(
        reaction,
        mol("[CH3]"),
    )

    assert reaction.is_radical is False
    assert reaction.radical_atom_idxs == ()


# =============================================================================
# Radical annotation before deduplication
# =============================================================================


def test_annotate_radicals_skips_inactive_reactions(
    monkeypatch,
):
    progression = make_progression(
        make_session()
    )

    reaction = make_reaction(
        activity_stats=False
    )

    calls = []

    monkeypatch.setattr(
        progression,
        "_sanitize_molecule",
        lambda mol: (
            calls.append(True)
            or (mol, True)
        ),
    )

    progression._annotate_radicals_before_deduplication(
        [
            reaction,
        ]
    )

    assert calls == []


def test_annotate_radicals_handles_missing_product():
    progression = make_progression(
        make_session()
    )

    reaction = make_reaction()

    reaction.product_combined_RDmol = None
    reaction.is_radical = True
    reaction.radical_atom_idxs = (
        10,
    )

    progression._annotate_radicals_before_deduplication(
        [
            reaction,
        ]
    )

    assert reaction.is_radical is False
    assert reaction.radical_atom_idxs == ()


def test_annotate_radicals_handles_failed_sanitization(
    monkeypatch,
):
    progression = make_progression(
        make_session()
    )

    reaction = make_reaction()

    reaction.is_radical = True
    reaction.radical_atom_idxs = (
        1,
    )

    monkeypatch.setattr(
        progression,
        "_sanitize_molecule",
        lambda mol: (
            mol,
            False,
        ),
    )

    progression._annotate_radicals_before_deduplication(
        [
            reaction,
        ]
    )

    assert reaction.is_radical is False
    assert reaction.radical_atom_idxs == ()


def test_annotate_radicals_delegates_metadata_on_success(
    monkeypatch,
):
    progression = make_progression(
        make_session()
    )

    reaction = make_reaction()

    sanitized = mol("CC")

    calls = []

    monkeypatch.setattr(
        progression,
        "_sanitize_molecule",
        lambda mol: (
            sanitized,
            True,
        ),
    )

    monkeypatch.setattr(
        progression,
        "_set_reaction_radical_metadata",
        lambda rxn, molecule:
        calls.append(
            (
                rxn,
                molecule,
            )
        ),
    )

    progression._annotate_radicals_before_deduplication(
        [
            reaction,
        ]
    )

    assert calls == [
        (
            reaction,
            sanitized,
        )
    ]


# =============================================================================
# Sanitization
# =============================================================================


def test_sanitize_simple_molecule_succeeds():
    progression = make_progression(
        make_session()
    )

    result, success = (
        progression
        ._sanitize_molecule(
            mol("CC")
        )
    )

    assert success is True

    assert isinstance(
        result,
        Chem.Mol,
    )


def test_sanitize_returns_new_molecule():
    progression = make_progression(
        make_session()
    )

    original = mol("CC")

    result, success = (
        progression
        ._sanitize_molecule(
            original
        )
    )

    assert success is True
    assert result is not original


def test_fix_radical_and_sanitize_returns_molecule():
    progression = make_progression(
        make_session()
    )

    result = (
        progression
        ._fix_radical_and_sanitize(
            mol("[CH3]")
        )
    )

    assert isinstance(
        result,
        Chem.Mol,
    )


# =============================================================================
# Product preparation for FG detection
# =============================================================================


def test_prepare_products_skips_inactive_reactions():
    inactive = make_reaction(
        activity_stats=False
    )

    session = make_session(
        reaction_metadata=[
            inactive
        ]
    )

    progression = make_progression(
        session
    )

    result = (
        progression
        ._prepare_products_for_idx_based_fg_detection()
    )

    assert result == []


def test_prepare_products_creates_role_for_active_reaction(
    monkeypatch,
):
    reaction = make_reaction(
        reaction_id=7,
        product_smiles="CC",
        delete_atom=False,
    )

    session = make_session(
        reaction_metadata=[
            reaction
        ]
    )

    progression = make_progression(
        session
    )

    sanitized = mol("CC")

    monkeypatch.setattr(
        progression,
        "_sanitize_molecule",
        lambda molecule: (
            sanitized,
            True,
        ),
    )

    result = (
        progression
        ._prepare_products_for_idx_based_fg_detection()
    )

    assert len(result) == 1

    role = result[0]

    assert (
        role.name
        == "new_7"
    )

    assert role.smiles == (
        Chem.MolToSmiles(
            sanitized
        )
    )

    assert (
        role.rdkit_mol
        is sanitized
    )


def test_prepare_products_uses_template_product_indices(
    monkeypatch,
):
    reaction = make_reaction(
        product_smiles="CCC",
        template_mapping={
            10: 0,
            11: 2,
        },
    )

    session = make_session(
        reaction_metadata=[
            reaction
        ]
    )

    progression = make_progression(
        session
    )

    monkeypatch.setattr(
        progression,
        "_sanitize_molecule",
        lambda molecule: (
            molecule,
            True,
        ),
    )

    role = (
        progression
        ._prepare_products_for_idx_based_fg_detection()[0]
    )

    assert (
        role.indexes_in_template
        == [
            0,
            2,
        ]
    )


def test_prepare_products_replaces_single_fragment_nondelete_product_with_sanitized_copy(
    monkeypatch,
):
    reaction = make_reaction(
        product_smiles="CC",
        delete_atom=False,
    )

    original = (
        reaction.product_combined_RDmol
    )

    sanitized = mol("CO")

    session = make_session(
        reaction_metadata=[
            reaction
        ]
    )

    progression = make_progression(
        session
    )

    monkeypatch.setattr(
        progression,
        "_sanitize_molecule",
        lambda molecule: (
            sanitized,
            True,
        ),
    )

    progression._prepare_products_for_idx_based_fg_detection()

    assert (
        reaction.product_combined_RDmol
        is not original
    )

    assert (
        Chem.MolToSmiles(
            reaction.product_combined_RDmol
        )
        == Chem.MolToSmiles(
            sanitized
        )
    )


def test_prepare_products_does_not_replace_delete_product(
    monkeypatch,
):
    reaction = make_reaction(
        product_smiles="CC",
        delete_atom=True,
    )

    original = (
        reaction.product_combined_RDmol
    )

    session = make_session(
        reaction_metadata=[
            reaction
        ]
    )

    progression = make_progression(
        session
    )

    monkeypatch.setattr(
        progression,
        "_sanitize_molecule",
        lambda molecule: (
            mol("CO"),
            True,
        ),
    )

    progression._prepare_products_for_idx_based_fg_detection()

    assert (
        reaction.product_combined_RDmol
        is original
    )


def test_prepare_products_failed_sanitization_clears_radical_metadata(
    monkeypatch,
    capsys,
):
    reaction = make_reaction()

    reaction.is_radical = True
    reaction.radical_atom_idxs = (
        1,
    )

    session = make_session(
        reaction_metadata=[
            reaction
        ]
    )

    progression = make_progression(
        session
    )

    best_effort = mol("CC")

    monkeypatch.setattr(
        progression,
        "_sanitize_molecule",
        lambda molecule: (
            best_effort,
            False,
        ),
    )

    result = (
        progression
        ._prepare_products_for_idx_based_fg_detection()
    )

    assert len(result) == 1

    assert reaction.is_radical is False
    assert reaction.radical_atom_idxs == ()

    assert (
        "sanitization failed"
        in capsys.readouterr().out
    )


# =============================================================================
# reaction_progression() control flow
# =============================================================================


def test_reaction_progression_zero_loop_returns_copy_of_existing_metadata():
    reaction = make_reaction()

    session = make_session(
        reaction_metadata=[
            reaction
        ]
    )

    progression = make_progression(
        session
    )

    result = progression.reaction_progression(
        max_loop=0
    )

    assert result == [
        reaction
    ]

    assert result is not (
        session.reaction_metadata
    )


def test_reaction_progression_negative_loop_returns_existing_metadata():
    reaction = make_reaction()

    session = make_session(
        reaction_metadata=[
            reaction
        ]
    )

    progression = make_progression(
        session
    )

    assert (
        progression.reaction_progression(
            max_loop=-1
        )
        == [
            reaction
        ]
    )


def test_reaction_progression_none_uses_session_iteration_depth(
    monkeypatch,
):
    session = make_session(
        reaction_iteration_depth=0
    )

    progression = make_progression(
        session
    )

    # If None correctly resolves to zero from the session,
    # no loop helpers should run.
    monkeypatch.setattr(
        progression,
        "_populate_monomer_roles",
        lambda: pytest.fail(
            "loop unexpectedly started"
        ),
    )

    assert (
        progression.reaction_progression(
            max_loop=None
        )
        == []
    )


def test_reaction_progression_first_iteration_populates_monomers(
    monkeypatch,
):
    session = make_session()

    progression = make_progression(
        session
    )

    calls = []

    monkeypatch.setattr(
        progression,
        "_populate_monomer_roles",
        lambda: calls.append(
            "populate"
        ),
    )

    monkeypatch.setattr(
        progression,
        "_set_is_looped_flag",
        lambda roles:
        calls.append(
            "looped"
        ),
    )

    monkeypatch.setattr(
        progression,
        "_prepare_products_for_idx_based_fg_detection",
        lambda: [],
    )

    progression.fg_detector = (
        SimpleNamespace(
            index_based_functional_groups_detector=
            lambda roles: []
        )
    )

    progression.reaction_progression(
        max_loop=1
    )

    assert calls == [
        "populate",
        "looped",
    ]

    assert (
        session
        .reaction_progression_session
        .iteration
        == 1
    )


def test_reaction_progression_stops_when_no_functional_groups(
    monkeypatch,
):
    session = make_session(
        reaction_metadata=[
            make_reaction()
        ]
    )

    progression = make_progression(
        session
    )

    monkeypatch.setattr(
        progression,
        "_populate_monomer_roles",
        lambda: None,
    )

    monkeypatch.setattr(
        progression,
        "_set_is_looped_flag",
        lambda roles: None,
    )

    prepared_roles = [
        object()
    ]

    monkeypatch.setattr(
        progression,
        "_prepare_products_for_idx_based_fg_detection",
        lambda: prepared_roles,
    )

    progression.fg_detector = (
        SimpleNamespace(
            index_based_functional_groups_detector=
            lambda roles: []
        )
    )

    result = progression.reaction_progression(
        max_loop=3
    )

    assert result == (
        session.reaction_metadata
    )

    assert (
        session
        .reaction_progression_session
        .iteration
        == 1
    )


def test_reaction_progression_stops_when_no_reactions(
    monkeypatch,
):
    initial_role = make_role()

    session = make_session(
        monomer_roles=[
            initial_role
        ]
    )

    progression = make_progression(
        session
    )

    monkeypatch.setattr(
        progression,
        "_populate_monomer_roles",
        lambda: None,
    )

    monkeypatch.setattr(
        progression,
        "_set_is_looped_flag",
        lambda roles: None,
    )

    monkeypatch.setattr(
        progression,
        "_prepare_products_for_idx_based_fg_detection",
        lambda: [
            object()
        ],
    )

    new_fg_role = make_role(
        name="new"
    )

    progression.fg_detector = (
        SimpleNamespace(
            index_based_functional_groups_detector=
            lambda roles: [
                new_fg_role
            ]
        )
    )

    progression.rxn_detector = (
        SimpleNamespace(
            index_based_reaction_detector=
            lambda roles: []
        )
    )

    result = progression.reaction_progression(
        max_loop=3
    )

    assert result == []

    assert (
        new_fg_role
        in session.monomer_roles
    )


def test_reaction_progression_converts_reaction_generator_to_list(
    monkeypatch,
):
    session = make_session()

    progression = make_progression(
        session
    )

    monkeypatch.setattr(
        progression,
        "_populate_monomer_roles",
        lambda: None,
    )

    monkeypatch.setattr(
        progression,
        "_set_is_looped_flag",
        lambda roles: None,
    )

    monkeypatch.setattr(
        progression,
        "_prepare_products_for_idx_based_fg_detection",
        lambda: [
            object()
        ],
    )

    progression.fg_detector = (
        SimpleNamespace(
            index_based_functional_groups_detector=
            lambda roles: [
                object()
            ]
        )
    )

    instance_1 = object()
    instance_2 = object()

    progression.rxn_detector = (
        SimpleNamespace(
            index_based_reaction_detector=
            lambda roles: (
                item
                for item in (
                    instance_1,
                    instance_2,
                )
            )
        )
    )

    received = []

    monkeypatch.setattr(
        progression,
        "_index_based_reaction_preparation",
        lambda reaction_instances:
        (
            received.append(
                reaction_instances
            )
            or []
        ),
    )

    monkeypatch.setattr(
        progression,
        "_annotate_radicals_before_deduplication",
        lambda reactions: None,
    )

    progression.deduplication_detector = (
        SimpleNamespace(
            compare_graphs_mol=
            lambda reactions, deep_check:
            reactions
        )
    )

    progression.reaction_progression(
        max_loop=1
    )

    assert received == [
        [
            instance_1,
            instance_2,
        ]
    ]


def test_reaction_progression_passes_deep_search_to_deduplicator(
    monkeypatch,
):
    session = make_session(
        deep_search=False
    )

    progression = make_progression(
        session
    )

    monkeypatch.setattr(
        progression,
        "_populate_monomer_roles",
        lambda: None,
    )

    monkeypatch.setattr(
        progression,
        "_set_is_looped_flag",
        lambda roles: None,
    )

    monkeypatch.setattr(
        progression,
        "_prepare_products_for_idx_based_fg_detection",
        lambda: [
            object()
        ],
    )

    progression.fg_detector = (
        SimpleNamespace(
            index_based_functional_groups_detector=
            lambda roles: [
                object()
            ]
        )
    )

    progression.rxn_detector = (
        SimpleNamespace(
            index_based_reaction_detector=
            lambda roles: [
                object()
            ]
        )
    )

    prepared = make_reaction()

    monkeypatch.setattr(
        progression,
        "_index_based_reaction_preparation",
        lambda reaction_instances: [
            prepared
        ],
    )

    monkeypatch.setattr(
        progression,
        "_annotate_radicals_before_deduplication",
        lambda reactions: None,
    )

    calls = []

    progression.deduplication_detector = (
        SimpleNamespace(
            compare_graphs_mol=
            lambda reactions, deep_check:
            (
                calls.append(
                    deep_check
                )
                or reactions
            )
        )
    )

    progression.reaction_progression(
        max_loop=1
    )

    assert calls == [
        False
    ]


def test_reaction_progression_annotates_radicals_before_deduplication(
    monkeypatch,
):
    session = make_session()

    progression = make_progression(
        session
    )

    order = []

    monkeypatch.setattr(
        progression,
        "_populate_monomer_roles",
        lambda: None,
    )

    monkeypatch.setattr(
        progression,
        "_set_is_looped_flag",
        lambda roles: None,
    )

    monkeypatch.setattr(
        progression,
        "_prepare_products_for_idx_based_fg_detection",
        lambda: [
            object()
        ],
    )

    progression.fg_detector = (
        SimpleNamespace(
            index_based_functional_groups_detector=
            lambda roles: [
                object()
            ]
        )
    )

    progression.rxn_detector = (
        SimpleNamespace(
            index_based_reaction_detector=
            lambda roles: [
                object()
            ]
        )
    )

    prepared = make_reaction()

    monkeypatch.setattr(
        progression,
        "_index_based_reaction_preparation",
        lambda reaction_instances: [
            prepared
        ],
    )

    monkeypatch.setattr(
        progression,
        "_annotate_radicals_before_deduplication",
        lambda reactions:
        order.append(
            "radicals"
        ),
    )

    progression.deduplication_detector = (
        SimpleNamespace(
            compare_graphs_mol=
            lambda reactions, deep_check:
            (
                order.append(
                    "dedup"
                )
                or reactions
            )
        )
    )

    progression.reaction_progression(
        max_loop=1
    )

    assert order == [
        "radicals",
        "dedup",
    ]


def test_reaction_progression_break_condition_stores_reactions(
    monkeypatch,
):
    initial = make_reaction(
        reaction_id=1
    )

    session = make_session(
        reaction_metadata=[
            initial
        ]
    )

    progression = make_progression(
        session
    )

    monkeypatch.setattr(
        progression,
        "_populate_monomer_roles",
        lambda: None,
    )

    monkeypatch.setattr(
        progression,
        "_set_is_looped_flag",
        lambda roles: None,
    )

    monkeypatch.setattr(
        progression,
        "_prepare_products_for_idx_based_fg_detection",
        lambda: [
            object()
        ],
    )

    progression.fg_detector = (
        SimpleNamespace(
            index_based_functional_groups_detector=
            lambda roles: [
                object()
            ]
        )
    )

    progression.rxn_detector = (
        SimpleNamespace(
            index_based_reaction_detector=
            lambda roles: [
                object()
            ]
        )
    )

    duplicate = make_reaction(
        reaction_id=2
    )

    monkeypatch.setattr(
        progression,
        "_index_based_reaction_preparation",
        lambda reaction_instances: [
            duplicate
        ],
    )

    monkeypatch.setattr(
        progression,
        "_annotate_radicals_before_deduplication",
        lambda reactions: None,
    )

    # Dedup leaves only the original reaction, so active pool size
    # stays 1 -> break.
    progression.deduplication_detector = (
        SimpleNamespace(
            compare_graphs_mol=
            lambda reactions, deep_check: [
                initial
            ]
        )
    )

    result = progression.reaction_progression(
        max_loop=5
    )

    assert result == [
        initial
    ]

    assert (
        session.reaction_metadata
        == [
            initial
        ]
    )