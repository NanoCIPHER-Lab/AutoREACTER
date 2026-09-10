from types import SimpleNamespace

import pandas as pd
import pytest
from rdkit import Chem

import AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions as prepare_module
from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import (
    MappingError,
    PrepareReactions,
    ReactionMetadata,
    SMARTSParsingError,
    ZeroActiveReactionsError,
)


# =============================================================================
# Helpers
# =============================================================================


def mol(smiles: str) -> Chem.Mol:
    molecule = Chem.MolFromSmiles(smiles)
    assert molecule is not None
    return molecule


def make_session(
    tmp_path,
    *,
    loop=False,
    deep_search=True,
    reaction_metadata=None,
    reaction_instances=None,
    reaction_id_counter=0,
):
    return SimpleNamespace(
        inputs=SimpleNamespace(
            loop=loop,
            deep_search=deep_search,
        ),
        staging_dir=tmp_path / "staging",
        reaction_metadata=list(
            reaction_metadata or []
        ),
        reaction_instances=list(
            reaction_instances or []
        ),
        reaction_id_counter=reaction_id_counter,
    )


def make_metadata(
    *,
    reaction_id=1,
    reactant_smiles="CC",
    product_smiles="CC",
    mapping=None,
    activity_stats=True,
    tmp_path=None,
):
    reactant = mol(reactant_smiles)
    product = mol(product_smiles)

    if mapping is None:
        mapping = {
            idx: idx
            for idx in range(
                min(
                    reactant.GetNumAtoms(),
                    product.GetNumAtoms(),
                )
            )
        }

    reverse_mapping = {
        product_idx: reactant_idx
        for reactant_idx, product_idx
        in mapping.items()
    }

    dataframe = pd.DataFrame(
        {
            "reactant_idx": list(
                mapping.keys()
            ),
            "product_idx": list(
                mapping.values()
            ),
        }
    )

    if len(dataframe):
        dataframe["first_shell"] = pd.Series(
            list(mapping.keys()),
            dtype="Int64",
        )

    csv_path = (
        tmp_path / f"reaction_{reaction_id}.csv"
        if tmp_path is not None
        else None
    )

    return ReactionMetadata(
        reaction_id=reaction_id,
        reactant_combined_RDmol=reactant,
        product_combined_RDmol=product,
        reactant_to_product_mapping=mapping,
        product_to_reactant_mapping=reverse_mapping,
        first_shell=list(
            mapping.keys()
        ),
        initiators=(
            list(mapping.keys())[:2]
        ),
        csv_path=csv_path,
        reaction_dataframe=dataframe,
        activity_stats=activity_stats,
    )


def make_fg(
    fg_1_indexes=None,
    fg_2_indexes=None,
):
    return SimpleNamespace(
        fg_1_indexes=fg_1_indexes,
        fg_2_indexes=fg_2_indexes,
    )


def make_monomer_role(
    *,
    smiles="C",
    rdkit_mol=None,
    is_monomer=True,
):
    return SimpleNamespace(
        smiles=smiles,
        rdkit_mol=rdkit_mol,
        is_monomer=is_monomer,
    )


def make_reaction_instance(
    *,
    monomer_1=None,
    monomer_2=None,
    same_reactants=False,
    reaction_smarts="[C:1].[O:2]>>[C:1][O:2]",
    delete_atom=False,
    functional_group_1=None,
    functional_group_2=None,
    reaction_name="test_reaction",
):
    if monomer_1 is None:
        monomer_1 = make_monomer_role(
            smiles="C",
            rdkit_mol=mol("C"),
        )

    if (
        monomer_2 is None
        and not same_reactants
    ):
        monomer_2 = make_monomer_role(
            smiles="O",
            rdkit_mol=mol("O"),
        )

    return SimpleNamespace(
        monomer_1=monomer_1,
        monomer_2=monomer_2,
        same_reactants=same_reactants,
        reaction_smarts=reaction_smarts,
        delete_atom=delete_atom,
        functional_group_1=functional_group_1,
        functional_group_2=functional_group_2,
        reaction_name=reaction_name,
    )


# =============================================================================
# Exceptions
# =============================================================================


@pytest.mark.parametrize(
    "exception_type",
    [
        MappingError,
        SMARTSParsingError,
        ZeroActiveReactionsError,
    ],
)
def test_custom_errors_are_exceptions(
    exception_type,
):
    assert issubclass(
        exception_type,
        Exception,
    )


# =============================================================================
# ReactionMetadata
# =============================================================================


def test_reaction_metadata_required_fields():
    reactant = mol("C")
    product = mol("C")

    metadata = ReactionMetadata(
        reaction_id=7,
        reactant_combined_RDmol=reactant,
        product_combined_RDmol=product,
        reactant_to_product_mapping={
            0: 0,
        },
        product_to_reactant_mapping={
            0: 0,
        },
    )

    assert metadata.reaction_id == 7
    assert (
        metadata.reactant_combined_RDmol
        is reactant
    )
    assert (
        metadata.product_combined_RDmol
        is product
    )


def test_reaction_metadata_defaults():
    metadata = ReactionMetadata(
        reaction_id=1,
        reactant_combined_RDmol=mol("C"),
        product_combined_RDmol=mol("C"),
        reactant_to_product_mapping={
            0: 0,
        },
        product_to_reactant_mapping={
            0: 0,
        },
    )

    assert (
        metadata
        .template_reactant_to_product_mapping
        is None
    )
    assert metadata.edge_atoms is None
    assert metadata.first_shell is None
    assert metadata.initiators is None
    assert metadata.byproduct_indices is None

    assert metadata.delete_atom is True
    assert metadata.delete_atom_idx is None

    assert metadata.is_radical is False
    assert metadata.radical_atom_idxs == ()
    assert metadata.activity_stats is True


def test_reaction_metadata_uses_slots():
    metadata = ReactionMetadata(
        reaction_id=1,
        reactant_combined_RDmol=mol("C"),
        product_combined_RDmol=mol("C"),
        reactant_to_product_mapping={
            0: 0,
        },
        product_to_reactant_mapping={
            0: 0,
        },
    )

    with pytest.raises(
        AttributeError
    ):
        metadata.unexpected = True


# =============================================================================
# Constructor
# =============================================================================


def test_constructor_sets_paths(
    tmp_path,
):
    session = make_session(
        tmp_path
    )

    preparer = PrepareReactions(
        session
    )

    assert preparer.session is session
    assert (
        preparer.inputs
        is session.inputs
    )

    assert preparer.staging_dir == (
        tmp_path / "staging"
    )

    assert preparer.cache == (
        tmp_path / "staging"
    )

    assert preparer.csv_cache == (
        tmp_path
        / "staging"
        / "csv_cache"
    )

    assert preparer.csv_cache.is_dir()


def test_constructor_preserves_existing_reaction_counter(
    tmp_path,
):
    session = make_session(
        tmp_path,
        reaction_id_counter=12,
    )

    PrepareReactions(
        session
    )

    assert (
        session.reaction_id_counter
        == 12
    )


def test_constructor_adds_missing_counter(
    tmp_path,
):
    session = SimpleNamespace(
        inputs=SimpleNamespace(),
        staging_dir=tmp_path,
    )

    PrepareReactions(
        session
    )

    assert (
        session.reaction_id_counter
        == 0
    )


# =============================================================================
# Zero-active-reaction guard
# =============================================================================


def test_zero_active_reactions_accepts_active_reaction(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    reaction = SimpleNamespace(
        activity_stats=True
    )

    preparer._zero_active_reactions_error(
        [
            reaction,
        ]
    )


def test_zero_active_reactions_accepts_mixed_list(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    preparer._zero_active_reactions_error(
        [
            SimpleNamespace(
                activity_stats=False
            ),
            SimpleNamespace(
                activity_stats=True
            ),
        ]
    )


def test_zero_active_reactions_rejects_all_inactive(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    with pytest.raises(
        ZeroActiveReactionsError,
        match="No active reactions",
    ):
        preparer._zero_active_reactions_error(
            [
                SimpleNamespace(
                    activity_stats=False
                )
            ]
        )


def test_zero_active_reactions_rejects_empty_list(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    with pytest.raises(
        ZeroActiveReactionsError
    ):
        preparer._zero_active_reactions_error(
            []
        )


# =============================================================================
# Functional-group index flattening
# =============================================================================


def test_flatten_fg_indexes_none(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    assert (
        preparer._flatten_fg_indexes(
            None
        )
        is None
    )


def test_flatten_fg_indexes_fg1(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    fg = make_fg(
        fg_1_indexes=[
            (
                1,
                2,
            ),
            (
                3,
                4,
            ),
        ]
    )

    assert (
        preparer._flatten_fg_indexes(
            fg
        )
        == {
            1,
            2,
            3,
            4,
        }
    )


def test_flatten_fg_indexes_both_groups(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    fg = make_fg(
        fg_1_indexes=[
            (
                1,
                2,
            )
        ],
        fg_2_indexes=[
            (
                2,
                5,
            )
        ],
    )

    assert (
        preparer._flatten_fg_indexes(
            fg
        )
        == {
            1,
            2,
            5,
        }
    )


def test_flatten_fg_indexes_empty_returns_none(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    fg = make_fg(
        fg_1_indexes=[],
        fg_2_indexes=[],
    )

    assert (
        preparer._flatten_fg_indexes(
            fg
        )
        is None
    )


# =============================================================================
# Forced initiator filtering
# =============================================================================


def test_initiators_no_restrictions(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    assert (
        preparer
        ._initiators_within_forced_indexes(
            [
                0,
                3,
            ],
            r1_atom_count=2,
            forced_indexes_1=None,
            forced_indexes_2=None,
        )
        is True
    )


def test_initiator_reactant1_allowed(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    assert (
        preparer
        ._initiators_within_forced_indexes(
            [
                1,
            ],
            2,
            {
                1,
            },
            None,
        )
    )


def test_initiator_reactant1_rejected(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    assert not (
        preparer
        ._initiators_within_forced_indexes(
            [
                1,
            ],
            2,
            {
                0,
            },
            None,
        )
    )


def test_initiator_reactant2_uses_local_index(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    # Combined index 4 with r1 size 3 -> local r2 index 1.
    assert (
        preparer
        ._initiators_within_forced_indexes(
            [
                4,
            ],
            3,
            None,
            {
                1,
            },
        )
    )


def test_initiator_reactant2_rejected(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    assert not (
        preparer
        ._initiators_within_forced_indexes(
            [
                4,
            ],
            3,
            None,
            {
                0,
            },
        )
    )


def test_initiator_at_r1_boundary_belongs_to_reactant2(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    assert (
        preparer
        ._initiators_within_forced_indexes(
            [
                3,
            ],
            3,
            None,
            {
                0,
            },
        )
    )


def test_empty_initiator_list_is_allowed(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    assert (
        preparer
        ._initiators_within_forced_indexes(
            [],
            2,
            set(),
            set(),
        )
    )


# =============================================================================
# Loop reactant copying
# =============================================================================


def test_copy_loop_reactant_missing_role_raises(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    with pytest.raises(
        SMARTSParsingError,
        match="missing its RDKit molecule",
    ):
        preparer._copy_loop_reactant_mol(
            None
        )


def test_copy_loop_reactant_missing_molecule_raises(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    role = make_monomer_role(
        rdkit_mol=None,
    )

    with pytest.raises(
        SMARTSParsingError
    ):
        preparer._copy_loop_reactant_mol(
            role
        )


def test_copy_loop_input_monomer_adds_hydrogens(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    original = mol("C")

    role = make_monomer_role(
        rdkit_mol=original,
        is_monomer=True,
    )

    result = (
        preparer
        ._copy_loop_reactant_mol(
            role
        )
    )

    assert (
        result.GetNumAtoms()
        == 5
    )

    assert (
        original.GetNumAtoms()
        == 1
    )


def test_copy_loop_generated_product_does_not_add_hydrogens(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    original = mol("C")

    role = make_monomer_role(
        rdkit_mol=original,
        is_monomer=False,
    )

    result = (
        preparer
        ._copy_loop_reactant_mol(
            role
        )
    )

    assert (
        result.GetNumAtoms()
        == 1
    )


def test_copy_loop_reactant_returns_copy(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    original = mol("CC")

    role = make_monomer_role(
        rdkit_mol=original,
        is_monomer=False,
    )

    result = (
        preparer
        ._copy_loop_reactant_mol(
            role
        )
    )

    assert result is not original


# =============================================================================
# Atom-map / isotope tracking
# =============================================================================


def test_assign_atom_maps_and_isotopes(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    r1 = mol("CC")
    r2 = mol("CO")

    preparer._assign_atom_map_numbers_and_set_isotopes(
        r1,
        r2,
    )

    assert [
        atom.GetAtomMapNum()
        for atom in r1.GetAtoms()
    ] == [
        1001,
        1002,
    ]

    assert [
        atom.GetIsotope()
        for atom in r1.GetAtoms()
    ] == [
        1001,
        1002,
    ]

    assert [
        atom.GetAtomMapNum()
        for atom in r2.GetAtoms()
    ] == [
        2001,
        2002,
    ]


def test_reassign_atom_map_numbers_from_isotope(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    molecule = mol("CC")

    molecule.GetAtomWithIdx(
        0
    ).SetIsotope(1234)

    molecule.GetAtomWithIdx(
        1
    ).SetIsotope(5678)

    preparer._reassign_atom_map_numbers_by_isotope(
        molecule
    )

    assert [
        atom.GetAtomMapNum()
        for atom in molecule.GetAtoms()
    ] == [
        1234,
        5678,
    ]

    assert [
        atom.GetIsotope()
        for atom in molecule.GetAtoms()
    ] == [
        0,
        0,
    ]


def test_reassign_zero_isotope_does_not_overwrite_map(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    molecule = mol("C")

    atom = molecule.GetAtomWithIdx(
        0
    )

    atom.SetAtomMapNum(77)
    atom.SetIsotope(0)

    preparer._reassign_atom_map_numbers_by_isotope(
        molecule
    )

    assert atom.GetAtomMapNum() == 77


def test_reveal_template_map_numbers(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    molecule = mol("CC")

    atom = molecule.GetAtomWithIdx(
        0
    )

    atom.SetIntProp(
        "old_mapno",
        2,
    )

    preparer._reveal_template_map_numbers(
        molecule
    )

    assert atom.GetAtomMapNum() == 2


def test_reveal_template_map_numbers_ignores_atom_without_property(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    molecule = mol("C")

    atom = molecule.GetAtomWithIdx(
        0
    )

    atom.SetAtomMapNum(88)

    preparer._reveal_template_map_numbers(
        molecule
    )

    assert atom.GetAtomMapNum() == 88


def test_clear_isotopes(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    first = mol("CC")
    second = mol("CO")

    for atom in first.GetAtoms():
        atom.SetIsotope(13)

    for atom in second.GetAtoms():
        atom.SetIsotope(18)

    preparer._clear_isotopes(
        first,
        second,
    )

    assert all(
        atom.GetIsotope() == 0
        for atom in first.GetAtoms()
    )

    assert all(
        atom.GetIsotope() == 0
        for atom in second.GetAtoms()
    )


# =============================================================================
# Reaction / reactant construction
# =============================================================================


def test_build_reaction_returns_rdkit_reaction(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    reaction = preparer._build_reaction(
        "[C:1].[O:2]>>[C:1][O:2]"
    )

    assert reaction is not None

    assert reaction.GetNumReactantTemplates() == 2
    assert reaction.GetNumProductTemplates() == 1


def test_build_reactants_valid_smiles_adds_explicit_hydrogens(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    first, second = (
        preparer._build_reactants(
            "C",
            "O",
        )
    )

    # CH4
    assert first.GetNumAtoms() == 5

    # H2O
    assert second.GetNumAtoms() == 3


def test_build_reactants_invalid_first_smiles(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    with pytest.raises(
        SMARTSParsingError,
        match="first reactant",
    ):
        preparer._build_reactants(
            "not-a-smiles",
            "O",
        )


def test_build_reactants_invalid_second_smiles(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    with pytest.raises(
        SMARTSParsingError,
        match="second reactant",
    ):
        preparer._build_reactants(
            "C",
            "not-a-smiles",
        )


def test_build_reaction_tuple_same_reactant(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    first = mol("C")
    second = mol("O")

    result = preparer._build_reaction_tuple(
        True,
        first,
        second,
    )

    assert len(result) == 1

    assert result[0][0] is first
    assert result[0][1] is first


def test_build_reaction_tuple_different_reactants(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    first = mol("C")
    second = mol("O")

    result = preparer._build_reaction_tuple(
        False,
        first,
        second,
    )

    assert result == [
        [
            first,
            second,
        ],
        [
            second,
            first,
        ],
    ]


# =============================================================================
# Consecutive-number helper
# =============================================================================


@pytest.mark.parametrize(
    "values",
    [
        [1],
        [1, 2],
        [3, 2, 1],
        [-2, -1, 0, 1],
    ],
)
def test_is_consecutive_true(
    tmp_path,
    values,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    assert preparer._is_consecutive(
        values
    )


@pytest.mark.parametrize(
    "values",
    [
        [],
        [1, 1],
        [1, 3],
        [0, 2, 3],
    ],
)
def test_is_consecutive_false(
    tmp_path,
    values,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    assert not preparer._is_consecutive(
        values
    )


# =============================================================================
# Atom-index mapping construction
# =============================================================================


def test_build_atom_index_mapping(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    reactant = mol("CC")
    product = mol("CC")

    reactant.GetAtomWithIdx(
        0
    ).SetAtomMapNum(1001)

    reactant.GetAtomWithIdx(
        1
    ).SetAtomMapNum(1002)

    product.GetAtomWithIdx(
        0
    ).SetAtomMapNum(1002)

    product.GetAtomWithIdx(
        1
    ).SetAtomMapNum(1001)

    mapping, dataframe = (
        preparer._build_atom_index_mapping(
            reactant,
            product,
        )
    )

    assert mapping == {
        0: 1,
        1: 0,
    }

    assert dataframe.to_dict(
        orient="records"
    ) == [
        {
            "reactant_idx": 0,
            "product_idx": 1,
        },
        {
            "reactant_idx": 1,
            "product_idx": 0,
        },
    ]


def test_build_atom_index_mapping_ignores_zero_map_number(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    reactant = mol("CC")
    product = mol("CC")

    reactant.GetAtomWithIdx(
        0
    ).SetAtomMapNum(0)

    reactant.GetAtomWithIdx(
        1
    ).SetAtomMapNum(10)

    product.GetAtomWithIdx(
        0
    ).SetAtomMapNum(10)

    mapping, _ = (
        preparer._build_atom_index_mapping(
            reactant,
            product,
        )
    )

    assert mapping == {
        1: 0,
    }


def test_build_atom_index_mapping_ignores_unmatched_reactant_tag(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    reactant = mol("C")
    product = mol("C")

    reactant.GetAtomWithIdx(
        0
    ).SetAtomMapNum(123)

    product.GetAtomWithIdx(
        0
    ).SetAtomMapNum(456)

    mapping, dataframe = (
        preparer._build_atom_index_mapping(
            reactant,
            product,
        )
    )

    assert mapping == {}
    assert dataframe.empty


# =============================================================================
# Mapping validation
# =============================================================================


def test_validate_mapping_accepts_complete_bijection(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    dataframe = pd.DataFrame(
        {
            "reactant_idx": [
                0,
                1,
            ],
            "product_idx": [
                1,
                0,
            ],
        }
    )

    preparer._validate_mapping(
        dataframe,
        mol("CC"),
        mol("CC"),
    )


def test_validate_mapping_rejects_none_dataframe(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    with pytest.raises(
        MappingError,
        match="empty dataframe",
    ):
        preparer._validate_mapping(
            None,
            mol("C"),
            mol("C"),
        )


def test_validate_mapping_rejects_empty_dataframe(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    with pytest.raises(
        MappingError,
        match="empty dataframe",
    ):
        preparer._validate_mapping(
            pd.DataFrame(),
            mol("C"),
            mol("C"),
        )


def test_validate_mapping_requires_columns(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    with pytest.raises(
        MappingError,
        match="required columns",
    ):
        preparer._validate_mapping(
            pd.DataFrame(
                {
                    "wrong": [
                        0
                    ]
                }
            ),
            mol("C"),
            mol("C"),
        )


def test_validate_mapping_rejects_count_mismatch(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    dataframe = pd.DataFrame(
        {
            "reactant_idx": pd.Series(
                [
                    0,
                    1,
                ],
                dtype="Int64",
            ),
            "product_idx": pd.Series(
                [
                    0,
                    pd.NA,
                ],
                dtype="Int64",
            ),
        }
    )

    with pytest.raises(
        MappingError,
        match="mismatch in atom counts",
    ):
        preparer._validate_mapping(
            dataframe,
            mol("CC"),
            mol("CC"),
        )


def test_validate_mapping_rejects_duplicate_reactant_idx(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    dataframe = pd.DataFrame(
        {
            "reactant_idx": [
                0,
                0,
            ],
            "product_idx": [
                0,
                1,
            ],
        }
    )

    with pytest.raises(
        MappingError,
        match="duplicate idxs.*reactant",
    ):
        preparer._validate_mapping(
            dataframe,
            mol("CC"),
            mol("CC"),
        )


def test_validate_mapping_rejects_duplicate_product_idx(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    dataframe = pd.DataFrame(
        {
            "reactant_idx": [
                0,
                1,
            ],
            "product_idx": [
                0,
                0,
            ],
        }
    )

    with pytest.raises(
        MappingError,
        match="duplicate idxs.*product",
    ):
        preparer._validate_mapping(
            dataframe,
            mol("CC"),
            mol("CC"),
        )


def test_validate_mapping_rejects_reactant_idx_above_bounds(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    dataframe = pd.DataFrame(
        {
            "reactant_idx": [
                0,
                2,
            ],
            "product_idx": [
                0,
                1,
            ],
        }
    )

    with pytest.raises(
        MappingError,
        match="reactant idx out of bounds",
    ):
        preparer._validate_mapping(
            dataframe,
            mol("CC"),
            mol("CC"),
        )


def test_validate_mapping_rejects_product_idx_above_bounds(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    dataframe = pd.DataFrame(
        {
            "reactant_idx": [
                0,
                1,
            ],
            "product_idx": [
                0,
                2,
            ],
        }
    )

    with pytest.raises(
        MappingError,
        match="product idx out of bounds",
    ):
        preparer._validate_mapping(
            dataframe,
            mol("CC"),
            mol("CC"),
        )


def test_validate_mapping_rejects_negative_reactant_idx(
    tmp_path,
):
    """
    A valid RDKit atom mapping cannot contain negative dataframe indices.

    This is an invariant test, not a change to reaction chemistry.
    """
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    dataframe = pd.DataFrame(
        {
            "reactant_idx": [
                -1,
                1,
            ],
            "product_idx": [
                0,
                1,
            ],
        }
    )

    with pytest.raises(
        MappingError,
        match="reactant idx out of bounds",
    ):
        preparer._validate_mapping(
            dataframe,
            mol("CC"),
            mol("CC"),
        )


def test_validate_mapping_rejects_negative_product_idx(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    dataframe = pd.DataFrame(
        {
            "reactant_idx": [
                0,
                1,
            ],
            "product_idx": [
                -1,
                1,
            ],
        }
    )

    with pytest.raises(
        MappingError,
        match="product idx out of bounds",
    ):
        preparer._validate_mapping(
            dataframe,
            mol("CC"),
            mol("CC"),
        )


def test_validate_mapping_rejects_incomplete_reactant(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    dataframe = pd.DataFrame(
        {
            "reactant_idx": [
                0,
                1,
            ],
            "product_idx": [
                0,
                1,
            ],
        }
    )

    with pytest.raises(
        MappingError,
        match="incomplete mapping for reactant",
    ):
        preparer._validate_mapping(
            dataframe,
            mol("CCC"),
            mol("CC"),
        )


def test_validate_mapping_rejects_incomplete_product(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    dataframe = pd.DataFrame(
        {
            "reactant_idx": [
                0,
                1,
            ],
            "product_idx": [
                0,
                1,
            ],
        }
    )

    with pytest.raises(
        MappingError,
        match="incomplete mapping for product",
    ):
        preparer._validate_mapping(
            dataframe,
            mol("CC"),
            mol("CCC"),
        )


# =============================================================================
# First shell / initiators
# =============================================================================


def test_assign_first_shell_and_initiators(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    reactant = mol("CCC")
    product = mol("CCC")

    product.GetAtomWithIdx(
        0
    ).SetAtomMapNum(1)

    product.GetAtomWithIdx(
        1
    ).SetAtomMapNum(2)

    product.GetAtomWithIdx(
        2
    ).SetAtomMapNum(999)

    first_shell, initiators = (
        preparer
        ._assign_first_shell_and_initiators(
            reactant,
            product,
            {
                0: 0,
                1: 1,
                2: 2,
            },
        )
    )

    assert first_shell == [
        0,
        1,
    ]

    assert initiators == [
        0,
        1,
    ]

    assert (
        reactant
        .GetAtomWithIdx(0)
        .GetAtomMapNum()
        == 1
    )

    assert (
        reactant
        .GetAtomWithIdx(1)
        .GetAtomMapNum()
        == 2
    )


def test_assign_first_shell_missing_reverse_mapping_raises(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    reactant = mol("CC")
    product = mol("CC")

    product.GetAtomWithIdx(
        0
    ).SetAtomMapNum(1)

    product.GetAtomWithIdx(
        1
    ).SetAtomMapNum(2)

    with pytest.raises(
        ValueError,
        match="not found in mapping",
    ):
        (
            preparer
            ._assign_first_shell_and_initiators(
                reactant,
                product,
                {
                    0: 0,
                },
            )
        )


def test_assign_first_shell_requires_exactly_two_initiators(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    reactant = mol("CC")
    product = mol("CC")

    product.GetAtomWithIdx(
        0
    ).SetAtomMapNum(1)

    product.GetAtomWithIdx(
        1
    ).SetAtomMapNum(3)

    with pytest.raises(
        ValueError,
        match="Expected 2 initiators",
    ):
        (
            preparer
            ._assign_first_shell_and_initiators(
                reactant,
                product,
                {
                    0: 0,
                    1: 1,
                },
            )
        )


# =============================================================================
# Byproduct detection
# =============================================================================


def test_detect_byproducts_disabled_returns_empty(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    assert (
        preparer._detect_byproducts(
            mol("CC.O"),
            {
                0: 0,
                1: 1,
                2: 2,
            },
            False,
        )
        == []
    )


def test_detect_byproducts_uses_smallest_fragment(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    product = mol(
        "CC.O"
    )

    result = preparer._detect_byproducts(
        product,
        {
            0: 10,
            1: 11,
            2: 12,
        },
        True,
    )

    assert result == [
        12,
    ]


def test_detect_byproducts_skips_unmapped_product_atoms(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    product = mol(
        "CC.O"
    )

    result = preparer._detect_byproducts(
        product,
        {
            0: 10,
            1: 11,
        },
        True,
    )

    assert result == []


# =============================================================================
# Duplicate handling
# =============================================================================


def test_detect_duplicates_keeps_first_unique(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    first = make_metadata(
        reaction_id=1
    )

    result = preparer._detect_duplicates(
        [
            first,
        ]
    )

    assert result == [
        first,
    ]

    assert first.activity_stats is True


def test_detect_duplicates_marks_later_duplicate_inactive_and_excludes_it(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    first = make_metadata(
        reaction_id=1,
        reactant_smiles="CC",
        product_smiles="CO",
    )

    second = make_metadata(
        reaction_id=2,
        reactant_smiles="CC",
        product_smiles="CO",
    )

    result = preparer._detect_duplicates(
        [
            first,
            second,
        ]
    )

    # Characterizes current runtime behavior:
    # duplicate is disabled AND omitted from the returned unique list.
    assert result == [
        first,
    ]

    assert first.activity_stats is True
    assert second.activity_stats is False


def test_detect_duplicates_same_reactant_different_product_is_unique(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    first = make_metadata(
        reaction_id=1,
        reactant_smiles="CC",
        product_smiles="CO",
    )

    second = make_metadata(
        reaction_id=2,
        reactant_smiles="CC",
        product_smiles="CN",
    )

    result = preparer._detect_duplicates(
        [
            first,
            second,
        ]
    )

    assert result == [
        first,
        second,
    ]


# =============================================================================
# Reaction instance processing - initial stage
# =============================================================================


def test_process_reaction_instances_initial_stage_uses_smiles(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    preparer = PrepareReactions(
        session
    )

    reaction = make_reaction_instance(
        monomer_1=make_monomer_role(
            smiles="C"
        ),
        monomer_2=make_monomer_role(
            smiles="O"
        ),
    )

    first = mol("C")
    second = mol("O")

    calls = []

    monkeypatch.setattr(
        preparer,
        "_build_reactants",
        lambda s1, s2: (
            calls.append(
                (
                    "reactants",
                    s1,
                    s2,
                )
            )
            or (
                first,
                second,
            )
        ),
    )

    reaction_tuple = [
        [
            first,
            second,
        ]
    ]

    monkeypatch.setattr(
        preparer,
        "_build_reaction_tuple",
        lambda same, r1, r2: (
            calls.append(
                (
                    "tuple",
                    same,
                )
            )
            or reaction_tuple
        ),
    )

    fake_rxn = object()

    monkeypatch.setattr(
        preparer,
        "_build_reaction",
        lambda smarts: (
            calls.append(
                (
                    "reaction",
                    smarts,
                )
            )
            or fake_rxn
        ),
    )

    monkeypatch.setattr(
        preparer,
        "_process_reaction_products",
        lambda **kwargs: (
            calls.append(
                (
                    "products",
                    kwargs,
                )
            )
            or [
                "metadata"
            ]
        ),
    )

    result = (
        preparer
        ._process_reaction_instances(
            [
                reaction,
            ],
            loop=False,
        )
    )

    assert result == [
        "metadata"
    ]

    assert (
        "reactants",
        "C",
        "O",
    ) in calls


def test_process_reaction_instances_same_reactants_uses_first_smiles_twice(
    tmp_path,
    monkeypatch,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    reaction = make_reaction_instance(
        monomer_1=make_monomer_role(
            smiles="CC"
        ),
        monomer_2=None,
        same_reactants=True,
    )

    captured = []

    monkeypatch.setattr(
        preparer,
        "_build_reactants",
        lambda first, second: (
            captured.append(
                (
                    first,
                    second,
                )
            )
            or (
                mol("CC"),
                mol("CC"),
            )
        ),
    )

    monkeypatch.setattr(
        preparer,
        "_build_reaction_tuple",
        lambda *args: [],
    )

    monkeypatch.setattr(
        preparer,
        "_build_reaction",
        lambda smarts: object(),
    )

    monkeypatch.setattr(
        preparer,
        "_process_reaction_products",
        lambda **kwargs:
        kwargs["reaction_metadata"],
    )

    preparer._process_reaction_instances(
        [
            reaction,
        ]
    )

    assert captured == [
        (
            "CC",
            "CC",
        )
    ]


# =============================================================================
# Reaction instance processing - loop mode
# =============================================================================


def test_process_reaction_instances_loop_skips_missing_monomer1_molecule(
    tmp_path,
    monkeypatch,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    reaction = make_reaction_instance(
        monomer_1=make_monomer_role(
            rdkit_mol=None,
        ),
        monomer_2=make_monomer_role(
            rdkit_mol=mol("O"),
        ),
    )

    monkeypatch.setattr(
        preparer,
        "_process_reaction_products",
        lambda **kwargs:
        pytest.fail(
            "reaction should have been skipped"
        ),
    )

    result = (
        preparer
        ._process_reaction_instances(
            [
                reaction,
            ],
            loop=True,
        )
    )

    assert result == []


def test_process_reaction_instances_loop_skips_missing_monomer2_molecule(
    tmp_path,
    monkeypatch,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    reaction = make_reaction_instance(
        monomer_1=make_monomer_role(
            rdkit_mol=mol("C"),
        ),
        monomer_2=make_monomer_role(
            rdkit_mol=None,
        ),
    )

    monkeypatch.setattr(
        preparer,
        "_process_reaction_products",
        lambda **kwargs:
        pytest.fail(
            "reaction should have been skipped"
        ),
    )

    result = (
        preparer
        ._process_reaction_instances(
            [
                reaction,
            ],
            loop=True,
        )
    )

    assert result == []


def test_process_reaction_instances_loop_passes_forced_fg_indices(
    tmp_path,
    monkeypatch,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    reaction = make_reaction_instance(
        monomer_1=make_monomer_role(
            rdkit_mol=mol("CC"),
            is_monomer=False,
        ),
        monomer_2=make_monomer_role(
            rdkit_mol=mol("CO"),
            is_monomer=False,
        ),
        functional_group_1=make_fg(
            fg_1_indexes=[
                (
                    1,
                    2,
                )
            ]
        ),
        functional_group_2=make_fg(
            fg_1_indexes=[
                (
                    3,
                    4,
                )
            ]
        ),
    )

    captured = []

    monkeypatch.setattr(
        preparer,
        "_build_reaction",
        lambda smarts: object(),
    )

    monkeypatch.setattr(
        preparer,
        "_process_reaction_products",
        lambda **kwargs: (
            captured.append(
                kwargs
            )
            or []
        ),
    )

    preparer._process_reaction_instances(
        [
            reaction,
        ],
        loop=True,
    )

    assert len(captured) == 1

    assert (
        captured[0][
            "forced_indexes_1"
        ]
        == {
            1,
            2,
        }
    )

    assert (
        captured[0][
            "forced_indexes_2"
        ]
        == {
            3,
            4,
        }
    )

    assert len(
        captured[0][
            "reaction_tuple"
        ]
    ) == 1


# =============================================================================
# Reaction product processing orchestration
# =============================================================================


class FakeReaction:
    def __init__(
        self,
        responses,
    ):
        self.responses = list(
            responses
        )
        self.calls = []

    def RunReactants(
        self,
        reactants,
    ):
        self.calls.append(
            reactants
        )

        if self.responses:
            return self.responses.pop(0)

        return ()


def configure_product_processing_helpers(
    preparer,
    monkeypatch,
    *,
    first_shell=None,
    initiators=None,
    byproducts=None,
):
    if first_shell is None:
        first_shell = [
            0,
            1,
        ]

    if initiators is None:
        initiators = [
            0,
            1,
        ]

    if byproducts is None:
        byproducts = []

    monkeypatch.setattr(
        preparer,
        "_reassign_atom_map_numbers_by_isotope",
        lambda molecule: None,
    )

    monkeypatch.setattr(
        preparer,
        "_build_atom_index_mapping",
        lambda reactant, product: (
            {
                0: 0,
                1: 1,
            },
            pd.DataFrame(
                {
                    "reactant_idx": [
                        0,
                        1,
                    ],
                    "product_idx": [
                        0,
                        1,
                    ],
                }
            ),
        ),
    )

    monkeypatch.setattr(
        preparer,
        "_reveal_template_map_numbers",
        lambda molecule: None,
    )

    monkeypatch.setattr(
        preparer,
        "_validate_mapping",
        lambda *args: None,
    )

    monkeypatch.setattr(
        preparer,
        "_assign_first_shell_and_initiators",
        lambda *args: (
            first_shell,
            initiators,
        ),
    )

    monkeypatch.setattr(
        preparer,
        "_detect_byproducts",
        lambda *args:
        byproducts,
    )


def test_process_reaction_products_creates_metadata_and_csv(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    preparer = PrepareReactions(
        session
    )

    product = Chem.CombineMols(
        mol("C"),
        mol("O"),
    )

    reaction = FakeReaction(
        [
            (
                (
                    product,
                ),
            )
        ]
    )

    configure_product_processing_helpers(
        preparer,
        monkeypatch,
    )

    result = (
        preparer
        ._process_reaction_products(
            rxn=reaction,
            csv_cache=preparer.csv_cache,
            reaction_tuple=[
                [
                    mol("C"),
                    mol("O"),
                ]
            ],
            delete_atoms=False,
        )
    )

    assert len(result) == 1

    metadata = result[0]

    assert metadata.reaction_id == 1

    assert (
        session.reaction_id_counter
        == 1
    )

    assert metadata.csv_path.is_file()

    assert metadata.first_shell == [
        0,
        1,
    ]

    assert metadata.initiators == [
        0,
        1,
    ]

    assert (
        metadata.activity_stats
        is True
    )


def test_process_reaction_products_appends_to_existing_list(
    tmp_path,
    monkeypatch,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    existing = [
        "existing"
    ]

    product = Chem.CombineMols(
        mol("C"),
        mol("O"),
    )

    reaction = FakeReaction(
        [
            (
                (
                    product,
                ),
            )
        ]
    )

    configure_product_processing_helpers(
        preparer,
        monkeypatch,
    )

    result = (
        preparer
        ._process_reaction_products(
            reaction,
            preparer.csv_cache,
            [
                [
                    mol("C"),
                    mol("O"),
                ]
            ],
            reaction_metadata=existing,
        )
    )

    assert result is existing
    assert result[0] == "existing"
    assert len(result) == 2


def test_process_reaction_products_tries_reverse_order_after_failure(
    tmp_path,
    monkeypatch,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    product = Chem.CombineMols(
        mol("C"),
        mol("O"),
    )

    reaction = FakeReaction(
        [
            (),
            (
                (
                    product,
                ),
            ),
        ]
    )

    configure_product_processing_helpers(
        preparer,
        monkeypatch,
    )

    result = (
        preparer
        ._process_reaction_products(
            reaction,
            preparer.csv_cache,
            [
                [
                    mol("C"),
                    mol("O"),
                ]
            ],
        )
    )

    assert len(reaction.calls) == 2
    assert len(result) == 1


def test_process_reaction_products_both_orders_fail_returns_no_metadata(
    tmp_path,
    capsys,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    # Use a real RDKit ChemicalReaction here because the failure-reporting
    # path calls ReactionToSmarts(rxn).
    #
    # This reaction requires two nitrogen reactants, while we provide
    # carbon and oxygen, so both reactant orders legitimately produce
    # no products.
    reaction = preparer._build_reaction(
        "[N:1].[N:2]>>[N:1][N:2]"
    )

    result = (
        preparer
        ._process_reaction_products(
            reaction,
            preparer.csv_cache,
            [
                [
                    mol("C"),
                    mol("O"),
                ]
            ],
        )
    )

    assert result == []

    assert (
        preparer.session
        .reaction_id_counter
        == 0
    )

    output = capsys.readouterr().out

    assert (
        "Reaction failed in both orders"
        in output
    )

    assert (
        "RDKit failed to react"
        in output
    )

    assert (
        "Reaction SMARTS:"
        in output
    )


def test_process_reaction_products_forced_indexes_can_reject_product(
    tmp_path,
    monkeypatch,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    product = Chem.CombineMols(
        mol("C"),
        mol("O"),
    )

    reaction = FakeReaction(
        [
            (
                (
                    product,
                ),
            )
        ]
    )

    configure_product_processing_helpers(
        preparer,
        monkeypatch,
        initiators=[
            0,
            1,
        ],
    )

    result = (
        preparer
        ._process_reaction_products(
            reaction,
            preparer.csv_cache,
            [
                [
                    mol("C"),
                    mol("O"),
                ]
            ],
            forced_indexes_1={
                999,
            },
            forced_indexes_2={
                999,
            },
        )
    )

    assert result == []

    assert (
        preparer.session
        .reaction_id_counter
        == 0
    )


def test_process_reaction_products_sets_delete_atom_idx(
    tmp_path,
    monkeypatch,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    product = Chem.CombineMols(
        mol("C"),
        mol("O"),
    )

    reaction = FakeReaction(
        [
            (
                (
                    product,
                ),
            )
        ]
    )

    configure_product_processing_helpers(
        preparer,
        monkeypatch,
        byproducts=[
            1,
        ],
    )

    result = (
        preparer
        ._process_reaction_products(
            reaction,
            preparer.csv_cache,
            [
                [
                    mol("C"),
                    mol("O"),
                ]
            ],
            delete_atoms=True,
        )
    )

    metadata = result[0]

    assert metadata.byproduct_indices == [
        1,
    ]

    assert metadata.delete_atom_idx == 1


# =============================================================================
# _prepare_reactions_stage
# =============================================================================


def test_prepare_reactions_stage_adds_template_mapping_and_edges(
    tmp_path,
    monkeypatch,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    reaction = make_metadata(
        reaction_id=1,
        reactant_smiles="CCC",
        product_smiles="CCC",
        mapping={
            0: 0,
            1: 1,
            2: 2,
        },
        tmp_path=tmp_path,
    )

    reaction.reaction_dataframe = (
        pd.DataFrame(
            {
                "reactant_idx": [
                    0,
                    1,
                    2,
                ],
                "product_idx": [
                    0,
                    1,
                    2,
                ],
                "first_shell": pd.Series(
                    [
                        1,
                        pd.NA,
                        pd.NA,
                    ],
                    dtype="Int64",
                ),
            }
        )
    )

    monkeypatch.setattr(
        preparer,
        "_process_reaction_instances",
        lambda detected, loop=False: [
            reaction
        ],
    )

    monkeypatch.setattr(
        preparer,
        "_detect_duplicates",
        lambda reactions:
        reactions,
    )

    monkeypatch.setattr(
        prepare_module,
        "reaction_atom_walker",
        lambda molecule, first_shell, mapping:
        (
            {
                0: 0,
                1: 1,
            },
            [
                2,
            ],
        ),
    )

    session = make_session(
        tmp_path,
        reaction_instances=[
            object()
        ],
    )

    result = (
        preparer
        ._prepare_reactions_stage(
            session
        )
    )

    assert result == [
        reaction
    ]

    assert (
        reaction
        .template_reactant_to_product_mapping
        == {
            0: 0,
            1: 1,
        }
    )

    assert reaction.edge_atoms == [
        2,
    ]

    assert (
        "template_reactant_idx"
        in reaction
        .reaction_dataframe
        .columns
    )

    assert (
        "edge_atoms"
        in reaction
        .reaction_dataframe
        .columns
    )

    assert reaction.csv_path.is_file()


def test_prepare_reactions_stage_skips_inactive_reactions(
    tmp_path,
    monkeypatch,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    reaction = make_metadata(
        activity_stats=False,
        tmp_path=tmp_path,
    )

    monkeypatch.setattr(
        preparer,
        "_process_reaction_instances",
        lambda detected, loop=False: [
            reaction
        ],
    )

    monkeypatch.setattr(
        preparer,
        "_detect_duplicates",
        lambda reactions:
        reactions,
    )

    monkeypatch.setattr(
        prepare_module,
        "reaction_atom_walker",
        lambda *args:
        pytest.fail(
            "inactive reaction should not be walked"
        ),
    )

    result = (
        preparer
        ._prepare_reactions_stage(
            []
        )
    )

    assert result == [
        reaction
    ]


def test_prepare_reactions_stage_accepts_direct_reaction_instance_list(
    tmp_path,
    monkeypatch,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    instances = [
        object(),
        object(),
    ]

    captured = []

    monkeypatch.setattr(
        preparer,
        "_process_reaction_instances",
        lambda detected, loop=False: (
            captured.append(
                (
                    detected,
                    loop,
                )
            )
            or []
        ),
    )

    monkeypatch.setattr(
        preparer,
        "_detect_duplicates",
        lambda reactions:
        reactions,
    )

    result = (
        preparer
        ._prepare_reactions_stage(
            instances,
            loop=True,
        )
    )

    assert result == []

    assert captured == [
        (
            instances,
            True,
        )
    ]


# =============================================================================
# _index_based_reaction_preparation
# =============================================================================


def test_index_based_reaction_preparation_creates_fresh_preparer_and_uses_loop_mode(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    preparer = PrepareReactions(
        session
    )

    calls = []

    class FakePreparer:
        def __init__(
            self,
            received_session,
        ):
            calls.append(
                (
                    "init",
                    received_session,
                )
            )

        def _prepare_reactions_stage(
            self,
            reaction_instances,
            loop=False,
        ):
            calls.append(
                (
                    "stage",
                    reaction_instances,
                    loop,
                )
            )

            return [
                "prepared"
            ]

    monkeypatch.setattr(
        prepare_module,
        "PrepareReactions",
        FakePreparer,
    )

    instances = [
        object()
    ]

    result = (
        preparer
        ._index_based_reaction_preparation(
            instances
        )
    )

    assert result == [
        "prepared"
    ]

    assert calls == [
        (
            "init",
            session,
        ),
        (
            "stage",
            instances,
            True,
        ),
    ]


# =============================================================================
# prepare_reactions() orchestration
# =============================================================================


def test_prepare_reactions_nonloop_deduplicates(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path,
        loop=False,
        deep_search=False,
    )

    preparer = PrepareReactions(
        session
    )

    prepared = [
        SimpleNamespace(
            activity_stats=True
        )
    ]

    final = [
        SimpleNamespace(
            activity_stats=True
        )
    ]

    monkeypatch.setattr(
        preparer,
        "_prepare_reactions_stage",
        lambda received_session:
        prepared,
    )

    calls = []

    class FakeDeduplicator:
        def compare_graphs_mol(
            self,
            reactions,
            deep_check=True,
        ):
            calls.append(
                (
                    reactions,
                    deep_check,
                )
            )

            return final

    monkeypatch.setattr(
        prepare_module,
        "DeduplicationDetector",
        FakeDeduplicator,
    )

    result = preparer.prepare_reactions(
        session
    )

    assert result is session

    assert (
        session.reaction_metadata
        is final
    )

    assert calls == [
        (
            prepared,
            False,
        )
    ]


def test_prepare_reactions_nonloop_rejects_zero_active_before_dedup(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path,
        loop=False,
    )

    preparer = PrepareReactions(
        session
    )

    monkeypatch.setattr(
        preparer,
        "_prepare_reactions_stage",
        lambda received_session: [
            SimpleNamespace(
                activity_stats=False
            )
        ],
    )

    with pytest.raises(
        ZeroActiveReactionsError
    ):
        preparer.prepare_reactions(
            session
        )


def test_prepare_reactions_loop_uses_progression(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path,
        loop=True,
    )

    preparer = PrepareReactions(
        session
    )

    initial = [
        SimpleNamespace(
            activity_stats=True
        )
    ]

    final = [
        SimpleNamespace(
            activity_stats=True
        )
    ]

    monkeypatch.setattr(
        preparer,
        "_prepare_reactions_stage",
        lambda received_session:
        initial,
    )

    calls = []

    class FakeProgression:
        def __init__(
            self,
            received_session,
            preparer=None,
        ):
            calls.append(
                (
                    "init",
                    received_session,
                    preparer,
                )
            )

        def reaction_progression(
            self,
        ):
            calls.append(
                (
                    "run",
                )
            )

            return final

    monkeypatch.setattr(
        prepare_module,
        "ReactionProgression",
        FakeProgression,
    )

    result = preparer.prepare_reactions(
        session
    )

    assert result is session

    assert (
        session.reaction_metadata
        is final
    )

    assert calls == [
        (
            "init",
            session,
            preparer,
        ),
        (
            "run",
        ),
    ]


def test_prepare_reactions_loop_rejects_zero_active_final_pool(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path,
        loop=True,
    )

    preparer = PrepareReactions(
        session
    )

    monkeypatch.setattr(
        preparer,
        "_prepare_reactions_stage",
        lambda received_session: [
            SimpleNamespace(
                activity_stats=True
            )
        ],
    )

    class FakeProgression:
        def __init__(
            self,
            session,
            preparer=None,
        ):
            pass

        def reaction_progression(
            self,
        ):
            return [
                SimpleNamespace(
                    activity_stats=False
                )
            ]

    monkeypatch.setattr(
        prepare_module,
        "ReactionProgression",
        FakeProgression,
    )

    with pytest.raises(
        ZeroActiveReactionsError
    ):
        preparer.prepare_reactions(
            session
        )


# =============================================================================
# Highlighted reaction-template images
# =============================================================================


def make_image_metadata():
    metadata = make_metadata(
        reaction_id=1,
        reactant_smiles="CCO",
        product_smiles="CCO",
        mapping={
            0: 0,
            1: 1,
            2: 2,
        },
    )

    metadata.template_reactant_to_product_mapping = {
        0: 0,
        1: 1,
    }

    metadata.edge_atoms = [
        2,
    ]

    metadata.byproduct_indices = [
        2,
    ]

    metadata.delete_atom = True

    metadata.reaction_dataframe = (
        pd.DataFrame(
            {
                "initiators": pd.Series(
                    [
                        0,
                        1,
                        pd.NA,
                    ],
                    dtype="Int64",
                )
            }
        )
    )

    return metadata


def test_reaction_templates_image_empty_metadata_returns_none(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    session = SimpleNamespace(
        reaction_metadata=[]
    )

    assert (
        preparer
        .reaction_templates_highlighted_image_grid(
            session
        )
        is None
    )


def test_reaction_templates_image_all_inactive_returns_none(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    metadata = make_image_metadata()
    metadata.activity_stats = False

    session = SimpleNamespace(
        reaction_metadata=[
            metadata
        ]
    )

    assert (
        preparer
        .reaction_templates_highlighted_image_grid(
            session
        )
        is None
    )


@pytest.mark.parametrize(
    "highlight_type",
    [
        "template",
        "edge",
        "initiators",
        "delete",
    ],
)
def test_reaction_templates_highlight_types_generate_image(
    tmp_path,
    highlight_type,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    metadata = make_image_metadata()

    session = SimpleNamespace(
        reaction_metadata=[
            metadata
        ]
    )

    image = (
        preparer
        .reaction_templates_highlighted_image_grid(
            session,
            highlight_type=highlight_type,
        )
    )

    assert image is not None
    assert image.size[0] > 0
    assert image.size[1] > 0


def test_reaction_templates_image_does_not_clear_original_atom_maps(
    tmp_path,
):
    preparer = PrepareReactions(
        make_session(tmp_path)
    )

    metadata = make_image_metadata()

    metadata.reactant_combined_RDmol.GetAtomWithIdx(
        0
    ).SetAtomMapNum(123)

    metadata.product_combined_RDmol.GetAtomWithIdx(
        0
    ).SetAtomMapNum(456)

    session = SimpleNamespace(
        reaction_metadata=[
            metadata
        ]
    )

    preparer.reaction_templates_highlighted_image_grid(
        session,
        highlight_type="template",
    )

    assert (
        metadata
        .reactant_combined_RDmol
        .GetAtomWithIdx(0)
        .GetAtomMapNum()
        == 123
    )

    assert (
        metadata
        .product_combined_RDmol
        .GetAtomWithIdx(0)
        .GetAtomMapNum()
        == 456
    )
