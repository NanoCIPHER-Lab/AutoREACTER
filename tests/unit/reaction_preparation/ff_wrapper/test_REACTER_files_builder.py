from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

import AutoREACTER.reaction_preparation.ff_wrapper.REACTER_files_builder as builder_module
from AutoREACTER.reaction_preparation.ff_wrapper.REACTER_files_builder import (
    REACTERFiles,
    REACTERFilesBuilder,
)


# =============================================================================
# Helpers
# =============================================================================


def make_monomer(
    *,
    monomer_id=1,
    data_id="data1",
    name="mma",
    lmp_molecule_file=None,
):
    return SimpleNamespace(
        id=monomer_id,
        data_id=data_id,
        name=name,
        lmp_molecule_file=lmp_molecule_file,
    )


def make_reaction_metadata(
    *,
    reaction_id=1,
    activity_stats=True,
    dataframe=None,
    delete_atom=False,
):
    if dataframe is None:
        dataframe = pd.DataFrame()

    return SimpleNamespace(
        reaction_id=reaction_id,
        activity_stats=activity_stats,
        reaction_dataframe=dataframe,
        delete_atom=delete_atom,
        map_file=None,
        map_file_with_delete_ids=None,
        pre_reaction_file=None,
        post_reaction_file=None,
    )


def make_session(
    tmp_path,
    *,
    force_field="PCFF",
    wildcards=False,
    deduplicate=True,
    monomers=None,
    reactions=None,
):
    inputs = SimpleNamespace(
        force_field=force_field,
        wildcards=wildcards,
        deduplicate_reaction_templates=deduplicate,
        monomers=list(
            monomers or []
        ),
    )

    return SimpleNamespace(
        inputs=inputs,
        reaction_metadata=list(
            reactions or []
        ),
        staging_dir=tmp_path / "staging",
        output_dir=tmp_path / "output",
        ff_files=None,
        reacter_files=None,
    )


def make_builder(
    tmp_path,
    **kwargs,
):
    session = make_session(
        tmp_path,
        **kwargs,
    )

    return (
        REACTERFilesBuilder(
            session
        ),
        session,
    )


def make_data_files(
    *,
    data_file,
    molecule_file,
):
    return SimpleNamespace(
        data_file=Path(
            data_file
        ),
        lmp_molecule_file=Path(
            molecule_file
        ),
    )


def make_template_ff_entry(
    *,
    reaction_id,
    pre_file,
    post_file,
):
    return SimpleNamespace(
        reaction_id=reaction_id,
        pre_reaction_file=make_data_files(
            data_file=Path(
                str(pre_file) + ".data"
            ),
            molecule_file=pre_file,
        ),
        post_reaction_file=make_data_files(
            data_file=Path(
                str(post_file) + ".data"
            ),
            molecule_file=post_file,
        ),
    )


# =============================================================================
# REACTERFiles
# =============================================================================


def test_reacter_files_stores_values(
    tmp_path,
):
    result = REACTERFiles(
        force_field_data=(
            tmp_path / "force_field.data"
        ),
        in_file=(
            tmp_path / "in.script"
        ),
        molecule_files=[
            "molecule"
        ],
        template_files=[
            "template"
        ],
    )

    assert result.force_field_data == (
        tmp_path / "force_field.data"
    )

    assert result.in_file == (
        tmp_path / "in.script"
    )

    assert result.molecule_files == [
        "molecule"
    ]

    assert result.template_files == [
        "template"
    ]


def test_reacter_files_uses_slots(
    tmp_path,
):
    result = REACTERFiles(
        force_field_data=(
            tmp_path / "ff.data"
        ),
        in_file=(
            tmp_path / "in.script"
        ),
        molecule_files=[],
        template_files=[],
    )

    with pytest.raises(
        AttributeError
    ):
        result.extra = 1


# =============================================================================
# Constructor
# =============================================================================


def test_constructor_stores_session(
    tmp_path,
):
    builder, session = (
        make_builder(
            tmp_path
        )
    )

    assert (
        builder.session
        is session
    )


def test_constructor_creates_reacter_cache(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    assert builder.cache_dir == (
        tmp_path
        / "staging"
        / "lunar"
        / "REACTER_files"
    )

    assert (
        builder.cache_dir.is_dir()
    )


def test_constructor_reads_force_field_options(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path,
            force_field="CVFF",
            wildcards=True,
            deduplicate=False,
        )
    )

    assert (
        builder.force_field
        == "CVFF"
    )

    assert (
        builder.wildcards
        is True
    )

    assert (
        builder
        .deduplicate_reaction_templates
        is False
    )


# =============================================================================
# _get_ending_integer
# =============================================================================


@pytest.mark.parametrize(
    "value, expected",
    [
        (
            "pre1",
            1,
        ),
        (
            "post27",
            27,
        ),
        (
            "data003",
            3,
        ),
        (
            "RXN_99",
            99,
        ),
        (
            "1.5",
            5,
        ),
    ],
)
def test_get_ending_integer(
    tmp_path,
    value,
    expected,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    assert (
        builder._get_ending_integer(
            value
        )
        == expected
    )


@pytest.mark.parametrize(
    "value",
    [
        "",
        "pre",
        "post_",
        "123abc",
        "reaction_test",
    ],
)
def test_get_ending_integer_none(
    tmp_path,
    value,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    assert (
        builder._get_ending_integer(
            value
        )
        is None
    )


# =============================================================================
# _ensure_dir
# =============================================================================


def test_ensure_dir_creates_directory(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    target = (
        tmp_path
        / "nested"
        / "directory"
    )

    result = builder._ensure_dir(
        str(target)
    )

    assert result == str(
        target
    )

    assert target.is_dir()


def test_ensure_dir_existing_directory_is_ok(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    target = (
        tmp_path / "existing"
    )

    target.mkdir()

    result = builder._ensure_dir(
        str(target)
    )

    assert result == str(
        target
    )


def test_ensure_dir_blocking_file_raises(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    target = (
        tmp_path / "blocked"
    )

    target.write_text(
        "file",
        encoding="utf-8",
    )

    with pytest.raises(
        FileExistsError,
        match="expected dir",
    ):
        builder._ensure_dir(
            str(target)
        )


def test_ensure_dir_can_remove_blocking_file(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    target = (
        tmp_path / "blocked"
    )

    target.write_text(
        "file",
        encoding="utf-8",
    )

    builder._ensure_dir(
        str(target),
        remove_blocking_file=True,
    )

    assert target.is_dir()


# =============================================================================
# _col_int_list
# =============================================================================


def test_col_int_list_missing_column_returns_empty(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    df = pd.DataFrame(
        {
            "other": [
                1,
                2,
            ]
        }
    )

    assert (
        builder._col_int_list(
            "initiators",
            df,
        )
        == []
    )


def test_col_int_list_empty_column_returns_empty(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    df = pd.DataFrame(
        {
            "initiators": [
                None,
                None,
            ]
        }
    )

    assert (
        builder._col_int_list(
            "initiators",
            df,
        )
        == []
    )


def test_col_int_list_converts_integer_like_values(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    df = pd.DataFrame(
        {
            "initiators": [
                1,
                2.0,
                "3",
                "4.0",
            ]
        }
    )

    assert (
        builder._col_int_list(
            "initiators",
            df,
        )
        == [
            1,
            2,
            3,
            4,
        ]
    )


def test_col_int_list_removes_duplicates_preserving_order(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    df = pd.DataFrame(
        {
            "edge_atoms": [
                5,
                2,
                5,
                3,
                2,
            ]
        }
    )

    assert (
        builder._col_int_list(
            "edge_atoms",
            df,
        )
        == [
            5,
            2,
            3,
        ]
    )


def test_col_int_list_ignores_nan(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    df = pd.DataFrame(
        {
            "byproduct_idx": [
                1,
                None,
                2,
            ]
        }
    )

    assert (
        builder._col_int_list(
            "byproduct_idx",
            df,
        )
        == [
            1,
            2,
        ]
    )


def test_col_int_list_non_numeric_raises(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    df = pd.DataFrame(
        {
            "initiators": [
                1,
                "bad",
            ]
        }
    )

    with pytest.raises(
        ValueError,
        match="contains non-numeric value",
    ):
        builder._col_int_list(
            "initiators",
            df,
        )


def test_col_int_list_non_integer_numeric_raises(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    df = pd.DataFrame(
        {
            "initiators": [
                1,
                2.5,
            ]
        }
    )

    with pytest.raises(
        ValueError,
        match="contains non-integer value",
    ):
        builder._col_int_list(
            "initiators",
            df,
        )


# =============================================================================
# _load_molecule_file
# =============================================================================


def test_load_molecule_file_finds_all_sections(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    path = (
        tmp_path / "template.molecule"
    )

    path.write_text(
        """Types

1 1 # c

Charges

1 0.0 # c

Coords

1 0.0 0.0 0.0 # c

Bonds

1 1 1 1 # c c

Angles

1 1 1 1 1 # c c c

Dihedrals

1 1 1 1 1 1 # c c c c

Impropers

1 1 1 1 1 1 # c c c c

""",
        encoding="utf-8",
    )

    result = (
        builder._load_molecule_file(
            path
        )
    )

    (
        lines,
        type_start,
        charge_start,
        coord_start,
        bond_start,
        angle_start,
        dihedral_start,
        improper_start,
    ) = result

    assert type_start == 2
    assert charge_start == 6
    assert coord_start == 10
    assert bond_start == 14
    assert angle_start == 18
    assert dihedral_start == 22
    assert improper_start == 26

    assert lines[0] == "Types"


def test_load_molecule_file_strips_line_whitespace(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    path = (
        tmp_path / "template.molecule"
    )

    path.write_text(
        """  Types  

  1 10 # c  

Charges

1 0.0 # c

Coords

1 0 0 0 # c

""",
        encoding="utf-8",
    )

    lines, *_ = (
        builder._load_molecule_file(
            path
        )
    )

    assert lines[0] == "Types"
    assert lines[2] == "1 10 # c"


def test_load_molecule_file_requires_types(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    path = (
        tmp_path / "bad.molecule"
    )

    path.write_text(
        """Charges

1 0 # c

Coords

1 0 0 0 # c
""",
        encoding="utf-8",
    )

    with pytest.raises(
        ValueError,
        match="Essential sections",
    ):
        builder._load_molecule_file(
            path
        )


def test_load_molecule_file_requires_charges(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    path = (
        tmp_path / "bad.molecule"
    )

    path.write_text(
        """Types

1 1 # c

Coords

1 0 0 0 # c
""",
        encoding="utf-8",
    )

    with pytest.raises(
        ValueError,
        match="Essential sections",
    ):
        builder._load_molecule_file(
            path
        )


def test_load_molecule_file_requires_coords(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    path = (
        tmp_path / "bad.molecule"
    )

    path.write_text(
        """Types

1 1 # c

Charges

1 0 # c
""",
        encoding="utf-8",
    )

    with pytest.raises(
        ValueError,
        match="Essential sections",
    ):
        builder._load_molecule_file(
            path
        )


def test_load_molecule_file_missing_optional_topology_sections_returns_none(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    path = (
        tmp_path / "minimal.molecule"
    )

    path.write_text(
        """Types

1 1 # c

Charges

1 0 # c

Coords

1 0 0 0 # c

""",
        encoding="utf-8",
    )

    result = (
        builder._load_molecule_file(
            path
        )
    )

    assert result[4] is None
    assert result[5] is None
    assert result[6] is None
    assert result[7] is None


# =============================================================================
# _molecule_file_format
# =============================================================================


def test_molecule_file_format_contains_counts(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path,
            force_field="PCFF",
        )
    )

    text = builder._molecule_file_format(
        5,
        4,
        3,
        2,
        1,
        "TYPES\n",
        "CHARGES\n",
        "COORDS\n",
        "BONDS\n",
        "ANGLES\n",
        "DIHEDRALS\n",
        "IMPROPERS\n",
        "template_pre_1.molecule",
    )

    assert "5 atoms" in text
    assert "4 bonds" in text
    assert "3 angles" in text
    assert "2 dihedrals" in text
    assert "1 impropers" in text


def test_molecule_file_format_contains_force_field_and_filename(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path,
            force_field="PCFF",
        )
    )

    text = builder._molecule_file_format(
        1,
        0,
        0,
        0,
        0,
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "template_pre_7.molecule",
    )

    assert (
        "template_pre_7.molecule"
        in text
    )

    assert (
        "PCFF forcefield"
        in text
    )


def test_molecule_file_format_section_order(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    text = builder._molecule_file_format(
        1,
        1,
        1,
        1,
        1,
        "TYPE_DATA\n",
        "CHARGE_DATA\n",
        "COORD_DATA\n",
        "BOND_DATA\n",
        "ANGLE_DATA\n",
        "DIHEDRAL_DATA\n",
        "IMPROPER_DATA\n",
        "test.molecule",
    )

    assert (
        text.index("Types")
        < text.index("Charges")
        < text.index("Coords")
        < text.index("Bonds")
        < text.index("Angles")
        < text.index("Dihedrals")
        < text.index("Impropers")
    )


# =============================================================================
# _molecule_file_preparation
# =============================================================================


def test_molecule_file_preparation_calls_modifier_pipeline(
    tmp_path,
    monkeypatch,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    fake_df = pd.DataFrame(
        {
            "atom_index": [
                1,
                2,
            ],
            "new_atom_index": [
                1,
                2,
            ],
        }
    )

    monkeypatch.setattr(
        builder,
        "_load_molecule_file",
        lambda path: (
            ["fake-lines"],
            1,
            2,
            3,
            4,
            5,
            6,
            7,
        ),
    )

    events = []

    def fake_types(
        lines,
        indexes,
        start,
    ):
        events.append(
            (
                "types",
                indexes,
                start,
            )
        )

        return (
            fake_df,
            "TYPES\n",
            2,
            {
                10: 1,
                20: 2,
            },
            False,
        )

    monkeypatch.setattr(
        builder_module,
        "modify_types",
        fake_types,
    )

    monkeypatch.setattr(
        builder_module,
        "modify_charges",
        lambda lines, df, start:
        (
            events.append(
                (
                    "charges",
                    start,
                )
            )
            or "CHARGES\n"
        ),
    )

    monkeypatch.setattr(
        builder_module,
        "modify_coords",
        lambda lines, df, start:
        (
            events.append(
                (
                    "coords",
                    start,
                )
            )
            or "COORDS\n"
        ),
    )

    monkeypatch.setattr(
        builder_module,
        "modify_bonds",
        lambda lines, df, start, legacy_mode:
        (
            events.append(
                (
                    "bonds",
                    start,
                    legacy_mode,
                )
            )
            or (
                "BONDS\n",
                1,
            )
        ),
    )

    monkeypatch.setattr(
        builder_module,
        "modify_angles",
        lambda lines, df, start, legacy_mode:
        (
            events.append(
                (
                    "angles",
                    start,
                    legacy_mode,
                )
            )
            or (
                "ANGLES\n",
                2,
            )
        ),
    )

    monkeypatch.setattr(
        builder_module,
        "modify_dihedrals",
        lambda lines, df, start, legacy_mode:
        (
            events.append(
                (
                    "dihedrals",
                    start,
                    legacy_mode,
                )
            )
            or (
                "DIHEDRALS\n",
                3,
            )
        ),
    )

    monkeypatch.setattr(
        builder_module,
        "modify_impropers",
        lambda lines, df, start, legacy_mode:
        (
            events.append(
                (
                    "impropers",
                    start,
                    legacy_mode,
                )
            )
            or (
                "IMPROPERS\n",
                4,
            )
        ),
    )

    monkeypatch.setattr(
        builder,
        "_molecule_file_format",
        lambda *args:
        "FINAL-MOLECULE",
    )

    result, mapping = (
        builder
        ._molecule_file_preparation(
            tmp_path / "input.molecule",
            [
                10,
                20,
            ],
            "template_pre_1.molecule",
        )
    )

    assert result == (
        "FINAL-MOLECULE"
    )

    assert mapping == {
        10: 1,
        20: 2,
    }

    assert events == [
        (
            "types",
            [
                10,
                20,
            ],
            1,
        ),
        (
            "charges",
            2,
        ),
        (
            "coords",
            3,
        ),
        (
            "bonds",
            4,
            False,
        ),
        (
            "angles",
            5,
            False,
        ),
        (
            "dihedrals",
            6,
            False,
        ),
        (
            "impropers",
            7,
            False,
        ),
    ]


# =============================================================================
# _map_file_write
# =============================================================================


def test_map_file_requires_exactly_two_initiators(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    with pytest.raises(
        ValueError,
        match="Expected exactly 2 initiator atoms",
    ):
        builder._map_file_write(
            reactant_to_product={
                0: 0,
            },
            initiator_atoms=[
                0,
            ],
            edge_atoms=[],
            delete_ids=[],
            file_name="RXN_1.map",
        )


def test_standard_map_has_no_delete_ids_section(
    tmp_path,
):
    """
    Important AutoREACTER contract:

    The normal RXN_N.map is the standard map and does NOT contain DeleteIDs.
    Delete IDs belong only in the optional supplementary map.
    """
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    text = builder._map_file_write(
        reactant_to_product={
            0: 1,
            1: 0,
        },
        initiator_atoms=[
            0,
            1,
        ],
        edge_atoms=[
            2,
        ],
        delete_ids=[],
        file_name="RXN_1.map",
    )

    assert "deleteIDs" not in text
    assert "\nDeleteIDs\n" not in text


def test_map_file_with_delete_ids_contains_optional_section(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    text = builder._map_file_write(
        reactant_to_product={
            0: 0,
            1: 1,
        },
        initiator_atoms=[
            0,
            1,
        ],
        edge_atoms=[],
        delete_ids=[
            3,
            2,
        ],
        file_name=(
            "RXN_1_with_delete_ids.map"
        ),
    )

    assert "2 deleteIDs" in text

    assert (
        "\nDeleteIDs\n\n"
        in text
    )

    # Written as 1-based and sorted.
    delete_block = (
        text.split(
            "DeleteIDs\n\n",
            1,
        )[1]
    )

    assert delete_block.splitlines() == [
        "3",
        "4",
    ]


def test_map_file_uses_one_based_indices(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    text = builder._map_file_write(
        reactant_to_product={
            0: 2,
            4: 1,
        },
        initiator_atoms=[
            0,
            4,
        ],
        edge_atoms=[
            3,
        ],
        delete_ids=[],
        file_name="RXN_5.map",
    )

    assert (
        "\nInitiatorIDs\n\n1\n5\n"
        in text
    )

    assert (
        "\nEdgeIDs\n\n4\n"
        in text
    )

    assert "1     3" in text
    assert "5     2" in text


def test_map_file_sorts_initiators_edges_and_equivalences(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    text = builder._map_file_write(
        reactant_to_product={
            5: 8,
            1: 2,
            3: 4,
        },
        initiator_atoms=[
            5,
            1,
        ],
        edge_atoms=[
            8,
            2,
            4,
        ],
        delete_ids=[],
        file_name="RXN_1.map",
    )

    initiator_block = (
        text.split(
            "InitiatorIDs\n\n",
            1,
        )[1]
        .split(
            "\nEdgeIDs",
            1,
        )[0]
        .strip()
        .splitlines()
    )

    assert initiator_block == [
        "2",
        "6",
    ]

    edge_block = (
        text.split(
            "EdgeIDs\n\n",
            1,
        )[1]
        .split(
            "\nEquivalences",
            1,
        )[0]
        .strip()
        .splitlines()
    )

    assert edge_block == [
        "3",
        "5",
        "9",
    ]


# =============================================================================
# _build_bond_react_templates
# =============================================================================


def test_build_templates_missing_post_entry_raises(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    pre = (
        tmp_path / "pre.molecule"
    )

    pre.write_text(
        "pre",
        encoding="utf-8",
    )

    with pytest.raises(
        ValueError,
        match="Corresponding post-reaction file",
    ):
        builder._build_bond_react_templates(
            file_dict={
                "pre_1": pre,
            },
            reactant_to_product={
                0: 0,
            },
            initiator_atoms=[
                0,
                1,
            ],
            edge_atoms=[],
            delete_ids=[],
        )


def test_build_templates_missing_pre_file_raises(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    post = (
        tmp_path / "post.molecule"
    )

    post.write_text(
        "post",
        encoding="utf-8",
    )

    with pytest.raises(
        FileNotFoundError,
        match="Missing pre file",
    ):
        builder._build_bond_react_templates(
            file_dict={
                "pre_1": (
                    tmp_path
                    / "missing.molecule"
                ),
                "post_1": post,
            },
            reactant_to_product={
                0: 0,
            },
            initiator_atoms=[
                0,
                1,
            ],
            edge_atoms=[],
            delete_ids=[],
        )


def test_build_templates_missing_post_file_raises(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    pre = (
        tmp_path / "pre.molecule"
    )

    pre.write_text(
        "pre",
        encoding="utf-8",
    )

    with pytest.raises(
        FileNotFoundError,
        match="Missing post file",
    ):
        builder._build_bond_react_templates(
            file_dict={
                "pre_1": pre,
                "post_1": (
                    tmp_path
                    / "missing.molecule"
                ),
            },
            reactant_to_product={
                0: 0,
            },
            initiator_atoms=[
                0,
                1,
            ],
            edge_atoms=[],
            delete_ids=[],
        )


def test_build_templates_writes_standard_pre_post_and_map(
    tmp_path,
    monkeypatch,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    pre = (
        tmp_path / "pre_source.molecule"
    )

    post = (
        tmp_path / "post_source.molecule"
    )

    pre.write_text(
        "pre",
        encoding="utf-8",
    )

    post.write_text(
        "post",
        encoding="utf-8",
    )

    def fake_prepare(
        path,
        indexes,
        file_name,
    ):
        if Path(path) == pre:
            return (
                "PRE-CONTENT",
                {
                    1: 1,
                    2: 2,
                },
            )

        return (
            "POST-CONTENT",
            {
                1: 1,
                2: 2,
            },
        )

    monkeypatch.setattr(
        builder,
        "_molecule_file_preparation",
        fake_prepare,
    )

    monkeypatch.setattr(
        builder,
        "_map_file_write",
        lambda *args, **kwargs:
        "STANDARD-MAP",
    )

    pre_out, post_out, map_path = (
        builder
        ._build_bond_react_templates(
            file_dict={
                "pre_1": pre,
                "post_1": post,
            },
            reactant_to_product={
                0: 0,
                1: 1,
            },
            initiator_atoms=[
                0,
                1,
            ],
            edge_atoms=[],
            delete_ids=[],
        )
    )

    assert Path(
        pre_out
    ).read_text(
        encoding="utf-8"
    ) == "PRE-CONTENT"

    assert Path(
        post_out
    ).read_text(
        encoding="utf-8"
    ) == "POST-CONTENT"

    assert Path(
        map_path
    ).read_text(
        encoding="utf-8"
    ) == "STANDARD-MAP"

    assert Path(
        pre_out
    ).name == (
        "template_pre_1.molecule"
    )

    assert Path(
        post_out
    ).name == (
        "template_post_1.molecule"
    )

    assert Path(
        map_path
    ).name == "RXN_1.map"


def test_build_templates_standard_map_suppresses_delete_ids(
    tmp_path,
    monkeypatch,
):
    """
    Even when delete atoms exist, the PRIMARY map must be generated
    with delete_ids=[].
    """
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    pre = tmp_path / "pre.molecule"
    post = tmp_path / "post.molecule"

    pre.write_text("pre")
    post.write_text("post")

    monkeypatch.setattr(
        builder,
        "_molecule_file_preparation",
        lambda path, indexes, file_name:
        (
            "CONTENT",
            {
                1: 1,
                2: 2,
                3: 3,
            },
        ),
    )

    calls = []

    def fake_map(
        mapping,
        initiators,
        edges,
        deletes,
        file_name,
    ):
        calls.append(
            (
                mapping,
                initiators,
                edges,
                deletes,
                file_name,
            )
        )

        return file_name

    monkeypatch.setattr(
        builder,
        "_map_file_write",
        fake_map,
    )

    builder._build_bond_react_templates(
        file_dict={
            "pre_1": pre,
            "post_1": post,
        },
        reactant_to_product={
            0: 0,
            1: 1,
            2: 2,
        },
        initiator_atoms=[
            0,
            1,
        ],
        edge_atoms=[
            2,
        ],
        delete_ids=[
            2,
        ],
    )

    # First map is always standard.
    assert calls[0][3] == []

    assert calls[0][4] == (
        "RXN_1.map"
    )


def test_build_templates_optional_delete_map_written_only_when_needed(
    tmp_path,
    monkeypatch,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    pre = tmp_path / "pre.molecule"
    post = tmp_path / "post.molecule"

    pre.write_text("pre")
    post.write_text("post")

    monkeypatch.setattr(
        builder,
        "_molecule_file_preparation",
        lambda path, indexes, file_name:
        (
            "CONTENT",
            {
                1: 1,
                2: 2,
                3: 3,
            },
        ),
    )

    calls = []

    def fake_map(
        mapping,
        initiators,
        edges,
        deletes,
        file_name,
    ):
        calls.append(
            (
                list(deletes),
                file_name,
            )
        )

        return file_name

    monkeypatch.setattr(
        builder,
        "_map_file_write",
        fake_map,
    )

    builder._build_bond_react_templates(
        file_dict={
            "pre_4": pre,
            "post_4": post,
        },
        reactant_to_product={
            0: 0,
            1: 1,
            2: 2,
        },
        initiator_atoms=[
            0,
            1,
        ],
        edge_atoms=[],
        delete_ids=[
            2,
        ],
    )

    assert calls == [
        (
            [],
            "RXN_4.map",
        ),
        (
            [
                2,
            ],
            "RXN_4_with_delete_ids.map",
        ),
    ]

    assert (
        builder.cache_dir
        / "RXN_4_with_delete_ids.map"
    ).is_file()


def test_build_templates_does_not_write_optional_delete_map_when_no_delete_ids(
    tmp_path,
    monkeypatch,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    pre = tmp_path / "pre.molecule"
    post = tmp_path / "post.molecule"

    pre.write_text("pre")
    post.write_text("post")

    monkeypatch.setattr(
        builder,
        "_molecule_file_preparation",
        lambda path, indexes, file_name:
        (
            "CONTENT",
            {
                1: 1,
                2: 2,
            },
        ),
    )

    monkeypatch.setattr(
        builder,
        "_map_file_write",
        lambda *args, **kwargs:
        "MAP",
    )

    builder._build_bond_react_templates(
        file_dict={
            "pre_2": pre,
            "post_2": post,
        },
        reactant_to_product={
            0: 0,
            1: 1,
        },
        initiator_atoms=[
            0,
            1,
        ],
        edge_atoms=[],
        delete_ids=[],
    )

    assert not (
        builder.cache_dir
        / "RXN_2_with_delete_ids.map"
    ).exists()


def test_build_templates_reindexes_mapping_into_trimmed_template_space(
    tmp_path,
    monkeypatch,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    pre = tmp_path / "pre.molecule"
    post = tmp_path / "post.molecule"

    pre.write_text("pre")
    post.write_text("post")

    def fake_prepare(
        path,
        indexes,
        file_name,
    ):
        if Path(path) == pre:
            return (
                "PRE",
                {
                    1: 2,
                    2: 1,
                },
            )

        return (
            "POST",
            {
                1: 2,
                3: 1,
            },
        )

    monkeypatch.setattr(
        builder,
        "_molecule_file_preparation",
        fake_prepare,
    )

    mappings = []

    def fake_map(
        mapping,
        initiators,
        edges,
        deletes,
        file_name,
    ):
        mappings.append(
            dict(mapping)
        )

        return "MAP"

    monkeypatch.setattr(
        builder,
        "_map_file_write",
        fake_map,
    )

    builder._build_bond_react_templates(
        file_dict={
            "pre_1": pre,
            "post_1": post,
        },
        reactant_to_product={
            0: 2,
            1: 0,
        },
        initiator_atoms=[
            0,
            1,
        ],
        edge_atoms=[],
        delete_ids=[],
    )

    # full 0 -> full 2
    # reactant old1 -> template2 -> zero-based1
    # product  old3 -> template1 -> zero-based0
    #
    # full 1 -> full 0
    # reactant old2 -> template1 -> zero-based0
    # product  old1 -> template2 -> zero-based1
    assert mappings[0] == {
        1: 0,
        0: 1,
    }


def test_build_templates_missing_filtered_mapping_raises(
    tmp_path,
    monkeypatch,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    pre = tmp_path / "pre.molecule"
    post = tmp_path / "post.molecule"

    pre.write_text("pre")
    post.write_text("post")

    monkeypatch.setattr(
        builder,
        "_molecule_file_preparation",
        lambda path, indexes, file_name:
        (
            "CONTENT",
            {
                1: 1,
            },
        ),
    )

    with pytest.raises(
        ValueError,
        match=(
            "Template mapping indices "
            "missing after filtering"
        ),
    ):
        builder._build_bond_react_templates(
            file_dict={
                "pre_1": pre,
                "post_1": post,
            },
            reactant_to_product={
                0: 0,
                1: 1,
            },
            initiator_atoms=[
                0,
                1,
            ],
            edge_atoms=[],
            delete_ids=[],
        )


# =============================================================================
# _copy_lunar_files_to_cache
# =============================================================================


def test_copy_lunar_files_copies_force_field_and_input(
    tmp_path,
):
    monomer = make_monomer()

    builder, session = (
        make_builder(
            tmp_path,
            monomers=[
                monomer
            ],
        )
    )

    ff_src = (
        tmp_path / "force_field.data"
    )

    in_src = (
        tmp_path
        / "in.create_atoms.script"
    )

    ff_src.write_text(
        "FF",
        encoding="utf-8",
    )

    in_src.write_text(
        "INPUT",
        encoding="utf-8",
    )

    ff_files = SimpleNamespace(
        force_field_data=ff_src,
        in_file=in_src,
        molecule_files=[],
    )

    ff_dest, in_dest = (
        builder
        ._copy_lunar_files_to_cache(
            ff_files
        )
    )

    assert ff_dest == (
        builder.cache_dir
        / "force_field.data"
    )

    assert in_dest == (
        builder.cache_dir
        / "in.create_atoms.script"
    )

    assert ff_dest.read_text() == "FF"
    assert in_dest.read_text() == "INPUT"


def test_copy_lunar_files_attaches_molecule_by_name(
    tmp_path,
):
    monomer = make_monomer(
        monomer_id=7,
        data_id="data7",
        name="mma",
    )

    builder, _ = (
        make_builder(
            tmp_path,
            monomers=[
                monomer
            ],
        )
    )

    ff_src = tmp_path / "ff.data"
    in_src = tmp_path / "in.script"
    mol_src = tmp_path / "mma.lmpmol"

    ff_src.write_text("ff")
    in_src.write_text("in")
    mol_src.write_text("mol")

    ff_files = SimpleNamespace(
        force_field_data=ff_src,
        in_file=in_src,
        molecule_files=[
            SimpleNamespace(
                id="mma",
                molecule_files=SimpleNamespace(
                    lmp_molecule_file=mol_src,
                ),
            )
        ],
    )

    builder._copy_lunar_files_to_cache(
        ff_files
    )

    assert (
        monomer.lmp_molecule_file
        == builder.cache_dir
        / "mma.molecule"
    )

    assert (
        monomer
        .lmp_molecule_file
        .read_text()
        == "mol"
    )


@pytest.mark.parametrize(
    "lookup_id",
    [
        "7",
        "data7",
        "mma",
    ],
)
def test_copy_lunar_files_matches_monomer_by_id_data_id_or_name(
    tmp_path,
    lookup_id,
):
    monomer = make_monomer(
        monomer_id=7,
        data_id="data7",
        name="mma",
    )

    builder, _ = (
        make_builder(
            tmp_path,
            monomers=[
                monomer
            ],
        )
    )

    ff_src = tmp_path / "ff.data"
    in_src = tmp_path / "in.script"
    mol_src = tmp_path / "mol.lmpmol"

    ff_src.write_text("ff")
    in_src.write_text("in")
    mol_src.write_text("mol")

    ff_files = SimpleNamespace(
        force_field_data=ff_src,
        in_file=in_src,
        molecule_files=[
            SimpleNamespace(
                id=lookup_id,
                molecule_files=SimpleNamespace(
                    lmp_molecule_file=mol_src,
                ),
            )
        ],
    )

    builder._copy_lunar_files_to_cache(
        ff_files
    )

    assert (
        monomer.lmp_molecule_file
        is not None
    )


def test_copy_lunar_force_field_failure_wrapped(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    ff_files = SimpleNamespace(
        force_field_data=(
            tmp_path / "missing.data"
        ),
        in_file=(
            tmp_path / "missing.script"
        ),
        molecule_files=[],
    )

    with pytest.raises(
        FileNotFoundError,
        match=(
            "Failed to copy force field data"
        ),
    ):
        builder._copy_lunar_files_to_cache(
            ff_files
        )


# =============================================================================
# _copy_path_to_final
# =============================================================================


def test_copy_path_to_final_none_returns_none(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    assert (
        builder._copy_path_to_final(
            None,
            tmp_path / "final",
        )
        is None
    )


def test_copy_path_to_final_missing_file_returns_none(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    result = (
        builder._copy_path_to_final(
            tmp_path / "missing.map",
            tmp_path / "final",
        )
    )

    assert result is None


def test_copy_path_to_final_copies_existing_file(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    src = tmp_path / "RXN_1.map"

    src.write_text(
        "MAP",
        encoding="utf-8",
    )

    final_dir = (
        tmp_path / "final"
    )

    result = (
        builder._copy_path_to_final(
            src,
            final_dir,
        )
    )

    assert result == (
        final_dir / "RXN_1.map"
    )

    assert result.read_text() == "MAP"


def test_copy_path_to_final_same_source_and_destination(
    tmp_path,
):
    builder, _ = (
        make_builder(
            tmp_path
        )
    )

    final_dir = (
        tmp_path / "final"
    )

    final_dir.mkdir()

    src = (
        final_dir / "RXN_1.map"
    )

    src.write_text(
        "MAP",
        encoding="utf-8",
    )

    result = (
        builder._copy_path_to_final(
            src,
            final_dir,
        )
    )

    assert result == src
    assert result.read_text() == "MAP"


# =============================================================================
# molecule_template_preparation - core orchestration
# =============================================================================


def test_molecule_template_preparation_builds_template_mapping_from_dataframe(
    tmp_path,
    monkeypatch,
):
    df = pd.DataFrame(
        {
            "template_reactant_idx": [
                0,
                1,
            ],
            "template_product_idx": [
                2,
                0,
            ],
            "initiators": [
                0,
                1,
            ],
            "edge_atoms": [
                3,
                None,
            ],
            "byproduct_idx": [
                None,
                None,
            ],
        }
    )

    metadata = (
        make_reaction_metadata(
            reaction_id=5,
            dataframe=df,
        )
    )

    monomer = make_monomer()

    builder, session = (
        make_builder(
            tmp_path,
            monomers=[
                monomer
            ],
            reactions=[
                metadata
            ],
        )
    )

    pre_source = (
        tmp_path / "pre5.lmpmol"
    )

    post_source = (
        tmp_path / "post5.lmpmol"
    )

    pre_source.write_text("pre")
    post_source.write_text("post")

    session.ff_files = (
        SimpleNamespace(
            template_files=[
                make_template_ff_entry(
                    reaction_id=5,
                    pre_file=pre_source,
                    post_file=post_source,
                )
            ],
        )
    )

    ff_cache = (
        builder.cache_dir
        / "force_field.data"
    )

    in_cache = (
        builder.cache_dir
        / "in.script"
    )

    ff_cache.write_text("ff")
    in_cache.write_text("in")

    monkeypatch.setattr(
        builder,
        "_copy_lunar_files_to_cache",
        lambda ff_files:
        (
            ff_cache,
            in_cache,
        ),
    )

    pre_out = (
        builder.cache_dir
        / "template_pre_5.molecule"
    )

    post_out = (
        builder.cache_dir
        / "template_post_5.molecule"
    )

    map_out = (
        builder.cache_dir
        / "RXN_5.map"
    )

    for path in (
        pre_out,
        post_out,
        map_out,
    ):
        path.write_text(
            path.name,
            encoding="utf-8",
        )

    calls = []

    def fake_build(
        *,
        file_dict,
        reactant_to_product,
        initiator_atoms,
        edge_atoms,
        delete_ids,
    ):
        calls.append(
            {
                "file_dict": file_dict,
                "mapping": reactant_to_product,
                "initiators": initiator_atoms,
                "edges": edge_atoms,
                "deletes": delete_ids,
            }
        )

        return (
            pre_out,
            post_out,
            map_out,
        )

    monkeypatch.setattr(
        builder,
        "_build_bond_react_templates",
        fake_build,
    )

    class FakeDetector:
        def compare_lammps_templates(
            self,
            template_files,
            wildcards,
        ):
            return template_files

    monkeypatch.setattr(
        builder_module,
        "DeduplicationDetector",
        FakeDetector,
    )

    builder.molecule_template_preparation(
        session
    )

    assert calls[0][
        "mapping"
    ] == {
        0: 2,
        1: 0,
    }

    assert calls[0][
        "initiators"
    ] == [
        0,
        1,
    ]

    assert calls[0][
        "edges"
    ] == [
        3,
    ]

    assert calls[0][
        "deletes"
    ] == []


def test_delete_atom_false_does_not_pass_byproducts_as_delete_ids(
    tmp_path,
    monkeypatch,
):
    df = pd.DataFrame(
        {
            "template_reactant_idx": [
                0,
                1,
            ],
            "template_product_idx": [
                0,
                1,
            ],
            "initiators": [
                0,
                1,
            ],
            "byproduct_idx": [
                5,
                6,
            ],
        }
    )

    metadata = make_reaction_metadata(
        reaction_id=1,
        dataframe=df,
        delete_atom=False,
    )

    builder, session = (
        make_builder(
            tmp_path,
            reactions=[
                metadata
            ],
        )
    )

    pre = tmp_path / "pre.lmpmol"
    post = tmp_path / "post.lmpmol"

    pre.write_text("pre")
    post.write_text("post")

    session.ff_files = (
        SimpleNamespace(
            template_files=[
                make_template_ff_entry(
                    reaction_id=1,
                    pre_file=pre,
                    post_file=post,
                )
            ],
        )
    )

    ff_cache = builder.cache_dir / "ff.data"
    in_cache = builder.cache_dir / "in.script"

    ff_cache.write_text("ff")
    in_cache.write_text("in")

    monkeypatch.setattr(
        builder,
        "_copy_lunar_files_to_cache",
        lambda files:
        (
            ff_cache,
            in_cache,
        ),
    )

    captured = []

    def fake_build(**kwargs):
        captured.append(
            kwargs["delete_ids"]
        )

        pre_out = (
            builder.cache_dir
            / "template_pre_1.molecule"
        )

        post_out = (
            builder.cache_dir
            / "template_post_1.molecule"
        )

        map_out = (
            builder.cache_dir
            / "RXN_1.map"
        )

        for path in (
            pre_out,
            post_out,
            map_out,
        ):
            path.write_text("x")

        return (
            pre_out,
            post_out,
            map_out,
        )

    monkeypatch.setattr(
        builder,
        "_build_bond_react_templates",
        fake_build,
    )

    monkeypatch.setattr(
        builder_module,
        "DeduplicationDetector",
        lambda:
        SimpleNamespace(
            compare_lammps_templates=(
                lambda template_files, wildcards:
                template_files
            )
        ),
    )

    builder.molecule_template_preparation(
        session
    )

    assert captured == [
        []
    ]


def test_delete_atom_true_passes_byproducts_as_optional_delete_ids(
    tmp_path,
    monkeypatch,
):
    df = pd.DataFrame(
        {
            "template_reactant_idx": [
                0,
                1,
            ],
            "template_product_idx": [
                0,
                1,
            ],
            "initiators": [
                0,
                1,
            ],
            "byproduct_idx": [
                5,
                6,
            ],
        }
    )

    metadata = make_reaction_metadata(
        reaction_id=1,
        dataframe=df,
        delete_atom=True,
    )

    builder, session = (
        make_builder(
            tmp_path,
            reactions=[
                metadata
            ],
        )
    )

    pre = tmp_path / "pre.lmpmol"
    post = tmp_path / "post.lmpmol"

    pre.write_text("pre")
    post.write_text("post")

    session.ff_files = (
        SimpleNamespace(
            template_files=[
                make_template_ff_entry(
                    reaction_id=1,
                    pre_file=pre,
                    post_file=post,
                )
            ],
        )
    )

    ff_cache = builder.cache_dir / "ff.data"
    in_cache = builder.cache_dir / "in.script"

    ff_cache.write_text("ff")
    in_cache.write_text("in")

    monkeypatch.setattr(
        builder,
        "_copy_lunar_files_to_cache",
        lambda files:
        (
            ff_cache,
            in_cache,
        ),
    )

    captured = []

    def fake_build(**kwargs):
        captured.append(
            kwargs["delete_ids"]
        )

        pre_out = (
            builder.cache_dir
            / "template_pre_1.molecule"
        )

        post_out = (
            builder.cache_dir
            / "template_post_1.molecule"
        )

        map_out = (
            builder.cache_dir
            / "RXN_1.map"
        )

        for path in (
            pre_out,
            post_out,
            map_out,
        ):
            path.write_text("x")

        return (
            pre_out,
            post_out,
            map_out,
        )

    monkeypatch.setattr(
        builder,
        "_build_bond_react_templates",
        fake_build,
    )

    monkeypatch.setattr(
        builder_module,
        "DeduplicationDetector",
        lambda:
        SimpleNamespace(
            compare_lammps_templates=(
                lambda template_files, wildcards:
                template_files
            )
        ),
    )

    builder.molecule_template_preparation(
        session
    )

    assert captured == [
        [
            5,
            6,
        ]
    ]


def test_optional_delete_map_is_not_required(
    tmp_path,
    monkeypatch,
):
    """
    Important contract:

    A reaction is valid with only RXN_N.map.
    RXN_N_with_delete_ids.map is supplementary and may remain None.
    """
    df = pd.DataFrame(
        {
            "template_reactant_idx": [
                0,
                1,
            ],
            "template_product_idx": [
                0,
                1,
            ],
            "initiators": [
                0,
                1,
            ],
        }
    )

    metadata = make_reaction_metadata(
        reaction_id=1,
        dataframe=df,
        delete_atom=False,
    )

    builder, session = (
        make_builder(
            tmp_path,
            reactions=[
                metadata
            ],
        )
    )

    pre_source = tmp_path / "pre.lmpmol"
    post_source = tmp_path / "post.lmpmol"

    pre_source.write_text("pre")
    post_source.write_text("post")

    session.ff_files = (
        SimpleNamespace(
            template_files=[
                make_template_ff_entry(
                    reaction_id=1,
                    pre_file=pre_source,
                    post_file=post_source,
                )
            ],
        )
    )

    ff_cache = builder.cache_dir / "ff.data"
    in_cache = builder.cache_dir / "in.script"
    pre_out = builder.cache_dir / "template_pre_1.molecule"
    post_out = builder.cache_dir / "template_post_1.molecule"
    map_out = builder.cache_dir / "RXN_1.map"

    for path in (
        ff_cache,
        in_cache,
        pre_out,
        post_out,
        map_out,
    ):
        path.write_text("x")

    monkeypatch.setattr(
        builder,
        "_copy_lunar_files_to_cache",
        lambda files:
        (
            ff_cache,
            in_cache,
        ),
    )

    monkeypatch.setattr(
        builder,
        "_build_bond_react_templates",
        lambda **kwargs:
        (
            pre_out,
            post_out,
            map_out,
        ),
    )

    monkeypatch.setattr(
        builder_module,
        "DeduplicationDetector",
        lambda:
        SimpleNamespace(
            compare_lammps_templates=(
                lambda template_files, wildcards:
                template_files
            )
        ),
    )

    builder.molecule_template_preparation(
        session
    )

    assert (
        metadata.map_file
        is not None
    )

    assert (
        metadata.pre_reaction_file
        is not None
    )

    assert (
        metadata.post_reaction_file
        is not None
    )

    assert (
        metadata.map_file_with_delete_ids
        is None
    )

    assert len(
        session
        .reacter_files
        .template_files
    ) == 1


# =============================================================================
# Deduplication routing
# =============================================================================


def test_template_deduplication_enabled_calls_compare_lammps_templates(
    tmp_path,
    monkeypatch,
):
    builder, session = (
        make_builder(
            tmp_path,
            deduplicate=True,
            wildcards=True,
        )
    )

    builder_reaction = (
        make_reaction_metadata(
            reaction_id=1,
            activity_stats=True,
        )
    )

    builder_reaction.map_file = (
        tmp_path / "RXN_1.map"
    )

    builder_reaction.pre_reaction_file = (
        tmp_path / "pre.molecule"
    )

    builder_reaction.post_reaction_file = (
        tmp_path / "post.molecule"
    )

    session.reaction_metadata = [
        builder_reaction
    ]

    session.ff_files = SimpleNamespace(
        template_files=[],
    )

    ff_cache = builder.cache_dir / "ff.data"
    in_cache = builder.cache_dir / "in.script"

    ff_cache.write_text("ff")
    in_cache.write_text("in")

    monkeypatch.setattr(
        builder,
        "_copy_lunar_files_to_cache",
        lambda files:
        (
            ff_cache,
            in_cache,
        ),
    )

    calls = []

    class FakeDetector:
        def compare_lammps_templates(
            self,
            template_files,
            wildcards,
        ):
            calls.append(
                (
                    list(template_files),
                    wildcards,
                )
            )

            return [
                "DEDUPED"
            ]

    monkeypatch.setattr(
        builder_module,
        "DeduplicationDetector",
        FakeDetector,
    )

    # Keep already-set reaction paths intact.
    monkeypatch.setattr(
        builder,
        "_copy_path_to_final",
        lambda src, final:
        Path(src)
        if src is not None
        else None,
    )

    builder.molecule_template_preparation(
        session
    )

    assert calls == [
        (
            [
                builder_reaction
            ],
            True,
        )
    ]

    assert (
        session
        .reacter_files
        .template_files
        == [
            "DEDUPED"
        ]
    )


def test_deduplication_disabled_skips_compare(
    tmp_path,
    monkeypatch,
    capsys,
):
    builder, session = (
        make_builder(
            tmp_path,
            deduplicate=False,
            wildcards=False,
        )
    )

    reaction = (
        make_reaction_metadata(
            reaction_id=1,
            activity_stats=True,
        )
    )

    reaction.map_file = (
        tmp_path / "RXN_1.map"
    )

    reaction.pre_reaction_file = (
        tmp_path / "pre.molecule"
    )

    reaction.post_reaction_file = (
        tmp_path / "post.molecule"
    )

    session.reaction_metadata = [
        reaction
    ]

    session.ff_files = (
        SimpleNamespace(
            template_files=[],
        )
    )

    ff_cache = builder.cache_dir / "ff.data"
    in_cache = builder.cache_dir / "in.script"

    ff_cache.write_text("ff")
    in_cache.write_text("in")

    monkeypatch.setattr(
        builder,
        "_copy_lunar_files_to_cache",
        lambda files:
        (
            ff_cache,
            in_cache,
        ),
    )

    class FakeDetector:
        def compare_lammps_templates(
            self,
            **kwargs,
        ):
            pytest.fail(
                "deduplication must not run"
            )

        def write_wildcard_maps(
            self,
            **kwargs,
        ):
            pytest.fail(
                "wildcard rewrite must not run"
            )

    monkeypatch.setattr(
        builder_module,
        "DeduplicationDetector",
        FakeDetector,
    )

    monkeypatch.setattr(
        builder,
        "_copy_path_to_final",
        lambda src, final:
        Path(src)
        if src is not None
        else None,
    )

    builder.molecule_template_preparation(
        session
    )

    assert (
        "Skipping LAMMPS reaction-template deduplication"
        in capsys.readouterr().out
    )

    assert (
        session
        .reacter_files
        .template_files
        == [
            reaction
        ]
    )


def test_wildcards_without_deduplication_calls_write_wildcard_maps(
    tmp_path,
    monkeypatch,
):
    builder, session = (
        make_builder(
            tmp_path,
            deduplicate=False,
            wildcards=True,
        )
    )

    reaction = (
        make_reaction_metadata(
            reaction_id=1,
            activity_stats=True,
        )
    )

    reaction.map_file = (
        tmp_path / "RXN_1.map"
    )

    reaction.pre_reaction_file = (
        tmp_path / "pre.molecule"
    )

    reaction.post_reaction_file = (
        tmp_path / "post.molecule"
    )

    session.reaction_metadata = [
        reaction
    ]

    session.ff_files = (
        SimpleNamespace(
            template_files=[],
        )
    )

    ff_cache = builder.cache_dir / "ff.data"
    in_cache = builder.cache_dir / "in.script"

    ff_cache.write_text("ff")
    in_cache.write_text("in")

    monkeypatch.setattr(
        builder,
        "_copy_lunar_files_to_cache",
        lambda files:
        (
            ff_cache,
            in_cache,
        ),
    )

    calls = []

    class FakeDetector:
        def write_wildcard_maps(
            self,
            template_files,
        ):
            calls.append(
                list(
                    template_files
                )
            )

            return [
                "WILDCARD"
            ]

    monkeypatch.setattr(
        builder_module,
        "DeduplicationDetector",
        FakeDetector,
    )

    monkeypatch.setattr(
        builder,
        "_copy_path_to_final",
        lambda src, final:
        Path(src)
        if src is not None
        else None,
    )

    builder.molecule_template_preparation(
        session
    )

    assert calls == [
        [
            reaction
        ]
    ]

    assert (
        session
        .reacter_files
        .template_files
        == [
            "WILDCARD"
        ]
    )


# =============================================================================
# Active-template filtering
# =============================================================================


def test_only_active_complete_reactions_enter_final_template_list(
    tmp_path,
    monkeypatch,
):
    complete = (
        make_reaction_metadata(
            reaction_id=1,
            activity_stats=True,
        )
    )

    inactive = (
        make_reaction_metadata(
            reaction_id=2,
            activity_stats=False,
        )
    )

    incomplete = (
        make_reaction_metadata(
            reaction_id=3,
            activity_stats=True,
        )
    )

    for reaction in (
        complete,
        inactive,
    ):
        reaction.map_file = (
            tmp_path
            / f"RXN_{reaction.reaction_id}.map"
        )

        reaction.pre_reaction_file = (
            tmp_path
            / f"pre{reaction.reaction_id}.molecule"
        )

        reaction.post_reaction_file = (
            tmp_path
            / f"post{reaction.reaction_id}.molecule"
        )

    incomplete.map_file = (
        tmp_path / "RXN_3.map"
    )

    incomplete.pre_reaction_file = (
        tmp_path / "pre3.molecule"
    )

    incomplete.post_reaction_file = None

    builder, session = (
        make_builder(
            tmp_path,
            deduplicate=False,
            reactions=[
                complete,
                inactive,
                incomplete,
            ],
        )
    )

    session.ff_files = (
        SimpleNamespace(
            template_files=[],
        )
    )

    ff_cache = builder.cache_dir / "ff.data"
    in_cache = builder.cache_dir / "in.script"

    ff_cache.write_text("ff")
    in_cache.write_text("in")

    monkeypatch.setattr(
        builder,
        "_copy_lunar_files_to_cache",
        lambda files:
        (
            ff_cache,
            in_cache,
        ),
    )

    monkeypatch.setattr(
        builder,
        "_copy_path_to_final",
        lambda src, final:
        Path(src)
        if src is not None
        else None,
    )

    monkeypatch.setattr(
        builder_module,
        "DeduplicationDetector",
        lambda:
        SimpleNamespace(),
    )

    builder.molecule_template_preparation(
        session
    )

    assert (
        session
        .reacter_files
        .template_files
        == [
            complete
        ]
    )


# =============================================================================
# Final output validation
# =============================================================================


def test_missing_final_force_field_copy_raises(
    tmp_path,
    monkeypatch,
):
    builder, session = (
        make_builder(
            tmp_path
        )
    )

    session.ff_files = (
        SimpleNamespace(
            template_files=[],
        )
    )

    ff_cache = (
        builder.cache_dir
        / "force_field.data"
    )

    in_cache = (
        builder.cache_dir
        / "in.script"
    )

    ff_cache.write_text("ff")
    in_cache.write_text("in")

    monkeypatch.setattr(
        builder,
        "_copy_lunar_files_to_cache",
        lambda files:
        (
            ff_cache,
            in_cache,
        ),
    )

    def fake_copy(
        src,
        final,
    ):
        if Path(src) == ff_cache:
            return None

        return Path(src)

    monkeypatch.setattr(
        builder,
        "_copy_path_to_final",
        fake_copy,
    )

    with pytest.raises(
        FileNotFoundError,
        match="force_field.data was not copied",
    ):
        builder.molecule_template_preparation(
            session
        )


def test_missing_final_input_copy_raises(
    tmp_path,
    monkeypatch,
):
    builder, session = (
        make_builder(
            tmp_path
        )
    )

    session.ff_files = (
        SimpleNamespace(
            template_files=[],
        )
    )

    ff_cache = (
        builder.cache_dir
        / "force_field.data"
    )

    in_cache = (
        builder.cache_dir
        / "in.script"
    )

    ff_cache.write_text("ff")
    in_cache.write_text("in")

    monkeypatch.setattr(
        builder,
        "_copy_lunar_files_to_cache",
        lambda files:
        (
            ff_cache,
            in_cache,
        ),
    )

    def fake_copy(
        src,
        final,
    ):
        if Path(src) == in_cache:
            return None

        return Path(src)

    monkeypatch.setattr(
        builder,
        "_copy_path_to_final",
        fake_copy,
    )

    with pytest.raises(
        FileNotFoundError,
        match="LAMMPS input file was not copied",
    ):
        builder.molecule_template_preparation(
            session
        )


def test_molecule_template_preparation_returns_none_and_sets_session_result(
    tmp_path,
    monkeypatch,
):
    monomer = make_monomer(
        lmp_molecule_file=(
            tmp_path
            / "mma.molecule"
        )
    )

    monomer.lmp_molecule_file.write_text(
        "mol"
    )

    builder, session = (
        make_builder(
            tmp_path,
            monomers=[
                monomer
            ],
            deduplicate=False,
        )
    )

    session.ff_files = (
        SimpleNamespace(
            template_files=[],
        )
    )

    ff_cache = (
        builder.cache_dir
        / "force_field.data"
    )

    in_cache = (
        builder.cache_dir
        / "in.script"
    )

    ff_cache.write_text("ff")
    in_cache.write_text("in")

    monkeypatch.setattr(
        builder,
        "_copy_lunar_files_to_cache",
        lambda files:
        (
            ff_cache,
            in_cache,
        ),
    )

    monkeypatch.setattr(
        builder_module,
        "DeduplicationDetector",
        lambda:
        SimpleNamespace(),
    )

    result = (
        builder
        .molecule_template_preparation(
            session
        )
    )

    assert result is None

    assert isinstance(
        session.reacter_files,
        REACTERFiles,
    )

    assert (
        session
        .reacter_files
        .force_field_data
        == Path(session.output_dir)
        / "force_field.data"
    )

    assert (
        session
        .reacter_files
        .molecule_files
        == [
            monomer
        ]
    )