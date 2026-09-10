from pathlib import Path
from types import SimpleNamespace

import pytest

import AutoREACTER.reaction_preparation.ff_wrapper.lunar_client.merge_builder as merge_builder
from AutoREACTER.reaction_preparation.ff_wrapper.lunar_client.merge_builder import (
    write_bond_react_merge_input,
)


# =============================================================================
# Helpers
# =============================================================================


def make_result(
    *,
    result_id,
    molecule,
    filename,
):
    return SimpleNamespace(
        id=result_id,
        molecule=molecule,
        all2lmp_data_file=filename,
    )


# =============================================================================
# Basic file creation
# =============================================================================


def test_write_merge_input_creates_output_directory(
    tmp_path,
):
    output_dir = (
        tmp_path
        / "nested"
        / "bond_react_merge"
    )

    all2lmp_dir = (
        tmp_path / "all2lmp"
    )

    result = (
        write_bond_react_merge_input(
            output_dir,
            all2lmp_dir,
            [],
        )
    )

    assert output_dir.is_dir()

    assert result == (
        output_dir
        / "merge_input.txt"
    )

    assert result.is_file()


def test_write_merge_input_returns_path_object(
    tmp_path,
):
    result = (
        write_bond_react_merge_input(
            tmp_path / "merge",
            tmp_path / "all2lmp",
            [],
        )
    )

    assert isinstance(
        result,
        Path,
    )


def test_write_merge_input_empty_results_still_writes_header(
    tmp_path,
):
    merge_file = (
        write_bond_react_merge_input(
            tmp_path / "merge",
            tmp_path / "all2lmp",
            [],
        )
    )

    text = merge_file.read_text()

    assert (
        '# anything following the "#" character will be ignored'
        in text
    )

    assert "# file-tag" in text

    assert "filename" in text

    assert (
        "comment (required)"
        in text
    )


def test_write_merge_input_writes_parent_directory_comment(
    tmp_path,
):
    merge_file = (
        write_bond_react_merge_input(
            tmp_path / "merge",
            tmp_path / "all2lmp",
            [],
        )
    )

    text = merge_file.read_text()

    assert (
        "Specify the parent_directory"
        in text
    )


# =============================================================================
# Molecule/data entries
# =============================================================================


def test_single_molecule_is_tagged_data1(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        merge_builder,
        "normalize_path",
        lambda path: str(path),
    )

    result = make_result(
        result_id="monomer1",
        molecule=True,
        filename="monomer.data",
    )

    merge_file = (
        write_bond_react_merge_input(
            tmp_path / "merge",
            tmp_path / "all2lmp",
            [
                result,
            ],
        )
    )

    text = merge_file.read_text()

    assert "data1" in text
    assert "monomer.data" in text

    assert (
        "# This datafile will have all coeffs in it"
        in text
    )


def test_multiple_molecules_are_numbered_sequentially(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        merge_builder,
        "normalize_path",
        lambda path: str(path),
    )

    results = [
        make_result(
            result_id="monomer_a",
            molecule=True,
            filename="a.data",
        ),
        make_result(
            result_id="monomer_b",
            molecule=True,
            filename="b.data",
        ),
        make_result(
            result_id="monomer_c",
            molecule=True,
            filename="c.data",
        ),
    ]

    merge_file = (
        write_bond_react_merge_input(
            tmp_path / "merge",
            tmp_path / "all2lmp",
            results,
        )
    )

    text = merge_file.read_text()

    assert "data1" in text
    assert "data2" in text
    assert "data3" in text

    assert (
        text.index("data1")
        < text.index("data2")
        < text.index("data3")
    )


def test_molecule_path_uses_all2lmp_cache_directory(
    tmp_path,
    monkeypatch,
):
    captured = []

    def fake_normalize(path):
        captured.append(
            Path(path)
        )

        return str(path)

    monkeypatch.setattr(
        merge_builder,
        "normalize_path",
        fake_normalize,
    )

    all2lmp_dir = (
        tmp_path / "all2lmp"
    )

    result = make_result(
        result_id="molecule1",
        molecule=True,
        filename="sample.data",
    )

    write_bond_react_merge_input(
        tmp_path / "merge",
        all2lmp_dir,
        [
            result,
        ],
    )

    assert captured == [
        all2lmp_dir
        / "sample.data"
    ]


def test_only_molecule_true_results_become_data_entries(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        merge_builder,
        "normalize_path",
        lambda path: str(path),
    )

    results = [
        make_result(
            result_id="monomer1",
            molecule=True,
            filename="monomer.data",
        ),
        make_result(
            result_id="reaction_without_valid_suffix",
            molecule=False,
            filename="reaction.data",
        ),
    ]

    merge_file = (
        write_bond_react_merge_input(
            tmp_path / "merge",
            tmp_path / "all2lmp",
            results,
        )
    )

    text = merge_file.read_text()

    assert "data1" in text
    assert "monomer.data" in text

    assert "data2" not in text


# =============================================================================
# Reaction pair handling
# =============================================================================


def test_single_pre_post_pair_is_written(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        merge_builder,
        "normalize_path",
        lambda path: str(path),
    )

    results = [
        make_result(
            result_id="pre1",
            molecule=False,
            filename="pre1.data",
        ),
        make_result(
            result_id="post1",
            molecule=False,
            filename="post1.data",
        ),
    ]

    merge_file = (
        write_bond_react_merge_input(
            tmp_path / "merge",
            tmp_path / "all2lmp",
            results,
        )
    )

    text = merge_file.read_text()

    assert "pre1" in text
    assert "post1" in text

    assert "pre1.data" in text
    assert "post1.data" in text

    assert "# for rxn1" in text


def test_reaction_pairs_are_sorted_by_reaction_id(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        merge_builder,
        "normalize_path",
        lambda path: str(path),
    )

    results = [
        make_result(
            result_id="pre10",
            molecule=False,
            filename="pre10.data",
        ),
        make_result(
            result_id="post10",
            molecule=False,
            filename="post10.data",
        ),
        make_result(
            result_id="pre2",
            molecule=False,
            filename="pre2.data",
        ),
        make_result(
            result_id="post2",
            molecule=False,
            filename="post2.data",
        ),
    ]

    merge_file = (
        write_bond_react_merge_input(
            tmp_path / "merge",
            tmp_path / "all2lmp",
            results,
        )
    )

    text = merge_file.read_text()

    # Reaction IDs are sorted numerically:
    # original ID 2 becomes output rxn1,
    # original ID 10 becomes output rxn2.
    assert (
        text.index("pre2.data")
        < text.index("pre10.data")
    )

    lines = text.splitlines()

    pre2_line = next(
        line
        for line in lines
        if "pre2.data" in line
    )

    pre10_line = next(
        line
        for line in lines
        if "pre10.data" in line
    )

    assert pre2_line.startswith(
        "pre1"
    )

    assert pre10_line.startswith(
        "pre2"
    )


def test_reaction_output_numbering_is_contiguous(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        merge_builder,
        "normalize_path",
        lambda path: str(path),
    )

    results = [
        make_result(
            result_id="pre5",
            molecule=False,
            filename="pre5.data",
        ),
        make_result(
            result_id="post5",
            molecule=False,
            filename="post5.data",
        ),
        make_result(
            result_id="pre99",
            molecule=False,
            filename="pre99.data",
        ),
        make_result(
            result_id="post99",
            molecule=False,
            filename="post99.data",
        ),
    ]

    merge_file = (
        write_bond_react_merge_input(
            tmp_path / "merge",
            tmp_path / "all2lmp",
            results,
        )
    )

    lines = merge_file.read_text().splitlines()

    pre_lines = [
        line
        for line in lines
        if line.startswith("pre")
    ]

    post_lines = [
        line
        for line in lines
        if line.startswith("post")
    ]

    assert pre_lines[0].startswith(
        "pre1"
    )

    assert pre_lines[1].startswith(
        "pre2"
    )

    assert post_lines[0].startswith(
        "post1"
    )

    assert post_lines[1].startswith(
        "post2"
    )


def test_pre_post_input_order_does_not_matter(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        merge_builder,
        "normalize_path",
        lambda path: str(path),
    )

    results = [
        make_result(
            result_id="post7",
            molecule=False,
            filename="post7.data",
        ),
        make_result(
            result_id="pre7",
            molecule=False,
            filename="pre7.data",
        ),
    ]

    merge_file = (
        write_bond_react_merge_input(
            tmp_path / "merge",
            tmp_path / "all2lmp",
            results,
        )
    )

    text = merge_file.read_text()

    assert "pre7.data" in text
    assert "post7.data" in text


def test_missing_post_reaction_raises(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        merge_builder,
        "normalize_path",
        lambda path: str(path),
    )

    results = [
        make_result(
            result_id="pre4",
            molecule=False,
            filename="pre4.data",
        )
    ]

    with pytest.raises(
        ValueError,
        match="Incomplete reaction pair for reaction ID 4",
    ):
        write_bond_react_merge_input(
            tmp_path / "merge",
            tmp_path / "all2lmp",
            results,
        )


def test_missing_pre_reaction_raises(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        merge_builder,
        "normalize_path",
        lambda path: str(path),
    )

    results = [
        make_result(
            result_id="post4",
            molecule=False,
            filename="post4.data",
        )
    ]

    with pytest.raises(
        ValueError,
        match="Incomplete reaction pair for reaction ID 4",
    ):
        write_bond_react_merge_input(
            tmp_path / "merge",
            tmp_path / "all2lmp",
            results,
        )


def test_nonmolecule_without_trailing_integer_is_ignored(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        merge_builder,
        "normalize_path",
        lambda path: str(path),
    )

    results = [
        make_result(
            result_id="reaction",
            molecule=False,
            filename="ignored.data",
        )
    ]

    merge_file = (
        write_bond_react_merge_input(
            tmp_path / "merge",
            tmp_path / "all2lmp",
            results,
        )
    )

    text = merge_file.read_text()

    assert "ignored.data" not in text
    assert "pre1" not in text
    assert "post1" not in text


def test_nonmolecule_with_integer_but_unknown_prefix_is_grouped_then_fails(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        merge_builder,
        "normalize_path",
        lambda path: str(path),
    )

    results = [
        make_result(
            result_id="reaction5",
            molecule=False,
            filename="reaction5.data",
        )
    ]

    # Current behavior:
    # reaction5 creates reaction_pairs[5], but does not populate
    # either pre or post, so the pair is incomplete.
    with pytest.raises(
        ValueError,
        match="Incomplete reaction pair for reaction ID 5",
    ):
        write_bond_react_merge_input(
            tmp_path / "merge",
            tmp_path / "all2lmp",
            results,
        )


# =============================================================================
# Duplicate IDs / current overwrite semantics
# =============================================================================


def test_later_pre_entry_with_same_reaction_id_replaces_earlier_one(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        merge_builder,
        "normalize_path",
        lambda path: str(path),
    )

    results = [
        make_result(
            result_id="pre1",
            molecule=False,
            filename="old_pre.data",
        ),
        make_result(
            result_id="pre1",
            molecule=False,
            filename="new_pre.data",
        ),
        make_result(
            result_id="post1",
            molecule=False,
            filename="post.data",
        ),
    ]

    merge_file = (
        write_bond_react_merge_input(
            tmp_path / "merge",
            tmp_path / "all2lmp",
            results,
        )
    )

    text = merge_file.read_text()

    assert "new_pre.data" in text
    assert "old_pre.data" not in text


def test_later_post_entry_with_same_reaction_id_replaces_earlier_one(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        merge_builder,
        "normalize_path",
        lambda path: str(path),
    )

    results = [
        make_result(
            result_id="pre1",
            molecule=False,
            filename="pre.data",
        ),
        make_result(
            result_id="post1",
            molecule=False,
            filename="old_post.data",
        ),
        make_result(
            result_id="post1",
            molecule=False,
            filename="new_post.data",
        ),
    ]

    merge_file = (
        write_bond_react_merge_input(
            tmp_path / "merge",
            tmp_path / "all2lmp",
            results,
        )
    )

    text = merge_file.read_text()

    assert "new_post.data" in text
    assert "old_post.data" not in text


# =============================================================================
# Mixed molecule / reaction output
# =============================================================================


def test_molecules_are_written_before_reactions(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        merge_builder,
        "normalize_path",
        lambda path: str(path),
    )

    results = [
        make_result(
            result_id="pre1",
            molecule=False,
            filename="pre.data",
        ),
        make_result(
            result_id="monomer",
            molecule=True,
            filename="monomer.data",
        ),
        make_result(
            result_id="post1",
            molecule=False,
            filename="post.data",
        ),
    ]

    merge_file = (
        write_bond_react_merge_input(
            tmp_path / "merge",
            tmp_path / "all2lmp",
            results,
        )
    )

    text = merge_file.read_text()

    assert (
        text.index("monomer.data")
        < text.index("pre.data")
    )

    assert (
        text.index("pre.data")
        < text.index("post.data")
    )


def test_data_and_reaction_numbering_are_independent(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        merge_builder,
        "normalize_path",
        lambda path: str(path),
    )

    results = [
        make_result(
            result_id="monomer1",
            molecule=True,
            filename="m1.data",
        ),
        make_result(
            result_id="monomer2",
            molecule=True,
            filename="m2.data",
        ),
        make_result(
            result_id="pre50",
            molecule=False,
            filename="pre50.data",
        ),
        make_result(
            result_id="post50",
            molecule=False,
            filename="post50.data",
        ),
    ]

    merge_file = (
        write_bond_react_merge_input(
            tmp_path / "merge",
            tmp_path / "all2lmp",
            results,
        )
    )

    lines = merge_file.read_text().splitlines()

    assert any(
        line.startswith("data1")
        and "m1.data" in line
        for line in lines
    )

    assert any(
        line.startswith("data2")
        and "m2.data" in line
        for line in lines
    )

    assert any(
        line.startswith("pre1")
        and "pre50.data" in line
        for line in lines
    )

    assert any(
        line.startswith("post1")
        and "post50.data" in line
        for line in lines
    )


# =============================================================================
# normalize_path integration boundary
# =============================================================================


def test_normalize_path_called_for_every_written_input_file(
    tmp_path,
    monkeypatch,
):
    calls = []

    def fake_normalize(path):
        calls.append(
            Path(path)
        )

        return f"NORMALIZED:{path}"

    monkeypatch.setattr(
        merge_builder,
        "normalize_path",
        fake_normalize,
    )

    all2lmp_dir = (
        tmp_path / "all2lmp"
    )

    results = [
        make_result(
            result_id="monomer",
            molecule=True,
            filename="molecule.data",
        ),
        make_result(
            result_id="pre3",
            molecule=False,
            filename="pre.data",
        ),
        make_result(
            result_id="post3",
            molecule=False,
            filename="post.data",
        ),
    ]

    write_bond_react_merge_input(
        tmp_path / "merge",
        all2lmp_dir,
        results,
    )

    assert calls == [
        all2lmp_dir
        / "molecule.data",
        all2lmp_dir
        / "pre.data",
        all2lmp_dir
        / "post.data",
    ]


def test_normalized_paths_are_used_in_output(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        merge_builder,
        "normalize_path",
        lambda path: (
            f"/normalized/{Path(path).name}"
        ),
    )

    results = [
        make_result(
            result_id="monomer",
            molecule=True,
            filename="molecule.data",
        ),
        make_result(
            result_id="pre1",
            molecule=False,
            filename="pre.data",
        ),
        make_result(
            result_id="post1",
            molecule=False,
            filename="post.data",
        ),
    ]

    merge_file = (
        write_bond_react_merge_input(
            tmp_path / "merge",
            tmp_path / "all2lmp",
            results,
        )
    )

    text = merge_file.read_text()

    assert (
        "/normalized/molecule.data"
        in text
    )

    assert (
        "/normalized/pre.data"
        in text
    )

    assert (
        "/normalized/post.data"
        in text
    )


# =============================================================================
# get_ending_integer integration boundary
# =============================================================================


def test_get_ending_integer_called_for_nonmolecule_results_only(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        merge_builder,
        "normalize_path",
        lambda path: str(path),
    )

    calls = []

    def fake_get_ending_integer(value):
        calls.append(value)

        if value == "preX":
            return 1

        if value == "postX":
            return 1

        return None

    monkeypatch.setattr(
        merge_builder,
        "get_ending_integer",
        fake_get_ending_integer,
    )

    results = [
        make_result(
            result_id="moleculeX",
            molecule=True,
            filename="molecule.data",
        ),
        make_result(
            result_id="preX",
            molecule=False,
            filename="pre.data",
        ),
        make_result(
            result_id="postX",
            molecule=False,
            filename="post.data",
        ),
    ]

    write_bond_react_merge_input(
        tmp_path / "merge",
        tmp_path / "all2lmp",
        results,
    )

    assert calls == [
        "preX",
        "postX",
    ]


# =============================================================================
# Input immutability
# =============================================================================


def test_function_does_not_mutate_results_list(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        merge_builder,
        "normalize_path",
        lambda path: str(path),
    )

    first = make_result(
        result_id="pre1",
        molecule=False,
        filename="pre.data",
    )

    second = make_result(
        result_id="post1",
        molecule=False,
        filename="post.data",
    )

    results = [
        first,
        second,
    ]

    original = list(
        results
    )

    write_bond_react_merge_input(
        tmp_path / "merge",
        tmp_path / "all2lmp",
        results,
    )

    assert results == original

    assert results[0] is first
    assert results[1] is second


def test_function_does_not_modify_result_objects(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        merge_builder,
        "normalize_path",
        lambda path: str(path),
    )

    pre = make_result(
        result_id="pre1",
        molecule=False,
        filename="pre.data",
    )

    post = make_result(
        result_id="post1",
        molecule=False,
        filename="post.data",
    )

    write_bond_react_merge_input(
        tmp_path / "merge",
        tmp_path / "all2lmp",
        [
            pre,
            post,
        ],
    )

    assert pre.id == "pre1"
    assert pre.molecule is False
    assert (
        pre.all2lmp_data_file
        == "pre.data"
    )

    assert post.id == "post1"
    assert post.molecule is False
    assert (
        post.all2lmp_data_file
        == "post.data"
    )