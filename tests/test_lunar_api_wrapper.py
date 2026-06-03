"""
Lightweight tests for lunar_api_wrapper.py.

Store this file at:
    tests/test_lunar_api_wrapper.py

These tests intentionally avoid running LUNAR itself. They test the wrapper's
path handling and merge_input.txt generation logic.
"""

from pathlib import Path

import pytest

from AutoREACTER.reaction_preparation.lunar_client.lunar_api_wrapper import (
    All2LMPResult,
    LunarAPIWrapper,
)


def make_wrapper_without_init(tmp_path: Path) -> LunarAPIWrapper:
    """Create a LunarAPIWrapper instance without calling __init__."""

    wrapper = LunarAPIWrapper.__new__(LunarAPIWrapper)
    wrapper.cache_dir = tmp_path / "lunar"
    wrapper.cache_atom_typing = wrapper.cache_dir / "atom_typing"
    wrapper.cache_all2lmp = wrapper.cache_dir / "all2lmp"
    wrapper.cache_bond_react_merge = wrapper.cache_dir / "bond_react_merge"
    wrapper.cache_mergeprep = wrapper.cache_dir / "all2lmp_2_bondreact_mergeprep"
    wrapper.LUNAR_LOCATION = tmp_path / "LUNAR"
    return wrapper


def test_get_ending_integer(tmp_path):
    wrapper = make_wrapper_without_init(tmp_path)

    assert wrapper._get_ending_integer("pre12") == 12
    assert wrapper._get_ending_integer("post3") == 3
    assert wrapper._get_ending_integer("data") is None


def test_normalize_windows_path_inside_wsl(tmp_path, monkeypatch):
    wrapper = make_wrapper_without_init(tmp_path)

    monkeypatch.setattr(wrapper, "_is_wsl", lambda: True)

    normalized = wrapper._normalize_path(r"C:\Users\janit\Documents\file.data")

    assert normalized == "/mnt/c/Users/janit/Documents/file.data"


def test_write_bond_react_merge_input(tmp_path, monkeypatch):
    wrapper = make_wrapper_without_init(tmp_path)

    wrapper.cache_all2lmp.mkdir(parents=True)
    wrapper.cache_bond_react_merge.mkdir(parents=True)

    # Keep paths deterministic for this unit test.
    monkeypatch.setattr(wrapper, "_normalize_path", lambda p: str(p))

    monomer = wrapper.cache_all2lmp / "BisphenolA_typed_IFF.data"
    pre = wrapper.cache_all2lmp / "pre1_typed_IFF.data"
    post = wrapper.cache_all2lmp / "post1_typed_IFF.data"

    results = [
        All2LMPResult(id="BisphenolA", molecule=True, all2lmp_data_file=monomer),
        All2LMPResult(id="pre1", molecule=False, all2lmp_data_file=pre),
        All2LMPResult(id="post1", molecule=False, all2lmp_data_file=post),
    ]

    merge_file = wrapper._write_bond_react_merge_input(results)
    text = merge_file.read_text(encoding="utf-8")

    assert merge_file.name == "merge_input.txt"
    assert "data1" in text
    assert "pre1" in text
    assert "post1" in text
    assert str(monomer) in text
    assert str(pre) in text
    assert str(post) in text


def test_write_bond_react_merge_input_rejects_missing_post(tmp_path):
    wrapper = make_wrapper_without_init(tmp_path)
    wrapper.cache_all2lmp.mkdir(parents=True)
    wrapper.cache_bond_react_merge.mkdir(parents=True)

    results = [
        All2LMPResult(
            id="pre1",
            molecule=False,
            all2lmp_data_file=wrapper.cache_all2lmp / "pre1_typed_IFF.data",
        ),
    ]

    with pytest.raises(ValueError, match="Incomplete reaction pair"):
        wrapper._write_bond_react_merge_input(results)
