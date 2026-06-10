"""
Lightweight tests for lunar_api_wrapper.py.

Store this file at:
    tests/test_lunar_api_wrapper.py

These tests intentionally avoid running LUNAR itself. They test the wrapper's
path handling and merge_input.txt generation logic.
"""

from pathlib import Path

import pytest

from AutoREACTER.reaction_preparation.ff_wrapper.lunar_client.lunar_executor import All2LMPResult
from AutoREACTER.reaction_preparation.ff_wrapper.lunar_client.lunar_utils import (
    get_ending_integer,
    normalize_path,
)
from AutoREACTER.reaction_preparation.ff_wrapper.lunar_client.merge_builder import write_bond_react_merge_input


def test_get_ending_integer(tmp_path):
    assert get_ending_integer("pre12") == 12
    assert get_ending_integer("post3") == 3
    assert get_ending_integer("data") is None


def test_normalize_windows_path_inside_wsl(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "AutoREACTER.reaction_preparation.ff_wrapper.lunar_client.lunar_utils.is_wsl",
        lambda: True,
    )
    normalized = normalize_path(r"C:\Users\janit\Documents\file.data")

    assert normalized == "/mnt/c/Users/janit/Documents/file.data"


def test_write_bond_react_merge_input(tmp_path, monkeypatch):
    cache_all2lmp = tmp_path / "lunar" / "all2lmp"
    cache_bond_react_merge = tmp_path / "lunar" / "bond_react_merge"
    cache_all2lmp.mkdir(parents=True)
    cache_bond_react_merge.mkdir(parents=True)

    # Keep paths deterministic for this unit test.
    monkeypatch.setattr(
        "AutoREACTER.reaction_preparation.ff_wrapper.lunar_client.merge_builder.normalize_path",
        lambda p: str(p),
    )

    monomer = Path("BisphenolA_typed_IFF.data")
    pre = Path("pre1_typed_IFF.data")
    post = Path("post1_typed_IFF.data")

    results = [
        All2LMPResult(id="BisphenolA", molecule=True, all2lmp_data_file=monomer),
        All2LMPResult(id="pre1", molecule=False, all2lmp_data_file=pre),
        All2LMPResult(id="post1", molecule=False, all2lmp_data_file=post),
    ]

    merge_file = write_bond_react_merge_input(
        cache_bond_react_merge=cache_bond_react_merge,
        cache_all2lmp=cache_all2lmp,
        all2lmp_results=results,
    )
    text = merge_file.read_text(encoding="utf-8")

    assert merge_file.name == "merge_input.txt"
    assert "data1" in text
    assert "pre1" in text
    assert "post1" in text
    assert str(cache_all2lmp / monomer) in text
    assert str(cache_all2lmp / pre) in text
    assert str(cache_all2lmp / post) in text


def test_write_bond_react_merge_input_rejects_missing_post(tmp_path):
    cache_all2lmp = tmp_path / "lunar" / "all2lmp"
    cache_bond_react_merge = tmp_path / "lunar" / "bond_react_merge"
    cache_all2lmp.mkdir(parents=True)
    cache_bond_react_merge.mkdir(parents=True)

    results = [
        All2LMPResult(
            id="pre1",
            molecule=False,
            all2lmp_data_file=Path("pre1_typed_IFF.data"),
        ),
    ]

    with pytest.raises(ValueError, match="Incomplete reaction pair"):
        write_bond_react_merge_input(
            cache_bond_react_merge=cache_bond_react_merge,
            cache_all2lmp=cache_all2lmp,
            all2lmp_results=results,
        )
