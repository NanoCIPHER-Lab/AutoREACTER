from pathlib import Path

import pytest

import AutoREACTER.reaction_preparation.ff_wrapper.lunar_client.lunar_utils as lunar_utils
from AutoREACTER.reaction_preparation.ff_wrapper.lunar_client.lunar_utils import (
    get_ending_integer,
    is_wsl,
    loading_screen,
    move_merge_outputs,
    normalize_path,
)


# =============================================================================
# is_wsl
# =============================================================================


def test_is_wsl_true_when_platform_release_contains_microsoft(
    monkeypatch,
):
    monkeypatch.setattr(
        lunar_utils.platform,
        "release",
        lambda: (
            "5.15.153.1-microsoft-standard-WSL2"
        ),
    )

    monkeypatch.delenv(
        "WSL_INTEROP",
        raising=False,
    )

    assert is_wsl() is True


def test_is_wsl_release_check_is_case_insensitive(
    monkeypatch,
):
    monkeypatch.setattr(
        lunar_utils.platform,
        "release",
        lambda: (
            "5.15-MICROSOFT-STANDARD"
        ),
    )

    monkeypatch.delenv(
        "WSL_INTEROP",
        raising=False,
    )

    assert is_wsl() is True


def test_is_wsl_true_when_wsl_interop_environment_exists(
    monkeypatch,
):
    monkeypatch.setattr(
        lunar_utils.platform,
        "release",
        lambda: "Linux",
    )

    monkeypatch.setenv(
        "WSL_INTEROP",
        "/run/WSL/123_interop",
    )

    assert is_wsl() is True


def test_is_wsl_false_on_normal_linux(
    monkeypatch,
):
    monkeypatch.setattr(
        lunar_utils.platform,
        "release",
        lambda: (
            "6.8.0-generic"
        ),
    )

    monkeypatch.delenv(
        "WSL_INTEROP",
        raising=False,
    )

    assert is_wsl() is False


# =============================================================================
# normalize_path
# =============================================================================


def test_normalize_path_accepts_path_object(
    monkeypatch,
):
    monkeypatch.setattr(
        lunar_utils,
        "is_wsl",
        lambda: False,
    )

    monkeypatch.setattr(
        lunar_utils.platform,
        "system",
        lambda: "Linux",
    )

    result = normalize_path(
        Path("/tmp/test/path")
    )

    assert result == (
        "/tmp/test/path"
    )


def test_normalize_path_strips_whitespace_and_double_quotes(
    monkeypatch,
):
    monkeypatch.setattr(
        lunar_utils,
        "is_wsl",
        lambda: False,
    )

    monkeypatch.setattr(
        lunar_utils.platform,
        "system",
        lambda: "Linux",
    )

    result = normalize_path(
        '  "/tmp/test/path"  '
    )

    assert result == (
        "/tmp/test/path"
    )


def test_normalize_path_strips_single_quotes(
    monkeypatch,
):
    monkeypatch.setattr(
        lunar_utils,
        "is_wsl",
        lambda: False,
    )

    monkeypatch.setattr(
        lunar_utils.platform,
        "system",
        lambda: "Linux",
    )

    result = normalize_path(
        "'/tmp/test/path'"
    )

    assert result == (
        "/tmp/test/path"
    )


def test_normalize_path_converts_backslashes_on_linux(
    monkeypatch,
):
    monkeypatch.setattr(
        lunar_utils,
        "is_wsl",
        lambda: False,
    )

    monkeypatch.setattr(
        lunar_utils.platform,
        "system",
        lambda: "Linux",
    )

    result = normalize_path(
        r"/tmp/test\folder\file"
    )

    assert result == (
        "/tmp/test/folder/file"
    )


def test_normalize_path_windows_drive_to_wsl(
    monkeypatch,
):
    monkeypatch.setattr(
        lunar_utils,
        "is_wsl",
        lambda: True,
    )

    result = normalize_path(
        "C:/Users/Janitha/LUNAR"
    )

    assert result == (
        "/mnt/c/Users/Janitha/LUNAR"
    )


def test_normalize_path_windows_backslashes_to_wsl(
    monkeypatch,
):
    monkeypatch.setattr(
        lunar_utils,
        "is_wsl",
        lambda: True,
    )

    result = normalize_path(
        r"D:\Research\LUNAR\data"
    )

    assert result == (
        "/mnt/d/Research/LUNAR/data"
    )


def test_normalize_path_wsl_drive_letter_is_lowercase(
    monkeypatch,
):
    monkeypatch.setattr(
        lunar_utils,
        "is_wsl",
        lambda: True,
    )

    result = normalize_path(
        "Z:/Some/Folder"
    )

    assert result == (
        "/mnt/z/Some/Folder"
    )


def test_normalize_path_existing_wsl_path_is_preserved(
    monkeypatch,
):
    monkeypatch.setattr(
        lunar_utils,
        "is_wsl",
        lambda: True,
    )

    result = normalize_path(
        "/mnt/c/Users/Janitha/LUNAR"
    )

    assert result == (
        "/mnt/c/Users/Janitha/LUNAR"
    )


def test_normalize_path_wsl_path_to_windows(
    monkeypatch,
):
    monkeypatch.setattr(
        lunar_utils,
        "is_wsl",
        lambda: False,
    )

    monkeypatch.setattr(
        lunar_utils.platform,
        "system",
        lambda: "Windows",
    )

    result = normalize_path(
        "/mnt/c/Users/Janitha/LUNAR"
    )

    assert result == (
        r"C:\Users\Janitha\LUNAR"
    )


def test_normalize_path_wsl_path_without_leading_slash_to_windows(
    monkeypatch,
):
    monkeypatch.setattr(
        lunar_utils,
        "is_wsl",
        lambda: False,
    )

    monkeypatch.setattr(
        lunar_utils.platform,
        "system",
        lambda: "Windows",
    )

    result = normalize_path(
        "mnt/d/Research/LUNAR"
    )

    assert result == (
        r"D:\Research\LUNAR"
    )


def test_normalize_path_linux_collapses_parent_segments(
    monkeypatch,
):
    monkeypatch.setattr(
        lunar_utils,
        "is_wsl",
        lambda: False,
    )

    monkeypatch.setattr(
        lunar_utils.platform,
        "system",
        lambda: "Linux",
    )

    result = normalize_path(
        "/tmp/a/../b"
    )

    assert result == (
        "/tmp/b"
    )


# =============================================================================
# get_ending_integer
# =============================================================================


@pytest.mark.parametrize(
    "value, expected",
    [
        (
            "pre12",
            12,
        ),
        (
            "reaction_123",
            123,
        ),
        (
            "RXN_1",
            1,
        ),
        (
            "abc0007",
            7,
        ),
        (
            "42",
            42,
        ),
        (
            "test0",
            0,
        ),
    ],
)
def test_get_ending_integer(
    value,
    expected,
):
    assert (
        get_ending_integer(
            value
        )
        == expected
    )


@pytest.mark.parametrize(
    "value, expected",
    [
        ("pre12", 12),
        ("reaction_123", 123),
        ("RXN_1", 1),
        ("abc0007", 7),
        ("42", 42),
        ("test0", 0),
        ("1.5", 5),
    ],
)
def test_get_ending_integer(
    value,
    expected,
):
    assert (
        get_ending_integer(value)
        == expected
    )


@pytest.mark.parametrize(
    "value",
    [
        "",
        "abc",
        "123abc",
        "rxn_1_test",
    ],
)
def test_get_ending_integer_none_when_no_trailing_integer(
    value,
):
    assert (
        get_ending_integer(value)
        is None
    )

# =============================================================================
# move_merge_outputs
# =============================================================================


def test_move_merge_outputs_creates_destination(
    tmp_path,
):
    src = tmp_path / "src"
    dst = tmp_path / "dst"

    src.mkdir()

    move_merge_outputs(
        src,
        dst,
    )

    assert dst.is_dir()


def test_move_merge_outputs_moves_expected_files(
    tmp_path,
):
    src = tmp_path / "src"
    dst = tmp_path / "dst"

    src.mkdir()

    filenames = [
        "system_merged.data",
        "system_merged.lmpmol",
        "force_field.data",
        "log.lammps",
        "all2lmp.log",
        "other.log",
    ]

    for filename in filenames:
        (
            src / filename
        ).write_text(
            filename,
            encoding="utf-8",
        )

    move_merge_outputs(
        src,
        dst,
    )

    for filename in filenames:
        assert (
            dst / filename
        ).is_file()

        assert not (
            src / filename
        ).exists()


def test_move_merge_outputs_leaves_irrelevant_files_in_source(
    tmp_path,
):
    src = tmp_path / "src"
    dst = tmp_path / "dst"

    src.mkdir()

    irrelevant = [
        "input.data",
        "notes.txt",
        "template.molecule",
        "random.json",
    ]

    for filename in irrelevant:
        (
            src / filename
        ).write_text(
            "keep",
            encoding="utf-8",
        )

    move_merge_outputs(
        src,
        dst,
    )

    for filename in irrelevant:
        assert (
            src / filename
        ).is_file()

        assert not (
            dst / filename
        ).exists()


def test_move_merge_outputs_preserves_file_contents(
    tmp_path,
):
    src = tmp_path / "src"
    dst = tmp_path / "dst"

    src.mkdir()

    source_file = (
        src / "sample_merged.data"
    )

    source_file.write_text(
        "LAMMPS DATA CONTENT",
        encoding="utf-8",
    )

    move_merge_outputs(
        src,
        dst,
    )

    assert (
        dst
        / "sample_merged.data"
    ).read_text(
        encoding="utf-8"
    ) == "LAMMPS DATA CONTENT"


def test_move_merge_outputs_accepts_string_paths(
    tmp_path,
):
    src = tmp_path / "src"
    dst = tmp_path / "dst"

    src.mkdir()

    (
        src / "force_field.data"
    ).write_text(
        "force field",
        encoding="utf-8",
    )

    move_merge_outputs(
        str(src),
        str(dst),
    )

    assert (
        dst / "force_field.data"
    ).is_file()


def test_move_merge_outputs_empty_source_is_noop_except_destination_creation(
    tmp_path,
):
    src = tmp_path / "src"
    dst = tmp_path / "dst"

    src.mkdir()

    move_merge_outputs(
        src,
        dst,
    )

    assert dst.is_dir()

    assert list(
        dst.iterdir()
    ) == []


# =============================================================================
# loading_screen
# =============================================================================


def test_loading_screen_prints_banner_and_ready(
    monkeypatch,
    capsys,
):
    monkeypatch.setattr(
        lunar_utils.time,
        "sleep",
        lambda seconds: None,
    )

    loading_screen()

    output = capsys.readouterr().out

    assert "Loading LUNAR" in output
    assert "Ready!" in output

    # Stable portion of the ASCII banner.
    assert "█████" in output


def test_loading_screen_uses_custom_name(
    monkeypatch,
    capsys,
):
    monkeypatch.setattr(
        lunar_utils.time,
        "sleep",
        lambda seconds: None,
    )

    loading_screen(
        "all2lmp"
    )

    output = capsys.readouterr().out

    assert (
        "Loading all2lmp"
        in output
    )


def test_loading_screen_runs_ten_spinner_steps(
    monkeypatch,
):
    sleep_calls = []

    monkeypatch.setattr(
        lunar_utils.time,
        "sleep",
        lambda seconds:
        sleep_calls.append(
            seconds
        ),
    )

    loading_screen(
        "test"
    )

    assert sleep_calls == [
        0.1,
        0.1,
        0.1,
        0.1,
        0.1,
        0.1,
        0.1,
        0.1,
        0.1,
        0.1,
    ]


def test_loading_screen_returns_none(
    monkeypatch,
):
    monkeypatch.setattr(
        lunar_utils.time,
        "sleep",
        lambda seconds: None,
    )

    result = loading_screen(
        "test"
    )

    assert result is None