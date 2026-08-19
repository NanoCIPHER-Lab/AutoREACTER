from pathlib import Path
from types import SimpleNamespace

import pytest

import AutoREACTER.reaction_preparation.ff_wrapper.lunar_client.locate_lunar as locate_lunar


# =============================================================================
# Helpers
# =============================================================================


def make_valid_lunar_root(
    root: Path,
) -> Path:
    root.mkdir(
        parents=True,
        exist_ok=True,
    )

    for filename in (
        "LUNAR.py",
        "atom_typing.py",
        "all2lmp.py",
        "bond_react_merge.py",
    ):
        (
            root / filename
        ).write_text(
            "# test\n",
            encoding="utf-8",
        )

    (
        root / "src"
    ).mkdir(
        exist_ok=True,
    )

    (
        root / "frc_files"
    ).mkdir(
        exist_ok=True,
    )

    return root


def fake_config(
    tmp_path,
    lunar_root=None,
):
    config_file = (
        tmp_path / "config.py"
    )

    config_file.write_text(
        "LUNAR_ROOT_DIR = None\n",
        encoding="utf-8",
    )

    return SimpleNamespace(
        __file__=str(config_file),
        LUNAR_ROOT_DIR=lunar_root,
    )


# =============================================================================
# _normalize_path
# =============================================================================


def test_normalize_path_none_returns_none():
    assert (
        locate_lunar._normalize_path(
            None
        )
        is None
    )


def test_normalize_path_empty_string_returns_none():
    assert (
        locate_lunar._normalize_path(
            ""
        )
        is None
    )


def test_normalize_path_strips_whitespace(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        locate_lunar.os.path,
        "exists",
        lambda path: (
            False
            if path == "/mnt"
            else Path(path).exists()
        ),
    )

    expected = (
        tmp_path / "folder"
    ).resolve()

    result = (
        locate_lunar._normalize_path(
            f"   {expected}   "
        )
    )

    assert result == str(
        expected
    )


def test_normalize_path_strips_double_quotes(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        locate_lunar.os.path,
        "exists",
        lambda path: (
            False
            if path == "/mnt"
            else Path(path).exists()
        ),
    )

    expected = (
        tmp_path / "folder"
    ).resolve()

    result = (
        locate_lunar._normalize_path(
            f'"{expected}"'
        )
    )

    assert result == str(
        expected
    )


def test_normalize_path_strips_single_quotes(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        locate_lunar.os.path,
        "exists",
        lambda path: (
            False
            if path == "/mnt"
            else Path(path).exists()
        ),
    )

    expected = (
        tmp_path / "folder"
    ).resolve()

    result = (
        locate_lunar._normalize_path(
            f"'{expected}'"
        )
    )

    assert result == str(
        expected
    )


def test_normalize_path_expands_environment_variable(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setenv(
        "AUTOREACTER_LUNAR_TEST",
        str(tmp_path),
    )

    monkeypatch.setattr(
        locate_lunar.os.path,
        "exists",
        lambda path: (
            False
            if path == "/mnt"
            else Path(path).exists()
        ),
    )

    result = (
        locate_lunar._normalize_path(
            "$AUTOREACTER_LUNAR_TEST/lunar"
        )
    )

    assert result == str(
        (
            tmp_path / "lunar"
        ).resolve()
    )


def test_normalize_path_expands_user_home(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setenv(
        "HOME",
        str(tmp_path),
    )

    monkeypatch.setattr(
        locate_lunar.os.path,
        "exists",
        lambda path: (
            False
            if path == "/mnt"
            else Path(path).exists()
        ),
    )

    result = (
        locate_lunar._normalize_path(
            "~/LUNAR"
        )
    )

    assert result == str(
        (
            tmp_path / "LUNAR"
        ).resolve()
    )


def test_normalize_path_converts_windows_path_to_wsl_when_mnt_exists(
    monkeypatch,
):
    monkeypatch.setattr(
        locate_lunar.os.path,
        "exists",
        lambda path: (
            path == "/mnt"
        ),
    )

    result = (
        locate_lunar._normalize_path(
            r"C:\Users\Janitha\LUNAR"
        )
    )

    assert result == (
        "/mnt/c/Users/Janitha/LUNAR"
    )


def test_normalize_path_wsl_drive_is_lowercase(
    monkeypatch,
):
    monkeypatch.setattr(
        locate_lunar.os.path,
        "exists",
        lambda path: (
            path == "/mnt"
        ),
    )

    result = (
        locate_lunar._normalize_path(
            r"D:\Research\LUNAR"
        )
    )

    assert result == (
        "/mnt/d/Research/LUNAR"
    )


def test_normalize_path_does_not_convert_windows_drive_without_mnt(
    tmp_path,
    monkeypatch,
):
    monkeypatch.chdir(
        tmp_path
    )

    monkeypatch.setattr(
        locate_lunar.os.path,
        "exists",
        lambda path: False,
    )

    result = (
        locate_lunar._normalize_path(
            r"C:\Users\Test"
        )
    )

    # Characterizes current non-WSL behavior.
    assert result == str(
        Path(
            r"C:\Users\Test"
        ).resolve()
    )


# =============================================================================
# _is_valid_dir
# =============================================================================


def test_is_valid_dir_none_is_false():
    assert (
        locate_lunar._is_valid_dir(
            None
        )
        is False
    )


def test_is_valid_dir_missing_path_is_false(
    tmp_path,
):
    assert (
        locate_lunar._is_valid_dir(
            tmp_path / "missing"
        )
        is False
    )


def test_is_valid_dir_file_is_false(
    tmp_path,
):
    path = (
        tmp_path / "not_directory"
    )

    path.write_text(
        "test",
        encoding="utf-8",
    )

    assert (
        locate_lunar._is_valid_dir(
            path
        )
        is False
    )


def test_is_valid_dir_accepts_complete_lunar_root(
    tmp_path,
):
    root = make_valid_lunar_root(
        tmp_path / "LUNAR"
    )

    assert (
        locate_lunar._is_valid_dir(
            root
        )
        is True
    )


@pytest.mark.parametrize(
    "missing_file",
    [
        "LUNAR.py",
        "atom_typing.py",
        "all2lmp.py",
        "bond_react_merge.py",
    ],
)
def test_is_valid_dir_requires_each_required_file(
    tmp_path,
    missing_file,
):
    root = make_valid_lunar_root(
        tmp_path / "LUNAR"
    )

    (
        root / missing_file
    ).unlink()

    assert (
        locate_lunar._is_valid_dir(
            root
        )
        is False
    )


@pytest.mark.parametrize(
    "missing_directory",
    [
        "src",
        "frc_files",
    ],
)
def test_is_valid_dir_requires_each_required_directory(
    tmp_path,
    missing_directory,
):
    root = make_valid_lunar_root(
        tmp_path / "LUNAR"
    )

    (
        root / missing_directory
    ).rmdir()

    assert (
        locate_lunar._is_valid_dir(
            root
        )
        is False
    )


# =============================================================================
# _write_config_py
# =============================================================================


def test_write_config_py_writes_string_path(
    tmp_path,
    monkeypatch,
):
    cfg = fake_config(
        tmp_path
    )

    monkeypatch.setattr(
        locate_lunar,
        "config",
        cfg,
    )

    locate_lunar._write_config_py(
        "/tmp/LUNAR"
    )

    text = Path(
        cfg.__file__
    ).read_text(
        encoding="utf-8",
    )

    assert text == (
        "LUNAR_ROOT_DIR = '/tmp/LUNAR'\n"
    )


def test_write_config_py_writes_none(
    tmp_path,
    monkeypatch,
):
    cfg = fake_config(
        tmp_path
    )

    monkeypatch.setattr(
        locate_lunar,
        "config",
        cfg,
    )

    locate_lunar._write_config_py(
        None
    )

    text = Path(
        cfg.__file__
    ).read_text(
        encoding="utf-8",
    )

    assert text == (
        "LUNAR_ROOT_DIR = None\n"
    )


def test_write_config_py_repr_handles_quotes(
    tmp_path,
    monkeypatch,
):
    cfg = fake_config(
        tmp_path
    )

    monkeypatch.setattr(
        locate_lunar,
        "config",
        cfg,
    )

    locate_lunar._write_config_py(
        "/tmp/user's/LUNAR"
    )

    text = Path(
        cfg.__file__
    ).read_text(
        encoding="utf-8",
    )

    namespace = {}

    exec(
        text,
        namespace,
    )

    assert (
        namespace[
            "LUNAR_ROOT_DIR"
        ]
        == "/tmp/user's/LUNAR"
    )


# =============================================================================
# set_LUNAR_loc
# =============================================================================


def test_set_lunar_loc_valid_path(
    tmp_path,
    monkeypatch,
):
    root = make_valid_lunar_root(
        tmp_path / "LUNAR"
    )

    cfg = fake_config(
        tmp_path
    )

    monkeypatch.setattr(
        locate_lunar,
        "config",
        cfg,
    )

    result = (
        locate_lunar.set_LUNAR_loc(
            str(root)
        )
    )

    assert result == str(
        root.resolve()
    )

    assert (
        cfg.LUNAR_ROOT_DIR
        == str(root.resolve())
    )

    text = Path(
        cfg.__file__
    ).read_text(
        encoding="utf-8",
    )

    assert str(
        root.resolve()
    ) in text


def test_set_lunar_loc_invalid_path_raises(
    tmp_path,
    monkeypatch,
):
    cfg = fake_config(
        tmp_path
    )

    monkeypatch.setattr(
        locate_lunar,
        "config",
        cfg,
    )

    with pytest.raises(
        ValueError,
        match="Not a valid directory",
    ):
        locate_lunar.set_LUNAR_loc(
            str(
                tmp_path / "missing"
            )
        )


def test_set_lunar_loc_invalid_path_does_not_change_config(
    tmp_path,
    monkeypatch,
):
    original = "/old/LUNAR"

    cfg = fake_config(
        tmp_path,
        lunar_root=original,
    )

    monkeypatch.setattr(
        locate_lunar,
        "config",
        cfg,
    )

    with pytest.raises(
        ValueError
    ):
        locate_lunar.set_LUNAR_loc(
            str(
                tmp_path / "missing"
            )
        )

    assert (
        cfg.LUNAR_ROOT_DIR
        == original
    )


# =============================================================================
# reset_LUNAR_loc
# =============================================================================


def test_reset_lunar_loc_sets_config_to_none(
    tmp_path,
    monkeypatch,
):
    cfg = fake_config(
        tmp_path,
        lunar_root="/old/LUNAR",
    )

    monkeypatch.setattr(
        locate_lunar,
        "config",
        cfg,
    )

    result = (
        locate_lunar.reset_LUNAR_loc()
    )

    assert result is None

    assert (
        cfg.LUNAR_ROOT_DIR
        is None
    )


def test_reset_lunar_loc_persists_none(
    tmp_path,
    monkeypatch,
):
    cfg = fake_config(
        tmp_path,
        lunar_root="/old/LUNAR",
    )

    monkeypatch.setattr(
        locate_lunar,
        "config",
        cfg,
    )

    locate_lunar.reset_LUNAR_loc()

    assert (
        Path(
            cfg.__file__
        ).read_text(
            encoding="utf-8",
        )
        == "LUNAR_ROOT_DIR = None\n"
    )


# =============================================================================
# _ask_cli
# =============================================================================


def test_ask_cli_returns_normalized_path(
    monkeypatch,
):
    monkeypatch.setattr(
        "builtins.input",
        lambda prompt: (
            "  /tmp/LUNAR  "
        ),
    )

    monkeypatch.setattr(
        locate_lunar,
        "_normalize_path",
        lambda path: (
            f"normalized:{path}"
        ),
    )

    result = (
        locate_lunar._ask_cli()
    )

    assert result == (
        "normalized:/tmp/LUNAR"
    )


def test_ask_cli_empty_input_returns_none(
    monkeypatch,
):
    monkeypatch.setattr(
        "builtins.input",
        lambda prompt: "   ",
    )

    assert (
        locate_lunar._ask_cli()
        is None
    )


# =============================================================================
# _ask_gui
# =============================================================================


def test_ask_gui_returns_normalized_selected_folder(
    monkeypatch,
):
    events = []

    class FakeRoot:
        def withdraw(self):
            events.append(
                "withdraw"
            )

        def destroy(self):
            events.append(
                "destroy"
            )

    monkeypatch.setattr(
        locate_lunar.tk,
        "Tk",
        lambda: FakeRoot(),
    )

    monkeypatch.setattr(
        locate_lunar.filedialog,
        "askdirectory",
        lambda title:
        "/tmp/LUNAR",
    )

    monkeypatch.setattr(
        locate_lunar,
        "_normalize_path",
        lambda path:
        "normalized",
    )

    result = (
        locate_lunar._ask_gui()
    )

    assert result == "normalized"

    assert events == [
        "withdraw",
        "destroy",
    ]


def test_ask_gui_cancel_returns_none(
    monkeypatch,
):
    events = []

    class FakeRoot:
        def withdraw(self):
            events.append(
                "withdraw"
            )

        def destroy(self):
            events.append(
                "destroy"
            )

    monkeypatch.setattr(
        locate_lunar.tk,
        "Tk",
        lambda: FakeRoot(),
    )

    monkeypatch.setattr(
        locate_lunar.filedialog,
        "askdirectory",
        lambda title: "",
    )

    result = (
        locate_lunar._ask_gui()
    )

    assert result is None

    assert events == [
        "withdraw",
        "destroy",
    ]


# =============================================================================
# get_LUNAR_loc - auto detection
# =============================================================================


def test_get_lunar_loc_auto_detects_parent(
    tmp_path,
    monkeypatch,
):
    lunar_root = make_valid_lunar_root(
        tmp_path / "LUNAR"
    )

    module_dir = (
        lunar_root
        / "some"
        / "nested"
    )

    module_dir.mkdir(
        parents=True,
    )

    fake_module_file = (
        module_dir
        / "locate_lunar.py"
    )

    fake_module_file.write_text(
        "",
        encoding="utf-8",
    )

    cfg = fake_config(
        tmp_path
    )

    monkeypatch.setattr(
        locate_lunar,
        "config",
        cfg,
    )

    monkeypatch.setattr(
        locate_lunar,
        "__file__",
        str(fake_module_file),
    )

    result = (
        locate_lunar.get_LUNAR_loc(
            force_prompt=False,
        )
    )

    assert result == str(
        lunar_root
    )

    assert (
        cfg.LUNAR_ROOT_DIR
        == str(lunar_root)
    )


def test_get_lunar_loc_auto_detection_precedes_saved_config(
    tmp_path,
    monkeypatch,
):
    auto_root = make_valid_lunar_root(
        tmp_path
        / "auto"
        / "LUNAR"
    )

    saved_root = make_valid_lunar_root(
        tmp_path
        / "saved"
        / "LUNAR"
    )

    nested = (
        auto_root / "pkg"
    )

    nested.mkdir()

    fake_module = (
        nested / "locate_lunar.py"
    )

    fake_module.write_text(
        "",
        encoding="utf-8",
    )

    cfg = fake_config(
        tmp_path,
        lunar_root=str(saved_root),
    )

    monkeypatch.setattr(
        locate_lunar,
        "config",
        cfg,
    )

    monkeypatch.setattr(
        locate_lunar,
        "__file__",
        str(fake_module),
    )

    result = (
        locate_lunar.get_LUNAR_loc()
    )

    # Characterizes current precedence:
    # auto-detection occurs before saved configuration.
    assert result == str(
        auto_root
    )

    assert (
        cfg.LUNAR_ROOT_DIR
        == str(auto_root)
    )


# =============================================================================
# get_LUNAR_loc - saved configuration
# =============================================================================


def test_get_lunar_loc_returns_saved_valid_path(
    tmp_path,
    monkeypatch,
    capsys,
):
    saved_root = make_valid_lunar_root(
        tmp_path / "saved_lunar"
    )

    cfg = fake_config(
        tmp_path,
        lunar_root=str(saved_root),
    )

    monkeypatch.setattr(
        locate_lunar,
        "config",
        cfg,
    )

    original_validator = (
        locate_lunar._is_valid_dir
    )

    def validator(path):
        return (
            Path(path) == saved_root
            and original_validator(path)
        )

    monkeypatch.setattr(
        locate_lunar,
        "_is_valid_dir",
        validator,
    )

    result = (
        locate_lunar.get_LUNAR_loc(
            force_prompt=False,
        )
    )

    assert result == str(
        saved_root
    )

    assert (
        "Using saved LUNAR root directory"
        in capsys.readouterr().out
    )


def test_get_lunar_loc_force_prompt_ignores_saved_config(
    tmp_path,
    monkeypatch,
):
    cfg = fake_config(
        tmp_path,
        lunar_root="/saved/LUNAR",
    )

    monkeypatch.setattr(
        locate_lunar,
        "config",
        cfg,
    )

    monkeypatch.setattr(
        locate_lunar,
        "_ask_cli",
        lambda: None,
    )

    result = (
        locate_lunar.get_LUNAR_loc(
            force_prompt=True,
            use_gui=False,
        )
    )

    assert result is None


# =============================================================================
# get_LUNAR_loc - CLI prompting
# =============================================================================


def test_get_lunar_loc_cli_accepts_valid_prompted_path(
    tmp_path,
    monkeypatch,
):
    lunar_root = make_valid_lunar_root(
        tmp_path / "LUNAR"
    )

    cfg = fake_config(
        tmp_path
    )

    monkeypatch.setattr(
        locate_lunar,
        "config",
        cfg,
    )

    monkeypatch.setattr(
        locate_lunar,
        "_ask_cli",
        lambda: str(
            lunar_root
        ),
    )

    result = (
        locate_lunar.get_LUNAR_loc(
            force_prompt=True,
            use_gui=False,
        )
    )

    assert result == str(
        lunar_root
    )

    assert (
        cfg.LUNAR_ROOT_DIR
        == str(lunar_root)
    )


def test_get_lunar_loc_cli_cancel_returns_none(
    tmp_path,
    monkeypatch,
):
    cfg = fake_config(
        tmp_path
    )

    monkeypatch.setattr(
        locate_lunar,
        "config",
        cfg,
    )

    monkeypatch.setattr(
        locate_lunar,
        "_ask_cli",
        lambda: None,
    )

    assert (
        locate_lunar.get_LUNAR_loc(
            force_prompt=True,
            use_gui=False,
        )
        is None
    )


def test_get_lunar_loc_cli_retries_after_invalid_path(
    tmp_path,
    monkeypatch,
    capsys,
):
    valid_root = make_valid_lunar_root(
        tmp_path / "valid_lunar"
    )

    invalid_root = (
        tmp_path / "invalid_lunar"
    )

    invalid_root.mkdir()

    cfg = fake_config(
        tmp_path
    )

    monkeypatch.setattr(
        locate_lunar,
        "config",
        cfg,
    )

    answers = iter(
        [
            str(invalid_root),
            str(valid_root),
        ]
    )

    monkeypatch.setattr(
        locate_lunar,
        "_ask_cli",
        lambda: next(
            answers
        ),
    )

    result = (
        locate_lunar.get_LUNAR_loc(
            force_prompt=True,
            use_gui=False,
        )
    )

    assert result == str(
        valid_root
    )

    output = (
        capsys.readouterr().out
    )

    assert (
        "Invalid LUNAR directory"
        in output
    )

    assert (
        "Expected files:"
        in output
    )

    assert (
        "Expected directory:"
        in output
    )


# =============================================================================
# get_LUNAR_loc - GUI prompting
# =============================================================================


def test_get_lunar_loc_gui_accepts_valid_path(
    tmp_path,
    monkeypatch,
):
    lunar_root = make_valid_lunar_root(
        tmp_path / "LUNAR"
    )

    cfg = fake_config(
        tmp_path
    )

    monkeypatch.setattr(
        locate_lunar,
        "config",
        cfg,
    )

    monkeypatch.setattr(
        locate_lunar,
        "_ask_gui",
        lambda: str(
            lunar_root
        ),
    )

    result = (
        locate_lunar.get_LUNAR_loc(
            force_prompt=True,
            use_gui=True,
        )
    )

    assert result == str(
        lunar_root
    )


def test_get_lunar_loc_gui_invalid_path_shows_error_then_cancel(
    tmp_path,
    monkeypatch,
):
    invalid_root = (
        tmp_path / "invalid"
    )

    invalid_root.mkdir()

    cfg = fake_config(
        tmp_path
    )

    monkeypatch.setattr(
        locate_lunar,
        "config",
        cfg,
    )

    answers = iter(
        [
            str(invalid_root),
            None,
        ]
    )

    monkeypatch.setattr(
        locate_lunar,
        "_ask_gui",
        lambda: next(
            answers
        ),
    )

    events = []

    class FakeRoot:
        def withdraw(self):
            events.append(
                "withdraw"
            )

        def destroy(self):
            events.append(
                "destroy"
            )

    monkeypatch.setattr(
        locate_lunar.tk,
        "Tk",
        lambda: FakeRoot(),
    )

    monkeypatch.setattr(
        locate_lunar.messagebox,
        "showerror",
        lambda title, message:
        events.append(
            (
                title,
                message,
            )
        ),
    )

    result = (
        locate_lunar.get_LUNAR_loc(
            force_prompt=True,
            use_gui=True,
        )
    )

    assert result is None

    assert any(
        isinstance(event, tuple)
        and event[0] == "Invalid Folder"
        for event in events
    )


def test_get_lunar_loc_none_use_gui_uses_global_setting(
    tmp_path,
    monkeypatch,
):
    lunar_root = make_valid_lunar_root(
        tmp_path / "LUNAR"
    )

    cfg = fake_config(
        tmp_path
    )

    monkeypatch.setattr(
        locate_lunar,
        "config",
        cfg,
    )

    monkeypatch.setattr(
        locate_lunar,
        "USE_GUI",
        True,
    )

    calls = []

    monkeypatch.setattr(
        locate_lunar,
        "_ask_gui",
        lambda: (
            calls.append(
                "gui"
            )
            or str(lunar_root)
        ),
    )

    monkeypatch.setattr(
        locate_lunar,
        "_ask_cli",
        lambda: pytest.fail(
            "CLI prompt should not be used"
        ),
    )

    result = (
        locate_lunar.get_LUNAR_loc(
            force_prompt=True,
            use_gui=None,
        )
    )

    assert result == str(
        lunar_root
    )

    assert calls == [
        "gui"
    ]