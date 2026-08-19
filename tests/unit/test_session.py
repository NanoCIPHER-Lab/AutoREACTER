import importlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest


session_module = importlib.import_module("AutoREACTER.session")
cache_module = importlib.import_module("AutoREACTER.cache")

Session = session_module.Session


# =============================================================================
# Helpers
# =============================================================================


def make_validated_inputs(simulation_name="test_simulation"):
    """
    Minimal object representing InputParser.validate_inputs() output.

    Session/read_input only requires simulation_name directly during this
    stage, so a lightweight namespace keeps these tests isolated from
    InputParser itself.
    """
    return SimpleNamespace(
        simulation_name=simulation_name
    )


def write_json(path: Path, data: dict) -> Path:
    path.write_text(
        json.dumps(data),
        encoding="utf-8",
    )
    return path


def configure_read_input_dependencies(
    monkeypatch,
    tmp_path,
    *,
    validated_inputs=None,
    validation_error=None,
):
    """
    Replace Initialization, InputParser, and GetCacheDir so read_input()
    can be tested without touching the real AutoREACTER environment.
    """
    if validated_inputs is None:
        validated_inputs = make_validated_inputs()

    calls = {
        "initialization": 0,
        "parser_instances": 0,
        "validate_calls": 0,
        "validated_data": None,
        "clear_staging": None,
        "cache_calls": 0,
    }

    staging_dir = tmp_path / "staging"

    def fake_initialization():
        calls["initialization"] += 1

    class FakeInputParser:
        def __init__(self):
            calls["parser_instances"] += 1

        def validate_inputs(self, data):
            calls["validate_calls"] += 1
            calls["validated_data"] = data

            if validation_error is not None:
                raise validation_error

            return validated_inputs

    def fake_get_cache_dir(clear_staging=True):
        calls["cache_calls"] += 1
        calls["clear_staging"] = clear_staging

        staging_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        return SimpleNamespace(
            staging_dir=staging_dir
        )

    monkeypatch.setattr(
        session_module,
        "Initialization",
        fake_initialization,
    )

    monkeypatch.setattr(
        session_module,
        "InputParser",
        FakeInputParser,
    )

    monkeypatch.setattr(
        cache_module,
        "GetCacheDir",
        fake_get_cache_dir,
    )

    return calls, staging_dir, validated_inputs


# =============================================================================
# Session dataclass
# =============================================================================


def test_session_stores_required_core_fields(tmp_path):
    inputs = make_validated_inputs()

    session = Session(
        inputs=inputs,
        staging_dir=tmp_path / "staging",
        output_dir=tmp_path / "output",
        images_dir=tmp_path / "output" / "images",
    )

    assert session.inputs is inputs
    assert session.staging_dir == tmp_path / "staging"
    assert session.output_dir == tmp_path / "output"

    assert (
        session.images_dir
        == tmp_path / "output" / "images"
    )


def test_session_pipeline_lists_default_to_empty(tmp_path):
    session = Session(
        inputs=make_validated_inputs(),
        staging_dir=tmp_path / "staging",
        output_dir=tmp_path / "output",
        images_dir=tmp_path / "images",
    )

    assert session.monomer_roles == []
    assert session.reaction_instances == []
    assert session.non_reactants == []
    assert session.reaction_metadata == []


def test_session_mutable_defaults_are_independent(tmp_path):
    session_1 = Session(
        inputs=make_validated_inputs("one"),
        staging_dir=tmp_path / "staging1",
        output_dir=tmp_path / "output1",
        images_dir=tmp_path / "images1",
    )

    session_2 = Session(
        inputs=make_validated_inputs("two"),
        staging_dir=tmp_path / "staging2",
        output_dir=tmp_path / "output2",
        images_dir=tmp_path / "images2",
    )

    session_1.monomer_roles.append("role")
    session_1.reaction_instances.append("reaction")
    session_1.non_reactants.append("nonreactant")
    session_1.reaction_metadata.append("metadata")

    assert session_2.monomer_roles == []
    assert session_2.reaction_instances == []
    assert session_2.non_reactants == []
    assert session_2.reaction_metadata == []


def test_session_file_bundles_default_to_none(tmp_path):
    session = Session(
        inputs=make_validated_inputs(),
        staging_dir=tmp_path / "staging",
        output_dir=tmp_path / "output",
        images_dir=tmp_path / "images",
    )

    assert session.ff_files is None
    assert session.reacter_files is None


def test_session_runtime_defaults(tmp_path):
    session = Session(
        inputs=make_validated_inputs(),
        staging_dir=tmp_path / "staging",
        output_dir=tmp_path / "output",
        images_dir=tmp_path / "images",
    )

    assert session.reaction_id_counter == 0

    assert (
        session.reaction_progression_session
        is None
    )

    assert session.deep_search is True
    assert session.reaction_iteration_depth == 5
    assert session.wildcards is True

    assert (
        session.deduplicate_reaction_templates
        is True
    )

    assert (
        session.write_second_reaction_stage
        is True
    )


def test_session_uses_slots(tmp_path):
    session = Session(
        inputs=make_validated_inputs(),
        staging_dir=tmp_path / "staging",
        output_dir=tmp_path / "output",
        images_dir=tmp_path / "images",
    )

    with pytest.raises(AttributeError):
        session.random_new_attribute = 123


# =============================================================================
# _resolve_input_path
# =============================================================================


def test_resolve_input_path_returns_absolute_json_path(
    tmp_path,
):
    input_file = tmp_path / "input.json"
    input_file.write_text(
        "{}",
        encoding="utf-8",
    )

    result = session_module._resolve_input_path(
        input_file
    )

    assert result == input_file.resolve()
    assert result.is_absolute()


def test_resolve_input_path_accepts_uppercase_json_suffix(
    tmp_path,
):
    input_file = tmp_path / "input.JSON"
    input_file.write_text(
        "{}",
        encoding="utf-8",
    )

    result = session_module._resolve_input_path(
        input_file
    )

    assert result == input_file.resolve()


def test_resolve_input_path_rejects_missing_file(
    tmp_path,
):
    missing = tmp_path / "missing.json"

    with pytest.raises(
        FileNotFoundError,
        match="Input file not found",
    ):
        session_module._resolve_input_path(
            missing
        )


def test_resolve_input_path_rejects_non_json_file(
    tmp_path,
):
    input_file = tmp_path / "input.txt"
    input_file.write_text(
        "{}",
        encoding="utf-8",
    )

    with pytest.raises(
        ValueError,
        match="Input file must be a JSON file",
    ):
        session_module._resolve_input_path(
            input_file
        )


def test_resolve_input_path_rejects_directory_named_json(
    tmp_path,
):
    directory = tmp_path / "fake.json"
    directory.mkdir()

    with pytest.raises(
        ValueError,
        match="Input file must be a JSON file",
    ):
        session_module._resolve_input_path(
            directory
        )


# =============================================================================
# _clear_directory
# =============================================================================


def test_clear_directory_removes_files_but_preserves_root(
    tmp_path,
):
    directory = tmp_path / "output"
    directory.mkdir()

    (directory / "one.txt").write_text("one")
    (directory / "two.txt").write_text("two")

    session_module._clear_directory(
        directory
    )

    assert directory.exists()
    assert directory.is_dir()
    assert list(directory.iterdir()) == []


def test_clear_directory_removes_nested_directories(
    tmp_path,
):
    directory = tmp_path / "output"
    nested = directory / "nested" / "deeper"

    nested.mkdir(parents=True)

    (nested / "data.txt").write_text(
        "content",
        encoding="utf-8",
    )

    session_module._clear_directory(
        directory
    )

    assert directory.exists()
    assert list(directory.iterdir()) == []


def test_clear_directory_nonexistent_path_is_noop(
    tmp_path,
):
    path = tmp_path / "does_not_exist"

    session_module._clear_directory(path)

    assert not path.exists()


def test_clear_directory_file_path_is_noop(
    tmp_path,
):
    path = tmp_path / "file.txt"
    path.write_text(
        "keep me",
        encoding="utf-8",
    )

    session_module._clear_directory(path)

    assert path.exists()

    assert (
        path.read_text(encoding="utf-8")
        == "keep me"
    )


# =============================================================================
# _resolve_output_dir
# =============================================================================


@pytest.mark.parametrize(
    "raw_output_dir",
    [
        None,
        "",
        "   ",
    ],
)
def test_resolve_output_dir_uses_default_when_missing(
    tmp_path,
    raw_output_dir,
):
    input_path = tmp_path / "input.json"

    result = session_module._resolve_output_dir(
        raw_output_dir=raw_output_dir,
        input_path=input_path,
        simulation_name="my_sim",
    )

    expected = (
        tmp_path
        / "AutoREACTER_outputs"
        / "my_sim"
    ).resolve()

    assert result == expected


def test_resolve_output_dir_relative_path_is_relative_to_input(
    tmp_path,
):
    input_path = tmp_path / "input.json"

    result = session_module._resolve_output_dir(
        raw_output_dir="custom/output",
        input_path=input_path,
        simulation_name="ignored_here",
    )

    assert result == (
        tmp_path / "custom" / "output"
    ).resolve()


def test_resolve_output_dir_absolute_path_is_preserved(
    tmp_path,
):
    input_path = tmp_path / "input.json"

    absolute_output = (
        tmp_path
        / "absolute_output"
    ).resolve()

    result = session_module._resolve_output_dir(
        raw_output_dir=str(absolute_output),
        input_path=input_path,
        simulation_name="sim",
    )

    assert result == absolute_output


def test_resolve_output_dir_windows_forward_slash_path():
    input_path = Path("/tmp/input.json")

    result = session_module._resolve_output_dir(
        raw_output_dir=(
            "C:/Users/Janitha/Documents/ARX"
        ),
        input_path=input_path,
        simulation_name="sim",
    )

    assert result == Path(
        "/mnt/c/Users/Janitha/Documents/ARX"
    ).resolve()


def test_resolve_output_dir_windows_backslash_path():
    input_path = Path("/tmp/input.json")

    result = session_module._resolve_output_dir(
        raw_output_dir=(
            r"D:\Projects\AutoREACTER\outputs"
        ),
        input_path=input_path,
        simulation_name="sim",
    )

    assert result == Path(
        "/mnt/d/Projects/AutoREACTER/outputs"
    ).resolve()


def test_resolve_output_dir_windows_drive_is_lowercased():
    result = session_module._resolve_output_dir(
        raw_output_dir=r"E:\Research\run",
        input_path=Path("/tmp/input.json"),
        simulation_name="sim",
    )

    assert result == Path(
        "/mnt/e/Research/run"
    ).resolve()


def test_resolve_output_dir_expands_user_home(
    tmp_path,
    monkeypatch,
):
    fake_home = tmp_path / "home"
    fake_home.mkdir()

    monkeypatch.setenv(
        "HOME",
        str(fake_home),
    )

    result = session_module._resolve_output_dir(
        raw_output_dir="~/arx_output",
        input_path=tmp_path / "input.json",
        simulation_name="sim",
    )

    assert result == (
        fake_home / "arx_output"
    ).resolve()


# =============================================================================
# read_input
# =============================================================================


def test_read_input_calls_initialization_once(
    tmp_path,
    monkeypatch,
):
    calls, _, _ = configure_read_input_dependencies(
        monkeypatch,
        tmp_path,
    )

    input_file = write_json(
        tmp_path / "input.json",
        {
            "simulation_name": "test_simulation",
        },
    )

    session_module.read_input(input_file)

    assert calls["initialization"] == 1


def test_read_input_constructs_input_parser_once(
    tmp_path,
    monkeypatch,
):
    calls, _, _ = configure_read_input_dependencies(
        monkeypatch,
        tmp_path,
    )

    input_file = write_json(
        tmp_path / "input.json",
        {
            "simulation_name": "test_simulation",
        },
    )

    session_module.read_input(input_file)

    assert calls["parser_instances"] == 1
    assert calls["validate_calls"] == 1


def test_read_input_passes_json_data_to_validator(
    tmp_path,
    monkeypatch,
):
    calls, _, _ = configure_read_input_dependencies(
        monkeypatch,
        tmp_path,
    )

    data = {
        "simulation_name": "demo",
        "deep_search": False,
        "wildcards": True,
    }

    input_file = write_json(
        tmp_path / "input.json",
        data,
    )

    session_module.read_input(input_file)

    assert calls["validated_data"] == data


def test_read_input_propagates_clear_staging_true(
    tmp_path,
    monkeypatch,
):
    calls, _, _ = configure_read_input_dependencies(
        monkeypatch,
        tmp_path,
    )

    input_file = write_json(
        tmp_path / "input.json",
        {
            "simulation_name": "demo",
        },
    )

    session_module.read_input(
        input_file,
        clear_staging=True,
    )

    assert calls["clear_staging"] is True


def test_read_input_propagates_clear_staging_false(
    tmp_path,
    monkeypatch,
):
    calls, _, _ = configure_read_input_dependencies(
        monkeypatch,
        tmp_path,
    )

    input_file = write_json(
        tmp_path / "input.json",
        {
            "simulation_name": "demo",
        },
    )

    session_module.read_input(
        input_file,
        clear_staging=False,
    )

    assert calls["clear_staging"] is False


def test_read_input_returns_session(
    tmp_path,
    monkeypatch,
):
    _, staging_dir, validated = (
        configure_read_input_dependencies(
            monkeypatch,
            tmp_path,
        )
    )

    input_file = write_json(
        tmp_path / "input.json",
        {
            "simulation_name": "test_simulation",
        },
    )

    result = session_module.read_input(
        input_file
    )

    assert isinstance(result, Session)
    assert result.inputs is validated
    assert result.staging_dir == staging_dir


def test_read_input_creates_default_output_directory(
    tmp_path,
    monkeypatch,
):
    configure_read_input_dependencies(
        monkeypatch,
        tmp_path,
        validated_inputs=make_validated_inputs(
            "my_simulation"
        ),
    )

    input_file = write_json(
        tmp_path / "input.json",
        {
            "simulation_name": "my_simulation",
        },
    )

    result = session_module.read_input(
        input_file
    )

    expected = (
        tmp_path
        / "AutoREACTER_outputs"
        / "my_simulation"
    ).resolve()

    assert result.output_dir == expected
    assert expected.is_dir()


def test_read_input_creates_images_directory(
    tmp_path,
    monkeypatch,
):
    configure_read_input_dependencies(
        monkeypatch,
        tmp_path,
    )

    input_file = write_json(
        tmp_path / "input.json",
        {
            "simulation_name": "test_simulation",
        },
    )

    result = session_module.read_input(
        input_file
    )

    assert (
        result.images_dir
        == result.output_dir / "images"
    )

    assert result.images_dir.is_dir()


def test_read_input_uses_relative_custom_output_dir(
    tmp_path,
    monkeypatch,
):
    configure_read_input_dependencies(
        monkeypatch,
        tmp_path,
        validated_inputs=make_validated_inputs(
            "demo"
        ),
    )

    input_file = write_json(
        tmp_path / "input.json",
        {
            "simulation_name": "demo",
            "output_dir": "custom_outputs",
        },
    )

    result = session_module.read_input(
        input_file
    )

    assert result.output_dir == (
        tmp_path / "custom_outputs"
    ).resolve()


def test_read_input_uses_absolute_custom_output_dir(
    tmp_path,
    monkeypatch,
):
    configure_read_input_dependencies(
        monkeypatch,
        tmp_path,
        validated_inputs=make_validated_inputs(
            "demo"
        ),
    )

    custom_output = (
        tmp_path / "absolute_custom"
    ).resolve()

    input_file = write_json(
        tmp_path / "input.json",
        {
            "simulation_name": "demo",
            "output_dir": str(custom_output),
        },
    )

    result = session_module.read_input(
        input_file
    )

    assert result.output_dir == custom_output
    assert custom_output.is_dir()


def test_read_input_clears_existing_output_directory(
    tmp_path,
    monkeypatch,
):
    configure_read_input_dependencies(
        monkeypatch,
        tmp_path,
        validated_inputs=make_validated_inputs(
            "demo"
        ),
    )

    output_dir = tmp_path / "existing_output"
    output_dir.mkdir()

    old_file = output_dir / "old.txt"
    old_file.write_text(
        "old data",
        encoding="utf-8",
    )

    nested = output_dir / "old_folder"
    nested.mkdir()

    (nested / "old_nested.txt").write_text(
        "old nested data",
        encoding="utf-8",
    )

    input_file = write_json(
        tmp_path / "input.json",
        {
            "simulation_name": "demo",
            "output_dir": str(output_dir),
        },
    )

    result = session_module.read_input(
        input_file
    )

    assert result.output_dir == output_dir.resolve()

    assert not old_file.exists()
    assert not nested.exists()

    assert (
        output_dir / "images"
    ).is_dir()


def test_read_input_rejects_output_path_that_is_file(
    tmp_path,
    monkeypatch,
):
    configure_read_input_dependencies(
        monkeypatch,
        tmp_path,
        validated_inputs=make_validated_inputs(
            "demo"
        ),
    )

    output_file = tmp_path / "not_a_directory"
    output_file.write_text(
        "existing file",
        encoding="utf-8",
    )

    input_file = write_json(
        tmp_path / "input.json",
        {
            "simulation_name": "demo",
            "output_dir": str(output_file),
        },
    )

    with pytest.raises(
        ValueError,
        match=(
            "Resolved output_dir exists "
            "but is not a directory"
        ),
    ):
        session_module.read_input(
            input_file
        )


@pytest.mark.parametrize(
    "simulation_name",
    [
        "../escape",
        "folder/simulation",
        ".",
        "..",
    ],
)
def test_read_input_rejects_unsafe_simulation_name(
    tmp_path,
    monkeypatch,
    simulation_name,
):
    configure_read_input_dependencies(
        monkeypatch,
        tmp_path,
        validated_inputs=make_validated_inputs(
            simulation_name
        ),
    )

    input_file = write_json(
        tmp_path / "input.json",
        {
            "simulation_name": simulation_name,
        },
    )

    with pytest.raises(
        ValueError,
        match=(
            "Invalid simulation_name "
            "for output directory"
        ),
    ):
        session_module.read_input(
            input_file
        )


def test_read_input_propagates_validation_error(
    tmp_path,
    monkeypatch,
):
    configure_read_input_dependencies(
        monkeypatch,
        tmp_path,
        validation_error=ValueError(
            "bad inputs"
        ),
    )

    input_file = write_json(
        tmp_path / "input.json",
        {
            "simulation_name": "demo",
        },
    )

    with pytest.raises(
        ValueError,
        match="bad inputs",
    ):
        session_module.read_input(
            input_file
        )


def test_read_input_invalid_json_propagates_json_error(
    tmp_path,
    monkeypatch,
):
    configure_read_input_dependencies(
        monkeypatch,
        tmp_path,
    )

    input_file = tmp_path / "input.json"

    input_file.write_text(
        "{ definitely not valid json",
        encoding="utf-8",
    )

    with pytest.raises(
        json.JSONDecodeError
    ):
        session_module.read_input(
            input_file
        )


def test_read_input_missing_file_fails_before_initialization(
    tmp_path,
    monkeypatch,
):
    calls, _, _ = configure_read_input_dependencies(
        monkeypatch,
        tmp_path,
    )

    missing = tmp_path / "missing.json"

    with pytest.raises(FileNotFoundError):
        session_module.read_input(
            missing
        )

    assert calls["initialization"] == 0
    assert calls["cache_calls"] == 0
    assert calls["parser_instances"] == 0


def test_read_input_prints_session_information(
    tmp_path,
    monkeypatch,
    capsys,
):
    configure_read_input_dependencies(
        monkeypatch,
        tmp_path,
        validated_inputs=make_validated_inputs(
            "print_test"
        ),
    )

    input_file = write_json(
        tmp_path / "input.json",
        {
            "simulation_name": "print_test",
        },
    )

    result = session_module.read_input(
        input_file
    )

    output = capsys.readouterr().out

    assert (
        "[INFO] Initialized AutoREACTER Session"
        in output
    )

    assert (
        "[INFO] Simulation Name: print_test"
        in output
    )

    assert str(input_file.resolve()) in output
    assert str(result.staging_dir) in output
    assert str(result.output_dir) in output