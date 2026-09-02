from pathlib import Path
from types import SimpleNamespace

import pytest

import AutoREACTER as arx


# =============================================================================
# Fixtures / helpers
# =============================================================================


@pytest.fixture(autouse=True)
def reset_active_workflow(monkeypatch):
    """
    Every test starts without an active global workflow.

    AutoREACTER intentionally keeps one package-level active workflow, so
    isolating this state prevents tests from affecting one another.
    """
    monkeypatch.setattr(
        arx,
        "_active_workflow",
        None,
    )


class FakeWorkflow:
    """
    Lightweight ARXCLI-like object used to verify public API delegation.
    """

    def __init__(self):
        self.session = object()
        self.calls = []

    def show_molecules(self):
        self.calls.append(
            ("show_molecules",)
        )
        return "molecules-image"

    def show_functional_groups(self):
        self.calls.append(
            ("show_functional_groups",)
        )
        return "functional-groups-image"

    def show_reactions(self):
        self.calls.append(
            ("show_reactions",)
        )
        return "reactions-image"

    def select_reactions(self):
        self.calls.append(
            ("select_reactions",)
        )
        return None

    def show_non_reactants(self):
        self.calls.append(
            ("show_non_reactants",)
        )
        return "non-reactants-image"

    def select_non_reactants(self):
        self.calls.append(
            ("select_non_reactants",)
        )
        return None

    def prepare_reactions(self):
        self.calls.append(
            ("prepare_reactions",)
        )
        return None

    def show_reaction_templates(
        self,
        highlight_type="template",
    ):
        self.calls.append(
            (
                "show_reaction_templates",
                highlight_type,
            )
        )

        return (
            f"templates-{highlight_type}"
        )

    def process(self):
        self.calls.append(
            ("process",)
        )
        return None


def install_fake_workflow(
    monkeypatch,
):
    workflow = FakeWorkflow()

    monkeypatch.setattr(
        arx,
        "_active_workflow",
        workflow,
    )

    return workflow


# =============================================================================
# Package metadata
# =============================================================================


def test_package_title():
    assert arx.__title__ == "AutoREACTER"


def test_package_version():
    assert arx.__version__ == "0.3"


def test_package_release_matches_version():
    assert (
        arx.__release__
        == arx.__version__
    )


def test_package_license():
    assert arx.__license__ == "MIT"


def test_package_authors():
    assert arx.__authors__ == [
        "Janitha Mahanthe",
        "Jacob Gissinger",
    ]


def test_package_author_string_matches_authors():
    assert (
        arx.__author__
        == ", ".join(arx.__authors__)
    )


# =============================================================================
# _ensure_workflow
# =============================================================================


def test_ensure_workflow_raises_without_active_session():
    with pytest.raises(
        RuntimeError,
        match="No active session",
    ):
        arx._ensure_workflow()


def test_ensure_workflow_error_tells_user_to_call_run():
    with pytest.raises(RuntimeError) as exc_info:
        arx._ensure_workflow()

    message = str(exc_info.value)

    assert "arx.run" in message
    assert "your_file.json" in message


def test_ensure_workflow_returns_active_workflow(
    monkeypatch,
):
    workflow = object()

    monkeypatch.setattr(
        arx,
        "_active_workflow",
        workflow,
    )

    assert (
        arx._ensure_workflow()
        is workflow
    )


# =============================================================================
# session()
# =============================================================================


def test_session_requires_active_workflow():
    with pytest.raises(
        RuntimeError,
        match="No active session",
    ):
        arx.session()


def test_session_returns_active_workflow_session(
    monkeypatch,
):
    workflow = FakeWorkflow()

    monkeypatch.setattr(
        arx,
        "_active_workflow",
        workflow,
    )

    assert (
        arx.session()
        is workflow.session
    )


# =============================================================================
# run()
# =============================================================================


def test_run_rejects_missing_input_file(
    tmp_path,
    monkeypatch,
):
    constructed = []

    class FakeARXCLI:
        def __init__(
            self,
            input_path,
        ):
            constructed.append(
                input_path
            )

    monkeypatch.setattr(
        arx,
        "ARXCLI",
        FakeARXCLI,
    )

    missing = (
        tmp_path
        / "missing.json"
    )

    with pytest.raises(
        FileNotFoundError,
        match="Input file not found",
    ):
        arx.run(missing)

    assert constructed == []

    assert (
        arx._active_workflow
        is None
    )


def test_run_constructs_arxcli_with_resolved_path(
    tmp_path,
    monkeypatch,
):
    input_file = (
        tmp_path / "input.json"
    )

    input_file.write_text(
        "{}",
        encoding="utf-8",
    )

    received = []

    class FakeARXCLI:
        def __init__(
            self,
            input_path,
        ):
            received.append(
                input_path
            )

    monkeypatch.setattr(
        arx,
        "ARXCLI",
        FakeARXCLI,
    )

    result = arx.run(
        input_file
    )

    assert received == [
        input_file.resolve()
    ]

    assert (
        result
        is arx._active_workflow
    )


def test_run_accepts_string_path(
    tmp_path,
    monkeypatch,
):
    input_file = (
        tmp_path / "input.json"
    )

    input_file.write_text(
        "{}",
        encoding="utf-8",
    )

    received = []

    class FakeARXCLI:
        def __init__(
            self,
            input_path,
        ):
            received.append(
                input_path
            )

    monkeypatch.setattr(
        arx,
        "ARXCLI",
        FakeARXCLI,
    )

    arx.run(
        str(input_file)
    )

    assert received == [
        input_file.resolve()
    ]

    assert isinstance(
        received[0],
        Path,
    )


def test_run_resolves_relative_path(
    tmp_path,
    monkeypatch,
):
    input_file = (
        tmp_path / "input.json"
    )

    input_file.write_text(
        "{}",
        encoding="utf-8",
    )

    monkeypatch.chdir(
        tmp_path
    )

    received = []

    class FakeARXCLI:
        def __init__(
            self,
            input_path,
        ):
            received.append(
                input_path
            )

    monkeypatch.setattr(
        arx,
        "ARXCLI",
        FakeARXCLI,
    )

    arx.run(
        "input.json"
    )

    assert received == [
        input_file.resolve()
    ]


def test_run_returns_created_workflow(
    tmp_path,
    monkeypatch,
):
    input_file = (
        tmp_path / "input.json"
    )

    input_file.write_text(
        "{}",
        encoding="utf-8",
    )

    created = object()

    monkeypatch.setattr(
        arx,
        "ARXCLI",
        lambda path: created,
    )

    result = arx.run(
        input_file
    )

    assert result is created

    assert (
        arx._active_workflow
        is created
    )


def test_run_replaces_existing_workflow(
    tmp_path,
    monkeypatch,
):
    first_input = (
        tmp_path / "first.json"
    )

    second_input = (
        tmp_path / "second.json"
    )

    first_input.write_text(
        "{}",
        encoding="utf-8",
    )

    second_input.write_text(
        "{}",
        encoding="utf-8",
    )

    created = []

    class FakeARXCLI:
        def __init__(
            self,
            input_path,
        ):
            self.input_path = (
                input_path
            )

            created.append(self)

    monkeypatch.setattr(
        arx,
        "ARXCLI",
        FakeARXCLI,
    )

    first = arx.run(
        first_input
    )

    second = arx.run(
        second_input
    )

    assert first is created[0]
    assert second is created[1]

    assert (
        arx._active_workflow
        is second
    )

    assert (
        first is not second
    )


def test_run_does_not_replace_existing_workflow_if_new_construction_fails(
    tmp_path,
    monkeypatch,
):
    input_file = (
        tmp_path / "input.json"
    )

    input_file.write_text(
        "{}",
        encoding="utf-8",
    )

    existing = object()

    monkeypatch.setattr(
        arx,
        "_active_workflow",
        existing,
    )

    class FailingARXCLI:
        def __init__(
            self,
            input_path,
        ):
            raise RuntimeError(
                "construction failed"
            )

    monkeypatch.setattr(
        arx,
        "ARXCLI",
        FailingARXCLI,
    )

    with pytest.raises(
        RuntimeError,
        match="construction failed",
    ):
        arx.run(input_file)

    # Assignment happens only after successful ARXCLI construction.
    assert (
        arx._active_workflow
        is existing
    )


# =============================================================================
# Delegation before run()
# =============================================================================


@pytest.mark.parametrize(
    "api_call",
    [
        lambda: arx.show_molecules(),
        lambda: arx.show_functional_groups(),
        lambda: arx.show_reactions(),
        lambda: arx.select_reactions(),
        lambda: arx.show_non_reactants(),
        lambda: arx.select_non_reactants(),
        lambda: arx.prepare_reactions(),
        lambda: arx.show_reaction_templates(),
        lambda: arx.process(),
    ],
)
def test_public_api_requires_run_first(
    api_call,
):
    with pytest.raises(
        RuntimeError,
        match="No active session",
    ):
        api_call()


# =============================================================================
# Public API delegation
# =============================================================================


def test_show_molecules_delegates(
    monkeypatch,
):
    workflow = install_fake_workflow(
        monkeypatch
    )

    result = arx.show_molecules()

    assert (
        result
        == "molecules-image"
    )

    assert workflow.calls == [
        ("show_molecules",)
    ]


def test_show_functional_groups_delegates(
    monkeypatch,
):
    workflow = install_fake_workflow(
        monkeypatch
    )

    result = (
        arx.show_functional_groups()
    )

    assert (
        result
        == "functional-groups-image"
    )

    assert workflow.calls == [
        ("show_functional_groups",)
    ]


def test_show_reactions_delegates(
    monkeypatch,
):
    workflow = install_fake_workflow(
        monkeypatch
    )

    result = arx.show_reactions()

    assert (
        result
        == "reactions-image"
    )

    assert workflow.calls == [
        ("show_reactions",)
    ]


def test_select_reactions_delegates(
    monkeypatch,
):
    workflow = install_fake_workflow(
        monkeypatch
    )

    result = arx.select_reactions()

    assert result is None

    assert workflow.calls == [
        ("select_reactions",)
    ]


def test_show_non_reactants_delegates(
    monkeypatch,
):
    workflow = install_fake_workflow(
        monkeypatch
    )

    result = (
        arx.show_non_reactants()
    )

    assert (
        result
        == "non-reactants-image"
    )

    assert workflow.calls == [
        ("show_non_reactants",)
    ]


def test_show_non_reactants_propagates_none(
    monkeypatch,
):
    workflow = FakeWorkflow()

    def return_none():
        workflow.calls.append(
            ("show_non_reactants",)
        )

        return None

    workflow.show_non_reactants = (
        return_none
    )

    monkeypatch.setattr(
        arx,
        "_active_workflow",
        workflow,
    )

    assert (
        arx.show_non_reactants()
        is None
    )


def test_select_non_reactants_delegates(
    monkeypatch,
):
    workflow = install_fake_workflow(
        monkeypatch
    )

    result = (
        arx.select_non_reactants()
    )

    assert result is None

    assert workflow.calls == [
        ("select_non_reactants",)
    ]


def test_prepare_reactions_delegates(
    monkeypatch,
):
    workflow = install_fake_workflow(
        monkeypatch
    )

    result = (
        arx.prepare_reactions()
    )

    assert result is None

    assert workflow.calls == [
        ("prepare_reactions",)
    ]


def test_process_delegates(
    monkeypatch,
):
    workflow = install_fake_workflow(
        monkeypatch
    )

    result = arx.process()

    assert result is None

    assert workflow.calls == [
        ("process",)
    ]


# =============================================================================
# show_reaction_templates()
# =============================================================================


@pytest.mark.parametrize(
    "highlight_type",
    [
        "template",
        "edge",
        "initiators",
        "delete",
    ],
)
def test_show_reaction_templates_accepts_supported_types(
    monkeypatch,
    highlight_type,
):
    workflow = install_fake_workflow(
        monkeypatch
    )

    result = (
        arx.show_reaction_templates(
            highlight_type
        )
    )

    assert result == (
        f"templates-{highlight_type}"
    )

    assert workflow.calls == [
        (
            "show_reaction_templates",
            highlight_type,
        )
    ]


@pytest.mark.parametrize(
    "highlight_type, expected",
    [
        ("TEMPLATE", "template"),
        ("Edge", "edge"),
        ("INITIATORS", "initiators"),
        ("Delete", "delete"),
    ],
)
def test_show_reaction_templates_is_case_insensitive(
    monkeypatch,
    highlight_type,
    expected,
):
    workflow = install_fake_workflow(
        monkeypatch
    )

    result = (
        arx.show_reaction_templates(
            highlight_type
        )
    )

    assert result == (
        f"templates-{expected}"
    )

    assert workflow.calls == [
        (
            "show_reaction_templates",
            expected,
        )
    ]


def test_show_reaction_templates_defaults_to_template(
    monkeypatch,
):
    workflow = install_fake_workflow(
        monkeypatch
    )

    result = (
        arx.show_reaction_templates()
    )

    assert (
        result
        == "templates-template"
    )

    assert workflow.calls == [
        (
            "show_reaction_templates",
            "template",
        )
    ]


def test_show_reaction_templates_none_defaults_to_template(
    monkeypatch,
):
    workflow = install_fake_workflow(
        monkeypatch
    )

    result = (
        arx.show_reaction_templates(
            None
        )
    )

    assert (
        result
        == "templates-template"
    )

    assert workflow.calls == [
        (
            "show_reaction_templates",
            "template",
        )
    ]


@pytest.mark.parametrize(
    "highlight_type",
    [
        "",
        "wrong",
        "reaction",
        "atoms",
        "templates",
    ],
)
def test_show_reaction_templates_rejects_invalid_type(
    monkeypatch,
    highlight_type,
):
    workflow = install_fake_workflow(
        monkeypatch
    )

    if highlight_type == "":
        # Current API intentionally treats false-like/empty input
        # as the default "template".
        result = (
            arx.show_reaction_templates(
                highlight_type
            )
        )

        assert (
            result
            == "templates-template"
        )

        return

    with pytest.raises(
        ValueError,
        match="Invalid highlight_type",
    ):
        arx.show_reaction_templates(
            highlight_type
        )

    # Validation occurs before delegation.
    assert workflow.calls == []


def test_show_reaction_templates_error_lists_allowed_values(
    monkeypatch,
):
    install_fake_workflow(
        monkeypatch
    )

    with pytest.raises(
        ValueError,
    ) as exc_info:
        arx.show_reaction_templates(
            "bad"
        )

    message = str(
        exc_info.value
    )

    assert "template" in message
    assert "edge" in message
    assert "initiators" in message
    assert "delete" in message


# =============================================================================
# __all__ / public export contract
# =============================================================================


def test_public_all_contains_package_metadata():
    expected = {
        "__title__",
        "__version__",
        "__release__",
        "__authors__",
        "__license__",
    }

    assert expected.issubset(
        set(arx.__all__)
    )


def test_public_all_contains_user_workflow_commands():
    expected = {
        "run",
        "show_molecules",
        "show_functional_groups",
        "show_reactions",
        "select_reactions",
        "show_non_reactants",
        "select_non_reactants",
        "prepare_reactions",
        "show_reaction_templates",
        "process",
    }

    assert expected.issubset(
        set(arx.__all__)
    )


def test_session_is_part_of_public_api():
    """
    session() is used directly by the documented/user-facing workflow:

        session = arx.session()

    Therefore it should be exported alongside the other public API helpers.
    """
    assert "session" in arx.__all__


def test_public_all_has_no_duplicates():
    assert len(arx.__all__) == len(
        set(arx.__all__)
    )


def test_every_name_in_public_all_exists():
    for name in arx.__all__:
        assert hasattr(
            arx,
            name,
        ), (
            f"{name!r} appears in "
            "__all__ but does not exist"
        )