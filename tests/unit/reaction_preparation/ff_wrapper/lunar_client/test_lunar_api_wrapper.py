from pathlib import Path
from types import SimpleNamespace

import pytest

import AutoREACTER.reaction_preparation.ff_wrapper.lunar_client.lunar_api_wrapper as lunar_api
from AutoREACTER.reaction_preparation.ff_wrapper.lunar_client.lunar_api_wrapper import (
    LunarAPIWrapper,
)


# =============================================================================
# Helpers
# =============================================================================


def make_session(
    tmp_path,
    *,
    force_field="PCFF",
):
    inputs = SimpleNamespace(
        force_field=force_field,
    )

    return SimpleNamespace(
        inputs=inputs,
        staging_dir=tmp_path / "staging",
    )


class FakeExecutor:
    def __init__(
        self,
        lunar_location,
        cache_dir,
    ):
        self.lunar_location = Path(
            lunar_location
        )

        self.cache_dir = Path(
            cache_dir
        )

        self.cache_atom_typing = (
            self.cache_dir
            / "atom_typing"
        )

        self.cache_all2lmp = (
            self.cache_dir
            / "all2lmp"
        )

        self.cache_bond_react_merge = (
            self.cache_dir
            / "bond_react_merge"
        )

        self.calls = []

    def run_atom_typing(
        self,
        updated_inputs,
        prepared_reactions,
        force_field,
    ):
        self.calls.append(
            (
                "atom_typing",
                updated_inputs,
                prepared_reactions,
                force_field,
            )
        )

        return [
            "atom-typing-result"
        ]

    def run_all2lmp(
        self,
        atom_typing_results,
        frc_file,
    ):
        self.calls.append(
            (
                "all2lmp",
                atom_typing_results,
                frc_file,
            )
        )

        return [
            "all2lmp-result"
        ]

    def run_bond_react_merge(
        self,
        merge_input_file_path,
        all2lmp_results,
    ):
        self.calls.append(
            (
                "bond_react_merge",
                merge_input_file_path,
                all2lmp_results,
            )
        )

        return "final-files"


# =============================================================================
# Constructor
# =============================================================================


def test_constructor_stores_session(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    monkeypatch.setattr(
        lunar_api,
        "get_LUNAR_loc",
        lambda use_gui=False:
        "/tmp/LUNAR",
    )

    monkeypatch.setattr(
        lunar_api,
        "LunarExecutor",
        FakeExecutor,
    )

    wrapper = LunarAPIWrapper(
        session
    )

    assert wrapper.session is session
    assert wrapper.inputs is session.inputs


def test_constructor_creates_lunar_cache_directory(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    monkeypatch.setattr(
        lunar_api,
        "get_LUNAR_loc",
        lambda use_gui=False:
        "/tmp/LUNAR",
    )

    monkeypatch.setattr(
        lunar_api,
        "LunarExecutor",
        FakeExecutor,
    )

    wrapper = LunarAPIWrapper(
        session
    )

    expected = (
        tmp_path
        / "staging"
        / "lunar"
    )

    assert wrapper.cache_dir == expected

    assert expected.is_dir()


def test_constructor_requests_lunar_location_without_gui(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    calls = []

    def fake_get_lunar_loc(
        use_gui=None,
    ):
        calls.append(
            use_gui
        )

        return "/tmp/LUNAR"

    monkeypatch.setattr(
        lunar_api,
        "get_LUNAR_loc",
        fake_get_lunar_loc,
    )

    monkeypatch.setattr(
        lunar_api,
        "LunarExecutor",
        FakeExecutor,
    )

    LunarAPIWrapper(
        session
    )

    assert calls == [
        False
    ]


def test_constructor_stores_lunar_location_as_path(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    monkeypatch.setattr(
        lunar_api,
        "get_LUNAR_loc",
        lambda use_gui=False:
        "/tmp/LUNAR",
    )

    monkeypatch.setattr(
        lunar_api,
        "LunarExecutor",
        FakeExecutor,
    )

    wrapper = LunarAPIWrapper(
        session
    )

    assert isinstance(
        wrapper.LUNAR_LOCATION,
        Path,
    )

    assert wrapper.LUNAR_LOCATION == (
        Path("/tmp/LUNAR")
    )


def test_constructor_initializes_executor_with_expected_paths(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    captured = []

    class CapturingExecutor:
        def __init__(
            self,
            lunar_location,
            cache_dir,
        ):
            captured.append(
                (
                    lunar_location,
                    cache_dir,
                )
            )

    monkeypatch.setattr(
        lunar_api,
        "get_LUNAR_loc",
        lambda use_gui=False:
        "/tmp/LUNAR",
    )

    monkeypatch.setattr(
        lunar_api,
        "LunarExecutor",
        CapturingExecutor,
    )

    wrapper = LunarAPIWrapper(
        session
    )

    assert captured == [
        (
            Path("/tmp/LUNAR"),
            tmp_path
            / "staging"
            / "lunar",
        )
    ]

    assert isinstance(
        wrapper.executor,
        CapturingExecutor,
    )


# =============================================================================
# lunar_workflow - loading screen
# =============================================================================


def test_lunar_workflow_calls_loading_screen(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    monkeypatch.setattr(
        lunar_api,
        "get_LUNAR_loc",
        lambda use_gui=False:
        "/tmp/LUNAR",
    )

    monkeypatch.setattr(
        lunar_api,
        "LunarExecutor",
        FakeExecutor,
    )

    loading_calls = []

    monkeypatch.setattr(
        lunar_api,
        "loading_screen",
        lambda name:
        loading_calls.append(
            name
        ),
    )

    monkeypatch.setattr(
        lunar_api,
        "get_force_field_file",
        lambda **kwargs:
        Path("/tmp/pcff.frc"),
    )

    monkeypatch.setattr(
        lunar_api,
        "write_bond_react_merge_input",
        lambda **kwargs:
        Path("/tmp/merge_input.txt"),
    )

    wrapper = LunarAPIWrapper(
        session
    )

    wrapper.lunar_workflow(
        session.inputs,
        [],
    )

    assert loading_calls == [
        "LUNAR Workflow"
    ]


# =============================================================================
# Force-field lookup
# =============================================================================


def test_lunar_workflow_uses_force_field_from_updated_inputs(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path,
        force_field="PCFF",
    )

    updated_inputs = SimpleNamespace(
        force_field="COMPASS",
    )

    monkeypatch.setattr(
        lunar_api,
        "get_LUNAR_loc",
        lambda use_gui=False:
        "/tmp/LUNAR",
    )

    monkeypatch.setattr(
        lunar_api,
        "LunarExecutor",
        FakeExecutor,
    )

    monkeypatch.setattr(
        lunar_api,
        "loading_screen",
        lambda name: None,
    )

    ff_calls = []

    def fake_get_force_field_file(
        *,
        force_field,
        lunar_location,
    ):
        ff_calls.append(
            (
                force_field,
                lunar_location,
            )
        )

        return Path(
            "/tmp/compass.frc"
        )

    monkeypatch.setattr(
        lunar_api,
        "get_force_field_file",
        fake_get_force_field_file,
    )

    monkeypatch.setattr(
        lunar_api,
        "write_bond_react_merge_input",
        lambda **kwargs:
        Path("/tmp/merge_input.txt"),
    )

    wrapper = LunarAPIWrapper(
        session
    )

    wrapper.lunar_workflow(
        updated_inputs,
        [],
    )

    assert ff_calls == [
        (
            "COMPASS",
            Path("/tmp/LUNAR"),
        )
    ]


# =============================================================================
# Stage 1 - atom typing
# =============================================================================


def test_lunar_workflow_passes_inputs_and_reactions_to_atom_typing(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    updated_inputs = SimpleNamespace(
        force_field="PCFF",
    )

    prepared_reactions = [
        object(),
        object(),
    ]

    monkeypatch.setattr(
        lunar_api,
        "get_LUNAR_loc",
        lambda use_gui=False:
        "/tmp/LUNAR",
    )

    monkeypatch.setattr(
        lunar_api,
        "LunarExecutor",
        FakeExecutor,
    )

    monkeypatch.setattr(
        lunar_api,
        "loading_screen",
        lambda name: None,
    )

    monkeypatch.setattr(
        lunar_api,
        "get_force_field_file",
        lambda **kwargs:
        Path("/tmp/pcff.frc"),
    )

    monkeypatch.setattr(
        lunar_api,
        "write_bond_react_merge_input",
        lambda **kwargs:
        Path("/tmp/merge_input.txt"),
    )

    wrapper = LunarAPIWrapper(
        session
    )

    wrapper.lunar_workflow(
        updated_inputs,
        prepared_reactions,
    )

    assert (
        wrapper.executor.calls[0]
        == (
            "atom_typing",
            updated_inputs,
            prepared_reactions,
            "PCFF",
        )
    )


# =============================================================================
# Stage 2 - all2lmp
# =============================================================================


def test_lunar_workflow_passes_atom_typing_results_to_all2lmp(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    monkeypatch.setattr(
        lunar_api,
        "get_LUNAR_loc",
        lambda use_gui=False:
        "/tmp/LUNAR",
    )

    monkeypatch.setattr(
        lunar_api,
        "LunarExecutor",
        FakeExecutor,
    )

    monkeypatch.setattr(
        lunar_api,
        "loading_screen",
        lambda name: None,
    )

    ff_file = Path(
        "/tmp/pcff.frc"
    )

    monkeypatch.setattr(
        lunar_api,
        "get_force_field_file",
        lambda **kwargs:
        ff_file,
    )

    monkeypatch.setattr(
        lunar_api,
        "write_bond_react_merge_input",
        lambda **kwargs:
        Path("/tmp/merge_input.txt"),
    )

    wrapper = LunarAPIWrapper(
        session
    )

    wrapper.lunar_workflow(
        session.inputs,
        [],
    )

    assert (
        wrapper.executor.calls[1]
        == (
            "all2lmp",
            [
                "atom-typing-result"
            ],
            ff_file,
        )
    )


# =============================================================================
# Stage 3 - merge input builder
# =============================================================================


def test_lunar_workflow_builds_merge_input_with_executor_caches(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    monkeypatch.setattr(
        lunar_api,
        "get_LUNAR_loc",
        lambda use_gui=False:
        "/tmp/LUNAR",
    )

    monkeypatch.setattr(
        lunar_api,
        "LunarExecutor",
        FakeExecutor,
    )

    monkeypatch.setattr(
        lunar_api,
        "loading_screen",
        lambda name: None,
    )

    monkeypatch.setattr(
        lunar_api,
        "get_force_field_file",
        lambda **kwargs:
        Path("/tmp/pcff.frc"),
    )

    merge_calls = []

    def fake_write_merge(
        *,
        cache_bond_react_merge,
        cache_all2lmp,
        all2lmp_results,
    ):
        merge_calls.append(
            (
                cache_bond_react_merge,
                cache_all2lmp,
                all2lmp_results,
            )
        )

        return Path(
            "/tmp/merge_input.txt"
        )

    monkeypatch.setattr(
        lunar_api,
        "write_bond_react_merge_input",
        fake_write_merge,
    )

    wrapper = LunarAPIWrapper(
        session
    )

    wrapper.lunar_workflow(
        session.inputs,
        [],
    )

    assert merge_calls == [
        (
            wrapper.executor.cache_bond_react_merge,
            wrapper.executor.cache_all2lmp,
            [
                "all2lmp-result"
            ],
        )
    ]


# =============================================================================
# Stage 4 - bond/react merge
# =============================================================================


def test_lunar_workflow_passes_merge_file_and_all2lmp_results_to_final_stage(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    monkeypatch.setattr(
        lunar_api,
        "get_LUNAR_loc",
        lambda use_gui=False:
        "/tmp/LUNAR",
    )

    monkeypatch.setattr(
        lunar_api,
        "LunarExecutor",
        FakeExecutor,
    )

    monkeypatch.setattr(
        lunar_api,
        "loading_screen",
        lambda name: None,
    )

    monkeypatch.setattr(
        lunar_api,
        "get_force_field_file",
        lambda **kwargs:
        Path("/tmp/pcff.frc"),
    )

    merge_input = Path(
        "/tmp/merge_input.txt"
    )

    monkeypatch.setattr(
        lunar_api,
        "write_bond_react_merge_input",
        lambda **kwargs:
        merge_input,
    )

    wrapper = LunarAPIWrapper(
        session
    )

    wrapper.lunar_workflow(
        session.inputs,
        [],
    )

    assert (
        wrapper.executor.calls[2]
        == (
            "bond_react_merge",
            merge_input,
            [
                "all2lmp-result"
            ],
        )
    )


# =============================================================================
# Workflow result
# =============================================================================


def test_lunar_workflow_returns_final_files(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    monkeypatch.setattr(
        lunar_api,
        "get_LUNAR_loc",
        lambda use_gui=False:
        "/tmp/LUNAR",
    )

    monkeypatch.setattr(
        lunar_api,
        "LunarExecutor",
        FakeExecutor,
    )

    monkeypatch.setattr(
        lunar_api,
        "loading_screen",
        lambda name: None,
    )

    monkeypatch.setattr(
        lunar_api,
        "get_force_field_file",
        lambda **kwargs:
        Path("/tmp/pcff.frc"),
    )

    monkeypatch.setattr(
        lunar_api,
        "write_bond_react_merge_input",
        lambda **kwargs:
        Path("/tmp/merge_input.txt"),
    )

    wrapper = LunarAPIWrapper(
        session
    )

    result = wrapper.lunar_workflow(
        session.inputs,
        [],
    )

    assert result == (
        "final-files"
    )


# =============================================================================
# Execution order
# =============================================================================


def test_lunar_workflow_executes_stages_in_order(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    events = []

    monkeypatch.setattr(
        lunar_api,
        "get_LUNAR_loc",
        lambda use_gui=False:
        "/tmp/LUNAR",
    )

    monkeypatch.setattr(
        lunar_api,
        "loading_screen",
        lambda name:
        events.append(
            "loading"
        ),
    )

    monkeypatch.setattr(
        lunar_api,
        "get_force_field_file",
        lambda **kwargs:
        (
            events.append(
                "force-field"
            )
            or Path(
                "/tmp/pcff.frc"
            )
        ),
    )

    class OrderedExecutor:
        def __init__(
            self,
            lunar_location,
            cache_dir,
        ):
            self.cache_all2lmp = (
                Path(cache_dir)
                / "all2lmp"
            )

            self.cache_bond_react_merge = (
                Path(cache_dir)
                / "bond_react_merge"
            )

        def run_atom_typing(
            self,
            **kwargs,
        ):
            events.append(
                "atom-typing"
            )

            return [
                "typed"
            ]

        def run_all2lmp(
            self,
            **kwargs,
        ):
            events.append(
                "all2lmp"
            )

            return [
                "converted"
            ]

        def run_bond_react_merge(
            self,
            **kwargs,
        ):
            events.append(
                "bond-react-merge"
            )

            return "final"

    monkeypatch.setattr(
        lunar_api,
        "LunarExecutor",
        OrderedExecutor,
    )

    monkeypatch.setattr(
        lunar_api,
        "write_bond_react_merge_input",
        lambda **kwargs:
        (
            events.append(
                "merge-builder"
            )
            or Path(
                "/tmp/merge_input.txt"
            )
        ),
    )

    wrapper = LunarAPIWrapper(
        session
    )

    result = wrapper.lunar_workflow(
        session.inputs,
        [],
    )

    assert result == "final"

    assert events == [
        "loading",
        "force-field",
        "atom-typing",
        "all2lmp",
        "merge-builder",
        "bond-react-merge",
    ]


# =============================================================================
# Failure propagation
# =============================================================================


def test_force_field_lookup_error_propagates(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    monkeypatch.setattr(
        lunar_api,
        "get_LUNAR_loc",
        lambda use_gui=False:
        "/tmp/LUNAR",
    )

    monkeypatch.setattr(
        lunar_api,
        "LunarExecutor",
        FakeExecutor,
    )

    monkeypatch.setattr(
        lunar_api,
        "loading_screen",
        lambda name: None,
    )

    def fail(**kwargs):
        raise FileNotFoundError(
            "force field missing"
        )

    monkeypatch.setattr(
        lunar_api,
        "get_force_field_file",
        fail,
    )

    wrapper = LunarAPIWrapper(
        session
    )

    with pytest.raises(
        FileNotFoundError,
        match="force field missing",
    ):
        wrapper.lunar_workflow(
            session.inputs,
            [],
        )


def test_atom_typing_error_stops_later_stages(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    monkeypatch.setattr(
        lunar_api,
        "get_LUNAR_loc",
        lambda use_gui=False:
        "/tmp/LUNAR",
    )

    monkeypatch.setattr(
        lunar_api,
        "loading_screen",
        lambda name: None,
    )

    monkeypatch.setattr(
        lunar_api,
        "get_force_field_file",
        lambda **kwargs:
        Path("/tmp/pcff.frc"),
    )

    class FailingExecutor:
        def __init__(
            self,
            lunar_location,
            cache_dir,
        ):
            self.cache_all2lmp = (
                Path(cache_dir)
                / "all2lmp"
            )

            self.cache_bond_react_merge = (
                Path(cache_dir)
                / "bond_react_merge"
            )

        def run_atom_typing(
            self,
            **kwargs,
        ):
            raise RuntimeError(
                "atom typing failed"
            )

        def run_all2lmp(
            self,
            **kwargs,
        ):
            pytest.fail(
                "all2lmp must not run"
            )

        def run_bond_react_merge(
            self,
            **kwargs,
        ):
            pytest.fail(
                "merge must not run"
            )

    monkeypatch.setattr(
        lunar_api,
        "LunarExecutor",
        FailingExecutor,
    )

    monkeypatch.setattr(
        lunar_api,
        "write_bond_react_merge_input",
        lambda **kwargs:
        pytest.fail(
            "merge builder must not run"
        ),
    )

    wrapper = LunarAPIWrapper(
        session
    )

    with pytest.raises(
        RuntimeError,
        match="atom typing failed",
    ):
        wrapper.lunar_workflow(
            session.inputs,
            [],
        )


def test_all2lmp_error_stops_merge_stages(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    monkeypatch.setattr(
        lunar_api,
        "get_LUNAR_loc",
        lambda use_gui=False:
        "/tmp/LUNAR",
    )

    monkeypatch.setattr(
        lunar_api,
        "loading_screen",
        lambda name: None,
    )

    monkeypatch.setattr(
        lunar_api,
        "get_force_field_file",
        lambda **kwargs:
        Path("/tmp/pcff.frc"),
    )

    class FailingExecutor:
        def __init__(
            self,
            lunar_location,
            cache_dir,
        ):
            self.cache_all2lmp = (
                Path(cache_dir)
                / "all2lmp"
            )

            self.cache_bond_react_merge = (
                Path(cache_dir)
                / "bond_react_merge"
            )

        def run_atom_typing(
            self,
            **kwargs,
        ):
            return [
                "typed"
            ]

        def run_all2lmp(
            self,
            **kwargs,
        ):
            raise RuntimeError(
                "all2lmp failed"
            )

        def run_bond_react_merge(
            self,
            **kwargs,
        ):
            pytest.fail(
                "bond merge must not run"
            )

    monkeypatch.setattr(
        lunar_api,
        "LunarExecutor",
        FailingExecutor,
    )

    monkeypatch.setattr(
        lunar_api,
        "write_bond_react_merge_input",
        lambda **kwargs:
        pytest.fail(
            "merge builder must not run"
        ),
    )

    wrapper = LunarAPIWrapper(
        session
    )

    with pytest.raises(
        RuntimeError,
        match="all2lmp failed",
    ):
        wrapper.lunar_workflow(
            session.inputs,
            [],
        )


def test_merge_builder_error_stops_final_merge(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    monkeypatch.setattr(
        lunar_api,
        "get_LUNAR_loc",
        lambda use_gui=False:
        "/tmp/LUNAR",
    )

    monkeypatch.setattr(
        lunar_api,
        "loading_screen",
        lambda name: None,
    )

    monkeypatch.setattr(
        lunar_api,
        "get_force_field_file",
        lambda **kwargs:
        Path("/tmp/pcff.frc"),
    )

    class Executor:
        def __init__(
            self,
            lunar_location,
            cache_dir,
        ):
            self.cache_all2lmp = (
                Path(cache_dir)
                / "all2lmp"
            )

            self.cache_bond_react_merge = (
                Path(cache_dir)
                / "bond_react_merge"
            )

        def run_atom_typing(
            self,
            **kwargs,
        ):
            return [
                "typed"
            ]

        def run_all2lmp(
            self,
            **kwargs,
        ):
            return [
                "converted"
            ]

        def run_bond_react_merge(
            self,
            **kwargs,
        ):
            pytest.fail(
                "final merge must not run"
            )

    monkeypatch.setattr(
        lunar_api,
        "LunarExecutor",
        Executor,
    )

    def fail_merge_builder(
        **kwargs,
    ):
        raise ValueError(
            "bad merge input"
        )

    monkeypatch.setattr(
        lunar_api,
        "write_bond_react_merge_input",
        fail_merge_builder,
    )

    wrapper = LunarAPIWrapper(
        session
    )

    with pytest.raises(
        ValueError,
        match="bad merge input",
    ):
        wrapper.lunar_workflow(
            session.inputs,
            [],
        )