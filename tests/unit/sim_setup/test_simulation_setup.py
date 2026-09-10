from types import SimpleNamespace

import pytest

import AutoREACTER.sim_setup.simulation_setup as sim_setup_module
from AutoREACTER.sim_setup.simulation_setup import (
    SimulationSetupManager,
)


# =============================================================================
# Helpers
# =============================================================================


def make_session(
    tmp_path,
    *,
    inputs=None,
    reacter_files=None,
):
    if inputs is None:
        inputs = SimpleNamespace(
            simulation_name="Polymer",
            simulations=[],
            write_second_reaction_stage=False,
        )

    if reacter_files is None:
        reacter_files = SimpleNamespace()

    return SimpleNamespace(
        inputs=inputs,
        reacter_files=reacter_files,
        output_dir=tmp_path,
    )


# =============================================================================
# setup_and_write_simulation
# =============================================================================


def test_setup_and_write_simulation_creates_property_calculator(
    tmp_path,
    monkeypatch,
):
    inputs = SimpleNamespace(
        simulation_name="Polymer",
    )

    session = make_session(
        tmp_path,
        inputs=inputs,
    )

    captured = {}

    class FakeCalculator:
        def __init__(self, setup):
            captured["setup"] = setup

        def process_all(self):
            return inputs

    monkeypatch.setattr(
        sim_setup_module,
        "SystemPropertyCalculations",
        FakeCalculator,
    )

    monkeypatch.setattr(
        SimulationSetupManager,
        "generate_input_files",
        lambda self, **kwargs: None,
    )

    manager = SimulationSetupManager()

    manager.setup_and_write_simulation(
        session
    )

    assert (
        captured["setup"]
        is inputs
    )


def test_setup_and_write_simulation_calls_process_all_once(
    tmp_path,
    monkeypatch,
):
    inputs = SimpleNamespace(
        simulation_name="Polymer",
    )

    session = make_session(
        tmp_path,
        inputs=inputs,
    )

    calls = []

    class FakeCalculator:
        def __init__(self, setup):
            self.setup = setup

        def process_all(self):
            calls.append(
                "process_all"
            )

            return self.setup

    monkeypatch.setattr(
        sim_setup_module,
        "SystemPropertyCalculations",
        FakeCalculator,
    )

    monkeypatch.setattr(
        SimulationSetupManager,
        "generate_input_files",
        lambda self, **kwargs: None,
    )

    manager = SimulationSetupManager()

    manager.setup_and_write_simulation(
        session
    )

    assert calls == [
        "process_all"
    ]


def test_setup_and_write_simulation_passes_updated_setup_to_writer_layer(
    tmp_path,
    monkeypatch,
):
    original_setup = SimpleNamespace(
        simulation_name="Original",
    )

    updated_setup = SimpleNamespace(
        simulation_name="Updated",
    )

    reacter_files = SimpleNamespace(
        template_files=[]
    )

    session = make_session(
        tmp_path,
        inputs=original_setup,
        reacter_files=reacter_files,
    )

    class FakeCalculator:
        def __init__(self, setup):
            assert (
                setup
                is original_setup
            )

        def process_all(self):
            return updated_setup

    monkeypatch.setattr(
        sim_setup_module,
        "SystemPropertyCalculations",
        FakeCalculator,
    )

    captured = {}

    def fake_generate(
        self,
        *,
        setup,
        reacter_files,
        run_dir,
    ):
        captured["setup"] = setup
        captured["reacter_files"] = reacter_files
        captured["run_dir"] = run_dir

    monkeypatch.setattr(
        SimulationSetupManager,
        "generate_input_files",
        fake_generate,
    )

    manager = SimulationSetupManager()

    manager.setup_and_write_simulation(
        session
    )

    assert (
        captured["setup"]
        is updated_setup
    )

    assert (
        captured["reacter_files"]
        is reacter_files
    )

    assert (
        captured["run_dir"]
        == tmp_path
    )


def test_setup_and_write_simulation_uses_session_reacter_files(
    tmp_path,
    monkeypatch,
):
    setup = SimpleNamespace(
        simulation_name="Polymer",
    )

    reacter_files = SimpleNamespace(
        marker="reaction-files"
    )

    session = make_session(
        tmp_path,
        inputs=setup,
        reacter_files=reacter_files,
    )

    class FakeCalculator:
        def __init__(self, setup):
            self.setup = setup

        def process_all(self):
            return self.setup

    monkeypatch.setattr(
        sim_setup_module,
        "SystemPropertyCalculations",
        FakeCalculator,
    )

    captured = {}

    def fake_generate(
        self,
        *,
        setup,
        reacter_files,
        run_dir,
    ):
        captured["reacter_files"] = (
            reacter_files
        )

    monkeypatch.setattr(
        SimulationSetupManager,
        "generate_input_files",
        fake_generate,
    )

    manager = SimulationSetupManager()

    manager.setup_and_write_simulation(
        session
    )

    assert (
        captured["reacter_files"]
        is reacter_files
    )


def test_setup_and_write_simulation_uses_session_output_dir(
    tmp_path,
    monkeypatch,
):
    setup = SimpleNamespace(
        simulation_name="Polymer",
    )

    session = make_session(
        tmp_path,
        inputs=setup,
    )

    class FakeCalculator:
        def __init__(self, setup):
            self.setup = setup

        def process_all(self):
            return self.setup

    monkeypatch.setattr(
        sim_setup_module,
        "SystemPropertyCalculations",
        FakeCalculator,
    )

    captured = {}

    def fake_generate(
        self,
        *,
        setup,
        reacter_files,
        run_dir,
    ):
        captured["run_dir"] = run_dir

    monkeypatch.setattr(
        SimulationSetupManager,
        "generate_input_files",
        fake_generate,
    )

    manager = SimulationSetupManager()

    manager.setup_and_write_simulation(
        session
    )

    assert (
        captured["run_dir"]
        == tmp_path
    )


def test_setup_and_write_simulation_calls_calculation_before_generation(
    tmp_path,
    monkeypatch,
):
    setup = SimpleNamespace(
        simulation_name="Polymer",
    )

    session = make_session(
        tmp_path,
        inputs=setup,
    )

    events = []

    class FakeCalculator:
        def __init__(self, setup):
            events.append(
                "calculator_init"
            )

            self.setup = setup

        def process_all(self):
            events.append(
                "process_all"
            )

            return self.setup

    monkeypatch.setattr(
        sim_setup_module,
        "SystemPropertyCalculations",
        FakeCalculator,
    )

    def fake_generate(
        self,
        **kwargs,
    ):
        events.append(
            "generate_input_files"
        )

    monkeypatch.setattr(
        SimulationSetupManager,
        "generate_input_files",
        fake_generate,
    )

    manager = SimulationSetupManager()

    manager.setup_and_write_simulation(
        session
    )

    assert events == [
        "calculator_init",
        "process_all",
        "generate_input_files",
    ]


def test_setup_and_write_simulation_returns_none(
    tmp_path,
    monkeypatch,
):
    setup = SimpleNamespace(
        simulation_name="Polymer",
    )

    session = make_session(
        tmp_path,
        inputs=setup,
    )

    class FakeCalculator:
        def __init__(self, setup):
            self.setup = setup

        def process_all(self):
            return self.setup

    monkeypatch.setattr(
        sim_setup_module,
        "SystemPropertyCalculations",
        FakeCalculator,
    )

    monkeypatch.setattr(
        SimulationSetupManager,
        "generate_input_files",
        lambda self, **kwargs: None,
    )

    manager = SimulationSetupManager()

    result = (
        manager
        .setup_and_write_simulation(
            session
        )
    )

    assert result is None


def test_setup_and_write_simulation_prints_success_message(
    tmp_path,
    monkeypatch,
    capsys,
):
    setup = SimpleNamespace(
        simulation_name="Polymer",
    )

    session = make_session(
        tmp_path,
        inputs=setup,
    )

    class FakeCalculator:
        def __init__(self, setup):
            self.setup = setup

        def process_all(self):
            return self.setup

    monkeypatch.setattr(
        sim_setup_module,
        "SystemPropertyCalculations",
        FakeCalculator,
    )

    monkeypatch.setattr(
        SimulationSetupManager,
        "generate_input_files",
        lambda self, **kwargs: None,
    )

    manager = SimulationSetupManager()

    manager.setup_and_write_simulation(
        session
    )

    captured = (
        capsys.readouterr()
    )

    expected_path = (
        tmp_path
        / "LAMMPS_input_files"
    )

    assert (
        "[SUCCESS]"
        in captured.out
    )

    assert (
        "All 5 simulation stages written"
        in captured.out
    )

    assert (
        str(expected_path)
        in captured.out
    )


# =============================================================================
# Failure propagation
# =============================================================================


def test_calculation_error_propagates_and_generation_is_not_called(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    class FakeCalculator:
        def __init__(self, setup):
            pass

        def process_all(self):
            raise ValueError(
                "calculation failed"
            )

    monkeypatch.setattr(
        sim_setup_module,
        "SystemPropertyCalculations",
        FakeCalculator,
    )

    called = []

    monkeypatch.setattr(
        SimulationSetupManager,
        "generate_input_files",
        lambda self, **kwargs:
        called.append(True),
    )

    manager = SimulationSetupManager()

    with pytest.raises(
        ValueError,
        match="calculation failed",
    ):
        manager.setup_and_write_simulation(
            session
        )

    assert called == []


def test_generation_error_propagates(
    tmp_path,
    monkeypatch,
):
    setup = SimpleNamespace(
        simulation_name="Polymer",
    )

    session = make_session(
        tmp_path,
        inputs=setup,
    )

    class FakeCalculator:
        def __init__(self, setup):
            self.setup = setup

        def process_all(self):
            return self.setup

    monkeypatch.setattr(
        sim_setup_module,
        "SystemPropertyCalculations",
        FakeCalculator,
    )

    def fake_generate(
        self,
        **kwargs,
    ):
        raise RuntimeError(
            "writer failed"
        )

    monkeypatch.setattr(
        SimulationSetupManager,
        "generate_input_files",
        fake_generate,
    )

    manager = SimulationSetupManager()

    with pytest.raises(
        RuntimeError,
        match="writer failed",
    ):
        manager.setup_and_write_simulation(
            session
        )


def test_success_message_not_printed_when_generation_fails(
    tmp_path,
    monkeypatch,
    capsys,
):
    setup = SimpleNamespace(
        simulation_name="Polymer",
    )

    session = make_session(
        tmp_path,
        inputs=setup,
    )

    class FakeCalculator:
        def __init__(self, setup):
            self.setup = setup

        def process_all(self):
            return self.setup

    monkeypatch.setattr(
        sim_setup_module,
        "SystemPropertyCalculations",
        FakeCalculator,
    )

    def fake_generate(
        self,
        **kwargs,
    ):
        raise RuntimeError(
            "writer failed"
        )

    monkeypatch.setattr(
        SimulationSetupManager,
        "generate_input_files",
        fake_generate,
    )

    manager = SimulationSetupManager()

    with pytest.raises(
        RuntimeError
    ):
        manager.setup_and_write_simulation(
            session
        )

    captured = (
        capsys.readouterr()
    )

    assert (
        "[SUCCESS]"
        not in captured.out
    )


# =============================================================================
# generate_input_files
# =============================================================================


def test_generate_input_files_constructs_writer_with_reacter_files(
    tmp_path,
    monkeypatch,
):
    reacter_files = SimpleNamespace(
        marker="rf"
    )

    setup = SimpleNamespace(
        simulation_name="Polymer",
    )

    captured = {}

    class FakeWriter:
        def __init__(
            self,
            *,
            reacter_files,
        ):
            captured[
                "reacter_files"
            ] = reacter_files

        def write_all_files(
            self,
            run_dir,
            setup,
        ):
            pass

    monkeypatch.setattr(
        sim_setup_module,
        "Writer",
        FakeWriter,
    )

    manager = SimulationSetupManager()

    manager.generate_input_files(
        setup=setup,
        reacter_files=reacter_files,
        run_dir=tmp_path,
    )

    assert (
        captured["reacter_files"]
        is reacter_files
    )


def test_generate_input_files_calls_write_all_files_with_expected_arguments(
    tmp_path,
    monkeypatch,
):
    reacter_files = SimpleNamespace(
        marker="rf"
    )

    setup = SimpleNamespace(
        simulation_name="Polymer",
        write_second_reaction_stage=False,
    )

    captured = {}

    class FakeWriter:
        def __init__(
            self,
            *,
            reacter_files,
        ):
            self.reacter_files = (
                reacter_files
            )

        def write_all_files(
            self,
            run_dir,
            passed_setup,
        ):
            captured["run_dir"] = (
                run_dir
            )

            captured["setup"] = (
                passed_setup
            )

    monkeypatch.setattr(
        sim_setup_module,
        "Writer",
        FakeWriter,
    )

    manager = SimulationSetupManager()

    manager.generate_input_files(
        setup=setup,
        reacter_files=reacter_files,
        run_dir=tmp_path,
    )

    assert (
        captured["run_dir"]
        == tmp_path
    )

    assert (
        captured["setup"]
        is setup
    )


def test_generate_input_files_does_not_modify_second_stage_flag(
    tmp_path,
    monkeypatch,
):
    """
    SimulationSetupManager passes the setup object through unchanged.

    Whether stage 2 is generated belongs to Writer.
    """
    reacter_files = SimpleNamespace()

    setup = SimpleNamespace(
        simulation_name="Polymer",
        write_second_reaction_stage=False,
    )

    observed = {}

    class FakeWriter:
        def __init__(
            self,
            *,
            reacter_files,
        ):
            pass

        def write_all_files(
            self,
            run_dir,
            passed_setup,
        ):
            observed["flag"] = (
                passed_setup
                .write_second_reaction_stage
            )

    monkeypatch.setattr(
        sim_setup_module,
        "Writer",
        FakeWriter,
    )

    manager = SimulationSetupManager()

    manager.generate_input_files(
        setup=setup,
        reacter_files=reacter_files,
        run_dir=tmp_path,
    )

    assert (
        observed["flag"]
        is False
    )

    assert (
        setup.write_second_reaction_stage
        is False
    )


def test_generate_input_files_returns_none(
    tmp_path,
    monkeypatch,
):
    class FakeWriter:
        def __init__(
            self,
            *,
            reacter_files,
        ):
            pass

        def write_all_files(
            self,
            run_dir,
            setup,
        ):
            return None

    monkeypatch.setattr(
        sim_setup_module,
        "Writer",
        FakeWriter,
    )

    manager = SimulationSetupManager()

    result = (
        manager
        .generate_input_files(
            setup=SimpleNamespace(),
            reacter_files=(
                SimpleNamespace()
            ),
            run_dir=tmp_path,
        )
    )

    assert result is None


def test_writer_error_propagates_from_generate_input_files(
    tmp_path,
    monkeypatch,
):
    class FakeWriter:
        def __init__(
            self,
            *,
            reacter_files,
        ):
            pass

        def write_all_files(
            self,
            run_dir,
            setup,
        ):
            raise RuntimeError(
                "LAMMPS writing failed"
            )

    monkeypatch.setattr(
        sim_setup_module,
        "Writer",
        FakeWriter,
    )

    manager = SimulationSetupManager()

    with pytest.raises(
        RuntimeError,
        match="LAMMPS writing failed",
    ):
        manager.generate_input_files(
            setup=SimpleNamespace(),
            reacter_files=(
                SimpleNamespace()
            ),
            run_dir=tmp_path,
        )