from pathlib import Path
from types import SimpleNamespace

import pytest

import AutoREACTER.sim_setup.writers.writer as writer_module
from AutoREACTER.sim_setup.writers.writer import Writer


# =============================================================================
# Helpers
# =============================================================================


def make_reacter_files():
    return SimpleNamespace(
        in_file=Path("in.create_atoms.script"),
        force_field_data=Path("force_field.data"),
        molecule_files=[],
        template_files=[],
    )


def make_simulation(
    *,
    tag="sim1",
    temperature=500.0,
):
    return SimpleNamespace(
        tag=tag,
        temperature=temperature,
    )


def make_setup(
    *,
    simulation_name="Polymer",
    simulations=None,
    write_second_reaction_stage=True,
):
    return SimpleNamespace(
        simulation_name=simulation_name,
        simulations=list(
            simulations or []
        ),
        write_second_reaction_stage=write_second_reaction_stage,
    )


class FakeInitialSettings:
    instances = []

    def __init__(self, reacter_files):
        self.reacter_files = reacter_files
        self.get_calls = 0

        self.settings = SimpleNamespace(
            units="real",
        )

        self.__class__.instances.append(
            self
        )

    def get_LUNAR_lammps_settings(self):
        self.get_calls += 1

        return self.settings


# =============================================================================
# Constructor
# =============================================================================


def test_constructor_stores_reacter_files(
    monkeypatch,
):
    FakeInitialSettings.instances.clear()

    monkeypatch.setattr(
        writer_module,
        "LammpsInitialSettings",
        FakeInitialSettings,
    )

    reacter_files = (
        make_reacter_files()
    )

    writer = Writer(
        reacter_files=reacter_files
    )

    assert (
        writer.reacter_files
        is reacter_files
    )


def test_constructor_creates_lammps_initial_settings(
    monkeypatch,
):
    FakeInitialSettings.instances.clear()

    monkeypatch.setattr(
        writer_module,
        "LammpsInitialSettings",
        FakeInitialSettings,
    )

    reacter_files = (
        make_reacter_files()
    )

    writer = Writer(
        reacter_files=reacter_files
    )

    assert (
        len(FakeInitialSettings.instances)
        == 1
    )

    instance = (
        FakeInitialSettings.instances[0]
    )

    assert (
        instance.reacter_files
        is reacter_files
    )

    assert (
        writer.lammps_initial_setup
        is instance
    )


def test_constructor_gets_lunar_settings_once(
    monkeypatch,
):
    FakeInitialSettings.instances.clear()

    monkeypatch.setattr(
        writer_module,
        "LammpsInitialSettings",
        FakeInitialSettings,
    )

    writer = Writer(
        reacter_files=(
            make_reacter_files()
        )
    )

    instance = (
        FakeInitialSettings.instances[0]
    )

    assert (
        instance.get_calls
        == 1
    )

    assert (
        writer.settings
        is instance.settings
    )


# =============================================================================
# Main output directories
# =============================================================================


def test_write_all_files_creates_lammps_input_files_directory(
    tmp_path,
    monkeypatch,
):
    FakeInitialSettings.instances.clear()

    monkeypatch.setattr(
        writer_module,
        "LammpsInitialSettings",
        FakeInitialSettings,
    )

    writer = Writer(
        make_reacter_files()
    )

    setup = make_setup(
        simulations=[]
    )

    writer.write_all_files(
        tmp_path,
        setup,
    )

    assert (
        tmp_path
        / "LAMMPS_input_files"
    ).is_dir()


def test_write_all_files_creates_simulation_subdirectory(
    tmp_path,
    monkeypatch,
):
    FakeInitialSettings.instances.clear()

    monkeypatch.setattr(
        writer_module,
        "LammpsInitialSettings",
        FakeInitialSettings,
    )

    # Prevent real writer execution.
    monkeypatch.setattr(
        writer_module,
        "DensificationWriter",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        writer_module,
        "PreEqWriter",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        writer_module,
        "RxnFirstStageWriter",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        writer_module,
        "RxnSecondStageWriter",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        writer_module,
        "PostEqWriter",
        lambda **kwargs: None,
    )

    writer = Writer(
        make_reacter_files()
    )

    setup = make_setup(
        simulation_name="Styrene",
        simulations=[
            make_simulation(
                tag="298K"
            )
        ],
    )

    writer.write_all_files(
        tmp_path,
        setup,
    )

    assert (
        tmp_path
        / "LAMMPS_input_files"
        / "Styrene_298K"
    ).is_dir()


def test_each_simulation_gets_own_subdirectory(
    tmp_path,
    monkeypatch,
):
    FakeInitialSettings.instances.clear()

    monkeypatch.setattr(
        writer_module,
        "LammpsInitialSettings",
        FakeInitialSettings,
    )

    monkeypatch.setattr(
        writer_module,
        "DensificationWriter",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        writer_module,
        "PreEqWriter",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        writer_module,
        "RxnFirstStageWriter",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        writer_module,
        "RxnSecondStageWriter",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        writer_module,
        "PostEqWriter",
        lambda **kwargs: None,
    )

    writer = Writer(
        make_reacter_files()
    )

    setup = make_setup(
        simulation_name="Styrene",
        simulations=[
            make_simulation(
                tag="298K"
            ),
            make_simulation(
                tag="373K"
            ),
            make_simulation(
                tag="523K"
            ),
        ],
    )

    writer.write_all_files(
        tmp_path,
        setup,
    )

    base = (
        tmp_path
        / "LAMMPS_input_files"
    )

    assert (
        base / "Styrene_298K"
    ).is_dir()

    assert (
        base / "Styrene_373K"
    ).is_dir()

    assert (
        base / "Styrene_523K"
    ).is_dir()


# =============================================================================
# Stage orchestration
# =============================================================================


def test_all_five_stages_are_called_when_second_stage_enabled(
    tmp_path,
    monkeypatch,
):
    FakeInitialSettings.instances.clear()

    monkeypatch.setattr(
        writer_module,
        "LammpsInitialSettings",
        FakeInitialSettings,
    )

    events = []

    def fake_stage(name):
        def factory(**kwargs):
            events.append(
                (
                    name,
                    kwargs,
                )
            )

            return SimpleNamespace()

        return factory

    monkeypatch.setattr(
        writer_module,
        "DensificationWriter",
        fake_stage("densification"),
    )

    monkeypatch.setattr(
        writer_module,
        "PreEqWriter",
        fake_stage("pre_eq"),
    )

    monkeypatch.setattr(
        writer_module,
        "RxnFirstStageWriter",
        fake_stage("rxn_first"),
    )

    monkeypatch.setattr(
        writer_module,
        "RxnSecondStageWriter",
        fake_stage("rxn_second"),
    )

    monkeypatch.setattr(
        writer_module,
        "PostEqWriter",
        fake_stage("post_eq"),
    )

    simulation = make_simulation()

    setup = make_setup(
        simulations=[
            simulation
        ],
        write_second_reaction_stage=True,
    )

    writer = Writer(
        make_reacter_files()
    )

    writer.write_all_files(
        tmp_path,
        setup,
    )

    assert [
        name
        for name, _
        in events
    ] == [
        "densification",
        "pre_eq",
        "rxn_first",
        "rxn_second",
        "post_eq",
    ]


def test_second_stage_is_skipped_when_disabled(
    tmp_path,
    monkeypatch,
):
    FakeInitialSettings.instances.clear()

    monkeypatch.setattr(
        writer_module,
        "LammpsInitialSettings",
        FakeInitialSettings,
    )

    events = []

    def fake_stage(name):
        def factory(**kwargs):
            events.append(name)
            return SimpleNamespace()

        return factory

    monkeypatch.setattr(
        writer_module,
        "DensificationWriter",
        fake_stage("densification"),
    )

    monkeypatch.setattr(
        writer_module,
        "PreEqWriter",
        fake_stage("pre_eq"),
    )

    monkeypatch.setattr(
        writer_module,
        "RxnFirstStageWriter",
        fake_stage("rxn_first"),
    )

    monkeypatch.setattr(
        writer_module,
        "RxnSecondStageWriter",
        fake_stage("rxn_second"),
    )

    monkeypatch.setattr(
        writer_module,
        "PostEqWriter",
        fake_stage("post_eq"),
    )

    writer = Writer(
        make_reacter_files()
    )

    writer.write_all_files(
        tmp_path,
        make_setup(
            simulations=[
                make_simulation()
            ],
            write_second_reaction_stage=False,
        ),
    )

    assert events == [
        "densification",
        "pre_eq",
        "rxn_first",
        "post_eq",
    ]

    assert (
        "rxn_second"
        not in events
    )


# =============================================================================
# Common arguments
# =============================================================================


def test_densification_receives_expected_arguments(
    tmp_path,
    monkeypatch,
):
    FakeInitialSettings.instances.clear()

    monkeypatch.setattr(
        writer_module,
        "LammpsInitialSettings",
        FakeInitialSettings,
    )

    captured = {}

    monkeypatch.setattr(
        writer_module,
        "DensificationWriter",
        lambda **kwargs:
        captured.update(kwargs),
    )

    monkeypatch.setattr(
        writer_module,
        "PreEqWriter",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        writer_module,
        "RxnFirstStageWriter",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        writer_module,
        "RxnSecondStageWriter",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        writer_module,
        "PostEqWriter",
        lambda **kwargs: None,
    )

    reacter_files = (
        make_reacter_files()
    )

    simulation = (
        make_simulation(
            tag="500K"
        )
    )

    writer = Writer(
        reacter_files
    )

    writer.write_all_files(
        tmp_path,
        make_setup(
            simulation_name="Polymer",
            simulations=[
                simulation
            ],
        ),
    )

    expected_dir = (
        tmp_path
        / "LAMMPS_input_files"
        / "Polymer_500K"
    )

    assert (
        captured["out_dir"]
        == expected_dir
    )

    assert (
        captured["settings"]
        is writer.settings
    )

    assert (
        captured["reacter_files"]
        is reacter_files
    )

    assert (
        captured["simulation"]
        is simulation
    )

    assert (
        captured["sim_name"]
        == "Polymer"
    )


def test_pre_eq_receives_expected_arguments(
    tmp_path,
    monkeypatch,
):
    FakeInitialSettings.instances.clear()

    monkeypatch.setattr(
        writer_module,
        "LammpsInitialSettings",
        FakeInitialSettings,
    )

    captured = {}

    monkeypatch.setattr(
        writer_module,
        "DensificationWriter",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        writer_module,
        "PreEqWriter",
        lambda **kwargs:
        captured.update(kwargs),
    )

    monkeypatch.setattr(
        writer_module,
        "RxnFirstStageWriter",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        writer_module,
        "RxnSecondStageWriter",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        writer_module,
        "PostEqWriter",
        lambda **kwargs: None,
    )

    simulation = make_simulation()

    writer = Writer(
        make_reacter_files()
    )

    writer.write_all_files(
        tmp_path,
        make_setup(
            simulation_name="Polymer",
            simulations=[
                simulation
            ],
        ),
    )

    assert (
        captured["settings"]
        is writer.settings
    )

    assert (
        captured["simulation"]
        is simulation
    )

    assert (
        captured["sim_name"]
        == "Polymer"
    )

    assert (
        "reacter_files"
        not in captured
    )


def test_first_reaction_stage_receives_reacter_files(
    tmp_path,
    monkeypatch,
):
    FakeInitialSettings.instances.clear()

    monkeypatch.setattr(
        writer_module,
        "LammpsInitialSettings",
        FakeInitialSettings,
    )

    captured = {}

    monkeypatch.setattr(
        writer_module,
        "DensificationWriter",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        writer_module,
        "PreEqWriter",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        writer_module,
        "RxnFirstStageWriter",
        lambda **kwargs:
        captured.update(kwargs),
    )

    monkeypatch.setattr(
        writer_module,
        "RxnSecondStageWriter",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        writer_module,
        "PostEqWriter",
        lambda **kwargs: None,
    )

    reacter_files = (
        make_reacter_files()
    )

    writer = Writer(
        reacter_files
    )

    writer.write_all_files(
        tmp_path,
        make_setup(
            simulations=[
                make_simulation()
            ],
        ),
    )

    assert (
        captured["reacter_files"]
        is reacter_files
    )


def test_second_reaction_stage_receives_reacter_files(
    tmp_path,
    monkeypatch,
):
    FakeInitialSettings.instances.clear()

    monkeypatch.setattr(
        writer_module,
        "LammpsInitialSettings",
        FakeInitialSettings,
    )

    captured = {}

    monkeypatch.setattr(
        writer_module,
        "DensificationWriter",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        writer_module,
        "PreEqWriter",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        writer_module,
        "RxnFirstStageWriter",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        writer_module,
        "RxnSecondStageWriter",
        lambda **kwargs:
        captured.update(kwargs),
    )

    monkeypatch.setattr(
        writer_module,
        "PostEqWriter",
        lambda **kwargs: None,
    )

    reacter_files = (
        make_reacter_files()
    )

    writer = Writer(
        reacter_files
    )

    writer.write_all_files(
        tmp_path,
        make_setup(
            simulations=[
                make_simulation()
            ],
            write_second_reaction_stage=True,
        ),
    )

    assert (
        captured["reacter_files"]
        is reacter_files
    )


# =============================================================================
# Post-equilibration flag forwarding
# =============================================================================


@pytest.mark.parametrize(
    "write_second_stage",
    [
        True,
        False,
    ],
)
def test_post_eq_receives_second_stage_flag(
    tmp_path,
    monkeypatch,
    write_second_stage,
):
    FakeInitialSettings.instances.clear()

    monkeypatch.setattr(
        writer_module,
        "LammpsInitialSettings",
        FakeInitialSettings,
    )

    captured = {}

    monkeypatch.setattr(
        writer_module,
        "DensificationWriter",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        writer_module,
        "PreEqWriter",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        writer_module,
        "RxnFirstStageWriter",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        writer_module,
        "RxnSecondStageWriter",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        writer_module,
        "PostEqWriter",
        lambda **kwargs:
        captured.update(kwargs),
    )

    writer = Writer(
        make_reacter_files()
    )

    writer.write_all_files(
        tmp_path,
        make_setup(
            simulations=[
                make_simulation()
            ],
            write_second_reaction_stage=(
                write_second_stage
            ),
        ),
    )

    assert (
        captured[
            "write_second_reaction_stage"
        ]
        is write_second_stage
    )


# =============================================================================
# Multiple simulations
# =============================================================================


def test_all_stages_run_for_each_simulation(
    tmp_path,
    monkeypatch,
):
    FakeInitialSettings.instances.clear()

    monkeypatch.setattr(
        writer_module,
        "LammpsInitialSettings",
        FakeInitialSettings,
    )

    calls = []

    def fake_stage(name):
        def factory(**kwargs):
            calls.append(
                (
                    name,
                    kwargs["simulation"].tag,
                )
            )

            return SimpleNamespace()

        return factory

    monkeypatch.setattr(
        writer_module,
        "DensificationWriter",
        fake_stage("dense"),
    )

    monkeypatch.setattr(
        writer_module,
        "PreEqWriter",
        fake_stage("pre"),
    )

    monkeypatch.setattr(
        writer_module,
        "RxnFirstStageWriter",
        fake_stage("rxn1"),
    )

    monkeypatch.setattr(
        writer_module,
        "RxnSecondStageWriter",
        fake_stage("rxn2"),
    )

    monkeypatch.setattr(
        writer_module,
        "PostEqWriter",
        fake_stage("post"),
    )

    writer = Writer(
        make_reacter_files()
    )

    writer.write_all_files(
        tmp_path,
        make_setup(
            simulations=[
                make_simulation(
                    tag="298K"
                ),
                make_simulation(
                    tag="500K"
                ),
            ],
            write_second_reaction_stage=True,
        ),
    )

    assert calls == [
        ("dense", "298K"),
        ("pre", "298K"),
        ("rxn1", "298K"),
        ("rxn2", "298K"),
        ("post", "298K"),
        ("dense", "500K"),
        ("pre", "500K"),
        ("rxn1", "500K"),
        ("rxn2", "500K"),
        ("post", "500K"),
    ]


def test_same_settings_object_is_shared_across_all_stages(
    tmp_path,
    monkeypatch,
):
    FakeInitialSettings.instances.clear()

    monkeypatch.setattr(
        writer_module,
        "LammpsInitialSettings",
        FakeInitialSettings,
    )

    received_settings = []

    def fake_stage(**kwargs):
        received_settings.append(
            kwargs["settings"]
        )

        return SimpleNamespace()

    monkeypatch.setattr(
        writer_module,
        "DensificationWriter",
        fake_stage,
    )

    monkeypatch.setattr(
        writer_module,
        "PreEqWriter",
        fake_stage,
    )

    monkeypatch.setattr(
        writer_module,
        "RxnFirstStageWriter",
        fake_stage,
    )

    monkeypatch.setattr(
        writer_module,
        "RxnSecondStageWriter",
        fake_stage,
    )

    monkeypatch.setattr(
        writer_module,
        "PostEqWriter",
        fake_stage,
    )

    writer = Writer(
        make_reacter_files()
    )

    writer.write_all_files(
        tmp_path,
        make_setup(
            simulations=[
                make_simulation()
            ],
            write_second_reaction_stage=True,
        ),
    )

    assert len(
        received_settings
    ) == 5

    assert all(
        settings
        is writer.settings
        for settings
        in received_settings
    )


# =============================================================================
# Empty setup
# =============================================================================


def test_empty_simulation_list_creates_only_root_lammps_directory(
    tmp_path,
    monkeypatch,
):
    FakeInitialSettings.instances.clear()

    monkeypatch.setattr(
        writer_module,
        "LammpsInitialSettings",
        FakeInitialSettings,
    )

    writer = Writer(
        make_reacter_files()
    )

    writer.write_all_files(
        tmp_path,
        make_setup(
            simulations=[]
        ),
    )

    root = (
        tmp_path
        / "LAMMPS_input_files"
    )

    assert root.is_dir()

    assert list(
        root.iterdir()
    ) == []


def test_write_all_files_returns_none(
    tmp_path,
    monkeypatch,
):
    FakeInitialSettings.instances.clear()

    monkeypatch.setattr(
        writer_module,
        "LammpsInitialSettings",
        FakeInitialSettings,
    )

    writer = Writer(
        make_reacter_files()
    )

    result = writer.write_all_files(
        tmp_path,
        make_setup(
            simulations=[]
        ),
    )

    assert result is None