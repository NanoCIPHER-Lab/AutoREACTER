from pathlib import Path
from types import SimpleNamespace

import pytest

from AutoREACTER.sim_setup.writers.pre_eq_writer import (
    PreEqWriter,
)


# =============================================================================
# Helpers
# =============================================================================


def make_settings(
    *,
    neighbor=None,
    neigh_modify=None,
):
    return SimpleNamespace(
        units="real",
        dimension="3",
        boundary="p p p",
        atom_style="full",
        bond_style="class2",
        angle_style="class2",
        dihedral_style="class2",
        improper_style="class2",
        special_bonds="lj/coul 0 0 1",
        pair_style="lj/class2/coul/long 12.0",
        kspace_style="pppm 1.0e-4",
        pair_modify="mix sixthpower",
        neighbor=neighbor,
        neigh_modify=neigh_modify,
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


def make_writer_without_init(
    tmp_path,
    *,
    settings=None,
    sim_name="Test",
):
    writer = object.__new__(
        PreEqWriter
    )

    writer.settings = (
        settings
        if settings is not None
        else make_settings()
    )

    writer.out_dir = tmp_path
    writer.sim_name = sim_name

    return writer


def read_script(
    tmp_path,
    filename,
):
    return (
        tmp_path
        / "2_pre_equilibration"
        / filename
    ).read_text(
        encoding="utf-8"
    )


def command_lines(text):
    """
    Return tokenized non-comment, non-empty LAMMPS commands.

    This avoids fragile assertions based on cosmetic column spacing.
    """
    result = []

    for line in text.splitlines():
        stripped = line.strip()

        if not stripped:
            continue

        if stripped.startswith("#"):
            continue

        result.append(
            stripped.split()
        )

    return result


# =============================================================================
# Constructor
# =============================================================================


def test_constructor_stores_configuration(
    tmp_path,
    monkeypatch,
):
    settings = make_settings()

    simulation = make_simulation()

    calls = []

    monkeypatch.setattr(
        PreEqWriter,
        "write_pre_eq_file",
        lambda self, simulation:
        (
            calls.append(simulation)
            or "in.test_pre_equilibration"
        ),
    )

    writer = PreEqWriter(
        out_dir=tmp_path,
        settings=settings,
        simulation=simulation,
        sim_name="Polymer",
    )

    assert writer.settings is settings
    assert writer.out_dir == tmp_path
    assert writer.sim_name == "Polymer"

    assert calls == [
        simulation
    ]

    assert (
        writer.in_pre_eq_file_name
        == "in.test_pre_equilibration"
    )


# =============================================================================
# Basic file generation
# =============================================================================


def test_write_pre_eq_creates_stage_directory(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    writer.write_pre_eq_file(
        make_simulation()
    )

    assert (
        tmp_path
        / "2_pre_equilibration"
    ).is_dir()


def test_write_pre_eq_creates_input_file(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path,
        sim_name="Polymer",
    )

    filename = writer.write_pre_eq_file(
        make_simulation(
            tag="500K"
        )
    )

    assert filename == (
        "in.Polymer_500K_pre_equilibration"
    )

    assert (
        tmp_path
        / "2_pre_equilibration"
        / filename
    ).is_file()


def test_write_pre_eq_returns_only_filename(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_pre_eq_file(
        make_simulation()
    )

    assert isinstance(
        filename,
        str,
    )

    assert "/" not in filename
    assert "\\" not in filename


def test_script_header_contains_tag_and_autoreacter(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path,
        sim_name="Epoxy",
    )

    filename = writer.write_pre_eq_file(
        make_simulation(
            tag="373K"
        )
    )

    text = read_script(
        tmp_path,
        filename,
    )

    assert (
        "# Epoxy_373K "
        "Pre-Equilibration Script"
        in text
    )

    assert (
        "Generated"
        in text
    )

    assert (
        "by AutoREACTER"
        in text
    )


# =============================================================================
# LAMMPS settings
# =============================================================================


def test_script_contains_core_lammps_settings(
    tmp_path,
):
    settings = make_settings()

    writer = make_writer_without_init(
        tmp_path,
        settings=settings,
    )

    filename = writer.write_pre_eq_file(
        make_simulation()
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "units",
        "real",
    ] in commands

    assert [
        "dimension",
        "3",
    ] in commands

    assert [
        "boundary",
        "p",
        "p",
        "p",
    ] in commands

    assert [
        "atom_style",
        "full",
    ] in commands


def test_script_contains_force_field_styles(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_pre_eq_file(
        make_simulation()
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "bond_style",
        "class2",
    ] in commands

    assert [
        "angle_style",
        "class2",
    ] in commands

    assert [
        "dihedral_style",
        "class2",
    ] in commands

    assert [
        "improper_style",
        "class2",
    ] in commands


def test_script_contains_pair_settings(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_pre_eq_file(
        make_simulation()
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "pair_style",
        "lj/class2/coul/long",
        "12.0",
    ] in commands

    assert [
        "kspace_style",
        "pppm",
        "1.0e-4",
    ] in commands

    assert [
        "pair_modify",
        "mix",
        "sixthpower",
    ] in commands


def test_optional_neighbor_setting_is_written(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path,
        settings=make_settings(
            neighbor="2.0 bin",
        ),
    )

    filename = writer.write_pre_eq_file(
        make_simulation()
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "neighbor",
        "2.0",
        "bin",
    ] in commands


def test_optional_neigh_modify_setting_is_written(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path,
        settings=make_settings(
            neigh_modify=(
                "delay 0 every 1 "
                "check yes one 5000 "
                "page 100000"
            ),
        ),
    )

    filename = writer.write_pre_eq_file(
        make_simulation()
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "neigh_modify",
        "delay",
        "0",
        "every",
        "1",
        "check",
        "yes",
        "one",
        "5000",
        "page",
        "100000",
    ] in commands


def test_optional_neighbor_commands_are_omitted_when_none(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path,
        settings=make_settings(
            neighbor=None,
            neigh_modify=None,
        ),
    )

    filename = writer.write_pre_eq_file(
        make_simulation()
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    keywords = [
        command[0]
        for command in commands
    ]

    assert "neighbor" not in keywords
    assert "neigh_modify" not in keywords


# =============================================================================
# Input from densification
# =============================================================================


def test_script_reads_shrinked_box_from_previous_stage(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path,
        sim_name="Polymer",
    )

    filename = writer.write_pre_eq_file(
        make_simulation(
            tag="500K"
        )
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "read_data",
        "Polymer_500K_shrinked_box.data",
        "&",
    ] in commands


def test_read_data_includes_extra_topology_capacity(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_pre_eq_file(
        make_simulation()
    )

    text = read_script(
        tmp_path,
        filename,
    )

    assert (
        "extra/bond/per/atom 50"
        in text
    )

    assert (
        "extra/angle/per/atom 50"
        in text
    )

    assert (
        "extra/dihedral/per/atom 50"
        in text
    )

    assert (
        "extra/improper/per/atom 50"
        in text
    )

    assert (
        "extra/special/per/atom 50"
        in text
    )


# =============================================================================
# Minimization / base simulation setup
# =============================================================================


def test_script_contains_initial_minimization(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_pre_eq_file(
        make_simulation()
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "minimize",
        "0.0",
        "1.0e-8",
        "1000",
        "100000",
    ] in commands


def test_script_resets_timestep(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_pre_eq_file(
        make_simulation()
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "reset_timestep",
        "0",
    ] in commands


def test_script_uses_one_fs_timestep(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_pre_eq_file(
        make_simulation()
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "timestep",
        "1",
    ] in commands


def test_script_uses_thermo_interval_1000(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_pre_eq_file(
        make_simulation()
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "thermo",
        "1000",
    ] in commands


def test_script_contains_expected_thermo_style(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_pre_eq_file(
        make_simulation()
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "thermo_style",
        "custom",
        "step",
        "time",
        "temp",
        "press",
        "density",
        "vol",
        "pe",
        "ke",
        "etotal",
    ] in commands


# =============================================================================
# Dump / restart
# =============================================================================


def test_script_contains_expected_dump(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path,
        sim_name="Polymer",
    )

    filename = writer.write_pre_eq_file(
        make_simulation(
            tag="298K"
        )
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "dump",
        "dump_2",
        "all",
        "xyz",
        "1000",
        "Polymer_298K_pre_equilibration.xyz",
    ] in commands


def test_script_sets_dump_type_labels(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_pre_eq_file(
        make_simulation()
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "dump_modify",
        "dump_2",
        "types",
        "labels",
    ] in commands


def test_script_contains_two_restart_files(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path,
        sim_name="Polymer",
    )

    filename = writer.write_pre_eq_file(
        make_simulation(
            tag="500K"
        )
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "restart",
        "100",
        "Polymer_500K_pre_equilibration_backup1.restart",
        "Polymer_500K_pre_equilibration_backup2.restart",
    ] in commands


# =============================================================================
# Stage 1 - NVT temperature ramp
# =============================================================================


def test_stage_1_ramps_from_298_15_to_target_temperature(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_pre_eq_file(
        make_simulation(
            temperature=523.0
        )
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "fix",
        "nvt_1",
        "all",
        "nvt",
        "temp",
        "298.15",
        "523.0",
        "100.0",
    ] in commands


def test_stage_1_runs_50000_steps(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_pre_eq_file(
        make_simulation()
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    run_commands = [
        command
        for command in commands
        if command[0] == "run"
    ]

    assert run_commands[0] == [
        "run",
        "50000",
    ]


def test_stage_1_unfixes_nvt_1(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_pre_eq_file(
        make_simulation()
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "unfix",
        "nvt_1",
    ] in commands


# =============================================================================
# Current Stage 2 behavior
# =============================================================================


def test_stage_2_is_final_nvt_at_target_temperature(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_pre_eq_file(
        make_simulation(
            temperature=373.0
        )
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "fix",
        "nvt_3",
        "all",
        "nvt",
        "temp",
        "373.0",
        "373.0",
        "100.0",
    ] in commands


def test_stage_2_runs_50000_steps(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_pre_eq_file(
        make_simulation()
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    run_commands = [
        command
        for command in commands
        if command[0] == "run"
    ]

    assert run_commands == [
        [
            "run",
            "50000",
        ],
        [
            "run",
            "50000",
        ],
    ]


def test_stage_2_unfixes_nvt_3(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_pre_eq_file(
        make_simulation()
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "unfix",
        "nvt_3",
    ] in commands


def test_current_script_contains_no_active_npt_fix(
    tmp_path,
):
    """
    Characterization of CURRENT runtime behavior.

    The docstring/comments still refer to an NPT stage, but the NPT code
    itself is commented out in production. Generated scripts therefore
    contain no active NPT command.
    """
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_pre_eq_file(
        make_simulation()
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    active_npt_commands = [
        command
        for command in commands
        if "npt" in command
    ]

    assert active_npt_commands == []


def test_current_script_has_exactly_two_active_fix_commands(
    tmp_path,
):
    """
    Current equilibration workflow contains:
        nvt_1
        nvt_3

    No active npt_2.
    """
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_pre_eq_file(
        make_simulation()
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    fix_commands = [
        command
        for command in commands
        if command[0] == "fix"
    ]

    assert len(fix_commands) == 2

    assert fix_commands[0][1] == "nvt_1"
    assert fix_commands[1][1] == "nvt_3"


# =============================================================================
# Cleanup / output
# =============================================================================


def test_script_undumps_dump_2(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_pre_eq_file(
        make_simulation()
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "undump",
        "dump_2",
    ] in commands


def test_script_writes_pre_equilibrated_data(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path,
        sim_name="Polymer",
    )

    filename = writer.write_pre_eq_file(
        make_simulation(
            tag="500K"
        )
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "write_data",
        "Polymer_500K_pre_equilibrated.data",
    ] in commands


# =============================================================================
# Ordering / workflow characterization
# =============================================================================


def test_major_commands_are_in_expected_order(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_pre_eq_file(
        make_simulation()
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    keywords = [
        command[0]
        for command in commands
    ]

    read_index = keywords.index(
        "read_data"
    )

    minimize_index = keywords.index(
        "minimize"
    )

    first_fix_index = keywords.index(
        "fix"
    )

    first_run_index = keywords.index(
        "run"
    )

    write_index = keywords.index(
        "write_data"
    )

    assert (
        read_index
        < minimize_index
        < first_fix_index
        < first_run_index
        < write_index
    )


def test_two_nvt_stages_use_simulation_temperature(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_pre_eq_file(
        make_simulation(
            temperature=600.0
        )
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    nvt_1 = next(
        command
        for command in commands
        if (
            command[0] == "fix"
            and command[1] == "nvt_1"
        )
    )

    nvt_3 = next(
        command
        for command in commands
        if (
            command[0] == "fix"
            and command[1] == "nvt_3"
        )
    )

    assert nvt_1 == [
        "fix",
        "nvt_1",
        "all",
        "nvt",
        "temp",
        "298.15",
        "600.0",
        "100.0",
    ]

    assert nvt_3 == [
        "fix",
        "nvt_3",
        "all",
        "nvt",
        "temp",
        "600.0",
        "600.0",
        "100.0",
    ]