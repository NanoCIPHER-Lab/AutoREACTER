from types import SimpleNamespace

import pytest

from AutoREACTER.sim_setup.writers.post_eq_writer import (
    PostEqWriter,
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
        PostEqWriter
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
        / "5_post_equilibration"
        / filename
    ).read_text(
        encoding="utf-8"
    )


def command_lines(text):
    """
    Return tokenized active LAMMPS commands.

    Comments and blank lines are ignored so tests do not depend on
    cosmetic alignment spacing.
    """
    commands = []

    for line in text.splitlines():
        stripped = line.strip()

        if not stripped:
            continue

        if stripped.startswith("#"):
            continue

        commands.append(
            stripped.split()
        )

    return commands


# =============================================================================
# Constructor
# =============================================================================


def test_constructor_stores_configuration_and_writes_file(
    tmp_path,
    monkeypatch,
):
    settings = make_settings()

    simulation = make_simulation()

    calls = []

    def fake_write(
        self,
        *,
        simulation,
        write_second_reaction_stage=True,
    ):
        calls.append(
            (
                simulation,
                write_second_reaction_stage,
            )
        )

        return "in.test_post_equilibration"

    monkeypatch.setattr(
        PostEqWriter,
        "write_post_eq_file",
        fake_write,
    )

    writer = PostEqWriter(
        out_dir=tmp_path,
        settings=settings,
        simulation=simulation,
        sim_name="Polymer",
    )

    assert writer.settings is settings
    assert writer.out_dir == tmp_path
    assert writer.sim_name == "Polymer"

    assert calls == [
        (
            simulation,
            True,
        )
    ]


def test_constructor_forwards_false_second_stage_flag(
    tmp_path,
    monkeypatch,
):
    calls = []

    def fake_write(
        self,
        *,
        simulation,
        write_second_reaction_stage=True,
    ):
        calls.append(
            write_second_reaction_stage
        )

        return "in.test"

    monkeypatch.setattr(
        PostEqWriter,
        "write_post_eq_file",
        fake_write,
    )

    PostEqWriter(
        out_dir=tmp_path,
        settings=make_settings(),
        simulation=make_simulation(),
        sim_name="Polymer",
        write_second_reaction_stage=False,
    )

    assert calls == [
        False
    ]


# =============================================================================
# File generation
# =============================================================================


def test_write_post_eq_creates_stage_directory(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    writer.write_post_eq_file(
        make_simulation()
    )

    assert (
        tmp_path
        / "5_post_equilibration"
    ).is_dir()


def test_write_post_eq_creates_expected_input_file(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path,
        sim_name="Polymer",
    )

    filename = writer.write_post_eq_file(
        make_simulation(
            tag="500K"
        )
    )

    assert filename == (
        "in.Polymer_500K_post_equilibration"
    )

    assert (
        tmp_path
        / "5_post_equilibration"
        / filename
    ).is_file()


def test_write_post_eq_returns_only_filename(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_post_eq_file(
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

    filename = writer.write_post_eq_file(
        make_simulation(
            tag="373K"
        )
    )

    text = read_script(
        tmp_path,
        filename,
    )

    assert (
        "#Epoxy_373K "
        "Post-Equilibration Script"
        in text
    )

    assert "Generated" in text

    assert (
        "by AutoREACTER"
        in text
    )


# =============================================================================
# Input reacted data selection
# =============================================================================


def test_default_reads_second_reaction_stage_output(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path,
        sim_name="Polymer",
    )

    filename = writer.write_post_eq_file(
        make_simulation(
            tag="500K"
        ),
        write_second_reaction_stage=True,
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "read_data",
        "Polymer_500K_reacted_1M-3.5_5.0A.data",
        "&",
    ] in commands


def test_false_second_stage_reads_first_reaction_output(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path,
        sim_name="Polymer",
    )

    filename = writer.write_post_eq_file(
        make_simulation(
            tag="500K"
        ),
        write_second_reaction_stage=False,
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "read_data",
        "Polymer_500K_reacted_0M-1M_3.5A.data",
        "&",
    ] in commands


def test_read_data_includes_extra_topology_capacity(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_post_eq_file(
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
# LAMMPS settings
# =============================================================================


def test_script_contains_core_lammps_settings(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_post_eq_file(
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

    filename = writer.write_post_eq_file(
        make_simulation()
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "angle_style",
        "class2",
    ] in commands

    assert [
        "bond_style",
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

    filename = writer.write_post_eq_file(
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


def test_neighbor_is_written_when_defined(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path,
        settings=make_settings(
            neighbor="2.0 bin",
        ),
    )

    filename = writer.write_post_eq_file(
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


def test_neigh_modify_is_written_when_defined(
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

    filename = writer.write_post_eq_file(
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


def test_neighbor_commands_omitted_when_none(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path,
        settings=make_settings(
            neighbor=None,
            neigh_modify=None,
        ),
    )

    filename = writer.write_post_eq_file(
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
# Setup / minimization
# =============================================================================


def test_script_contains_minimization(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_post_eq_file(
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

    filename = writer.write_post_eq_file(
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

    filename = writer.write_post_eq_file(
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


def test_script_uses_thermo_1000(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_post_eq_file(
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

    filename = writer.write_post_eq_file(
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

    filename = writer.write_post_eq_file(
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
        "dump",
        "dump_1",
        "all",
        "xyz",
        "100",
        "Polymer_500K_post_equilibration.xyz",
    ] in commands


def test_script_sets_dump_type_labels(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_post_eq_file(
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
        "dump_1",
        "types",
        "labels",
    ] in commands


def test_script_contains_rotating_restart_files(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path,
        sim_name="Polymer",
    )

    filename = writer.write_post_eq_file(
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
        "Polymer_500K_post_equilibration_backup1.restart",
        "Polymer_500K_post_equilibration_backup2.restart",
    ] in commands


# =============================================================================
# Stage 1 - NVT cool to 300 K
# =============================================================================


def test_stage_1_nvt_cools_from_sim_temperature_to_300(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_post_eq_file(
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
        "523.0",
        "300.0",
        "100.0",
    ] in commands


def test_stage_1_runs_100000_steps(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_post_eq_file(
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
        "100000",
    ]


def test_stage_1_unfixes_nvt_1(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_post_eq_file(
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
# Stage 2 - NPT
# =============================================================================


def test_stage_2_uses_npt_at_300_k_and_zero_pressure(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_post_eq_file(
        make_simulation()
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "fix",
        "npt_2",
        "all",
        "npt",
        "temp",
        "300.0",
        "300.0",
        "100.0",
        "iso",
        "0.0",
        "0.0",
        "1000.0",
    ] in commands


def test_stage_2_unfixes_npt_2(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_post_eq_file(
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
        "npt_2",
    ] in commands


# =============================================================================
# Stage 3 - final NVT
# =============================================================================


def test_stage_3_final_nvt_is_at_300_k(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_post_eq_file(
        make_simulation()
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
        "300.0",
        "300.0",
        "100.0",
    ] in commands


def test_stage_3_unfixes_nvt_3(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_post_eq_file(
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


def test_exactly_three_run_commands_each_100000_steps(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_post_eq_file(
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
            "100000",
        ],
        [
            "run",
            "100000",
        ],
        [
            "run",
            "100000",
        ],
    ]


def test_exactly_three_active_fix_commands(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_post_eq_file(
        make_simulation()
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    fixes = [
        command
        for command in commands
        if command[0] == "fix"
    ]

    assert len(fixes) == 3

    assert fixes[0][1] == "nvt_1"
    assert fixes[1][1] == "npt_2"
    assert fixes[2][1] == "nvt_3"


# =============================================================================
# Output / cleanup
# =============================================================================


def test_script_undumps_trajectory(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_post_eq_file(
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
        "dump_1",
    ] in commands


def test_script_writes_post_equilibrated_data(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path,
        sim_name="Polymer",
    )

    filename = writer.write_post_eq_file(
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
        "Polymer_500K_post_equilibrated.data",
    ] in commands


# =============================================================================
# Workflow order
# =============================================================================


def test_major_commands_are_in_expected_order(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_post_eq_file(
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


def test_fix_and_unfix_order_is_preserved(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_post_eq_file(
        make_simulation()
    )

    commands = command_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    fix_sequence = [
        command[:2]
        for command in commands
        if command[0] in {
            "fix",
            "unfix",
        }
    ]

    assert fix_sequence == [
        [
            "fix",
            "nvt_1",
        ],
        [
            "unfix",
            "nvt_1",
        ],
        [
            "fix",
            "npt_2",
        ],
        [
            "unfix",
            "npt_2",
        ],
        [
            "fix",
            "nvt_3",
        ],
        [
            "unfix",
            "nvt_3",
        ],
    ]


# =============================================================================
# Reaction-map separation
# =============================================================================


def test_post_equilibration_script_contains_no_reaction_map_commands(
    tmp_path,
):
    """
    Post-equilibration consumes an already-reacted data file.

    Standard RXN_N.map files and optional RXN_N_with_delete_ids.map files
    are reaction-stage inputs and must not be referenced here.
    """
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_post_eq_file(
        make_simulation()
    )

    text = read_script(
        tmp_path,
        filename,
    )

    assert "RXN_" not in text

    assert ".map" not in text


def test_post_equilibration_contains_no_bond_react_fix(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = writer.write_post_eq_file(
        make_simulation()
    )

    text = read_script(
        tmp_path,
        filename,
    )

    assert "bond/react" not in text
