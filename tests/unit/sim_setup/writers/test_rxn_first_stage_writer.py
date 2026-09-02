from pathlib import Path
from types import SimpleNamespace

import pytest

import AutoREACTER.sim_setup.writers.rxn_first_stage_writer as rxn_module
from AutoREACTER.sim_setup.writers.rxn_first_stage_writer import (
    RxnFirstStageWriter,
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


def make_template(
    *,
    reaction_id=1,
    activity_stats=True,
    pre_reaction_file=None,
    post_reaction_file=None,
    map_file=None,
    map_file_with_delete_ids=None,
):
    return SimpleNamespace(
        reaction_id=reaction_id,
        activity_stats=activity_stats,
        pre_reaction_file=pre_reaction_file,
        post_reaction_file=post_reaction_file,
        map_file=map_file,
        map_file_with_delete_ids=map_file_with_delete_ids,
    )


def make_reacter_files(
    *,
    template_files=None,
):
    return SimpleNamespace(
        template_files=list(
            template_files or []
        ),
    )


def make_real_template(
    tmp_path,
    *,
    reaction_id=1,
    activity_stats=True,
    with_delete_map=False,
):
    pre = (
        tmp_path
        / f"template_pre_{reaction_id}.molecule"
    )

    post = (
        tmp_path
        / f"template_post_{reaction_id}.molecule"
    )

    standard_map = (
        tmp_path
        / f"RXN_{reaction_id}.map"
    )

    pre.write_text(
        "PRE",
        encoding="utf-8",
    )

    post.write_text(
        "POST",
        encoding="utf-8",
    )

    standard_map.write_text(
        "STANDARD MAP",
        encoding="utf-8",
    )

    delete_map = None

    if with_delete_map:
        delete_map = (
            tmp_path
            / (
                f"RXN_{reaction_id}"
                "_with_delete_ids.map"
            )
        )

        delete_map.write_text(
            "DELETE MAP",
            encoding="utf-8",
        )

    template = make_template(
        reaction_id=reaction_id,
        activity_stats=activity_stats,
        pre_reaction_file=pre,
        post_reaction_file=post,
        map_file=standard_map,
        map_file_with_delete_ids=delete_map,
    )

    return (
        template,
        pre,
        post,
        standard_map,
        delete_map,
    )


def make_writer_without_init(
    tmp_path,
    *,
    settings=None,
    reacter_files=None,
    sim_name="Test",
):
    writer = object.__new__(
        RxnFirstStageWriter
    )

    writer.settings = (
        settings
        if settings is not None
        else make_settings()
    )

    writer.out_dir = tmp_path
    writer.sim_name = sim_name

    writer.reacter_files = (
        reacter_files
        if reacter_files is not None
        else make_reacter_files()
    )

    return writer


def read_script(
    tmp_path,
    filename,
):
    return (
        tmp_path
        / "3_reaction_first_stage"
        / filename
    ).read_text(
        encoding="utf-8"
    )


def active_lines(text):
    """
    Return non-empty, non-comment lines.

    Important here because comments deliberately mention the optional
    _with_delete_ids.map filename.
    """
    result = []

    for line in text.splitlines():
        stripped = line.strip()

        if not stripped:
            continue

        if stripped.startswith("#"):
            continue

        result.append(
            stripped
        )

    return result


def command_tokens(text):
    return [
        line.split()
        for line in active_lines(text)
    ]


# =============================================================================
# Constructor
# =============================================================================


def test_constructor_stores_configuration(
    tmp_path,
    monkeypatch,
):
    settings = make_settings()

    reacter_files = (
        make_reacter_files()
    )

    simulation = (
        make_simulation()
    )

    calls = []

    monkeypatch.setattr(
        RxnFirstStageWriter,
        "write_first_stage_reaction_files",
        lambda self, simulation:
        (
            calls.append(simulation)
            or "in.test_reaction"
        ),
    )

    writer = RxnFirstStageWriter(
        out_dir=tmp_path,
        settings=settings,
        reacter_files=reacter_files,
        simulation=simulation,
        sim_name="Polymer",
    )

    assert writer.settings is settings

    assert (
        writer.reacter_files
        is reacter_files
    )

    assert writer.out_dir == tmp_path

    assert (
        writer.sim_name
        == "Polymer"
    )

    assert calls == [
        simulation
    ]

    assert (
        writer.first_stage_file_name
        == "in.test_reaction"
    )


# =============================================================================
# Basic file generation
# =============================================================================


def test_writer_creates_first_reaction_stage_directory(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    writer.write_first_stage_reaction_files(
        make_simulation()
    )

    assert (
        tmp_path
        / "3_reaction_first_stage"
    ).is_dir()


def test_writer_creates_expected_input_file(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path,
        sim_name="Polymer",
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation(
                tag="500K"
            )
        )
    )

    assert filename == (
        "in.Polymer_500K_reaction"
    )

    assert (
        tmp_path
        / "3_reaction_first_stage"
        / filename
    ).is_file()


def test_writer_returns_filename_only(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
    )

    assert isinstance(
        filename,
        str,
    )

    assert "/" not in filename
    assert "\\" not in filename


def test_script_header_contains_tag(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path,
        sim_name="Epoxy",
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation(
                tag="373K"
            )
        )
    )

    text = read_script(
        tmp_path,
        filename,
    )

    assert (
        "# Epoxy_373K "
        "First Reaction Stage Script"
        in text
    )

    assert (
        "by AutoREACTER"
        in text
    )


# =============================================================================
# LAMMPS settings
# =============================================================================


def test_script_contains_initialization_settings(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
    )

    commands = command_tokens(
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

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
    )

    commands = command_tokens(
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

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
    )

    commands = command_tokens(
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


def test_optional_neighbor_settings_are_written(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path,
        settings=make_settings(
            neighbor="2.0 bin",
            neigh_modify=(
                "delay 0 every 1 "
                "check yes"
            ),
        ),
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
    )

    commands = command_tokens(
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

    assert [
        "neigh_modify",
        "delay",
        "0",
        "every",
        "1",
        "check",
        "yes",
    ] in commands


def test_optional_neighbor_settings_omitted_when_none(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path,
        settings=make_settings(
            neighbor=None,
            neigh_modify=None,
        ),
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
    )

    commands = command_tokens(
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
# Input data
# =============================================================================


def test_reads_pre_equilibrated_data(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path,
        sim_name="Polymer",
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation(
                tag="500K"
            )
        )
    )

    commands = command_tokens(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "read_data",
        "Polymer_500K_pre_equilibrated.data",
        "&",
    ] in commands


def test_read_data_reserves_extra_topology_capacity(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
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
# Minimization / velocity / timestep
# =============================================================================


def test_script_contains_expected_minimization(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
    )

    commands = command_tokens(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "minimize",
        "1.0e-4",
        "1.0e-6",
        "1000",
        "10000",
    ] in commands


def test_velocity_uses_temperature_and_random_seed(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        rxn_module.random,
        "randint",
        lambda a, b: 123456,
    )

    writer = make_writer_without_init(
        tmp_path
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation(
                temperature=523.0
            )
        )
    )

    commands = command_tokens(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "velocity",
        "all",
        "create",
        "523.0",
        "123456",
        "dist",
        "gaussian",
    ] in commands


def test_script_uses_one_fs_timestep(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
    )

    commands = command_tokens(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "timestep",
        "1.0",
    ] in commands


def test_script_uses_thermo_interval_100(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
    )

    commands = command_tokens(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "thermo",
        "100",
    ] in commands


def test_script_resets_timestep(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
    )

    commands = command_tokens(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "reset_timestep",
        "0",
    ] in commands


# =============================================================================
# Reaction template definitions
# =============================================================================


def test_active_template_writes_pre_and_post_molecule_commands(
    tmp_path,
):
    (
        template,
        pre,
        post,
        _,
        _,
    ) = make_real_template(
        tmp_path,
        reaction_id=7,
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=[
                    template
                ]
            )
        ),
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
    )

    commands = command_tokens(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "molecule",
        "mol_pre_7",
        pre.name,
    ] in commands

    assert [
        "molecule",
        "mol_post_7",
        post.name,
    ] in commands


def test_inactive_template_is_not_written(
    tmp_path,
):
    (
        template,
        _,
        _,
        _,
        _,
    ) = make_real_template(
        tmp_path,
        reaction_id=9,
        activity_stats=False,
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=[
                    template
                ]
            )
        ),
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
    )

    text = read_script(
        tmp_path,
        filename,
    )

    assert "mol_pre_9" not in text
    assert "mol_post_9" not in text
    assert "rxn_stp_9" not in text


def test_missing_activity_stats_defaults_to_active(
    tmp_path,
):
    (
        template,
        _,
        _,
        _,
        _,
    ) = make_real_template(
        tmp_path,
        reaction_id=3,
    )

    del template.activity_stats

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=[
                    template
                ]
            )
        ),
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
    )

    text = read_script(
        tmp_path,
        filename,
    )

    assert "mol_pre_3" in text
    assert "rxn_stp_3" in text


def test_multiple_active_templates_are_written(
    tmp_path,
):
    template_1, *_ = (
        make_real_template(
            tmp_path,
            reaction_id=1,
        )
    )

    template_2, *_ = (
        make_real_template(
            tmp_path,
            reaction_id=2,
        )
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=[
                    template_1,
                    template_2,
                ]
            )
        ),
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
    )

    text = read_script(
        tmp_path,
        filename,
    )

    assert "rxn_stp_1" in text
    assert "rxn_stp_2" in text


# =============================================================================
# Standard map behavior
# =============================================================================


def test_standard_map_is_used_when_no_delete_map_exists(
    tmp_path,
):
    (
        template,
        _,
        _,
        standard_map,
        _,
    ) = make_real_template(
        tmp_path,
        reaction_id=4,
        with_delete_map=False,
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=[
                    template
                ]
            )
        ),
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
    )

    lines = active_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    react_line = next(
        line
        for line in lines
        if "rxn_stp_4" in line
    )

    assert (
        standard_map.name
        in react_line
    )


def test_standard_map_remains_default_even_when_delete_map_exists(
    tmp_path,
):
    """
    Core AutoREACTER behavior.

    The supplementary DeleteIDs map is supplied to the user, but the generated
    LAMMPS script MUST continue to use RXN_N.map by default.
    """
    (
        template,
        _,
        _,
        standard_map,
        delete_map,
    ) = make_real_template(
        tmp_path,
        reaction_id=5,
        with_delete_map=True,
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=[
                    template
                ]
            )
        ),
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
    )

    lines = active_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    react_line = next(
        line
        for line in lines
        if "rxn_stp_5" in line
    )

    assert (
        standard_map.name
        in react_line
    )

    assert (
        delete_map.name
        not in react_line
    )


def test_delete_map_is_never_automatically_selected(
    tmp_path,
):
    (
        template,
        _,
        _,
        _,
        delete_map,
    ) = make_real_template(
        tmp_path,
        reaction_id=8,
        with_delete_map=True,
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=[
                    template
                ]
            )
        ),
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
    )

    active = "\n".join(
        active_lines(
            read_script(
                tmp_path,
                filename,
            )
        )
    )

    assert (
        delete_map.name
        not in active
    )


def test_reaction_command_contains_expected_options(
    tmp_path,
):
    (
        template,
        _,
        _,
        standard_map,
        _,
    ) = make_real_template(
        tmp_path,
        reaction_id=2,
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=[
                    template
                ]
            )
        ),
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
    )

    lines = active_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    react_line = next(
        line
        for line in lines
        if "rxn_stp_2" in line
    )

    tokens = react_line.split()

    assert tokens[:2] == [
        "react",
        "rxn_stp_2",
    ]

    assert "all" in tokens
    assert "1" in tokens
    assert "0.0" in tokens
    assert "3.5" in tokens

    assert (
        "mol_pre_2"
        in tokens
    )

    assert (
        "mol_post_2"
        in tokens
    )

    assert (
        standard_map.name
        in tokens
    )

    assert tokens[-4:] == [
        "stabilize_steps",
        "60",
        "rescale_charges",
        "yes",
    ]


# =============================================================================
# fix bond/react
# =============================================================================


def test_fix_bond_react_header_is_written(
    tmp_path,
):
    template, *_ = (
        make_real_template(
            tmp_path,
            reaction_id=1,
        )
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=[
                    template
                ]
            )
        ),
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
    )

    text = read_script(
        tmp_path,
        filename,
    )

    assert (
        "rxns all bond/react "
        "stabilization yes "
        "statted_grp 0.03"
        in text
    )


def test_multiple_reactions_are_joined_into_same_fix(
    tmp_path,
):
    template_1, *_ = (
        make_real_template(
            tmp_path,
            reaction_id=1,
        )
    )

    template_2, *_ = (
        make_real_template(
            tmp_path,
            reaction_id=2,
        )
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=[
                    template_1,
                    template_2,
                ]
            )
        ),
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
    )

    text = read_script(
        tmp_path,
        filename,
    )

    assert "rxn_stp_1" in text
    assert "rxn_stp_2" in text

    assert (
        " & "
        in text
    )


# =============================================================================
# Thermostat / output
# =============================================================================


def test_active_nvt_uses_simulation_temperature(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation(
                temperature=373.0
            )
        )
    )

    commands = command_tokens(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "fix",
        "1",
        "statted_grp_REACT",
        "nvt",
        "temp",
        "373.0",
        "373.0",
        "100.0",
    ] in commands


def test_npt_command_is_currently_commented_out(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
    )

    active = active_lines(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert not any(
        " npt " in f" {line} "
        for line in active
    )


def test_thermo_style_contains_reaction_output(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
    )

    commands = command_tokens(
        read_script(
            tmp_path,
            filename,
        )
    )

    thermo_style = next(
        command
        for command in commands
        if command[0] == "thermo_style"
    )

    assert "custom" in thermo_style
    assert "step" in thermo_style
    assert "time" in thermo_style
    assert "temp" in thermo_style
    assert "f_rxns[*]" in thermo_style
    assert "press" in thermo_style
    assert "density" in thermo_style
    assert "vol" in thermo_style
    assert "pe" in thermo_style
    assert "ke" in thermo_style
    assert "etotal" in thermo_style


def test_reaction_stage_runs_one_million_steps(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
    )

    commands = command_tokens(
        read_script(
            tmp_path,
            filename,
        )
    )

    assert [
        "run",
        "1000000",
    ] in commands


def test_output_names_use_first_stage_range(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path,
        sim_name="Polymer",
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation(
                tag="500K"
            )
        )
    )

    text = read_script(
        tmp_path,
        filename,
    )

    output_base = (
        "Polymer_500K"
        "_reacted_0M-1M_3.5A"
    )

    assert (
        f"{output_base}.xyz"
        in text
    )

    assert (
        f"{output_base}_backup1.restart"
        in text
    )

    assert (
        f"{output_base}_backup2.restart"
        in text
    )

    assert (
        f"{output_base}.restart"
        in text
    )

    assert (
        f"{output_base}.data"
        in text
    )


def test_write_data_uses_nofix(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation()
        )
    )

    commands = command_tokens(
        read_script(
            tmp_path,
            filename,
        )
    )

    write_commands = [
        command
        for command in commands
        if command[0] == "write_data"
    ]

    assert len(
        write_commands
    ) == 1

    assert (
        write_commands[0][-1]
        == "nofix"
    )


# =============================================================================
# _copy_required_files
# =============================================================================


def test_copy_required_files_copies_standard_map_and_templates(
    tmp_path,
):
    (
        template,
        pre,
        post,
        standard_map,
        _,
    ) = make_real_template(
        tmp_path,
        reaction_id=1,
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=[
                    template
                ]
            )
        ),
    )

    dest = (
        tmp_path / "reaction"
    )

    dest.mkdir()

    writer._copy_required_files(
        dest
    )

    assert (
        dest / standard_map.name
    ).is_file()

    assert (
        dest / pre.name
    ).is_file()

    assert (
        dest / post.name
    ).is_file()


def test_copy_required_files_copies_optional_delete_map_when_present(
    tmp_path,
):
    """
    Supplementary DeleteIDs map is copied for the user.

    It is NOT automatically used by the generated LAMMPS script.
    """
    (
        template,
        _,
        _,
        standard_map,
        delete_map,
    ) = make_real_template(
        tmp_path,
        reaction_id=2,
        with_delete_map=True,
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=[
                    template
                ]
            )
        ),
    )

    dest = (
        tmp_path / "reaction"
    )

    dest.mkdir()

    writer._copy_required_files(
        dest
    )

    assert (
        dest / standard_map.name
    ).is_file()

    assert (
        dest / delete_map.name
    ).is_file()


def test_copy_required_files_without_delete_map_still_succeeds(
    tmp_path,
):
    (
        template,
        pre,
        post,
        standard_map,
        delete_map,
    ) = make_real_template(
        tmp_path,
        reaction_id=3,
        with_delete_map=False,
    )

    assert delete_map is None

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=[
                    template
                ]
            )
        ),
    )

    dest = (
        tmp_path / "reaction"
    )

    dest.mkdir()

    writer._copy_required_files(
        dest
    )

    assert (
        dest / standard_map.name
    ).is_file()

    assert (
        dest / pre.name
    ).is_file()

    assert (
        dest / post.name
    ).is_file()


def test_copy_required_files_missing_delete_map_attribute_still_succeeds(
    tmp_path,
):
    (
        template,
        _,
        _,
        standard_map,
        _,
    ) = make_real_template(
        tmp_path,
        reaction_id=4,
    )

    del (
        template
        .map_file_with_delete_ids
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=[
                    template
                ]
            )
        ),
    )

    dest = (
        tmp_path / "reaction"
    )

    dest.mkdir()

    writer._copy_required_files(
        dest
    )

    assert (
        dest / standard_map.name
    ).is_file()


def test_copy_required_files_missing_standard_map_raises(
    tmp_path,
):
    pre = (
        tmp_path / "pre.molecule"
    )

    post = (
        tmp_path / "post.molecule"
    )

    pre.write_text(
        "PRE",
        encoding="utf-8",
    )

    post.write_text(
        "POST",
        encoding="utf-8",
    )

    missing_map = (
        tmp_path / "missing.map"
    )

    template = make_template(
        reaction_id=1,
        pre_reaction_file=pre,
        post_reaction_file=post,
        map_file=missing_map,
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=[
                    template
                ]
            )
        ),
    )

    dest = (
        tmp_path / "dest"
    )

    dest.mkdir()

    with pytest.raises(
        FileNotFoundError,
        match=(
            "Required reaction file "
            "not found"
        ),
    ):
        writer._copy_required_files(
            dest
        )


def test_copy_required_files_missing_pre_template_raises(
    tmp_path,
):
    post = (
        tmp_path / "post.molecule"
    )

    map_file = (
        tmp_path / "RXN_1.map"
    )

    post.write_text(
        "POST",
        encoding="utf-8",
    )

    map_file.write_text(
        "MAP",
        encoding="utf-8",
    )

    template = make_template(
        reaction_id=1,
        pre_reaction_file=(
            tmp_path
            / "missing_pre.molecule"
        ),
        post_reaction_file=post,
        map_file=map_file,
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=[
                    template
                ]
            )
        ),
    )

    dest = (
        tmp_path / "dest"
    )

    dest.mkdir()

    with pytest.raises(
        FileNotFoundError,
    ):
        writer._copy_required_files(
            dest
        )


def test_copy_required_files_missing_post_template_raises(
    tmp_path,
):
    pre = (
        tmp_path / "pre.molecule"
    )

    map_file = (
        tmp_path / "RXN_1.map"
    )

    pre.write_text(
        "PRE",
        encoding="utf-8",
    )

    map_file.write_text(
        "MAP",
        encoding="utf-8",
    )

    template = make_template(
        reaction_id=1,
        pre_reaction_file=pre,
        post_reaction_file=(
            tmp_path
            / "missing_post.molecule"
        ),
        map_file=map_file,
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=[
                    template
                ]
            )
        ),
    )

    dest = (
        tmp_path / "dest"
    )

    dest.mkdir()

    with pytest.raises(
        FileNotFoundError,
    ):
        writer._copy_required_files(
            dest
        )


def test_copy_required_files_missing_optional_delete_map_path_raises_when_declared(
    tmp_path,
):
    """
    Current behavior:

    The delete map is optional in the sense that the attribute may be None.

    But if metadata explicitly points to a delete map, that declared file
    must actually exist for it to be copied.
    """
    (
        template,
        _,
        _,
        _,
        _,
    ) = make_real_template(
        tmp_path,
        reaction_id=6,
    )

    template.map_file_with_delete_ids = (
        tmp_path
        / "RXN_6_with_delete_ids.map"
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=[
                    template
                ]
            )
        ),
    )

    dest = (
        tmp_path / "dest"
    )

    dest.mkdir()

    with pytest.raises(
        FileNotFoundError,
        match=(
            "Required reaction file "
            "not found"
        ),
    ):
        writer._copy_required_files(
            dest
        )


def test_copy_required_files_ignores_inactive_templates(
    tmp_path,
):
    template = make_template(
        reaction_id=99,
        activity_stats=False,
        pre_reaction_file=(
            tmp_path
            / "missing_pre.molecule"
        ),
        post_reaction_file=(
            tmp_path
            / "missing_post.molecule"
        ),
        map_file=(
            tmp_path
            / "missing.map"
        ),
        map_file_with_delete_ids=(
            tmp_path
            / "missing_delete.map"
        ),
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=[
                    template
                ]
            )
        ),
    )

    dest = (
        tmp_path / "dest"
    )

    dest.mkdir()

    writer._copy_required_files(
        dest
    )

    assert list(
        dest.iterdir()
    ) == []


# =============================================================================
# Full delete-map integration
# =============================================================================


def test_delete_map_is_supplied_but_standard_map_remains_in_generated_script(
    tmp_path,
):
    """
    End-to-end characterization of the intended map behavior.

    Both maps are copied:
        RXN_10.map
        RXN_10_with_delete_ids.map

    But the generated fix bond/react command references:
        RXN_10.map
    """
    (
        template,
        pre,
        post,
        standard_map,
        delete_map,
    ) = make_real_template(
        tmp_path,
        reaction_id=10,
        with_delete_map=True,
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=[
                    template
                ]
            )
        ),
        sim_name="Polymer",
    )

    filename = (
        writer
        .write_first_stage_reaction_files(
            make_simulation(
                tag="500K"
            )
        )
    )

    stage_dir = (
        tmp_path
        / "3_reaction_first_stage"
    )

    assert (
        stage_dir / pre.name
    ).is_file()

    assert (
        stage_dir / post.name
    ).is_file()

    assert (
        stage_dir / standard_map.name
    ).is_file()

    assert (
        stage_dir / delete_map.name
    ).is_file()

    text = read_script(
        tmp_path,
        filename,
    )

    active = "\n".join(
        active_lines(text)
    )

    assert (
        standard_map.name
        in active
    )

    assert (
        delete_map.name
        not in active
    )