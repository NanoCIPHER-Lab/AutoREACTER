from types import SimpleNamespace

import pytest

import AutoREACTER.sim_setup.writers.rxn_second_stage_writer as rxn_module
from AutoREACTER.sim_setup.writers.rxn_second_stage_writer import (
    RxnSecondStageWriter,
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
        RxnSecondStageWriter
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
        / "4_reaction_second_stage"
        / filename
    ).read_text(
        encoding="utf-8"
    )


def active_lines(text):
    """
    Return non-empty, non-comment LAMMPS lines.

    Comments are ignored because the generated file intentionally
    mentions RXN_i_with_delete_ids.map in explanatory comments.
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
        RxnSecondStageWriter,
        "write_second_stage_reaction_files",
        lambda self, simulation:
        (
            calls.append(simulation)
            or "in.test_stage_2"
        ),
    )

    writer = RxnSecondStageWriter(
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
        writer.second_stage_file_name
        == "in.test_stage_2"
    )


# =============================================================================
# Active templates
# =============================================================================


def test_no_active_templates_raises(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=[]
            )
        ),
    )

    with pytest.raises(
        ValueError,
        match=(
            "No active reaction templates "
            "available for the second reaction stage"
        ),
    ):
        writer.write_second_stage_reaction_files(
            make_simulation()
        )


def test_all_inactive_templates_raise(
    tmp_path,
):
    template = make_template(
        reaction_id=1,
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

    with pytest.raises(
        ValueError,
        match="No active reaction templates",
    ):
        writer.write_second_stage_reaction_files(
            make_simulation()
        )


def test_missing_activity_stats_defaults_to_active(
    tmp_path,
):
    template, *_ = (
        make_real_template(
            tmp_path,
            reaction_id=3,
        )
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
        .write_second_stage_reaction_files(
            make_simulation()
        )
    )

    text = read_script(
        tmp_path,
        filename,
    )

    assert "rxn_stp_3" in text


# =============================================================================
# File generation
# =============================================================================


def test_writer_creates_second_stage_directory(
    tmp_path,
):
    template, *_ = (
        make_real_template(
            tmp_path
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

    writer.write_second_stage_reaction_files(
        make_simulation()
    )

    assert (
        tmp_path
        / "4_reaction_second_stage"
    ).is_dir()


def test_writer_creates_expected_input_file(
    tmp_path,
):
    template, *_ = (
        make_real_template(
            tmp_path
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
        sim_name="Polymer",
    )

    filename = (
        writer
        .write_second_stage_reaction_files(
            make_simulation(
                tag="500K"
            )
        )
    )

    assert filename == (
        "in.Polymer_500K_reaction_stage_2"
    )

    assert (
        tmp_path
        / "4_reaction_second_stage"
        / filename
    ).is_file()


def test_writer_returns_filename_only(
    tmp_path,
):
    template, *_ = (
        make_real_template(
            tmp_path
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
        .write_second_stage_reaction_files(
            make_simulation()
        )
    )

    assert isinstance(
        filename,
        str,
    )

    assert "/" not in filename
    assert "\\" not in filename


def test_header_contains_tag_and_autoreacter(
    tmp_path,
):
    template, *_ = (
        make_real_template(
            tmp_path
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
        sim_name="Epoxy",
    )

    filename = (
        writer
        .write_second_stage_reaction_files(
            make_simulation(
                tag="500K"
            )
        )
    )

    text = read_script(
        tmp_path,
        filename,
    )

    assert (
        "# Epoxy_500K "
        "Second Reaction Stage Script"
        in text
    )

    assert (
        "by AutoREACTER"
        in text
    )


# =============================================================================
# LAMMPS initialization
# =============================================================================


def test_script_contains_initialization_settings(
    tmp_path,
):
    template, *_ = (
        make_real_template(
            tmp_path
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
        .write_second_stage_reaction_files(
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
    template, *_ = (
        make_real_template(
            tmp_path
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
        .write_second_stage_reaction_files(
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
    template, *_ = (
        make_real_template(
            tmp_path
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
        .write_second_stage_reaction_files(
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
    template, *_ = (
        make_real_template(
            tmp_path
        )
    )

    writer = make_writer_without_init(
        tmp_path,
        settings=make_settings(
            neighbor="2.0 bin",
            neigh_modify=(
                "delay 0 every 1 "
                "check yes"
            ),
        ),
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
        .write_second_stage_reaction_files(
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
    template, *_ = (
        make_real_template(
            tmp_path
        )
    )

    writer = make_writer_without_init(
        tmp_path,
        settings=make_settings(
            neighbor=None,
            neigh_modify=None,
        ),
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
        .write_second_stage_reaction_files(
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
# Input structure
# =============================================================================


def test_reads_first_stage_reacted_data(
    tmp_path,
):
    template, *_ = (
        make_real_template(
            tmp_path
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
        sim_name="Polymer",
    )

    filename = (
        writer
        .write_second_stage_reaction_files(
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
        "Polymer_500K_reacted_0M-1M_3.5A.data",
        "&",
    ] in commands


def test_read_data_reserves_extra_topology_capacity(
    tmp_path,
):
    template, *_ = (
        make_real_template(
            tmp_path
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
        .write_second_stage_reaction_files(
            make_simulation()
        )
    )

    text = read_script(
        tmp_path,
        filename,
    )

    assert "extra/bond/per/atom 50" in text
    assert "extra/angle/per/atom 50" in text
    assert "extra/dihedral/per/atom 50" in text
    assert "extra/improper/per/atom 50" in text
    assert "extra/special/per/atom 50" in text


# =============================================================================
# Minimization / velocity / timestep
# =============================================================================


def test_script_contains_minimization(
    tmp_path,
):
    template, *_ = (
        make_real_template(
            tmp_path
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
        .write_second_stage_reaction_files(
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


def test_velocity_uses_simulation_temperature_and_seed(
    tmp_path,
    monkeypatch,
):
    template, *_ = (
        make_real_template(
            tmp_path
        )
    )

    monkeypatch.setattr(
        rxn_module.random,
        "randint",
        lambda a, b: 123456,
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
        .write_second_stage_reaction_files(
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


def test_second_stage_timestep_thermo_and_reset(
    tmp_path,
):
    template, *_ = (
        make_real_template(
            tmp_path
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
        .write_second_stage_reaction_files(
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

    assert [
        "thermo",
        "100",
    ] in commands

    assert [
        "reset_timestep",
        "1000000",
    ] in commands


# =============================================================================
# Reaction definitions
# =============================================================================


def test_active_template_writes_pre_and_post_molecules(
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
        .write_second_stage_reaction_files(
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


def test_inactive_template_not_written_when_active_template_exists(
    tmp_path,
):
    active, *_ = (
        make_real_template(
            tmp_path,
            reaction_id=1,
        )
    )

    inactive, *_ = (
        make_real_template(
            tmp_path,
            reaction_id=2,
            activity_stats=False,
        )
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=[
                    active,
                    inactive,
                ]
            )
        ),
    )

    filename = (
        writer
        .write_second_stage_reaction_files(
            make_simulation()
        )
    )

    text = read_script(
        tmp_path,
        filename,
    )

    assert "rxn_stp_1" in text
    assert "rxn_stp_2" not in text


def test_reaction_uses_five_angstrom_cutoff(
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
        .write_second_stage_reaction_files(
            make_simulation()
        )
    )

    react_line = next(
        line
        for line in active_lines(
            read_script(
                tmp_path,
                filename,
            )
        )
        if "rxn_stp_4" in line
    )

    tokens = (
        react_line.split()
    )

    assert "0.0" in tokens
    assert "5.0" in tokens

    assert (
        standard_map.name
        in tokens
    )


def test_reaction_uses_stabilize_steps_200(
    tmp_path,
):
    template, *_ = (
        make_real_template(
            tmp_path,
            reaction_id=5,
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
        .write_second_stage_reaction_files(
            make_simulation()
        )
    )

    react_line = next(
        line
        for line in active_lines(
            read_script(
                tmp_path,
                filename,
            )
        )
        if "rxn_stp_5" in line
    )

    assert (
        react_line.split()[-4:]
        == [
            "stabilize_steps",
            "200",
            "rescale_charges",
            "yes",
        ]
    )


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
        reaction_id=6,
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
        .write_second_stage_reaction_files(
            make_simulation()
        )
    )

    react_line = next(
        line
        for line in active_lines(
            read_script(
                tmp_path,
                filename,
            )
        )
        if "rxn_stp_6" in line
    )

    assert (
        standard_map.name
        in react_line
    )


def test_standard_map_remains_default_when_delete_map_exists(
    tmp_path,
):
    """
    AutoREACTER supplies the optional delete-ID map to the user,
    but the generated LAMMPS command must still use RXN_N.map.
    """
    (
        template,
        _,
        _,
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
    )

    filename = (
        writer
        .write_second_stage_reaction_files(
            make_simulation()
        )
    )

    react_line = next(
        line
        for line in active_lines(
            read_script(
                tmp_path,
                filename,
            )
        )
        if "rxn_stp_10" in line
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
        reaction_id=12,
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
        .write_second_stage_reaction_files(
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


# =============================================================================
# fix bond/react thermo vector
# =============================================================================


def test_single_reaction_uses_explicit_f_rxns_1(
    tmp_path,
):
    template, *_ = (
        make_real_template(
            tmp_path,
            reaction_id=17,
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
        .write_second_stage_reaction_files(
            make_simulation()
        )
    )

    commands = command_tokens(
        read_script(
            tmp_path,
            filename,
        )
    )

    thermo = next(
        command
        for command in commands
        if command[0] == "thermo_style"
    )

    assert (
        "f_rxns[1]"
        in thermo
    )

    assert (
        "f_rxns[*]"
        not in thermo
    )


def test_multiple_reactions_use_contiguous_thermo_indices(
    tmp_path,
):
    templates = [
        make_real_template(
            tmp_path,
            reaction_id=reaction_id,
        )[0]
        for reaction_id in (
            4,
            17,
            68,
        )
    ]

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=templates
            )
        ),
    )

    filename = (
        writer
        .write_second_stage_reaction_files(
            make_simulation()
        )
    )

    commands = command_tokens(
        read_script(
            tmp_path,
            filename,
        )
    )

    thermo = next(
        command
        for command in commands
        if command[0] == "thermo_style"
    )

    assert "f_rxns[1]" in thermo
    assert "f_rxns[2]" in thermo
    assert "f_rxns[3]" in thermo

    assert "f_rxns[4]" not in thermo


def test_thermo_indices_do_not_use_reaction_ids(
    tmp_path,
):
    template_17, *_ = (
        make_real_template(
            tmp_path,
            reaction_id=17,
        )
    )

    template_68, *_ = (
        make_real_template(
            tmp_path,
            reaction_id=68,
        )
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=[
                    template_17,
                    template_68,
                ]
            )
        ),
    )

    filename = (
        writer
        .write_second_stage_reaction_files(
            make_simulation()
        )
    )

    text = read_script(
        tmp_path,
        filename,
    )

    assert "f_rxns[1]" in text
    assert "f_rxns[2]" in text

    assert "f_rxns[17]" not in text
    assert "f_rxns[68]" not in text


# =============================================================================
# Thermostat / run
# =============================================================================


def test_active_nvt_uses_simulation_temperature(
    tmp_path,
):
    template, *_ = (
        make_real_template(
            tmp_path
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
        .write_second_stage_reaction_files(
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


def test_npt_command_remains_commented_out(
    tmp_path,
):
    template, *_ = (
        make_real_template(
            tmp_path
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
        .write_second_stage_reaction_files(
            make_simulation()
        )
    )

    assert not any(
        " npt " in f" {line} "
        for line in active_lines(
            read_script(
                tmp_path,
                filename,
            )
        )
    )


def test_second_stage_runs_2500000_steps(
    tmp_path,
):
    template, *_ = (
        make_real_template(
            tmp_path
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
        .write_second_stage_reaction_files(
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
        "2500000",
    ] in commands


# =============================================================================
# Output naming
# =============================================================================


def test_output_names_use_second_stage_range(
    tmp_path,
):
    template, *_ = (
        make_real_template(
            tmp_path
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
        sim_name="Polymer",
    )

    filename = (
        writer
        .write_second_stage_reaction_files(
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
        "_reacted_1M-3.5_5.0A"
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


def test_final_write_data_uses_nofix(
    tmp_path,
):
    template, *_ = (
        make_real_template(
            tmp_path
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
        .write_second_stage_reaction_files(
            make_simulation()
        )
    )

    commands = command_tokens(
        read_script(
            tmp_path,
            filename,
        )
    )

    write_data = [
        command
        for command in commands
        if command[0] == "write_data"
    ]

    assert len(
        write_data
    ) == 1

    assert (
        write_data[0][-1]
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
        tmp_path / "stage2"
    )

    dest.mkdir()

    writer._copy_required_files(
        dest
    )

    assert (
        dest / pre.name
    ).is_file()

    assert (
        dest / post.name
    ).is_file()

    assert (
        dest / standard_map.name
    ).is_file()


def test_copy_required_files_copies_optional_delete_map(
    tmp_path,
):
    """
    DeleteIDs map is copied for the user, but is not automatically used.
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
        tmp_path / "stage2"
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


def test_no_delete_map_is_valid(
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
        tmp_path / "stage2"
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


def test_missing_delete_map_attribute_is_valid(
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
        tmp_path / "stage2"
    )

    dest.mkdir()

    writer._copy_required_files(
        dest
    )

    assert (
        dest / standard_map.name
    ).is_file()


def test_missing_standard_map_raises(
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

    template = make_template(
        reaction_id=1,
        pre_reaction_file=pre,
        post_reaction_file=post,
        map_file=(
            tmp_path
            / "missing.map"
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
        tmp_path / "stage2"
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


def test_missing_pre_template_raises(
    tmp_path,
):
    post = (
        tmp_path / "post.molecule"
    )

    standard_map = (
        tmp_path / "RXN_1.map"
    )

    post.write_text(
        "POST",
        encoding="utf-8",
    )

    standard_map.write_text(
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
        map_file=standard_map,
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
        tmp_path / "stage2"
    )

    dest.mkdir()

    with pytest.raises(
        FileNotFoundError,
    ):
        writer._copy_required_files(
            dest
        )


def test_missing_post_template_raises(
    tmp_path,
):
    pre = (
        tmp_path / "pre.molecule"
    )

    standard_map = (
        tmp_path / "RXN_1.map"
    )

    pre.write_text(
        "PRE",
        encoding="utf-8",
    )

    standard_map.write_text(
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
        map_file=standard_map,
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
        tmp_path / "stage2"
    )

    dest.mkdir()

    with pytest.raises(
        FileNotFoundError,
    ):
        writer._copy_required_files(
            dest
        )


def test_declared_but_missing_delete_map_raises(
    tmp_path,
):
    template, *_ = (
        make_real_template(
            tmp_path,
            reaction_id=5,
        )
    )

    template.map_file_with_delete_ids = (
        tmp_path
        / "RXN_5_with_delete_ids.map"
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
        tmp_path / "stage2"
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


def test_inactive_templates_are_not_copied(
    tmp_path,
):
    active, *_ = (
        make_real_template(
            tmp_path,
            reaction_id=1,
        )
    )

    inactive = make_template(
        reaction_id=2,
        activity_stats=False,
        map_file=(
            tmp_path
            / "missing.map"
        ),
        pre_reaction_file=(
            tmp_path
            / "missing_pre"
        ),
        post_reaction_file=(
            tmp_path
            / "missing_post"
        ),
        map_file_with_delete_ids=(
            tmp_path
            / "missing_delete"
        ),
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                template_files=[
                    active,
                    inactive,
                ]
            )
        ),
    )

    dest = (
        tmp_path / "stage2"
    )

    dest.mkdir()

    writer._copy_required_files(
        dest
    )

    assert (
        dest / "RXN_1.map"
    ).is_file()

    assert not (
        dest / "missing.map"
    ).exists()


# =============================================================================
# Delete-map full integration
# =============================================================================


def test_delete_map_is_copied_but_standard_map_remains_used(
    tmp_path,
):
    """
    Full contract:

    Files supplied:
        RXN_20.map
        RXN_20_with_delete_ids.map

    Generated LAMMPS script:
        uses RXN_20.map
    """
    (
        template,
        pre,
        post,
        standard_map,
        delete_map,
    ) = make_real_template(
        tmp_path,
        reaction_id=20,
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
        .write_second_stage_reaction_files(
            make_simulation(
                tag="500K"
            )
        )
    )

    stage_dir = (
        tmp_path
        / "4_reaction_second_stage"
    )

    assert (
        stage_dir / pre.name
    ).is_file()

    assert (
        stage_dir / post.name
    ).is_file()

    assert (
        stage_dir
        / standard_map.name
    ).is_file()

    assert (
        stage_dir
        / delete_map.name
    ).is_file()

    active = "\n".join(
        active_lines(
            read_script(
                tmp_path,
                filename,
            )
        )
    )

    assert (
        standard_map.name
        in active
    )

    assert (
        delete_map.name
        not in active
    )
