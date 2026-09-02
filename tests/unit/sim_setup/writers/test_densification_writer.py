from pathlib import Path
from types import SimpleNamespace

import pytest

import AutoREACTER.sim_setup.writers.densification_writer as dens_module
from AutoREACTER.sim_setup.writers.densification_writer import (
    DensificationWriter,
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


def make_monomer(
    *,
    name="mma",
    data_id="data1",
    molecule_file=None,
):
    return SimpleNamespace(
        name=name,
        data_id=data_id,
        lmp_molecule_file=molecule_file,
    )


def make_reacter_files(
    *,
    force_field_data,
    molecule_files=None,
    template_files=None,
):
    return SimpleNamespace(
        force_field_data=Path(force_field_data),
        molecule_files=list(
            molecule_files or []
        ),
        template_files=list(
            template_files or []
        ),
    )


def make_simulation(
    *,
    tag="sim1",
    initial_box_length=40.0,
    temperature=500.0,
    density=1.2,
    monomer_counts=None,
):
    return SimpleNamespace(
        tag=tag,
        initial_box_length=initial_box_length,
        temperature=temperature,
        density=density,
        monomer_counts=(
            {"mma": 10}
            if monomer_counts is None
            else monomer_counts
        ),
    )


def make_writer_without_init(
    tmp_path,
    *,
    settings=None,
    reacter_files=None,
    sim_name="Test",
):
    writer = object.__new__(
        DensificationWriter
    )

    writer.settings = (
        settings
        if settings is not None
        else make_settings()
    )

    if reacter_files is None:
        ff_file = (
            tmp_path / "force_field.data"
        )

        ff_file.write_text(
            "ff",
            encoding="utf-8",
        )

        reacter_files = make_reacter_files(
            force_field_data=ff_file,
        )

    writer.reacter_files = (
        reacter_files
    )

    writer.erate = -0.001
    writer.timestep = 0.01
    writer.out_dir = tmp_path
    writer.sim_name = sim_name

    return writer


# =============================================================================
# Constructor
# =============================================================================


def test_constructor_stores_configuration(
    tmp_path,
    monkeypatch,
):
    settings = make_settings()

    ff_file = (
        tmp_path / "force_field.data"
    )

    ff_file.write_text(
        "ff",
        encoding="utf-8",
    )

    reacter_files = make_reacter_files(
        force_field_data=ff_file,
    )

    simulation = make_simulation()

    calls = []

    monkeypatch.setattr(
        DensificationWriter,
        "write_lammps_densification_file",
        lambda self, simulation:
        (
            calls.append(simulation)
            or "in.test_densification"
        ),
    )

    writer = DensificationWriter(
        out_dir=tmp_path,
        settings=settings,
        reacter_files=reacter_files,
        simulation=simulation,
        sim_name="Test",
    )

    assert writer.settings is settings

    assert (
        writer.reacter_files
        is reacter_files
    )

    assert writer.out_dir == tmp_path

    assert writer.sim_name == "Test"

    assert writer.erate == pytest.approx(
        -0.001
    )

    assert writer.timestep == pytest.approx(
        0.01
    )

    assert calls == [
        simulation
    ]

    assert (
        writer.in_dense_file_name
        == "in.test_densification"
    )


# =============================================================================
# _write_empty_box_data
# =============================================================================


def test_write_empty_box_creates_file(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    output_dir = (
        tmp_path / "box"
    )

    output_dir.mkdir()

    simulation = make_simulation(
        initial_box_length=40.0,
    )

    writer._write_empty_box_data(
        simulation,
        output_dir,
    )

    path = (
        output_dir
        / "empty_box.data"
    )

    assert path.is_file()


def test_write_empty_box_is_centered_at_origin(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    output_dir = (
        tmp_path / "box"
    )

    output_dir.mkdir()

    writer._write_empty_box_data(
        make_simulation(
            initial_box_length=40.0,
        ),
        output_dir,
    )

    text = (
        output_dir
        / "empty_box.data"
    ).read_text(
        encoding="utf-8"
    )

    assert (
        "-20.00 20.00 xlo xhi"
        in text
    )

    assert (
        "-20.00 20.00 ylo yhi"
        in text
    )

    assert (
        "-20.00 20.00 zlo zhi"
        in text
    )


def test_write_empty_box_contains_zero_topology_counts(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    output_dir = (
        tmp_path / "box"
    )

    output_dir.mkdir()

    writer._write_empty_box_data(
        make_simulation(),
        output_dir,
    )

    text = (
        output_dir
        / "empty_box.data"
    ).read_text(
        encoding="utf-8"
    )

    assert "0 atoms" in text
    assert "0 bonds" in text
    assert "0 angles" in text
    assert "0 dihedrals" in text
    assert "0 impropers" in text


def test_write_empty_box_rounds_bounds_to_two_decimals(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    output_dir = (
        tmp_path / "box"
    )

    output_dir.mkdir()

    writer._write_empty_box_data(
        make_simulation(
            initial_box_length=33.333,
        ),
        output_dir,
    )

    text = (
        output_dir
        / "empty_box.data"
    ).read_text(
        encoding="utf-8"
    )

    assert (
        "-16.67 16.67 xlo xhi"
        in text
    )


# =============================================================================
# _get_force_field_types
# =============================================================================


def test_get_force_field_types_reads_all_counts(
    tmp_path,
):
    ff_file = (
        tmp_path / "ff.data"
    )

    ff_file.write_text(
        """LAMMPS data file

12 atom types
8 bond types
6 angle types
4 dihedral types
3 improper types

Masses
""",
        encoding="utf-8",
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                force_field_data=ff_file,
            )
        ),
    )

    result = (
        writer._get_force_field_types()
    )

    assert result == {
        "atom": 12,
        "bond": 8,
        "angle": 6,
        "dihedral": 4,
        "improper": 3,
    }


def test_get_force_field_types_defaults_improper_to_zero(
    tmp_path,
):
    ff_file = (
        tmp_path / "ff.data"
    )

    ff_file.write_text(
        """10 atom types
5 bond types
4 angle types
3 dihedral types

Masses
""",
        encoding="utf-8",
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                force_field_data=ff_file,
            )
        ),
    )

    result = (
        writer._get_force_field_types()
    )

    assert (
        result["improper"]
        == 0
    )


def test_get_force_field_types_missing_file_raises(
    tmp_path,
):
    missing = (
        tmp_path / "missing.data"
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                force_field_data=missing,
            )
        ),
    )

    with pytest.raises(
        FileNotFoundError,
        match=(
            "Force field data file "
            "not found"
        ),
    ):
        writer._get_force_field_types()


@pytest.mark.parametrize(
    "missing_key",
    [
        "atom",
        "bond",
        "angle",
        "dihedral",
    ],
)
def test_get_force_field_types_requires_core_headers(
    tmp_path,
    missing_key,
):
    values = {
        "atom": 10,
        "bond": 8,
        "angle": 6,
        "dihedral": 4,
    }

    lines = []

    for key, value in (
        values.items()
    ):
        if key != missing_key:
            lines.append(
                f"{value} {key} types"
            )

    lines.append("")
    lines.append("Masses")

    ff_file = (
        tmp_path / "ff.data"
    )

    ff_file.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                force_field_data=ff_file,
            )
        ),
    )

    with pytest.raises(
        ValueError,
        match=(
            f"Required LAMMPS header "
            f"'{missing_key} types'"
        ),
    ):
        writer._get_force_field_types()


def test_get_force_field_types_stops_at_mass_section(
    tmp_path,
):
    ff_file = (
        tmp_path / "ff.data"
    )

    ff_file.write_text(
        """10 atom types

Masses

8 bond types
6 angle types
4 dihedral types
""",
        encoding="utf-8",
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                force_field_data=ff_file,
            )
        ),
    )

    with pytest.raises(
        ValueError,
        match="bond types",
    ):
        writer._get_force_field_types()


# =============================================================================
# _run_calculation
# =============================================================================


def test_default_run_calculation(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    total_steps, thermo = (
        writer._run_calculation()
    )

    assert total_steps == 50000
    assert thermo == 1000


@pytest.mark.parametrize(
    (
        "timestep, "
        "expected_steps, "
        "expected_thermo"
    ),
    [
        (
            0.01,
            50000,
            1000,
        ),
        (
            0.1,
            5000,
            100,
        ),
        (
            1.0,
            500,
            10,
        ),
        (
            10.0,
            50,
            1,
        ),
        (
            100.0,
            5,
            1,
        ),
    ],
)
def test_run_calculation_rounding_branches(
    tmp_path,
    timestep,
    expected_steps,
    expected_thermo,
):
    writer = make_writer_without_init(
        tmp_path
    )

    writer.erate = -0.001
    writer.timestep = timestep

    total_steps, thermo = (
        writer._run_calculation()
    )

    assert (
        total_steps
        == expected_steps
    )

    assert (
        thermo
        == expected_thermo
    )


def test_run_calculation_positive_erate_raises(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path
    )

    writer.erate = 0.001

    with pytest.raises(
        ValueError,
        match=(
            "Calculated negative steps"
        ),
    ):
        writer._run_calculation()


# =============================================================================
# write_lammps_densification_file
# =============================================================================


def test_write_densification_creates_stage_directory_and_script(
    tmp_path,
    monkeypatch,
):
    ff_file = (
        tmp_path / "force_field.data"
    )

    ff_file.write_text(
        "ff",
        encoding="utf-8",
    )

    monomer_file = (
        tmp_path / "mma.molecule"
    )

    monomer_file.write_text(
        "mol",
        encoding="utf-8",
    )

    rf = make_reacter_files(
        force_field_data=ff_file,
        molecule_files=[
            make_monomer(
                name="mma",
                molecule_file=monomer_file,
            )
        ],
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=rf,
        sim_name="Polymer",
    )

    monkeypatch.setattr(
        writer,
        "_get_force_field_types",
        lambda: {
            "atom": 10,
            "bond": 20,
            "angle": 30,
            "dihedral": 40,
            "improper": 5,
        },
    )

    monkeypatch.setattr(
        writer,
        "_run_calculation",
        lambda: (
            12345,
            250,
        ),
    )

    monkeypatch.setattr(
        writer,
        "_copy_required_files",
        lambda dest_dir: None,
    )

    monkeypatch.setattr(
        dens_module.random,
        "randint",
        lambda a, b: 12345,
    )

    simulation = make_simulation(
        tag="500K",
        temperature=500.0,
        density=1.23,
        monomer_counts={
            "mma": 25,
        },
    )

    filename = (
        writer
        .write_lammps_densification_file(
            simulation
        )
    )

    stage_dir = (
        tmp_path
        / "1_densification"
    )

    assert stage_dir.is_dir()

    assert filename == (
        "in.Polymer_500K_densification"
    )

    assert (
        stage_dir / filename
    ).is_file()

    assert (
        stage_dir
        / "empty_box.data"
    ).is_file()


def test_write_densification_contains_force_field_settings(
    tmp_path,
    monkeypatch,
):
    ff_file = (
        tmp_path / "force_field.data"
    )

    ff_file.write_text(
        "ff",
        encoding="utf-8",
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
                force_field_data=ff_file,
            )
        ),
    )

    monkeypatch.setattr(
        writer,
        "_get_force_field_types",
        lambda: {
            "atom": 1,
            "bond": 2,
            "angle": 3,
            "dihedral": 4,
            "improper": 5,
        },
    )

    monkeypatch.setattr(
        writer,
        "_run_calculation",
        lambda: (
            1000,
            100,
        ),
    )

    monkeypatch.setattr(
        writer,
        "_copy_required_files",
        lambda dest_dir: None,
    )

    monkeypatch.setattr(
        dens_module.random,
        "randint",
        lambda a, b: 11111,
    )

    filename = (
        writer
        .write_lammps_densification_file(
            make_simulation()
        )
    )

    text = (
        tmp_path
        / "1_densification"
        / filename
    ).read_text(
        encoding="utf-8"
    )

    assert (
        "units            real"
        in text
    )

    assert (
        "dimension        3"
        in text
    )

    assert (
        "boundary         p p p"
        in text
    )

    assert (
        "atom_style       full"
        in text
    )

    assert (
        "bond_style       class2"
        in text
    )

    assert (
        "angle_style      class2"
        in text
    )

    assert (
        "dihedral_style   class2"
        in text
    )

    assert (
        "improper_style   class2"
        in text
    )

    assert (
        "pair_style       "
        "lj/class2/coul/long 12.0"
        in text
    )

    assert (
        "neighbor         "
        "2.0 bin"
        in text
    )

    assert (
        "neigh_modify     "
        "delay 0 every 1 check yes"
        in text
    )


def test_write_densification_omits_optional_neighbor_commands_when_none(
    tmp_path,
    monkeypatch,
):
    ff_file = (
        tmp_path / "force_field.data"
    )

    ff_file.write_text(
        "ff",
        encoding="utf-8",
    )

    writer = make_writer_without_init(
        tmp_path,
        settings=make_settings(
            neighbor=None,
            neigh_modify=None,
        ),
        reacter_files=(
            make_reacter_files(
                force_field_data=ff_file,
            )
        ),
    )

    monkeypatch.setattr(
        writer,
        "_get_force_field_types",
        lambda: {
            "atom": 1,
            "bond": 1,
            "angle": 1,
            "dihedral": 1,
            "improper": 0,
        },
    )

    monkeypatch.setattr(
        writer,
        "_run_calculation",
        lambda: (
            1000,
            100,
        ),
    )

    monkeypatch.setattr(
        writer,
        "_copy_required_files",
        lambda dest_dir: None,
    )

    monkeypatch.setattr(
        dens_module.random,
        "randint",
        lambda a, b: 11111,
    )

    filename = (
        writer
        .write_lammps_densification_file(
            make_simulation()
        )
    )

    text = (
        tmp_path
        / "1_densification"
        / filename
    ).read_text(
        encoding="utf-8"
    )

    assert (
        "\nneighbor "
        not in text
    )

    assert (
        "\nneigh_modify "
        not in text
    )


def test_write_densification_uses_force_field_type_counts(
    tmp_path,
    monkeypatch,
):
    ff_file = (
        tmp_path / "force_field.data"
    )

    ff_file.write_text(
        "ff",
        encoding="utf-8",
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                force_field_data=ff_file,
            )
        ),
    )

    monkeypatch.setattr(
        writer,
        "_get_force_field_types",
        lambda: {
            "atom": 11,
            "bond": 22,
            "angle": 33,
            "dihedral": 44,
            "improper": 55,
        },
    )

    monkeypatch.setattr(
        writer,
        "_run_calculation",
        lambda: (
            1000,
            100,
        ),
    )

    monkeypatch.setattr(
        writer,
        "_copy_required_files",
        lambda dest_dir: None,
    )

    monkeypatch.setattr(
        dens_module.random,
        "randint",
        lambda a, b: 11111,
    )

    filename = (
        writer
        .write_lammps_densification_file(
            make_simulation()
        )
    )

    text = (
        tmp_path
        / "1_densification"
        / filename
    ).read_text(
        encoding="utf-8"
    )

    assert (
        "extra/atom/types 11"
        in text
    )

    assert (
        "extra/bond/types 22"
        in text
    )

    assert (
        "extra/angle/types 33"
        in text
    )

    assert (
        "extra/dihedral/types 44"
        in text
    )

    assert (
        "extra/improper/types 55"
        in text
    )


def test_write_densification_defines_molecules_and_inserts_counts(
    tmp_path,
    monkeypatch,
):
    ff_file = (
        tmp_path / "force_field.data"
    )

    ff_file.write_text(
        "ff",
        encoding="utf-8",
    )

    mma_file = (
        tmp_path / "mma.molecule"
    )

    tegdma_file = (
        tmp_path / "tegdma.molecule"
    )

    mma_file.write_text(
        "mma",
        encoding="utf-8",
    )

    tegdma_file.write_text(
        "tegdma",
        encoding="utf-8",
    )

    rf = make_reacter_files(
        force_field_data=ff_file,
        molecule_files=[
            make_monomer(
                name="mma",
                molecule_file=mma_file,
            ),
            make_monomer(
                name="tegdma",
                molecule_file=tegdma_file,
            ),
        ],
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=rf,
    )

    monkeypatch.setattr(
        writer,
        "_get_force_field_types",
        lambda: {
            "atom": 1,
            "bond": 1,
            "angle": 1,
            "dihedral": 1,
            "improper": 0,
        },
    )

    monkeypatch.setattr(
        writer,
        "_run_calculation",
        lambda: (
            5000,
            100,
        ),
    )

    monkeypatch.setattr(
        writer,
        "_copy_required_files",
        lambda dest_dir: None,
    )

    seeds = iter(
        [
            11111,
            22222,
            33333,
        ]
    )

    monkeypatch.setattr(
        dens_module.random,
        "randint",
        lambda a, b:
        next(seeds),
    )

    simulation = make_simulation(
        monomer_counts={
            "mma": 25,
            "tegdma": 5,
            "unknown": 999,
        },
    )

    filename = (
        writer
        .write_lammps_densification_file(
            simulation
        )
    )

    text = (
        tmp_path
        / "1_densification"
        / filename
    ).read_text(
        encoding="utf-8"
    )

    assert (
        "mol_1 mma.molecule"
        in text
    )

    assert (
        "mol_2 tegdma.molecule"
        in text
    )

    assert (
        "random 25 11111 NULL "
        "overlap 2.0 maxtry 100 "
        "mol mol_1 11111"
        in text
    )

    assert (
        "random 5 22222 NULL "
        "overlap 2.0 maxtry 100 "
        "mol mol_2 22222"
        in text
    )

    assert (
        "random 999"
        not in text
    )

    assert (
        "all create 500.0 33333 "
        "dist gaussian"
        in text
    )


def test_molecule_name_falls_back_to_data_id(
    tmp_path,
    monkeypatch,
):
    ff_file = (
        tmp_path / "force_field.data"
    )

    ff_file.write_text(
        "ff",
        encoding="utf-8",
    )

    mol_file = (
        tmp_path / "data42.molecule"
    )

    mol_file.write_text(
        "mol",
        encoding="utf-8",
    )

    rf = make_reacter_files(
        force_field_data=ff_file,
        molecule_files=[
            make_monomer(
                name=None,
                data_id="data42",
                molecule_file=mol_file,
            )
        ],
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=rf,
    )

    monkeypatch.setattr(
        writer,
        "_get_force_field_types",
        lambda: {
            "atom": 1,
            "bond": 1,
            "angle": 1,
            "dihedral": 1,
            "improper": 0,
        },
    )

    monkeypatch.setattr(
        writer,
        "_run_calculation",
        lambda: (
            100,
            10,
        ),
    )

    monkeypatch.setattr(
        writer,
        "_copy_required_files",
        lambda dest_dir: None,
    )

    monkeypatch.setattr(
        dens_module.random,
        "randint",
        lambda a, b: 12345,
    )

    filename = (
        writer
        .write_lammps_densification_file(
            make_simulation(
                monomer_counts={
                    "data42": 7
                }
            )
        )
    )

    text = (
        tmp_path
        / "1_densification"
        / filename
    ).read_text(
        encoding="utf-8"
    )

    assert (
        "mol_1 data42.molecule"
        in text
    )

    assert (
        "random 7 12345"
        in text
    )


def test_write_densification_contains_run_conditions(
    tmp_path,
    monkeypatch,
):
    ff_file = (
        tmp_path / "force_field.data"
    )

    ff_file.write_text(
        "ff",
        encoding="utf-8",
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                force_field_data=ff_file,
            )
        ),
    )

    monkeypatch.setattr(
        writer,
        "_get_force_field_types",
        lambda: {
            "atom": 1,
            "bond": 1,
            "angle": 1,
            "dihedral": 1,
            "improper": 0,
        },
    )

    monkeypatch.setattr(
        writer,
        "_run_calculation",
        lambda: (
            45678,
            321,
        ),
    )

    monkeypatch.setattr(
        writer,
        "_copy_required_files",
        lambda dest_dir: None,
    )

    monkeypatch.setattr(
        dens_module.random,
        "randint",
        lambda a, b: 55555,
    )

    simulation = make_simulation(
        temperature=373.0,
        density=1.25,
    )

    filename = (
        writer
        .write_lammps_densification_file(
            simulation
        )
    )

    text = (
        tmp_path
        / "1_densification"
        / filename
    ).read_text(
        encoding="utf-8"
    )

    assert (
        "nvt temp 373.0 373.0 100.0"
        in text
    )

    assert (
        "x erate -0.001 "
        "y erate -0.001 "
        "z erate -0.001"
        in text
    )

    assert (
        "v_my_den > 1.25 "
        "error continue"
        in text
    )

    # Do not make these assertions depend on cosmetic alignment spacing.
    script_lines = (
        text.splitlines()
    )

    assert any(
        line.split()
        == [
            "thermo",
            "321",
        ]
        for line in script_lines
    )

    assert any(
        line.split()
        == [
            "run",
            "45678",
        ]
        for line in script_lines
    )

    assert any(
        line.split()
        == [
            "timestep",
            "0.01",
        ]
        for line in script_lines
    )


def test_write_densification_calls_copy_required_files_with_stage_dir(
    tmp_path,
    monkeypatch,
):
    ff_file = (
        tmp_path / "force_field.data"
    )

    ff_file.write_text(
        "ff",
        encoding="utf-8",
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                force_field_data=ff_file,
            )
        ),
    )

    monkeypatch.setattr(
        writer,
        "_get_force_field_types",
        lambda: {
            "atom": 1,
            "bond": 1,
            "angle": 1,
            "dihedral": 1,
            "improper": 0,
        },
    )

    monkeypatch.setattr(
        writer,
        "_run_calculation",
        lambda: (
            100,
            10,
        ),
    )

    monkeypatch.setattr(
        dens_module.random,
        "randint",
        lambda a, b: 12345,
    )

    calls = []

    monkeypatch.setattr(
        writer,
        "_copy_required_files",
        lambda dest_dir:
        calls.append(
            dest_dir
        ),
    )

    writer.write_lammps_densification_file(
        make_simulation()
    )

    assert calls == [
        tmp_path
        / "1_densification"
    ]


# =============================================================================
# _copy_required_files
# =============================================================================


def test_copy_required_files_copies_force_field(
    tmp_path,
):
    ff_file = (
        tmp_path / "force_field.data"
    )

    ff_file.write_text(
        "FF",
        encoding="utf-8",
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                force_field_data=ff_file,
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

    copied = (
        dest
        / "force_field.data"
    )

    assert copied.is_file()

    assert (
        copied.read_text(
            encoding="utf-8"
        )
        == "FF"
    )


def test_copy_required_files_copies_monomer_molecules(
    tmp_path,
):
    ff_file = (
        tmp_path / "force_field.data"
    )

    ff_file.write_text(
        "FF",
        encoding="utf-8",
    )

    mma = (
        tmp_path / "mma.molecule"
    )

    tegdma = (
        tmp_path / "tegdma.molecule"
    )

    mma.write_text(
        "MMA",
        encoding="utf-8",
    )

    tegdma.write_text(
        "TEGDMA",
        encoding="utf-8",
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                force_field_data=ff_file,
                molecule_files=[
                    make_monomer(
                        name="mma",
                        molecule_file=mma,
                    ),
                    make_monomer(
                        name="tegdma",
                        molecule_file=tegdma,
                    ),
                ],
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

    assert (
        dest
        / "mma.molecule"
    ).read_text(
        encoding="utf-8"
    ) == "MMA"

    assert (
        dest
        / "tegdma.molecule"
    ).read_text(
        encoding="utf-8"
    ) == "TEGDMA"


def test_copy_required_files_ignores_missing_force_field(
    tmp_path,
):
    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                force_field_data=(
                    tmp_path
                    / "missing.data"
                ),
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


def test_copy_required_files_ignores_missing_molecule_file(
    tmp_path,
):
    ff_file = (
        tmp_path / "ff.data"
    )

    ff_file.write_text(
        "ff",
        encoding="utf-8",
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                force_field_data=ff_file,
                molecule_files=[
                    make_monomer(
                        molecule_file=(
                            tmp_path
                            / "missing.molecule"
                        )
                    )
                ],
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

    assert not (
        dest
        / "missing.molecule"
    ).exists()


def test_copy_required_files_ignores_none_molecule_path(
    tmp_path,
):
    ff_file = (
        tmp_path / "ff.data"
    )

    ff_file.write_text(
        "ff",
        encoding="utf-8",
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                force_field_data=ff_file,
                molecule_files=[
                    make_monomer(
                        molecule_file=None
                    )
                ],
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

    assert (
        dest / "ff.data"
    ).is_file()


def test_densification_does_not_copy_reaction_maps(
    tmp_path,
):
    """
    Reaction maps do NOT belong in 1_densification.

    Both:
        RXN_N.map
    and the optional:
        RXN_N_with_delete_ids.map

    belong to the REACTER/reaction workflow.
    """
    ff_file = (
        tmp_path
        / "force_field.data"
    )

    ff_file.write_text(
        "ff",
        encoding="utf-8",
    )

    standard_map = (
        tmp_path
        / "RXN_1.map"
    )

    delete_map = (
        tmp_path
        / "RXN_1_with_delete_ids.map"
    )

    standard_map.write_text(
        "STANDARD",
        encoding="utf-8",
    )

    delete_map.write_text(
        "DELETE",
        encoding="utf-8",
    )

    template = SimpleNamespace(
        map_file=standard_map,
        map_file_with_delete_ids=delete_map,
    )

    writer = make_writer_without_init(
        tmp_path,
        reacter_files=(
            make_reacter_files(
                force_field_data=ff_file,
                template_files=[
                    template
                ],
            )
        ),
    )

    dest = (
        tmp_path
        / "densification"
    )

    dest.mkdir()

    writer._copy_required_files(
        dest
    )

    assert not (
        dest
        / "RXN_1.map"
    ).exists()

    assert not (
        dest
        / "RXN_1_with_delete_ids.map"
    ).exists()