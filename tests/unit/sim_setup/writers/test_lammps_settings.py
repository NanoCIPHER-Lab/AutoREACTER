from pathlib import Path
from types import SimpleNamespace

import pytest

from AutoREACTER.sim_setup.writers.lammps_settings import (
    LammpsInitialSettings,
    LammpsSettings,
)


# =============================================================================
# Helpers
# =============================================================================


def make_reacter_files(
    in_file,
):
    return SimpleNamespace(
        in_file=Path(in_file),
    )


# =============================================================================
# LammpsSettings dataclass
# =============================================================================


def test_lammps_settings_stores_all_values():
    settings = LammpsSettings(
        units="real",
        dimension="3",
        boundary="p p p",
        atom_style="full",
        bond_style="class2",
        angle_style="class2",
        dihedral_style="class2",
        improper_style="class2",
        special_bonds="lj/coul 0 0 1",
        pair_style="lj/class2/coul/long 8.5",
        kspace_style="pppm 1e-4",
        pair_modify="tail yes mix sixthpower",
        neighbor="2.0 bin",
        neigh_modify="delay 0 every 1 check yes",
    )

    assert settings.units == "real"
    assert settings.dimension == "3"
    assert settings.boundary == "p p p"
    assert settings.atom_style == "full"
    assert settings.bond_style == "class2"
    assert settings.angle_style == "class2"
    assert settings.dihedral_style == "class2"
    assert settings.improper_style == "class2"
    assert settings.special_bonds == "lj/coul 0 0 1"
    assert settings.pair_style == "lj/class2/coul/long 8.5"
    assert settings.kspace_style == "pppm 1e-4"
    assert settings.pair_modify == "tail yes mix sixthpower"
    assert settings.neighbor == "2.0 bin"
    assert settings.neigh_modify == "delay 0 every 1 check yes"


def test_lammps_settings_uses_slots():
    settings = LammpsSettings(
        units="real",
        dimension="3",
        boundary="p p p",
        atom_style="full",
        bond_style="class2",
        angle_style="class2",
        dihedral_style="class2",
        improper_style="class2",
        special_bonds="lj/coul 0 0 1",
        pair_style="lj/class2/coul/long 8.5",
        kspace_style="pppm 1e-4",
        pair_modify="tail yes mix sixthpower",
        neighbor=None,
        neigh_modify=None,
    )

    with pytest.raises(AttributeError):
        settings.extra = "not allowed"


# =============================================================================
# Constructor
# =============================================================================


def test_constructor_stores_reacter_files(
    tmp_path,
):
    in_file = tmp_path / "in.create_atoms.script"

    reacter_files = make_reacter_files(
        in_file
    )

    result = LammpsInitialSettings(
        reacter_files
    )

    assert result.reacter_files is reacter_files


def test_constructor_stores_lunar_input_file_location(
    tmp_path,
):
    in_file = tmp_path / "in.create_atoms.script"

    reacter_files = make_reacter_files(
        in_file
    )

    result = LammpsInitialSettings(
        reacter_files
    )

    assert (
        result.lunar_in_file_location
        == in_file
    )


# =============================================================================
# _defults
# =============================================================================


def test_defaults_returns_lammps_settings(
    tmp_path,
):
    initial = LammpsInitialSettings(
        make_reacter_files(
            tmp_path / "unused.script"
        )
    )

    result = initial._defults()

    assert isinstance(
        result,
        LammpsSettings,
    )


def test_defaults_exact_values(
    tmp_path,
):
    initial = LammpsInitialSettings(
        make_reacter_files(
            tmp_path / "unused.script"
        )
    )

    result = initial._defults()

    assert result.units == "real"
    assert result.dimension == "3"
    assert result.boundary == "p p p"
    assert result.atom_style == "full"

    assert result.bond_style == "class2"
    assert result.angle_style == "class2"
    assert result.dihedral_style == "class2"
    assert result.improper_style == "class2"

    assert result.special_bonds == (
        "lj/coul 0 0 1"
    )

    assert result.pair_style == (
        "lj/class2/coul/long 8.5"
    )

    assert result.kspace_style == (
        "pppm 1e-4"
    )

    assert result.pair_modify == (
        "tail yes mix sixthpower"
    )

    assert result.neighbor is None
    assert result.neigh_modify is None


def test_defaults_returns_new_object_each_time(
    tmp_path,
):
    initial = LammpsInitialSettings(
        make_reacter_files(
            tmp_path / "unused.script"
        )
    )

    first = initial._defults()
    second = initial._defults()

    assert first is not second


# =============================================================================
# get_LUNAR_lammps_settings - defaults
# =============================================================================


def test_missing_lunar_input_file_uses_defaults(
    tmp_path,
):
    initial = LammpsInitialSettings(
        make_reacter_files(
            tmp_path / "missing.script"
        )
    )

    result = (
        initial
        .get_LUNAR_lammps_settings()
    )

    assert result.units == "real"
    assert result.dimension == "3"
    assert result.boundary == "p p p"
    assert result.atom_style == "full"

    assert result.bond_style == "class2"
    assert result.angle_style == "class2"
    assert result.dihedral_style == "class2"
    assert result.improper_style == "class2"

    assert result.special_bonds == (
        "lj/coul 0 0 1"
    )

    assert result.pair_style == (
        "lj/class2/coul/long 8.5"
    )

    assert result.kspace_style == (
        "pppm 1e-4"
    )

    assert result.pair_modify == (
        "tail yes mix sixthpower"
    )

    assert result.neighbor is None
    assert result.neigh_modify is None


def test_missing_lunar_input_file_prints_warning(
    tmp_path,
    capsys,
):
    missing = (
        tmp_path / "missing.script"
    )

    initial = LammpsInitialSettings(
        make_reacter_files(
            missing
        )
    )

    initial.get_LUNAR_lammps_settings()

    output = (
        capsys.readouterr().out
    )

    assert (
        "Could not read LAMMPS input file"
        in output
    )

    assert str(missing) in output

    assert (
        "Proceeding with default LAMMPS settings"
        in output
    )


def test_empty_lunar_input_file_uses_defaults(
    tmp_path,
):
    in_file = (
        tmp_path / "in.script"
    )

    in_file.write_text(
        "",
        encoding="utf-8",
    )

    initial = LammpsInitialSettings(
        make_reacter_files(
            in_file
        )
    )

    result = (
        initial
        .get_LUNAR_lammps_settings()
    )

    assert result.units == "real"

    assert result.pair_style == (
        "lj/class2/coul/long 8.5"
    )

    assert result.neighbor is None


# =============================================================================
# Override parsing
# =============================================================================


def test_lunar_file_overrides_single_setting(
    tmp_path,
):
    in_file = tmp_path / "in.script"

    in_file.write_text(
        "units metal\n",
        encoding="utf-8",
    )

    initial = LammpsInitialSettings(
        make_reacter_files(
            in_file
        )
    )

    result = (
        initial
        .get_LUNAR_lammps_settings()
    )

    assert result.units == "metal"

    # Everything else remains default.
    assert result.dimension == "3"

    assert result.atom_style == "full"


def test_lunar_file_overrides_multiple_settings(
    tmp_path,
):
    in_file = tmp_path / "in.script"

    in_file.write_text(
        """units real
dimension 2
boundary f f p
atom_style molecular
bond_style harmonic
angle_style harmonic
dihedral_style opls
improper_style cvff
special_bonds lj/coul 0.0 0.0 0.5
pair_style lj/cut 10.0
kspace_style ewald 1e-5
pair_modify mix arithmetic
neighbor 3.0 bin
neigh_modify delay 5 every 2 check no
""",
        encoding="utf-8",
    )

    initial = LammpsInitialSettings(
        make_reacter_files(
            in_file
        )
    )

    result = (
        initial
        .get_LUNAR_lammps_settings()
    )

    assert result.units == "real"
    assert result.dimension == "2"
    assert result.boundary == "f f p"

    assert (
        result.atom_style
        == "molecular"
    )

    assert (
        result.bond_style
        == "harmonic"
    )

    assert (
        result.angle_style
        == "harmonic"
    )

    assert (
        result.dihedral_style
        == "opls"
    )

    assert (
        result.improper_style
        == "cvff"
    )

    assert (
        result.special_bonds
        == "lj/coul 0.0 0.0 0.5"
    )

    assert (
        result.pair_style
        == "lj/cut 10.0"
    )

    assert (
        result.kspace_style
        == "ewald 1e-5"
    )

    assert (
        result.pair_modify
        == "mix arithmetic"
    )

    assert (
        result.neighbor
        == "3.0 bin"
    )

    assert (
        result.neigh_modify
        == "delay 5 every 2 check no"
    )


def test_inline_comments_are_removed(
    tmp_path,
):
    in_file = tmp_path / "in.script"

    in_file.write_text(
        """
pair_style lj/class2/coul/long 12.0 # cutoff
kspace_style pppm 1.0e-4 # electrostatics
neighbor 2.0 bin # neighbor list
""",
        encoding="utf-8",
    )

    initial = LammpsInitialSettings(
        make_reacter_files(
            in_file
        )
    )

    result = (
        initial
        .get_LUNAR_lammps_settings()
    )

    assert (
        result.pair_style
        == "lj/class2/coul/long 12.0"
    )

    assert (
        result.kspace_style
        == "pppm 1.0e-4"
    )

    assert (
        result.neighbor
        == "2.0 bin"
    )


def test_full_comment_lines_are_ignored(
    tmp_path,
):
    in_file = tmp_path / "in.script"

    in_file.write_text(
        """# units metal
# pair_style lj/cut 15.0
units real
""",
        encoding="utf-8",
    )

    initial = LammpsInitialSettings(
        make_reacter_files(
            in_file
        )
    )

    result = (
        initial
        .get_LUNAR_lammps_settings()
    )

    assert result.units == "real"

    assert (
        result.pair_style
        == "lj/class2/coul/long 8.5"
    )


def test_blank_lines_are_ignored(
    tmp_path,
):
    in_file = tmp_path / "in.script"

    in_file.write_text(
        """

units metal


pair_style lj/cut 9.0


""",
        encoding="utf-8",
    )

    initial = LammpsInitialSettings(
        make_reacter_files(
            in_file
        )
    )

    result = (
        initial
        .get_LUNAR_lammps_settings()
    )

    assert result.units == "metal"

    assert (
        result.pair_style
        == "lj/cut 9.0"
    )


# =============================================================================
# Prefix matching behavior
# =============================================================================


def test_keyword_must_begin_clean_line(
    tmp_path,
):
    in_file = tmp_path / "in.script"

    in_file.write_text(
        """variable units string metal
units real
""",
        encoding="utf-8",
    )

    initial = LammpsInitialSettings(
        make_reacter_files(
            in_file
        )
    )

    result = (
        initial
        .get_LUNAR_lammps_settings()
    )

    assert result.units == "real"


def test_leading_whitespace_before_keyword_is_allowed(
    tmp_path,
):
    in_file = tmp_path / "in.script"

    in_file.write_text(
        "    units metal\n",
        encoding="utf-8",
    )

    initial = LammpsInitialSettings(
        make_reacter_files(
            in_file
        )
    )

    result = (
        initial
        .get_LUNAR_lammps_settings()
    )

    assert result.units == "metal"


def test_empty_value_does_not_override_default(
    tmp_path,
):
    in_file = tmp_path / "in.script"

    in_file.write_text(
        "units\n",
        encoding="utf-8",
    )

    initial = LammpsInitialSettings(
        make_reacter_files(
            in_file
        )
    )

    result = (
        initial
        .get_LUNAR_lammps_settings()
    )

    assert result.units == "real"


# =============================================================================
# Current duplicate-line behavior
# =============================================================================


def test_last_matching_line_wins(
    tmp_path,
):
    """
    Characterization of current nested-loop behavior.

    Every matching line updates the same setting, so the last one wins.
    """
    in_file = tmp_path / "in.script"

    in_file.write_text(
        """pair_style lj/cut 8.0
pair_style lj/cut 10.0
pair_style lj/class2/coul/long 12.0
""",
        encoding="utf-8",
    )

    initial = LammpsInitialSettings(
        make_reacter_files(
            in_file
        )
    )

    result = (
        initial
        .get_LUNAR_lammps_settings()
    )

    assert (
        result.pair_style
        == "lj/class2/coul/long 12.0"
    )


def test_last_neighbor_line_wins(
    tmp_path,
):
    in_file = tmp_path / "in.script"

    in_file.write_text(
        """neighbor 1.0 bin
neighbor 2.0 bin
neighbor 3.0 bin
""",
        encoding="utf-8",
    )

    initial = LammpsInitialSettings(
        make_reacter_files(
            in_file
        )
    )

    result = (
        initial
        .get_LUNAR_lammps_settings()
    )

    assert (
        result.neighbor
        == "3.0 bin"
    )


# =============================================================================
# Current startswith semantics
# =============================================================================


def test_keyword_matching_currently_uses_startswith(
    tmp_path,
):
    """
    Characterization test.

    Current code uses:

        clean_line.startswith(keyword)

    rather than requiring a whitespace boundary after the keyword.

    Therefore a line such as "units_extra metal" is interpreted as a
    units setting with value "_extra metal".

    This is odd, but do not change production during characterization.
    """
    in_file = tmp_path / "in.script"

    in_file.write_text(
        "units_extra metal\n",
        encoding="utf-8",
    )

    initial = LammpsInitialSettings(
        make_reacter_files(
            in_file
        )
    )

    result = (
        initial
        .get_LUNAR_lammps_settings()
    )

    assert (
        result.units
        == "_extra metal"
    )


def test_pair_style_does_not_override_pair_modify(
    tmp_path,
):
    in_file = tmp_path / "in.script"

    in_file.write_text(
        """pair_style lj/cut 10.0
pair_modify mix arithmetic
""",
        encoding="utf-8",
    )

    initial = LammpsInitialSettings(
        make_reacter_files(
            in_file
        )
    )

    result = (
        initial
        .get_LUNAR_lammps_settings()
    )

    assert (
        result.pair_style
        == "lj/cut 10.0"
    )

    assert (
        result.pair_modify
        == "mix arithmetic"
    )


# =============================================================================
# Realistic LUNAR-style input
# =============================================================================


def test_realistic_lunar_settings_file(
    tmp_path,
):
    in_file = (
        tmp_path
        / "in.create_atoms.script"
    )

    in_file.write_text(
        """# LUNAR generated settings

units            real
dimension        3
boundary         p p p
atom_style       full

bond_style       class2
angle_style      class2
dihedral_style   class2
improper_style   class2

pair_style       lj/class2/coul/long 12.0
kspace_style     pppm 1.0e-4
pair_modify      mix sixthpower

neighbor         2.0 bin
neigh_modify     delay 0 every 1 check yes one 5000 page 100000
""",
        encoding="utf-8",
    )

    initial = LammpsInitialSettings(
        make_reacter_files(
            in_file
        )
    )

    result = (
        initial
        .get_LUNAR_lammps_settings()
    )

    assert result.units == "real"
    assert result.dimension == "3"
    assert result.boundary == "p p p"
    assert result.atom_style == "full"

    assert result.bond_style == "class2"
    assert result.angle_style == "class2"
    assert result.dihedral_style == "class2"
    assert result.improper_style == "class2"

    assert (
        result.pair_style
        == "lj/class2/coul/long 12.0"
    )

    assert (
        result.kspace_style
        == "pppm 1.0e-4"
    )

    assert (
        result.pair_modify
        == "mix sixthpower"
    )

    assert (
        result.neighbor
        == "2.0 bin"
    )

    assert (
        result.neigh_modify
        == (
            "delay 0 every 1 check yes "
            "one 5000 page 100000"
        )
    )