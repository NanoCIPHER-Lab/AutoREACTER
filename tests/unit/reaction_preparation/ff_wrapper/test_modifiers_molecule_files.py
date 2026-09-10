import pandas as pd
import pytest

from AutoREACTER.reaction_preparation.ff_wrapper.modifiers_molecule_files import (
    modify_angles,
    modify_bonds,
    modify_charges,
    modify_coords,
    modify_dihedrals,
    modify_impropers,
    modify_types,
)


# =============================================================================
# Helpers
# =============================================================================


def make_type_df(
    mapping,
):
    """
    Build the minimal type DataFrame required by topology modifiers.

    mapping:
        {
            old_atom_index: new_atom_index,
            ...
        }
    """
    return pd.DataFrame(
        {
            "atom_index": list(
                mapping.keys()
            ),
            "new_atom_index": list(
                mapping.values()
            ),
        }
    )


# =============================================================================
# modify_types
# =============================================================================


def test_modify_types_filters_atoms_and_reindexes():
    lines = [
        "1 10 # C\n",
        "2 20 # H\n",
        "3 30 # O\n",
        "4 40 # N\n",
        "\n",
    ]

    template_indexes = [
        1,
        3,
    ]

    (
        df,
        section,
        count,
        mapping,
        legacy_mode,
    ) = modify_types(
        lines,
        template_indexes,
        0,
    )

    assert list(
        df["atom_index"]
    ) == [
        1,
        3,
    ]

    assert list(
        df["new_atom_index"]
    ) == [
        1,
        2,
    ]

    assert count == 2

    assert mapping == {
        1: 1,
        3: 2,
    }

    assert legacy_mode is True

    assert "10" in section
    assert "30" in section

    assert "20" not in section
    assert "40" not in section


def test_modify_types_numeric_atom_types_are_ints():
    lines = [
        "1 12 # c1\n",
        "2 7 # hc\n",
        "\n",
    ]

    (
        df,
        _,
        _,
        _,
        legacy_mode,
    ) = modify_types(
        lines,
        [
            1,
            2,
        ],
        0,
    )

    assert df.loc[
        0,
        "atom_type",
    ] == 12

    assert df.loc[
        1,
        "atom_type",
    ] == 7

    assert legacy_mode is True


def test_modify_types_string_atom_type_disables_legacy_mode():
    lines = [
        "1 c1 # c1\n",
        "2 hc # hc\n",
        "\n",
    ]

    (
        df,
        _,
        _,
        _,
        legacy_mode,
    ) = modify_types(
        lines,
        [
            1,
            2,
        ],
        0,
    )

    assert (
        df.loc[
            0,
            "atom_type",
        ]
        == "c1"
    )

    assert (
        df.loc[
            1,
            "atom_type",
        ]
        == "hc"
    )

    assert legacy_mode is False


def test_modify_types_one_string_type_switches_global_legacy_mode_false():
    lines = [
        "1 10 # C\n",
        "2 hc # H\n",
        "\n",
    ]

    (
        df,
        _,
        _,
        _,
        legacy_mode,
    ) = modify_types(
        lines,
        [
            1,
            2,
        ],
        0,
    )

    # Characterizes current mixed parsing behavior:
    # first row remains integer, second row is string,
    # while legacy_mode becomes False globally.
    assert df.loc[
        0,
        "atom_type",
    ] == 10

    assert (
        df.loc[
            1,
            "atom_type",
        ]
        == "hc"
    )

    assert legacy_mode is False


def test_modify_types_stops_at_first_blank_line():
    lines = [
        "1 10 # C\n",
        "2 20 # H\n",
        "\n",
        "3 30 # O\n",
    ]

    (
        df,
        _,
        count,
        _,
        _,
    ) = modify_types(
        lines,
        [
            1,
            2,
            3,
        ],
        0,
    )

    assert count == 2

    assert list(
        df["atom_index"]
    ) == [
        1,
        2,
    ]


def test_modify_types_respects_start_index():
    lines = [
        "Types\n",
        "1 10 # C\n",
        "2 20 # H\n",
        "\n",
    ]

    (
        df,
        _,
        count,
        _,
        _,
    ) = modify_types(
        lines,
        [
            1,
            2,
        ],
        1,
    )

    assert count == 2

    assert list(
        df["atom_index"]
    ) == [
        1,
        2,
    ]


def test_modify_types_ignores_lines_with_fewer_than_four_fields():
    lines = [
        "this is short\n",
        "1 10 # C\n",
        "\n",
    ]

    (
        df,
        _,
        count,
        _,
        _,
    ) = modify_types(
        lines,
        [
            1,
        ],
        0,
    )

    assert count == 1

    assert list(
        df["atom_index"]
    ) == [
        1,
    ]


def test_modify_types_template_indexes_are_sorted_in_place():
    lines = [
        "1 10 # C\n",
        "2 20 # H\n",
        "3 30 # O\n",
        "\n",
    ]

    template_indexes = [
        3,
        1,
        2,
    ]

    modify_types(
        lines,
        template_indexes,
        0,
    )

    # Characterization:
    # production currently mutates the supplied list.
    assert template_indexes == [
        1,
        2,
        3,
    ]


def test_modify_types_output_order_follows_input_rows_not_template_list():
    lines = [
        "3 30 # O\n",
        "1 10 # C\n",
        "2 20 # H\n",
        "\n",
    ]

    template_indexes = [
        2,
        1,
        3,
    ]

    (
        df,
        _,
        _,
        mapping,
        _,
    ) = modify_types(
        lines,
        template_indexes,
        0,
    )

    # Even though template_indexes gets sorted,
    # surviving rows retain their order from the source section.
    assert list(
        df["atom_index"]
    ) == [
        3,
        1,
        2,
    ]

    assert mapping == {
        3: 1,
        1: 2,
        2: 3,
    }


def test_modify_types_all_selected():
    lines = [
        "10 4 # ca\n",
        "20 8 # hc\n",
        "\n",
    ]

    (
        df,
        _,
        count,
        mapping,
        _,
    ) = modify_types(
        lines,
        [
            10,
            20,
        ],
        0,
    )

    assert count == 2

    assert mapping == {
        10: 1,
        20: 2,
    }

    assert list(
        df["new_atom_index"]
    ) == [
        1,
        2,
    ]


def test_modify_types_none_selected_returns_empty_filtered_dataframe():
    lines = [
        "1 10 # C\n",
        "2 20 # H\n",
        "\n",
    ]

    (
        df,
        section,
        count,
        mapping,
        legacy_mode,
    ) = modify_types(
        lines,
        [
            100,
        ],
        0,
    )

    assert df.empty

    assert section == ""

    assert count == 0

    assert mapping == {}

    assert legacy_mode is True


def test_modify_types_preserves_hash_and_real_type():
    lines = [
        "5 12 HASHVALUE c2\n",
        "\n",
    ]

    (
        df,
        section,
        _,
        _,
        _,
    ) = modify_types(
        lines,
        [
            5,
        ],
        0,
    )

    assert (
        df.loc[
            0,
            "hash",
        ]
        == "HASHVALUE"
    )

    assert (
        df.loc[
            0,
            "atom_type_real",
        ]
        == "c2"
    )

    assert "HASHVALUE" in section
    assert "c2" in section


# =============================================================================
# modify_charges
# =============================================================================


def test_modify_charges_filters_and_reindexes_atoms():
    type_df = make_type_df(
        {
            2: 1,
            4: 2,
        }
    )

    lines = [
        "1 -0.1 # c1\n",
        "2 0.2 # c2\n",
        "3 -0.3 # h\n",
        "4 0.4 # o\n",
        "\n",
    ]

    section = modify_charges(
        lines,
        type_df,
        0,
    )

    assert "0.200000" in section
    assert "0.400000" in section

    assert "-0.100000" not in section
    assert "-0.300000" not in section

    output_lines = (
        section.strip().splitlines()
    )

    assert len(
        output_lines
    ) == 2

    assert (
        output_lines[0]
        .split()[0]
        == "1"
    )

    assert (
        output_lines[1]
        .split()[0]
        == "2"
    )


def test_modify_charges_formats_six_decimal_places():
    type_df = make_type_df(
        {
            1: 1,
        }
    )

    section = modify_charges(
        [
            "1 0.123456789 # c1\n",
            "\n",
        ],
        type_df,
        0,
    )

    assert (
        "0.123457"
        in section
    )


def test_modify_charges_stops_at_blank_line():
    type_df = make_type_df(
        {
            1: 1,
            2: 2,
        }
    )

    section = modify_charges(
        [
            "1 0.1 # c1\n",
            "\n",
            "2 0.2 # c2\n",
        ],
        type_df,
        0,
    )

    assert "0.100000" in section
    assert "0.200000" not in section


def test_modify_charges_respects_start_index():
    type_df = make_type_df(
        {
            1: 1,
        }
    )

    lines = [
        "Charges\n",
        "1 -0.25 # c1\n",
        "\n",
    ]

    section = modify_charges(
        lines,
        type_df,
        1,
    )

    assert "-0.250000" in section


def test_modify_charges_no_matching_atoms_returns_empty_string():
    type_df = make_type_df(
        {
            10: 1,
        }
    )

    section = modify_charges(
        [
            "1 0.1 # c1\n",
            "2 0.2 # c2\n",
            "\n",
        ],
        type_df,
        0,
    )

    assert section == ""


def test_modify_charges_preserves_hash_and_real_type():
    type_df = make_type_df(
        {
            1: 1,
        }
    )

    section = modify_charges(
        [
            "1 -0.5 HASH c2\n",
            "\n",
        ],
        type_df,
        0,
    )

    assert "HASH" in section
    assert "c2" in section


# =============================================================================
# modify_coords
# =============================================================================


def test_modify_coords_filters_and_reindexes():
    type_df = make_type_df(
        {
            2: 1,
            3: 2,
        }
    )

    lines = [
        "1 1.0 2.0 3.0 # c1\n",
        "2 4.0 5.0 6.0 # c2\n",
        "3 7.0 8.0 9.0 # o\n",
        "\n",
    ]

    section = modify_coords(
        lines,
        type_df,
        0,
    )

    assert "4.000000" in section
    assert "5.000000" in section
    assert "6.000000" in section

    assert "7.000000" in section
    assert "8.000000" in section
    assert "9.000000" in section

    assert "1.000000" not in section


def test_modify_coords_formats_six_decimal_places():
    type_df = make_type_df(
        {
            1: 1,
        }
    )

    section = modify_coords(
        [
            "1 1.1234567 -2.3456789 3.1 # c1\n",
            "\n",
        ],
        type_df,
        0,
    )

    assert "1.123457" in section
    assert "-2.345679" in section
    assert "3.100000" in section


def test_modify_coords_no_matching_atoms_returns_empty_string():
    type_df = make_type_df(
        {
            99: 1,
        }
    )

    section = modify_coords(
        [
            "1 1.0 2.0 3.0 # c\n",
            "\n",
        ],
        type_df,
        0,
    )

    assert section == ""


def test_modify_coords_stops_at_blank_line():
    type_df = make_type_df(
        {
            1: 1,
            2: 2,
        }
    )

    section = modify_coords(
        [
            "1 1.0 2.0 3.0 # c\n",
            "\n",
            "2 4.0 5.0 6.0 # o\n",
        ],
        type_df,
        0,
    )

    assert "1.000000" in section
    assert "4.000000" not in section


# =============================================================================
# modify_bonds
# =============================================================================


def test_modify_bonds_keeps_bond_when_both_atoms_survive():
    type_df = make_type_df(
        {
            5: 1,
            9: 2,
        }
    )

    section, count = modify_bonds(
        [
            "7 12 5 9 # c1 c2\n",
            "\n",
        ],
        type_df,
        0,
        True,
    )

    assert count == 1

    parts = section.split()

    assert parts[:4] == [
        "1",
        "12",
        "1",
        "2",
    ]


def test_modify_bonds_drops_bond_if_first_atom_missing():
    type_df = make_type_df(
        {
            2: 1,
        }
    )

    section, count = modify_bonds(
        [
            "1 5 1 2 # c1 c2\n",
            "\n",
        ],
        type_df,
        0,
        True,
    )

    assert section == ""
    assert count == 0


def test_modify_bonds_drops_bond_if_second_atom_missing():
    type_df = make_type_df(
        {
            1: 1,
        }
    )

    section, count = modify_bonds(
        [
            "1 5 1 2 # c1 c2\n",
            "\n",
        ],
        type_df,
        0,
        True,
    )

    assert section == ""
    assert count == 0


def test_modify_bonds_reindexes_bond_ids_sequentially():
    type_df = make_type_df(
        {
            1: 1,
            2: 2,
            3: 3,
        }
    )

    section, count = modify_bonds(
        [
            "50 5 1 2 # c1 c2\n",
            "99 6 2 3 # c2 c3\n",
            "\n",
        ],
        type_df,
        0,
        True,
    )

    rows = section.strip().splitlines()

    assert count == 2

    assert rows[0].split()[0] == "1"
    assert rows[1].split()[0] == "2"


def test_modify_bonds_legacy_mode_uses_integer_type():
    type_df = make_type_df(
        {
            1: 1,
            2: 2,
        }
    )

    section, count = modify_bonds(
        [
            "1 25 1 2 # c1 c2\n",
            "\n",
        ],
        type_df,
        0,
        True,
    )

    assert count == 1
    assert section.split()[1] == "25"


def test_modify_bonds_nonlegacy_mode_accepts_string_type():
    type_df = make_type_df(
        {
            1: 1,
            2: 2,
        }
    )

    section, count = modify_bonds(
        [
            "1 c1-c2 1 2 # c1 c2\n",
            "\n",
        ],
        type_df,
        0,
        False,
    )

    assert count == 1

    assert (
        section.split()[1]
        == "c1-c2"
    )


def test_modify_bonds_preserves_atom_mapping_from_type_df():
    type_df = make_type_df(
        {
            10: 2,
            20: 1,
        }
    )

    section, _ = modify_bonds(
        [
            "1 5 10 20 # a b\n",
            "\n",
        ],
        type_df,
        0,
        True,
    )

    parts = section.split()

    assert parts[2:4] == [
        "2",
        "1",
    ]


# =============================================================================
# modify_angles
# =============================================================================


def test_modify_angles_keeps_complete_angle():
    type_df = make_type_df(
        {
            10: 1,
            20: 2,
            30: 3,
        }
    )

    section, count = modify_angles(
        [
            "9 7 10 20 30 # c1 c2 c3\n",
            "\n",
        ],
        type_df,
        0,
        True,
    )

    assert count == 1

    parts = section.split()

    assert parts[:5] == [
        "1",
        "7",
        "1",
        "2",
        "3",
    ]


def test_modify_angles_drops_angle_when_any_atom_missing():
    type_df = make_type_df(
        {
            10: 1,
            20: 2,
        }
    )

    section, count = modify_angles(
        [
            "1 7 10 20 30 # c1 c2 c3\n",
            "\n",
        ],
        type_df,
        0,
        True,
    )

    assert section == ""
    assert count == 0


def test_modify_angles_reindexes_ids_sequentially():
    type_df = make_type_df(
        {
            1: 1,
            2: 2,
            3: 3,
            4: 4,
        }
    )

    section, count = modify_angles(
        [
            "50 4 1 2 3 # a b c\n",
            "80 5 2 3 4 # b c d\n",
            "\n",
        ],
        type_df,
        0,
        True,
    )

    rows = section.strip().splitlines()

    assert count == 2

    assert rows[0].split()[0] == "1"
    assert rows[1].split()[0] == "2"


def test_modify_angles_nonlegacy_accepts_string_type():
    type_df = make_type_df(
        {
            1: 1,
            2: 2,
            3: 3,
        }
    )

    section, count = modify_angles(
        [
            "1 c1-c2-c3 1 2 3 # c1 c2 c3\n",
            "\n",
        ],
        type_df,
        0,
        False,
    )

    assert count == 1

    assert (
        section.split()[1]
        == "c1-c2-c3"
    )


# =============================================================================
# modify_dihedrals
# =============================================================================


def test_modify_dihedrals_keeps_complete_dihedral():
    type_df = make_type_df(
        {
            10: 1,
            20: 2,
            30: 3,
            40: 4,
        }
    )

    section, count = modify_dihedrals(
        [
            "15 8 10 20 30 40 # a b c d\n",
            "\n",
        ],
        type_df,
        0,
        True,
    )

    assert count == 1

    parts = section.split()

    assert parts[:6] == [
        "1",
        "8",
        "1",
        "2",
        "3",
        "4",
    ]


def test_modify_dihedrals_drops_when_any_atom_missing():
    type_df = make_type_df(
        {
            10: 1,
            20: 2,
            30: 3,
        }
    )

    section, count = modify_dihedrals(
        [
            "1 8 10 20 30 40 # a b c d\n",
            "\n",
        ],
        type_df,
        0,
        True,
    )

    assert section == ""
    assert count == 0


def test_modify_dihedrals_reindexes_ids():
    type_df = make_type_df(
        {
            1: 1,
            2: 2,
            3: 3,
            4: 4,
            5: 5,
        }
    )

    section, count = modify_dihedrals(
        [
            "90 8 1 2 3 4 # a b c d\n",
            "92 9 2 3 4 5 # b c d e\n",
            "\n",
        ],
        type_df,
        0,
        True,
    )

    rows = section.strip().splitlines()

    assert count == 2

    assert rows[0].split()[0] == "1"
    assert rows[1].split()[0] == "2"


def test_modify_dihedrals_nonlegacy_accepts_string_type():
    type_df = make_type_df(
        {
            1: 1,
            2: 2,
            3: 3,
            4: 4,
        }
    )

    section, count = modify_dihedrals(
        [
            "1 torsion_name 1 2 3 4 # a b c d\n",
            "\n",
        ],
        type_df,
        0,
        False,
    )

    assert count == 1

    assert (
        section.split()[1]
        == "torsion_name"
    )


# =============================================================================
# modify_impropers
# =============================================================================


def test_modify_impropers_keeps_complete_improper():
    type_df = make_type_df(
        {
            10: 1,
            20: 2,
            30: 3,
            40: 4,
        }
    )

    section, count = modify_impropers(
        [
            "25 11 10 20 30 40 # a b c d\n",
            "\n",
        ],
        type_df,
        0,
        True,
    )

    assert count == 1

    parts = section.split()

    assert parts[:6] == [
        "1",
        "11",
        "1",
        "2",
        "3",
        "4",
    ]


def test_modify_impropers_drops_when_any_atom_missing():
    type_df = make_type_df(
        {
            10: 1,
            20: 2,
            30: 3,
        }
    )

    section, count = modify_impropers(
        [
            "1 11 10 20 30 40 # a b c d\n",
            "\n",
        ],
        type_df,
        0,
        True,
    )

    assert section == ""
    assert count == 0


def test_modify_impropers_reindexes_ids_sequentially():
    type_df = make_type_df(
        {
            1: 1,
            2: 2,
            3: 3,
            4: 4,
            5: 5,
        }
    )

    section, count = modify_impropers(
        [
            "40 11 1 2 3 4 # a b c d\n",
            "70 12 2 3 4 5 # b c d e\n",
            "\n",
        ],
        type_df,
        0,
        True,
    )

    rows = section.strip().splitlines()

    assert count == 2

    assert rows[0].split()[0] == "1"
    assert rows[1].split()[0] == "2"


def test_modify_impropers_nonlegacy_accepts_string_type():
    type_df = make_type_df(
        {
            1: 1,
            2: 2,
            3: 3,
            4: 4,
        }
    )

    section, count = modify_impropers(
        [
            "1 improper_name 1 2 3 4 # a b c d\n",
            "\n",
        ],
        type_df,
        0,
        False,
    )

    assert count == 1

    assert (
        section.split()[1]
        == "improper_name"
    )


# =============================================================================
# Filtering invariants across topology types
# =============================================================================


@pytest.mark.parametrize(
    "function,line",
    [
        (
            modify_bonds,
            "1 1 1 99 # a z\n",
        ),
        (
            modify_angles,
            "1 1 1 2 99 # a b z\n",
        ),
        (
            modify_dihedrals,
            "1 1 1 2 3 99 # a b c z\n",
        ),
        (
            modify_impropers,
            "1 1 1 2 3 99 # a b c z\n",
        ),
    ],
)
def test_topology_entry_removed_when_it_references_atom_outside_template(
    function,
    line,
):
    type_df = make_type_df(
        {
            1: 1,
            2: 2,
            3: 3,
        }
    )

    section, count = function(
        [
            line,
            "\n",
        ],
        type_df,
        0,
        True,
    )

    assert section == ""
    assert count == 0


@pytest.mark.parametrize(
    "function,line",
    [
        (
            modify_bonds,
            "1 10 5 10 # a b\n",
        ),
        (
            modify_angles,
            "1 10 5 10 15 # a b c\n",
        ),
        (
            modify_dihedrals,
            "1 10 5 10 15 20 # a b c d\n",
        ),
        (
            modify_impropers,
            "1 10 5 10 15 20 # a b c d\n",
        ),
    ],
)
def test_topology_functions_preserve_nontrivial_atom_remapping(
    function,
    line,
):
    type_df = make_type_df(
        {
            5: 4,
            10: 3,
            15: 2,
            20: 1,
        }
    )

    section, count = function(
        [
            line,
            "\n",
        ],
        type_df,
        0,
        True,
    )

    assert count == 1

    parts = section.split()

    if function is modify_bonds:
        assert parts[2:4] == [
            "4",
            "3",
        ]

    elif function is modify_angles:
        assert parts[2:5] == [
            "4",
            "3",
            "2",
        ]

    else:
        assert parts[2:6] == [
            "4",
            "3",
            "2",
            "1",
        ]


# =============================================================================
# Malformed / short rows
# =============================================================================


def test_modify_bonds_ignores_short_rows():
    type_df = make_type_df(
        {
            1: 1,
            2: 2,
        }
    )

    section, count = modify_bonds(
        [
            "1 2 1\n",
            "2 3 1 2 # a b\n",
            "\n",
        ],
        type_df,
        0,
        True,
    )

    assert count == 1
    assert section.split()[0] == "1"


def test_modify_angles_ignores_short_rows():
    type_df = make_type_df(
        {
            1: 1,
            2: 2,
            3: 3,
        }
    )

    section, count = modify_angles(
        [
            "1 2 1 2\n",
            "2 3 1 2 3 # a b c\n",
            "\n",
        ],
        type_df,
        0,
        True,
    )

    assert count == 1


def test_modify_dihedrals_ignores_short_rows():
    type_df = make_type_df(
        {
            1: 1,
            2: 2,
            3: 3,
            4: 4,
        }
    )

    section, count = modify_dihedrals(
        [
            "1 2 1 2 3\n",
            "2 3 1 2 3 4 # a b c d\n",
            "\n",
        ],
        type_df,
        0,
        True,
    )

    assert count == 1


def test_modify_impropers_ignores_short_rows():
    type_df = make_type_df(
        {
            1: 1,
            2: 2,
            3: 3,
            4: 4,
        }
    )

    section, count = modify_impropers(
        [
            "1 2 1 2 3\n",
            "2 3 1 2 3 4 # a b c d\n",
            "\n",
        ],
        type_df,
        0,
        True,
    )

    assert count == 1


# =============================================================================
# Cross-function integration
# =============================================================================


def test_modify_types_mapping_can_drive_all_topology_modifiers():
    type_lines = [
        "10 5 # c1\n",
        "20 6 # c2\n",
        "30 7 # c3\n",
        "40 8 # c4\n",
        "\n",
    ]

    (
        type_df,
        _,
        number_of_types,
        index_mapping,
        legacy_mode,
    ) = modify_types(
        type_lines,
        [
            10,
            20,
            30,
            40,
        ],
        0,
    )

    assert number_of_types == 4

    assert index_mapping == {
        10: 1,
        20: 2,
        30: 3,
        40: 4,
    }

    assert legacy_mode is True

    bond_section, bonds = modify_bonds(
        [
            "1 1 10 20 # c1 c2\n",
            "\n",
        ],
        type_df,
        0,
        legacy_mode,
    )

    angle_section, angles = modify_angles(
        [
            "1 1 10 20 30 # c1 c2 c3\n",
            "\n",
        ],
        type_df,
        0,
        legacy_mode,
    )

    dihedral_section, dihedrals = (
        modify_dihedrals(
            [
                "1 1 10 20 30 40 # c1 c2 c3 c4\n",
                "\n",
            ],
            type_df,
            0,
            legacy_mode,
        )
    )

    improper_section, impropers = (
        modify_impropers(
            [
                "1 1 10 20 30 40 # c1 c2 c3 c4\n",
                "\n",
            ],
            type_df,
            0,
            legacy_mode,
        )
    )

    assert bonds == 1
    assert angles == 1
    assert dihedrals == 1
    assert impropers == 1

    assert bond_section.split()[2:4] == [
        "1",
        "2",
    ]

    assert angle_section.split()[2:5] == [
        "1",
        "2",
        "3",
    ]

    assert dihedral_section.split()[2:6] == [
        "1",
        "2",
        "3",
        "4",
    ]

    assert improper_section.split()[2:6] == [
        "1",
        "2",
        "3",
        "4",
    ]


def test_string_type_mode_propagates_to_topology_modifiers():
    type_lines = [
        "1 c1 # c1\n",
        "2 c2 # c2\n",
        "3 c3 # c3\n",
        "4 c4 # c4\n",
        "\n",
    ]

    (
        type_df,
        _,
        _,
        _,
        legacy_mode,
    ) = modify_types(
        type_lines,
        [
            1,
            2,
            3,
            4,
        ],
        0,
    )

    assert legacy_mode is False

    bond_section, _ = modify_bonds(
        [
            "1 c1-c2 1 2 # c1 c2\n",
            "\n",
        ],
        type_df,
        0,
        legacy_mode,
    )

    angle_section, _ = modify_angles(
        [
            "1 c1-c2-c3 1 2 3 # c1 c2 c3\n",
            "\n",
        ],
        type_df,
        0,
        legacy_mode,
    )

    dihedral_section, _ = modify_dihedrals(
        [
            "1 c1-c2-c3-c4 1 2 3 4 # c1 c2 c3 c4\n",
            "\n",
        ],
        type_df,
        0,
        legacy_mode,
    )

    improper_section, _ = modify_impropers(
        [
            "1 improper-c1 1 2 3 4 # c1 c2 c3 c4\n",
            "\n",
        ],
        type_df,
        0,
        legacy_mode,
    )

    assert (
        bond_section.split()[1]
        == "c1-c2"
    )

    assert (
        angle_section.split()[1]
        == "c1-c2-c3"
    )

    assert (
        dihedral_section.split()[1]
        == "c1-c2-c3-c4"
    )

    assert (
        improper_section.split()[1]
        == "improper-c1"
    )
