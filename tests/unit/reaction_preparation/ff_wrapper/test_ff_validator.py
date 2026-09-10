from pathlib import Path
from types import SimpleNamespace

import pytest

import AutoREACTER.reaction_preparation.ff_wrapper.ff_validator as ff_validator
from AutoREACTER.reaction_preparation.ff_wrapper.ff_validator import (
    FFValidator,
    ForceFieldValidationError,
)


# =============================================================================
# Helpers
# =============================================================================


def make_ff_files(
    ff_file: Path,
):
    return SimpleNamespace(
        force_field_data=ff_file,
    )


def write_force_field(
    path: Path,
    text: str,
) -> Path:
    path.write_text(
        text,
        encoding="utf-8",
    )

    return path


def make_uninitialized_validator(
    ff_file=None,
):
    """
    Construct FFValidator without calling __init__/validate.

    Useful for testing the helper methods independently.
    """
    validator = object.__new__(
        FFValidator
    )

    validator.ff_file = ff_file

    if ff_file is not None:
        validator.ff_files = (
            make_ff_files(
                ff_file
            )
        )

    return validator


# =============================================================================
# Exception type
# =============================================================================


def test_force_field_validation_error_is_exception():
    assert issubclass(
        ForceFieldValidationError,
        Exception,
    )


# =============================================================================
# Constructor
# =============================================================================


def test_constructor_stores_ff_files_and_force_field_path(
    tmp_path,
    monkeypatch,
):
    ff_file = (
        tmp_path / "force_field.data"
    )

    ff_files = make_ff_files(
        ff_file
    )

    calls = []

    monkeypatch.setattr(
        FFValidator,
        "validate",
        lambda self:
        (
            calls.append(self)
            or True
        ),
    )

    validator = FFValidator(
        ff_files
    )

    assert (
        validator.ff_files
        is ff_files
    )

    assert (
        validator.ff_file
        == ff_file
    )

    assert calls == [
        validator
    ]


def test_constructor_validates_immediately(
    tmp_path,
):
    missing = (
        tmp_path
        / "missing.data"
    )

    with pytest.raises(
        ForceFieldValidationError,
        match="File not found",
    ):
        FFValidator(
            make_ff_files(
                missing
            )
        )


# =============================================================================
# validate - file existence
# =============================================================================


def test_validate_missing_file_raises(
    tmp_path,
):
    missing = (
        tmp_path
        / "missing.data"
    )

    validator = (
        make_uninitialized_validator(
            missing
        )
    )

    with pytest.raises(
        ForceFieldValidationError,
        match="File not found",
    ) as exc_info:
        validator.validate()

    assert str(
        missing
    ) in str(
        exc_info.value
    )


def test_validate_existing_valid_file_returns_true(
    tmp_path,
):
    ff_file = write_force_field(
        tmp_path / "force_field.data",
        """
Pair Coeffs

1 0.100 3.500
2 0.200 4.000

""",
    )

    validator = FFValidator(
        make_ff_files(
            ff_file
        )
    )

    assert (
        validator.validate()
        is True
    )


# =============================================================================
# validate - current section-presence behavior
# =============================================================================


def test_validate_empty_file_currently_returns_true(
    tmp_path,
):
    """
    Characterization test.

    Although the class documentation says required section labels
    are checked, the current implementation skips any section that
    is absent.
    """
    ff_file = write_force_field(
        tmp_path / "empty.data",
        "",
    )

    validator = FFValidator(
        make_ff_files(
            ff_file
        )
    )

    assert (
        validator.validate()
        is True
    )


def test_validate_file_with_no_recognized_sections_currently_returns_true(
    tmp_path,
):
    """
    Characterizes the current implementation:
    unrecognized/missing coefficient sections are ignored.
    """
    ff_file = write_force_field(
        tmp_path / "other.data",
        """
Masses

1 12.011
2 1.008

Atoms

1 1 1 0.0 0.0 0.0 0.0
""",
    )

    validator = FFValidator(
        make_ff_files(
            ff_file
        )
    )

    assert (
        validator.validate()
        is True
    )


def test_validate_missing_some_coeff_sections_are_skipped(
    tmp_path,
):
    ff_file = write_force_field(
        tmp_path / "partial.data",
        """
Bond Coeffs

1 100.0 1.5
2 200.0 1.4

""",
    )

    validator = FFValidator(
        make_ff_files(
            ff_file
        )
    )

    assert (
        validator.validate()
        is True
    )


@pytest.mark.parametrize(
    "section_name",
    [
        "Pair Coeffs",
        "Bond Coeffs",
        "Angle Coeffs",
        "Dihedral Coeffs",
    ],
)
def test_validate_recognizes_each_supported_section(
    tmp_path,
    section_name,
):
    ff_file = write_force_field(
        tmp_path / "force_field.data",
        f"""
{section_name}

1 1.0 2.0
2 3.0 4.0

""",
    )

    validator = FFValidator(
        make_ff_files(
            ff_file
        )
    )

    assert (
        validator.validate()
        is True
    )


# =============================================================================
# find_section_start
# =============================================================================


def test_find_section_start_returns_correct_index():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "LAMMPS data file\n",
        "\n",
        "Masses\n",
        "\n",
        "Pair Coeffs\n",
        "\n",
    ]

    assert (
        validator.find_section_start(
            lines,
            "Pair Coeffs",
        )
        == 4
    )


def test_find_section_start_returns_none_when_absent():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "Masses\n",
        "Atoms\n",
        "Bonds\n",
    ]

    assert (
        validator.find_section_start(
            lines,
            "Pair Coeffs",
        )
        is None
    )


def test_find_section_start_uses_substring_matching():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "Pair Coeffs  # class2\n",
    ]

    assert (
        validator.find_section_start(
            lines,
            "Pair Coeffs",
        )
        == 0
    )


def test_find_section_start_returns_first_match():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "Pair Coeffs\n",
        "\n",
        "Pair Coeffs\n",
    ]

    assert (
        validator.find_section_start(
            lines,
            "Pair Coeffs",
        )
        == 0
    )


def test_find_section_start_is_case_sensitive():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "pair coeffs\n",
    ]

    assert (
        validator.find_section_start(
            lines,
            "Pair Coeffs",
        )
        is None
    )


# =============================================================================
# find_data_block
# =============================================================================


def test_find_data_block_skips_blank_lines_after_header():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "Pair Coeffs\n",
        "\n",
        "\n",
        "1 0.1 3.5\n",
        "2 0.2 4.0\n",
        "\n",
        "Bond Coeffs\n",
    ]

    result = (
        validator.find_data_block(
            lines,
            0,
        )
    )

    assert result == (
        3,
        5,
    )


def test_find_data_block_stops_at_first_blank_line():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "Pair Coeffs\n",
        "\n",
        "1 0.1 3.5\n",
        "2 0.2 4.0\n",
        "\n",
        "3 0.3 5.0\n",
    ]

    result = (
        validator.find_data_block(
            lines,
            0,
        )
    )

    assert result == (
        2,
        4,
    )


def test_find_data_block_can_end_at_eof():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "Bond Coeffs\n",
        "\n",
        "1 100.0 1.5\n",
        "2 200.0 1.4\n",
    ]

    result = (
        validator.find_data_block(
            lines,
            0,
        )
    )

    assert result == (
        2,
        4,
    )


def test_find_data_block_empty_after_header_returns_eof_indexes():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "Angle Coeffs\n",
        "\n",
        "\n",
    ]

    result = (
        validator.find_data_block(
            lines,
            0,
        )
    )

    assert result == (
        3,
        3,
    )


# =============================================================================
# check_coefficients - valid input
# =============================================================================


def test_check_coefficients_single_valid_row():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "1 0.100 3.500\n",
    ]

    assert (
        validator.check_coefficients(
            start_line=0,
            end_line=1,
            lines=lines,
            section_name="Pair Coeffs",
        )
        is True
    )


def test_check_coefficients_multiple_rows_same_count():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "1 0.100 3.500\n",
        "2 0.200 4.000\n",
        "3 0.300 4.500\n",
    ]

    assert (
        validator.check_coefficients(
            start_line=0,
            end_line=3,
            lines=lines,
            section_name="Pair Coeffs",
        )
        is True
    )


def test_check_coefficients_accepts_negative_values():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "1 -1.25 2.50 -0.75\n",
        "2 1.25 -2.50 0.75\n",
    ]

    assert (
        validator.check_coefficients(
            0,
            2,
            lines,
            "Dihedral Coeffs",
        )
        is True
    )


def test_check_coefficients_accepts_scientific_notation():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "1 1.0e-4 2.5E+2\n",
        "2 -3.0e-5 4.0E-1\n",
    ]

    assert (
        validator.check_coefficients(
            0,
            2,
            lines,
            "Pair Coeffs",
        )
        is True
    )


def test_check_coefficients_accepts_zero_values_when_row_not_all_zero():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "1 0.0 0.0 1.5 0.0\n",
    ]

    assert (
        validator.check_coefficients(
            0,
            1,
            lines,
            "Dihedral Coeffs",
        )
        is True
    )


def test_check_coefficients_ignores_inline_comments():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "1 0.100 3.500 # carbon\n",
        "2 0.200 4.000 # hydrogen\n",
    ]

    assert (
        validator.check_coefficients(
            0,
            2,
            lines,
            "Pair Coeffs",
        )
        is True
    )


def test_check_coefficients_ignores_full_comment_lines_inside_block():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "# comment only\n",
        "1 10.0 1.5\n",
        "# another comment\n",
        "2 20.0 1.4\n",
    ]

    assert (
        validator.check_coefficients(
            0,
            4,
            lines,
            "Bond Coeffs",
        )
        is True
    )


def test_check_coefficients_currently_does_not_validate_type_identifier():
    """
    Characterization test.

    values[0] is treated as the type identifier, but current production
    code only converts values[1:] to floats.
    """
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "not_an_integer 1.0 2.0\n",
        "another_id 3.0 4.0\n",
    ]

    assert (
        validator.check_coefficients(
            0,
            2,
            lines,
            "Pair Coeffs",
        )
        is True
    )


# =============================================================================
# check_coefficients - dynamic coefficient count
# =============================================================================


def test_first_data_row_sets_expected_coefficient_count():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "1 1.0 2.0 3.0 4.0\n",
        "2 5.0 6.0 7.0 8.0\n",
    ]

    assert (
        validator.check_coefficients(
            0,
            2,
            lines,
            "Dihedral Coeffs",
        )
        is True
    )


def test_inconsistent_coefficient_count_raises():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "1 1.0 2.0 3.0\n",
        "2 4.0 5.0\n",
    ]

    with pytest.raises(
        ForceFieldValidationError,
        match="Inconsistent data",
    ) as exc_info:
        validator.check_coefficients(
            0,
            2,
            lines,
            "Angle Coeffs",
        )

    message = str(
        exc_info.value
    )

    assert (
        "Expected 3 coefficients"
        in message
    )

    assert (
        "found 2"
        in message
    )


def test_inconsistent_coefficient_error_includes_section_name():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "1 1.0 2.0\n",
        "2 3.0\n",
    ]

    with pytest.raises(
        ForceFieldValidationError
    ) as exc_info:
        validator.check_coefficients(
            0,
            2,
            lines,
            "Bond Coeffs",
        )

    assert (
        "Bond Coeffs"
        in str(
            exc_info.value
        )
    )


def test_comment_does_not_affect_dynamic_coefficient_count():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "# ignored comment\n",
        "1 1.0 2.0 3.0\n",
        "2 4.0 5.0 6.0\n",
    ]

    assert (
        validator.check_coefficients(
            0,
            3,
            lines,
            "Angle Coeffs",
        )
        is True
    )


# =============================================================================
# check_coefficients - non-numeric values
# =============================================================================


def test_non_numeric_coefficient_raises():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "1 1.0 not-a-number 3.0\n",
    ]

    with pytest.raises(
        ForceFieldValidationError,
        match="Non-numeric data",
    ):
        validator.check_coefficients(
            0,
            1,
            lines,
            "Angle Coeffs",
        )


def test_non_numeric_error_includes_section_name():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "1 1.0 BAD\n",
    ]

    with pytest.raises(
        ForceFieldValidationError
    ) as exc_info:
        validator.check_coefficients(
            0,
            1,
            lines,
            "Pair Coeffs",
        )

    assert (
        "Pair Coeffs"
        in str(
            exc_info.value
        )
    )


def test_non_numeric_error_reports_one_based_line_number():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "header\n",
        "ignored\n",
        "1 1.0 2.0\n",
        "2 BAD 4.0\n",
    ]

    with pytest.raises(
        ForceFieldValidationError
    ) as exc_info:
        validator.check_coefficients(
            start_line=2,
            end_line=4,
            lines=lines,
            section_name="Pair Coeffs",
        )

    assert (
        "line 4"
        in str(
            exc_info.value
        )
    )


# =============================================================================
# check_coefficients - all-zero rows
# =============================================================================


@pytest.mark.parametrize(
    "line",
    [
        "1 0.0\n",
        "1 0 0\n",
        "1 0.000000 0.000000 0.000000\n",
        "1 -0.0 0.0\n",
    ],
)
def test_all_zero_coefficients_raise(
    line,
):
    validator = (
        make_uninitialized_validator()
    )

    with pytest.raises(
        ForceFieldValidationError,
        match="All zeros",
    ):
        validator.check_coefficients(
            0,
            1,
            [
                line
            ],
            "Pair Coeffs",
        )


def test_all_zero_error_includes_section_name():
    validator = (
        make_uninitialized_validator()
    )

    with pytest.raises(
        ForceFieldValidationError
    ) as exc_info:
        validator.check_coefficients(
            0,
            1,
            [
                "1 0.0 0.0\n"
            ],
            "Dihedral Coeffs",
        )

    assert (
        "Dihedral Coeffs"
        in str(
            exc_info.value
        )
    )


def test_all_zero_error_includes_original_row():
    validator = (
        make_uninitialized_validator()
    )

    line = (
        "7 0.000 0.000 0.000 # test type\n"
    )

    with pytest.raises(
        ForceFieldValidationError
    ) as exc_info:
        validator.check_coefficients(
            0,
            1,
            [
                line
            ],
            "Angle Coeffs",
        )

    assert (
        "7 0.000 0.000 0.000"
        in str(
            exc_info.value
        )
    )


# =============================================================================
# check_coefficients - no usable data
# =============================================================================


def test_empty_data_block_raises_no_data():
    validator = (
        make_uninitialized_validator()
    )

    with pytest.raises(
        ForceFieldValidationError,
        match="No data found in section Pair Coeffs",
    ):
        validator.check_coefficients(
            0,
            0,
            [],
            "Pair Coeffs",
        )


def test_comment_only_block_raises_no_data():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "# comment one\n",
        "# comment two\n",
    ]

    with pytest.raises(
        ForceFieldValidationError,
        match="No data found in section Bond Coeffs",
    ):
        validator.check_coefficients(
            0,
            2,
            lines,
            "Bond Coeffs",
        )


def test_whitespace_and_comments_only_raise_no_data():
    validator = (
        make_uninitialized_validator()
    )

    lines = [
        "   # comment\n",
        "      \n",
        "# another\n",
    ]

    with pytest.raises(
        ForceFieldValidationError,
        match="No data found",
    ):
        validator.check_coefficients(
            0,
            3,
            lines,
            "Angle Coeffs",
        )


# =============================================================================
# validate - full section integration
# =============================================================================


def test_validate_checks_multiple_coefficient_sections(
    tmp_path,
):
    ff_file = write_force_field(
        tmp_path / "force_field.data",
        """
Pair Coeffs

1 0.10 3.50
2 0.20 4.00

Bond Coeffs

1 300.0 1.50
2 250.0 1.40

Angle Coeffs

1 50.0 109.5
2 60.0 120.0

Dihedral Coeffs

1 1.0 2.0 3.0
2 4.0 5.0 6.0

""",
    )

    validator = FFValidator(
        make_ff_files(
            ff_file
        )
    )

    assert (
        validator.validate()
        is True
    )


def test_validate_invalid_later_section_still_raises(
    tmp_path,
):
    ff_file = write_force_field(
        tmp_path / "force_field.data",
        """
Pair Coeffs

1 0.10 3.50
2 0.20 4.00

Bond Coeffs

1 300.0 1.50
2 0.0 0.0

""",
    )

    with pytest.raises(
        ForceFieldValidationError,
        match="All zeros",
    ):
        FFValidator(
            make_ff_files(
                ff_file
            )
        )


def test_validate_section_with_no_data_raises(
    tmp_path,
):
    ff_file = write_force_field(
        tmp_path / "force_field.data",
        """
Pair Coeffs


""",
    )

    with pytest.raises(
        ForceFieldValidationError,
        match="No data found in section Pair Coeffs",
    ):
        FFValidator(
            make_ff_files(
                ff_file
            )
        )


def test_validate_inconsistent_rows_in_file_raises(
    tmp_path,
):
    ff_file = write_force_field(
        tmp_path / "force_field.data",
        """
Angle Coeffs

1 50.0 109.5 2.0
2 60.0 120.0

""",
    )

    with pytest.raises(
        ForceFieldValidationError,
        match="Inconsistent data",
    ):
        FFValidator(
            make_ff_files(
                ff_file
            )
        )


def test_validate_non_numeric_file_data_raises(
    tmp_path,
):
    ff_file = write_force_field(
        tmp_path / "force_field.data",
        """
Bond Coeffs

1 300.0 1.50
2 WRONG 1.40

""",
    )

    with pytest.raises(
        ForceFieldValidationError,
        match="Non-numeric data",
    ):
        FFValidator(
            make_ff_files(
                ff_file
            )
        )


# =============================================================================
# Comments / LAMMPS-style headers
# =============================================================================


def test_validate_class2_section_header_comment(
    tmp_path,
):
    ff_file = write_force_field(
        tmp_path / "force_field.data",
        """
Dihedral Coeffs  # class2

1 0.0 0.0 0.0514 0.0 -0.143 0.0
2 0.0 0.0 0.0316 0.0 -0.1681 0.0

""",
    )

    validator = FFValidator(
        make_ff_files(
            ff_file
        )
    )

    assert (
        validator.validate()
        is True
    )


def test_validate_inline_atom_type_comments_do_not_count_as_coefficients(
    tmp_path,
):
    ff_file = write_force_field(
        tmp_path / "force_field.data",
        """
Pair Coeffs

1 0.0540 4.0100 # c1
2 0.0640 3.8540 # c2

""",
    )

    validator = FFValidator(
        make_ff_files(
            ff_file
        )
    )

    assert (
        validator.validate()
        is True
    )