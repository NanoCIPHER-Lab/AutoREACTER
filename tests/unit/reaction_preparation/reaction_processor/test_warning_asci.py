from AutoREACTER.reaction_preparation.reaction_processor.warning_asci import (
    ascii_art,
    print_warning,
)


# =============================================================================
# ascii_art
# =============================================================================


def test_ascii_art_uppercases_message(capsys):
    ascii_art(
        "reaction progression warning"
    )

    output = capsys.readouterr().out

    assert (
        "WARNING: REACTION PROGRESSION WARNING"
        in output
    )


def test_ascii_art_prints_warning_prefix(capsys):
    ascii_art(
        "test message"
    )

    output = capsys.readouterr().out

    assert output.startswith(
        "WARNING: TEST MESSAGE"
    )


def test_ascii_art_prints_banner(capsys):
    ascii_art(
        "test"
    )

    output = capsys.readouterr().out

    # Stable fragments from the ASCII banner.
    assert "____" in output
    assert "|_  _|" in output
    assert "WARNING:" in output


def test_ascii_art_prints_multiline_output(capsys):
    ascii_art(
        "hello"
    )

    output = capsys.readouterr().out

    nonempty_lines = [
        line
        for line in output.splitlines()
        if line.strip()
    ]

    # One warning line + several ASCII-art lines.
    assert len(nonempty_lines) >= 6


def test_ascii_art_handles_empty_message(capsys):
    ascii_art("")

    output = capsys.readouterr().out

    assert output.startswith(
        "WARNING:"
    )

    assert "____" in output


def test_ascii_art_does_not_return_value():
    result = ascii_art(
        "warning"
    )

    assert result is None


# =============================================================================
# print_warning
# =============================================================================


def test_print_warning_prints_beta_phase_message(
    capsys,
):
    print_warning()

    output = capsys.readouterr().out

    assert (
        "ENTERING THE REACTION PROGRESSION LOOP "
        "IS STILL IN THE BETA PHASE."
        in output
    )


def test_print_warning_mentions_chemical_accuracy(
    capsys,
):
    print_warning()

    output = capsys.readouterr().out

    assert (
        "RESULTS CAN BE CHEMICALLY INACCURATE"
        in output
    )


def test_print_warning_uses_warning_prefix(
    capsys,
):
    print_warning()

    output = capsys.readouterr().out

    assert output.startswith(
        "WARNING:"
    )


def test_print_warning_prints_ascii_banner(
    capsys,
):
    print_warning()

    output = capsys.readouterr().out

    assert "____" in output
    assert "|_  _|" in output


def test_print_warning_returns_none():
    result = print_warning()

    assert result is None