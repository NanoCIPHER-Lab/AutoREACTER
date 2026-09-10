import importlib

import pytest

from AutoREACTER.initialization import Initialization


# =============================================================================
# Constructor
# =============================================================================


def test_initialization_runs_all_steps_in_order(monkeypatch):
    calls = []

    monkeypatch.setattr(
        Initialization,
        "moldule_imports",
        classmethod(
            lambda cls: calls.append("modules")
        ),
    )

    monkeypatch.setattr(
        Initialization,
        "ASCII_Mupt_reaction_LAMMPS",
        classmethod(
            lambda cls: calls.append("ascii")
        ),
    )

    monkeypatch.setattr(
        Initialization,
        "print_version",
        classmethod(
            lambda cls: calls.append("version")
        ),
    )

    Initialization()

    assert calls == [
        "modules",
        "ascii",
        "version",
    ]


# =============================================================================
# Required module imports
# =============================================================================


def test_moldule_imports_imports_all_required_modules(
    monkeypatch,
    capsys,
):
    imported_modules = []

    def fake_import_module(name):
        imported_modules.append(name)
        return object()

    monkeypatch.setattr(
        importlib,
        "import_module",
        fake_import_module,
    )

    Initialization.moldule_imports()

    assert imported_modules == [
        "rdkit",
        "pandas",
        "numpy",
        "networkx",
    ]

    output = capsys.readouterr().out

    assert (
        "All required modules are successfully imported."
        in output
    )


def test_moldule_imports_stops_when_module_is_missing(
    monkeypatch,
):
    imported_modules = []

    def fake_import_module(name):
        imported_modules.append(name)

        if name == "numpy":
            error = ModuleNotFoundError(
                "No module named 'numpy'"
            )
            error.name = "numpy"
            raise error

        return object()

    monkeypatch.setattr(
        importlib,
        "import_module",
        fake_import_module,
    )

    with pytest.raises(RuntimeError):
        Initialization.moldule_imports()

    assert imported_modules == [
        "rdkit",
        "pandas",
        "numpy",
    ]


def test_moldule_imports_missing_module_error_mentions_module_name(
    monkeypatch,
):
    def fake_import_module(name):
        if name == "pandas":
            error = ModuleNotFoundError(
                "No module named 'pandas'"
            )
            error.name = "pandas"
            raise error

        return object()

    monkeypatch.setattr(
        importlib,
        "import_module",
        fake_import_module,
    )

    with pytest.raises(
        RuntimeError,
        match="pandas",
    ):
        Initialization.moldule_imports()


def test_moldule_imports_missing_module_error_has_helpful_message(
    monkeypatch,
):
    def fake_import_module(name):
        error = ModuleNotFoundError(
            f"No module named '{name}'"
        )
        error.name = name
        raise error

    monkeypatch.setattr(
        importlib,
        "import_module",
        fake_import_module,
    )

    with pytest.raises(RuntimeError) as exc_info:
        Initialization.moldule_imports()

    message = str(exc_info.value)

    assert "Required module not found" in message
    assert "Please install the missing module" in message
    assert "Exiting program" in message


def test_moldule_imports_does_not_print_success_when_import_fails(
    monkeypatch,
    capsys,
):
    def fake_import_module(name):
        error = ModuleNotFoundError(
            f"No module named '{name}'"
        )
        error.name = name
        raise error

    monkeypatch.setattr(
        importlib,
        "import_module",
        fake_import_module,
    )

    with pytest.raises(RuntimeError):
        Initialization.moldule_imports()

    output = capsys.readouterr().out

    assert (
        "All required modules are successfully imported."
        not in output
    )


# =============================================================================
# ASCII banner
# =============================================================================


def test_ascii_banner_prints_autoreacter_name(
    capsys,
):
    Initialization.ASCII_Mupt_reaction_LAMMPS()

    output = capsys.readouterr().out

    assert output.strip() != ""

    # The banner is ASCII art, so use stable fragments rather
    # than depending on every whitespace character.
    assert "oooo" in output
    assert "888" in output


def test_ascii_banner_contains_multiple_lines(
    capsys,
):
    Initialization.ASCII_Mupt_reaction_LAMMPS()

    output = capsys.readouterr().out

    lines = [
        line
        for line in output.splitlines()
        if line.strip()
    ]

    assert len(lines) >= 7


# =============================================================================
# Version printing
# =============================================================================


def test_print_version_uses_autoreacter_version(
    monkeypatch,
    capsys,
):
    import AutoREACTER

    monkeypatch.setattr(
        AutoREACTER,
        "__version__",
        "9.9.9-test",
    )

    Initialization.print_version()

    output = capsys.readouterr().out

    assert (
        "AutoREACTER version: 9.9.9-test"
        in output
    )


def test_print_version_prints_exact_prefix(
    capsys,
):
    import AutoREACTER

    Initialization.print_version()

    output = capsys.readouterr().out.strip()

    assert output.startswith(
        "AutoREACTER version:"
    )

    assert AutoREACTER.__version__ in output