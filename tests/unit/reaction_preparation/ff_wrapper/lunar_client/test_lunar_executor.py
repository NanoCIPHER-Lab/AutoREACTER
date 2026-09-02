from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
import subprocess
import sys

import pytest

import AutoREACTER.reaction_preparation.ff_wrapper.lunar_client.lunar_executor as lunar_executor
from AutoREACTER.reaction_preparation.ff_wrapper.lunar_client.lunar_executor import (
    All2LMPResult,
    AtomTypingResult,
    LunarExecutor,
)


# =============================================================================
# Helpers
# =============================================================================


def make_monomer(
    *,
    name="mma",
    molecule_3Dmol_path=None,
    status=True,
):
    return SimpleNamespace(
        name=name,
        molecule_3Dmol_path=molecule_3Dmol_path,
        status=status,
    )


def make_inputs(
    monomers=None,
):
    return SimpleNamespace(
        monomers=list(
            monomers or []
        )
    )


def make_reaction(
    *,
    reaction_id=1,
    pre_path=None,
    post_path=None,
):
    return SimpleNamespace(
        reaction_id=reaction_id,
        reactant_combined_3Dmol_path=pre_path,
        product_combined_3Dmol_path=post_path,
    )


def make_executor(
    tmp_path,
):
    lunar_root = (
        tmp_path / "LUNAR"
    )

    lunar_root.mkdir(
        parents=True,
        exist_ok=True,
    )

    return LunarExecutor(
        lunar_location=lunar_root,
        cache_dir=tmp_path / "cache",
    )


def install_fake_atom_typing_subprocess(
    monkeypatch,
):
    """
    Fake atom_typing.py execution and create exactly the files that
    LunarExecutor expects.
    """
    calls = []

    def fake_run(
        command,
        check=False,
        **kwargs,
    ):
        calls.append(
            {
                "command": list(command),
                "check": check,
                **kwargs,
            }
        )

        topo = Path(
            command[
                command.index("-topo") + 1
            ]
        )

        output_dir = Path(
            command[
                command.index("-dir") + 1
            ]
        )

        output_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        # Monomer output naming is based on the input stem.
        stem = topo.stem

        (
            output_dir
            / f"{stem}_typed.data"
        ).write_text(
            "typed data",
            encoding="utf-8",
        )

        (
            output_dir
            / f"{stem}_typed.nta"
        ).write_text(
            "nta",
            encoding="utf-8",
        )

        return SimpleNamespace(
            returncode=0
        )

    monkeypatch.setattr(
        lunar_executor.subprocess,
        "run",
        fake_run,
    )

    return calls


# =============================================================================
# Fake FF-wrapper objects for bond_react_merge tests
# =============================================================================


@dataclass
class FakeDataFiles:
    data_file: Path
    lmp_molecule_file: Path


@dataclass
class FakeMoleculeFile:
    id: str
    molecule_files: FakeDataFiles


@dataclass
class FakeTemplateFile:
    reaction_id: int | None
    pre_reaction_file: FakeDataFiles
    post_reaction_file: FakeDataFiles


@dataclass
class FakeFFFiles:
    force_field_data: Path
    in_file: Path
    molecule_files: list
    template_files: list


def install_fake_ff_structures(
    monkeypatch,
):
    monkeypatch.setattr(
        lunar_executor,
        "DataFiles",
        FakeDataFiles,
    )

    monkeypatch.setattr(
        lunar_executor,
        "MoleculeFile",
        FakeMoleculeFile,
    )

    monkeypatch.setattr(
        lunar_executor,
        "TemplateFile",
        FakeTemplateFile,
    )

    monkeypatch.setattr(
        lunar_executor,
        "FFFiles",
        FakeFFFiles,
    )


# =============================================================================
# Dataclasses
# =============================================================================


def test_atom_typing_result_fields(
    tmp_path,
):
    result = AtomTypingResult(
        id="data1",
        molecule=True,
        typed_data_file=(
            tmp_path / "typed.data"
        ),
        nta_file=(
            tmp_path / "typed.nta"
        ),
    )

    assert result.id == "data1"
    assert result.molecule is True

    assert result.typed_data_file == (
        tmp_path / "typed.data"
    )

    assert result.nta_file == (
        tmp_path / "typed.nta"
    )


def test_atom_typing_result_uses_slots(
    tmp_path,
):
    result = AtomTypingResult(
        id="data1",
        molecule=True,
        typed_data_file=(
            tmp_path / "typed.data"
        ),
        nta_file=(
            tmp_path / "typed.nta"
        ),
    )

    with pytest.raises(
        AttributeError
    ):
        result.extra = 1


def test_all2lmp_result_fields():
    result = All2LMPResult(
        id="pre1",
        molecule=False,
        all2lmp_data_file=Path(
            "pre1_typed_IFF.data"
        ),
    )

    assert result.id == "pre1"
    assert result.molecule is False

    assert (
        result.all2lmp_data_file
        == Path(
            "pre1_typed_IFF.data"
        )
    )


def test_all2lmp_result_uses_slots():
    result = All2LMPResult(
        id="data1",
        molecule=True,
        all2lmp_data_file=Path(
            "data1_typed_IFF.data"
        ),
    )

    with pytest.raises(
        AttributeError
    ):
        result.extra = 1


# =============================================================================
# Constructor
# =============================================================================


def test_constructor_sets_lunar_script_paths(
    tmp_path,
):
    lunar_root = (
        tmp_path / "LUNAR"
    )

    executor = LunarExecutor(
        lunar_root,
        tmp_path / "cache",
    )

    assert executor.lunar_location == (
        lunar_root
    )

    assert executor.atom_typing_py == (
        lunar_root
        / "atom_typing.py"
    )

    assert executor.all2lmp_py == (
        lunar_root
        / "all2lmp.py"
    )

    assert executor.bond_react_merge_py == (
        lunar_root
        / "bond_react_merge.py"
    )


def test_constructor_creates_stage_caches(
    tmp_path,
):
    cache = (
        tmp_path / "cache"
    )

    executor = LunarExecutor(
        tmp_path / "LUNAR",
        cache,
    )

    assert executor.cache_atom_typing == (
        cache / "atom_typing"
    )

    assert executor.cache_all2lmp == (
        cache / "all2lmp"
    )

    assert (
        executor.cache_bond_react_merge
        == cache / "bond_react_merge"
    )

    assert (
        executor.cache_atom_typing.is_dir()
    )

    assert (
        executor.cache_all2lmp.is_dir()
    )

    assert (
        executor
        .cache_bond_react_merge
        .is_dir()
    )


def test_constructor_accepts_string_paths(
    tmp_path,
):
    executor = LunarExecutor(
        str(tmp_path / "LUNAR"),
        str(tmp_path / "cache"),
    )

    assert isinstance(
        executor.lunar_location,
        Path,
    )

    assert isinstance(
        executor.cache_atom_typing,
        Path,
    )


# =============================================================================
# run_atom_typing - monomers
# =============================================================================


def test_run_atom_typing_active_monomer(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    molecule_file = (
        tmp_path / "mma.mol"
    )

    molecule_file.write_text(
        "molecule",
        encoding="utf-8",
    )

    inputs = make_inputs(
        [
            make_monomer(
                name="mma",
                molecule_3Dmol_path=molecule_file,
            )
        ]
    )

    calls = (
        install_fake_atom_typing_subprocess(
            monkeypatch
        )
    )

    monkeypatch.setattr(
        lunar_executor.time,
        "sleep",
        lambda seconds: None,
    )

    results = executor.run_atom_typing(
        inputs,
        prepared_reactions=[],
        force_field="PCFF",
    )

    assert len(results) == 1

    result = results[0]

    assert result.id == "mma"
    assert result.molecule is True

    assert result.typed_data_file == (
        executor.cache_atom_typing
        / "mma_typed.data"
    )

    assert result.nta_file == (
        executor.cache_atom_typing
        / "mma_typed.nta"
    )

    assert len(calls) == 1


def test_run_atom_typing_command_for_monomer(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    molecule_file = (
        tmp_path / "styrene.mol"
    )

    molecule_file.write_text(
        "",
        encoding="utf-8",
    )

    inputs = make_inputs(
        [
            make_monomer(
                name="styrene",
                molecule_3Dmol_path=molecule_file,
            )
        ]
    )

    calls = (
        install_fake_atom_typing_subprocess(
            monkeypatch
        )
    )

    monkeypatch.setattr(
        lunar_executor.time,
        "sleep",
        lambda seconds: None,
    )

    executor.run_atom_typing(
        inputs,
        [],
        "PCFF",
    )

    assert calls[0][
        "command"
    ] == [
        sys.executable,
        str(
            executor.atom_typing_py
        ),
        "-topo",
        str(molecule_file),
        "-dir",
        str(
            executor.cache_atom_typing
        ),
        "-ff",
        "PCFF",
        "-del-method",
        "mass",
        "-del-crit",
        "0",
    ]

    assert (
        calls[0]["check"]
        is True
    )


def test_run_atom_typing_skips_inactive_monomer(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    calls = []

    monkeypatch.setattr(
        lunar_executor.subprocess,
        "run",
        lambda *args, **kwargs:
        calls.append(
            args
        ),
    )

    monkeypatch.setattr(
        lunar_executor.time,
        "sleep",
        lambda seconds: None,
    )

    inputs = make_inputs(
        [
            make_monomer(
                name="inactive",
                molecule_3Dmol_path=None,
                status=False,
            )
        ]
    )

    results = executor.run_atom_typing(
        inputs,
        [],
        "PCFF",
    )

    assert results == []
    assert calls == []


def test_run_atom_typing_missing_active_monomer_3d_path_raises(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    monkeypatch.setattr(
        lunar_executor.subprocess,
        "run",
        lambda *args, **kwargs:
        pytest.fail(
            "subprocess must not run"
        ),
    )

    inputs = make_inputs(
        [
            make_monomer(
                name="mma",
                molecule_3Dmol_path=None,
                status=True,
            )
        ]
    )

    with pytest.raises(
        ValueError,
        match="missing molecule_3Dmol_path",
    ):
        executor.run_atom_typing(
            inputs,
            [],
            "PCFF",
        )


def test_run_atom_typing_status_missing_defaults_active(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    molecule_file = (
        tmp_path / "mma.mol"
    )

    molecule_file.write_text(
        "",
        encoding="utf-8",
    )

    monomer = SimpleNamespace(
        name="mma",
        molecule_3Dmol_path=molecule_file,
    )

    inputs = make_inputs(
        [
            monomer
        ]
    )

    install_fake_atom_typing_subprocess(
        monkeypatch
    )

    monkeypatch.setattr(
        lunar_executor.time,
        "sleep",
        lambda seconds: None,
    )

    results = executor.run_atom_typing(
        inputs,
        [],
        "PCFF",
    )

    assert len(results) == 1
    assert results[0].id == "mma"


# =============================================================================
# run_atom_typing - reaction templates
# =============================================================================


def test_run_atom_typing_processes_pre_and_post_templates(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    pre = (
        tmp_path / "reaction_pre.mol"
    )

    post = (
        tmp_path / "reaction_post.mol"
    )

    pre.write_text(
        "",
        encoding="utf-8",
    )

    post.write_text(
        "",
        encoding="utf-8",
    )

    reaction = make_reaction(
        reaction_id=7,
        pre_path=pre,
        post_path=post,
    )

    install_fake_atom_typing_subprocess(
        monkeypatch
    )

    monkeypatch.setattr(
        lunar_executor.time,
        "sleep",
        lambda seconds: None,
    )

    results = executor.run_atom_typing(
        make_inputs(),
        [
            reaction
        ],
        "PCFF",
    )

    assert [
        result.id
        for result in results
    ] == [
        "pre7",
        "post7",
    ]

    assert all(
        result.molecule is False
        for result in results
    )


def test_run_atom_typing_template_outputs_use_separate_directories(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    pre = (
        tmp_path / "foo_pre.mol"
    )

    post = (
        tmp_path / "foo_post.mol"
    )

    pre.write_text("")
    post.write_text("")

    reaction = make_reaction(
        reaction_id=3,
        pre_path=pre,
        post_path=post,
    )

    install_fake_atom_typing_subprocess(
        monkeypatch
    )

    monkeypatch.setattr(
        lunar_executor.time,
        "sleep",
        lambda seconds: None,
    )

    results = executor.run_atom_typing(
        make_inputs(),
        [
            reaction
        ],
        "PCFF",
    )

    assert (
        results[0]
        .typed_data_file
        == executor.cache_atom_typing
        / "pre3"
        / "foo_pre_typed.data"
    )

    assert (
        results[1]
        .typed_data_file
        == executor.cache_atom_typing
        / "post3"
        / "foo_post_typed.data"
    )


def test_run_atom_typing_result_order_is_monomers_then_templates(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    monomer_file = (
        tmp_path / "mma.mol"
    )

    pre = (
        tmp_path / "pre.mol"
    )

    post = (
        tmp_path / "post.mol"
    )

    for path in (
        monomer_file,
        pre,
        post,
    ):
        path.write_text("")

    inputs = make_inputs(
        [
            make_monomer(
                name="mma",
                molecule_3Dmol_path=monomer_file,
            )
        ]
    )

    reaction = make_reaction(
        reaction_id=1,
        pre_path=pre,
        post_path=post,
    )

    install_fake_atom_typing_subprocess(
        monkeypatch
    )

    monkeypatch.setattr(
        lunar_executor.time,
        "sleep",
        lambda seconds: None,
    )

    results = executor.run_atom_typing(
        inputs,
        [
            reaction
        ],
        "PCFF",
    )

    assert [
        result.id
        for result in results
    ] == [
        "mma",
        "pre1",
        "post1",
    ]


def test_run_atom_typing_sleeps_after_every_executed_item(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    monomer_file = (
        tmp_path / "mma.mol"
    )

    pre = (
        tmp_path / "pre.mol"
    )

    post = (
        tmp_path / "post.mol"
    )

    for path in (
        monomer_file,
        pre,
        post,
    ):
        path.write_text("")

    inputs = make_inputs(
        [
            make_monomer(
                name="mma",
                molecule_3Dmol_path=monomer_file,
            )
        ]
    )

    reaction = make_reaction(
        reaction_id=1,
        pre_path=pre,
        post_path=post,
    )

    install_fake_atom_typing_subprocess(
        monkeypatch
    )

    sleeps = []

    monkeypatch.setattr(
        lunar_executor.time,
        "sleep",
        lambda seconds:
        sleeps.append(seconds),
    )

    executor.run_atom_typing(
        inputs,
        [
            reaction
        ],
        "PCFF",
    )

    assert sleeps == [
        0.1,
        0.1,
        0.1,
    ]


# =============================================================================
# run_atom_typing - output verification
# =============================================================================


def test_run_atom_typing_missing_expected_output_raises(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    molecule_file = (
        tmp_path / "mma.mol"
    )

    molecule_file.write_text("")

    monkeypatch.setattr(
        lunar_executor.subprocess,
        "run",
        lambda *args, **kwargs:
        SimpleNamespace(
            returncode=0
        ),
    )

    monkeypatch.setattr(
        lunar_executor.time,
        "sleep",
        lambda seconds: None,
    )

    inputs = make_inputs(
        [
            make_monomer(
                name="mma",
                molecule_3Dmol_path=molecule_file,
            )
        ]
    )

    with pytest.raises(
        FileNotFoundError,
        match="Expected LUNAR output not found for mma",
    ):
        executor.run_atom_typing(
            inputs,
            [],
            "PCFF",
        )


def test_run_atom_typing_prints_generated_message(
    tmp_path,
    monkeypatch,
    capsys,
):
    executor = make_executor(
        tmp_path
    )

    molecule_file = (
        tmp_path / "mma.mol"
    )

    molecule_file.write_text("")

    install_fake_atom_typing_subprocess(
        monkeypatch
    )

    monkeypatch.setattr(
        lunar_executor.time,
        "sleep",
        lambda seconds: None,
    )

    executor.run_atom_typing(
        make_inputs(
            [
                make_monomer(
                    name="mma",
                    molecule_3Dmol_path=molecule_file,
                )
            ]
        ),
        [],
        "PCFF",
    )

    assert (
        "[LUNAR atom_typing] Generated files for mma"
        in capsys.readouterr().out
    )


def test_run_atom_typing_subprocess_failure_propagates(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    molecule_file = (
        tmp_path / "mma.mol"
    )

    molecule_file.write_text("")

    def fail(*args, **kwargs):
        raise subprocess.CalledProcessError(
            1,
            "atom_typing",
        )

    monkeypatch.setattr(
        lunar_executor.subprocess,
        "run",
        fail,
    )

    with pytest.raises(
        subprocess.CalledProcessError
    ):
        executor.run_atom_typing(
            make_inputs(
                [
                    make_monomer(
                        name="mma",
                        molecule_3Dmol_path=molecule_file,
                    )
                ]
            ),
            [],
            "PCFF",
        )


# =============================================================================
# run_all2lmp
# =============================================================================


def install_fake_all2lmp_subprocess(
    executor,
    monkeypatch,
):
    calls = []

    def fake_run(
        command,
        check=False,
        **kwargs,
    ):
        calls.append(
            {
                "command": list(command),
                "check": check,
                **kwargs,
            }
        )

        topo = Path(
            command[
                command.index("-topo") + 1
            ]
        )

        # Executor output naming is derived from the identifier.
        # The typed input file stem is e.g. "pre1_typed".
        stem = topo.stem

        if stem.endswith(
            "_typed"
        ):
            result_id = stem[
                :-len("_typed")
            ]
        else:
            result_id = stem

        (
            executor.cache_all2lmp
            / f"{result_id}_typed_IFF.data"
        ).write_text(
            "IFF",
            encoding="utf-8",
        )

        return SimpleNamespace(
            returncode=0
        )

    monkeypatch.setattr(
        lunar_executor.subprocess,
        "run",
        fake_run,
    )

    return calls


def test_run_all2lmp_converts_single_result(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    typed = (
        tmp_path / "data1_typed.data"
    )

    nta = (
        tmp_path / "data1_typed.nta"
    )

    typed.write_text("")
    nta.write_text("")

    input_result = AtomTypingResult(
        id="data1",
        molecule=True,
        typed_data_file=typed,
        nta_file=nta,
    )

    install_fake_all2lmp_subprocess(
        executor,
        monkeypatch,
    )

    monkeypatch.setattr(
        lunar_executor.time,
        "sleep",
        lambda seconds: None,
    )

    results = executor.run_all2lmp(
        [
            input_result
        ],
        tmp_path / "pcff.frc",
    )

    assert len(results) == 1

    assert results[0] == (
        All2LMPResult(
            id="data1",
            molecule=True,
            all2lmp_data_file=Path(
                "data1_typed_IFF.data"
            ),
        )
    )


def test_run_all2lmp_exact_command(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    typed = (
        tmp_path / "pre2_typed.data"
    )

    nta = (
        tmp_path / "pre2_typed.nta"
    )

    typed.write_text("")
    nta.write_text("")

    frc = (
        tmp_path / "pcff.frc"
    )

    calls = (
        install_fake_all2lmp_subprocess(
            executor,
            monkeypatch,
        )
    )

    monkeypatch.setattr(
        lunar_executor.time,
        "sleep",
        lambda seconds: None,
    )

    executor.run_all2lmp(
        [
            AtomTypingResult(
                id="pre2",
                molecule=False,
                typed_data_file=typed,
                nta_file=nta,
            )
        ],
        frc,
    )

    assert calls[0][
        "command"
    ] == [
        sys.executable,
        str(
            executor.all2lmp_py
        ),
        "-topo",
        str(typed),
        "-nta",
        str(nta),
        "-frc",
        str(frc),
        "-asm",
        "T",
        "-dir",
        str(
            executor.cache_all2lmp
        ),
    ]

    assert (
        calls[0]["check"]
        is True
    )


def test_run_all2lmp_preserves_id_and_molecule_flag(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    typed = (
        tmp_path / "post9_typed.data"
    )

    nta = (
        tmp_path / "post9_typed.nta"
    )

    typed.write_text("")
    nta.write_text("")

    install_fake_all2lmp_subprocess(
        executor,
        monkeypatch,
    )

    monkeypatch.setattr(
        lunar_executor.time,
        "sleep",
        lambda seconds: None,
    )

    result = executor.run_all2lmp(
        [
            AtomTypingResult(
                id="post9",
                molecule=False,
                typed_data_file=typed,
                nta_file=nta,
            )
        ],
        tmp_path / "pcff.frc",
    )[0]

    assert result.id == "post9"
    assert result.molecule is False


def test_run_all2lmp_processes_results_in_input_order(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    entries = []

    for identifier, molecule in [
        ("data1", True),
        ("pre1", False),
        ("post1", False),
    ]:
        typed = (
            tmp_path
            / f"{identifier}_typed.data"
        )

        nta = (
            tmp_path
            / f"{identifier}_typed.nta"
        )

        typed.write_text("")
        nta.write_text("")

        entries.append(
            AtomTypingResult(
                id=identifier,
                molecule=molecule,
                typed_data_file=typed,
                nta_file=nta,
            )
        )

    install_fake_all2lmp_subprocess(
        executor,
        monkeypatch,
    )

    monkeypatch.setattr(
        lunar_executor.time,
        "sleep",
        lambda seconds: None,
    )

    results = executor.run_all2lmp(
        entries,
        tmp_path / "pcff.frc",
    )

    assert [
        result.id
        for result in results
    ] == [
        "data1",
        "pre1",
        "post1",
    ]


def test_run_all2lmp_sleeps_once_after_all_commands(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    typed = (
        tmp_path / "data1_typed.data"
    )

    nta = (
        tmp_path / "data1_typed.nta"
    )

    typed.write_text("")
    nta.write_text("")

    install_fake_all2lmp_subprocess(
        executor,
        monkeypatch,
    )

    sleeps = []

    monkeypatch.setattr(
        lunar_executor.time,
        "sleep",
        lambda seconds:
        sleeps.append(seconds),
    )

    executor.run_all2lmp(
        [
            AtomTypingResult(
                id="data1",
                molecule=True,
                typed_data_file=typed,
                nta_file=nta,
            )
        ],
        tmp_path / "pcff.frc",
    )

    assert sleeps == [
        0.1
    ]


def test_run_all2lmp_empty_input_still_sleeps_once(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    sleeps = []

    monkeypatch.setattr(
        lunar_executor.time,
        "sleep",
        lambda seconds:
        sleeps.append(seconds),
    )

    results = executor.run_all2lmp(
        [],
        tmp_path / "pcff.frc",
    )

    assert results == []

    assert sleeps == [
        0.1
    ]


def test_run_all2lmp_missing_output_raises(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    typed = (
        tmp_path / "data1_typed.data"
    )

    nta = (
        tmp_path / "data1_typed.nta"
    )

    typed.write_text("")
    nta.write_text("")

    monkeypatch.setattr(
        lunar_executor.subprocess,
        "run",
        lambda *args, **kwargs:
        SimpleNamespace(
            returncode=0
        ),
    )

    monkeypatch.setattr(
        lunar_executor.time,
        "sleep",
        lambda seconds: None,
    )

    with pytest.raises(
        FileNotFoundError,
        match="Expected output file not found",
    ):
        executor.run_all2lmp(
            [
                AtomTypingResult(
                    id="data1",
                    molecule=True,
                    typed_data_file=typed,
                    nta_file=nta,
                )
            ],
            tmp_path / "pcff.frc",
        )


def test_run_all2lmp_prints_generated_message(
    tmp_path,
    monkeypatch,
    capsys,
):
    executor = make_executor(
        tmp_path
    )

    typed = (
        tmp_path / "data1_typed.data"
    )

    nta = (
        tmp_path / "data1_typed.nta"
    )

    typed.write_text("")
    nta.write_text("")

    install_fake_all2lmp_subprocess(
        executor,
        monkeypatch,
    )

    monkeypatch.setattr(
        lunar_executor.time,
        "sleep",
        lambda seconds: None,
    )

    executor.run_all2lmp(
        [
            AtomTypingResult(
                id="data1",
                molecule=True,
                typed_data_file=typed,
                nta_file=nta,
            )
        ],
        tmp_path / "pcff.frc",
    )

    assert (
        "[LUNAR all2lmp] Generated file for data1"
        in capsys.readouterr().out
    )


# =============================================================================
# run_bond_react_merge - subprocess
# =============================================================================


def test_run_bond_react_merge_exact_command(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    merge_input = (
        executor.cache_bond_react_merge
        / "merge_input.txt"
    )

    merge_input.write_text("")

    calls = []

    def fake_run(
        command,
        **kwargs,
    ):
        calls.append(
            (
                list(command),
                kwargs,
            )
        )

        return SimpleNamespace(
            returncode=0
        )

    monkeypatch.setattr(
        lunar_executor.subprocess,
        "run",
        fake_run,
    )

    monkeypatch.setattr(
        lunar_executor,
        "move_merge_outputs",
        lambda src, dst: None,
    )

    install_fake_ff_structures(
        monkeypatch
    )

    monkeypatch.setattr(
        lunar_executor,
        "FFValidator",
        lambda files: None,
    )

    executor.run_bond_react_merge(
        merge_input,
        [],
    )

    command, kwargs = calls[0]

    assert command == [
        sys.executable,
        str(
            executor.bond_react_merge_py
        ),
        "-files",
        "infile:merge_input.txt",
        "-atomstyle",
        "full",
        "-tl",
        "T",
        "-wrd",
        "F",
        "-map",
        "F",
    ]

    assert kwargs[
        "cwd"
    ] == str(
        executor.cache_bond_react_merge
    )

    assert (
        kwargs["check"]
        is True
    )


def test_run_bond_react_merge_sets_qt_offscreen_environment(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    merge_input = (
        executor.cache_bond_react_merge
        / "merge_input.txt"
    )

    merge_input.write_text("")

    captured_env = []

    def fake_run(
        command,
        **kwargs,
    ):
        captured_env.append(
            kwargs["env"]
        )

        return SimpleNamespace(
            returncode=0
        )

    monkeypatch.setenv(
        "AUTOREACTER_TEST_ENV",
        "preserved",
    )

    monkeypatch.setattr(
        lunar_executor.subprocess,
        "run",
        fake_run,
    )

    monkeypatch.setattr(
        lunar_executor,
        "move_merge_outputs",
        lambda src, dst: None,
    )

    install_fake_ff_structures(
        monkeypatch
    )

    monkeypatch.setattr(
        lunar_executor,
        "FFValidator",
        lambda files: None,
    )

    executor.run_bond_react_merge(
        merge_input,
        [],
    )

    env = captured_env[0]

    assert (
        env["QT_QPA_PLATFORM"]
        == "offscreen"
    )

    assert (
        env["AUTOREACTER_TEST_ENV"]
        == "preserved"
    )


def test_run_bond_react_merge_subprocess_failure_propagates(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    merge_input = (
        tmp_path / "merge_input.txt"
    )

    merge_input.write_text("")

    def fail(*args, **kwargs):
        raise subprocess.CalledProcessError(
            1,
            "bond_react_merge",
        )

    monkeypatch.setattr(
        lunar_executor.subprocess,
        "run",
        fail,
    )

    monkeypatch.setattr(
        lunar_executor,
        "move_merge_outputs",
        lambda *args:
        pytest.fail(
            "outputs must not move after subprocess failure"
        ),
    )

    with pytest.raises(
        subprocess.CalledProcessError
    ):
        executor.run_bond_react_merge(
            merge_input,
            [],
        )


# =============================================================================
# run_bond_react_merge - moving outputs
# =============================================================================


def test_run_bond_react_merge_moves_all2lmp_outputs_after_execution(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    merge_input = (
        tmp_path / "merge_input.txt"
    )

    merge_input.write_text("")

    order = []

    monkeypatch.setattr(
        lunar_executor.subprocess,
        "run",
        lambda *args, **kwargs:
        (
            order.append(
                "subprocess"
            )
            or SimpleNamespace(
                returncode=0
            )
        ),
    )

    def fake_move(
        src,
        dst,
    ):
        order.append(
            "move"
        )

        assert src == (
            executor.cache_all2lmp
        )

        assert dst == (
            executor
            .cache_bond_react_merge
        )

    monkeypatch.setattr(
        lunar_executor,
        "move_merge_outputs",
        fake_move,
    )

    install_fake_ff_structures(
        monkeypatch
    )

    monkeypatch.setattr(
        lunar_executor,
        "FFValidator",
        lambda files:
        order.append(
            "validate"
        ),
    )

    executor.run_bond_react_merge(
        merge_input,
        [],
    )

    assert order == [
        "subprocess",
        "move",
        "validate",
    ]


# =============================================================================
# run_bond_react_merge - resulting FFFiles
# =============================================================================


def test_run_bond_react_merge_builds_molecule_files(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    merge_input = (
        tmp_path / "merge_input.txt"
    )

    merge_input.write_text("")

    monkeypatch.setattr(
        lunar_executor.subprocess,
        "run",
        lambda *args, **kwargs:
        SimpleNamespace(
            returncode=0
        ),
    )

    monkeypatch.setattr(
        lunar_executor,
        "move_merge_outputs",
        lambda src, dst: None,
    )

    install_fake_ff_structures(
        monkeypatch
    )

    monkeypatch.setattr(
        lunar_executor,
        "FFValidator",
        lambda files: None,
    )

    result = (
        executor.run_bond_react_merge(
            merge_input,
            [
                All2LMPResult(
                    id="mma",
                    molecule=True,
                    all2lmp_data_file=Path(
                        "mma_typed_IFF.data"
                    ),
                )
            ],
        )
    )

    assert len(
        result.molecule_files
    ) == 1

    molecule_file = (
        result.molecule_files[0]
    )

    assert molecule_file.id == "mma"

    assert (
        molecule_file
        .molecule_files
        .data_file
        == executor.cache_bond_react_merge
        / "mma_typed_IFF_merged.data"
    )

    assert (
        molecule_file
        .molecule_files
        .lmp_molecule_file
        == executor.cache_bond_react_merge
        / "mma_typed_IFF_merged.lmpmol"
    )


def test_run_bond_react_merge_builds_template_from_pre_result(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    merge_input = (
        tmp_path / "merge_input.txt"
    )

    merge_input.write_text("")

    monkeypatch.setattr(
        lunar_executor.subprocess,
        "run",
        lambda *args, **kwargs:
        SimpleNamespace(
            returncode=0
        ),
    )

    monkeypatch.setattr(
        lunar_executor,
        "move_merge_outputs",
        lambda src, dst: None,
    )

    install_fake_ff_structures(
        monkeypatch
    )

    monkeypatch.setattr(
        lunar_executor,
        "FFValidator",
        lambda files: None,
    )

    results = [
        All2LMPResult(
            id="pre12",
            molecule=False,
            all2lmp_data_file=Path(
                "pre12_typed_IFF.data"
            ),
        ),
        All2LMPResult(
            id="post12",
            molecule=False,
            all2lmp_data_file=Path(
                "post12_typed_IFF.data"
            ),
        ),
    ]

    result = (
        executor.run_bond_react_merge(
            merge_input,
            results,
        )
    )

    assert len(
        result.template_files
    ) == 1

    template = (
        result.template_files[0]
    )

    assert (
        template.reaction_id
        == 12
    )

    assert (
        template
        .pre_reaction_file
        .data_file
        == executor.cache_bond_react_merge
        / "pre12_typed_IFF_merged.data"
    )

    assert (
        template
        .post_reaction_file
        .data_file
        == executor.cache_bond_react_merge
        / "post12_typed_IFF_merged.data"
    )


def test_run_bond_react_merge_post_result_does_not_create_second_template(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    merge_input = (
        tmp_path / "merge_input.txt"
    )

    merge_input.write_text("")

    monkeypatch.setattr(
        lunar_executor.subprocess,
        "run",
        lambda *args, **kwargs:
        SimpleNamespace(
            returncode=0
        ),
    )

    monkeypatch.setattr(
        lunar_executor,
        "move_merge_outputs",
        lambda src, dst: None,
    )

    install_fake_ff_structures(
        monkeypatch
    )

    monkeypatch.setattr(
        lunar_executor,
        "FFValidator",
        lambda files: None,
    )

    result = executor.run_bond_react_merge(
        merge_input,
        [
            All2LMPResult(
                id="pre1",
                molecule=False,
                all2lmp_data_file=Path(
                    "pre1_typed_IFF.data"
                ),
            ),
            All2LMPResult(
                id="post1",
                molecule=False,
                all2lmp_data_file=Path(
                    "post1_typed_IFF.data"
                ),
            ),
        ],
    )

    assert len(
        result.template_files
    ) == 1


def test_run_bond_react_merge_sets_force_field_and_input_files(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    merge_input = (
        tmp_path / "merge_input.txt"
    )

    merge_input.write_text("")

    monkeypatch.setattr(
        lunar_executor.subprocess,
        "run",
        lambda *args, **kwargs:
        SimpleNamespace(
            returncode=0
        ),
    )

    monkeypatch.setattr(
        lunar_executor,
        "move_merge_outputs",
        lambda src, dst: None,
    )

    install_fake_ff_structures(
        monkeypatch
    )

    monkeypatch.setattr(
        lunar_executor,
        "FFValidator",
        lambda files: None,
    )

    result = (
        executor.run_bond_react_merge(
            merge_input,
            [],
        )
    )

    assert (
        result.force_field_data
        == executor.cache_bond_react_merge
        / "force_field.data"
    )

    # Characterizes current implementation:
    # in.create_atoms.script remains referenced in all2lmp cache.
    assert (
        result.in_file
        == executor.cache_all2lmp
        / "in.create_atoms.script"
    )


def test_run_bond_react_merge_validates_final_files(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    merge_input = (
        tmp_path / "merge_input.txt"
    )

    merge_input.write_text("")

    monkeypatch.setattr(
        lunar_executor.subprocess,
        "run",
        lambda *args, **kwargs:
        SimpleNamespace(
            returncode=0
        ),
    )

    monkeypatch.setattr(
        lunar_executor,
        "move_merge_outputs",
        lambda src, dst: None,
    )

    install_fake_ff_structures(
        monkeypatch
    )

    validated = []

    monkeypatch.setattr(
        lunar_executor,
        "FFValidator",
        lambda files:
        validated.append(files),
    )

    result = (
        executor.run_bond_react_merge(
            merge_input,
            [],
        )
    )

    assert validated == [
        result
    ]


def test_run_bond_react_merge_returns_fffiles(
    tmp_path,
    monkeypatch,
):
    executor = make_executor(
        tmp_path
    )

    merge_input = (
        tmp_path / "merge_input.txt"
    )

    merge_input.write_text("")

    monkeypatch.setattr(
        lunar_executor.subprocess,
        "run",
        lambda *args, **kwargs:
        SimpleNamespace(
            returncode=0
        ),
    )

    monkeypatch.setattr(
        lunar_executor,
        "move_merge_outputs",
        lambda src, dst: None,
    )

    install_fake_ff_structures(
        monkeypatch
    )

    monkeypatch.setattr(
        lunar_executor,
        "FFValidator",
        lambda files: None,
    )

    result = (
        executor.run_bond_react_merge(
            merge_input,
            [],
        )
    )

    assert isinstance(
        result,
        FakeFFFiles,
    )