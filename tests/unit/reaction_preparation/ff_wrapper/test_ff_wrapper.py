import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace

import pytest

from AutoREACTER.reaction_preparation.ff_wrapper.ff_wrapper import (
    DataFiles,
    FFFiles,
    FFWrapper,
    MoleculeFile,
    TemplateFile,
)


# =============================================================================
# Helpers
# =============================================================================


def make_session(
    *,
    force_field="PCFF",
    reaction_metadata=None,
):
    inputs = SimpleNamespace(
        force_field=force_field,
    )

    return SimpleNamespace(
        inputs=inputs,
        reaction_metadata=(
            []
            if reaction_metadata is None
            else reaction_metadata
        ),
        ff_files=None,
    )


def install_fake_lunar_wrapper(
    monkeypatch,
    *,
    final_files="lunar-final-files",
    events=None,
):
    """
    Install a fake module for the lazy import performed inside
    FFWrapper.generate_force_field_files().
    """
    module_name = (
        "AutoREACTER.reaction_preparation."
        "ff_wrapper.lunar_client.lunar_api_wrapper"
    )

    fake_module = ModuleType(
        module_name
    )

    class FakeLunarAPIWrapper:
        def __init__(
            self,
            ARX,
        ):
            self.ARX = ARX

            if events is not None:
                events.append(
                    (
                        "lunar-init",
                        ARX,
                    )
                )

        def lunar_workflow(
            self,
            updated_inputs,
            prepared_reactions,
        ):
            if events is not None:
                events.append(
                    (
                        "lunar-workflow",
                        updated_inputs,
                        prepared_reactions,
                    )
                )

            return final_files

    fake_module.LunarAPIWrapper = (
        FakeLunarAPIWrapper
    )

    monkeypatch.setitem(
        sys.modules,
        module_name,
        fake_module,
    )

    return FakeLunarAPIWrapper


def install_fake_foyer_wrapper(
    monkeypatch,
    *,
    final_files="foyer-final-files",
    events=None,
):
    """
    Install a fake module for the current lazy Foyer routing branch.

    This does NOT claim Foyer support is functional. It only
    characterizes the current FFWrapper routing behavior.
    """
    module_name = (
        "AutoREACTER.reaction_preparation."
        "ff_wrapper.foyer_client.foyer_api_wrapper"
    )

    fake_module = ModuleType(
        module_name
    )

    class FakeFoyerAPIWrapper:
        def __init__(
            self,
            ARX,
            prepared_reactions_with_3d_mols,
        ):
            self.ARX = ARX
            self.prepared_reactions = (
                prepared_reactions_with_3d_mols
            )

            self.final_foyer_files = (
                final_files
            )

            if events is not None:
                events.append(
                    (
                        "foyer-init",
                        ARX,
                        prepared_reactions_with_3d_mols,
                    )
                )

    fake_module.FoyerAPIWrapper = (
        FakeFoyerAPIWrapper
    )

    monkeypatch.setitem(
        sys.modules,
        module_name,
        fake_module,
    )

    return FakeFoyerAPIWrapper


# =============================================================================
# DataFiles
# =============================================================================


def test_data_files_stores_paths(
    tmp_path,
):
    data_file = (
        tmp_path / "system.data"
    )

    molecule_file = (
        tmp_path / "system.lmpmol"
    )

    result = DataFiles(
        data_file=data_file,
        lmp_molecule_file=molecule_file,
    )

    assert result.data_file == (
        data_file
    )

    assert (
        result.lmp_molecule_file
        == molecule_file
    )


def test_data_files_uses_slots(
    tmp_path,
):
    result = DataFiles(
        data_file=(
            tmp_path / "a.data"
        ),
        lmp_molecule_file=(
            tmp_path / "a.lmpmol"
        ),
    )

    with pytest.raises(
        AttributeError
    ):
        result.extra = 1


# =============================================================================
# MoleculeFile
# =============================================================================


def test_molecule_file_stores_data_files(
    tmp_path,
):
    files = DataFiles(
        data_file=(
            tmp_path / "mma.data"
        ),
        lmp_molecule_file=(
            tmp_path / "mma.lmpmol"
        ),
    )

    result = MoleculeFile(
        id="mma",
        molecule_files=files,
    )

    assert result.id == "mma"

    assert (
        result.molecule_files
        is files
    )


def test_molecule_file_accepts_none():
    result = MoleculeFile(
        id="mma",
        molecule_files=None,
    )

    assert (
        result.molecule_files
        is None
    )


def test_molecule_file_uses_slots():
    result = MoleculeFile(
        id="mma",
        molecule_files=None,
    )

    with pytest.raises(
        AttributeError
    ):
        result.extra = 1


# =============================================================================
# TemplateFile
# =============================================================================


def test_template_file_stores_pre_and_post(
    tmp_path,
):
    pre = DataFiles(
        data_file=(
            tmp_path / "pre.data"
        ),
        lmp_molecule_file=(
            tmp_path / "pre.lmpmol"
        ),
    )

    post = DataFiles(
        data_file=(
            tmp_path / "post.data"
        ),
        lmp_molecule_file=(
            tmp_path / "post.lmpmol"
        ),
    )

    result = TemplateFile(
        reaction_id=7,
        pre_reaction_file=pre,
        post_reaction_file=post,
    )

    assert result.reaction_id == 7

    assert (
        result.pre_reaction_file
        is pre
    )

    assert (
        result.post_reaction_file
        is post
    )


def test_template_file_accepts_optional_values():
    result = TemplateFile(
        reaction_id=None,
        pre_reaction_file=None,
        post_reaction_file=None,
    )

    assert result.reaction_id is None
    assert result.pre_reaction_file is None
    assert result.post_reaction_file is None


def test_template_file_uses_slots():
    result = TemplateFile(
        reaction_id=None,
        pre_reaction_file=None,
        post_reaction_file=None,
    )

    with pytest.raises(
        AttributeError
    ):
        result.extra = 1


# =============================================================================
# FFFiles
# =============================================================================


def test_ff_files_stores_complete_output(
    tmp_path,
):
    molecule = MoleculeFile(
        id="mma",
        molecule_files=None,
    )

    template = TemplateFile(
        reaction_id=1,
        pre_reaction_file=None,
        post_reaction_file=None,
    )

    result = FFFiles(
        molecule_files=[
            molecule
        ],
        template_files=[
            template
        ],
        force_field_data=(
            tmp_path
            / "force_field.data"
        ),
        in_file=(
            tmp_path
            / "in.create_atoms.script"
        ),
    )

    assert result.molecule_files == [
        molecule
    ]

    assert result.template_files == [
        template
    ]

    assert result.force_field_data == (
        tmp_path
        / "force_field.data"
    )

    assert result.in_file == (
        tmp_path
        / "in.create_atoms.script"
    )


def test_ff_files_in_file_defaults_none(
    tmp_path,
):
    result = FFFiles(
        molecule_files=[],
        template_files=[],
        force_field_data=(
            tmp_path / "ff.data"
        ),
    )

    assert result.in_file is None


def test_ff_files_uses_slots(
    tmp_path,
):
    result = FFFiles(
        molecule_files=[],
        template_files=[],
        force_field_data=(
            tmp_path / "ff.data"
        ),
    )

    with pytest.raises(
        AttributeError
    ):
        result.extra = 1


# =============================================================================
# FFWrapper constructor
# =============================================================================


def test_ff_wrapper_constructor_stores_session():
    session = make_session()

    wrapper = FFWrapper(
        session
    )

    assert wrapper.session is session


def test_ff_wrapper_constructor_stores_inputs():
    session = make_session(
        force_field="CVFF"
    )

    wrapper = FFWrapper(
        session
    )

    assert (
        wrapper.inputs
        is session.inputs
    )


# =============================================================================
# LUNAR routing
# =============================================================================


@pytest.mark.parametrize(
    "force_field",
    [
        "PCFF",
        "PCFF-IFF",
        "compass",
        "CVFF",
        "CVFF-IFF",
        "DREIDING",
        "Clay-FF",
    ],
)
def test_supported_lunar_names_route_to_lunar(
    monkeypatch,
    capsys,
    force_field,
):
    session = make_session(
        force_field=force_field,
    )

    events = []

    install_fake_lunar_wrapper(
        monkeypatch,
        final_files="LUNAR_RESULT",
        events=events,
    )

    wrapper = FFWrapper(
        session
    )

    result = (
        wrapper
        .generate_force_field_files(
            session
        )
    )

    assert result is None

    assert (
        session.ff_files
        == "LUNAR_RESULT"
    )

    assert events[0] == (
        "lunar-init",
        session,
    )

    assert events[1] == (
        "lunar-workflow",
        session.inputs,
        session.reaction_metadata,
    )

    assert (
        f"Routing to LUNAR for force field: {force_field}"
        in capsys.readouterr().out
    )


def test_lunar_wrapper_receives_original_wrapper_session(
    monkeypatch,
):
    original_session = (
        make_session(
            force_field="PCFF"
        )
    )

    passed_session = (
        make_session(
            force_field="CVFF"
        )
    )

    events = []

    install_fake_lunar_wrapper(
        monkeypatch,
        events=events,
    )

    wrapper = FFWrapper(
        original_session
    )

    wrapper.generate_force_field_files(
        passed_session
    )

    # Characterizes current implementation:
    #
    # routing decisions and workflow arguments come from
    # the method's session, but LunarAPIWrapper is
    # constructed using self.session from FFWrapper.__init__.
    assert events[0] == (
        "lunar-init",
        original_session,
    )

    assert events[1] == (
        "lunar-workflow",
        passed_session.inputs,
        passed_session.reaction_metadata,
    )


def test_lunar_result_is_written_to_passed_session(
    monkeypatch,
):
    original_session = (
        make_session(
            force_field="PCFF"
        )
    )

    passed_session = (
        make_session(
            force_field="CVFF"
        )
    )

    original_session.ff_files = (
        "ORIGINAL"
    )

    install_fake_lunar_wrapper(
        monkeypatch,
        final_files="NEW_FILES",
    )

    wrapper = FFWrapper(
        original_session
    )

    wrapper.generate_force_field_files(
        passed_session
    )

    assert (
        passed_session.ff_files
        == "NEW_FILES"
    )

    # Current behavior: result is assigned to the
    # session passed to generate_force_field_files().
    assert (
        original_session.ff_files
        == "ORIGINAL"
    )


# =============================================================================
# Force-field fallback
# =============================================================================


def test_none_force_field_falls_back_to_pcff(
    monkeypatch,
    capsys,
):
    session = make_session(
        force_field=None,
    )

    install_fake_lunar_wrapper(
        monkeypatch,
        final_files="PCFF_RESULT",
    )

    wrapper = FFWrapper(
        session
    )

    wrapper.generate_force_field_files(
        session
    )

    assert (
        session.ff_files
        == "PCFF_RESULT"
    )

    assert (
        "Routing to LUNAR for force field: PCFF"
        in capsys.readouterr().out
    )


def test_empty_force_field_falls_back_to_pcff(
    monkeypatch,
):
    session = make_session(
        force_field="",
    )

    install_fake_lunar_wrapper(
        monkeypatch,
        final_files="PCFF_RESULT",
    )

    wrapper = FFWrapper(
        session
    )

    wrapper.generate_force_field_files(
        session
    )

    assert (
        session.ff_files
        == "PCFF_RESULT"
    )


# =============================================================================
# session.inputs fallback
# =============================================================================


def test_missing_passed_session_inputs_falls_back_to_wrapper_inputs(
    monkeypatch,
):
    original_session = (
        make_session(
            force_field="PCFF"
        )
    )

    passed_session = SimpleNamespace(
        inputs=None,
        reaction_metadata=[],
        ff_files=None,
    )

    events = []

    install_fake_lunar_wrapper(
        monkeypatch,
        final_files="RESULT",
        events=events,
    )

    wrapper = FFWrapper(
        original_session
    )

    wrapper.generate_force_field_files(
        passed_session
    )

    assert (
        passed_session.ff_files
        == "RESULT"
    )

    assert events[1] == (
        "lunar-workflow",
        original_session.inputs,
        passed_session.reaction_metadata,
    )


def test_both_input_sources_none_raises():
    original_session = SimpleNamespace(
        inputs=None,
    )

    passed_session = SimpleNamespace(
        inputs=None,
        reaction_metadata=[],
        ff_files=None,
    )

    wrapper = FFWrapper(
        original_session
    )

    with pytest.raises(
        ValueError,
        match=(
            "No inputs provided to FFWrapper "
            "and session inputs are None"
        ),
    ):
        wrapper.generate_force_field_files(
            passed_session
        )


# =============================================================================
# Reaction metadata forwarding
# =============================================================================


def test_reaction_metadata_is_forwarded_to_lunar(
    monkeypatch,
):
    reactions = [
        object(),
        object(),
        object(),
    ]

    session = make_session(
        force_field="PCFF",
        reaction_metadata=reactions,
    )

    events = []

    install_fake_lunar_wrapper(
        monkeypatch,
        events=events,
    )

    wrapper = FFWrapper(
        session
    )

    wrapper.generate_force_field_files(
        session
    )

    assert events[1][2] is reactions


# =============================================================================
# Foyer routing - current behavior only
# =============================================================================


@pytest.mark.parametrize(
    "force_field",
    [
        "OPLSAA",
        "GAFF",
    ],
)
def test_current_foyer_names_route_to_foyer(
    monkeypatch,
    capsys,
    force_field,
):
    """
    Characterization only.

    AutoREACTER currently marks Foyer force-field lookup as
    not implemented, but FFWrapper still contains a routing
    branch for OPLSAA and GAFF.

    This test does NOT claim Foyer works.
    """
    reactions = [
        object(),
    ]

    session = make_session(
        force_field=force_field,
        reaction_metadata=reactions,
    )

    events = []

    install_fake_foyer_wrapper(
        monkeypatch,
        final_files="FOYER_RESULT",
        events=events,
    )

    wrapper = FFWrapper(
        session
    )

    result = (
        wrapper
        .generate_force_field_files(
            session
        )
    )

    assert result is None

    assert (
        session.ff_files
        == "FOYER_RESULT"
    )

    assert events == [
        (
            "foyer-init",
            session,
            reactions,
        )
    ]

    assert (
        f"Routing to Foyer for force field: {force_field}"
        in capsys.readouterr().out
    )


# =============================================================================
# Unsupported force fields
# =============================================================================


@pytest.mark.parametrize(
    "force_field",
    [
        "AMBER",
        "CHARMM",
        "UFF",
        "MMFF94",
        "UNKNOWN",
        "pcff",
        "cvff",
        "COMPASS",
    ],
)
def test_unrecognized_force_field_raises(
    force_field,
):
    session = make_session(
        force_field=force_field,
    )

    wrapper = FFWrapper(
        session
    )

    with pytest.raises(
        ValueError,
        match=(
            "Unsupported or unrecognized "
            "force field requested"
        ),
    ) as exc_info:
        wrapper.generate_force_field_files(
            session
        )

    assert force_field in str(
        exc_info.value
    )


def test_unsupported_force_field_does_not_change_ff_files():
    session = make_session(
        force_field="UNKNOWN",
    )

    session.ff_files = (
        "existing-files"
    )

    wrapper = FFWrapper(
        session
    )

    with pytest.raises(
        ValueError
    ):
        wrapper.generate_force_field_files(
            session
        )

    assert (
        session.ff_files
        == "existing-files"
    )


# =============================================================================
# Backend failure propagation
# =============================================================================


def test_lunar_constructor_failure_propagates(
    monkeypatch,
):
    session = make_session(
        force_field="PCFF"
    )

    module_name = (
        "AutoREACTER.reaction_preparation."
        "ff_wrapper.lunar_client.lunar_api_wrapper"
    )

    fake_module = ModuleType(
        module_name
    )

    class FailingLunarAPIWrapper:
        def __init__(
            self,
            ARX,
        ):
            raise RuntimeError(
                "LUNAR init failed"
            )

    fake_module.LunarAPIWrapper = (
        FailingLunarAPIWrapper
    )

    monkeypatch.setitem(
        sys.modules,
        module_name,
        fake_module,
    )

    wrapper = FFWrapper(
        session
    )

    with pytest.raises(
        RuntimeError,
        match="LUNAR init failed",
    ):
        wrapper.generate_force_field_files(
            session
        )


def test_lunar_workflow_failure_propagates(
    monkeypatch,
):
    session = make_session(
        force_field="PCFF"
    )

    module_name = (
        "AutoREACTER.reaction_preparation."
        "ff_wrapper.lunar_client.lunar_api_wrapper"
    )

    fake_module = ModuleType(
        module_name
    )

    class FailingLunarAPIWrapper:
        def __init__(
            self,
            ARX,
        ):
            pass

        def lunar_workflow(
            self,
            updated_inputs,
            prepared_reactions,
        ):
            raise RuntimeError(
                "LUNAR workflow failed"
            )

    fake_module.LunarAPIWrapper = (
        FailingLunarAPIWrapper
    )

    monkeypatch.setitem(
        sys.modules,
        module_name,
        fake_module,
    )

    wrapper = FFWrapper(
        session
    )

    with pytest.raises(
        RuntimeError,
        match="LUNAR workflow failed",
    ):
        wrapper.generate_force_field_files(
            session
        )


def test_backend_failure_does_not_replace_existing_ff_files(
    monkeypatch,
):
    session = make_session(
        force_field="PCFF"
    )

    session.ff_files = (
        "old-files"
    )

    module_name = (
        "AutoREACTER.reaction_preparation."
        "ff_wrapper.lunar_client.lunar_api_wrapper"
    )

    fake_module = ModuleType(
        module_name
    )

    class FailingLunarAPIWrapper:
        def __init__(
            self,
            ARX,
        ):
            pass

        def lunar_workflow(
            self,
            updated_inputs,
            prepared_reactions,
        ):
            raise RuntimeError(
                "failure"
            )

    fake_module.LunarAPIWrapper = (
        FailingLunarAPIWrapper
    )

    monkeypatch.setitem(
        sys.modules,
        module_name,
        fake_module,
    )

    wrapper = FFWrapper(
        session
    )

    with pytest.raises(
        RuntimeError
    ):
        wrapper.generate_force_field_files(
            session
        )

    assert (
        session.ff_files
        == "old-files"
    )


# =============================================================================
# Return contract
# =============================================================================


def test_generate_force_field_files_returns_none_on_lunar_success(
    monkeypatch,
):
    session = make_session(
        force_field="PCFF"
    )

    install_fake_lunar_wrapper(
        monkeypatch
    )

    wrapper = FFWrapper(
        session
    )

    result = (
        wrapper
        .generate_force_field_files(
            session
        )
    )

    assert result is None


def test_generate_force_field_files_returns_none_on_current_foyer_branch(
    monkeypatch,
):
    session = make_session(
        force_field="GAFF"
    )

    install_fake_foyer_wrapper(
        monkeypatch
    )

    wrapper = FFWrapper(
        session
    )

    result = (
        wrapper
        .generate_force_field_files(
            session
        )
    )

    assert result is None