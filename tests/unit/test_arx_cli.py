from contextlib import contextmanager
import importlib
import os
from pathlib import Path
from types import SimpleNamespace

import pytest
from PIL import Image as PILImage


cli_module = importlib.import_module(
    "AutoREACTER.arx_cli"
)

ARXCLI = cli_module.ARXCLI
ErrorHandler = cli_module.ErrorHandler
NoReactionGenerated = (
    cli_module.NoReactionGenerated
)


# =============================================================================
# Helpers
# =============================================================================


def make_cli(tmp_path):
    """
    Create an ARXCLI instance without running its real constructor.

    This is useful for unit-testing individual controller methods without
    triggering read_input(), RDKit detection, force-field generation, or
    any other heavy pipeline component.
    """
    cli = ARXCLI.__new__(ARXCLI)

    output_dir = tmp_path / "output"
    images_dir = output_dir / "images"

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )
    images_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    cli.input = tmp_path / "input.json"

    cli.session = SimpleNamespace(
        output_dir=output_dir,
        images_dir=images_dir,
        reaction_instances=[],
        non_reactants=[],
    )

    cli.img_dir = images_dir

    cli.error_handler = (
        ErrorHandler().waterfall_order()
    )

    cli._fg_detected = False
    cli._reactions_detected = False
    cli._reactions_selected = False
    cli._non_reactants_detected = False
    cli._non_reactants_selected = False

    return cli


# =============================================================================
# ErrorHandler
# =============================================================================


def test_error_handler_waterfall_order():
    handler = ErrorHandler()

    result = handler.waterfall_order()

    assert result == {
        "select_reactions": False,
        "select_non_reactants": False,
        "process": False,
    }


def test_error_handler_returns_independent_dicts():
    handler = ErrorHandler()

    first = handler.waterfall_order()
    second = handler.waterfall_order()

    first["select_reactions"] = True

    assert (
        second["select_reactions"]
        is False
    )


# =============================================================================
# ARXCLI constructor
# =============================================================================


def test_arxcli_constructor_initializes_expected_state(
    tmp_path,
    monkeypatch,
):
    input_path = tmp_path / "input.json"

    session = SimpleNamespace(
        output_dir=tmp_path / "output",
        images_dir=tmp_path / "output" / "images",
        reaction_instances=[],
        non_reactants=[],
    )

    read_calls = []
    saved_input_paths = []
    saved_images = []
    ensure_calls = []

    fake_monomer_image = object()

    def fake_read_input(path):
        read_calls.append(path)
        return session

    class FakeInputParser:
        def initial_molecules_image_grid(
            self,
            passed_session,
        ):
            assert passed_session is session
            return fake_monomer_image

    def fake_save_input(
        self,
        path,
    ):
        saved_input_paths.append(path)

    def fake_save_image(
        self,
        image,
        path,
        is_non_reactant=False,
    ):
        saved_images.append(
            (
                image,
                Path(path),
                is_non_reactant,
            )
        )

    def fake_ensure_reactions(self):
        ensure_calls.append(True)

        # Verify flags are initialized before bootstrap.
        assert self._fg_detected is False
        assert self._reactions_detected is False
        assert self._reactions_selected is False

        assert (
            self._non_reactants_detected
            is False
        )

        assert (
            self._non_reactants_selected
            is False
        )

    monkeypatch.setattr(
        cli_module,
        "read_input",
        fake_read_input,
    )

    monkeypatch.setattr(
        cli_module,
        "InputParser",
        FakeInputParser,
    )

    monkeypatch.setattr(
        ARXCLI,
        "_save_input_json",
        fake_save_input,
    )

    monkeypatch.setattr(
        ARXCLI,
        "_save_rdkit_img",
        fake_save_image,
    )

    monkeypatch.setattr(
        ARXCLI,
        "_ensure_reactions_detected",
        fake_ensure_reactions,
    )

    cli = ARXCLI(input_path)

    assert cli.input == input_path
    assert cli.session is session

    assert (
        cli.img_dir
        == session.images_dir
    )

    assert read_calls == [
        input_path.resolve()
    ]

    assert saved_input_paths == [
        input_path.resolve()
    ]

    assert saved_images == [
        (
            fake_monomer_image,
            session.images_dir
            / "monomers.png",
            False,
        )
    ]

    assert ensure_calls == [
        True
    ]

    assert cli.error_handler == {
        "select_reactions": False,
        "select_non_reactants": False,
        "process": False,
    }


def test_arxcli_instances_do_not_share_waterfall_state(
    tmp_path,
    monkeypatch,
):
    input_1 = tmp_path / "one.json"
    input_2 = tmp_path / "two.json"

    session_1 = SimpleNamespace(
        output_dir=tmp_path / "out1",
        images_dir=tmp_path / "out1" / "images",
    )

    session_2 = SimpleNamespace(
        output_dir=tmp_path / "out2",
        images_dir=tmp_path / "out2" / "images",
    )

    sessions = iter(
        [
            session_1,
            session_2,
        ]
    )

    monkeypatch.setattr(
        cli_module,
        "read_input",
        lambda path: next(sessions),
    )

    class FakeInputParser:
        def initial_molecules_image_grid(
            self,
            session,
        ):
            return object()

    monkeypatch.setattr(
        cli_module,
        "InputParser",
        FakeInputParser,
    )

    monkeypatch.setattr(
        ARXCLI,
        "_save_input_json",
        lambda self, path: None,
    )

    monkeypatch.setattr(
        ARXCLI,
        "_save_rdkit_img",
        lambda self, image, path,
        is_non_reactant=False: None,
    )

    monkeypatch.setattr(
        ARXCLI,
        "_ensure_reactions_detected",
        lambda self: None,
    )

    cli_1 = ARXCLI(input_1)
    cli_2 = ARXCLI(input_2)

    cli_1.error_handler[
        "select_reactions"
    ] = True

    assert (
        cli_2.error_handler[
            "select_reactions"
        ]
        is False
    )


# =============================================================================
# show_molecules
# =============================================================================


def test_show_molecules_delegates_to_input_parser(
    tmp_path,
    monkeypatch,
):
    cli = make_cli(tmp_path)

    expected = object()
    calls = []

    class FakeInputParser:
        def initial_molecules_image_grid(
            self,
            session,
        ):
            calls.append(session)
            return expected

    monkeypatch.setattr(
        cli_module,
        "InputParser",
        FakeInputParser,
    )

    result = cli.show_molecules()

    assert result is expected
    assert calls == [
        cli.session
    ]


# =============================================================================
# show_functional_groups
# =============================================================================


def test_show_functional_groups_triggers_detection_when_needed(
    tmp_path,
    monkeypatch,
):
    cli = make_cli(tmp_path)

    expected = object()
    ensure_calls = []

    def fake_ensure():
        ensure_calls.append(True)
        cli._fg_detected = True

    monkeypatch.setattr(
        cli,
        "_ensure_fg_detected",
        fake_ensure,
    )

    class FakeDetector:
        def functional_group_highlighted_molecules_image_grid(
            self,
            session,
        ):
            assert session is cli.session
            return expected

    monkeypatch.setattr(
        cli_module,
        "FunctionalGroupsDetector",
        FakeDetector,
    )

    result = cli.show_functional_groups()

    assert result is expected
    assert ensure_calls == [
        True
    ]


def test_show_functional_groups_does_not_redetect_when_already_detected(
    tmp_path,
    monkeypatch,
):
    cli = make_cli(tmp_path)

    cli._fg_detected = True

    ensure_calls = []

    monkeypatch.setattr(
        cli,
        "_ensure_fg_detected",
        lambda: ensure_calls.append(True),
    )

    expected = object()

    class FakeDetector:
        def functional_group_highlighted_molecules_image_grid(
            self,
            session,
        ):
            return expected

    monkeypatch.setattr(
        cli_module,
        "FunctionalGroupsDetector",
        FakeDetector,
    )

    result = cli.show_functional_groups()

    assert result is expected
    assert ensure_calls == []


# =============================================================================
# show_reactions
# =============================================================================


def test_show_reactions_triggers_detection_when_needed(
    tmp_path,
    monkeypatch,
):
    cli = make_cli(tmp_path)

    expected = object()
    ensure_calls = []

    def fake_ensure():
        ensure_calls.append(True)
        cli._reactions_detected = True

    monkeypatch.setattr(
        cli,
        "_ensure_reactions_detected",
        fake_ensure,
    )

    class FakeDetector:
        def available_reaction_image_grid(
            self,
            session,
        ):
            assert session is cli.session
            return expected

    monkeypatch.setattr(
        cli_module,
        "ReactionDetector",
        FakeDetector,
    )

    result = cli.show_reactions()

    assert result is expected
    assert ensure_calls == [
        True
    ]


def test_show_reactions_does_not_redetect_when_already_detected(
    tmp_path,
    monkeypatch,
):
    cli = make_cli(tmp_path)

    cli._reactions_detected = True

    ensure_calls = []

    monkeypatch.setattr(
        cli,
        "_ensure_reactions_detected",
        lambda: ensure_calls.append(True),
    )

    expected = object()

    class FakeDetector:
        def available_reaction_image_grid(
            self,
            session,
        ):
            return expected

    monkeypatch.setattr(
        cli_module,
        "ReactionDetector",
        FakeDetector,
    )

    result = cli.show_reactions()

    assert result is expected
    assert ensure_calls == []


# =============================================================================
# select_reactions
# =============================================================================


def test_select_reactions_ensures_detection_first(
    tmp_path,
    monkeypatch,
):
    cli = make_cli(tmp_path)

    cli.session.reaction_instances = [
        object()
    ]

    calls = []

    def fake_ensure():
        calls.append("detect")
        cli._reactions_detected = True

    monkeypatch.setattr(
        cli,
        "_ensure_reactions_detected",
        fake_ensure,
    )

    class FakeDetector:
        def reaction_selection(
            self,
            session,
        ):
            assert session is cli.session
            calls.append("select")

    monkeypatch.setattr(
        cli_module,
        "ReactionDetector",
        FakeDetector,
    )

    cli.select_reactions()

    assert calls == [
        "detect",
        "select",
    ]

    assert cli._reactions_selected is True

    assert (
        cli.error_handler[
            "select_reactions"
        ]
        is True
    )


def test_select_reactions_raises_when_none_detected(
    tmp_path,
):
    cli = make_cli(tmp_path)

    cli._reactions_detected = True
    cli.session.reaction_instances = []

    with pytest.raises(
        RuntimeError,
        match="No reactions detected",
    ):
        cli.select_reactions()

    assert (
        cli._reactions_selected
        is False
    )

    assert (
        cli.error_handler[
            "select_reactions"
        ]
        is False
    )


def test_select_reactions_is_idempotent(
    tmp_path,
    monkeypatch,
):
    cli = make_cli(tmp_path)

    cli._reactions_detected = True
    cli._reactions_selected = True

    cli.session.reaction_instances = [
        object()
    ]

    calls = []

    class FakeDetector:
        def reaction_selection(
            self,
            session,
        ):
            calls.append(True)

    monkeypatch.setattr(
        cli_module,
        "ReactionDetector",
        FakeDetector,
    )

    cli.select_reactions()
    cli.select_reactions()

    assert calls == []


# =============================================================================
# show_non_reactants
# =============================================================================


def test_show_non_reactants_returns_none_when_no_visualization(
    tmp_path,
    monkeypatch,
):
    cli = make_cli(tmp_path)

    ensure_calls = []
    save_calls = []

    monkeypatch.setattr(
        cli,
        "_ensure_non_reactants_detected",
        lambda: ensure_calls.append(True),
    )

    class FakeDetector:
        def non_reactants_to_visualization(
            self,
            session,
        ):
            assert session is cli.session
            return None

    monkeypatch.setattr(
        cli_module,
        "NonReactantsDetector",
        FakeDetector,
    )

    monkeypatch.setattr(
        cli,
        "_save_rdkit_img",
        lambda *args, **kwargs:
        save_calls.append(True),
    )

    result = cli.show_non_reactants()

    assert result is None

    assert ensure_calls == [
        True
    ]

    assert save_calls == []


def test_show_non_reactants_saves_visualization(
    tmp_path,
    monkeypatch,
):
    cli = make_cli(tmp_path)

    image = object()
    saved = []

    monkeypatch.setattr(
        cli,
        "_ensure_non_reactants_detected",
        lambda: None,
    )

    class FakeDetector:
        def non_reactants_to_visualization(
            self,
            session,
        ):
            return image

    monkeypatch.setattr(
        cli_module,
        "NonReactantsDetector",
        FakeDetector,
    )

    def fake_save(
        img,
        path,
        is_non_reactant=False,
    ):
        saved.append(
            (
                img,
                Path(path),
                is_non_reactant,
            )
        )

    monkeypatch.setattr(
        cli,
        "_save_rdkit_img",
        fake_save,
    )

    result = cli.show_non_reactants()

    assert result is image

    assert saved == [
        (
            image,
            cli.img_dir
            / "non_reactants.png",
            True,
        )
    ]


# =============================================================================
# select_non_reactants
# =============================================================================


def test_select_non_reactants_calls_selection_when_candidates_exist(
    tmp_path,
    monkeypatch,
):
    cli = make_cli(tmp_path)

    cli.session.non_reactants = [
        object()
    ]

    calls = []

    monkeypatch.setattr(
        cli,
        "_ensure_non_reactants_detected",
        lambda: calls.append("detect"),
    )

    class FakeDetector:
        def non_reactant_selection(
            self,
            session,
        ):
            assert session is cli.session
            calls.append("select")

    monkeypatch.setattr(
        cli_module,
        "NonReactantsDetector",
        FakeDetector,
    )

    cli.select_non_reactants()

    assert calls == [
        "detect",
        "select",
    ]

    assert (
        cli._non_reactants_selected
        is True
    )

    assert (
        cli.error_handler[
            "select_non_reactants"
        ]
        is True
    )


def test_select_non_reactants_completes_when_none_exist(
    tmp_path,
    monkeypatch,
):
    cli = make_cli(tmp_path)

    cli.session.non_reactants = []

    selection_calls = []

    monkeypatch.setattr(
        cli,
        "_ensure_non_reactants_detected",
        lambda: None,
    )

    class FakeDetector:
        def non_reactant_selection(
            self,
            session,
        ):
            selection_calls.append(True)

    monkeypatch.setattr(
        cli_module,
        "NonReactantsDetector",
        FakeDetector,
    )

    cli.select_non_reactants()

    assert selection_calls == []

    assert (
        cli._non_reactants_selected
        is True
    )

    assert (
        cli.error_handler[
            "select_non_reactants"
        ]
        is True
    )


def test_select_non_reactants_is_idempotent(
    tmp_path,
    monkeypatch,
):
    cli = make_cli(tmp_path)

    cli._non_reactants_selected = True
    cli.session.non_reactants = [
        object()
    ]

    calls = []

    monkeypatch.setattr(
        cli,
        "_ensure_non_reactants_detected",
        lambda: calls.append("detect"),
    )

    class FakeDetector:
        def non_reactant_selection(
            self,
            session,
        ):
            calls.append("select")

    monkeypatch.setattr(
        cli_module,
        "NonReactantsDetector",
        FakeDetector,
    )

    cli.select_non_reactants()

    # Detection helper is called first by the current implementation,
    # but selection itself must not happen again.
    assert calls == [
        "detect"
    ]


# =============================================================================
# prepare_reactions
# =============================================================================


def test_prepare_reactions_delegates_and_marks_process_complete(
    tmp_path,
    monkeypatch,
):
    cli = make_cli(tmp_path)

    calls = []

    class FakePrepareReactions:
        def __init__(
            self,
            session,
        ):
            assert session is cli.session
            calls.append("init")

        def prepare_reactions(
            self,
            session,
        ):
            assert session is cli.session
            calls.append("prepare")

    monkeypatch.setattr(
        cli_module,
        "PrepareReactions",
        FakePrepareReactions,
    )

    result = cli.prepare_reactions()

    assert result is None

    assert calls == [
        "init",
        "prepare",
    ]

    assert (
        cli.error_handler["process"]
        is True
    )


def test_show_reaction_templates_delegates_highlight_type(
    tmp_path,
    monkeypatch,
):
    cli = make_cli(tmp_path)

    expected = object()
    calls = []

    class FakePrepareReactions:
        def __init__(
            self,
            session,
        ):
            assert session is cli.session

        def reaction_templates_highlighted_image_grid(
            self,
            session,
            highlight_type,
        ):
            calls.append(
                (
                    session,
                    highlight_type,
                )
            )
            return expected

    monkeypatch.setattr(
        cli_module,
        "PrepareReactions",
        FakePrepareReactions,
    )

    result = cli.show_reaction_templates(
        highlight_type="edge"
    )

    assert result is expected

    assert calls == [
        (
            cli.session,
            "edge",
        )
    ]


# =============================================================================
# process guards
# =============================================================================


def test_process_requires_reaction_selection(
    tmp_path,
):
    cli = make_cli(tmp_path)

    with pytest.raises(
        RuntimeError,
        match="Reactions have not been selected",
    ):
        cli.process()


def test_process_requires_non_reactant_selection(
    tmp_path,
):
    cli = make_cli(tmp_path)

    cli.error_handler[
        "select_reactions"
    ] = True

    with pytest.raises(
        RuntimeError,
        match="Non-reactants have not been selected",
    ):
        cli.process()


def test_process_requires_reaction_preparation(
    tmp_path,
):
    cli = make_cli(tmp_path)

    cli.error_handler[
        "select_reactions"
    ] = True

    cli.error_handler[
        "select_non_reactants"
    ] = True

    with pytest.raises(
        RuntimeError,
        match="Processing has not been completed",
    ):
        cli.process()


# =============================================================================
# process pipeline
# =============================================================================


def test_process_runs_pipeline_in_correct_order(
    tmp_path,
    monkeypatch,
):
    cli = make_cli(tmp_path)

    cli.error_handler[
        "select_reactions"
    ] = True

    cli.error_handler[
        "select_non_reactants"
    ] = True

    cli.error_handler[
        "process"
    ] = True

    calls = []

    @contextmanager
    def fake_writer(
        filename="AutoREACTER.log",
    ):
        calls.append("writer_enter")
        yield
        calls.append("writer_exit")

    monkeypatch.setattr(
        cli,
        "_writer",
        fake_writer,
    )

    class FakeMolecule3DPreparation:
        def __init__(
            self,
            session,
        ):
            assert session is cli.session

        def prepare_molecule_3d_geometry(
            self,
            session,
        ):
            assert session is cli.session
            calls.append("3d")

    class FakeFFWrapper:
        def __init__(
            self,
            session,
        ):
            assert session is cli.session

        def generate_force_field_files(
            self,
            session,
        ):
            assert session is cli.session
            calls.append("ff")

    class FakeREACTERFilesBuilder:
        def __init__(
            self,
            session,
        ):
            assert session is cli.session

        def molecule_template_preparation(
            self,
            session,
        ):
            assert session is cli.session
            calls.append("reacter")

    class FakeSimulationSetupManager:
        def setup_and_write_simulation(
            self,
            session,
        ):
            assert session is cli.session
            calls.append("simulation")

    class FakePrepareReactions:
        def __init__(
            self,
            session,
        ):
            assert session is cli.session

        def reaction_templates_highlighted_image_grid(
            self,
            session,
            highlight_type,
        ):
            assert session is cli.session

            calls.append(
                f"highlight:{highlight_type}"
            )

            return (
                f"image-{highlight_type}"
            )

    monkeypatch.setattr(
        cli_module,
        "Molecule3DPreparation",
        FakeMolecule3DPreparation,
    )

    monkeypatch.setattr(
        cli_module,
        "FFWrapper",
        FakeFFWrapper,
    )

    monkeypatch.setattr(
        cli_module,
        "REACTERFilesBuilder",
        FakeREACTERFilesBuilder,
    )

    monkeypatch.setattr(
        cli_module,
        "SimulationSetupManager",
        FakeSimulationSetupManager,
    )

    monkeypatch.setattr(
        cli_module,
        "PrepareReactions",
        FakePrepareReactions,
    )

    def fake_save(
        image,
        path,
        is_non_reactant=False,
    ):
        calls.append(
            f"save:{Path(path).name}"
        )

    monkeypatch.setattr(
        cli,
        "_save_rdkit_img",
        fake_save,
    )

    cli.process()

    assert calls == [
        "writer_enter",
        "3d",
        "ff",
        "reacter",
        "simulation",
        "highlight:template",
        "save:templates_template.png",
        "highlight:edge",
        "save:templates_edge.png",
        "highlight:initiators",
        "save:templates_initiators.png",
        "highlight:delete",
        "save:templates_delete.png",
        "writer_exit",
    ]

    assert (
        cli.error_handler["process"]
        is True
    )


# =============================================================================
# _save_input_json
# =============================================================================


def test_save_input_json_copies_to_output_directory(
    tmp_path,
):
    cli = make_cli(tmp_path)

    source = tmp_path / "original.json"

    source.write_text(
        '{"hello": "world"}',
        encoding="utf-8",
    )

    cli._save_input_json(source)

    destination = (
        cli.session.output_dir
        / "input.json"
    )

    assert destination.exists()

    assert (
        destination.read_text(
            encoding="utf-8"
        )
        == '{"hello": "world"}'
    )


def test_save_input_json_overwrites_existing_copy(
    tmp_path,
):
    cli = make_cli(tmp_path)

    source = tmp_path / "source.json"

    source.write_text(
        "new content",
        encoding="utf-8",
    )

    destination = (
        cli.session.output_dir
        / "input.json"
    )

    destination.write_text(
        "old content",
        encoding="utf-8",
    )

    cli._save_input_json(source)

    assert (
        destination.read_text(
            encoding="utf-8"
        )
        == "new content"
    )


# =============================================================================
# _ensure_fg_detected
# =============================================================================


def test_ensure_fg_detected_runs_detector_once(
    tmp_path,
    monkeypatch,
):
    cli = make_cli(tmp_path)

    calls = []

    class FakeDetector:
        def functional_groups_detector(
            self,
            session,
        ):
            assert session is cli.session
            calls.append("detect")

    monkeypatch.setattr(
        cli_module,
        "FunctionalGroupsDetector",
        FakeDetector,
    )

    cli._ensure_fg_detected()
    cli._ensure_fg_detected()

    assert calls == [
        "detect"
    ]

    assert cli._fg_detected is True


def test_ensure_fg_detected_sets_flag_only_after_success(
    tmp_path,
    monkeypatch,
):
    cli = make_cli(tmp_path)

    class FakeDetector:
        def functional_groups_detector(
            self,
            session,
        ):
            raise RuntimeError(
                "detector failed"
            )

    monkeypatch.setattr(
        cli_module,
        "FunctionalGroupsDetector",
        FakeDetector,
    )

    with pytest.raises(
        RuntimeError,
        match="detector failed",
    ):
        cli._ensure_fg_detected()

    assert cli._fg_detected is False


# =============================================================================
# _ensure_reactions_detected
# =============================================================================


def test_ensure_reactions_detected_runs_reaction_detection_once(
    tmp_path,
    monkeypatch,
):
    cli = make_cli(tmp_path)

    calls = []
    saved = []

    fg_image = object()
    reaction_image = object()

    class FakeFGDetector:
        def functional_groups_detector(
            self,
            session,
        ):
            calls.append("fg_detect")

        def functional_group_highlighted_molecules_image_grid(
            self,
            session,
        ):
            calls.append("fg_image")
            return fg_image

    class FakeReactionDetector:
        def reaction_detector(
            self,
            session,
        ):
            calls.append("reaction_detect")

        def available_reaction_image_grid(
            self,
            session,
        ):
            calls.append("reaction_image")
            return reaction_image

    monkeypatch.setattr(
        cli_module,
        "FunctionalGroupsDetector",
        FakeFGDetector,
    )

    monkeypatch.setattr(
        cli_module,
        "ReactionDetector",
        FakeReactionDetector,
    )

    def fake_save(
        image,
        path,
        is_non_reactant=False,
    ):
        saved.append(
            (
                image,
                Path(path).name,
            )
        )

    monkeypatch.setattr(
        cli,
        "_save_rdkit_img",
        fake_save,
    )

    cli._ensure_reactions_detected()
    cli._ensure_reactions_detected()

    assert calls.count(
        "fg_detect"
    ) == 1

    assert calls.count(
        "reaction_detect"
    ) == 1

    assert calls.count(
        "reaction_image"
    ) == 1

    # Current implementation regenerates/saves the FG visualization
    # whenever this private helper is called, even after reaction detection
    # has already completed. Reaction detection itself remains idempotent.
    assert calls.count(
        "fg_image"
    ) == 2

    assert cli._fg_detected is True

    assert (
        cli._reactions_detected
        is True
    )

    assert (
        (
            reaction_image,
            "reactions.png",
        )
        in saved
    )


def test_ensure_reactions_detected_does_not_mark_success_after_detector_error(
    tmp_path,
    monkeypatch,
):
    cli = make_cli(tmp_path)

    class FakeFGDetector:
        def functional_groups_detector(
            self,
            session,
        ):
            pass

        def functional_group_highlighted_molecules_image_grid(
            self,
            session,
        ):
            return object()

    class FakeReactionDetector:
        def reaction_detector(
            self,
            session,
        ):
            raise RuntimeError(
                "reaction failure"
            )

    monkeypatch.setattr(
        cli_module,
        "FunctionalGroupsDetector",
        FakeFGDetector,
    )

    monkeypatch.setattr(
        cli_module,
        "ReactionDetector",
        FakeReactionDetector,
    )

    monkeypatch.setattr(
        cli,
        "_save_rdkit_img",
        lambda *args, **kwargs: None,
    )

    with pytest.raises(
        RuntimeError,
        match="reaction failure",
    ):
        cli._ensure_reactions_detected()

    assert cli._fg_detected is True

    assert (
        cli._reactions_detected
        is False
    )


# =============================================================================
# _ensure_non_reactants_detected
# =============================================================================


def test_ensure_non_reactants_detected_selects_reactions_first(
    tmp_path,
    monkeypatch,
):
    cli = make_cli(tmp_path)

    calls = []

    def fake_select_reactions():
        calls.append("select_reactions")
        cli._reactions_selected = True

    monkeypatch.setattr(
        cli,
        "select_reactions",
        fake_select_reactions,
    )

    class FakeDetector:
        def non_monomer_detector(
            self,
            session,
        ):
            assert session is cli.session
            calls.append(
                "detect_non_reactants"
            )

    monkeypatch.setattr(
        cli_module,
        "NonReactantsDetector",
        FakeDetector,
    )

    cli._ensure_non_reactants_detected()

    assert calls == [
        "select_reactions",
        "detect_non_reactants",
    ]

    assert (
        cli._non_reactants_detected
        is True
    )


def test_ensure_non_reactants_detected_is_idempotent(
    tmp_path,
    monkeypatch,
):
    cli = make_cli(tmp_path)

    cli._reactions_selected = True

    calls = []

    class FakeDetector:
        def non_monomer_detector(
            self,
            session,
        ):
            calls.append(True)

    monkeypatch.setattr(
        cli_module,
        "NonReactantsDetector",
        FakeDetector,
    )

    cli._ensure_non_reactants_detected()
    cli._ensure_non_reactants_detected()

    assert calls == [
        True
    ]


# =============================================================================
# _save_rdkit_img
# =============================================================================


def test_save_rdkit_img_creates_parent_directories(
    tmp_path,
):
    cli = make_cli(tmp_path)

    image = PILImage.new(
        "RGB",
        (2, 2),
    )

    path = (
        tmp_path
        / "nested"
        / "more"
        / "image.png"
    )

    cli._save_rdkit_img(
        image,
        path,
    )

    assert path.exists()


def test_save_rdkit_img_saves_pil_image(
    tmp_path,
):
    cli = make_cli(tmp_path)

    image = PILImage.new(
        "RGB",
        (3, 3),
    )

    path = tmp_path / "image.png"

    cli._save_rdkit_img(
        image,
        path,
    )

    assert path.exists()

    loaded = PILImage.open(path)

    assert loaded.size == (
        3,
        3,
    )


@pytest.mark.parametrize(
    "data",
    [
        b"\x89PNGtest",
        bytearray(b"\x89PNGtest"),
    ],
)
def test_save_rdkit_img_writes_raw_bytes(
    tmp_path,
    data,
):
    cli = make_cli(tmp_path)

    path = tmp_path / "raw.bin"

    cli._save_rdkit_img(
        data,
        path,
    )

    assert (
        path.read_bytes()
        == bytes(data)
    )


def test_save_rdkit_img_writes_image_data_bytes(
    tmp_path,
):
    cli = make_cli(tmp_path)

    image = SimpleNamespace(
        data=b"image-data"
    )

    path = tmp_path / "image.bin"

    cli._save_rdkit_img(
        image,
        path,
    )

    assert (
        path.read_bytes()
        == b"image-data"
    )


def test_save_rdkit_img_writes_image_data_string(
    tmp_path,
):
    cli = make_cli(tmp_path)

    image = SimpleNamespace(
        data="<svg>hello</svg>"
    )

    path = tmp_path / "image.svg"

    cli._save_rdkit_img(
        image,
        path,
    )

    assert (
        path.read_text(
            encoding="utf-8"
        )
        == "<svg>hello</svg>"
    )


def test_save_rdkit_img_writes_direct_svg_string(
    tmp_path,
):
    cli = make_cli(tmp_path)

    path = tmp_path / "image.svg"

    cli._save_rdkit_img(
        "<svg>direct</svg>",
        path,
    )

    assert (
        path.read_text(
            encoding="utf-8"
        )
        == "<svg>direct</svg>"
    )


def test_save_rdkit_img_allows_none_for_non_reactants(
    tmp_path,
):
    cli = make_cli(tmp_path)

    path = tmp_path / "none.png"

    result = cli._save_rdkit_img(
        None,
        path,
        is_non_reactant=True,
    )

    assert result is None
    assert not path.exists()


def test_save_rdkit_img_rejects_none_for_required_image(
    tmp_path,
):
    cli = make_cli(tmp_path)

    with pytest.raises(
        NoReactionGenerated,
        match="No reaction was generated",
    ):
        cli._save_rdkit_img(
            None,
            tmp_path / "missing.png",
        )


def test_save_rdkit_img_rejects_unsupported_type(
    tmp_path,
):
    cli = make_cli(tmp_path)

    with pytest.raises(
        TypeError,
        match="Unsupported image type",
    ):
        cli._save_rdkit_img(
            12345,
            tmp_path / "bad.png",
        )


# =============================================================================
# _writer
# =============================================================================


def test_writer_writes_python_and_os_output_to_log(
    tmp_path,
    capfd,
):
    cli = make_cli(tmp_path)

    with cli._writer():
        print(
            "python-output",
            flush=True,
        )

        os.write(
            1,
            b"os-output\n",
        )

    log_path = (
        cli.session.output_dir
        / "AutoREACTER.log"
    )

    log_text = log_path.read_text(
        encoding="utf-8"
    )

    assert "python-output" in log_text
    assert "os-output" in log_text

    terminal_output = (
        capfd.readouterr().out
    )

    assert "python-output" in terminal_output
    assert "os-output" in terminal_output


def test_writer_supports_custom_filename(
    tmp_path,
):
    cli = make_cli(tmp_path)

    with cli._writer(
        "custom.log"
    ):
        print(
            "custom-output",
            flush=True,
        )

    log_path = (
        cli.session.output_dir
        / "custom.log"
    )

    assert log_path.exists()

    assert (
        "custom-output"
        in log_path.read_text(
            encoding="utf-8"
        )
    )


def test_writer_restores_stdout_after_exception(
    tmp_path,
    capfd,
):
    cli = make_cli(tmp_path)

    with pytest.raises(
        RuntimeError,
        match="boom",
    ):
        with cli._writer(
            "exception.log"
        ):
            print(
                "before-exception",
                flush=True,
            )

            raise RuntimeError(
                "boom"
            )

    print(
        "after-exception",
        flush=True,
    )

    terminal_output = (
        capfd.readouterr().out
    )

    assert (
        "after-exception"
        in terminal_output
    )

    log_path = (
        cli.session.output_dir
        / "exception.log"
    )

    log_text = log_path.read_text(
        encoding="utf-8"
    )

    assert (
        "before-exception"
        in log_text
    )

    assert (
        "after-exception"
        not in log_text
    )


def test_writer_appends_to_existing_log(
    tmp_path,
):
    cli = make_cli(tmp_path)

    with cli._writer():
        print(
            "first-run",
            flush=True,
        )

    with cli._writer():
        print(
            "second-run",
            flush=True,
        )

    log_path = (
        cli.session.output_dir
        / "AutoREACTER.log"
    )

    log_text = log_path.read_text(
        encoding="utf-8"
    )

    assert "first-run" in log_text
    assert "second-run" in log_text


# =============================================================================
# __repr__
# =============================================================================


def test_repr_with_input_path(
    tmp_path,
):
    cli = ARXCLI.__new__(
        ARXCLI
    )

    cli.input = (
        tmp_path / "input.json"
    )

    assert repr(cli) == (
        f"ARXCLI(input={cli.input})"
    )


def test_repr_without_input():
    cli = ARXCLI.__new__(
        ARXCLI
    )

    assert repr(cli) == "ARXCLI()"