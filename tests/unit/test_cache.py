from pathlib import Path
from types import SimpleNamespace

import pytest

import AutoREACTER.cache as cache_module
from AutoREACTER.cache import (
    GetCacheDir,
    RunDirectoryManager,
)


# =============================================================================
# Helpers
# =============================================================================


def make_fake_reacter_files(old_base: Path):
    """
    Build a lightweight REACTERFiles-like object.

    The production function only requires these attributes, so using
    SimpleNamespace keeps cache tests independent from the force-field and
    reaction-template builders.
    """
    force_field_data = old_base / "force_field.data"
    in_file = old_base / "in.create_atoms.script"
    molecule_file = old_base / "monomer.molecule"
    map_file = old_base / "RXN_1.map"
    pre_file = old_base / "template_pre_1.molecule"
    post_file = old_base / "template_post_1.molecule"

    for path in [
        force_field_data,
        in_file,
        molecule_file,
        map_file,
        pre_file,
        post_file,
    ]:
        path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )
        path.write_text(
            path.name,
            encoding="utf-8",
        )

    molecule = SimpleNamespace(
        molecule_files=SimpleNamespace(
            lmp_molecule_file=molecule_file,
        )
    )

    template = SimpleNamespace(
        map_file=map_file,
        pre_reaction_file=SimpleNamespace(
            lmp_molecule_file=pre_file,
        ),
        post_reaction_file=SimpleNamespace(
            lmp_molecule_file=post_file,
        ),
    )

    return SimpleNamespace(
        force_field_data=force_field_data,
        in_file=in_file,
        molecule_files=[molecule],
        template_files=[template],
    )


# =============================================================================
# GetCacheDir construction
# =============================================================================


def test_get_cache_dir_creates_staging_directory(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        cache_module.tempfile,
        "gettempdir",
        lambda: str(tmp_path),
    )

    cache = GetCacheDir()

    assert cache.staging_dir.exists()
    assert cache.staging_dir.is_dir()
    assert cache.staging_dir.parent == tmp_path


def test_get_cache_dir_default_clears_existing_staging_contents(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        cache_module.tempfile,
        "gettempdir",
        lambda: str(tmp_path),
    )

    staging = (
        tmp_path
        / "AutoREACTER_staging"
    )
    staging.mkdir()

    old_file = staging / "old.txt"
    old_file.write_text(
        "old",
        encoding="utf-8",
    )

    old_directory = staging / "nested"
    old_directory.mkdir()

    (
        old_directory / "nested.txt"
    ).write_text(
        "nested",
        encoding="utf-8",
    )

    cache = GetCacheDir()

    assert cache.staging_dir == staging
    assert staging.exists()
    assert list(staging.iterdir()) == []


def test_get_cache_dir_clear_staging_true_clears_contents(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        cache_module.tempfile,
        "gettempdir",
        lambda: str(tmp_path),
    )

    staging = (
        tmp_path
        / "AutoREACTER_staging"
    )
    staging.mkdir()

    old_file = staging / "old.txt"
    old_file.write_text(
        "old",
        encoding="utf-8",
    )

    GetCacheDir(
        clear_staging=True
    )

    assert not old_file.exists()


def test_get_cache_dir_clear_staging_false_preserves_contents(
    tmp_path,
    monkeypatch,
):
    """
    The constructor documents clear_staging=False as disabling automatic
    staging cleanup. Existing temporary data must therefore remain untouched.
    """
    monkeypatch.setattr(
        cache_module.tempfile,
        "gettempdir",
        lambda: str(tmp_path),
    )

    staging = (
        tmp_path
        / "AutoREACTER_staging"
    )
    staging.mkdir()

    old_file = staging / "keep.txt"
    old_file.write_text(
        "keep me",
        encoding="utf-8",
    )

    cache = GetCacheDir(
        clear_staging=False
    )

    assert cache.staging_dir == staging
    assert old_file.exists()

    assert (
        old_file.read_text(
            encoding="utf-8"
        )
        == "keep me"
    )


# =============================================================================
# clear_staging_dir
# =============================================================================


def test_clear_staging_dir_removes_files(
    tmp_path,
):
    cache = GetCacheDir.__new__(
        GetCacheDir
    )

    cache.staging_dir = (
        tmp_path / "staging"
    )
    cache.staging_dir.mkdir()

    file_path = (
        cache.staging_dir / "file.txt"
    )
    file_path.write_text(
        "data",
        encoding="utf-8",
    )

    cache.clear_staging_dir()

    assert not file_path.exists()
    assert cache.staging_dir.exists()


def test_clear_staging_dir_removes_nested_directories(
    tmp_path,
):
    cache = GetCacheDir.__new__(
        GetCacheDir
    )

    cache.staging_dir = (
        tmp_path / "staging"
    )

    nested = (
        cache.staging_dir
        / "a"
        / "b"
        / "c"
    )
    nested.mkdir(
        parents=True
    )

    (
        nested / "data.txt"
    ).write_text(
        "data",
        encoding="utf-8",
    )

    cache.clear_staging_dir()

    assert cache.staging_dir.exists()
    assert list(
        cache.staging_dir.iterdir()
    ) == []


def test_clear_staging_dir_creates_missing_root(
    tmp_path,
):
    cache = GetCacheDir.__new__(
        GetCacheDir
    )

    cache.staging_dir = (
        tmp_path
        / "missing_staging"
    )

    assert not cache.staging_dir.exists()

    cache.clear_staging_dir()

    assert cache.staging_dir.is_dir()


def test_clear_staging_dir_prints_success(
    tmp_path,
    capsys,
):
    cache = GetCacheDir.__new__(
        GetCacheDir
    )

    cache.staging_dir = (
        tmp_path / "staging"
    )
    cache.staging_dir.mkdir()

    cache.clear_staging_dir()

    output = capsys.readouterr().out

    assert (
        "[OK] Cleared staging cache:"
        in output
    )

    assert str(
        cache.staging_dir
    ) in output


def test_clear_staging_dir_continues_after_failed_item(
    tmp_path,
    monkeypatch,
    capsys,
):
    cache = GetCacheDir.__new__(
        GetCacheDir
    )

    cache.staging_dir = (
        tmp_path / "staging"
    )
    cache.staging_dir.mkdir()

    failing_dir = (
        cache.staging_dir
        / "cannot_remove"
    )
    failing_dir.mkdir()

    removable_file = (
        cache.staging_dir
        / "remove_me.txt"
    )
    removable_file.write_text(
        "remove",
        encoding="utf-8",
    )

    real_rmtree = (
        cache_module.shutil.rmtree
    )

    def fake_rmtree(path, *args, **kwargs):
        if Path(path) == failing_dir:
            raise PermissionError(
                "simulated failure"
            )

        return real_rmtree(
            path,
            *args,
            **kwargs,
        )

    monkeypatch.setattr(
        cache_module.shutil,
        "rmtree",
        fake_rmtree,
    )

    cache.clear_staging_dir()

    assert failing_dir.exists()
    assert not removable_file.exists()

    output = capsys.readouterr().out

    assert (
        "[WARN] Failed to remove staging cache item"
        in output
    )

    assert (
        "[WARN] Staging cache partially cleared"
        in output
    )


# =============================================================================
# RunDirectoryManager construction
# =============================================================================


def test_run_directory_manager_creates_base_directory(
    tmp_path,
):
    base_dir = (
        tmp_path
        / "runs"
        / "nested"
    )

    assert not base_dir.exists()

    manager = RunDirectoryManager(
        base_dir
    )

    assert manager.base_dir == base_dir
    assert base_dir.is_dir()


def test_run_directory_manager_accepts_string_path(
    tmp_path,
):
    base_dir = tmp_path / "runs"

    manager = RunDirectoryManager(
        str(base_dir)
    )

    assert isinstance(
        manager.base_dir,
        Path,
    )

    assert (
        manager.base_dir
        == base_dir
    )


# =============================================================================
# remove_path
# =============================================================================


def test_remove_path_removes_file(
    tmp_path,
):
    manager = RunDirectoryManager(
        tmp_path / "runs"
    )

    path = tmp_path / "file.txt"
    path.write_text(
        "data",
        encoding="utf-8",
    )

    manager.remove_path(path)

    assert not path.exists()


def test_remove_path_removes_directory_recursively(
    tmp_path,
):
    manager = RunDirectoryManager(
        tmp_path / "runs"
    )

    directory = tmp_path / "folder"
    nested = directory / "nested"
    nested.mkdir(
        parents=True
    )

    (
        nested / "data.txt"
    ).write_text(
        "data",
        encoding="utf-8",
    )

    manager.remove_path(
        directory
    )

    assert not directory.exists()


def test_remove_path_nonexistent_path_is_noop(
    tmp_path,
):
    manager = RunDirectoryManager(
        tmp_path / "runs"
    )

    missing = (
        tmp_path / "missing"
    )

    manager.remove_path(missing)

    assert not missing.exists()


# =============================================================================
# move_into_run
# =============================================================================


def test_move_into_run_moves_all_contents(
    tmp_path,
):
    manager = RunDirectoryManager(
        tmp_path / "runs"
    )

    source = tmp_path / "source"
    destination = (
        tmp_path / "destination"
    )

    source.mkdir()
    destination.mkdir()

    (
        source / "one.txt"
    ).write_text(
        "one",
        encoding="utf-8",
    )

    (
        source / "two.txt"
    ).write_text(
        "two",
        encoding="utf-8",
    )

    result = manager.move_into_run(
        source,
        destination,
    )

    assert result == destination

    assert (
        destination / "one.txt"
    ).read_text(
        encoding="utf-8"
    ) == "one"

    assert (
        destination / "two.txt"
    ).read_text(
        encoding="utf-8"
    ) == "two"

    assert list(source.iterdir()) == []


def test_move_into_run_moves_directories(
    tmp_path,
):
    manager = RunDirectoryManager(
        tmp_path / "runs"
    )

    source = tmp_path / "source"
    destination = (
        tmp_path / "destination"
    )

    source.mkdir()
    destination.mkdir()

    nested = source / "nested"
    nested.mkdir()

    (
        nested / "data.txt"
    ).write_text(
        "nested data",
        encoding="utf-8",
    )

    manager.move_into_run(
        source,
        destination,
    )

    moved_file = (
        destination
        / "nested"
        / "data.txt"
    )

    assert moved_file.exists()

    assert (
        moved_file.read_text(
            encoding="utf-8"
        )
        == "nested data"
    )


def test_move_into_run_overwrites_existing_file(
    tmp_path,
):
    manager = RunDirectoryManager(
        tmp_path / "runs"
    )

    source = tmp_path / "source"
    destination = (
        tmp_path / "destination"
    )

    source.mkdir()
    destination.mkdir()

    source_file = (
        source / "same.txt"
    )
    source_file.write_text(
        "new",
        encoding="utf-8",
    )

    existing = (
        destination / "same.txt"
    )
    existing.write_text(
        "old",
        encoding="utf-8",
    )

    manager.move_into_run(
        source,
        destination,
    )

    assert (
        existing.read_text(
            encoding="utf-8"
        )
        == "new"
    )


def test_move_into_run_overwrites_existing_directory(
    tmp_path,
):
    manager = RunDirectoryManager(
        tmp_path / "runs"
    )

    source = tmp_path / "source"
    destination = (
        tmp_path / "destination"
    )

    source.mkdir()
    destination.mkdir()

    incoming = source / "folder"
    incoming.mkdir()

    (
        incoming / "new.txt"
    ).write_text(
        "new",
        encoding="utf-8",
    )

    existing = (
        destination / "folder"
    )
    existing.mkdir()

    (
        existing / "old.txt"
    ).write_text(
        "old",
        encoding="utf-8",
    )

    manager.move_into_run(
        source,
        destination,
    )

    assert not (
        destination
        / "folder"
        / "old.txt"
    ).exists()

    assert (
        destination
        / "folder"
        / "new.txt"
    ).exists()


def test_move_into_run_preserves_unrelated_destination_files(
    tmp_path,
):
    manager = RunDirectoryManager(
        tmp_path / "runs"
    )

    source = tmp_path / "source"
    destination = (
        tmp_path / "destination"
    )

    source.mkdir()
    destination.mkdir()

    (
        source / "incoming.txt"
    ).write_text(
        "incoming",
        encoding="utf-8",
    )

    unrelated = (
        destination / "keep.txt"
    )
    unrelated.write_text(
        "keep",
        encoding="utf-8",
    )

    manager.move_into_run(
        source,
        destination,
    )

    assert unrelated.exists()

    assert (
        unrelated.read_text(
            encoding="utf-8"
        )
        == "keep"
    )


def test_move_into_run_prints_destination(
    tmp_path,
    capsys,
):
    manager = RunDirectoryManager(
        tmp_path / "runs"
    )

    source = tmp_path / "source"
    destination = (
        tmp_path / "destination"
    )

    source.mkdir()
    destination.mkdir()

    manager.move_into_run(
        source,
        destination,
    )

    output = capsys.readouterr().out

    assert "[OK] Moved files" in output
    assert str(destination) in output


# =============================================================================
# move_reacter_files
# =============================================================================


def test_move_reacter_files_requires_reacter_directory(
    tmp_path,
):
    manager = RunDirectoryManager(
        tmp_path / "runs"
    )

    staging = tmp_path / "staging"
    final = tmp_path / "final"

    staging.mkdir()
    final.mkdir()

    fake_reacter_files = (
        SimpleNamespace()
    )

    expected = (
        staging
        / "lunar"
        / "REACTER_files"
    )

    with pytest.raises(
        FileNotFoundError,
        match="REACTER files not found",
    ):
        manager.move_reacter_files(
            fake_reacter_files,
            staging,
            final,
        )

    assert not expected.exists()


def test_move_reacter_files_moves_physical_files(
    tmp_path,
):
    manager = RunDirectoryManager(
        tmp_path / "runs"
    )

    staging = tmp_path / "staging"
    final = tmp_path / "final"

    old_base = (
        staging
        / "lunar"
        / "REACTER_files"
    )

    old_base.mkdir(
        parents=True
    )
    final.mkdir()

    reacter_files = (
        make_fake_reacter_files(
            old_base
        )
    )

    manager.move_reacter_files(
        reacter_files,
        staging,
        final,
    )

    expected_names = {
        "force_field.data",
        "in.create_atoms.script",
        "monomer.molecule",
        "RXN_1.map",
        "template_pre_1.molecule",
        "template_post_1.molecule",
    }

    assert {
        item.name
        for item in final.iterdir()
    } == expected_names

    assert (
        list(old_base.iterdir())
        == []
    )


def test_move_reacter_files_remaps_run_level_paths(
    tmp_path,
):
    manager = RunDirectoryManager(
        tmp_path / "runs"
    )

    staging = tmp_path / "staging"
    final = tmp_path / "final"

    old_base = (
        staging
        / "lunar"
        / "REACTER_files"
    )

    old_base.mkdir(
        parents=True
    )
    final.mkdir()

    reacter_files = (
        make_fake_reacter_files(
            old_base
        )
    )

    result = manager.move_reacter_files(
        reacter_files,
        staging,
        final,
    )

    assert result is reacter_files

    assert (
        result.force_field_data
        == final / "force_field.data"
    )

    assert (
        result.in_file
        == final / "in.create_atoms.script"
    )


def test_move_reacter_files_remaps_molecule_file_paths(
    tmp_path,
):
    manager = RunDirectoryManager(
        tmp_path / "runs"
    )

    staging = tmp_path / "staging"
    final = tmp_path / "final"

    old_base = (
        staging
        / "lunar"
        / "REACTER_files"
    )

    old_base.mkdir(
        parents=True
    )
    final.mkdir()

    reacter_files = (
        make_fake_reacter_files(
            old_base
        )
    )

    manager.move_reacter_files(
        reacter_files,
        staging,
        final,
    )

    molecule_path = (
        reacter_files
        .molecule_files[0]
        .molecule_files
        .lmp_molecule_file
    )

    assert (
        molecule_path
        == final / "monomer.molecule"
    )


def test_move_reacter_files_remaps_template_paths(
    tmp_path,
):
    manager = RunDirectoryManager(
        tmp_path / "runs"
    )

    staging = tmp_path / "staging"
    final = tmp_path / "final"

    old_base = (
        staging
        / "lunar"
        / "REACTER_files"
    )

    old_base.mkdir(
        parents=True
    )
    final.mkdir()

    reacter_files = (
        make_fake_reacter_files(
            old_base
        )
    )

    manager.move_reacter_files(
        reacter_files,
        staging,
        final,
    )

    template = (
        reacter_files.template_files[0]
    )

    assert (
        template.map_file
        == final / "RXN_1.map"
    )

    assert (
        template
        .pre_reaction_file
        .lmp_molecule_file
        == final / "template_pre_1.molecule"
    )

    assert (
        template
        .post_reaction_file
        .lmp_molecule_file
        == final / "template_post_1.molecule"
    )


def test_move_reacter_files_preserves_paths_outside_old_base(
    tmp_path,
):
    manager = RunDirectoryManager(
        tmp_path / "runs"
    )

    staging = tmp_path / "staging"
    final = tmp_path / "final"

    old_base = (
        staging
        / "lunar"
        / "REACTER_files"
    )

    old_base.mkdir(
        parents=True
    )
    final.mkdir()

    reacter_files = (
        make_fake_reacter_files(
            old_base
        )
    )

    external_file = (
        tmp_path / "external.data"
    )
    external_file.write_text(
        "external",
        encoding="utf-8",
    )

    reacter_files.force_field_data = (
        external_file
    )

    manager.move_reacter_files(
        reacter_files,
        staging,
        final,
    )

    assert (
        reacter_files.force_field_data
        == external_file
    )


def test_move_reacter_files_prints_success(
    tmp_path,
    capsys,
):
    manager = RunDirectoryManager(
        tmp_path / "runs"
    )

    staging = tmp_path / "staging"
    final = tmp_path / "final"

    old_base = (
        staging
        / "lunar"
        / "REACTER_files"
    )

    old_base.mkdir(
        parents=True
    )
    final.mkdir()

    reacter_files = (
        make_fake_reacter_files(
            old_base
        )
    )

    manager.move_reacter_files(
        reacter_files,
        staging,
        final,
    )

    output = capsys.readouterr().out

    assert (
        "[OK] REACTER files moved"
        in output
    )

    assert str(final) in output