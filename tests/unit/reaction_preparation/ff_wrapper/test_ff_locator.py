from pathlib import Path

import pytest

import AutoREACTER.reaction_preparation.ff_wrapper.ff_locator as ff_locator
from AutoREACTER.reaction_preparation.ff_wrapper.ff_locator import (
    get_force_field_file,
)


# =============================================================================
# Helpers
# =============================================================================


def make_fake_internal_ff_tree(
    tmp_path: Path,
    *,
    include_pcff: bool = True,
) -> Path:
    """
    Build the directory that pkg_resources.files("AutoREACTER")
    is expected to expose.
    """
    package_root = (
        tmp_path / "AutoREACTER"
    )

    ff_dir = (
        package_root
        / "reaction_preparation"
        / "ff_wrapper"
        / "FF_files"
    )

    ff_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if include_pcff:
        (
            ff_dir / "pcff.frc"
        ).write_text(
            "fake pcff",
            encoding="utf-8",
        )

    return package_root


def make_fake_lunar(
    tmp_path: Path,
    *,
    files=(),
) -> Path:
    lunar_root = (
        tmp_path / "LUNAR"
    )

    frc_dir = (
        lunar_root / "frc_files"
    )

    frc_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    for filename in files:
        (
            frc_dir / filename
        ).write_text(
            "fake frc",
            encoding="utf-8",
        )

    return lunar_root


def patch_internal_package_root(
    monkeypatch,
    package_root: Path,
):
    monkeypatch.setattr(
        ff_locator.pkg_resources,
        "files",
        lambda package_name: package_root,
    )


# =============================================================================
# Bundled PCFF
# =============================================================================


@pytest.mark.parametrize(
    "force_field",
    [
        "PCFF",
        "PCFF-IFF",
    ],
)
def test_internal_pcff_variants_resolve_to_same_file(
    tmp_path,
    monkeypatch,
    force_field,
):
    package_root = (
        make_fake_internal_ff_tree(
            tmp_path
        )
    )

    patch_internal_package_root(
        monkeypatch,
        package_root,
    )

    result = get_force_field_file(
        force_field
    )

    expected = (
        package_root
        / "reaction_preparation"
        / "ff_wrapper"
        / "FF_files"
        / "pcff.frc"
    )

    assert result == expected


@pytest.mark.parametrize(
    "force_field",
    [
        "PCFF",
        "PCFF-IFF",
    ],
)
def test_internal_pcff_does_not_require_lunar_location(
    tmp_path,
    monkeypatch,
    force_field,
):
    package_root = (
        make_fake_internal_ff_tree(
            tmp_path
        )
    )

    patch_internal_package_root(
        monkeypatch,
        package_root,
    )

    result = get_force_field_file(
        force_field=force_field,
        lunar_location=None,
    )

    assert result.is_file()


def test_internal_pcff_returns_path_object(
    tmp_path,
    monkeypatch,
):
    package_root = (
        make_fake_internal_ff_tree(
            tmp_path
        )
    )

    patch_internal_package_root(
        monkeypatch,
        package_root,
    )

    result = get_force_field_file(
        "PCFF"
    )

    assert isinstance(
        result,
        Path,
    )


def test_internal_pcff_missing_file_raises_value_error(
    tmp_path,
    monkeypatch,
):
    package_root = (
        make_fake_internal_ff_tree(
            tmp_path,
            include_pcff=False,
        )
    )

    patch_internal_package_root(
        monkeypatch,
        package_root,
    )

    expected = (
        package_root
        / "reaction_preparation"
        / "ff_wrapper"
        / "FF_files"
        / "pcff.frc"
    )

    with pytest.raises(
        ValueError,
        match="Force field file not found",
    ) as exc_info:
        get_force_field_file(
            "PCFF"
        )

    assert str(
        expected
    ) in str(
        exc_info.value
    )


def test_real_packaged_pcff_file_is_available():
    """
    Small packaging regression check.

    This intentionally uses the real package resource rather than
    monkeypatching pkg_resources. If pcff.frc disappears from the
    installation/package, AutoREACTER's supported PCFF backend is broken.
    """
    result = get_force_field_file(
        "PCFF"
    )

    assert isinstance(
        result,
        Path,
    )

    assert result.is_file()

    assert result.name == (
        "pcff.frc"
    )


def test_pcff_prefers_internal_file_even_when_lunar_location_given(
    tmp_path,
    monkeypatch,
):
    package_root = (
        make_fake_internal_ff_tree(
            tmp_path
        )
    )

    patch_internal_package_root(
        monkeypatch,
        package_root,
    )

    lunar_root = (
        make_fake_lunar(
            tmp_path,
            files=[
                "compass_published.frc",
            ],
        )
    )

    result = get_force_field_file(
        force_field="PCFF",
        lunar_location=lunar_root,
    )

    assert result == (
        package_root
        / "reaction_preparation"
        / "ff_wrapper"
        / "FF_files"
        / "pcff.frc"
    )


# =============================================================================
# pkg_resources failure
# =============================================================================


def test_pkg_resources_failure_is_wrapped_as_runtime_error(
    monkeypatch,
):
    def fail(package_name):
        raise RuntimeError(
            "resource lookup exploded"
        )

    monkeypatch.setattr(
        ff_locator.pkg_resources,
        "files",
        fail,
    )

    with pytest.raises(
        RuntimeError,
        match="Error locating FF_files directory using pkg_resources",
    ) as exc_info:
        get_force_field_file(
            "PCFF"
        )

    assert (
        "resource lookup exploded"
        in str(
            exc_info.value
        )
    )


@pytest.mark.parametrize(
    "force_field",
    [
        "PCFF",
        "PCFF-IFF",
        "compass",
        "CVFF",
        "CVFF-IFF",
        "DREIDING",
    ],
)
def test_lunar_force_field_branch_attempts_internal_resource_lookup_first(
    monkeypatch,
    force_field,
):
    calls = []

    def fail(package_name):
        calls.append(
            package_name
        )

        raise RuntimeError(
            "boom"
        )

    monkeypatch.setattr(
        ff_locator.pkg_resources,
        "files",
        fail,
    )

    with pytest.raises(
        RuntimeError,
        match="Error locating FF_files",
    ):
        get_force_field_file(
            force_field,
            lunar_location="/tmp/LUNAR",
        )

    assert calls == [
        "AutoREACTER"
    ]


# =============================================================================
# External LUNAR force fields
# =============================================================================


@pytest.mark.parametrize(
    "force_field, expected_filename",
    [
        (
            "compass",
            "compass_published.frc",
        ),
        (
            "CVFF-IFF",
            "cvff_aug.frc",
        ),
        (
            "CVFF",
            "cvff.frc",
        ),
        (
            "DREIDING",
            "all2lmp_dreiding.frc",
        ),
    ],
)
def test_external_lunar_force_field_mapping(
    tmp_path,
    monkeypatch,
    force_field,
    expected_filename,
):
    package_root = (
        make_fake_internal_ff_tree(
            tmp_path
        )
    )

    patch_internal_package_root(
        monkeypatch,
        package_root,
    )

    lunar_root = (
        make_fake_lunar(
            tmp_path,
            files=[
                expected_filename
            ],
        )
    )

    result = get_force_field_file(
        force_field=force_field,
        lunar_location=lunar_root,
    )

    assert result == (
        lunar_root
        / "frc_files"
        / expected_filename
    )

    assert result.is_file()


@pytest.mark.parametrize(
    "force_field",
    [
        "compass",
        "CVFF-IFF",
        "CVFF",
        "DREIDING",
    ],
)
def test_external_lunar_force_field_without_lunar_location_raises(
    tmp_path,
    monkeypatch,
    force_field,
):
    package_root = (
        make_fake_internal_ff_tree(
            tmp_path
        )
    )

    patch_internal_package_root(
        monkeypatch,
        package_root,
    )

    with pytest.raises(
        ValueError,
        match="Required LUNAR installation not found to resolve",
    ) as exc_info:
        get_force_field_file(
            force_field=force_field,
            lunar_location=None,
        )

    assert force_field in str(
        exc_info.value
    )


@pytest.mark.parametrize(
    "force_field, expected_filename",
    [
        (
            "compass",
            "compass_published.frc",
        ),
        (
            "CVFF-IFF",
            "cvff_aug.frc",
        ),
        (
            "CVFF",
            "cvff.frc",
        ),
        (
            "DREIDING",
            "all2lmp_dreiding.frc",
        ),
    ],
)
def test_external_lunar_force_field_missing_file_raises(
    tmp_path,
    monkeypatch,
    force_field,
    expected_filename,
):
    package_root = (
        make_fake_internal_ff_tree(
            tmp_path
        )
    )

    patch_internal_package_root(
        monkeypatch,
        package_root,
    )

    lunar_root = (
        make_fake_lunar(
            tmp_path
        )
    )

    expected = (
        lunar_root
        / "frc_files"
        / expected_filename
    )

    with pytest.raises(
        ValueError,
        match="Force field file not found",
    ) as exc_info:
        get_force_field_file(
            force_field=force_field,
            lunar_location=lunar_root,
        )

    assert str(
        expected
    ) in str(
        exc_info.value
    )


def test_lunar_location_accepts_string_path(
    tmp_path,
    monkeypatch,
):
    package_root = (
        make_fake_internal_ff_tree(
            tmp_path
        )
    )

    patch_internal_package_root(
        monkeypatch,
        package_root,
    )

    lunar_root = (
        make_fake_lunar(
            tmp_path,
            files=[
                "cvff.frc"
            ],
        )
    )

    result = get_force_field_file(
        force_field="CVFF",
        lunar_location=str(
            lunar_root
        ),
    )

    assert result == (
        lunar_root
        / "frc_files"
        / "cvff.frc"
    )


# =============================================================================
# Exact current force-field naming contract
# =============================================================================


def test_compass_lowercase_is_supported(
    tmp_path,
    monkeypatch,
):
    package_root = (
        make_fake_internal_ff_tree(
            tmp_path
        )
    )

    patch_internal_package_root(
        monkeypatch,
        package_root,
    )

    lunar_root = (
        make_fake_lunar(
            tmp_path,
            files=[
                "compass_published.frc"
            ],
        )
    )

    result = get_force_field_file(
        "compass",
        lunar_location=lunar_root,
    )

    assert result.name == (
        "compass_published.frc"
    )


@pytest.mark.parametrize(
    "force_field",
    [
        "COMPASS",
        "Compass",
        "pcff",
        "Pcff",
        "cvff",
        "dreiding",
    ],
)
def test_force_field_names_are_currently_case_sensitive(
    force_field,
):
    """
    Characterizes the current locator contract.

    Do not silently change this test into case-insensitive behavior unless
    force-field normalization is deliberately moved into this function.
    """
    with pytest.raises(
        ValueError,
        match="Unsupported force field requested",
    ):
        get_force_field_file(
            force_field
        )


def test_leading_or_trailing_whitespace_is_not_normalized():
    with pytest.raises(
        ValueError,
        match="Unsupported force field requested",
    ):
        get_force_field_file(
            " PCFF "
        )


# =============================================================================
# Foyer placeholders
# =============================================================================


@pytest.mark.parametrize(
    "force_field",
    [
        "OPLSAA",
        "GAFF",
    ],
)
def test_foyer_force_fields_are_explicitly_not_implemented(
    force_field,
):
    with pytest.raises(
        NotImplementedError,
        match="support is currently in development",
    ) as exc_info:
        get_force_field_file(
            force_field
        )

    assert force_field in str(
        exc_info.value
    )


def test_foyer_placeholder_does_not_need_lunar_location():
    with pytest.raises(
        NotImplementedError
    ):
        get_force_field_file(
            force_field="OPLSAA",
            lunar_location=None,
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
        "",
        "unknown",
    ],
)
def test_unsupported_force_field_raises_value_error(
    force_field,
):
    with pytest.raises(
        ValueError,
        match="Unsupported force field requested",
    ) as exc_info:
        get_force_field_file(
            force_field
        )

    assert repr(
        force_field
    ).strip("'") in str(
        exc_info.value
    ) or force_field == ""


def test_unsupported_error_includes_requested_force_field():
    with pytest.raises(
        ValueError
    ) as exc_info:
        get_force_field_file(
            "NOT_A_FORCE_FIELD"
        )

    assert (
        "NOT_A_FORCE_FIELD"
        in str(
            exc_info.value
        )
    )


# =============================================================================
# Path/layout contract
# =============================================================================


def test_external_force_fields_are_expected_under_frc_files(
    tmp_path,
    monkeypatch,
):
    package_root = (
        make_fake_internal_ff_tree(
            tmp_path
        )
    )

    patch_internal_package_root(
        monkeypatch,
        package_root,
    )

    lunar_root = (
        tmp_path / "LUNAR"
    )

    lunar_root.mkdir()

    # Correct filename, deliberately placed in the WRONG directory.
    (
        lunar_root
        / "cvff.frc"
    ).write_text(
        "wrong location",
        encoding="utf-8",
    )

    with pytest.raises(
        ValueError,
        match="Force field file not found",
    ):
        get_force_field_file(
            force_field="CVFF",
            lunar_location=lunar_root,
        )


def test_directory_named_like_force_field_file_is_rejected(
    tmp_path,
    monkeypatch,
):
    package_root = (
        make_fake_internal_ff_tree(
            tmp_path
        )
    )

    patch_internal_package_root(
        monkeypatch,
        package_root,
    )

    lunar_root = (
        make_fake_lunar(
            tmp_path
        )
    )

    fake_file_as_directory = (
        lunar_root
        / "frc_files"
        / "cvff.frc"
    )

    fake_file_as_directory.mkdir()

    with pytest.raises(
        ValueError,
        match="Force field file not found",
    ):
        get_force_field_file(
            force_field="CVFF",
            lunar_location=lunar_root,
        )