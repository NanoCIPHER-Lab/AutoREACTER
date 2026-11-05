from pathlib import Path

from AutoREACTER.reaction_preparation.ff_wrapper.lunar_client.config import (
    LUNAR_ROOT_DIR,
)


def test_lunar_root_dir_is_path_string_path_or_none():
    assert (
        LUNAR_ROOT_DIR is None
        or isinstance(LUNAR_ROOT_DIR, (str, Path))
    )


def test_lunar_root_dir_is_not_empty_when_configured():
    if LUNAR_ROOT_DIR is None:
        return

    assert str(LUNAR_ROOT_DIR).strip()
