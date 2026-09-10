from AutoREACTER.reaction_preparation.ff_wrapper.lunar_client.config import (
    LUNAR_ROOT_DIR,
)


def test_lunar_root_dir_is_string():
    assert isinstance(
        LUNAR_ROOT_DIR,
        str,
    )


def test_lunar_root_dir_is_not_empty():
    assert LUNAR_ROOT_DIR.strip()