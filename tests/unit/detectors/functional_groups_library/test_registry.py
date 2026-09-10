import pytest

from AutoREACTER.detectors.functional_groups_library import registry


def test_load_functional_groups_returns_dict():
    """Registry should return a flat dictionary."""
    groups = registry.load_functional_groups()

    assert isinstance(groups, dict)


def test_load_functional_groups_is_not_empty():
    """The active functional-group registry should contain entries."""
    groups = registry.load_functional_groups()

    assert groups


def test_load_functional_groups_merges_modules(monkeypatch):
    """Entries from multiple modules should be merged into one dictionary."""
    module_1 = {
        "group_a": {
            "functionality_type": "mono",
            "smarts_1": "[O]",
            "group_name": "group_a_name",
            "comments": None,
        }
    }

    module_2 = {
        "group_b": {
            "functionality_type": "mono",
            "smarts_1": "[N]",
            "group_name": "group_b_name",
            "comments": None,
        }
    }

    monkeypatch.setattr(
        registry,
        "_FUNCTIONAL_GROUP_MODULES",
        [module_1, module_2],
    )

    groups = registry.load_functional_groups()

    assert groups == {
        "group_a": module_1["group_a"],
        "group_b": module_2["group_b"],
    }


def test_duplicate_entry_key_raises_value_error(monkeypatch):
    """Duplicate dictionary keys across modules must be rejected."""
    module_1 = {
        "duplicate_key": {
            "functionality_type": "mono",
            "smarts_1": "[O]",
            "group_name": "oxygen_group",
            "comments": None,
        }
    }

    module_2 = {
        "duplicate_key": {
            "functionality_type": "mono",
            "smarts_1": "[N]",
            "group_name": "nitrogen_group",
            "comments": None,
        }
    }

    monkeypatch.setattr(
        registry,
        "_FUNCTIONAL_GROUP_MODULES",
        [module_1, module_2],
    )

    with pytest.raises(
        ValueError,
        match="Duplicate functional-group key: duplicate_key",
    ):
        registry.load_functional_groups()


def test_duplicate_group_name_raises_value_error(monkeypatch):
    """Different keys must not share the same group_name."""
    module_1 = {
        "key_a": {
            "functionality_type": "mono",
            "smarts_1": "[O]",
            "group_name": "same_group",
            "comments": None,
        }
    }

    module_2 = {
        "key_b": {
            "functionality_type": "mono",
            "smarts_1": "[N]",
            "group_name": "same_group",
            "comments": None,
        }
    }

    monkeypatch.setattr(
        registry,
        "_FUNCTIONAL_GROUP_MODULES",
        [module_1, module_2],
    )

    with pytest.raises(
        ValueError,
        match="Duplicate group_name",
    ):
        registry.load_functional_groups()


def test_duplicate_group_name_error_identifies_both_keys(monkeypatch):
    """Duplicate-name error should identify both conflicting entries."""
    module_1 = {
        "first_key": {
            "functionality_type": "mono",
            "smarts_1": "[O]",
            "group_name": "duplicate_name",
            "comments": None,
        }
    }

    module_2 = {
        "second_key": {
            "functionality_type": "mono",
            "smarts_1": "[N]",
            "group_name": "duplicate_name",
            "comments": None,
        }
    }

    monkeypatch.setattr(
        registry,
        "_FUNCTIONAL_GROUP_MODULES",
        [module_1, module_2],
    )

    with pytest.raises(ValueError) as exc_info:
        registry.load_functional_groups()

    message = str(exc_info.value)

    assert "duplicate_name" in message
    assert "first_key" in message
    assert "second_key" in message


def test_empty_module_list_returns_empty_dict(monkeypatch):
    """An empty module collection should produce an empty registry."""
    monkeypatch.setattr(
        registry,
        "_FUNCTIONAL_GROUP_MODULES",
        [],
    )

    assert registry.load_functional_groups() == {}


def test_module_level_functional_groups_matches_loader():
    """The exported FUNCTIONAL_GROUPS registry should match a fresh load."""
    assert registry.FUNCTIONAL_GROUPS == registry.load_functional_groups()


def test_functional_groups_library_exposes_monomer_types():
    """Backward-compatible class should expose the merged registry."""
    library = registry.FunctionalGroupsLibrary()

    assert hasattr(library, "monomer_types")
    assert isinstance(library.monomer_types, dict)
    assert library.monomer_types == registry.load_functional_groups()


def test_real_registry_has_unique_group_names():
    """Active production registry should contain no duplicate group_name values."""
    groups = registry.load_functional_groups()

    group_names = [
        entry["group_name"]
        for entry in groups.values()
    ]

    assert len(group_names) == len(set(group_names))