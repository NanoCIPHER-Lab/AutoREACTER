import pytest
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from AutoREACTER.detectors.reactions_library import registry


# =============================================================================
# Helpers
# =============================================================================

def make_valid_reaction(
    *,
    smarts="[C:1].[O:2]>>[C:1]-[O:2]",
    **extra,
):
    """Create a minimal valid AutoREACTER reaction-library entry."""
    reaction = {
        "reaction": smarts,
    }

    reaction.update(extra)

    return reaction


def templates_from_smarts(smarts):
    """
    Convert reaction SMARTS into independent reactant/product Mol copies.

    RDKit's GetReactantTemplate() and GetProductTemplate() return objects
    owned by the ChemicalReaction. Therefore, the returned templates must
    be copied before the ChemicalReaction goes out of scope.

    Without the Chem.Mol copies below, tests may retain dangling C++ object
    references and can eventually cause a segmentation fault.
    """
    reaction = rdChemReactions.ReactionFromSmarts(smarts)

    assert reaction is not None

    reactants = [
        Chem.Mol(
            reaction.GetReactantTemplate(i)
        )
        for i in range(
            reaction.GetNumReactantTemplates()
        )
    ]

    products = [
        Chem.Mol(
            reaction.GetProductTemplate(i)
        )
        for i in range(
            reaction.GetNumProductTemplates()
        )
    ]

    return reactants, products


# =============================================================================
# ReactionLibraryValidationError
# =============================================================================

def test_validation_error_is_value_error():
    assert issubclass(
        registry.ReactionLibraryValidationError,
        ValueError,
    )


# =============================================================================
# _atom_maps_in_templates()
# =============================================================================

def test_atom_maps_in_templates_collects_all_nonzero_maps():
    reactants, _ = templates_from_smarts(
        "[C:1].[O:2].[N:7]>>[C:1]-[O:2]-[N:7]"
    )

    result = registry._atom_maps_in_templates(
        reactants
    )

    assert result == {
        1,
        2,
        7,
    }


def test_atom_maps_in_templates_ignores_unmapped_atoms():
    reactants, _ = templates_from_smarts(
        "[C:1].O>>[C:1]-O"
    )

    result = registry._atom_maps_in_templates(
        reactants
    )

    assert result == {
        1,
    }


def test_atom_maps_in_templates_empty_templates_returns_empty_set():
    result = registry._atom_maps_in_templates(
        []
    )

    assert result == set()


def test_atom_maps_in_templates_deduplicates_map_numbers():
    reactants, _ = templates_from_smarts(
        "[C:1].[O:1]>>[C:1]-[O:1]"
    )

    result = registry._atom_maps_in_templates(
        reactants
    )

    assert result == {
        1,
    }


# =============================================================================
# _has_bond_between_atom_maps()
# =============================================================================

def test_has_bond_between_atom_maps_returns_true_when_bond_exists():
    _, products = templates_from_smarts(
        "[C:1].[O:2]>>[C:1]-[O:2]"
    )

    result = registry._has_bond_between_atom_maps(
        products,
        1,
        2,
    )

    assert result is True


def test_has_bond_between_atom_maps_is_order_independent():
    _, products = templates_from_smarts(
        "[C:1].[O:2]>>[C:1]-[O:2]"
    )

    result = registry._has_bond_between_atom_maps(
        products,
        2,
        1,
    )

    assert result is True


def test_has_bond_between_atom_maps_returns_false_without_bond():
    _, products = templates_from_smarts(
        "[C:1].[O:2]>>[C:1].[O:2]"
    )

    result = registry._has_bond_between_atom_maps(
        products,
        1,
        2,
    )

    assert result is False


def test_has_bond_between_atom_maps_returns_false_for_missing_map():
    _, products = templates_from_smarts(
        "[C:1].[O:2]>>[C:1]-[O:2]"
    )

    result = registry._has_bond_between_atom_maps(
        products,
        1,
        99,
    )

    assert result is False


def test_has_bond_between_atom_maps_empty_templates_returns_false():
    result = registry._has_bond_between_atom_maps(
        [],
        1,
        2,
    )

    assert result is False


# =============================================================================
# _validate_reaction_smarts()
# =============================================================================

def test_valid_reaction_smarts_has_no_errors():
    reaction = make_valid_reaction()

    errors = registry._validate_reaction_smarts(
        "test_reaction",
        reaction,
    )

    assert errors == []


def test_missing_reaction_key_is_reported():
    errors = registry._validate_reaction_smarts(
        "test_reaction",
        {},
    )

    assert errors == [
        "test_reaction: missing required key 'reaction'"
    ]


def test_empty_reaction_smarts_is_reported_as_missing():
    errors = registry._validate_reaction_smarts(
        "test_reaction",
        {
            "reaction": "",
        },
    )

    assert errors == [
        "test_reaction: missing required key 'reaction'"
    ]


def test_validation_can_be_explicitly_disabled():
    """
    Special reactions can intentionally bypass the initiator-bond check.
    """
    reaction = {
        "reaction": "THIS DOES NOT NEED TO BE VALIDATED",
        "validate_initiator_bond": False,
    }

    errors = registry._validate_reaction_smarts(
        "special_reaction",
        reaction,
    )

    assert errors == []


def test_missing_rdkit_is_reported(monkeypatch):
    monkeypatch.setattr(
        registry,
        "rdChemReactions",
        None,
    )

    errors = registry._validate_reaction_smarts(
        "test_reaction",
        make_valid_reaction(),
    )

    assert errors == [
        "test_reaction: RDKit is required to validate reaction SMARTS"
    ]


@pytest.mark.parametrize(
    "initiator_maps",
    [
        (),
        (1,),
        (1, 2, 3),
        [1],
        [1, 2, 3],
    ],
)
def test_initiator_atom_maps_must_contain_exactly_two_maps(
    initiator_maps,
):
    reaction = make_valid_reaction(
        initiator_atom_maps=initiator_maps,
    )

    errors = registry._validate_reaction_smarts(
        "test_reaction",
        reaction,
    )

    assert errors == [
        "test_reaction: initiator_atom_maps must contain exactly two atom maps"
    ]


def test_invalid_reaction_smarts_is_reported():
    reaction = {
        "reaction": "not_a_reaction_smarts",
    }

    errors = registry._validate_reaction_smarts(
        "bad_reaction",
        reaction,
    )

    assert errors

    assert errors[0].startswith(
        "bad_reaction:"
    )

    assert (
        "invalid reaction SMARTS"
        in errors[0]
        or
        "could not parse reaction SMARTS"
        in errors[0]
    )


def test_missing_initiator_map_from_reactants_is_reported():
    reaction = {
        "reaction": "[C:1].O>>[C:1]-[O:2]",
    }

    errors = registry._validate_reaction_smarts(
        "test_reaction",
        reaction,
    )

    assert any(
        "initiator atom maps missing from reactants: [2]"
        in error
        for error in errors
    )


def test_missing_initiator_map_from_products_is_reported():
    reaction = {
        "reaction": "[C:1].[O:2]>>[C:1]-O",
    }

    errors = registry._validate_reaction_smarts(
        "test_reaction",
        reaction,
    )

    assert any(
        "initiator atom maps missing from products: [2]"
        in error
        for error in errors
    )


def test_missing_initiator_bond_is_reported():
    reaction = {
        "reaction": "[C:1].[O:2]>>[C:1].[O:2]",
    }

    errors = registry._validate_reaction_smarts(
        "test_reaction",
        reaction,
    )

    assert (
        "test_reaction: product does not contain required "
        "AutoREACTER initiator bond between atom maps 1 and 2"
        in errors
    )


def test_missing_product_map_can_generate_multiple_validation_errors():
    """
    If a mapped initiator disappears from the product, both the missing-map
    error and missing-bond error should be preserved.
    """
    reaction = {
        "reaction": "[C:1].[O:2]>>[C:1]-O",
    }

    errors = registry._validate_reaction_smarts(
        "test_reaction",
        reaction,
    )

    assert any(
        "missing from products"
        in error
        for error in errors
    )

    assert any(
        "product does not contain required"
        in error
        for error in errors
    )


def test_custom_initiator_atom_maps_are_supported():
    reaction = {
        "reaction": "[C:7].[O:9]>>[C:7]-[O:9]",
        "initiator_atom_maps": (
            7,
            9,
        ),
    }

    errors = registry._validate_reaction_smarts(
        "custom_maps",
        reaction,
    )

    assert errors == []


def test_custom_maps_still_require_product_bond():
    reaction = {
        "reaction": "[C:7].[O:9]>>[C:7].[O:9]",
        "initiator_atom_maps": (
            7,
            9,
        ),
    }

    errors = registry._validate_reaction_smarts(
        "custom_maps",
        reaction,
    )

    assert any(
        "between atom maps 7 and 9"
        in error
        for error in errors
    )


def test_initiator_atom_maps_are_converted_to_int():
    """
    String map numbers are accepted because registry converts them using int().
    """
    reaction = {
        "reaction": "[C:7].[O:9]>>[C:7]-[O:9]",
        "initiator_atom_maps": (
            "7",
            "9",
        ),
    }

    errors = registry._validate_reaction_smarts(
        "custom_maps",
        reaction,
    )

    assert errors == []


# =============================================================================
# validate_reactions()
# =============================================================================

def test_validate_reactions_accepts_valid_dictionary():
    reactions = {
        "reaction_1": make_valid_reaction(),
        "reaction_2": make_valid_reaction(
            smarts="[N:1].[C:2]>>[N:1]-[C:2]"
        ),
    }

    result = registry.validate_reactions(
        reactions
    )

    assert result is None


def test_validate_reactions_rejects_non_dictionary_entry():
    reactions = {
        "bad_reaction": "not a dictionary",
    }

    with pytest.raises(
        registry.ReactionLibraryValidationError,
        match="reaction entry must be a dictionary",
    ):
        registry.validate_reactions(
            reactions
        )


def test_validate_reactions_rejects_missing_reaction_smarts():
    reactions = {
        "bad_reaction": {},
    }

    with pytest.raises(
        registry.ReactionLibraryValidationError,
        match="missing required key 'reaction'",
    ):
        registry.validate_reactions(
            reactions
        )


def test_validate_reactions_aggregates_multiple_reaction_errors():
    reactions = {
        "bad_1": {},
        "bad_2": {
            "reaction": "[C:1].[O:2]>>[C:1].[O:2]"
        },
    }

    with pytest.raises(
        registry.ReactionLibraryValidationError
    ) as exc_info:
        registry.validate_reactions(
            reactions
        )

    message = str(
        exc_info.value
    )

    assert "bad_1" in message
    assert "bad_2" in message
    assert "missing required key 'reaction'" in message
    assert "product does not contain required" in message


def test_validation_error_message_has_expected_header():
    reactions = {
        "bad": {},
    }

    with pytest.raises(
        registry.ReactionLibraryValidationError
    ) as exc_info:
        registry.validate_reactions(
            reactions
        )

    assert str(
        exc_info.value
    ).startswith(
        "Reaction library validation failed:"
    )


def test_empty_reaction_dictionary_is_valid():
    result = registry.validate_reactions(
        {}
    )

    assert result is None


# =============================================================================
# load_reactions()
# =============================================================================

def test_load_reactions_merges_modules(monkeypatch):
    module_1 = {
        "reaction_a": make_valid_reaction(),
    }

    module_2 = {
        "reaction_b": make_valid_reaction(
            smarts="[N:1].[C:2]>>[N:1]-[C:2]"
        ),
    }

    monkeypatch.setattr(
        registry,
        "_REACTION_MODULES",
        [
            module_1,
            module_2,
        ],
    )

    result = registry.load_reactions()

    assert set(
        result
    ) == {
        "reaction_a",
        "reaction_b",
    }

    assert (
        result["reaction_a"]
        is module_1["reaction_a"]
    )

    assert (
        result["reaction_b"]
        is module_2["reaction_b"]
    )


def test_load_reactions_rejects_duplicate_reaction_names(
    monkeypatch,
):
    module_1 = {
        "duplicate": make_valid_reaction(),
    }

    module_2 = {
        "duplicate": make_valid_reaction(),
    }

    monkeypatch.setattr(
        registry,
        "_REACTION_MODULES",
        [
            module_1,
            module_2,
        ],
    )

    with pytest.raises(
        ValueError,
        match="Duplicate reaction name: duplicate",
    ):
        registry.load_reactions()


def test_load_reactions_empty_modules_returns_empty_dict(
    monkeypatch,
):
    monkeypatch.setattr(
        registry,
        "_REACTION_MODULES",
        [],
    )

    result = registry.load_reactions()

    assert result == {}


def test_load_reactions_calls_validation(monkeypatch):
    module = {
        "reaction_a": {
            "reaction": "anything"
        }
    }

    monkeypatch.setattr(
        registry,
        "_REACTION_MODULES",
        [
            module,
        ],
    )

    captured = {}

    def fake_validate(reactions):
        captured["reactions"] = reactions

    monkeypatch.setattr(
        registry,
        "validate_reactions",
        fake_validate,
    )

    result = registry.load_reactions()

    assert (
        captured["reactions"]
        == result
    )

    assert "reaction_a" in result


def test_load_reactions_propagates_validation_failure(
    monkeypatch,
):
    module = {
        "bad": {},
    }

    monkeypatch.setattr(
        registry,
        "_REACTION_MODULES",
        [
            module,
        ],
    )

    with pytest.raises(
        registry.ReactionLibraryValidationError
    ):
        registry.load_reactions()


# =============================================================================
# ReactionLibrary backward compatibility
# =============================================================================

def test_reaction_library_exposes_reactions(
    monkeypatch,
):
    expected = {
        "test": make_valid_reaction(),
    }

    monkeypatch.setattr(
        registry,
        "load_reactions",
        lambda: expected,
    )

    library = registry.ReactionLibrary()

    assert library.reactions is expected


# =============================================================================
# Production registry sanity checks
# =============================================================================

def test_module_level_reactions_is_dictionary():
    assert isinstance(
        registry.REACTIONS,
        dict,
    )


def test_module_level_reactions_is_not_empty():
    assert registry.REACTIONS


def test_production_reaction_names_are_unique():
    names = list(
        registry.REACTIONS.keys()
    )

    assert (
        len(names)
        == len(set(names))
    )


def test_production_registry_passes_validation():
    """
    All currently enabled production reactions must satisfy the
    registry validation rules.
    """
    result = registry.validate_reactions(
        registry.REACTIONS
    )

    assert result is None


def test_fresh_production_load_matches_exported_registry():
    freshly_loaded = registry.load_reactions()

    assert (
        freshly_loaded
        == registry.REACTIONS
    )