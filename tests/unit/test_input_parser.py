import math
from types import SimpleNamespace

import pytest
from rdkit import Chem

from AutoREACTER.input_parser import (
    CompatibilityError,
    DuplicateMonomerError,
    InputConflictError,
    InputError,
    InputParser,
    InputSchemaError,
    MonomerEntry,
    NumericFieldError,
    Simulation,
    SimulationSetup,
    SmilesValidationError,
)


# =============================================================================
# Fixtures / helpers
# =============================================================================


@pytest.fixture
def parser():
    return InputParser()


def make_counts_input():
    return {
        "simulation_name": "counts_demo",
        "simulations": [
            {
                "tag": "small",
                "temperature": 300,
                "density": 0.8,
                "monomer_counts": {
                    "ethanol": 10,
                    "ethylamine": 20,
                },
            },
            {
                "tag": "large",
                "temperature": 400,
                "density": 1.0,
                "monomer_counts": {
                    "ethanol": 100,
                    "ethylamine": 200,
                },
            },
        ],
        "monomers": [
            {
                "name": "ethanol",
                "smiles": "CCO",
            },
            {
                "name": "ethylamine",
                "smiles": "CCN",
            },
        ],
    }


def make_ratio_input():
    return {
        "simulation_name": "ratio_demo",
        "simulations": [
            {
                "tag": "small",
                "temperature": 300,
                "density": 0.8,
                "total_atoms": 10000,
                "monomer_ratios": {
                    "ethanol": 1.0,
                    "ethylamine": 2.0,
                },
            },
            {
                "tag": "large",
                "temperature": 400,
                "density": 1.0,
                "total_atoms": 100000,
                "monomer_ratios": {
                    "ethanol": 1.0,
                    "ethylamine": 2.0,
                },
            },
        ],
        "monomers": [
            {
                "name": "ethanol",
                "smiles": "CCO",
            },
            {
                "name": "ethylamine",
                "smiles": "CCN",
            },
        ],
    }


# =============================================================================
# Exception hierarchy
# =============================================================================


@pytest.mark.parametrize(
    "exception_class",
    [
        InputSchemaError,
        InputConflictError,
        NumericFieldError,
        SmilesValidationError,
        DuplicateMonomerError,
        CompatibilityError,
    ],
)
def test_input_exceptions_inherit_from_input_error(
    exception_class,
):
    assert issubclass(
        exception_class,
        InputError,
    )


# =============================================================================
# Dataclasses
# =============================================================================


def test_monomer_entry_defaults():
    monomer = MonomerEntry(
        id=1,
        data_id="data_1",
        name="ethanol",
        smiles="CCO",
        count={"sim": 10},
        ratio=None,
    )

    assert monomer.rdkit_mol is None
    assert monomer.molecule_3Dmol_path is None
    assert monomer.lmp_molecule_file is None
    assert monomer.num_atoms is None
    assert monomer.molecular_weight is None
    assert monomer.status is True


def test_monomer_entry_uses_slots():
    monomer = MonomerEntry(
        id=1,
        data_id="data_1",
        name="ethanol",
        smiles="CCO",
        count={"sim": 10},
        ratio=None,
    )

    with pytest.raises(AttributeError):
        monomer.random_attribute = 10


def test_simulation_defaults():
    simulation = Simulation(
        tag="sim",
        temperature=300.0,
        density=1.0,
    )

    assert simulation.monomer_counts is None
    assert simulation.monomer_ratios is None
    assert simulation.total_atoms is None
    assert simulation.initial_box_volume is None
    assert simulation.initial_box_length is None


def test_simulation_setup_workflow_defaults():
    setup = SimulationSetup(
        simulation_name="demo",
        temperature=[300.0],
        density=[1.0],
        force_field="PCFF",
        monomers=[],
    )

    assert setup.deep_search is True
    assert setup.loop is True
    assert setup.reaction_iteration_depth == 5
    assert setup.wildcards is True

    assert (
        setup.deduplicate_reaction_templates
        is True
    )

    assert (
        setup.write_second_reaction_stage
        is False
    )


# =============================================================================
# validate_basic_format
# =============================================================================


def test_validate_basic_format_accepts_required_keys(
    parser,
):
    parser.validate_basic_format(
        {
            "simulation_name": "demo",
            "simulations": [],
            "monomers": [],
        }
    )


@pytest.mark.parametrize(
    "value",
    [
        None,
        [],
        "string",
        123,
    ],
)
def test_validate_basic_format_requires_dictionary(
    parser,
    value,
):
    with pytest.raises(
        InputSchemaError,
        match="Expected input to be a dictionary",
    ):
        parser.validate_basic_format(value)


@pytest.mark.parametrize(
    "missing_key",
    [
        "simulation_name",
        "simulations",
        "monomers",
    ],
)
def test_validate_basic_format_requires_each_core_key(
    parser,
    missing_key,
):
    data = {
        "simulation_name": "demo",
        "simulations": [],
        "monomers": [],
    }

    del data[missing_key]

    with pytest.raises(
        InputSchemaError,
        match="Missing required key",
    ):
        parser.validate_basic_format(data)


# =============================================================================
# Composition-mode detection
# =============================================================================


def test_get_inputs_mode_detects_counts(parser):
    simulations = [
        {
            "monomer_counts": {
                "a": 1,
            }
        }
    ]

    assert (
        parser._get_inputs_mode(simulations)
        == "counts"
    )


def test_get_inputs_mode_detects_ratio(parser):
    simulations = [
        {
            "monomer_ratios": {
                "a": 1.0,
            }
        }
    ]

    assert (
        parser._get_inputs_mode(simulations)
        == "ratio"
    )


@pytest.mark.parametrize(
    "simulations",
    [
        [],
        None,
        {},
        "bad",
    ],
)
def test_get_inputs_mode_requires_nonempty_list(
    parser,
    simulations,
):
    with pytest.raises(InputSchemaError):
        parser._get_inputs_mode(simulations)


def test_get_inputs_mode_requires_dict_entries(
    parser,
):
    with pytest.raises(InputSchemaError):
        parser._get_inputs_mode(
            ["not a dictionary"]
        )


def test_get_inputs_mode_rejects_counts_and_ratios_together(
    parser,
):
    with pytest.raises(
        InputConflictError,
        match="not both",
    ):
        parser._get_inputs_mode(
            [
                {
                    "monomer_counts": {
                        "a": 1,
                    },
                    "monomer_ratios": {
                        "a": 1.0,
                    },
                }
            ]
        )


def test_get_inputs_mode_rejects_missing_composition(
    parser,
):
    with pytest.raises(InputSchemaError):
        parser._get_inputs_mode(
            [
                {
                    "tag": "sim",
                }
            ]
        )


def test_get_inputs_mode_rejects_mixed_modes(
    parser,
):
    with pytest.raises(
        InputConflictError,
        match="same composition method",
    ):
        parser._get_inputs_mode(
            [
                {
                    "monomer_counts": {
                        "a": 1,
                    }
                },
                {
                    "monomer_ratios": {
                        "a": 1.0,
                    }
                },
            ]
        )


# =============================================================================
# Temperature / density validation
# =============================================================================


@pytest.mark.parametrize(
    "value, expected",
    [
        (1, 1.0),
        (300, 300.0),
        (298.15, 298.15),
    ],
)
def test_validate_temperature_accepts_positive_numbers(
    parser,
    value,
    expected,
):
    assert (
        parser._validate_temperature(value)
        == expected
    )


@pytest.mark.parametrize(
    "value",
    [
        0,
        -1,
        -300.0,
        True,
        False,
        "300",
        None,
    ],
)
def test_validate_temperature_rejects_invalid_values(
    parser,
    value,
):
    with pytest.raises(NumericFieldError):
        parser._validate_temperature(value)


@pytest.mark.parametrize(
    "value, expected",
    [
        (1, 1.0),
        (0.8, 0.8),
        (1.25, 1.25),
    ],
)
def test_validate_density_accepts_positive_numbers(
    parser,
    value,
    expected,
):
    assert parser._validate_density(value) == expected


@pytest.mark.parametrize(
    "value",
    [
        0,
        -1,
        True,
        False,
        "0.8",
        None,
    ],
)
def test_validate_density_rejects_invalid_values(
    parser,
    value,
):
    with pytest.raises(NumericFieldError):
        parser._validate_density(value)


# =============================================================================
# Force-field validation
# =============================================================================


def test_force_field_defaults_to_pcff(parser):
    assert (
        parser._validate_force_field(None)
        == "PCFF"
    )


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("PCFF", "PCFF"),
        ("pcff", "PCFF"),
        ("PCFF-IFF", "PCFF-IFF"),
        ("pcff-iff", "PCFF-IFF"),
        ("CVFF", "CVFF"),
        ("cvff", "CVFF"),
        ("CVFF-IFF", "CVFF-IFF"),
        ("clayff", "Clay-FF"),
        ("Clay-FF", "Clay-FF"),
        ("dreiding", "DREIDING"),
        ("drieding", "DREIDING"),
    ],
)
def test_force_field_aliases_are_normalized(
    parser,
    raw,
    expected,
):
    assert (
        parser._validate_force_field(raw)
        == expected
    )


def test_compass_normalizes_to_declared_canonical_name(
    parser,
):
    """
    ForceFieldType declares the canonical spelling as 'Compass'.

    This test intentionally checks that the runtime normalizer agrees
    with that public contract.
    """
    assert (
        parser._validate_force_field("compass")
        == "Compass"
    )


@pytest.mark.parametrize(
    "raw",
    [
        "",
        "   ",
        123,
        [],
        {},
    ],
)
def test_force_field_rejects_invalid_schema(
    parser,
    raw,
):
    with pytest.raises(InputSchemaError):
        parser._validate_force_field(raw)


def test_force_field_rejects_unknown_name(parser):
    with pytest.raises(
        InputSchemaError,
        match="Unsupported force field",
    ):
        parser._validate_force_field(
            "not-a-force-field"
        )


@pytest.mark.parametrize(
    "raw",
    [
        "OPLSAA",
        "opls",
        "opls-aa",
        "GAFF",
        "gaff",
    ],
)
def test_current_lunar_workflow_rejects_incompatible_force_fields(
    parser,
    raw,
):
    with pytest.raises(CompatibilityError):
        parser._validate_force_field(raw)


# =============================================================================
# Workflow option normalization
# =============================================================================


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("deep_search", "deepsearch"),
        ("deep-search", "deepsearch"),
        ("Deep Search", "deepsearch"),
        (" DEEP_SEARCH ", "deepsearch"),
        ("write_second_reaction_stage", "writesecondreactionstage"),
    ],
)
def test_normalized_option_key(
    parser,
    raw,
    expected,
):
    assert (
        parser._normalized_option_key(raw)
        == expected
    )


def test_get_workflow_option_returns_default(
    parser,
):
    assert (
        parser._get_workflow_option(
            inputs={},
            key="wildcards",
            default=True,
        )
        is True
    )


def test_get_workflow_option_accepts_alias(
    parser,
):
    value = parser._get_workflow_option(
        inputs={
            "use-wildcards": False,
        },
        key="wildcards",
        default=True,
        aliases=[
            "use_wildcards",
        ],
    )

    assert value is False


def test_get_workflow_option_accepts_equivalent_duplicate_alias_values(
    parser,
):
    value = parser._get_workflow_option(
        inputs={
            "wildcards": False,
            "use-wildcards": False,
        },
        key="wildcards",
        default=True,
        aliases=[
            "use_wildcards",
        ],
    )

    assert value is False


def test_get_workflow_option_rejects_conflicting_alias_values(
    parser,
):
    with pytest.raises(
        InputConflictError,
        match="Conflicting values",
    ):
        parser._get_workflow_option(
            inputs={
                "wildcards": True,
                "use-wildcards": False,
            },
            key="wildcards",
            default=True,
            aliases=[
                "use_wildcards",
            ],
        )


# =============================================================================
# Boolean workflow options
# =============================================================================


@pytest.mark.parametrize(
    "value",
    [
        True,
        "true",
        "TRUE",
        "yes",
        "Y",
        "on",
        "1",
    ],
)
def test_validate_bool_option_true_values(
    parser,
    value,
):
    result = parser._validate_bool_option(
        inputs={
            "feature": value,
        },
        key="feature",
        default=False,
    )

    assert result is True


@pytest.mark.parametrize(
    "value",
    [
        False,
        "false",
        "FALSE",
        "no",
        "N",
        "off",
        "0",
    ],
)
def test_validate_bool_option_false_values(
    parser,
    value,
):
    result = parser._validate_bool_option(
        inputs={
            "feature": value,
        },
        key="feature",
        default=True,
    )

    assert result is False


@pytest.mark.parametrize(
    "value",
    [
        1,
        0,
        2,
        None,
        [],
        "maybe",
    ],
)
def test_validate_bool_option_rejects_invalid_values(
    parser,
    value,
):
    with pytest.raises(InputSchemaError):
        parser._validate_bool_option(
            inputs={
                "feature": value,
            },
            key="feature",
            default=True,
        )


# =============================================================================
# Reaction iteration depth
# =============================================================================


def test_reaction_iteration_depth_defaults_to_five(
    parser,
):
    assert (
        parser._validate_reaction_iteration_depth(
            {}
        )
        == 5
    )


@pytest.mark.parametrize(
    "value, expected",
    [
        (0, 0),
        (1, 1),
        (5, 5),
        (20, 20),
        (False, 0),
        (True, 5),
        ("0", 0),
        ("1", 1),
        ("10", 10),
        ("false", 0),
        ("no", 0),
        ("off", 0),
        ("none", 0),
        ("true", 5),
        ("yes", 5),
        ("on", 5),
    ],
)
def test_reaction_iteration_depth_accepts_supported_values(
    parser,
    value,
    expected,
):
    result = parser._validate_reaction_iteration_depth(
        {
            "reaction_iteration_depth": value,
        }
    )

    assert result == expected


@pytest.mark.parametrize(
    "alias",
    [
        "rxn_iteration_depth",
        "reaction_depth",
        "iteration_depth",
        "max_loop_count",
        "max_iterations",
        "iterations",
        "loop",
    ],
)
def test_reaction_iteration_depth_accepts_aliases(
    parser,
    alias,
):
    result = parser._validate_reaction_iteration_depth(
        {
            alias: 7,
        }
    )

    assert result == 7


@pytest.mark.parametrize(
    "value",
    [
        -1,
        -10,
        1.5,
        None,
        [],
        {},
        "hello",
        "2.5",
    ],
)
def test_reaction_iteration_depth_rejects_invalid_values(
    parser,
    value,
):
    with pytest.raises(InputSchemaError):
        parser._validate_reaction_iteration_depth(
            {
                "reaction_iteration_depth": value,
            }
        )


def test_reaction_iteration_depth_rejects_conflicting_aliases(
    parser,
):
    with pytest.raises(InputConflictError):
        parser._validate_reaction_iteration_depth(
            {
                "reaction_iteration_depth": 5,
                "loop": 3,
            }
        )


def test_validate_loop_zero_disables_loop(parser):
    assert (
        parser._validate_loop(
            {
                "reaction_iteration_depth": 0,
            }
        )
        == (False, None)
    )


def test_validate_loop_positive_depth_enables_loop(
    parser,
):
    assert (
        parser._validate_loop(
            {
                "reaction_iteration_depth": 8,
            }
        )
        == (True, 8)
    )


# =============================================================================
# SMILES validation
# =============================================================================


def test_validate_smiles_returns_canonical_smiles(
    parser,
):
    smiles, mol = parser._validate_smiles(
        "OCC"
    )

    assert smiles == "CCO"
    assert isinstance(mol, Chem.Mol)


def test_validate_smiles_strips_whitespace(
    parser,
):
    smiles, _ = parser._validate_smiles(
        "  CCO  "
    )

    assert smiles == "CCO"


@pytest.mark.parametrize(
    "value",
    [
        None,
        "",
        "   ",
        123,
        [],
    ],
)
def test_validate_smiles_requires_nonempty_string(
    parser,
    value,
):
    with pytest.raises(SmilesValidationError):
        parser._validate_smiles(value)


def test_validate_smiles_rejects_invalid_rdkit_smiles(
    parser,
):
    with pytest.raises(
        SmilesValidationError,
        match="Invalid SMILES",
    ):
        parser._validate_smiles(
            "C1CC"
        )


# =============================================================================
# Duplicate SMILES
# =============================================================================


def test_validate_no_duplicate_smiles_adds_new_smiles(
    parser,
):
    seen = []

    returned = parser.validate_no_duplicate_smiles(
        "CCO",
        seen,
    )

    assert returned is seen
    assert seen == ["CCO"]


def test_validate_no_duplicate_smiles_rejects_duplicate(
    parser,
):
    seen = ["CCO"]

    with pytest.raises(
        DuplicateMonomerError,
        match="Duplicate monomer",
    ):
        parser.validate_no_duplicate_smiles(
            "CCO",
            seen,
        )


def test_validate_inputs_detects_canonical_smiles_duplicates(
    parser,
):
    data = make_counts_input()

    data["monomers"] = [
        {
            "name": "one",
            "smiles": "CCO",
        },
        {
            "name": "two",
            "smiles": "OCC",
        },
    ]

    for simulation in data["simulations"]:
        simulation["monomer_counts"] = {
            "one": 1,
            "two": 1,
        }

    with pytest.raises(DuplicateMonomerError):
        parser.validate_inputs(data)


# =============================================================================
# Derived molecular properties
# =============================================================================


def test_derive_molecule_properties_includes_hydrogens(
    parser,
):
    mol = Chem.MolFromSmiles("C")

    num_atoms, molecular_weight = (
        parser._derive_molecule_properties(
            mol
        )
    )

    # methane = 1 C + 4 H
    assert num_atoms == 5

    assert molecular_weight == pytest.approx(
        16.043,
        rel=1e-3,
    )


def test_int_to_dict(parser):
    assert parser._int_to_dict(42) == {
        "_": 42
    }


# =============================================================================
# Simulation validation
# =============================================================================


def test_validate_simulations_counts_mode(
    parser,
):
    systems = make_counts_input()[
        "simulations"
    ]

    result = parser._validate_simulations(
        systems,
        "counts",
    )

    assert result["method"] == "counts"

    assert result["temperatures"] == [
        300.0,
        400.0,
    ]

    assert result["density"] == [
        0.8,
        1.0,
    ]

    assert len(result["simulations"]) == 2

    assert all(
        isinstance(simulation, Simulation)
        for simulation in result["simulations"]
    )


def test_validate_simulations_ratio_mode(
    parser,
):
    systems = make_ratio_input()[
        "simulations"
    ]

    result = parser._validate_simulations(
        systems,
        "ratio",
    )

    assert result["method"] == "ratio"

    assert (
        result["simulations"][0].total_atoms
        == 10000
    )

    assert (
        result["simulations"][1].total_atoms
        == 100000
    )


def test_validate_simulations_rejects_duplicate_tags(
    parser,
):
    systems = make_counts_input()[
        "simulations"
    ]

    systems[1]["tag"] = "small"

    with pytest.raises(
        InputSchemaError,
        match="Duplicate system tag",
    ):
        parser._validate_simulations(
            systems,
            "counts",
        )


@pytest.mark.parametrize(
    "tag",
    [
        "",
        "   ",
        None,
        123,
    ],
)
def test_validate_simulations_requires_valid_tag(
    parser,
    tag,
):
    systems = make_counts_input()[
        "simulations"
    ]

    systems[0]["tag"] = tag

    with pytest.raises(InputSchemaError):
        parser._validate_simulations(
            systems,
            "counts",
        )


def test_counts_mode_rejects_total_atoms(
    parser,
):
    systems = make_counts_input()[
        "simulations"
    ]

    systems[0]["total_atoms"] = 10000

    with pytest.raises(
        InputSchemaError,
        match="total_atoms",
    ):
        parser._validate_simulations(
            systems,
            "counts",
        )


@pytest.mark.parametrize(
    "value",
    [
        -1,
        1.5,
        True,
        "10",
    ],
)
def test_counts_mode_rejects_invalid_counts(
    parser,
    value,
):
    systems = make_counts_input()[
        "simulations"
    ]

    systems[0]["monomer_counts"][
        "ethanol"
    ] = value

    with pytest.raises(NumericFieldError):
        parser._validate_simulations(
            systems,
            "counts",
        )


def test_counts_mode_allows_zero_count(
    parser,
):
    systems = make_counts_input()[
        "simulations"
    ]

    systems[0]["monomer_counts"][
        "ethanol"
    ] = 0

    result = parser._validate_simulations(
        systems,
        "counts",
    )

    assert (
        result["simulations"][0]
        .monomer_counts["ethanol"]
        == 0
    )


@pytest.mark.parametrize(
    "value",
    [
        0,
        -1,
        True,
        1.5,
        "10000",
        None,
    ],
)
def test_ratio_mode_requires_positive_integer_total_atoms(
    parser,
    value,
):
    systems = make_ratio_input()[
        "simulations"
    ]

    systems[0]["total_atoms"] = value

    with pytest.raises(NumericFieldError):
        parser._validate_simulations(
            systems,
            "ratio",
        )


@pytest.mark.parametrize(
    "value",
    [
        -1,
        True,
        "1.0",
        None,
    ],
)
def test_ratio_mode_rejects_invalid_ratio_values(
    parser,
    value,
):
    systems = make_ratio_input()[
        "simulations"
    ]

    systems[0]["monomer_ratios"][
        "ethanol"
    ] = value

    with pytest.raises(NumericFieldError):
        parser._validate_simulations(
            systems,
            "ratio",
        )


def test_ratio_mode_allows_zero_ratio(
    parser,
):
    systems = make_ratio_input()[
        "simulations"
    ]

    for system in systems:
        system["monomer_ratios"][
            "ethanol"
        ] = 0

    result = parser._validate_simulations(
        systems,
        "ratio",
    )

    assert (
        result["simulations"][0]
        .monomer_ratios["ethanol"]
        == 0
    )


def test_ratio_mode_requires_identical_ratios_between_systems(
    parser,
):
    systems = make_ratio_input()[
        "simulations"
    ]

    systems[1]["monomer_ratios"][
        "ethanol"
    ] = 3.0

    with pytest.raises(
        InputSchemaError,
        match="identical",
    ):
        parser._validate_simulations(
            systems,
            "ratio",
        )


# =============================================================================
# System/monomer key consistency
# =============================================================================


def test_system_monomer_keys_accept_exact_match(
    parser,
):
    data = make_counts_input()

    parser._validate_system_monomer_keys(
        data,
        data["simulations"],
        "counts",
    )


def test_system_monomer_keys_reject_unknown_name(
    parser,
):
    data = make_counts_input()

    data["simulations"][0][
        "monomer_counts"
    ]["unknown"] = 10

    with pytest.raises(
        InputSchemaError,
        match="unknown monomer",
    ):
        parser._validate_system_monomer_keys(
            data,
            data["simulations"],
            "counts",
        )


def test_system_monomer_keys_reject_missing_name(
    parser,
):
    data = make_counts_input()

    del data["simulations"][0][
        "monomer_counts"
    ]["ethanol"]

    with pytest.raises(
        InputSchemaError,
        match="missing monomer",
    ):
        parser._validate_system_monomer_keys(
            data,
            data["simulations"],
            "counts",
        )


def test_missing_monomer_name_receives_data_id_name(
    parser,
):
    data = {
        "monomers": [
            {
                "smiles": "CCO",
            }
        ]
    }

    systems = [
        {
            "tag": "sim",
            "monomer_counts": {
                "data_1": 2,
            },
        }
    ]

    parser._validate_system_monomer_keys(
        data,
        systems,
        "counts",
    )

    assert (
        data["monomers"][0]["name"]
        == "data_1"
    )


# =============================================================================
# Monomer-entry validation
# =============================================================================


def test_validate_monomer_entry_counts_mode(
    parser,
):
    data = make_counts_input()

    monomers = parser._validate_monomer_entry(
        data,
        "counts",
        data["simulations"],
    )

    assert len(monomers) == 2

    ethanol = monomers[0]

    assert ethanol.id == 1
    assert ethanol.data_id == "data_1"
    assert ethanol.name == "ethanol"
    assert ethanol.smiles == "CCO"

    assert ethanol.count == {
        "small": 10,
        "large": 100,
    }

    assert ethanol.ratio is None
    assert isinstance(ethanol.rdkit_mol, Chem.Mol)
    assert ethanol.num_atoms is not None
    assert ethanol.molecular_weight is not None


def test_validate_monomer_entry_ratio_mode(
    parser,
):
    data = make_ratio_input()

    monomers = parser._validate_monomer_entry(
        data,
        "ratio",
        data["simulations"],
    )

    ethanol = monomers[0]

    assert ethanol.count is None
    assert ethanol.ratio == 1.0


def test_validate_monomer_entry_requires_list(
    parser,
):
    with pytest.raises(InputSchemaError):
        parser._validate_monomer_entry(
            {
                "monomers": {},
            },
            "counts",
            [],
        )


def test_validate_monomer_entry_rejects_non_dictionary_entry(
    parser,
):
    with pytest.raises(InputSchemaError):
        parser._validate_monomer_entry(
            {
                "monomers": [
                    "bad entry"
                ],
            },
            "counts",
            [],
        )


# =============================================================================
# Legacy composition validator
# =============================================================================


def test_validate_composition_accepts_counts_targets(
    parser,
):
    composition = {
        "targets": [
            {
                "tag": "one",
            },
            {
                "tag": "two",
            },
        ]
    }

    assert (
        parser._validate_composition(
            composition,
            "counts",
        )
        is composition
    )


def test_validate_composition_requires_targets(
    parser,
):
    with pytest.raises(InputSchemaError):
        parser._validate_composition(
            {},
            "counts",
        )


def test_validate_composition_rejects_duplicate_tags(
    parser,
):
    with pytest.raises(InputSchemaError):
        parser._validate_composition(
            {
                "targets": [
                    {
                        "tag": "same",
                    },
                    {
                        "tag": "same",
                    },
                ]
            },
            "counts",
        )


def test_validate_composition_ratio_requires_total_atoms(
    parser,
):
    with pytest.raises(NumericFieldError):
        parser._validate_composition(
            {
                "targets": [
                    {
                        "tag": "one",
                    }
                ]
            },
            "ratio",
        )


def test_validate_composition_counts_rejects_total_atoms(
    parser,
):
    with pytest.raises(InputSchemaError):
        parser._validate_composition(
            {
                "targets": [
                    {
                        "tag": "one",
                        "total_atoms": 1000,
                    }
                ]
            },
            "counts",
        )


# =============================================================================
# Legacy numeric validator
# =============================================================================


def test_legacy_numeric_validator_accepts_valid_input(
    parser,
):
    parser._validate_numeric_fields(
        {
            "density": 1.0,
            "temperature": [
                300,
                400,
            ],
            "number_of_monomers": {
                "a": 5,
                "b": 10,
            },
        }
    )


def test_legacy_numeric_validator_rejects_bad_density(
    parser,
):
    with pytest.raises(NumericFieldError):
        parser._validate_numeric_fields(
            {
                "density": 0,
                "temperature": 300,
                "number_of_monomers": {
                    "a": 1,
                },
            }
        )


def test_legacy_numeric_validator_rejects_bad_temperature(
    parser,
):
    with pytest.raises(NumericFieldError):
        parser._validate_numeric_fields(
            {
                "density": 1.0,
                "temperature": -10,
                "number_of_monomers": {
                    "a": 1,
                },
            }
        )


def test_legacy_numeric_validator_rejects_bad_monomer_count(
    parser,
):
    with pytest.raises(NumericFieldError):
        parser._validate_numeric_fields(
            {
                "density": 1.0,
                "temperature": 300,
                "number_of_monomers": {
                    "a": 0,
                },
            }
        )


# =============================================================================
# Full validate_inputs(): counts mode
# =============================================================================


def test_validate_inputs_counts_mode_end_to_end(
    parser,
):
    data = make_counts_input()

    result = parser.validate_inputs(data)

    assert isinstance(result, SimulationSetup)

    assert (
        result.simulation_name
        == "counts_demo"
    )

    assert result.composition_method == "counts"

    assert result.temperature == [
        300.0,
        400.0,
    ]

    assert result.density == [
        0.8,
        1.0,
    ]

    assert result.force_field == "PCFF"
    assert len(result.monomers) == 2
    assert len(result.simulations) == 2

    assert result.loop is True
    assert result.max_loop_count == 5
    assert result.reaction_iteration_depth == 5

    assert result.deep_search is True
    assert result.wildcards is True

    assert (
        result.deduplicate_reaction_templates
        is True
    )

    assert (
        result.write_second_reaction_stage
        is True
    )


def test_validate_inputs_preserves_raw_input_reference(
    parser,
):
    data = make_counts_input()

    result = parser.validate_inputs(data)

    assert result.input_json is data


def test_validate_inputs_zero_iteration_depth_disables_loop(
    parser,
):
    data = make_counts_input()

    data["reaction_iteration_depth"] = 0

    result = parser.validate_inputs(data)

    assert result.loop is False
    assert result.max_loop_count is None
    assert result.reaction_iteration_depth == 0


def test_validate_inputs_workflow_options_can_be_disabled(
    parser,
):
    data = make_counts_input()

    data.update(
        {
            "deep_search": False,
            "wildcards": False,
            "deduplicate_reaction_templates": False,
            "write_second_reaction_stage": False,
        }
    )

    result = parser.validate_inputs(data)

    assert result.deep_search is False
    assert result.wildcards is False

    assert (
        result.deduplicate_reaction_templates
        is False
    )

    assert (
        result.write_second_reaction_stage
        is False
    )


def test_validate_inputs_accepts_workflow_aliases(
    parser,
):
    data = make_counts_input()

    data["deep-search"] = "no"
    data["use-wildcards"] = "false"
    data["template_dedup"] = "off"
    data["stage_2"] = "0"
    data["iterations"] = "3"

    result = parser.validate_inputs(data)

    assert result.deep_search is False
    assert result.wildcards is False

    assert (
        result.deduplicate_reaction_templates
        is False
    )

    assert (
        result.write_second_reaction_stage
        is False
    )

    assert result.reaction_iteration_depth == 3
    assert result.max_loop_count == 3


# =============================================================================
# Full validate_inputs(): ratio mode
# =============================================================================


def test_validate_inputs_ratio_mode_end_to_end(
    parser,
):
    data = make_ratio_input()

    result = parser.validate_inputs(data)

    assert isinstance(result, SimulationSetup)

    assert result.composition_method == "ratio"

    assert result.temperature == [
        300.0,
        400.0,
    ]

    assert result.density == [
        0.8,
        1.0,
    ]

    assert len(result.monomers) == 2
    assert len(result.simulations) == 2

    assert result.monomers[0].count is None
    assert result.monomers[0].ratio == 1.0

    assert (
        result.simulations[0].total_atoms
        == 10000
    )


# =============================================================================
# Public-path malformed input checks
# =============================================================================


def test_validate_inputs_rejects_non_dictionary_monomer_with_schema_error(
    parser,
):
    """
    Public validation should translate malformed monomer entries into the
    parser's own schema exception rather than leaking AttributeError.
    """
    data = make_counts_input()

    data["monomers"] = [
        "not-a-dictionary"
    ]

    data["simulations"] = [
        {
            "tag": "sim",
            "temperature": 300,
            "density": 1.0,
            "monomer_counts": {
                "data_1": 1,
            },
        }
    ]

    with pytest.raises(InputSchemaError):
        parser.validate_inputs(data)


@pytest.mark.parametrize(
    "simulation_name",
    [
        "",
        "   ",
        None,
        123,
    ],
)
def test_validate_inputs_requires_nonempty_string_simulation_name(
    parser,
    simulation_name,
):
    """
    simulation_name is part of the public schema and later becomes an
    output-directory name, so invalid values should fail here.
    """
    data = make_counts_input()

    data["simulation_name"] = simulation_name

    with pytest.raises(InputSchemaError):
        parser.validate_inputs(data)


# =============================================================================
# Molecule representation / image generation
# =============================================================================


def test_molecule_representation_returns_molecules_and_legends(
    parser,
):
    mol_1 = Chem.MolFromSmiles("CCO")
    mol_2 = Chem.MolFromSmiles("CCN")

    setup = SimulationSetup(
        simulation_name="demo",
        temperature=[300.0],
        density=[1.0],
        force_field="PCFF",
        monomers=[
            MonomerEntry(
                id=1,
                data_id="data_1",
                name="ethanol",
                smiles="CCO",
                count={"sim": 1},
                ratio=None,
                rdkit_mol=mol_1,
            ),
            MonomerEntry(
                id=2,
                data_id="data_2",
                name="ethylamine",
                smiles="CCN",
                count={"sim": 1},
                ratio=None,
                rdkit_mol=mol_2,
            ),
        ],
    )

    molecules, legends = (
        parser.molecule_representation_of_initial_molecules(
            setup
        )
    )

    assert molecules == [
        mol_1,
        mol_2,
    ]

    assert legends == [
        "ethanol",
        "ethylamine",
    ]


def test_molecule_representation_uses_data_id_when_name_missing(
    parser,
):
    mol = Chem.MolFromSmiles("CCO")

    setup = SimulationSetup(
        simulation_name="demo",
        temperature=[300.0],
        density=[1.0],
        force_field="PCFF",
        monomers=[
            MonomerEntry(
                id=1,
                data_id="data_1",
                name=None,
                smiles="CCO",
                count={"sim": 1},
                ratio=None,
                rdkit_mol=mol,
            ),
        ],
    )

    _, legends = (
        parser.molecule_representation_of_initial_molecules(
            setup
        )
    )

    assert legends == [
        "data_1"
    ]


def test_initial_molecules_image_grid_uses_expected_draw_settings(
    parser,
    monkeypatch,
):
    mol_1 = Chem.MolFromSmiles("CCO")
    mol_2 = Chem.MolFromSmiles("CCN")

    setup = SimulationSetup(
        simulation_name="demo",
        temperature=[300.0],
        density=[1.0],
        force_field="PCFF",
        monomers=[
            MonomerEntry(
                id=1,
                data_id="data_1",
                name="ethanol",
                smiles="CCO",
                count={"sim": 1},
                ratio=None,
                rdkit_mol=mol_1,
            ),
            MonomerEntry(
                id=2,
                data_id="data_2",
                name="ethylamine",
                smiles="CCN",
                count={"sim": 1},
                ratio=None,
                rdkit_mol=mol_2,
            ),
        ],
    )

    session = SimpleNamespace(
        inputs=setup
    )

    captured = {}

    fake_image = object()

    def fake_grid(
        molecules,
        *,
        molsPerRow,
        subImgSize,
        legends,
    ):
        captured["molecules"] = molecules
        captured["molsPerRow"] = molsPerRow
        captured["subImgSize"] = subImgSize
        captured["legends"] = legends

        return fake_image

    monkeypatch.setattr(
        "AutoREACTER.input_parser.Draw.MolsToGridImage",
        fake_grid,
    )

    result = parser.initial_molecules_image_grid(
        session
    )

    assert result is fake_image

    assert captured["molecules"] == [
        mol_1,
        mol_2,
    ]

    assert captured["molsPerRow"] == 3

    assert captured["subImgSize"] == (
        400,
        400,
    )

    assert captured["legends"] == [
        "ethanol",
        "ethylamine",
    ]