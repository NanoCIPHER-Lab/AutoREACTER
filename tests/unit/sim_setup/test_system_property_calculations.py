import math
from types import SimpleNamespace

import pytest
from rdkit import Chem

from AutoREACTER.sim_setup.system_property_calculations import (
    CM_2_A3,
    N_A,
    NoneMonomerError,
    SystemPropertyCalculations,
)


# =============================================================================
# Helpers
# =============================================================================


def make_monomer(
    *,
    name="monomer",
    monomer_id=1,
    status=True,
    rdkit_mol=None,
    num_atoms=None,
    molecular_weight=None,
    count=None,
):
    return SimpleNamespace(
        id=monomer_id,
        name=name,
        status=status,
        rdkit_mol=rdkit_mol,
        num_atoms=num_atoms,
        molecular_weight=molecular_weight,
        count=count,
    )


def make_simulation(
    *,
    tag="sim1",
    total_atoms=100,
    density=1.0,
    monomer_ratios=None,
    monomer_counts=None,
):
    return SimpleNamespace(
        tag=tag,
        total_atoms=total_atoms,
        density=density,
        monomer_ratios=(
            {}
            if monomer_ratios is None
            else monomer_ratios
        ),
        monomer_counts=monomer_counts,
        initial_box_volume=None,
        initial_box_length=None,
    )


def make_setup(
    *,
    monomers=None,
    simulations=None,
    composition_method="ratio",
):
    return SimpleNamespace(
        monomers=list(
            monomers or []
        ),
        simulations=list(
            simulations or []
        ),
        composition_method=composition_method,
    )


# =============================================================================
# Constants
# =============================================================================


def test_avogadro_constant():
    assert N_A == pytest.approx(
        6.02214076e23
    )


def test_cm3_to_angstrom3_conversion():
    assert CM_2_A3 == pytest.approx(
        1.0e24
    )


# =============================================================================
# Exception
# =============================================================================


def test_none_monomer_error_is_exception():
    assert issubclass(
        NoneMonomerError,
        Exception,
    )


# =============================================================================
# Constructor
# =============================================================================


def test_constructor_stores_simulation_setup():
    setup = make_setup()

    calculator = (
        SystemPropertyCalculations(
            setup
        )
    )

    assert (
        calculator.simulation_setup
        is setup
    )


# =============================================================================
# process_all
# =============================================================================


def test_process_all_calls_steps_in_order(
    monkeypatch,
):
    setup = make_setup()

    calculator = (
        SystemPropertyCalculations(
            setup
        )
    )

    events = []

    monkeypatch.setattr(
        calculator,
        "_populate_monomer_properties",
        lambda:
        events.append(
            "populate"
        ),
    )

    monkeypatch.setattr(
        calculator,
        "_calculate_replica_properties",
        lambda:
        events.append(
            "calculate"
        ),
    )

    result = (
        calculator.process_all()
    )

    assert events == [
        "populate",
        "calculate",
    ]

    assert result is setup


def test_process_all_returns_same_setup_object(
    monkeypatch,
):
    setup = make_setup()

    calculator = (
        SystemPropertyCalculations(
            setup
        )
    )

    monkeypatch.setattr(
        calculator,
        "_populate_monomer_properties",
        lambda: None,
    )

    monkeypatch.setattr(
        calculator,
        "_calculate_replica_properties",
        lambda: None,
    )

    result = (
        calculator.process_all()
    )

    assert result is setup


# =============================================================================
# _populate_monomer_properties
# =============================================================================


def test_populate_active_monomer_counts_explicit_hydrogens():
    methane = (
        Chem.MolFromSmiles(
            "C"
        )
    )

    monomer = make_monomer(
        name="methane",
        rdkit_mol=methane,
    )

    setup = make_setup(
        monomers=[
            monomer
        ]
    )

    calculator = (
        SystemPropertyCalculations(
            setup
        )
    )

    calculator._populate_monomer_properties()

    # CH4 = one carbon + four hydrogens.
    assert monomer.num_atoms == 5


def test_populate_does_not_modify_original_rdkit_molecule():
    methane = (
        Chem.MolFromSmiles(
            "C"
        )
    )

    assert (
        methane.GetNumAtoms()
        == 1
    )

    monomer = make_monomer(
        name="methane",
        rdkit_mol=methane,
    )

    calculator = (
        SystemPropertyCalculations(
            make_setup(
                monomers=[
                    monomer
                ]
            )
        )
    )

    calculator._populate_monomer_properties()

    # AddHs() is applied to a temporary molecule.
    assert (
        methane.GetNumAtoms()
        == 1
    )

    assert monomer.num_atoms == 5


def test_populate_sets_molecular_weight():
    monomer = make_monomer(
        name="methane",
        rdkit_mol=(
            Chem.MolFromSmiles(
                "C"
            )
        ),
    )

    calculator = (
        SystemPropertyCalculations(
            make_setup(
                monomers=[
                    monomer
                ]
            )
        )
    )

    calculator._populate_monomer_properties()

    assert monomer.molecular_weight == pytest.approx(
        16.043,
        rel=1.0e-3,
    )


def test_populate_ethanol_full_atom_count():
    # C2H6O = 9 atoms total.
    monomer = make_monomer(
        name="ethanol",
        rdkit_mol=(
            Chem.MolFromSmiles(
                "CCO"
            )
        ),
    )

    calculator = (
        SystemPropertyCalculations(
            make_setup(
                monomers=[
                    monomer
                ]
            )
        )
    )

    calculator._populate_monomer_properties()

    assert monomer.num_atoms == 9


def test_populate_skips_inactive_monomer():
    monomer = make_monomer(
        name="inactive",
        status=False,
        rdkit_mol=None,
        num_atoms=999,
        molecular_weight=888,
    )

    calculator = (
        SystemPropertyCalculations(
            make_setup(
                monomers=[
                    monomer
                ]
            )
        )
    )

    calculator._populate_monomer_properties()

    assert monomer.num_atoms == 999

    assert (
        monomer.molecular_weight
        == 888
    )


def test_populate_active_none_molecule_raises():
    monomer = make_monomer(
        name="broken",
        monomer_id=42,
        status=True,
        rdkit_mol=None,
    )

    calculator = (
        SystemPropertyCalculations(
            make_setup(
                monomers=[
                    monomer
                ]
            )
        )
    )

    with pytest.raises(
        NoneMonomerError,
        match=(
            "Monomer with ID 42 "
            "has no RDKit Mol object"
        ),
    ):
        calculator._populate_monomer_properties()


def test_populate_prints_info_message(
    capsys,
):
    calculator = (
        SystemPropertyCalculations(
            make_setup()
        )
    )

    calculator._populate_monomer_properties()

    assert (
        "Populating monomer properties"
        in capsys.readouterr().out
    )


def test_populate_multiple_active_monomers():
    methane = make_monomer(
        name="methane",
        rdkit_mol=(
            Chem.MolFromSmiles(
                "C"
            )
        ),
    )

    ethane = make_monomer(
        name="ethane",
        rdkit_mol=(
            Chem.MolFromSmiles(
                "CC"
            )
        ),
    )

    calculator = (
        SystemPropertyCalculations(
            make_setup(
                monomers=[
                    methane,
                    ethane,
                ]
            )
        )
    )

    calculator._populate_monomer_properties()

    assert methane.num_atoms == 5

    # C2H6
    assert ethane.num_atoms == 8


# =============================================================================
# _get_monomer_counts
# =============================================================================


def test_get_monomer_counts_requires_total_atoms():
    monomer = make_monomer(
        name="A",
        num_atoms=10,
    )

    simulation = make_simulation(
        tag="test",
        total_atoms=None,
        monomer_ratios={
            "A": 1,
        },
    )

    calculator = (
        SystemPropertyCalculations(
            make_setup()
        )
    )

    with pytest.raises(
        ValueError,
        match=(
            "Simulation 'test' is "
            "missing 'total_atoms'"
        ),
    ):
        calculator._get_monomer_counts(
            simulation,
            {
                "A": monomer
            },
        )


def test_get_monomer_counts_basic_ratio_calculation():
    """
    atoms per ratio unit:

        A: ratio 2 * 10 atoms = 20
        B: ratio 1 * 20 atoms = 20
                               ----
                                40

    total_atoms = 100
    multiplier = 100 / 40 = 2.5

    A count = ceil(2 * 2.5) = 5
    B count = ceil(1 * 2.5) = 3
    """
    monomer_a = make_monomer(
        name="A",
        num_atoms=10,
    )

    monomer_b = make_monomer(
        name="B",
        num_atoms=20,
    )

    simulation = make_simulation(
        tag="sim",
        total_atoms=100,
        monomer_ratios={
            "A": 2,
            "B": 1,
        },
    )

    calculator = (
        SystemPropertyCalculations(
            make_setup()
        )
    )

    calculator._get_monomer_counts(
        simulation,
        {
            "A": monomer_a,
            "B": monomer_b,
        },
    )

    assert simulation.monomer_counts == {
        "A": 5,
        "B": 3,
    }


def test_ratio_counting_can_overshoot_requested_total_atoms():
    """
    Characterization of current ceil-based behavior.

    5 A molecules * 10 atoms = 50
    3 B molecules * 20 atoms = 60
    total actual atoms = 110

    Requested total_atoms was 100.
    """
    monomer_a = make_monomer(
        name="A",
        num_atoms=10,
    )

    monomer_b = make_monomer(
        name="B",
        num_atoms=20,
    )

    simulation = make_simulation(
        total_atoms=100,
        monomer_ratios={
            "A": 2,
            "B": 1,
        },
    )

    calculator = (
        SystemPropertyCalculations(
            make_setup()
        )
    )

    calculator._get_monomer_counts(
        simulation,
        {
            "A": monomer_a,
            "B": monomer_b,
        },
    )

    actual_atoms = (
        simulation.monomer_counts["A"]
        * monomer_a.num_atoms
        + simulation.monomer_counts["B"]
        * monomer_b.num_atoms
    )

    assert actual_atoms == 110

    assert actual_atoms >= (
        simulation.total_atoms
    )


def test_get_monomer_counts_initializes_none_dictionary():
    monomer = make_monomer(
        name="A",
        num_atoms=10,
    )

    simulation = make_simulation(
        total_atoms=100,
        monomer_ratios={
            "A": 1,
        },
        monomer_counts=None,
    )

    calculator = (
        SystemPropertyCalculations(
            make_setup()
        )
    )

    calculator._get_monomer_counts(
        simulation,
        {
            "A": monomer
        },
    )

    assert isinstance(
        simulation.monomer_counts,
        dict,
    )

    assert (
        simulation.monomer_counts["A"]
        == 10
    )


def test_get_monomer_counts_preserves_existing_dictionary_entries():
    monomer = make_monomer(
        name="A",
        num_atoms=10,
    )

    simulation = make_simulation(
        total_atoms=100,
        monomer_ratios={
            "A": 1,
        },
        monomer_counts={
            "existing": 99
        },
    )

    calculator = (
        SystemPropertyCalculations(
            make_setup()
        )
    )

    calculator._get_monomer_counts(
        simulation,
        {
            "A": monomer
        },
    )

    assert simulation.monomer_counts == {
        "existing": 99,
        "A": 10,
    }


def test_get_monomer_counts_records_count_on_monomer():
    monomer = make_monomer(
        name="A",
        num_atoms=10,
        count=None,
    )

    simulation = make_simulation(
        tag="500K",
        total_atoms=100,
        monomer_ratios={
            "A": 1,
        },
    )

    calculator = (
        SystemPropertyCalculations(
            make_setup()
        )
    )

    calculator._get_monomer_counts(
        simulation,
        {
            "A": monomer
        },
    )

    assert monomer.count == {
        "500K": 10
    }


def test_get_monomer_counts_preserves_other_simulation_counts():
    monomer = make_monomer(
        name="A",
        num_atoms=10,
        count={
            "old_sim": 4
        },
    )

    simulation = make_simulation(
        tag="new_sim",
        total_atoms=100,
        monomer_ratios={
            "A": 1,
        },
    )

    calculator = (
        SystemPropertyCalculations(
            make_setup()
        )
    )

    calculator._get_monomer_counts(
        simulation,
        {
            "A": monomer
        },
    )

    assert monomer.count == {
        "old_sim": 4,
        "new_sim": 10,
    }


def test_get_monomer_counts_ignores_ratio_names_not_active():
    monomer = make_monomer(
        name="A",
        num_atoms=10,
    )

    simulation = make_simulation(
        total_atoms=100,
        monomer_ratios={
            "A": 1,
            "inactive_or_unknown": 999,
        },
    )

    calculator = (
        SystemPropertyCalculations(
            make_setup()
        )
    )

    calculator._get_monomer_counts(
        simulation,
        {
            "A": monomer
        },
    )

    assert simulation.monomer_counts == {
        "A": 10
    }


def test_get_monomer_counts_no_active_ratio_monomers_raises():
    simulation = make_simulation(
        tag="empty",
        total_atoms=100,
        monomer_ratios={
            "missing": 1,
        },
    )

    calculator = (
        SystemPropertyCalculations(
            make_setup()
        )
    )

    with pytest.raises(
        ValueError,
        match=(
            "has no active monomers "
            "to calculate ratios"
        ),
    ):
        calculator._get_monomer_counts(
            simulation,
            {},
        )


def test_get_monomer_counts_minimum_is_one():
    monomer_a = make_monomer(
        name="A",
        num_atoms=1000,
    )

    monomer_b = make_monomer(
        name="B",
        num_atoms=1,
    )

    simulation = make_simulation(
        total_atoms=1,
        monomer_ratios={
            "A": 1,
            "B": 1,
        },
    )

    calculator = (
        SystemPropertyCalculations(
            make_setup()
        )
    )

    calculator._get_monomer_counts(
        simulation,
        {
            "A": monomer_a,
            "B": monomer_b,
        },
    )

    assert (
        simulation.monomer_counts["A"]
        >= 1
    )

    assert (
        simulation.monomer_counts["B"]
        >= 1
    )


# =============================================================================
# _calculate_box_dimensions
# =============================================================================


def test_calculate_box_dimensions_requires_counts():
    simulation = make_simulation(
        tag="missing_counts",
        monomer_counts=None,
    )

    calculator = (
        SystemPropertyCalculations(
            make_setup()
        )
    )

    with pytest.raises(
        ValueError,
        match=(
            "Simulation 'missing_counts' "
            "lacks monomer counts"
        ),
    ):
        calculator._calculate_box_dimensions(
            simulation,
            {},
        )


def test_calculate_box_dimensions_exact_formula():
    monomer = make_monomer(
        name="A",
        molecular_weight=100.0,
    )

    simulation = make_simulation(
        density=1.0,
        monomer_counts={
            "A": 10,
        },
    )

    calculator = (
        SystemPropertyCalculations(
            make_setup()
        )
    )

    calculator._calculate_box_dimensions(
        simulation,
        {
            "A": monomer
        },
    )

    total_mass_g = (
        10 * 100.0
    ) / N_A

    initial_density = (
        1.0 / 4.0
    )

    expected_volume = (
        total_mass_g
        / initial_density
        * CM_2_A3
    )

    expected_length = round(
        math.pow(
            expected_volume,
            1.0 / 3.0,
        ),
        2,
    )

    assert (
        simulation.initial_box_volume
        == pytest.approx(
            expected_volume
        )
    )

    assert (
        simulation.initial_box_length
        == expected_length
    )


def test_initial_box_uses_quarter_target_density():
    monomer = make_monomer(
        name="A",
        molecular_weight=50.0,
    )

    simulation = make_simulation(
        density=2.0,
        monomer_counts={
            "A": 4,
        },
    )

    calculator = (
        SystemPropertyCalculations(
            make_setup()
        )
    )

    calculator._calculate_box_dimensions(
        simulation,
        {
            "A": monomer
        },
    )

    total_mass_g = (
        4 * 50.0
    ) / N_A

    target_volume = (
        total_mass_g
        / 2.0
        * CM_2_A3
    )

    # Initial density = target density / 4,
    # therefore initial volume = target volume * 4.
    assert (
        simulation.initial_box_volume
        == pytest.approx(
            target_volume * 4.0
        )
    )


def test_box_length_is_cube_root_of_initial_volume():
    monomer = make_monomer(
        name="A",
        molecular_weight=75.0,
    )

    simulation = make_simulation(
        density=1.2,
        monomer_counts={
            "A": 15,
        },
    )

    calculator = (
        SystemPropertyCalculations(
            make_setup()
        )
    )

    calculator._calculate_box_dimensions(
        simulation,
        {
            "A": monomer
        },
    )

    expected = round(
        simulation.initial_box_volume
        ** (
            1.0 / 3.0
        ),
        2,
    )

    assert (
        simulation.initial_box_length
        == expected
    )


def test_box_length_is_rounded_to_two_decimal_places():
    monomer = make_monomer(
        name="A",
        molecular_weight=123.456,
    )

    simulation = make_simulation(
        density=1.234,
        monomer_counts={
            "A": 17,
        },
    )

    calculator = (
        SystemPropertyCalculations(
            make_setup()
        )
    )

    calculator._calculate_box_dimensions(
        simulation,
        {
            "A": monomer
        },
    )

    assert (
        simulation.initial_box_length
        == round(
            simulation.initial_box_length,
            2,
        )
    )


def test_box_mass_sums_multiple_active_monomers():
    monomer_a = make_monomer(
        name="A",
        molecular_weight=100.0,
    )

    monomer_b = make_monomer(
        name="B",
        molecular_weight=50.0,
    )

    simulation = make_simulation(
        density=1.0,
        monomer_counts={
            "A": 2,
            "B": 4,
        },
    )

    calculator = (
        SystemPropertyCalculations(
            make_setup()
        )
    )

    calculator._calculate_box_dimensions(
        simulation,
        {
            "A": monomer_a,
            "B": monomer_b,
        },
    )

    expected_mass = (
        (2 * 100.0)
        + (4 * 50.0)
    ) / N_A

    expected_volume = (
        expected_mass
        / 0.25
        * CM_2_A3
    )

    assert (
        simulation.initial_box_volume
        == pytest.approx(
            expected_volume
        )
    )


def test_box_calculation_ignores_counts_for_unknown_or_inactive_monomers():
    monomer = make_monomer(
        name="A",
        molecular_weight=100.0,
    )

    simulation = make_simulation(
        density=1.0,
        monomer_counts={
            "A": 2,
            "unknown": 100000,
        },
    )

    calculator = (
        SystemPropertyCalculations(
            make_setup()
        )
    )

    calculator._calculate_box_dimensions(
        simulation,
        {
            "A": monomer
        },
    )

    expected_mass = (
        2 * 100.0
    ) / N_A

    expected_volume = (
        expected_mass
        / 0.25
        * CM_2_A3
    )

    assert (
        simulation.initial_box_volume
        == pytest.approx(
            expected_volume
        )
    )


def test_box_calculation_empty_active_mass_produces_zero_box():
    """
    Characterization of current behavior.

    If monomer_counts exists but none of its names correspond to an
    active monomer, total mass remains zero.
    """
    simulation = make_simulation(
        density=1.0,
        monomer_counts={
            "unknown": 5,
        },
    )

    calculator = (
        SystemPropertyCalculations(
            make_setup()
        )
    )

    calculator._calculate_box_dimensions(
        simulation,
        {},
    )

    assert (
        simulation.initial_box_volume
        == 0.0
    )

    assert (
        simulation.initial_box_length
        == 0.0
    )


# =============================================================================
# _calculate_replica_properties
# =============================================================================


def test_calculate_replica_properties_builds_only_active_monomer_lookup(
    monkeypatch,
):
    active = make_monomer(
        name="active",
        status=True,
    )

    inactive = make_monomer(
        name="inactive",
        status=False,
    )

    simulation = make_simulation()

    setup = make_setup(
        monomers=[
            active,
            inactive,
        ],
        simulations=[
            simulation
        ],
        composition_method="count",
    )

    calculator = (
        SystemPropertyCalculations(
            setup
        )
    )

    captured = []

    monkeypatch.setattr(
        calculator,
        "_calculate_box_dimensions",
        lambda sim, monomers:
        captured.append(
            monomers
        ),
    )

    calculator._calculate_replica_properties()

    assert list(
        captured[0].keys()
    ) == [
        "active"
    ]

    assert (
        captured[0]["active"]
        is active
    )


def test_ratio_mode_calculates_counts_before_box(
    monkeypatch,
):
    simulation = make_simulation()

    setup = make_setup(
        simulations=[
            simulation
        ],
        composition_method="ratio",
    )

    calculator = (
        SystemPropertyCalculations(
            setup
        )
    )

    events = []

    monkeypatch.setattr(
        calculator,
        "_get_monomer_counts",
        lambda sim, monomers:
        events.append(
            "counts"
        ),
    )

    monkeypatch.setattr(
        calculator,
        "_calculate_box_dimensions",
        lambda sim, monomers:
        events.append(
            "box"
        ),
    )

    calculator._calculate_replica_properties()

    assert events == [
        "counts",
        "box",
    ]


def test_non_ratio_mode_skips_count_calculation(
    monkeypatch,
):
    simulation = make_simulation(
        monomer_counts={
            "A": 1
        },
    )

    setup = make_setup(
        simulations=[
            simulation
        ],
        composition_method="count",
    )

    calculator = (
        SystemPropertyCalculations(
            setup
        )
    )

    monkeypatch.setattr(
        calculator,
        "_get_monomer_counts",
        lambda *args:
        pytest.fail(
            "_get_monomer_counts must not "
            "run in count mode"
        ),
    )

    calls = []

    monkeypatch.setattr(
        calculator,
        "_calculate_box_dimensions",
        lambda sim, monomers:
        calls.append(
            sim
        ),
    )

    calculator._calculate_replica_properties()

    assert calls == [
        simulation
    ]


def test_calculate_replica_properties_processes_all_simulations(
    monkeypatch,
):
    sim1 = make_simulation(
        tag="sim1"
    )

    sim2 = make_simulation(
        tag="sim2"
    )

    sim3 = make_simulation(
        tag="sim3"
    )

    setup = make_setup(
        simulations=[
            sim1,
            sim2,
            sim3,
        ],
        composition_method="count",
    )

    calculator = (
        SystemPropertyCalculations(
            setup
        )
    )

    processed = []

    monkeypatch.setattr(
        calculator,
        "_calculate_box_dimensions",
        lambda sim, monomers:
        processed.append(
            sim.tag
        ),
    )

    calculator._calculate_replica_properties()

    assert processed == [
        "sim1",
        "sim2",
        "sim3",
    ]


# =============================================================================
# End-to-end property calculation
# =============================================================================


def test_process_all_real_small_ratio_system():
    """
    Small integration test across:

        RDKit properties
            ↓
        ratio counts
            ↓
        mass / box calculation
    """
    methane = make_monomer(
        name="methane",
        rdkit_mol=(
            Chem.MolFromSmiles(
                "C"
            )
        ),
    )

    ethane = make_monomer(
        name="ethane",
        rdkit_mol=(
            Chem.MolFromSmiles(
                "CC"
            )
        ),
    )

    simulation = make_simulation(
        tag="small",
        total_atoms=100,
        density=1.0,
        monomer_ratios={
            "methane": 1.0,
            "ethane": 1.0,
        },
    )

    setup = make_setup(
        monomers=[
            methane,
            ethane,
        ],
        simulations=[
            simulation
        ],
        composition_method="ratio",
    )

    calculator = (
        SystemPropertyCalculations(
            setup
        )
    )

    result = (
        calculator.process_all()
    )

    assert result is setup

    assert methane.num_atoms == 5
    assert ethane.num_atoms == 8

    assert (
        simulation.monomer_counts[
            "methane"
        ]
        > 0
    )

    assert (
        simulation.monomer_counts[
            "ethane"
        ]
        > 0
    )

    assert (
        simulation.initial_box_volume
        > 0
    )

    assert (
        simulation.initial_box_length
        > 0
    )

    assert (
        methane.count["small"]
        == simulation
        .monomer_counts[
            "methane"
        ]
    )

    assert (
        ethane.count["small"]
        == simulation
        .monomer_counts[
            "ethane"
        ]
    )