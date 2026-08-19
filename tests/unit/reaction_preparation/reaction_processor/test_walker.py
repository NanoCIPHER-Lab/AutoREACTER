from rdkit import Chem

from AutoREACTER.reaction_preparation.reaction_processor.walker import (
    get_new_neighbors,
    product_atom_walker,
    reactant_atom_walker,
    reaction_atom_walker,
)


# =============================================================================
# Helpers
# =============================================================================


def mol(smiles: str) -> Chem.Mol:
    molecule = Chem.MolFromSmiles(smiles)

    assert molecule is not None

    return molecule


# =============================================================================
# get_new_neighbors
# =============================================================================


def test_get_new_neighbors_returns_unvisited_neighbors():
    # 0 - 1 - 2 - 3
    molecule = mol("CCCC")

    visited = {
        1: [1],
        2: [],
    }

    result = get_new_neighbors(
        molecule,
        1,
        visited,
    )

    assert result == [0, 2]


def test_get_new_neighbors_excludes_atoms_from_previous_shells():
    molecule = mol("CCCC")

    visited = {
        1: [0],
        2: [1],
        3: [],
    }

    result = get_new_neighbors(
        molecule,
        1,
        visited,
    )

    assert result == [2]


def test_get_new_neighbors_excludes_atoms_from_any_visited_shell():
    molecule = mol("CCCCC")

    visited = {
        1: [0],
        2: [4],
        3: [2],
        4: [],
    }

    result = get_new_neighbors(
        molecule,
        2,
        visited,
    )

    # Atom 2 is connected to atoms 1 and 3.
    assert result == [1, 3]


def test_get_new_neighbors_returns_empty_for_isolated_atom():
    molecule = mol("[He]")

    visited = {
        1: [0],
    }

    result = get_new_neighbors(
        molecule,
        0,
        visited,
    )

    assert result == []


def test_get_new_neighbors_does_not_revisit_ring_atoms():
    # Cyclohexane:
    #
    # 0 connected to 1 and 5.
    molecule = mol("C1CCCCC1")

    visited = {
        1: [0],
        2: [1],
    }

    result = get_new_neighbors(
        molecule,
        1,
        visited,
    )

    assert result == [2]


def test_get_new_neighbors_does_not_modify_visited():
    molecule = mol("CCCC")

    visited = {
        1: [1],
        2: [],
    }

    original = {
        key: value.copy()
        for key, value
        in visited.items()
    }

    get_new_neighbors(
        molecule,
        1,
        visited,
    )

    assert visited == original


# =============================================================================
# reactant_atom_walker
# =============================================================================


def test_reactant_atom_walker_accepts_single_integer_start():
    molecule = mol("CCC")

    template_atoms, edge_atoms = (
        reactant_atom_walker(
            molecule,
            1,
            max_bonds=1,
        )
    )

    # max_bonds=1 means start atom plus atoms exactly one bond away.
    assert template_atoms == [
        1,
        0,
        2,
    ]

    assert edge_atoms == [
        0,
        2,
    ]


def test_reactant_atom_walker_accepts_list_start():
    molecule = mol("CCCCC")

    template_atoms, edge_atoms = (
        reactant_atom_walker(
            molecule,
            [1, 3],
            max_bonds=1,
        )
    )

    assert template_atoms == [
        1,
        3,
        0,
        2,
        4,
    ]

    assert edge_atoms == [
        0,
        2,
        4,
    ]


def test_reactant_atom_walker_accepts_tuple_start():
    molecule = mol("CCCCC")

    template_atoms, edge_atoms = (
        reactant_atom_walker(
            molecule,
            (1, 3),
            max_bonds=1,
        )
    )

    assert template_atoms == [
        1,
        3,
        0,
        2,
        4,
    ]

    assert edge_atoms == [
        0,
        2,
        4,
    ]


def test_reactant_atom_walker_one_bond_distance():
    # Branched molecule:
    #
    #       2
    #       |
    #   0 - 1 - 3
    #
    molecule = mol("CC(C)C")

    template_atoms, edge_atoms = (
        reactant_atom_walker(
            molecule,
            1,
            max_bonds=1,
        )
    )

    assert template_atoms == [
        1,
        0,
        2,
        3,
    ]

    assert edge_atoms == [
        0,
        2,
        3,
    ]


def test_reactant_atom_walker_two_bond_distance():
    # 0 - 1 - 2 - 3 - 4
    molecule = mol("CCCCC")

    template_atoms, edge_atoms = (
        reactant_atom_walker(
            molecule,
            2,
            max_bonds=2,
        )
    )

    assert template_atoms == [
        2,
        1,
        3,
        0,
        4,
    ]

    assert edge_atoms == [
        0,
        4,
    ]


def test_reactant_atom_walker_three_bond_distance():
    # 0 - 1 - 2 - 3 - 4 - 5 - 6
    molecule = mol("CCCCCCC")

    template_atoms, edge_atoms = (
        reactant_atom_walker(
            molecule,
            3,
            max_bonds=3,
        )
    )

    assert template_atoms == [
        3,
        2,
        4,
        1,
        5,
        0,
        6,
    ]

    assert edge_atoms == [
        0,
        6,
    ]


def test_reactant_atom_walker_ring_does_not_revisit_atoms():
    # Six-membered ring.
    molecule = mol("C1CCCCC1")

    template_atoms, edge_atoms = (
        reactant_atom_walker(
            molecule,
            0,
            max_bonds=2,
        )
    )

    # Distance 0: 0
    # Distance 1: 1, 5
    # Distance 2: 2, 4
    assert template_atoms == [
        0,
        1,
        5,
        2,
        4,
    ]

    assert edge_atoms == [
        2,
        4,
    ]

    assert len(
        template_atoms
    ) == len(
        set(template_atoms)
    )


def test_reactant_atom_walker_stops_at_requested_distance():
    molecule = mol(
        "CCCCCCC"
    )

    template_atoms, edge_atoms = (
        reactant_atom_walker(
            molecule,
            3,
            max_bonds=2,
        )
    )

    # Atoms 0 and 6 are three bonds away and must not be included.
    assert 0 not in template_atoms
    assert 6 not in template_atoms

    assert edge_atoms == [
        1,
        5,
    ]


def test_reactant_atom_walker_does_not_cross_disconnected_components():
    molecule = mol(
        "CCC.CCC"
    )

    template_atoms, edge_atoms = (
        reactant_atom_walker(
            molecule,
            0,
            max_bonds=4,
        )
    )

    assert set(
        template_atoms
    ) == {
        0,
        1,
        2,
    }

    # No atoms exist exactly four bonds from the seed.
    assert edge_atoms == []


def test_reactant_atom_walker_multiple_seeds_do_not_duplicate_shared_neighbor():
    # 0 - 1 - 2 - 3 - 4
    #
    # Seeds 1 and 3 both reach atom 2.
    molecule = mol(
        "CCCCC"
    )

    template_atoms, edge_atoms = (
        reactant_atom_walker(
            molecule,
            [1, 3],
            max_bonds=1,
        )
    )

    assert template_atoms.count(
        2
    ) == 1

    assert edge_atoms.count(
        2
    ) == 1


def test_reactant_atom_walker_preserves_seed_order():
    molecule = mol(
        "CCCCC"
    )

    template_atoms, _ = (
        reactant_atom_walker(
            molecule,
            [3, 1],
            max_bonds=1,
        )
    )

    assert template_atoms[:2] == [
        3,
        1,
    ]


# =============================================================================
# product_atom_walker
# =============================================================================


def test_product_atom_walker_maps_all_available_atoms():
    template_indexes = [
        0,
        2,
        4,
    ]

    mapping = {
        0: 10,
        2: 12,
        4: 14,
    }

    result = product_atom_walker(
        template_indexes,
        mapping,
    )

    assert result == {
        0: 10,
        2: 12,
        4: 14,
    }


def test_product_atom_walker_skips_unmapped_atoms():
    template_indexes = [
        0,
        1,
        2,
    ]

    mapping = {
        0: 10,
        2: 12,
    }

    result = product_atom_walker(
        template_indexes,
        mapping,
    )

    assert result == {
        0: 10,
        2: 12,
    }


def test_product_atom_walker_accepts_product_index_zero():
    mapping = {
        5: 0,
    }

    result = product_atom_walker(
        [5],
        mapping,
    )

    assert result == {
        5: 0,
    }


def test_product_atom_walker_skips_none_mapping_value():
    mapping = {
        1: None,
        2: 5,
    }

    result = product_atom_walker(
        [1, 2],
        mapping,
    )

    assert result == {
        2: 5,
    }


def test_product_atom_walker_empty_template_returns_empty_dict():
    result = product_atom_walker(
        [],
        {
            0: 10,
        },
    )

    assert result == {}


def test_product_atom_walker_empty_mapping_returns_empty_dict():
    result = product_atom_walker(
        [
            0,
            1,
        ],
        {},
    )

    assert result == {}


def test_product_atom_walker_preserves_template_order():
    result = product_atom_walker(
        [
            3,
            1,
            2,
        ],
        {
            1: 10,
            2: 20,
            3: 30,
        },
    )

    assert list(
        result.keys()
    ) == [
        3,
        1,
        2,
    ]


def test_product_atom_walker_does_not_modify_mapping():
    mapping = {
        0: 10,
        1: 11,
    }

    original = mapping.copy()

    product_atom_walker(
        [
            0,
            1,
        ],
        mapping,
    )

    assert mapping == original


# =============================================================================
# reaction_atom_walker
# =============================================================================


def test_reaction_atom_walker_combines_graph_walk_and_mapping():
    # 0 - 1 - 2 - 3 - 4
    molecule = mol(
        "CCCCC"
    )

    mapping = {
        0: 100,
        1: 101,
        2: 102,
        3: 103,
        4: 104,
    }

    mapped, edge_atoms = (
        reaction_atom_walker(
            molecule,
            2,
            mapping,
            max_bonds=2,
        )
    )

    assert mapped == {
        2: 102,
        1: 101,
        3: 103,
        0: 100,
        4: 104,
    }

    assert edge_atoms == [
        0,
        4,
    ]


def test_reaction_atom_walker_keeps_edge_atoms_in_reactant_index_space():
    molecule = mol(
        "CCCCC"
    )

    mapping = {
        0: 50,
        1: 51,
        2: 52,
        3: 53,
        4: 54,
    }

    _, edge_atoms = (
        reaction_atom_walker(
            molecule,
            2,
            mapping,
            max_bonds=2,
        )
    )

    assert edge_atoms == [
        0,
        4,
    ]

    assert 50 not in edge_atoms
    assert 54 not in edge_atoms


def test_reaction_atom_walker_skips_template_atoms_missing_from_mapping():
    molecule = mol(
        "CCCCC"
    )

    mapping = {
        1: 11,
        2: 12,
        3: 13,
    }

    mapped, edge_atoms = (
        reaction_atom_walker(
            molecule,
            2,
            mapping,
            max_bonds=2,
        )
    )

    assert mapped == {
        2: 12,
        1: 11,
        3: 13,
    }

    # Edge detection is based on the reactant graph, independent of
    # whether those atoms happen to have product mappings.
    assert edge_atoms == [
        0,
        4,
    ]


def test_reaction_atom_walker_multiple_start_atoms():
    molecule = mol(
        "CCCCC"
    )

    mapping = {
        idx: idx + 100
        for idx in range(5)
    }

    mapped, edge_atoms = (
        reaction_atom_walker(
            molecule,
            [1, 3],
            mapping,
            max_bonds=1,
        )
    )

    assert mapped == {
        1: 101,
        3: 103,
        0: 100,
        2: 102,
        4: 104,
    }

    assert edge_atoms == [
        0,
        2,
        4,
    ]


def test_reaction_atom_walker_disconnected_components():
    molecule = mol(
        "CC.CC"
    )

    mapping = {
        0: 10,
        1: 11,
        2: 12,
        3: 13,
    }

    mapped, edge_atoms = (
        reaction_atom_walker(
            molecule,
            0,
            mapping,
            max_bonds=3,
        )
    )

    assert mapped == {
        0: 10,
        1: 11,
    }

    assert edge_atoms == []


def test_reaction_atom_walker_does_not_modify_mapping():
    molecule = mol(
        "CCC"
    )

    mapping = {
        0: 10,
        1: 11,
        2: 12,
    }

    original = mapping.copy()

    reaction_atom_walker(
        molecule,
        1,
        mapping,
        max_bonds=1,
    )

    assert mapping == original