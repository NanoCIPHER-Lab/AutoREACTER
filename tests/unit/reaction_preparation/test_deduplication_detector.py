from pathlib import Path
from types import SimpleNamespace

import networkx as nx
import pytest
from rdkit import Chem

from AutoREACTER.reaction_preparation.deduplication_detector import (
    DeduplicationDetector,
)


# =============================================================================
# Helpers
# =============================================================================


def mol(smiles: str) -> Chem.Mol:
    molecule = Chem.MolFromSmiles(smiles)

    assert molecule is not None

    return molecule


def make_phase_graph(
    detector,
    node_labels,
    edges=(),
):
    graph = nx.Graph()

    for node_id, label in node_labels.items():
        graph.add_node(
            node_id,
            **{
                detector.NODE_ATTRIBUTE:
                    label,
            },
        )

    for atom1, atom2, bond_label in edges:
        graph.add_edge(
            atom1,
            atom2,
            **{
                detector.EDGE_ATTRIBUTE:
                    bond_label,
            },
        )

    return graph


def make_metadata(
    reactant_smiles="CC",
    product_smiles="CC",
    template_mapping=None,
    full_mapping=None,
    first_shell=None,
    activity_stats=True,
    reaction_id=1,
    is_radical=False,
):
    reactant = mol(
        reactant_smiles
    )

    product = mol(
        product_smiles
    )

    if full_mapping is None:
        full_mapping = {
            idx: idx
            for idx in range(
                min(
                    reactant.GetNumAtoms(),
                    product.GetNumAtoms(),
                )
            )
        }

    if template_mapping is None:
        template_mapping = (
            full_mapping.copy()
        )

    if first_shell is None:
        first_shell = list(
            template_mapping.keys()
        )

    return SimpleNamespace(
        reaction_id=reaction_id,
        reactant_combined_RDmol=reactant,
        product_combined_RDmol=product,
        template_reactant_to_product_mapping=template_mapping,
        reactant_to_product_mapping=full_mapping,
        first_shell=first_shell,
        activity_stats=activity_stats,
        is_radical=is_radical,
        pre_reaction_file=None,
        post_reaction_file=None,
        map_file=None,
    )


def write_lammps_molecule(
    path: Path,
    atom_types,
    bonds=(),
):
    lines = [
        "# test molecule\n",
        "\n",
        "Coords\n",
        "\n",
    ]

    for atom_id in atom_types:
        lines.append(
            f"{atom_id} 0.0 0.0 0.0\n"
        )

    lines.extend(
        [
            "\n",
            "Types\n",
            "\n",
        ]
    )

    for atom_id, atom_type in (
        atom_types.items()
    ):
        lines.append(
            f"{atom_id} {atom_type}\n"
        )

    lines.extend(
        [
            "\n",
            "Charges\n",
            "\n",
        ]
    )

    for atom_id in atom_types:
        lines.append(
            f"{atom_id} 0.0\n"
        )

    if bonds:
        lines.extend(
            [
                "\n",
                "Bonds\n",
                "\n",
            ]
        )

        for (
            bond_id,
            bond_type,
            atom1,
            atom2,
        ) in bonds:
            lines.append(
                f"{bond_id} "
                f"{bond_type} "
                f"{atom1} "
                f"{atom2}\n"
            )

    path.write_text(
        "".join(lines),
        encoding="utf-8",
    )

    return path


def write_map_file(
    path: Path,
    equivalences=(),
    edge_ids=(),
    initiator_ids=(),
    delete_ids=(),
    wildcards=(),
):
    lines = [
        "# AutoREACTER test map\n",
        "\n",
        f"{len(edge_ids)} edgeIDs\n",
        f"{len(equivalences)} equivalences\n",
    ]

    if delete_ids:
        lines.append(
            f"{len(delete_ids)} deleteIDs\n"
        )

    lines.append(
        f"{len(wildcards)} wildcards\n"
    )

    lines.extend(
        [
            "\n",
            "InitiatorIDs\n",
            "\n",
        ]
    )

    for atom_id in initiator_ids:
        lines.append(
            f"{atom_id}\n"
        )

    lines.extend(
        [
            "\n",
            "EdgeIDs\n",
            "\n",
        ]
    )

    for atom_id in edge_ids:
        lines.append(
            f"{atom_id}\n"
        )

    lines.extend(
        [
            "\n",
            "Equivalences\n",
            "\n",
        ]
    )

    for pre_id, post_id in equivalences:
        lines.append(
            f"{pre_id} {post_id}\n"
        )

    if delete_ids:
        lines.extend(
            [
                "\n",
                "DeleteIDs\n",
                "\n",
            ]
        )

        for atom_id in delete_ids:
            lines.append(
                f"{atom_id}\n"
            )

    lines.extend(
        [
            "\n",
            "Wildcards\n",
            "\n",
        ]
    )

    for atom_id in wildcards:
        lines.append(
            f"{atom_id}\n"
        )

    path.write_text(
        "".join(lines),
        encoding="utf-8",
    )

    return path


# =============================================================================
# Construction / cache isolation
# =============================================================================


def test_constructor_initializes_expected_cache_groups():
    detector = (
        DeduplicationDetector()
    )

    assert detector.seen_reactions == {
        detector.LAMMPS_COMPARISON_GROUP: [],
        detector.RDKIT_COMPARISON_GROUP: [],
    }

    assert (
        detector.seen_reaction_pairs
        == {
            detector.LAMMPS_COMPARISON_GROUP: [],
            detector.RDKIT_COMPARISON_GROUP: [],
        }
    )


def test_detector_instances_have_independent_caches():
    first = DeduplicationDetector()
    second = DeduplicationDetector()

    first.seen_reactions[
        first.RDKIT_COMPARISON_GROUP
    ].append(
        nx.Graph()
    )

    assert (
        second.seen_reactions[
            second.RDKIT_COMPARISON_GROUP
        ]
        == []
    )


# =============================================================================
# clear_cache
# =============================================================================


def test_clear_cache_clears_all_groups():
    detector = (
        DeduplicationDetector()
    )

    detector.seen_reactions[
        "rdkit"
    ].append(
        nx.Graph()
    )

    detector.seen_reaction_pairs[
        "lammps"
    ].append(
        (
            nx.Graph(),
            nx.Graph(),
        )
    )

    detector.clear_cache()

    assert all(
        not value
        for value
        in detector.seen_reactions.values()
    )

    assert all(
        not value
        for value
        in detector.seen_reaction_pairs.values()
    )


def test_clear_cache_specific_group_only():
    detector = (
        DeduplicationDetector()
    )

    detector.seen_reactions[
        "rdkit"
    ].append(
        nx.Graph()
    )

    detector.seen_reactions[
        "lammps"
    ].append(
        nx.Graph()
    )

    detector.clear_cache(
        "rdkit"
    )

    assert (
        detector.seen_reactions[
            "rdkit"
        ]
        == []
    )

    assert len(
        detector.seen_reactions[
            "lammps"
        ]
    ) == 1


def test_clear_cache_creates_unknown_group():
    detector = (
        DeduplicationDetector()
    )

    detector.clear_cache(
        "custom"
    )

    assert (
        detector.seen_reactions[
            "custom"
        ]
        == []
    )

    assert (
        detector.seen_reaction_pairs[
            "custom"
        ]
        == []
    )


# =============================================================================
# is_duplicate_pair
# =============================================================================


def test_duplicate_pair_first_pair_is_unique():
    detector = (
        DeduplicationDetector()
    )

    pre = make_phase_graph(
        detector,
        {
            0: "C",
            1: "O",
        },
        [
            (
                0,
                1,
                "SINGLE",
            )
        ],
    )

    post = make_phase_graph(
        detector,
        {
            0: "C",
            1: "O",
        },
        [
            (
                0,
                1,
                "DOUBLE",
            )
        ],
    )

    assert (
        detector.is_duplicate_pair(
            pre,
            post,
            "test",
        )
        is False
    )


def test_duplicate_pair_second_identical_pair_is_duplicate():
    detector = (
        DeduplicationDetector()
    )

    pre = make_phase_graph(
        detector,
        {
            0: "C",
            1: "O",
        },
        [
            (
                0,
                1,
                "SINGLE",
            )
        ],
    )

    post = make_phase_graph(
        detector,
        {
            0: "C",
            1: "O",
        },
        [
            (
                0,
                1,
                "DOUBLE",
            )
        ],
    )

    assert not detector.is_duplicate_pair(
        pre,
        post,
        "test",
    )

    assert detector.is_duplicate_pair(
        pre,
        post,
        "test",
    )


def test_duplicate_pair_requires_matching_pre_and_post_pair():
    detector = (
        DeduplicationDetector()
    )

    pre = make_phase_graph(
        detector,
        {
            0: "C",
        },
    )

    post_1 = make_phase_graph(
        detector,
        {
            0: "O",
        },
    )

    post_2 = make_phase_graph(
        detector,
        {
            0: "N",
        },
    )

    assert not detector.is_duplicate_pair(
        pre,
        post_1,
        "test",
    )

    assert not detector.is_duplicate_pair(
        pre,
        post_2,
        "test",
    )


def test_duplicate_pair_is_structure_based_not_node_id_based():
    detector = (
        DeduplicationDetector()
    )

    pre_1 = make_phase_graph(
        detector,
        {
            1: "C",
            2: "O",
        },
        [
            (
                1,
                2,
                "SINGLE",
            )
        ],
    )

    post_1 = make_phase_graph(
        detector,
        {
            1: "C",
            2: "O",
        },
        [
            (
                1,
                2,
                "DOUBLE",
            )
        ],
    )

    pre_2 = make_phase_graph(
        detector,
        {
            100: "C",
            200: "O",
        },
        [
            (
                100,
                200,
                "SINGLE",
            )
        ],
    )

    post_2 = make_phase_graph(
        detector,
        {
            100: "C",
            200: "O",
        },
        [
            (
                100,
                200,
                "DOUBLE",
            )
        ],
    )

    assert not detector.is_duplicate_pair(
        pre_1,
        post_1,
        "test",
    )

    assert detector.is_duplicate_pair(
        pre_2,
        post_2,
        "test",
    )


def test_duplicate_pair_cache_uses_copies():
    detector = (
        DeduplicationDetector()
    )

    pre = make_phase_graph(
        detector,
        {
            0: "C",
        },
    )

    post = make_phase_graph(
        detector,
        {
            0: "O",
        },
    )

    detector.is_duplicate_pair(
        pre,
        post,
        "test",
    )

    pre.add_node(
        99,
        atom_label="X",
    )

    cached_pre, _ = (
        detector
        .seen_reaction_pairs[
            "test"
        ][0]
    )

    assert 99 not in cached_pre


# =============================================================================
# Coupled RDKit graph handling
# =============================================================================


def test_couple_graphs_builds_pre_post_nodes_and_correspondence_edges():
    detector = (
        DeduplicationDetector()
    )

    pre = make_phase_graph(
        detector,
        {
            0: "C",
            1: "O",
        },
        [
            (
                0,
                1,
                "SINGLE",
            )
        ],
    )

    post = make_phase_graph(
        detector,
        {
            0: "C",
            1: "O",
        },
        [
            (
                0,
                1,
                "DOUBLE",
            )
        ],
    )

    coupled = detector._couple_graphs(
        pre,
        post,
    )

    assert set(
        coupled.nodes
    ) == {
        ("pre", 0),
        ("pre", 1),
        ("post", 0),
        ("post", 1),
    }

    assert coupled.has_edge(
        ("pre", 0),
        ("post", 0),
    )

    assert (
        coupled.edges[
            ("pre", 0),
            ("post", 0),
        ]["relationship"]
        == "atom_correspondence"
    )


def test_couple_graphs_requires_matching_atom_ids():
    detector = (
        DeduplicationDetector()
    )

    pre = make_phase_graph(
        detector,
        {
            0: "C",
            1: "O",
        },
    )

    post = make_phase_graph(
        detector,
        {
            0: "C",
            2: "O",
        },
    )

    with pytest.raises(
        ValueError,
        match="matching atom IDs",
    ):
        detector._couple_graphs(
            pre,
            post,
        )


def test_couple_graphs_preserves_radical_signature():
    detector = (
        DeduplicationDetector()
    )

    pre = make_phase_graph(
        detector,
        {
            0: "C",
        },
    )

    post = make_phase_graph(
        detector,
        {
            0: "C",
        },
    )

    pre.graph[
        detector.RADICAL_COUNT_ATTRIBUTE
    ] = 1

    pre.graph[
        detector.RADICAL_PRESENT_ATTRIBUTE
    ] = True

    post.graph[
        detector.RADICAL_COUNT_ATTRIBUTE
    ] = 0

    post.graph[
        detector.RADICAL_PRESENT_ATTRIBUTE
    ] = False

    coupled = detector._couple_graphs(
        pre,
        post,
    )

    assert coupled.graph[
        detector.RADICAL_SIGNATURE_ATTRIBUTE
    ] == (
        1,
        0,
        True,
        False,
    )


def test_add_phase_requires_node_label():
    detector = (
        DeduplicationDetector()
    )

    source = nx.Graph()
    source.add_node(0)

    target = nx.Graph()

    with pytest.raises(
        ValueError,
        match="atom_label",
    ):
        detector._add_phase_to_coupled_graph(
            source,
            target,
            "pre",
        )


def test_add_phase_requires_edge_label():
    detector = (
        DeduplicationDetector()
    )

    source = nx.Graph()

    source.add_node(
        0,
        atom_label="C",
    )

    source.add_node(
        1,
        atom_label="O",
    )

    source.add_edge(
        0,
        1,
    )

    target = nx.Graph()

    with pytest.raises(
        ValueError,
        match="bond_label",
    ):
        detector._add_phase_to_coupled_graph(
            source,
            target,
            "pre",
        )


# =============================================================================
# is_duplicate / coupled cache
# =============================================================================


def test_is_duplicate_first_coupled_reaction_is_unique():
    detector = (
        DeduplicationDetector()
    )

    pre = make_phase_graph(
        detector,
        {
            0: "C",
        },
    )

    post = make_phase_graph(
        detector,
        {
            0: "O",
        },
    )

    assert (
        detector.is_duplicate(
            pre,
            post,
            "test",
        )
        is False
    )


def test_is_duplicate_second_equivalent_coupled_reaction_is_duplicate():
    detector = (
        DeduplicationDetector()
    )

    pre = make_phase_graph(
        detector,
        {
            0: "C",
        },
    )

    post = make_phase_graph(
        detector,
        {
            0: "O",
        },
    )

    detector.is_duplicate(
        pre,
        post,
        "test",
    )

    assert detector.is_duplicate(
        pre,
        post,
        "test",
    )


def test_coupled_duplicate_radical_signature_must_match():
    detector = (
        DeduplicationDetector()
    )

    pre_1 = make_phase_graph(
        detector,
        {0: "C"},
    )

    post_1 = make_phase_graph(
        detector,
        {0: "C"},
    )

    pre_1.graph[
        detector.RADICAL_PRESENT_ATTRIBUTE
    ] = False

    post_1.graph[
        detector.RADICAL_PRESENT_ATTRIBUTE
    ] = False

    pre_2 = make_phase_graph(
        detector,
        {0: "C"},
    )

    post_2 = make_phase_graph(
        detector,
        {0: "C"},
    )

    pre_2.graph[
        detector.RADICAL_PRESENT_ATTRIBUTE
    ] = False

    post_2.graph[
        detector.RADICAL_PRESENT_ATTRIBUTE
    ] = True

    assert not detector.is_duplicate(
        pre_1,
        post_1,
        "test",
    )

    assert not detector.is_duplicate(
        pre_2,
        post_2,
        "test",
    )


# =============================================================================
# RDKit graph conversion
# =============================================================================


def test_rdkit_mol_to_networkx_rejects_none():
    detector = (
        DeduplicationDetector()
    )

    with pytest.raises(
        ValueError,
        match="None RDKit molecule",
    ):
        detector.rdkit_mol_to_networkx(
            None
        )


def test_rdkit_mol_to_networkx_builds_all_atoms():
    detector = (
        DeduplicationDetector()
    )

    molecule = mol(
        "CCO"
    )

    graph = (
        detector.rdkit_mol_to_networkx(
            molecule,
            deep_check=False,
        )
    )

    assert graph.number_of_nodes() == 3

    assert graph.number_of_edges() == 2


def test_rdkit_graph_node_label_contains_element():
    detector = (
        DeduplicationDetector()
    )

    graph = (
        detector.rdkit_mol_to_networkx(
            mol("CO"),
            deep_check=False,
        )
    )

    assert graph.nodes[
        0
    ][
        detector.NODE_ATTRIBUTE
    ][0] == "C"

    assert graph.nodes[
        1
    ][
        detector.NODE_ATTRIBUTE
    ][0] == "O"


def test_rdkit_graph_edges_include_bond_type():
    detector = (
        DeduplicationDetector()
    )

    graph = (
        detector.rdkit_mol_to_networkx(
            mol("C=C"),
            deep_check=False,
        )
    )

    assert (
        graph.edges[
            0,
            1,
        ][detector.EDGE_ATTRIBUTE]
        == "DOUBLE"
    )


def test_rdkit_graph_can_be_restricted_to_atom_indices():
    detector = (
        DeduplicationDetector()
    )

    graph = (
        detector.rdkit_mol_to_networkx(
            mol("CCCC"),
            atom_idxs={
                1,
                2,
            },
            deep_check=False,
        )
    )

    assert set(
        graph.nodes
    ) == {
        1,
        2,
    }

    assert graph.number_of_edges() == 1


def test_rdkit_graph_relabels_nodes():
    detector = (
        DeduplicationDetector()
    )

    graph = (
        detector.rdkit_mol_to_networkx(
            mol("CO"),
            atom_idxs={
                0,
                1,
            },
            idx_relabel={
                0: 10,
                1: 20,
            },
            deep_check=False,
        )
    )

    assert set(
        graph.nodes
    ) == {
        10,
        20,
    }

    assert graph.has_edge(
        10,
        20,
    )


def test_rdkit_graph_relabel_requires_all_included_indices():
    detector = (
        DeduplicationDetector()
    )

    with pytest.raises(
        ValueError,
        match="does not contain entries",
    ):
        detector.rdkit_mol_to_networkx(
            mol("CO"),
            atom_idxs={
                0,
                1,
            },
            idx_relabel={
                0: 10,
            },
        )


def test_rdkit_deep_check_includes_external_neighbor_signature():
    detector = (
        DeduplicationDetector()
    )

    molecule = mol(
        "CCO"
    )

    graph = (
        detector.rdkit_mol_to_networkx(
            molecule,
            atom_idxs={
                1,
            },
            deep_check=True,
        )
    )

    label = graph.nodes[
        1
    ][detector.NODE_ATTRIBUTE]

    assert label[0] == "C"

    signature = label[2]

    external_symbols = {
        item[0]
        for item in signature
    }

    assert external_symbols == {
        "C",
        "O",
    }


def test_rdkit_shallow_check_omits_external_environment_signature():
    detector = (
        DeduplicationDetector()
    )

    graph = (
        detector.rdkit_mol_to_networkx(
            mol("CCO"),
            atom_idxs={
                1,
            },
            deep_check=False,
        )
    )

    label = graph.nodes[
        1
    ][detector.NODE_ATTRIBUTE]

    assert len(label) == 2


# =============================================================================
# Radical helpers
# =============================================================================


def test_explicit_radical_atom_is_detected():
    detector = (
        DeduplicationDetector()
    )

    radical = mol(
        "[CH3]"
    )

    atom = radical.GetAtomWithIdx(
        0
    )

    assert (
        detector._is_radical_atom(
            atom
        )
        is True
    )


def test_normal_carbon_is_not_radical():
    detector = (
        DeduplicationDetector()
    )

    atom = (
        mol("C")
        .GetAtomWithIdx(0)
    )

    assert (
        detector._is_radical_atom(
            atom
        )
        is False
    )


def test_count_radical_atoms():
    detector = (
        DeduplicationDetector()
    )

    molecule = mol(
        "[CH3].[CH3]"
    )

    assert (
        detector._count_radical_atoms(
            molecule
        )
        == 2
    )


# =============================================================================
# Mapping selection
# =============================================================================


def test_select_template_mapping():
    detector = (
        DeduplicationDetector()
    )

    metadata = make_metadata(
        template_mapping={
            0: 1,
        },
    )

    result = (
        detector
        ._select_reactant_to_product_mapping(
            metadata,
            reaction_index=1,
            index_source="template",
        )
    )

    assert result == {
        0: 1,
    }


def test_select_template_mapping_requires_mapping():
    detector = (
        DeduplicationDetector()
    )

    metadata = make_metadata()

    metadata.template_reactant_to_product_mapping = None

    with pytest.raises(
        ValueError,
        match="template_reactant_to_product_mapping",
    ):
        detector._select_reactant_to_product_mapping(
            metadata,
            1,
            "template",
        )


def test_select_first_shell_mapping():
    detector = (
        DeduplicationDetector()
    )

    metadata = make_metadata(
        full_mapping={
            0: 2,
            1: 1,
            2: 0,
        },
        first_shell=[
            0,
            2,
        ],
    )

    result = (
        detector
        ._select_reactant_to_product_mapping(
            metadata,
            1,
            "first_shell",
        )
    )

    assert result == {
        0: 2,
        2: 0,
    }


def test_select_first_shell_requires_indices():
    detector = (
        DeduplicationDetector()
    )

    metadata = make_metadata()

    metadata.first_shell = None

    with pytest.raises(
        ValueError,
        match="first_shell",
    ):
        detector._select_reactant_to_product_mapping(
            metadata,
            1,
            "first_shell",
        )


def test_select_first_shell_requires_full_mapping():
    detector = (
        DeduplicationDetector()
    )

    metadata = make_metadata()

    metadata.first_shell = [
        0,
    ]

    metadata.reactant_to_product_mapping = None

    with pytest.raises(
        ValueError,
        match="reactant_to_product_mapping",
    ):
        detector._select_reactant_to_product_mapping(
            metadata,
            1,
            "first_shell",
        )


def test_select_mapping_rejects_unknown_source():
    detector = (
        DeduplicationDetector()
    )

    metadata = make_metadata()

    with pytest.raises(
        ValueError,
        match="Unsupported index_source",
    ):
        detector._select_reactant_to_product_mapping(
            metadata,
            1,
            "wrong",
        )


# =============================================================================
# compare_graphs_mol
# =============================================================================


def test_compare_graphs_mol_keeps_first_unique_reaction():
    detector = (
        DeduplicationDetector()
    )

    metadata = make_metadata()

    result = (
        detector.compare_graphs_mol(
            [
                metadata,
            ]
        )
    )

    assert result == [
        metadata,
    ]

    assert (
        metadata.activity_stats
        is True
    )


def test_compare_graphs_mol_disables_later_duplicate():
    detector = (
        DeduplicationDetector()
    )

    first = make_metadata(
        reaction_id=1,
    )

    second = make_metadata(
        reaction_id=2,
    )

    result = (
        detector.compare_graphs_mol(
            [
                first,
                second,
            ]
        )
    )

    assert result == [
        first,
    ]

    assert (
        first.activity_stats
        is True
    )

    assert (
        second.activity_stats
        is False
    )


def test_compare_graphs_mol_repeated_same_object_is_not_disabled():
    detector = (
        DeduplicationDetector()
    )

    metadata = make_metadata()

    result = (
        detector.compare_graphs_mol(
            [
                metadata,
                metadata,
            ]
        )
    )

    assert result == [
        metadata,
    ]

    assert (
        metadata.activity_stats
        is True
    )


def test_compare_graphs_mol_skips_already_inactive_reactions():
    detector = (
        DeduplicationDetector()
    )

    metadata = make_metadata(
        activity_stats=False
    )

    result = (
        detector.compare_graphs_mol(
            [
                metadata,
            ]
        )
    )

    assert result == []


def test_compare_graphs_mol_requires_reactant_molecule():
    detector = (
        DeduplicationDetector()
    )

    metadata = make_metadata()

    metadata.reactant_combined_RDmol = None

    with pytest.raises(
        ValueError,
        match="combined reactant",
    ):
        detector.compare_graphs_mol(
            [
                metadata,
            ]
        )


def test_compare_graphs_mol_requires_product_molecule():
    detector = (
        DeduplicationDetector()
    )

    metadata = make_metadata()

    metadata.product_combined_RDmol = None

    with pytest.raises(
        ValueError,
        match="combined product",
    ):
        detector.compare_graphs_mol(
            [
                metadata,
            ]
        )


def test_compare_graphs_mol_rejects_nonbijective_mapping():
    detector = (
        DeduplicationDetector()
    )

    metadata = make_metadata(
        template_mapping={
            0: 0,
            1: 0,
        },
    )

    with pytest.raises(
        ValueError,
        match="non-bijective",
    ):
        detector.compare_graphs_mol(
            [
                metadata,
            ]
        )


def test_compare_graphs_mol_supports_first_shell_source():
    detector = (
        DeduplicationDetector()
    )

    metadata = make_metadata(
        template_mapping=None,
        full_mapping={
            0: 0,
            1: 1,
        },
        first_shell=[
            0,
        ],
    )

    metadata.template_reactant_to_product_mapping = None

    result = detector.compare_graphs_mol(
        [
            metadata,
        ],
        index_source="first_shell",
    )

    assert result == [
        metadata,
    ]


def test_compare_graphs_mol_clears_rdkit_cache_each_pass():
    detector = (
        DeduplicationDetector()
    )

    first = make_metadata(
        reaction_id=1,
    )

    second = make_metadata(
        reaction_id=2,
    )

    assert detector.compare_graphs_mol(
        [
            first,
        ]
    ) == [
        first,
    ]

    # A fresh deduplication pass must not treat an equivalent reaction
    # as a duplicate merely because it appeared in a previous call.
    assert detector.compare_graphs_mol(
        [
            second,
        ]
    ) == [
        second,
    ]


def test_compare_graphs_mol_radical_metadata_changes_signature():
    detector = (
        DeduplicationDetector()
    )

    non_radical = make_metadata(
        reaction_id=1,
        is_radical=False,
    )

    radical = make_metadata(
        reaction_id=2,
        is_radical=True,
    )

    result = detector.compare_graphs_mol(
        [
            non_radical,
            radical,
        ]
    )

    assert result == [
        non_radical,
        radical,
    ]


def test_compare_graphs_mol_forwards_deep_check(
    monkeypatch,
):
    detector = (
        DeduplicationDetector()
    )

    metadata = make_metadata()

    calls = []

    original = (
        detector.rdkit_mol_to_networkx
    )

    def wrapper(
        molecule,
        atom_idxs=None,
        idx_relabel=None,
        deep_check=True,
    ):
        calls.append(
            deep_check
        )

        return original(
            molecule=molecule,
            atom_idxs=atom_idxs,
            idx_relabel=idx_relabel,
            deep_check=deep_check,
        )

    monkeypatch.setattr(
        detector,
        "rdkit_mol_to_networkx",
        wrapper,
    )

    detector.compare_graphs_mol(
        [
            metadata,
        ],
        deep_check=False,
    )

    assert calls == [
        False,
        False,
    ]


# =============================================================================
# LAMMPS template filename helpers
# =============================================================================


@pytest.mark.parametrize(
    "name, expected",
    [
        (
            "template_pre_1.molecule",
            "1",
        ),
        (
            "template_post_22.molecule",
            "22",
        ),
        (
            "template_pre_1_homo2.molecule",
            "1_homo2",
        ),
    ],
)
def test_lammps_reaction_id_from_template_path(
    name,
    expected,
):
    result = (
        DeduplicationDetector
        ._lammps_reaction_id_from_template_path(
            name
        )
    )

    assert result == expected


@pytest.mark.parametrize(
    "name",
    [
        "RXN_1.map",
        "pre_1.molecule",
        "template_1.molecule",
        "template_pre_1.txt",
    ],
)
def test_lammps_reaction_id_rejects_bad_filename(
    name,
):
    with pytest.raises(
        ValueError,
        match="Could not infer reaction ID",
    ):
        (
            DeduplicationDetector
            ._lammps_reaction_id_from_template_path(
                name
            )
        )


def test_lammps_map_path_from_template_path(
    tmp_path,
):
    detector = (
        DeduplicationDetector()
    )

    template = (
        tmp_path
        / "template_pre_42.molecule"
    )

    assert (
        detector
        ._lammps_map_path_from_template_path(
            template
        )
        == tmp_path / "RXN_42.map"
    )


# =============================================================================
# LAMMPS map parsing
# =============================================================================


def test_is_lammps_map_section_header():
    detector = (
        DeduplicationDetector
    )

    assert detector._is_lammps_map_section_header(
        "EdgeIDs"
    )

    assert detector._is_lammps_map_section_header(
        "Equivalences extra"
    )

    assert not detector._is_lammps_map_section_header(
        "1 2"
    )


@pytest.mark.parametrize(
    "line",
    [
        "11 edgeIDs",
        "2 equivalences",
        "0 wildcards",
        "  5 EDGEIDS  ",
    ],
)
def test_is_lammps_count_line_true(
    line,
):
    assert (
        DeduplicationDetector
        ._is_lammps_count_line(
            line
        )
    )


def test_is_lammps_count_line_false():
    assert not (
        DeduplicationDetector
        ._is_lammps_count_line(
            "2 deleteIDs"
        )
    )


def test_read_lammps_equivalences(
    tmp_path,
):
    path = write_map_file(
        tmp_path / "RXN_1.map",
        equivalences=[
            (
                1,
                10,
            ),
            (
                2,
                20,
            ),
        ],
    )

    result = (
        DeduplicationDetector
        ._read_lammps_equivalences(
            path
        )
    )

    assert result == {
        1: 10,
        2: 20,
    }


def test_read_lammps_equivalences_rejects_conflicting_pre_mapping(
    tmp_path,
):
    path = write_map_file(
        tmp_path / "RXN_1.map",
        equivalences=[
            (
                1,
                10,
            ),
            (
                1,
                20,
            ),
        ],
    )

    with pytest.raises(
        ValueError,
        match="Conflicting Equivalences",
    ):
        (
            DeduplicationDetector
            ._read_lammps_equivalences(
                path
            )
        )


def test_read_lammps_equivalences_rejects_nonbijective_post_mapping(
    tmp_path,
):
    path = write_map_file(
        tmp_path / "RXN_1.map",
        equivalences=[
            (
                1,
                10,
            ),
            (
                2,
                10,
            ),
        ],
    )

    with pytest.raises(
        ValueError,
        match="Non-bijective",
    ):
        (
            DeduplicationDetector
            ._read_lammps_equivalences(
                path
            )
        )


def test_read_lammps_edge_ids(
    tmp_path,
):
    path = write_map_file(
        tmp_path / "RXN_1.map",
        edge_ids=[
            5,
            8,
            11,
        ],
    )

    assert (
        DeduplicationDetector
        ._read_lammps_edge_ids(
            path
        )
        == [
            5,
            8,
            11,
        ]
    )


def test_read_lammps_single_column_section(
    tmp_path,
):
    path = write_map_file(
        tmp_path / "RXN_1.map",
        initiator_ids=[
            3,
            7,
        ],
    )

    assert (
        DeduplicationDetector
        ._read_lammps_single_column_section(
            path,
            "InitiatorIDs",
        )
        == [
            3,
            7,
        ]
    )


def test_read_lammps_integer_section_ignores_non_integer_rows(
    tmp_path,
):
    path = tmp_path / "RXN_1.map"

    path.write_text(
        """
Equivalences

1 10
bad row
2 20
""",
        encoding="utf-8",
    )

    rows = (
        DeduplicationDetector
        ._read_lammps_map_integer_section(
            path,
            "Equivalences",
        )
    )

    assert rows == [
        [
            1,
            10,
        ],
        [
            2,
            20,
        ],
    ]


# =============================================================================
# Map text manipulation
# =============================================================================


def test_remove_lammps_map_section():
    lines = [
        "# header\n",
        "EdgeIDs\n",
        "\n",
        "1\n",
        "2\n",
        "\n",
        "Equivalences\n",
        "\n",
        "1 1\n",
    ]

    result = (
        DeduplicationDetector
        ._remove_lammps_map_section(
            lines,
            "EdgeIDs",
        )
    )

    joined = "".join(result)

    assert "EdgeIDs" not in joined
    assert "Equivalences" in joined
    assert "1 1" in joined


def test_replace_lammps_count_lines():
    detector = (
        DeduplicationDetector()
    )

    lines = [
        "# header\n",
        "\n",
        "2 edgeIDs\n",
        "5 equivalences\n",
        "0 wildcards\n",
        "\n",
        "InitiatorIDs\n",
    ]

    result = (
        detector
        ._replace_lammps_count_lines(
            lines,
            edge_count=4,
            equivalence_count=8,
            wildcard_count=2,
        )
    )

    joined = "".join(result)

    assert "4 edgeIDs" in joined
    assert "8 equivalences" in joined
    assert "2 wildcards" in joined

    assert "5 equivalences" not in joined


# =============================================================================
# Wildcard-map writing
# =============================================================================


def test_write_lammps_wildcard_map_file(
    tmp_path,
):
    detector = (
        DeduplicationDetector()
    )

    path = write_map_file(
        tmp_path / "RXN_1.map",
        equivalences=[
            (
                1,
                10,
            ),
            (
                2,
                20,
            ),
        ],
        edge_ids=[
            1,
            2,
        ],
        initiator_ids=[
            1,
            2,
        ],
        delete_ids=[
            3,
        ],
    )

    detector._write_lammps_wildcard_map_file(
        path,
        wildcard_ids=[
            2,
            1,
            2,
        ],
    )

    text = path.read_text(
        encoding="utf-8"
    )

    assert "2 edgeIDs" in text
    assert "2 equivalences" in text
    assert "1 deleteIDs" in text
    assert "2 wildcards" in text

    assert "Wildcards" in text

    wildcards = (
        detector
        ._read_lammps_single_column_section(
            path,
            "Wildcards",
        )
    )

    assert wildcards == [
        2,
        1,
    ]


# =============================================================================
# Wildcard graph handling
# =============================================================================


def test_remove_nodes_if_present_returns_copy():
    graph = nx.Graph()

    graph.add_edges_from(
        [
            (
                1,
                2,
            ),
            (
                2,
                3,
            ),
        ]
    )

    result = (
        DeduplicationDetector
        ._remove_nodes_if_present(
            graph,
            [
                2,
                999,
            ],
        )
    )

    assert 2 not in result

    assert 2 in graph


def test_apply_lammps_wildcards_removes_pre_and_equivalent_post_nodes(
    tmp_path,
):
    detector = (
        DeduplicationDetector()
    )

    pre = make_phase_graph(
        detector,
        {
            1: "A",
            2: "B",
        },
        [
            (
                1,
                2,
                "1",
            )
        ],
    )

    post = make_phase_graph(
        detector,
        {
            10: "A",
            20: "B",
        },
        [
            (
                10,
                20,
                "1",
            )
        ],
    )

    map_path = write_map_file(
        tmp_path / "RXN_1.map",
        equivalences=[
            (
                1,
                10,
            ),
            (
                2,
                20,
            ),
        ],
        edge_ids=[
            2,
        ],
    )

    new_pre, new_post, mapping = (
        detector._apply_lammps_wildcards(
            pre,
            post,
            {
                1: 10,
                2: 20,
            },
            map_path,
        )
    )

    assert set(
        new_pre.nodes
    ) == {
        1,
    }

    assert set(
        new_post.nodes
    ) == {
        10,
    }

    assert mapping == {
        1: 10,
    }


# =============================================================================
# LAMMPS graph coupling
# =============================================================================


def test_couple_lammps_graphs_uses_explicit_equivalences():
    detector = (
        DeduplicationDetector()
    )

    pre = make_phase_graph(
        detector,
        {
            1: "A",
            2: "B",
        },
    )

    post = make_phase_graph(
        detector,
        {
            10: "A",
            20: "B",
        },
    )

    coupled = (
        detector
        ._couple_lammps_graphs(
            pre,
            post,
            {
                1: 20,
                2: 10,
            },
            "test.map",
        )
    )

    assert coupled.has_edge(
        ("pre", 1),
        ("post", 20),
    )

    assert coupled.has_edge(
        ("pre", 2),
        ("post", 10),
    )


def test_couple_lammps_graphs_rejects_missing_pre_atom():
    detector = (
        DeduplicationDetector()
    )

    pre = make_phase_graph(
        detector,
        {
            1: "A",
        },
    )

    post = make_phase_graph(
        detector,
        {
            10: "A",
        },
    )

    with pytest.raises(
        ValueError,
        match="references pre atom",
    ):
        detector._couple_lammps_graphs(
            pre,
            post,
            {
                2: 10,
            },
            "test.map",
        )


def test_couple_lammps_graphs_rejects_missing_post_atom():
    detector = (
        DeduplicationDetector()
    )

    pre = make_phase_graph(
        detector,
        {
            1: "A",
        },
    )

    post = make_phase_graph(
        detector,
        {
            10: "A",
        },
    )

    with pytest.raises(
        ValueError,
        match="references post atom",
    ):
        detector._couple_lammps_graphs(
            pre,
            post,
            {
                1: 20,
            },
            "test.map",
        )


# =============================================================================
# LAMMPS molecule parsing
# =============================================================================


def test_lammps_molecule_to_networkx_rejects_missing_file(
    tmp_path,
):
    detector = (
        DeduplicationDetector()
    )

    with pytest.raises(
        FileNotFoundError
    ):
        detector.lammps_molecule_to_networkx(
            tmp_path
            / "missing.molecule"
        )


def test_lammps_molecule_to_networkx_requires_types_section(
    tmp_path,
):
    detector = (
        DeduplicationDetector()
    )

    path = (
        tmp_path / "bad.molecule"
    )

    path.write_text(
        """
Coords

1 0 0 0
""",
        encoding="utf-8",
    )

    with pytest.raises(
        ValueError,
        match="Types section",
    ):
        detector.lammps_molecule_to_networkx(
            path
        )


def test_lammps_molecule_to_networkx_builds_nodes_and_bonds(
    tmp_path,
):
    detector = (
        DeduplicationDetector()
    )

    path = write_lammps_molecule(
        tmp_path
        / "template_pre_1.molecule",
        {
            1: "3",
            2: "7",
        },
        [
            (
                1,
                "2",
                1,
                2,
            )
        ],
    )

    graph = (
        detector
        .lammps_molecule_to_networkx(
            path
        )
    )

    assert graph.nodes[
        1
    ][detector.NODE_ATTRIBUTE] == "3"

    assert graph.nodes[
        2
    ][detector.NODE_ATTRIBUTE] == "7"

    assert graph.edges[
        1,
        2,
    ][detector.EDGE_ATTRIBUTE] == "2"


def test_lammps_graph_ignores_coordinates_and_charges(
    tmp_path,
):
    detector = (
        DeduplicationDetector()
    )

    path = write_lammps_molecule(
        tmp_path
        / "template_pre_1.molecule",
        {
            1: "5",
        },
    )

    graph = (
        detector
        .lammps_molecule_to_networkx(
            path
        )
    )

    assert graph.number_of_nodes() == 1

    assert set(
        graph.nodes[1]
    ) == {
        detector.NODE_ATTRIBUTE,
    }


def test_add_lammps_atoms_rejects_bad_line(
    tmp_path,
):
    detector = (
        DeduplicationDetector()
    )

    graph = nx.Graph()

    with pytest.raises(
        ValueError,
        match="Invalid Types line",
    ):
        detector._add_lammps_atoms(
            graph,
            [
                "1",
            ],
            tmp_path
            / "test.molecule",
        )


def test_add_lammps_atoms_rejects_invalid_id(
    tmp_path,
):
    detector = (
        DeduplicationDetector()
    )

    with pytest.raises(
        ValueError,
        match="Invalid atom ID",
    ):
        detector._add_lammps_atoms(
            nx.Graph(),
            [
                "abc 4",
            ],
            tmp_path
            / "test.molecule",
        )


def test_add_lammps_atoms_rejects_duplicate_atom_id(
    tmp_path,
):
    detector = (
        DeduplicationDetector()
    )

    with pytest.raises(
        ValueError,
        match="Duplicate atom ID",
    ):
        detector._add_lammps_atoms(
            nx.Graph(),
            [
                "1 4",
                "1 5",
            ],
            tmp_path
            / "test.molecule",
        )


def test_add_lammps_bonds_rejects_short_line(
    tmp_path,
):
    detector = (
        DeduplicationDetector()
    )

    graph = nx.Graph()

    graph.add_node(
        1,
        atom_label="1",
    )

    with pytest.raises(
        ValueError,
        match="Invalid Bonds line",
    ):
        detector._add_lammps_bonds(
            graph,
            [
                "1 2 1",
            ],
            tmp_path
            / "test.molecule",
        )


def test_add_lammps_bonds_rejects_undefined_atom(
    tmp_path,
):
    detector = (
        DeduplicationDetector()
    )

    graph = nx.Graph()

    graph.add_node(
        1,
        atom_label="1",
    )

    with pytest.raises(
        ValueError,
        match="undefined atom IDs",
    ):
        detector._add_lammps_bonds(
            graph,
            [
                "1 1 1 2",
            ],
            tmp_path
            / "test.molecule",
        )


# =============================================================================
# Complete LAMMPS template deduplication
# =============================================================================


def test_lammps_template_pair_without_map_uses_uncoupled_pair_cache(
    tmp_path,
):
    detector = (
        DeduplicationDetector()
    )

    pre = write_lammps_molecule(
        tmp_path
        / "template_pre_1.molecule",
        {
            1: "1",
            2: "2",
        },
        [
            (
                1,
                "1",
                1,
                2,
            )
        ],
    )

    post = write_lammps_molecule(
        tmp_path
        / "template_post_1.molecule",
        {
            1: "1",
            2: "2",
        },
        [
            (
                1,
                "2",
                1,
                2,
            )
        ],
    )

    assert not (
        detector
        .is_duplicate_lammps_template_pair(
            pre,
            post,
        )
    )

    assert (
        detector
        .is_duplicate_lammps_template_pair(
            pre,
            post,
        )
    )


def test_lammps_template_pair_with_map_uses_coupled_cache(
    tmp_path,
):
    detector = (
        DeduplicationDetector()
    )

    pre = write_lammps_molecule(
        tmp_path
        / "template_pre_1.molecule",
        {
            1: "1",
            2: "2",
        },
        [
            (
                1,
                "1",
                1,
                2,
            )
        ],
    )

    post = write_lammps_molecule(
        tmp_path
        / "template_post_1.molecule",
        {
            10: "1",
            20: "2",
        },
        [
            (
                1,
                "2",
                10,
                20,
            )
        ],
    )

    map_path = write_map_file(
        tmp_path
        / "RXN_1.map",
        equivalences=[
            (
                1,
                10,
            ),
            (
                2,
                20,
            ),
        ],
    )

    assert not (
        detector
        .is_duplicate_lammps_template_pair(
            pre,
            post,
            map_file_path=map_path,
        )
    )

    assert (
        detector
        .is_duplicate_lammps_template_pair(
            pre,
            post,
            map_file_path=map_path,
        )
    )


def test_lammps_template_pair_existing_map_requires_equivalences(
    tmp_path,
):
    detector = (
        DeduplicationDetector()
    )

    pre = write_lammps_molecule(
        tmp_path
        / "template_pre_1.molecule",
        {
            1: "1",
        },
    )

    post = write_lammps_molecule(
        tmp_path
        / "template_post_1.molecule",
        {
            1: "1",
        },
    )

    map_path = (
        tmp_path / "RXN_1.map"
    )

    map_path.write_text(
        """
EdgeIDs

1
""",
        encoding="utf-8",
    )

    with pytest.raises(
        ValueError,
        match="No Equivalences",
    ):
        (
            detector
            .is_duplicate_lammps_template_pair(
                pre,
                post,
                map_file_path=map_path,
            )
        )


# =============================================================================
# compare_graphs
# =============================================================================


def test_compare_graphs_skips_non_pre_files(
    tmp_path,
):
    detector = (
        DeduplicationDetector()
    )

    path = (
        tmp_path
        / "template_post_1.molecule"
    )

    path.write_text(
        "",
        encoding="utf-8",
    )

    assert (
        detector.compare_graphs(
            [
                path,
            ]
        )
        == {}
    )


def test_compare_graphs_skips_missing_post_file(
    tmp_path,
    capsys,
):
    detector = (
        DeduplicationDetector()
    )

    pre = (
        tmp_path
        / "template_pre_1.molecule"
    )

    pre.write_text(
        "",
        encoding="utf-8",
    )

    result = (
        detector.compare_graphs(
            [
                pre,
            ]
        )
    )

    assert result == {}

    assert (
        "does not exist"
        in capsys.readouterr().out
    )


def test_compare_graphs_records_result(
    tmp_path,
    monkeypatch,
):
    detector = (
        DeduplicationDetector()
    )

    pre = (
        tmp_path
        / "template_pre_1.molecule"
    )

    post = (
        tmp_path
        / "template_post_1.molecule"
    )

    pre.write_text(
        "",
        encoding="utf-8",
    )

    post.write_text(
        "",
        encoding="utf-8",
    )

    monkeypatch.setattr(
        detector,
        "is_duplicate_lammps_template_pair",
        lambda **kwargs: True,
    )

    result = detector.compare_graphs(
        [
            pre,
        ]
    )

    assert result == {
        str(pre): True,
    }


# =============================================================================
# write_wildcard_maps
# =============================================================================


def test_write_wildcard_maps_skips_inactive():
    detector = (
        DeduplicationDetector()
    )

    metadata = make_metadata(
        activity_stats=False
    )

    assert (
        detector.write_wildcard_maps(
            [
                metadata,
            ]
        )
        == []
    )


def test_write_wildcard_maps_disables_missing_map_definition():
    detector = (
        DeduplicationDetector()
    )

    metadata = make_metadata()

    metadata.map_file = None

    result = (
        detector.write_wildcard_maps(
            [
                metadata,
            ]
        )
    )

    assert result == []

    assert (
        metadata.activity_stats
        is False
    )


def test_write_wildcard_maps_disables_missing_map_file(
    tmp_path,
):
    detector = (
        DeduplicationDetector()
    )

    metadata = make_metadata()

    metadata.map_file = (
        tmp_path
        / "missing.map"
    )

    result = (
        detector.write_wildcard_maps(
            [
                metadata,
            ]
        )
    )

    assert result == []

    assert (
        metadata.activity_stats
        is False
    )


def test_write_wildcard_maps_writes_active_map(
    tmp_path,
):
    detector = (
        DeduplicationDetector()
    )

    map_path = write_map_file(
        tmp_path
        / "RXN_1.map",
        equivalences=[
            (
                1,
                1,
            )
        ],
        edge_ids=[
            1,
        ],
    )

    metadata = make_metadata()

    metadata.map_file = map_path

    result = (
        detector.write_wildcard_maps(
            [
                metadata,
            ]
        )
    )

    assert result == [
        metadata,
    ]

    wildcards = (
        detector
        ._read_lammps_single_column_section(
            map_path,
            "Wildcards",
        )
    )

    assert wildcards == [
        1,
    ]


# =============================================================================
# compare_lammps_templates orchestration
# =============================================================================


def test_compare_lammps_templates_skips_inactive():
    detector = (
        DeduplicationDetector()
    )

    metadata = make_metadata(
        activity_stats=False
    )

    assert (
        detector.compare_lammps_templates(
            [
                metadata,
            ]
        )
        == []
    )


def test_compare_lammps_templates_disables_missing_file_definitions():
    detector = (
        DeduplicationDetector()
    )

    metadata = make_metadata()

    result = (
        detector.compare_lammps_templates(
            [
                metadata,
            ]
        )
    )

    assert result == []

    assert (
        metadata.activity_stats
        is False
    )


def test_compare_lammps_templates_disables_missing_pre_file(
    tmp_path,
):
    detector = (
        DeduplicationDetector()
    )

    metadata = make_metadata()

    metadata.pre_reaction_file = (
        tmp_path
        / "missing_pre.molecule"
    )

    post = (
        tmp_path
        / "template_post_1.molecule"
    )

    post.write_text(
        "",
        encoding="utf-8",
    )

    metadata.post_reaction_file = post

    result = (
        detector.compare_lammps_templates(
            [
                metadata,
            ]
        )
    )

    assert result == []

    assert (
        metadata.activity_stats
        is False
    )


def test_compare_lammps_templates_disables_duplicate(
    tmp_path,
    monkeypatch,
):
    detector = (
        DeduplicationDetector()
    )

    pre = (
        tmp_path
        / "template_pre_1.molecule"
    )

    post = (
        tmp_path
        / "template_post_1.molecule"
    )

    pre.write_text(
        "",
        encoding="utf-8",
    )

    post.write_text(
        "",
        encoding="utf-8",
    )

    metadata = make_metadata()

    metadata.pre_reaction_file = pre
    metadata.post_reaction_file = post

    monkeypatch.setattr(
        detector,
        "is_duplicate_lammps_template_pair",
        lambda **kwargs: True,
    )

    result = (
        detector.compare_lammps_templates(
            [
                metadata,
            ]
        )
    )

    assert result == []

    assert (
        metadata.activity_stats
        is False
    )


def test_compare_lammps_templates_retains_unique(
    tmp_path,
    monkeypatch,
):
    detector = (
        DeduplicationDetector()
    )

    pre = (
        tmp_path
        / "template_pre_1.molecule"
    )

    post = (
        tmp_path
        / "template_post_1.molecule"
    )

    pre.write_text(
        "",
        encoding="utf-8",
    )

    post.write_text(
        "",
        encoding="utf-8",
    )

    metadata = make_metadata()

    metadata.pre_reaction_file = pre
    metadata.post_reaction_file = post

    monkeypatch.setattr(
        detector,
        "is_duplicate_lammps_template_pair",
        lambda **kwargs: False,
    )

    result = (
        detector.compare_lammps_templates(
            [
                metadata,
            ]
        )
    )

    assert result == [
        metadata,
    ]

    assert (
        metadata.activity_stats
        is True
    )
