from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from rdkit import Chem

from AutoREACTER.reaction_preparation.reaction_processor.utils import (
    add_column_safe,
    add_dict_as_new_columns,
    compare_rdkit_molecules_canonical,
    compare_set,
    extract_unique_references,
    prep_for_3d_molecule_generation,
    prepare_paths,
)


# =============================================================================
# prepare_paths
# =============================================================================


def test_prepare_paths_creates_directory(tmp_path):
    cache = tmp_path / "cache"

    result = prepare_paths(
        cache,
        "csv_cache",
    )

    assert result == (
        cache / "csv_cache"
    )

    assert result.exists()
    assert result.is_dir()


def test_prepare_paths_creates_parent_directories(tmp_path):
    cache = (
        tmp_path
        / "a"
        / "b"
        / "c"
    )

    result = prepare_paths(
        cache,
        "csv_cache",
    )

    assert result.is_dir()


def test_prepare_paths_accepts_string_cache_path(tmp_path):
    cache = tmp_path / "cache"

    result = prepare_paths(
        str(cache),
        "csv_cache",
    )

    assert isinstance(
        result,
        Path,
    )

    assert result == (
        cache / "csv_cache"
    )


def test_prepare_paths_accepts_nested_subdirectory(tmp_path):
    cache = tmp_path / "cache"

    result = prepare_paths(
        cache,
        "one/two/three",
    )

    assert result == (
        cache
        / "one"
        / "two"
        / "three"
    )

    assert result.is_dir()


def test_prepare_paths_is_idempotent(tmp_path):
    cache = tmp_path / "cache"

    first = prepare_paths(
        cache,
        "csv_cache",
    )

    marker = first / "keep.txt"
    marker.write_text(
        "keep",
        encoding="utf-8",
    )

    second = prepare_paths(
        cache,
        "csv_cache",
    )

    assert first == second
    assert marker.exists()

    assert (
        marker.read_text(
            encoding="utf-8"
        )
        == "keep"
    )


# =============================================================================
# add_dict_as_new_columns
# =============================================================================


def test_add_dict_as_new_columns_adds_default_columns():
    df = pd.DataFrame(
        {
            "reactant_idx": [0, 1, 2],
        }
    )

    mapping = {
        0: 2,
        1: 1,
        2: 0,
    }

    result = add_dict_as_new_columns(
        df,
        mapping,
    )

    assert result is df

    assert (
        result[
            "template_reactant_idx"
        ].tolist()
        == [0, 1, 2]
    )

    assert (
        result[
            "template_product_idx"
        ].tolist()
        == [2, 1, 0]
    )


def test_add_dict_as_new_columns_uses_custom_titles():
    df = pd.DataFrame(
        {"base": [10, 20]}
    )

    result = add_dict_as_new_columns(
        df,
        {
            5: 7,
            6: 8,
        },
        titles=[
            "left_idx",
            "right_idx",
        ],
    )

    assert (
        result["left_idx"].tolist()
        == [5, 6]
    )

    assert (
        result["right_idx"].tolist()
        == [7, 8]
    )


def test_add_dict_as_new_columns_uses_nullable_integer_dtype():
    df = pd.DataFrame(
        {"base": [10, 20, 30]}
    )

    result = add_dict_as_new_columns(
        df,
        {
            1: 4,
            2: 5,
        },
    )

    assert str(
        result[
            "template_reactant_idx"
        ].dtype
    ) == "Int64"

    assert str(
        result[
            "template_product_idx"
        ].dtype
    ) == "Int64"


def test_add_dict_as_new_columns_short_mapping_fills_remaining_rows_with_na():
    df = pd.DataFrame(
        {
            "base": [
                "a",
                "b",
                "c",
            ]
        }
    )

    result = add_dict_as_new_columns(
        df,
        {
            10: 20,
        },
    )

    values = result[
        "template_reactant_idx"
    ].tolist()

    assert values[0] == 10
    assert pd.isna(values[1])
    assert pd.isna(values[2])


def test_add_dict_as_new_columns_empty_mapping_creates_nullable_columns():
    df = pd.DataFrame(
        {
            "base": [1, 2],
        }
    )

    result = add_dict_as_new_columns(
        df,
        {},
    )

    assert (
        "template_reactant_idx"
        in result.columns
    )

    assert (
        "template_product_idx"
        in result.columns
    )

    assert result[
        "template_reactant_idx"
    ].isna().all()

    assert result[
        "template_product_idx"
    ].isna().all()


def test_add_dict_as_new_columns_preserves_existing_columns():
    df = pd.DataFrame(
        {
            "reactant_idx": [0, 1],
            "product_idx": [1, 0],
        }
    )

    add_dict_as_new_columns(
        df,
        {
            0: 1,
        },
    )

    assert (
        df["reactant_idx"].tolist()
        == [0, 1]
    )

    assert (
        df["product_idx"].tolist()
        == [1, 0]
    )


def test_add_dict_as_new_columns_is_position_safe_with_nondefault_dataframe_index():
    """
    Utility columns describe row-position data, not pandas index labels.

    A dataframe with a non-default index must therefore receive the values
    in row order rather than turning them into NA through index alignment.
    """
    df = pd.DataFrame(
        {
            "base": ["a", "b"],
        },
        index=[10, 20],
    )

    result = add_dict_as_new_columns(
        df,
        {
            4: 8,
            5: 9,
        },
    )

    assert (
        result[
            "template_reactant_idx"
        ].tolist()
        == [4, 5]
    )

    assert (
        result[
            "template_product_idx"
        ].tolist()
        == [8, 9]
    )


# =============================================================================
# add_column_safe
# =============================================================================


def test_add_column_safe_adds_column():
    df = pd.DataFrame(
        {
            "base": [1, 2, 3],
        }
    )

    result = add_column_safe(
        df,
        [10, 11, 12],
        "edge_atoms",
    )

    assert result is df

    assert (
        result["edge_atoms"].tolist()
        == [10, 11, 12]
    )


def test_add_column_safe_uses_nullable_integer_dtype():
    df = pd.DataFrame(
        {
            "base": [1, 2, 3],
        }
    )

    result = add_column_safe(
        df,
        [10],
        "edge_atoms",
    )

    assert (
        str(
            result[
                "edge_atoms"
            ].dtype
        )
        == "Int64"
    )


def test_add_column_safe_short_list_fills_remaining_rows_with_na():
    df = pd.DataFrame(
        {
            "base": [1, 2, 3],
        }
    )

    result = add_column_safe(
        df,
        [7],
        "edge_atoms",
    )

    values = result[
        "edge_atoms"
    ].tolist()

    assert values[0] == 7
    assert pd.isna(values[1])
    assert pd.isna(values[2])


def test_add_column_safe_empty_list_creates_empty_nullable_column():
    df = pd.DataFrame(
        {
            "base": [1, 2],
        }
    )

    result = add_column_safe(
        df,
        [],
        "edge_atoms",
    )

    assert (
        "edge_atoms"
        in result.columns
    )

    assert result[
        "edge_atoms"
    ].isna().all()


def test_add_column_safe_overwrites_existing_column():
    df = pd.DataFrame(
        {
            "edge_atoms": [100, 200],
        }
    )

    result = add_column_safe(
        df,
        [1, 2],
        "edge_atoms",
    )

    assert (
        result[
            "edge_atoms"
        ].tolist()
        == [1, 2]
    )


def test_add_column_safe_is_position_safe_with_nondefault_dataframe_index():
    """
    Values should be assigned by dataframe row position, not pandas index
    label alignment.
    """
    df = pd.DataFrame(
        {
            "base": ["a", "b"],
        },
        index=[50, 100],
    )

    result = add_column_safe(
        df,
        [3, 4],
        "edge_atoms",
    )

    assert (
        result[
            "edge_atoms"
        ].tolist()
        == [3, 4]
    )


# =============================================================================
# extract_unique_references
# =============================================================================


def test_extract_unique_references_empty_input():
    assert (
        extract_unique_references([])
        == []
    )


def test_extract_unique_references_extracts_string_values():
    reaction = SimpleNamespace(
        references={
            "paper": "doi:one",
            "source": "doi:two",
        }
    )

    result = extract_unique_references(
        [reaction]
    )

    assert result == [
        "doi:one",
        "doi:two",
    ]


def test_extract_unique_references_extracts_list_values():
    reaction = SimpleNamespace(
        references={
            "papers": [
                "doi:one",
                "doi:two",
            ]
        }
    )

    result = extract_unique_references(
        [reaction]
    )

    assert result == [
        "doi:one",
        "doi:two",
    ]


def test_extract_unique_references_extracts_tuple_values():
    reaction = SimpleNamespace(
        references={
            "papers": (
                "doi:one",
                "doi:two",
            )
        }
    )

    result = extract_unique_references(
        [reaction]
    )

    assert result == [
        "doi:one",
        "doi:two",
    ]


def test_extract_unique_references_extracts_nested_dict_string_values():
    reaction = SimpleNamespace(
        references={
            "papers": {
                "primary": "doi:one",
                "secondary": "doi:two",
            }
        }
    )

    result = extract_unique_references(
        [reaction]
    )

    assert result == [
        "doi:one",
        "doi:two",
    ]


def test_extract_unique_references_removes_duplicates_preserving_first_order():
    reaction_1 = SimpleNamespace(
        references={
            "a": "doi:one",
            "b": [
                "doi:two",
                "doi:one",
            ],
        }
    )

    reaction_2 = SimpleNamespace(
        references={
            "c": {
                "x": "doi:two",
                "y": "doi:three",
            }
        }
    )

    result = extract_unique_references(
        [
            reaction_1,
            reaction_2,
        ]
    )

    assert result == [
        "doi:one",
        "doi:two",
        "doi:three",
    ]


def test_extract_unique_references_ignores_non_string_values():
    reaction = SimpleNamespace(
        references={
            "number": 123,
            "none": None,
            "list": [
                1,
                "doi:one",
                None,
            ],
            "dict": {
                "a": 10,
                "b": "doi:two",
            },
        }
    )

    result = extract_unique_references(
        [reaction]
    )

    assert result == [
        "doi:one",
        "doi:two",
    ]


def test_extract_unique_references_uses_legacy_reference_fallback():
    reaction = SimpleNamespace(
        references=None,
        reference={
            "paper": "legacy-doi",
        },
    )

    result = extract_unique_references(
        [reaction]
    )

    assert result == [
        "legacy-doi",
    ]


def test_extract_unique_references_uses_legacy_reference_when_references_attribute_missing():
    reaction = SimpleNamespace(
        reference={
            "paper": "legacy-doi",
        }
    )

    result = extract_unique_references(
        [reaction]
    )

    assert result == [
        "legacy-doi",
    ]


def test_extract_unique_references_empty_references_does_not_use_legacy_fallback():
    """
    An existing references dictionary is the preferred current API.
    Legacy reference is only a fallback when references is unavailable/None.
    """
    reaction = SimpleNamespace(
        references={},
        reference={
            "paper": "legacy-doi",
        },
    )

    result = extract_unique_references(
        [reaction]
    )

    assert result == []


# =============================================================================
# compare_set
# =============================================================================


def make_metadata(
    reactant_smiles,
    product_smiles,
):
    return SimpleNamespace(
        reactant_combined_RDmol=Chem.MolFromSmiles(
            reactant_smiles
        ),
        product_combined_RDmol=Chem.MolFromSmiles(
            product_smiles
        ),
    )


def test_compare_set_returns_true_for_empty_existing_list():
    reactant = Chem.MolFromSmiles(
        "CC"
    )
    product = Chem.MolFromSmiles(
        "CO"
    )

    assert (
        compare_set(
            [],
            reactant,
            product,
        )
        is True
    )


def test_compare_set_returns_false_for_identical_pair():
    existing = [
        make_metadata(
            "CC",
            "CO",
        )
    ]

    result = compare_set(
        existing,
        Chem.MolFromSmiles("CC"),
        Chem.MolFromSmiles("CO"),
    )

    assert result is False


def test_compare_set_uses_canonical_structure_not_smiles_order():
    existing = [
        make_metadata(
            "CCO",
            "CCN",
        )
    ]

    result = compare_set(
        existing,
        Chem.MolFromSmiles("OCC"),
        Chem.MolFromSmiles("NCC"),
    )

    assert result is False


def test_compare_set_returns_true_when_reactant_differs():
    existing = [
        make_metadata(
            "CC",
            "CO",
        )
    ]

    result = compare_set(
        existing,
        Chem.MolFromSmiles("CCC"),
        Chem.MolFromSmiles("CO"),
    )

    assert result is True


def test_compare_set_returns_true_when_product_differs():
    existing = [
        make_metadata(
            "CC",
            "CO",
        )
    ]

    result = compare_set(
        existing,
        Chem.MolFromSmiles("CC"),
        Chem.MolFromSmiles("CN"),
    )

    assert result is True


def test_compare_set_ignores_atom_mapping_numbers():
    existing_reactant = Chem.MolFromSmiles(
        "CC"
    )
    existing_product = Chem.MolFromSmiles(
        "CO"
    )

    existing_reactant.GetAtomWithIdx(
        0
    ).SetAtomMapNum(11)

    existing_product.GetAtomWithIdx(
        0
    ).SetAtomMapNum(22)

    metadata = SimpleNamespace(
        reactant_combined_RDmol=existing_reactant,
        product_combined_RDmol=existing_product,
    )

    new_reactant = Chem.MolFromSmiles(
        "CC"
    )
    new_product = Chem.MolFromSmiles(
        "CO"
    )

    new_reactant.GetAtomWithIdx(
        1
    ).SetAtomMapNum(500)

    new_product.GetAtomWithIdx(
        1
    ).SetAtomMapNum(600)

    assert (
        compare_set(
            [metadata],
            new_reactant,
            new_product,
        )
        is False
    )


def test_compare_set_does_not_mutate_new_molecules_atom_maps():
    reactant = Chem.MolFromSmiles(
        "CC"
    )
    product = Chem.MolFromSmiles(
        "CO"
    )

    reactant.GetAtomWithIdx(
        0
    ).SetAtomMapNum(123)

    product.GetAtomWithIdx(
        0
    ).SetAtomMapNum(456)

    compare_set(
        [],
        reactant,
        product,
    )

    assert (
        reactant.GetAtomWithIdx(
            0
        ).GetAtomMapNum()
        == 123
    )

    assert (
        product.GetAtomWithIdx(
            0
        ).GetAtomMapNum()
        == 456
    )


def test_compare_set_does_not_mutate_existing_metadata_molecules():
    metadata = make_metadata(
        "CC",
        "CO",
    )

    metadata.reactant_combined_RDmol.GetAtomWithIdx(
        0
    ).SetAtomMapNum(101)

    metadata.product_combined_RDmol.GetAtomWithIdx(
        0
    ).SetAtomMapNum(202)

    compare_set(
        [metadata],
        Chem.MolFromSmiles("CC"),
        Chem.MolFromSmiles("CO"),
    )

    assert (
        metadata
        .reactant_combined_RDmol
        .GetAtomWithIdx(0)
        .GetAtomMapNum()
        == 101
    )

    assert (
        metadata
        .product_combined_RDmol
        .GetAtomWithIdx(0)
        .GetAtomMapNum()
        == 202
    )


# =============================================================================
# compare_rdkit_molecules_canonical
# =============================================================================


def test_compare_rdkit_molecules_canonical_finds_identical_structure():
    smiles = [
        "CCO",
    ]

    result_list, found = (
        compare_rdkit_molecules_canonical(
            smiles,
            "OCC",
        )
    )

    assert found is True
    assert result_list is smiles

    assert smiles == [
        "CCO",
    ]


def test_compare_rdkit_molecules_canonical_appends_new_structure():
    smiles = [
        "CCO",
    ]

    result_list, found = (
        compare_rdkit_molecules_canonical(
            smiles,
            "CCN",
        )
    )

    assert found is False
    assert result_list is smiles

    assert smiles == [
        "CCO",
        "CCN",
    ]


def test_compare_rdkit_molecules_canonical_appends_to_empty_list():
    smiles = []

    result_list, found = (
        compare_rdkit_molecules_canonical(
            smiles,
            "CC",
        )
    )

    assert found is False
    assert result_list is smiles

    assert smiles == [
        "CC",
    ]


@pytest.mark.parametrize(
    "candidate",
    [
        None,
        "",
    ],
)
def test_compare_rdkit_molecules_canonical_rejects_empty_candidate_without_appending(
    candidate,
):
    smiles = [
        "CC",
    ]

    result_list, found = (
        compare_rdkit_molecules_canonical(
            smiles,
            candidate,
        )
    )

    assert found is False
    assert result_list is smiles

    assert smiles == [
        "CC",
    ]


def test_compare_rdkit_molecules_canonical_rejects_invalid_candidate_without_appending():
    smiles = [
        "CC",
    ]

    result_list, found = (
        compare_rdkit_molecules_canonical(
            smiles,
            "this-is-not-smiles",
        )
    )

    assert found is False
    assert result_list is smiles

    assert smiles == [
        "CC",
    ]


def test_compare_rdkit_molecules_canonical_invalid_existing_entry_does_not_hide_later_match():
    """
    One malformed cached entry should not prevent the function from checking
    later valid entries for a chemical duplicate.
    """
    smiles = [
        "this-is-not-smiles",
        "CCO",
    ]

    result_list, found = (
        compare_rdkit_molecules_canonical(
            smiles,
            "OCC",
        )
    )

    assert found is True
    assert result_list is smiles

    assert smiles == [
        "this-is-not-smiles",
        "CCO",
    ]


# =============================================================================
# prep_for_3d_molecule_generation
# =============================================================================


def test_prep_for_3d_molecule_generation_builds_numbered_data_molecules():
    result = (
        prep_for_3d_molecule_generation(
            [
                "C",
                "CC",
            ],
            {},
        )
    )

    assert set(result) == {
        "data_1",
        "data_2",
    }

    assert isinstance(
        result["data_1"],
        Chem.Mol,
    )

    assert isinstance(
        result["data_2"],
        Chem.Mol,
    )


def test_prep_for_3d_molecule_generation_adds_explicit_hydrogens():
    result = (
        prep_for_3d_molecule_generation(
            [
                "C",
            ],
            {},
        )
    )

    methane = result[
        "data_1"
    ]

    # CH4 = 1 carbon + 4 explicit hydrogen atoms.
    assert (
        methane.GetNumAtoms()
        == 5
    )

    hydrogen_count = sum(
        atom.GetAtomicNum() == 1
        for atom in methane.GetAtoms()
    )

    assert hydrogen_count == 4


def test_prep_for_3d_molecule_generation_rejects_invalid_smiles():
    with pytest.raises(
        ValueError,
        match="Invalid SMILES",
    ):
        prep_for_3d_molecule_generation(
            [
                "CC",
                "not-a-smiles",
            ],
            {},
        )


def test_prep_for_3d_molecule_generation_error_mentions_invalid_smiles():
    invalid = "definitely-not-smiles"

    with pytest.raises(
        ValueError,
    ) as exc_info:
        prep_for_3d_molecule_generation(
            [
                invalid,
            ],
            {},
        )

    assert (
        invalid
        in str(exc_info.value)
    )


def test_prep_for_3d_molecule_generation_adds_pre_and_post_reaction_molecules():
    reactant = Chem.MolFromSmiles(
        "CC"
    )
    product = Chem.MolFromSmiles(
        "CO"
    )

    reaction_data = {
        "reaction_1": {
            "reactant": reactant,
            "product": product,
        }
    }

    result = (
        prep_for_3d_molecule_generation(
            [],
            reaction_data,
        )
    )

    assert (
        result["pre_1"]
        is reactant
    )

    assert (
        result["post_1"]
        is product
    )


def test_prep_for_3d_molecule_generation_skips_entry_missing_reactant():
    product = Chem.MolFromSmiles(
        "CO"
    )

    result = (
        prep_for_3d_molecule_generation(
            [],
            {
                "reaction": {
                    "product": product,
                }
            },
        )
    )

    assert result == {}


def test_prep_for_3d_molecule_generation_skips_entry_missing_product():
    reactant = Chem.MolFromSmiles(
        "CC"
    )

    result = (
        prep_for_3d_molecule_generation(
            [],
            {
                "reaction": {
                    "reactant": reactant,
                }
            },
        )
    )

    assert result == {}


def test_prep_for_3d_molecule_generation_skips_entry_with_none_values():
    result = (
        prep_for_3d_molecule_generation(
            [],
            {
                "reaction": {
                    "reactant": None,
                    "product": None,
                }
            },
        )
    )

    assert result == {}


def test_prep_for_3d_molecule_generation_combines_input_and_reaction_molecules():
    reactant = Chem.MolFromSmiles(
        "CO"
    )
    product = Chem.MolFromSmiles(
        "CN"
    )

    result = (
        prep_for_3d_molecule_generation(
            [
                "C",
                "CC",
            ],
            {
                "rxn": {
                    "reactant": reactant,
                    "product": product,
                }
            },
        )
    )

    assert set(result) == {
        "data_1",
        "data_2",
        "pre_1",
        "post_1",
    }


def test_prep_for_3d_molecule_generation_preserves_reaction_dictionary_order_numbering():
    reactant_1 = Chem.MolFromSmiles(
        "CC"
    )
    product_1 = Chem.MolFromSmiles(
        "CO"
    )

    reactant_2 = Chem.MolFromSmiles(
        "CCC"
    )
    product_2 = Chem.MolFromSmiles(
        "CCO"
    )

    result = (
        prep_for_3d_molecule_generation(
            [],
            {
                "first": {
                    "reactant": reactant_1,
                    "product": product_1,
                },
                "second": {
                    "reactant": reactant_2,
                    "product": product_2,
                },
            },
        )
    )

    assert result[
        "pre_1"
    ] is reactant_1

    assert result[
        "post_1"
    ] is product_1

    assert result[
        "pre_2"
    ] is reactant_2

    assert result[
        "post_2"
    ] is product_2