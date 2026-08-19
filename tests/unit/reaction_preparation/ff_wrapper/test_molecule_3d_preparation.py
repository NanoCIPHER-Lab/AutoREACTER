from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from rdkit import Chem

import AutoREACTER.reaction_preparation.ff_wrapper.molecule_3d_preparation as molecule_3d
from AutoREACTER.reaction_preparation.ff_wrapper.molecule_3d_preparation import (
    FragmentSeparationError,
    Molecule3DPreparation,
    Molecule3DPreparationError,
    OptimizationError,
)


# =============================================================================
# Helpers
# =============================================================================


def make_session(
    tmp_path,
    *,
    monomers=None,
    reactions=None,
):
    inputs = SimpleNamespace(
        monomers=list(
            monomers or []
        )
    )

    return SimpleNamespace(
        inputs=inputs,
        reaction_metadata=list(
            reactions or []
        ),
        staging_dir=tmp_path / "staging",
    )


def make_monomer(
    *,
    name="mma",
    status=True,
    rdkit_mol=None,
):
    return SimpleNamespace(
        name=name,
        status=status,
        rdkit_mol=rdkit_mol,
        molecule_3Dmol_path=None,
    )


def make_reaction(
    *,
    reaction_id=1,
    activity_stats=True,
    reactant=None,
    product=None,
):
    return SimpleNamespace(
        reaction_id=reaction_id,
        activity_stats=activity_stats,
        reactant_combined_RDmol=reactant,
        product_combined_RDmol=product,
        reactant_combined_3Dmol_path=None,
        product_combined_3Dmol_path=None,
    )


def make_two_fragment_mol():
    """
    Create two disconnected carbon atoms with a manually assigned conformer.
    No embedding is required for the separation tests.
    """
    mol = Chem.MolFromSmiles(
        "[C].[C]"
    )

    conf = Chem.Conformer(
        mol.GetNumAtoms()
    )

    conf.SetAtomPosition(
        0,
        (
            1.0,
            2.0,
            3.0,
        ),
    )

    conf.SetAtomPosition(
        1,
        (
            10.0,
            20.0,
            30.0,
        ),
    )

    mol.AddConformer(
        conf
    )

    return mol


def make_three_fragment_mol():
    mol = Chem.MolFromSmiles(
        "[C].[C].[C]"
    )

    conf = Chem.Conformer(
        mol.GetNumAtoms()
    )

    for i in range(
        mol.GetNumAtoms()
    ):
        conf.SetAtomPosition(
            i,
            (
                float(i),
                0.0,
                0.0,
            ),
        )

    mol.AddConformer(
        conf
    )

    return mol


def make_explicit_ch3_with_duplicate_h_property():
    """
    Build a carbon with three real H neighbors plus an erroneous
    NumExplicitHs property of 1.

    This reproduces the kind of duplicate hydrogen bookkeeping that
    _repair_reaction_molecule_for_3d() is designed to repair.
    """
    editable = Chem.RWMol()

    carbon = Chem.Atom(6)
    carbon.SetNumExplicitHs(1)

    carbon_idx = editable.AddAtom(
        carbon
    )

    for _ in range(3):
        hydrogen_idx = (
            editable.AddAtom(
                Chem.Atom(1)
            )
        )

        editable.AddBond(
            carbon_idx,
            hydrogen_idx,
            Chem.BondType.SINGLE,
        )

    mol = editable.GetMol()

    mol.UpdatePropertyCache(
        strict=False
    )

    return mol


def make_explicit_ch4_with_radical():
    """
    Build carbon with four explicit H neighbors but incorrectly mark it
    as a radical. Repair should clear that radical because bond valence
    is already four.
    """
    editable = Chem.RWMol()

    carbon = Chem.Atom(6)
    carbon.SetNumRadicalElectrons(
        1
    )

    carbon_idx = editable.AddAtom(
        carbon
    )

    for _ in range(4):
        hydrogen_idx = (
            editable.AddAtom(
                Chem.Atom(1)
            )
        )

        editable.AddBond(
            carbon_idx,
            hydrogen_idx,
            Chem.BondType.SINGLE,
        )

    mol = editable.GetMol()

    mol.UpdatePropertyCache(
        strict=False
    )

    return mol


def patch_basic_optimization_stack(
    preparer,
    monkeypatch,
    *,
    embed_result=0,
    mmff_available=True,
    mmff_result=0,
    uff_available=False,
    uff_result=0,
):
    """
    Patch the expensive/variable RDKit 3D operations while preserving
    _optimization() control flow.
    """
    events = []

    params = SimpleNamespace(
        randomSeed=None,
        useRandomCoords=False,
        ignoreSmoothingFailures=False,
    )

    monkeypatch.setattr(
        preparer,
        "_repair_reaction_molecule_for_3d",
        lambda mol: mol,
    )

    monkeypatch.setattr(
        molecule_3d.AllChem,
        "ETKDGv3",
        lambda: params,
    )

    def fake_embed(
        mol,
        passed_params,
    ):
        events.append(
            (
                "embed",
                mol,
                passed_params,
            )
        )

        return embed_result

    monkeypatch.setattr(
        molecule_3d.AllChem,
        "EmbedMolecule",
        fake_embed,
    )

    monkeypatch.setattr(
        molecule_3d.AllChem,
        "MMFFHasAllMoleculeParams",
        lambda mol:
        mmff_available,
    )

    def fake_mmff(
        mol,
        maxIters,
    ):
        events.append(
            (
                "mmff",
                maxIters,
            )
        )

        return mmff_result

    monkeypatch.setattr(
        molecule_3d.AllChem,
        "MMFFOptimizeMolecule",
        fake_mmff,
    )

    monkeypatch.setattr(
        molecule_3d.AllChem,
        "UFFHasAllMoleculeParams",
        lambda mol:
        uff_available,
    )

    def fake_uff(
        mol,
        maxIters,
    ):
        events.append(
            (
                "uff",
                maxIters,
            )
        )

        return uff_result

    monkeypatch.setattr(
        molecule_3d.AllChem,
        "UFFOptimizeMolecule",
        fake_uff,
    )

    def fake_write(
        mol,
        filename,
        **kwargs,
    ):
        events.append(
            (
                "write",
                Path(filename),
                kwargs,
            )
        )

        Path(
            filename
        ).write_text(
            "fake mol file",
            encoding="utf-8",
        )

    monkeypatch.setattr(
        molecule_3d.Chem,
        "MolToMolFile",
        fake_write,
    )

    return events, params


# =============================================================================
# Exception hierarchy
# =============================================================================


def test_molecule_3d_preparation_error_is_exception():
    assert issubclass(
        Molecule3DPreparationError,
        Exception,
    )


def test_fragment_separation_error_inherits_base_error():
    assert issubclass(
        FragmentSeparationError,
        Molecule3DPreparationError,
    )


def test_optimization_error_inherits_base_error():
    assert issubclass(
        OptimizationError,
        Molecule3DPreparationError,
    )


# =============================================================================
# Constructor / cache
# =============================================================================


def test_constructor_stores_session(
    tmp_path,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    assert (
        preparer.session
        is session
    )

    assert (
        preparer.inputs
        is session.inputs
    )


def test_constructor_creates_cache_directories(
    tmp_path,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    assert (
        preparer.cache_dir
        == tmp_path
        / "staging"
        / "3D_molecules"
    )

    assert (
        preparer.molecule_3d_path
        == preparer.cache_dir
        / "molecules_3Dmol"
    )

    assert (
        preparer.full_templates_path
        == preparer.cache_dir
        / "full_templates_3Dmol"
    )

    assert (
        preparer.cache_dir.is_dir()
    )

    assert (
        preparer
        .molecule_3d_path
        .is_dir()
    )

    assert (
        preparer
        .full_templates_path
        .is_dir()
    )


def test_cache_property_returns_main_cache(
    tmp_path,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    assert (
        preparer.cache
        == preparer.cache_dir
    )


# =============================================================================
# prepare_molecule_3d_geometry - monomers
# =============================================================================


def test_prepare_active_monomer_calls_add_hs_and_optimization(
    tmp_path,
    monkeypatch,
):
    original_mol = (
        Chem.MolFromSmiles(
            "CC"
        )
    )

    monomer = make_monomer(
        name="ethane",
        status=True,
        rdkit_mol=original_mol,
    )

    session = make_session(
        tmp_path,
        monomers=[
            monomer
        ],
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    added_hs = object()

    add_h_calls = []

    monkeypatch.setattr(
        molecule_3d.Chem,
        "AddHs",
        lambda mol:
        (
            add_h_calls.append(
                mol
            )
            or added_hs
        ),
    )

    optimization_calls = []

    output_path = (
        preparer.molecule_3d_path
        / "ethane.mol"
    )

    def fake_optimization(
        *,
        molecule_name,
        mol,
        cache_dir,
        separate_fragments=False,
    ):
        optimization_calls.append(
            (
                molecule_name,
                mol,
                cache_dir,
                separate_fragments,
            )
        )

        return output_path

    monkeypatch.setattr(
        preparer,
        "_optimization",
        fake_optimization,
    )

    result = (
        preparer
        .prepare_molecule_3d_geometry(
            session
        )
    )

    assert result is None

    assert add_h_calls == [
        original_mol
    ]

    assert optimization_calls == [
        (
            "ethane",
            added_hs,
            preparer.molecule_3d_path,
            False,
        )
    ]

    assert (
        monomer.molecule_3Dmol_path
        == output_path
    )


def test_prepare_skips_inactive_monomer(
    tmp_path,
    monkeypatch,
):
    monomer = make_monomer(
        name="inactive",
        status=False,
        rdkit_mol=None,
    )

    session = make_session(
        tmp_path,
        monomers=[
            monomer
        ],
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    monkeypatch.setattr(
        preparer,
        "_optimization",
        lambda **kwargs:
        pytest.fail(
            "inactive monomer must not be optimized"
        ),
    )

    preparer.prepare_molecule_3d_geometry(
        session
    )

    assert (
        monomer.molecule_3Dmol_path
        is None
    )


def test_prepare_active_monomer_without_rdkit_mol_raises(
    tmp_path,
):
    monomer = make_monomer(
        name="mma",
        status=True,
        rdkit_mol=None,
    )

    session = make_session(
        tmp_path,
        monomers=[
            monomer
        ],
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    with pytest.raises(
        Molecule3DPreparationError,
        match="RDKit Mol object is missing for molecule mma",
    ):
        preparer.prepare_molecule_3d_geometry(
            session
        )


def test_prepare_monomer_optimization_error_is_wrapped(
    tmp_path,
    monkeypatch,
):
    monomer = make_monomer(
        name="styrene",
        status=True,
        rdkit_mol=(
            Chem.MolFromSmiles(
                "C=C"
            )
        ),
    )

    session = make_session(
        tmp_path,
        monomers=[
            monomer
        ],
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    monkeypatch.setattr(
        preparer,
        "_optimization",
        lambda **kwargs:
        (_ for _ in ()).throw(
            RuntimeError(
                "boom"
            )
        ),
    )

    with pytest.raises(
        OptimizationError,
        match="Error optimizing molecule styrene: boom",
    ) as exc_info:
        preparer.prepare_molecule_3d_geometry(
            session
        )

    assert isinstance(
        exc_info.value.__cause__,
        RuntimeError,
    )


# =============================================================================
# prepare_molecule_3d_geometry - reactions
# =============================================================================


def test_prepare_skips_inactive_reaction(
    tmp_path,
    monkeypatch,
):
    reaction = make_reaction(
        reaction_id=5,
        activity_stats=False,
        reactant=object(),
        product=object(),
    )

    session = make_session(
        tmp_path,
        reactions=[
            reaction
        ],
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    monkeypatch.setattr(
        preparer,
        "_optimization",
        lambda **kwargs:
        pytest.fail(
            "inactive reaction must not be optimized"
        ),
    )

    preparer.prepare_molecule_3d_geometry(
        session
    )

    assert (
        reaction
        .reactant_combined_3Dmol_path
        is None
    )

    assert (
        reaction
        .product_combined_3Dmol_path
        is None
    )


def test_prepare_reaction_pre_and_post(
    tmp_path,
    monkeypatch,
):
    pre_mol = object()
    post_mol = object()

    reaction = make_reaction(
        reaction_id=12,
        activity_stats=True,
        reactant=pre_mol,
        product=post_mol,
    )

    session = make_session(
        tmp_path,
        reactions=[
            reaction
        ],
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    calls = []

    def fake_optimization(
        *,
        molecule_name,
        mol,
        cache_dir,
        separate_fragments=False,
    ):
        calls.append(
            (
                molecule_name,
                mol,
                cache_dir,
                separate_fragments,
            )
        )

        return (
            Path(cache_dir)
            / f"{molecule_name}.mol"
        )

    monkeypatch.setattr(
        preparer,
        "_optimization",
        fake_optimization,
    )

    preparer.prepare_molecule_3d_geometry(
        session
    )

    assert calls == [
        (
            "pre12",
            pre_mol,
            preparer.full_templates_path,
            True,
        ),
        (
            "post12",
            post_mol,
            preparer.full_templates_path,
            True,
        ),
    ]

    assert (
        reaction
        .reactant_combined_3Dmol_path
        == preparer.full_templates_path
        / "pre12.mol"
    )

    assert (
        reaction
        .product_combined_3Dmol_path
        == preparer.full_templates_path
        / "post12.mol"
    )


def test_prepare_reaction_missing_pre_molecule_is_skipped(
    tmp_path,
    monkeypatch,
):
    product = object()

    reaction = make_reaction(
        reaction_id=2,
        reactant=None,
        product=product,
    )

    session = make_session(
        tmp_path,
        reactions=[
            reaction
        ],
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    calls = []

    monkeypatch.setattr(
        preparer,
        "_optimization",
        lambda **kwargs:
        (
            calls.append(
                kwargs["molecule_name"]
            )
            or Path(
                f"/tmp/{kwargs['molecule_name']}.mol"
            )
        ),
    )

    preparer.prepare_molecule_3d_geometry(
        session
    )

    assert calls == [
        "post2"
    ]


def test_prepare_reaction_missing_post_molecule_is_skipped(
    tmp_path,
    monkeypatch,
):
    reactant = object()

    reaction = make_reaction(
        reaction_id=2,
        reactant=reactant,
        product=None,
    )

    session = make_session(
        tmp_path,
        reactions=[
            reaction
        ],
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    calls = []

    monkeypatch.setattr(
        preparer,
        "_optimization",
        lambda **kwargs:
        (
            calls.append(
                kwargs["molecule_name"]
            )
            or Path(
                f"/tmp/{kwargs['molecule_name']}.mol"
            )
        ),
    )

    preparer.prepare_molecule_3d_geometry(
        session
    )

    assert calls == [
        "pre2"
    ]


def test_prepare_reactant_optimization_error_is_wrapped(
    tmp_path,
    monkeypatch,
):
    reaction = make_reaction(
        reaction_id=8,
        reactant=object(),
        product=None,
    )

    session = make_session(
        tmp_path,
        reactions=[
            reaction
        ],
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    monkeypatch.setattr(
        preparer,
        "_optimization",
        lambda **kwargs:
        (_ for _ in ()).throw(
            RuntimeError(
                "pre failure"
            )
        ),
    )

    with pytest.raises(
        OptimizationError,
        match=(
            "Error optimizing reactant complex "
            "for reaction 8: pre failure"
        ),
    ):
        preparer.prepare_molecule_3d_geometry(
            session
        )


def test_prepare_product_optimization_error_is_wrapped(
    tmp_path,
    monkeypatch,
):
    reaction = make_reaction(
        reaction_id=9,
        reactant=None,
        product=object(),
    )

    session = make_session(
        tmp_path,
        reactions=[
            reaction
        ],
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    monkeypatch.setattr(
        preparer,
        "_optimization",
        lambda **kwargs:
        (_ for _ in ()).throw(
            RuntimeError(
                "post failure"
            )
        ),
    )

    with pytest.raises(
        OptimizationError,
        match=(
            "Error optimizing product complex "
            "for reaction 9: post failure"
        ),
    ):
        preparer.prepare_molecule_3d_geometry(
            session
        )


# =============================================================================
# _separate_fragments_3d
# =============================================================================


def test_separate_fragments_single_fragment_returns_same_object(
    tmp_path,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    mol = Chem.MolFromSmiles(
        "CC"
    )

    result = (
        preparer
        ._separate_fragments_3d(
            mol
        )
    )

    assert result is mol


def test_separate_fragments_more_than_two_raises(
    tmp_path,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    mol = (
        make_three_fragment_mol()
    )

    with pytest.raises(
        FragmentSeparationError,
        match="Expected 2 fragments",
    ):
        preparer._separate_fragments_3d(
            mol
        )


def test_separate_fragments_moves_only_second_fragment(
    tmp_path,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    mol = make_two_fragment_mol()

    conf = mol.GetConformer()

    first_before = np.array(
        conf.GetAtomPosition(0)
    )

    second_before = np.array(
        conf.GetAtomPosition(1)
    )

    result = (
        preparer
        ._separate_fragments_3d(
            mol
        )
    )

    conf_after = (
        result.GetConformer()
    )

    first_after = np.array(
        conf_after.GetAtomPosition(0)
    )

    second_after = np.array(
        conf_after.GetAtomPosition(1)
    )

    np.testing.assert_allclose(
        first_after,
        first_before,
    )

    # [C].[C] has molecular weight ~24,
    # round(24/100) == 0, therefore shift = 4 Å.
    np.testing.assert_allclose(
        second_after,
        second_before
        + np.array(
            [
                4.0,
                0.0,
                0.0,
            ]
        ),
    )


def test_separate_fragments_returns_same_molecule_instance(
    tmp_path,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    mol = make_two_fragment_mol()

    result = (
        preparer
        ._separate_fragments_3d(
            mol
        )
    )

    assert result is mol


# =============================================================================
# _repair_reaction_molecule_for_3d
# =============================================================================


def test_repair_does_not_change_atom_count(
    tmp_path,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    mol = (
        make_explicit_ch3_with_duplicate_h_property()
    )

    count_before = (
        mol.GetNumAtoms()
    )

    repaired = (
        preparer
        ._repair_reaction_molecule_for_3d(
            mol
        )
    )

    assert (
        repaired.GetNumAtoms()
        == count_before
    )


def test_repair_clears_duplicate_explicit_h_property(
    tmp_path,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    mol = (
        make_explicit_ch3_with_duplicate_h_property()
    )

    assert (
        mol.GetAtomWithIdx(0)
        .GetNumExplicitHs()
        == 1
    )

    repaired = (
        preparer
        ._repair_reaction_molecule_for_3d(
            mol
        )
    )

    carbon = (
        repaired.GetAtomWithIdx(
            0
        )
    )

    assert (
        carbon.GetNumExplicitHs()
        == 0
    )


def test_repair_sets_no_implicit_when_real_h_neighbors_exist(
    tmp_path,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    repaired = (
        preparer
        ._repair_reaction_molecule_for_3d(
            make_explicit_ch3_with_duplicate_h_property()
        )
    )

    carbon = (
        repaired.GetAtomWithIdx(
            0
        )
    )

    assert (
        carbon.GetNoImplicit()
        is True
    )


def test_repair_valence_three_neutral_carbon_becomes_radical(
    tmp_path,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    repaired = (
        preparer
        ._repair_reaction_molecule_for_3d(
            make_explicit_ch3_with_duplicate_h_property()
        )
    )

    carbon = (
        repaired.GetAtomWithIdx(
            0
        )
    )

    assert (
        carbon.GetNumRadicalElectrons()
        == 1
    )


def test_repair_valence_four_carbon_clears_radical(
    tmp_path,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    mol = (
        make_explicit_ch4_with_radical()
    )

    assert (
        mol.GetAtomWithIdx(0)
        .GetNumRadicalElectrons()
        == 1
    )

    repaired = (
        preparer
        ._repair_reaction_molecule_for_3d(
            mol
        )
    )

    carbon = (
        repaired.GetAtomWithIdx(
            0
        )
    )

    assert (
        carbon.GetNumRadicalElectrons()
        == 0
    )


def test_repair_does_not_mutate_original_molecule(
    tmp_path,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    original = (
        make_explicit_ch3_with_duplicate_h_property()
    )

    original_explicit_h = (
        original.GetAtomWithIdx(0)
        .GetNumExplicitHs()
    )

    repaired = (
        preparer
        ._repair_reaction_molecule_for_3d(
            original
        )
    )

    assert repaired is not original

    assert (
        original.GetAtomWithIdx(0)
        .GetNumExplicitHs()
        == original_explicit_h
    )

    assert (
        repaired.GetAtomWithIdx(0)
        .GetNumExplicitHs()
        == 0
    )


# =============================================================================
# _optimization - setup / embedding
# =============================================================================


def test_optimization_returns_expected_output_path(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    events, _ = (
        patch_basic_optimization_stack(
            preparer,
            monkeypatch,
        )
    )

    output_dir = (
        tmp_path / "output"
    )

    result = preparer._optimization(
        molecule_name="styrene",
        mol=Chem.MolFromSmiles(
            "CC"
        ),
        cache_dir=output_dir,
    )

    assert result == (
        output_dir
        / "styrene.mol"
    )

    assert result.is_file()


def test_optimization_creates_output_directory(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    patch_basic_optimization_stack(
        preparer,
        monkeypatch,
    )

    output_dir = (
        tmp_path
        / "deep"
        / "nested"
        / "output"
    )

    assert not output_dir.exists()

    preparer._optimization(
        "test",
        Chem.MolFromSmiles(
            "CC"
        ),
        output_dir,
    )

    assert output_dir.is_dir()


def test_optimization_sets_deterministic_etkdg_parameters(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    _, params = (
        patch_basic_optimization_stack(
            preparer,
            monkeypatch,
        )
    )

    preparer._optimization(
        "test",
        Chem.MolFromSmiles(
            "CC"
        ),
        tmp_path,
    )

    assert (
        params.randomSeed
        == 0xF00D
    )

    assert (
        params.useRandomCoords
        is True
    )

    assert (
        params.ignoreSmoothingFailures
        is True
    )


def test_optimization_embedding_failure_raises(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    patch_basic_optimization_stack(
        preparer,
        monkeypatch,
        embed_result=-1,
    )

    with pytest.raises(
        OptimizationError,
        match=(
            "Failed to embed molecule test "
            "in 3D"
        ),
    ):
        preparer._optimization(
            "test",
            Chem.MolFromSmiles(
                "CC"
            ),
            tmp_path,
        )


def test_optimization_repair_failure_is_wrapped(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    def fail_repair(
        mol,
    ):
        raise ValueError(
            "bad valence"
        )

    monkeypatch.setattr(
        preparer,
        "_repair_reaction_molecule_for_3d",
        fail_repair,
    )

    with pytest.raises(
        OptimizationError,
        match=(
            "Failed to repair molecule test "
            "before 3D embedding"
        ),
    ) as exc_info:
        preparer._optimization(
            "test",
            Chem.MolFromSmiles(
                "CC"
            ),
            tmp_path,
        )

    assert isinstance(
        exc_info.value.__cause__,
        ValueError,
    )


# =============================================================================
# _optimization - fragment separation
# =============================================================================


def test_optimization_does_not_separate_by_default(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    patch_basic_optimization_stack(
        preparer,
        monkeypatch,
    )

    monkeypatch.setattr(
        preparer,
        "_separate_fragments_3d",
        lambda mol:
        pytest.fail(
            "fragment separation should not run"
        ),
    )

    preparer._optimization(
        "test",
        Chem.MolFromSmiles(
            "CC"
        ),
        tmp_path,
        separate_fragments=False,
    )


def test_optimization_separates_when_requested(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    patch_basic_optimization_stack(
        preparer,
        monkeypatch,
    )

    calls = []

    monkeypatch.setattr(
        preparer,
        "_separate_fragments_3d",
        lambda mol:
        (
            calls.append(mol)
            or mol
        ),
    )

    preparer._optimization(
        "test",
        Chem.MolFromSmiles(
            "C.C"
        ),
        tmp_path,
        separate_fragments=True,
    )

    assert len(calls) == 1


# =============================================================================
# _optimization - force-field selection
# =============================================================================


def test_optimization_prefers_mmff_when_available(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    events, _ = (
        patch_basic_optimization_stack(
            preparer,
            monkeypatch,
            mmff_available=True,
            uff_available=True,
        )
    )

    preparer._optimization(
        "test",
        Chem.MolFromSmiles(
            "CC"
        ),
        tmp_path,
    )

    assert (
        "mmff",
        1000,
    ) in events

    assert not any(
        event[0] == "uff"
        for event in events
    )


def test_optimization_uses_uff_when_mmff_unavailable(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    events, _ = (
        patch_basic_optimization_stack(
            preparer,
            monkeypatch,
            mmff_available=False,
            uff_available=True,
        )
    )

    preparer._optimization(
        "test",
        Chem.MolFromSmiles(
            "CC"
        ),
        tmp_path,
    )

    assert (
        "uff",
        1000,
    ) in events

    assert not any(
        event[0] == "mmff"
        for event in events
    )


def test_optimization_without_mmff_or_uff_still_saves(
    tmp_path,
    monkeypatch,
    capsys,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    events, _ = (
        patch_basic_optimization_stack(
            preparer,
            monkeypatch,
            mmff_available=False,
            uff_available=False,
        )
    )

    result = preparer._optimization(
        "unsupported",
        Chem.MolFromSmiles(
            "CC"
        ),
        tmp_path,
    )

    assert result.is_file()

    assert (
        "no MMFF or UFF parameters are available"
        in capsys.readouterr().out
    )

    assert not any(
        event[0] in {
            "mmff",
            "uff",
        }
        for event in events
    )


def test_mmff_failure_minus_one_raises(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    patch_basic_optimization_stack(
        preparer,
        monkeypatch,
        mmff_available=True,
        mmff_result=-1,
    )

    with pytest.raises(
        OptimizationError,
        match=(
            "MMFF optimization failed "
            "for test"
        ),
    ):
        preparer._optimization(
            "test",
            Chem.MolFromSmiles(
                "CC"
            ),
            tmp_path,
        )


def test_uff_failure_minus_one_raises(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    patch_basic_optimization_stack(
        preparer,
        monkeypatch,
        mmff_available=False,
        uff_available=True,
        uff_result=-1,
    )

    with pytest.raises(
        OptimizationError,
        match=(
            "UFF optimization failed "
            "for test"
        ),
    ):
        preparer._optimization(
            "test",
            Chem.MolFromSmiles(
                "CC"
            ),
            tmp_path,
        )


@pytest.mark.parametrize(
    "use_mmff",
    [
        True,
        False,
    ],
)
def test_nonconverged_optimization_warns_but_saves(
    tmp_path,
    monkeypatch,
    capsys,
    use_mmff,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    patch_basic_optimization_stack(
        preparer,
        monkeypatch,
        mmff_available=use_mmff,
        mmff_result=1,
        uff_available=not use_mmff,
        uff_result=1,
    )

    result = preparer._optimization(
        "test",
        Chem.MolFromSmiles(
            "CC"
        ),
        tmp_path,
    )

    assert result.is_file()

    output = (
        capsys.readouterr().out
    )

    if use_mmff:
        assert (
            "MMFF optimization did not converge"
            in output
        )
    else:
        assert (
            "UFF optimization did not converge"
            in output
        )


# =============================================================================
# _optimization - atom-count invariant
# =============================================================================


def test_optimization_rejects_atom_count_change_during_repair(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    original = (
        Chem.MolFromSmiles(
            "C"
        )
    )

    monkeypatch.setattr(
        preparer,
        "_repair_reaction_molecule_for_3d",
        lambda mol:
        Chem.AddHs(
            mol
        ),
    )

    params = SimpleNamespace(
        randomSeed=None,
        useRandomCoords=False,
        ignoreSmoothingFailures=False,
    )

    monkeypatch.setattr(
        molecule_3d.AllChem,
        "ETKDGv3",
        lambda: params,
    )

    monkeypatch.setattr(
        molecule_3d.AllChem,
        "EmbedMolecule",
        lambda mol, params: 0,
    )

    monkeypatch.setattr(
        molecule_3d.AllChem,
        "MMFFHasAllMoleculeParams",
        lambda mol: False,
    )

    monkeypatch.setattr(
        molecule_3d.AllChem,
        "UFFHasAllMoleculeParams",
        lambda mol: False,
    )

    monkeypatch.setattr(
        molecule_3d.Chem,
        "MolToMolFile",
        lambda *args, **kwargs:
        pytest.fail(
            "file must not be written "
            "after atom-count mismatch"
        ),
    )

    with pytest.raises(
        OptimizationError,
        match="Atom count mismatch for test",
    ):
        preparer._optimization(
            "test",
            original,
            tmp_path,
        )


# =============================================================================
# _optimization - input molecule isolation
# =============================================================================


def test_optimization_works_on_copy_of_input_molecule(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    original = (
        Chem.MolFromSmiles(
            "CC"
        )
    )

    original.SetProp(
        "ORIGINAL_MARKER",
        "yes",
    )

    seen = []

    def fake_repair(
        mol,
    ):
        seen.append(
            mol
        )

        mol.SetProp(
            "MODIFIED_COPY",
            "yes",
        )

        return mol

    monkeypatch.setattr(
        preparer,
        "_repair_reaction_molecule_for_3d",
        fake_repair,
    )

    params = SimpleNamespace(
        randomSeed=None,
        useRandomCoords=False,
        ignoreSmoothingFailures=False,
    )

    monkeypatch.setattr(
        molecule_3d.AllChem,
        "ETKDGv3",
        lambda: params,
    )

    monkeypatch.setattr(
        molecule_3d.AllChem,
        "EmbedMolecule",
        lambda mol, params: 0,
    )

    monkeypatch.setattr(
        molecule_3d.AllChem,
        "MMFFHasAllMoleculeParams",
        lambda mol: False,
    )

    monkeypatch.setattr(
        molecule_3d.AllChem,
        "UFFHasAllMoleculeParams",
        lambda mol: False,
    )

    monkeypatch.setattr(
        molecule_3d.Chem,
        "MolToMolFile",
        lambda mol, filename, **kwargs:
        Path(filename).write_text(
            "fake",
            encoding="utf-8",
        ),
    )

    preparer._optimization(
        "test",
        original,
        tmp_path,
    )

    assert len(seen) == 1

    assert seen[0] is not original

    assert not original.HasProp(
        "MODIFIED_COPY"
    )

    assert (
        original.GetProp(
            "ORIGINAL_MARKER"
        )
        == "yes"
    )


# =============================================================================
# _optimization - output writer contract
# =============================================================================


def test_optimization_writes_with_stereo_and_without_kekulization(
    tmp_path,
    monkeypatch,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    events, _ = (
        patch_basic_optimization_stack(
            preparer,
            monkeypatch,
        )
    )

    preparer._optimization(
        "test",
        Chem.MolFromSmiles(
            "CC"
        ),
        tmp_path,
    )

    write_event = next(
        event
        for event in events
        if event[0] == "write"
    )

    _, path, kwargs = (
        write_event
    )

    assert path == (
        tmp_path / "test.mol"
    )

    assert (
        kwargs["includeStereo"]
        is True
    )

    assert (
        kwargs["kekulize"]
        is False
    )


def test_optimization_prints_save_location(
    tmp_path,
    monkeypatch,
    capsys,
):
    session = make_session(
        tmp_path
    )

    preparer = (
        Molecule3DPreparation(
            session
        )
    )

    patch_basic_optimization_stack(
        preparer,
        monkeypatch,
    )

    expected = (
        tmp_path / "styrene.mol"
    )

    preparer._optimization(
        "styrene",
        Chem.MolFromSmiles(
            "CC"
        ),
        tmp_path,
    )

    assert (
        f"Saving optimized styrene to {expected}"
        in capsys.readouterr().out
    )
