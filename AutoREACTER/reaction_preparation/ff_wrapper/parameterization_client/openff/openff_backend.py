"""OpenFF Toolkit + Interchange parameterization backend for AutoREACTER.

Uses optimized in-memory RDKit molecules, applies an OpenFF SMIRNOFF
force field, and exports native LAMMPS data files.

The shared LAMMPS merger must separately support OpenFF's exported
styles and parameter identities before full AutoREACTER integration.
"""

from pathlib import Path

import numpy as np
from rdkit import Chem

from openff.toolkit import ForceField, Molecule, Topology
from openff.interchange import Interchange
from openff.units import unit


class OpenFFBackend:
    name = "openff"

    def __init__(
        self,
        ff_spec: str | Path,
        *,
        include_type_labels: bool = False,
    ) -> None:
        if not ff_spec:
            raise ValueError("An OpenFF .offxml force field must be specified.")

        self._forcefield = ForceField(str(ff_spec))
        self._used_ids: set[str] = set()
        self.include_type_labels = include_type_labels

    def parameterize(
        self,
        mol: Chem.Mol,
        molecule_id: str,
        *,
        output_dir: Path,
    ) -> Path:
        """Parameterize one RDKit molecule and write <id>_typed.data."""

        self._validate_molecule_id(molecule_id)

        if molecule_id in self._used_ids:
            raise ValueError(
                f"Duplicate OpenFF molecule/template ID '{molecule_id}'."
            )

        if not isinstance(mol, Chem.Mol):
            raise TypeError(f"{molecule_id}: expected an RDKit Chem.Mol.")

        if mol.GetNumAtoms() == 0:
            raise ValueError(f"{molecule_id}: RDKit molecule is empty.")

        if mol.GetNumConformers() == 0:
            raise ValueError(
                f"{molecule_id}: no 3D conformer in RDKit molecule."
            )

        # AutoREACTER's 3D preparation supplies explicit hydrogens.
        # Do not allow OpenFF conversion to silently change atom count.
        if any(atom.GetNumImplicitHs() for atom in mol.GetAtoms()):
            raise ValueError(
                f"{molecule_id}: molecule has implicit hydrogens. "
                "Explicit hydrogens are required to preserve atom indices."
            )

        input_xyz = np.asarray(
            mol.GetConformer().GetPositions(),
            dtype=float,
        )

        if (
            input_xyz.shape != (mol.GetNumAtoms(), 3)
            or not np.isfinite(input_xyz).all()
        ):
            raise ValueError(
                f"{molecule_id}: invalid RDKit coordinates."
            )

        # Convert the existing in-memory RDKit molecule.
        # Undefined stereochemistry is rejected rather than guessed.
        off_mol = Molecule.from_rdkit(
            mol,
            allow_undefined_stereo=False,
            hydrogens_are_explicit=True,
        )

        if off_mol.n_atoms != mol.GetNumAtoms():
            raise RuntimeError(
                f"{molecule_id}: atom count changed during OpenFF conversion."
            )

        expected_elements = [
            atom.GetAtomicNum() for atom in mol.GetAtoms()
        ]
        actual_elements = [
            atom.atomic_number for atom in off_mol.atoms
        ]

        if expected_elements != actual_elements:
            raise RuntimeError(
                f"{molecule_id}: atom order/elements changed during conversion."
            )

        if not off_mol.conformers:
            raise RuntimeError(
                f"{molecule_id}: OpenFF conversion lost the 3D conformer."
            )

        off_xyz = np.asarray(
            off_mol.conformers[0].m_as(unit.angstrom),
            dtype=float,
        )

        if not np.allclose(
            off_xyz, input_xyz, rtol=0.0, atol=1e-4
        ):
            raise RuntimeError(
                f"{molecule_id}: coordinates or atom ordering changed."
            )

        topology = Topology.from_molecules([off_mol])

        # Rectangular simulation box with 10 Å padding on each side.
        lengths = np.ptp(input_xyz, axis=0) + 20.0
        box = unit.Quantity(
            np.diag(lengths),
            unit.angstrom,
        )

        # Assign SMIRNOFF parameters and partial charges.
        interchange = Interchange.from_smirnoff(
            force_field=self._forcefield,
            topology=topology,
            box=box,
        )

        if interchange.topology.n_atoms != mol.GetNumAtoms():
            raise RuntimeError(
                f"{molecule_id}: OpenFF parameterization changed atom count."
            )

        if interchange.positions is None:
            raise RuntimeError(
                f"{molecule_id}: Interchange has no coordinates."
            )

        interchange_xyz = np.asarray(
            interchange.positions.m_as(unit.angstrom),
            dtype=float,
        )

        if (
            interchange_xyz.shape != input_xyz.shape
            or not np.allclose(
                interchange_xyz, input_xyz, rtol=0.0, atol=1e-4
            )
        ):
            raise RuntimeError(
                f"{molecule_id}: Interchange changed atom ordering or coordinates."
            )

        # Native OpenFF Interchange LAMMPS output.
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        out_path = output_dir / f"{molecule_id}_typed.data"

        self._used_ids.add(molecule_id)

        try:
            interchange.to_lammps_datafile(
                str(out_path),
                include_type_labels=self.include_type_labels,
            )
        except Exception:
            out_path.unlink(missing_ok=True)
            raise

        if not out_path.is_file() or out_path.stat().st_size == 0:
            raise RuntimeError(
                f"{molecule_id}: OpenFF produced no usable LAMMPS file."
            )

        print(f"[OpenFF] Parameterized {molecule_id} -> {out_path}")

        return out_path

    @staticmethod
    def _validate_molecule_id(molecule_id: str) -> None:
        """Reject unsafe IDs that could corrupt output paths."""
        if (
            not isinstance(molecule_id, str)
            or not molecule_id
            or molecule_id in {".", ".."}
            or "/" in molecule_id
            or "\\" in molecule_id
            or "\x00" in molecule_id
            or Path(molecule_id).name != molecule_id
        ):
            raise ValueError(
                f"Unsafe molecule/template ID: {molecule_id!r}"
            )
