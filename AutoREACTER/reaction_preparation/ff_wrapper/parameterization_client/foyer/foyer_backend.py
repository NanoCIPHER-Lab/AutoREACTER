"""Foyer parameterization backend for AutoREACTER."""

from pathlib import Path
from typing import Sequence
import re

import numpy as np
import parmed as pmd
from foyer import Forcefield
from mbuild.formats.lammpsdata import write_lammpsdata
from rdkit import Chem


FF_FILES_DIR = Path(__file__).resolve().parents[2] / "FF_files"
SUPPLEMENT_XML = FF_FILES_DIR / "foyer_autoreacter_supplement.xml"

CROSS_FILE_OVERRIDES = {"arx_HCl_Cl": {"opls_401"}}


class FoyerBackend:
    name = "foyer"

    def __init__(self, ff_spec):
        self._forcefield = self._load_forcefield(ff_spec)
        self._used_ids: set[str] = set()

    def parameterize(
        self,
        mol: Chem.Mol,
        molecule_id: str,
        *,
        output_dir: Path,
    ) -> Path:
        return self._run_foyer_atom_typing(
            mol,
            molecule_id=molecule_id,
            output_dir=output_dir,
        )

    @staticmethod
    def _load_forcefield(
        ff_spec: str | Path | Sequence[str | Path],
        include_supplement: bool = True,
    ) -> Forcefield:
        """Load a built-in Foyer FF or XML file(s), plus optional supplement."""
        if ff_spec is None:
            raise ValueError("A Foyer force field must be specified.")

        specs = (
            [ff_spec] if isinstance(ff_spec, (str, Path)) else list(ff_spec)
        )
        if not specs:
            raise ValueError("The Foyer force-field specification is empty.")

        files: list[str] = []
        built_in_names: list[str] = []
        for spec in specs:
            path = Path(spec).expanduser()
            if path.is_file():
                files.append(str(path.resolve()))
            elif path.suffix.lower() == ".xml" or path.is_absolute():
                raise FileNotFoundError(f"Force-field XML file not found: {path}")
            else:
                name = str(spec).lower()
                if name in {"gaff", "gaff2"}:
                    raise ValueError(
                        f"'{name}' is not a Foyer built-in; provide a GAFF XML path."
                    )
                built_in_names.append(name)

        if len(built_in_names) > 1:
            raise ValueError(
                "Use at most one built-in Foyer force field; "
                f"got {built_in_names!r}."
            )

        if include_supplement:
            if not SUPPLEMENT_XML.is_file():
                raise FileNotFoundError(
                    f"AutoREACTER supplemental force field not found: {SUPPLEMENT_XML}"
                )
            files.append(str(SUPPLEMENT_XML))

        ff = Forcefield(
            forcefield_files=files or None,
            name=built_in_names[0] if built_in_names else None,
        )

        for type_name, overridden_names in CROSS_FILE_OVERRIDES.items():
            if type_name not in ff.atomTypeDefinitions:
                continue
            present = {
                name for name in overridden_names if name in ff.atomTypeDefinitions
            }
            ff.atomTypeOverrides[type_name] = set(
                ff.atomTypeOverrides.get(type_name, set())
            ) | present
        return ff

    def _run_foyer_atom_typing(
        self, mol: Chem.Mol, molecule_id: str, *, output_dir: Path
    ) -> Path:
        """Apply Foyer to one RDKit Mol and write ``<id>_typed.data``.

        Only RDKit objects are used as molecular inputs; .mol files are not
        reloaded. Explicit hydrogens/3D conformers must already be present
        where needed, and the supplement must actually provide parameters.
        """
        self._validate_molecule_id(molecule_id)
        if molecule_id in self._used_ids:
            raise ValueError(
                f"Duplicate Foyer molecule/template ID '{molecule_id}': "
                "would overwrite another output."
            )
        self._used_ids.add(molecule_id)

        if not isinstance(mol, Chem.Mol):
            raise TypeError(f"{molecule_id}: expected an RDKit Chem.Mol.")
        if mol.GetNumAtoms() == 0:
            raise ValueError(f"{molecule_id}: RDKit molecule is empty.")
        if mol.GetNumConformers() == 0:
            raise ValueError(
                f"{molecule_id}: no 3D conformer in the RDKit molecule "
                f"(SMILES: {Chem.MolToSmiles(mol)}). Ensure the 3D preparation "
                "stage saves optimized molecules back into Session."
            )

        structure = pmd.load_rdkit(mol)
        self._restore_rdkit_connectivity(mol, structure, molecule_id)

        try:
            typed_structure = self._forcefield.apply(
                structure,
                assert_dihedral_params=False,
            )
        except Exception as exc:
            exc.add_note(self._describe_typing_failure(mol, molecule_id, exc))
            raise

        if len(typed_structure.atoms) != mol.GetNumAtoms():
            raise RuntimeError(
                f"{molecule_id}: Foyer changed atom count "
                f"({mol.GetNumAtoms()} -> {len(typed_structure.atoms)})."
            )
        expected_elements = [a.GetAtomicNum() for a in mol.GetAtoms()]
        actual_elements = [a.atomic_number for a in typed_structure.atoms]
        if expected_elements != actual_elements:
            raise RuntimeError(
                f"{molecule_id}: atom ordering/elements changed during Foyer typing."
            )

        xyz = np.asarray(typed_structure.coordinates, dtype=float)
        if xyz.shape != (mol.GetNumAtoms(), 3) or not np.isfinite(xyz).all():
            raise RuntimeError(f"{molecule_id}: invalid 3D coordinates after typing.")
        extent = np.ptp(xyz, axis=0)
        padded = extent + 20.0  # Å: 10 Å padding on each side.
        typed_structure.box = np.array([*padded, 90.0, 90.0, 90.0])

        out_path = output_dir / f"{molecule_id}_typed.data"
        write_lammpsdata(
            typed_structure,
            str(out_path),
            atom_style="full",
            unit_style="real",
        )
        if not out_path.is_file() or out_path.stat().st_size == 0:
            raise RuntimeError(f"{molecule_id}: LAMMPS writer made no usable file.")
        print(f"[Foyer] Typed {molecule_id} -> {out_path}")
        return out_path

    @staticmethod
    def _validate_molecule_id(molecule_id: str) -> None:
        """Preserve original names while preventing unsafe/ambiguous paths."""
        if (
            not isinstance(molecule_id, str)
            or not molecule_id
            or molecule_id in {".", ".."}
            or "/" in molecule_id
            or "\\" in molecule_id
            or "\x00" in molecule_id
            or Path(molecule_id).name != molecule_id
        ):
            raise ValueError(f"Unsafe molecule/template ID: {molecule_id!r}")

    @staticmethod
    def _restore_rdkit_connectivity(
        mol: Chem.Mol,
        structure: pmd.Structure,
        molecule_id: str,
    ) -> None:
        """Remove spurious PDB-inferred bonds and restore RDKit bond orders."""
        if len(structure.atoms) != mol.GetNumAtoms():
            raise RuntimeError(
                f"{molecule_id}: ParmEd/RDKit atom count mismatch: "
                f"{len(structure.atoms)} vs {mol.GetNumAtoms()}."
            )
        expected_elements = [a.GetAtomicNum() for a in mol.GetAtoms()]
        actual_elements = [a.atomic_number for a in structure.atoms]
        if expected_elements != actual_elements:
            raise RuntimeError(
                f"{molecule_id}: ParmEd changed atom ordering/elements."
            )

        original_bonds = {
            frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())):
                bond.GetBondTypeAsDouble()
            for bond in mol.GetBonds()
        }
        for bond in list(structure.bonds):
            pair = frozenset((bond.atom1.idx, bond.atom2.idx))
            if pair not in original_bonds:
                bond.delete()
                structure.bonds.remove(bond)
            else:
                bond.order = original_bonds[pair]

        converted_pairs = [
            frozenset((bond.atom1.idx, bond.atom2.idx))
            for bond in structure.bonds
        ]
        if (
            len(converted_pairs) != len(original_bonds)
            or set(converted_pairs) != set(original_bonds)
        ):
            missing = sorted(
                tuple(sorted(pair))
                for pair in set(original_bonds) - set(converted_pairs)
            )
            raise RuntimeError(
                f"{molecule_id}: RDKit/ParmEd bond connectivity mismatch. "
                f"Missing bonds: {missing}."
            )

    @staticmethod
    def _describe_typing_failure(
        mol: Chem.Mol, molecule_id: str, err: Exception
    ) -> str:
        """Add molecule/atom context without hiding the original Foyer error."""
        lines = [
            f"While typing '{molecule_id}': {Chem.MolToSmiles(mol)} "
            f"({mol.GetNumAtoms()} RDKit atoms)."
        ]
        match = re.search(r"atom numbered (\d+)", str(err), re.IGNORECASE)
        if match:
            idx = int(match.group(1))
            if 0 <= idx < mol.GetNumAtoms():
                atom = mol.GetAtomWithIdx(idx)
                neighbors = ", ".join(
                    f"{n.GetSymbol()}{n.GetIdx()}"
                    for n in atom.GetNeighbors()
                )
                lines.append(
                    f"Reported atom {idx}: {atom.GetSymbol()}, neighbors: {neighbors}."
                )
                acyl_halide = Chem.MolFromSmarts("[CX3](=O)[F,Cl,Br,I]")
                if any(idx in atoms for atoms in mol.GetSubstructMatches(acyl_halide)):
                    lines.append(
                        "This atom is part of an acyl halide; verify that the "
                        "selected XML files include chemically appropriate "
                        "atom types AND bonded/nonbonded parameters."
                    )
        return "\n".join(lines)
