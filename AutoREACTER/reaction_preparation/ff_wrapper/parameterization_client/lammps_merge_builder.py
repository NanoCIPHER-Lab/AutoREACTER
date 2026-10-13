"""Merge independently parameterized LAMMPS data files for AutoREACTER.

Works on existing .data files, WITHOUT re-atom-typing or reading .mol files.
Maintains global atom/bond/angle/dihedral/improper IDs, writes per-structure
.data/.molecule files and a coefficient-only force_field.data with official
LAMMPS Type Labels sections. Does NOT generate an input script.
Currently supports mBuild-compatible full atom style, units real.

Unsupported/missing or conflicting parameters cause an explicit error rather
than silently creating an unparameterized reactive force field.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import math
import re

from AutoREACTER.reaction_preparation.ff_wrapper.ff_wrapper import (
    DataFiles, FFFiles, MoleculeFile, TemplateFile,
)


CATEGORIES = ("atom", "bond", "angle", "dihedral", "improper")
BONDED = ("bond", "angle", "dihedral", "improper")
COEFF_SECTIONS = {
    "atom": "Pair Coeffs", "bond": "Bond Coeffs", "angle": "Angle Coeffs",
    "dihedral": "Dihedral Coeffs", "improper": "Improper Coeffs",
}
TOPOLOGY_SECTIONS = {
    "bond": "Bonds", "angle": "Angles", "dihedral": "Dihedrals",
    "improper": "Impropers",
}
COUNTS = {"atoms": "atom", "bonds": "bond", "angles": "angle",
          "dihedrals": "dihedral", "impropers": "improper"}
ARITY = {"bond": 2, "angle": 3, "dihedral": 4, "improper": 4}
SUPPORTED_STYLES = {
    "bond": {"harmonic": 2},
    "angle": {"harmonic": 2},
    "dihedral": {"opls": 4},
    "improper": {"harmonic": 2, "cvff": 3},
}
TYPE_LABEL_SECTIONS = {
    "atom": "Atom Type Labels", "bond": "Bond Type Labels",
    "angle": "Angle Type Labels", "dihedral": "Dihedral Type Labels",
    "improper": "Improper Type Labels",
}
SECTIONS = {
    "Masses", "Atoms", *COEFF_SECTIONS.values(),
    *TOPOLOGY_SECTIONS.values(), *TYPE_LABEL_SECTIONS.values(),
}
EXCLUDED_SECTIONS = {
    "PairIJ Coeffs", "BondBond Coeffs", "BondAngle Coeffs",
    "MiddleBondTorsion Coeffs", "EndBondTorsion Coeffs",
    "AngleTorsion Coeffs", "AngleAngleTorsion Coeffs",
    "BondBond13 Coeffs", "AngleAngle Coeffs", "Velocities",
}


def _numbers(words: list[str], context: str) -> tuple[float, ...]:
    try:
        values = tuple(float(v) for v in words)
    except ValueError as exc:
        raise ValueError(f"{context}: malformed numerical values: {words!r}") from exc
    if not all(math.isfinite(v) for v in values):
        raise ValueError(f"{context}: NaN/Infinity in numbers")
    return values


def _equivalent(a: tuple[float, ...], b: tuple[float, ...]) -> bool:
    return len(a) == len(b) and all(
        math.isclose(x, y, rel_tol=1e-8, abs_tol=1e-9) for x, y in zip(a, b)
    )


def _identity(category: str, labels: tuple[str, ...]) -> tuple[str, ...]:
    if category in ("bond", "angle", "dihedral"):
        return min(labels, labels[::-1])
    # Improper order/convention can affect the potential, don't reorder it.
    return labels


@dataclass(frozen=True)
class Parameter:
    values: tuple[float, ...]
    labels: tuple[str, ...]


@dataclass
class ParsedData:
    path: Path
    counts: dict[str, int]
    atom_labels: dict[int, str]
    masses: dict[int, float]
    params: dict[str, dict[int, Parameter]]
    styles: dict[str, str]
    atoms: dict[int, list[str]]
    topology: dict[str, list[list[int]]]


def parse_lammps_data(path: Path) -> ParsedData:
    """Read the mBuild subset of LAMMPS .data; reject unsupported styles."""
    path = Path(path)
    lines = path.read_text().splitlines()
    if not lines or "units = real" not in lines[0].lower():
        raise ValueError(f"{path}: expected an mBuild LAMMPS data file in real units")

    contents: dict[str, list[str]] = {}
    styles: dict[str, str] = {}
    header: list[str] = []
    section: str | None = None
    for raw in lines[1:]:
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        clean = s.split("#", 1)[0].strip()
        if clean in EXCLUDED_SECTIONS:
            raise NotImplementedError(f"{path}: {clean} cannot be merged safely")
        if clean in SECTIONS:
            if clean in contents:
                raise ValueError(f"{path}: duplicate section {clean}")
            section = clean
            contents[clean] = []
            styles[clean] = s.split("#", 1)[1].strip().lower() if "#" in s else ""
            continue
        if section is None:
            header.append(clean)
        else:
            contents[section].append(s)

    if "Masses" not in contents or "Pair Coeffs" not in contents or "Atoms" not in contents:
        raise ValueError(f"{path}: missing Masses, Pair Coeffs, or Atoms section")
    if styles["Atoms"] != "full":
        raise ValueError(f"{path}: only Atoms # full is supported")

    counts = {key: 0 for key in COUNTS}
    ntypes = {key: 0 for key in CATEGORIES}
    for line in header:
        fields = line.split()
        if len(fields) == 2 and fields[0].isdigit() and fields[1] in counts:
            counts[fields[1]] = int(fields[0])
        elif len(fields) == 3 and fields[0].isdigit() and fields[-1] == "types" and fields[1] in ntypes:
            ntypes[fields[1]] = int(fields[0])
        elif len(fields) == 4 and fields[2:] in (["xlo", "xhi"], ["ylo", "yhi"], ["zlo", "zhi"]):
            _numbers(fields[:2], str(path))
        elif len(fields) == 6 and fields[3:] == ["xy", "xz", "yz"]:
            if any(abs(v) > 1e-10 for v in _numbers(fields[:3], str(path))):
                raise NotImplementedError(f"{path}: triclinic box not supported")
        else:
            raise ValueError(f"{path}: unexpected LAMMPS data header line {line!r}")

    # Prefer LAMMPS's explicit Atom Type Labels section, with the mBuild
    # Masses comments as the fallback for original Foyer-typed files.
    declared_labels: dict[int, str] = {}
    for raw in contents.get(TYPE_LABEL_SECTIONS["atom"], []):
        fields = raw.split("#", 1)[0].split()
        if len(fields) != 2:
            raise ValueError(f"{path}: malformed Atom Type Labels entry {raw!r}")
        idx = int(fields[0])
        if idx in declared_labels or fields[1] in declared_labels.values():
            raise ValueError(f"{path}: duplicate Atom Type Labels ID or label")
        declared_labels[idx] = fields[1]
    atom_labels: dict[int, str] = {}
    masses: dict[int, float] = {}
    for raw in contents["Masses"]:
        fields = raw.split("#", 1)
        pair = fields[0].split()
        label = fields[1].strip().split()[0] if len(fields) > 1 and fields[1].strip() else ""
        if len(pair) != 2:
            raise ValueError(f"{path}: malformed Masses row: {raw}")
        idx = int(pair[0])
        declared = declared_labels.get(idx)
        if declared and label and declared != label:
            raise ValueError(f"{path}: inconsistent Masses and Atom Type Labels at ID {idx}")
        label = declared or label
        if not label:
            raise ValueError(f"{path}: atom type {idx} lacks Atom Type Labels and # comment")
        mass = _numbers([pair[1]], str(path))[0]
        if idx in masses or mass <= 0:
            raise ValueError(f"{path}: repeated ID / invalid atom mass at {idx}")
        atom_labels[idx] = label
        masses[idx] = mass
    if set(masses) != set(range(1, ntypes["atom"] + 1)):
        raise ValueError(f"{path}: Masses do not cover 1..{ntypes['atom']} atom types")
    if declared_labels and set(declared_labels) != set(masses):
        raise ValueError(f"{path}: Atom Type Labels are not complete")

    params: dict[str, dict[int, Parameter]] = {cat: {} for cat in CATEGORIES}
    for category, name in COEFF_SECTIONS.items():
        if name not in contents:
            if ntypes[category]:
                raise ValueError(f"{path}: {name} missing for declared {category} types")
            continue
        style = styles[name].split()[0] if styles[name] else ""
        if category != "atom" and style not in SUPPORTED_STYLES[category]:
            raise NotImplementedError(f"{path}: unsupported {name} style: {style!r}")
        if category == "atom" and style not in ("", "lj", "lj/cut", "lj/cut/coul/long"):
            raise NotImplementedError(f"{path}: unsupported pair style annotation {style!r}")
        for raw in contents[name]:
            values_comment = raw.split("#", 1)
            tokens = values_comment[0].split()
            if not tokens:
                continue
            idx = int(tokens[0])
            coeff = _numbers(tokens[1:], f"{path}/{name}/{idx}")
            needed = 2 if category == "atom" else SUPPORTED_STYLES[category][style]
            if len(coeff) != needed or idx in params[category]:
                raise ValueError(f"{path}: malformed or duplicate {name} ID {idx}")
            labels = tuple(values_comment[1].split()) if len(values_comment) > 1 else ()
            if category == "atom" and labels and labels[0] != atom_labels.get(idx):
                raise ValueError(f"{path}: Pair Coeffs type {idx} label disagrees with Masses")
            params[category][idx] = Parameter(coeff, labels)
        if set(params[category]) != set(range(1, ntypes[category] + 1)):
            raise ValueError(f"{path}: incomplete {name} coefficient list")

    atoms: dict[int, list[str]] = {}
    for raw in contents["Atoms"]:
        fields = raw.split("#", 1)[0].split()
        if len(fields) not in (7, 10):
            raise ValueError(f"{path}: unsupported atom row {raw!r}")
        idx, typ = int(fields[0]), int(fields[2])
        if idx in atoms or typ not in masses or typ not in params["atom"]:
            raise ValueError(f"{path}: duplicate atom or undefined atom type {idx}/{typ}")
        _numbers(fields[3:7], f"{path}/atom {idx}")
        atoms[idx] = fields
    if set(atoms) != set(range(1, counts["atoms"] + 1)):
        raise ValueError(f"{path}: Atoms count/indices invalid")

    topology: dict[str, list[list[int]]] = {}
    for category, name in TOPOLOGY_SECTIONS.items():
        topology[category] = []
        for raw in contents.get(name, []):
            fields = raw.split("#", 1)[0].split()
            if len(fields) != 2 + ARITY[category]:
                raise ValueError(f"{path}: malformed {name} row {raw!r}")
            row = list(map(int, fields))
            if row[1] not in params[category] or any(a not in atoms for a in row[2:]):
                raise ValueError(f"{path}: {name} references undefined atom or coefficient")
            topology[category].append(row)
        if len(topology[category]) != counts[name.lower()]:
            raise ValueError(f"{path}: {name} header count mismatch")
        if {row[0] for row in topology[category]} != set(range(1, len(topology[category]) + 1)):
            raise ValueError(f"{path}: {name} IDs must be 1..N")

        # Original mBuild coefficient comments carry atom-label signatures.
        # If missing, infer signature from the topology using Atom Type Labels.
        for typ, par in list(params[category].items()):
            if len(par.labels) == ARITY[category]:
                continue
            used = {
                _identity(category, tuple(atom_labels[int(atoms[atomid][2])] for atomid in row[2:]))
                for row in topology[category] if row[1] == typ
            }
            if len(used) != 1:
                raise ValueError(
                    f"{path}: {name} type {typ} lacks usable # atom-type labels "
                    f"and topology cannot infer one unique signature ({used!r})"
                )
            params[category][typ] = Parameter(par.values, next(iter(used)))
    return ParsedData(path, counts, atom_labels, masses, params, styles, atoms, topology)


class LammpsMergeBuilder:
    """Compile globally consistent force-field and molecule/template file IDs."""

    def __init__(self, merge_dir: Path):
        self.merge_dir = Path(merge_dir)
        self.molecule_dir = self.merge_dir / "molecules"
        self.template_dir = self.merge_dir / "templates"
        for directory in (self.merge_dir, self.molecule_dir, self.template_dir):
            directory.mkdir(parents=True, exist_ok=True)

    def build(self, typing_results: Iterable[object]) -> FFFiles:
        results = list(typing_results)
        if not results:
            raise ValueError("No parameterized structures to merge")
        parsed = [(result, parse_lammps_data(result.typed_data_file)) for result in results]
        atom_ids: dict[str, int] = {}
        atom_info: dict[str, tuple[float, tuple[float, ...]]] = {}
        bonded_ids: dict[str, dict[tuple[str, tuple[str, ...]], int]] = {
            cat: {} for cat in BONDED
        }
        bonded_params: dict[str, list[Parameter]] = {cat: [] for cat in BONDED}
        styles: dict[str, str] = {}
        per_molecule_maps: dict[str, dict[str, dict[int, int]]] = {}

        for result, doc in parsed:
            if result.id in per_molecule_maps:
                raise ValueError(f"Duplicate molecule/template identifier {result.id!r}")
            maps: dict[str, dict[int, int]] = {cat: {} for cat in CATEGORIES}
            for local, label in sorted(doc.atom_labels.items()):
                details = (doc.masses[local], doc.params["atom"][local].values)
                if label not in atom_ids:
                    atom_ids[label] = len(atom_ids) + 1
                    atom_info[label] = details
                else:
                    old = atom_info[label]
                    if not math.isclose(old[0], details[0], rel_tol=1e-8, abs_tol=1e-7) or not _equivalent(old[1], details[1]):
                        raise ValueError(f"{doc.path}: conflicting mass/pair params for atom label {label!r}")
                maps["atom"][local] = atom_ids[label]

            for cat in BONDED:
                section = COEFF_SECTIONS[cat]
                if not doc.params[cat]:
                    continue
                style = doc.styles[section].split()[0]
                if cat in styles and styles[cat] != style:
                    raise NotImplementedError(
                        f"Conflicting {cat} styles: {styles[cat]} vs {style}"
                    )
                styles[cat] = style
                for local, par in sorted(doc.params[cat].items()):
                    labels = _identity(cat, par.labels)
                    identity = (style, labels)
                    if identity not in bonded_ids[cat]:
                        bonded_ids[cat][identity] = len(bonded_ids[cat]) + 1
                        bonded_params[cat].append(Parameter(par.values, labels))
                    else:
                        gid = bonded_ids[cat][identity]
                        if not _equivalent(bonded_params[cat][gid - 1].values, par.values):
                            raise ValueError(
                                f"{doc.path}: {cat} {labels} has incompatible parameters "
                                "between typed files"
                            )
                    maps[cat][local] = bonded_ids[cat][identity]
            per_molecule_maps[result.id] = maps

        global_atoms = {atom_ids[label]: (label, *atom_info[label]) for label in atom_ids}
        force_field_data = self.merge_dir / "force_field.data"
        self._write_data(force_field_data, {}, {cat: [] for cat in BONDED}, global_atoms,
                         bonded_params, styles, coefficient_only=True)

        molecule_files: list[MoleculeFile] = []
        reaction_pairs: dict[int, dict[str, DataFiles]] = defaultdict(dict)
        for result, doc in parsed:
            maps = per_molecule_maps[result.id]
            atoms: dict[int, list[str]] = {}
            for i, local in doc.atoms.items():
                fields = list(local)
                fields[2] = str(maps["atom"][int(fields[2])])
                atoms[i] = fields
            interactions: dict[str, list[list[int]]] = {cat: [] for cat in BONDED}
            for cat in BONDED:
                for row in doc.topology[cat]:
                    r = list(row)
                    r[1] = maps[cat][r[1]]
                    interactions[cat].append(r)

            target = self.molecule_dir if result.molecule else self.template_dir
            data_file = target / f"{result.id}_typed_merged.data"
            molecule_path = target / f"{result.id}.molecule"
            self._write_data(data_file, atoms, interactions, global_atoms, bonded_params, styles)
            self._write_molecule(molecule_path, result.id, atoms, interactions)
            bundle = DataFiles(data_file=data_file, lmp_molecule_file=molecule_path)
            if result.molecule:
                molecule_files.append(MoleculeFile(id=result.id, molecule_files=bundle))
            else:
                match = re.fullmatch(r"(pre|post)(\d+)", str(result.id))
                if not match:
                    raise ValueError(f"Invalid reaction template ID {result.id!r}; expected preN/postN")
                prefix, reaction_id = match.group(1), int(match.group(2))
                if prefix in reaction_pairs[reaction_id]:
                    raise ValueError(f"Duplicate {prefix} template for reaction {reaction_id}")
                reaction_pairs[reaction_id][prefix] = bundle

        templates = []
        for reaction_id, entries in sorted(reaction_pairs.items()):
            if set(entries) != {"pre", "post"}:
                raise ValueError(f"Reaction {reaction_id}: requires both pre and post templates")
            templates.append(TemplateFile(
                reaction_id=reaction_id,
                pre_reaction_file=entries["pre"],
                post_reaction_file=entries["post"],
            ))

        # No in.create_atoms.script: the project handles LAMMPS runtime inputs.
        # Remove a stale script from an earlier version of THIS builder.
        stale_script = self.merge_dir / "in.create_atoms.script"
        if stale_script.is_file():
            stale_script.unlink()
        final = FFFiles(
            molecule_files=molecule_files, template_files=templates,
            force_field_data=force_field_data, in_file=None,
        )
        print("[LAMMPS merge] Global types: " + ", ".join(
            f"{cat}={len(atom_ids) if cat == 'atom' else len(bonded_params[cat])}"
            for cat in CATEGORIES
        ))
        print(f"[LAMMPS merge] Wrote {len(molecule_files)} monomer(s), {len(templates)} reaction(s) -> {self.merge_dir}")
        return final

    @staticmethod
    def _write_data(path: Path, atoms: dict[int, list[str]],
                    topology: dict[str, list[list[int]]],
                    atom_params: dict[int, tuple[str, float, tuple[float, ...]]],
                    bonded_params: dict[str, list[Parameter]], styles: dict[str, str],
                    coefficient_only: bool = False) -> None:
        lines = [f"AutoREACTER global types: {path.name}; units = real", ""]
        lines.append(f"{len(atoms)} atoms")
        for cat in BONDED:
            lines.append(f"{len(topology[cat])} {TOPOLOGY_SECTIONS[cat].lower()}")
        lines.append("")
        lines.append(f"{len(atom_params)} atom types")
        for cat in BONDED:
            lines.append(f"{len(bonded_params[cat])} {cat} types")
        lines.append("")

        if coefficient_only:
            coords = [(-100.0, 100.0)] * 3
        else:
            positions = [tuple(float(f) for f in row[4:7]) for row in atoms.values()]
            coords = [(min(p[i] for p in positions) - 10.0,
                       max(p[i] for p in positions) + 10.0) for i in range(3)]
        for (lo, hi), axis in zip(coords, "xyz"):
            lines.append(f"{lo:.8f} {hi:.8f} {axis}lo {axis}hi")
        # LAMMPS type-label maps for all global type IDs, including bonded types.
        # Labels are metadata; numeric IDs are still used in .lmpmol files.
        lines.extend(["", TYPE_LABEL_SECTIONS["atom"], ""])
        for idx, (label, _, _) in sorted(atom_params.items()):
            if not label or label[0].isdigit() or any(c.isspace() for c in label) or "#" in label:
                raise ValueError(f"Invalid LAMMPS atom type label: {label!r}")
            lines.append(f"{idx} {label}")
        for cat in BONDED:
            params = bonded_params[cat]
            if params:
                lines.extend(["", TYPE_LABEL_SECTIONS[cat], ""])
                used_type_labels = set()
                for idx, par in enumerate(params, 1):
                    # Use the actual atom-type signature as a readable label.
                    # Keep numeric IDs in topology, to avoid relying on label maps
                    # when a molecule template is loaded in isolation.
                    type_label = "-".join(par.labels)
                    if (not type_label or type_label[0].isdigit()
                        or any(c.isspace() for c in type_label) or "#" in type_label):
                        raise ValueError(f"Invalid {cat} type label: {type_label!r}")
                    if type_label in used_type_labels:
                        raise ValueError(f"Colliding {cat} labels: {type_label!r}")
                    used_type_labels.add(type_label)
                    lines.append(f"{idx} {type_label}")
        lines.extend(["", "Masses", ""])
        for idx, (label, mass, _) in sorted(atom_params.items()):
            lines.append(f"{idx} {mass:.9f} # {label}")
        lines.extend(["", "Pair Coeffs", ""])
        for idx, (label, _, pair) in sorted(atom_params.items()):
            lines.append(f"{idx} {' '.join(f'{v:.10g}' for v in pair)} # {label}")
        for cat in BONDED:
            params = bonded_params[cat]
            if not params:
                continue
            lines.extend(["", f"{COEFF_SECTIONS[cat]} # {styles[cat]}", ""])
            for idx, p in enumerate(params, 1):
                lines.append(f"{idx} {' '.join(f'{v:.10g}' for v in p.values)} # {' '.join(p.labels)}")
        if atoms:
            lines.extend(["", "Atoms # full", ""])
            lines.extend(" ".join(atoms[idx]) for idx in sorted(atoms))
        for cat in BONDED:
            if topology[cat]:
                lines.extend(["", TOPOLOGY_SECTIONS[cat], ""])
                lines.extend(" ".join(map(str, row)) for row in topology[cat])
        path.write_text("\n".join(lines) + "\n")

    @staticmethod
    def _write_molecule(path: Path, label: str, atoms: dict[int, list[str]],
                        topology: dict[str, list[list[int]]]) -> None:
        lines = [f"# AutoREACTER {label}: global numeric type IDs", "",
                 f"{len(atoms)} atoms"]
        for cat in BONDED:
            lines.append(f"{len(topology[cat])} {TOPOLOGY_SECTIONS[cat].lower()}")
        for section, indices in (("Types", (2,)), ("Charges", (3,)),
                                 ("Coords", (4, 5, 6))):
            lines.extend(["", section, ""])
            for idx, fields in sorted(atoms.items()):
                lines.append(" ".join([str(idx)] + [fields[i] for i in indices]))
        # REACTERFilesBuilder expects all section headings even with zero items.
        for cat in BONDED:
            lines.extend(["", TOPOLOGY_SECTIONS[cat], ""])
            lines.extend(" ".join(map(str, row)) for row in topology[cat])
        path.write_text("\n".join(lines) + "\n")


# Backward compatibility for the existing Foyer backend.
FoyerMergeBuilder = LammpsMergeBuilder
