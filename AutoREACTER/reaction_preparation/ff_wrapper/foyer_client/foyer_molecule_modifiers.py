"""Standard LAMMPS molecule-file modifiers for AutoREACTER's Foyer pathway.

Companion to foyer_merge_builder.py.  The LUNAR modifiers expect extended
nonstandard rows (``ID type # hash real_type``), whereas LAMMPS .molecule
files use standard ``ID type``, ``ID charge``, and so on.  This module exposes
the SAME six functions and return signatures as modifiers_molecule_files.py,
without relying on synthetic comment fields or changing the LUNAR parser.

All atom/template index filtering is 1-based.  Bonded interaction numbers
are global type IDs allocated by FoyerMergeBuilder, and are not renumbered.
"""
from __future__ import annotations

import math
import pandas as pd


SECTION_HEADINGS = {
    "Types", "Charges", "Coords", "Bonds", "Angles", "Dihedrals", "Impropers",
    "Molecules", "Fragments", "Masses", "Special Bond Counts", "Special Bonds",
}


def _rows(lines: list[str], start: int | None, ncols: int, section: str):
    if start is None:
        return []
    result = []
    for raw in lines[start:]:
        clean = raw.split('#', 1)[0].strip()
        if not clean:
            break
        # A zero-count section has no data lines. When the caller begins at
        # heading_index+2, the next section heading can be the first line.
        if clean in SECTION_HEADINGS:
            break
        fields = clean.split()
        if len(fields) != ncols:
            raise ValueError(f"{section}: expected {ncols} fields; got {raw!r}")
        result.append(fields)
    return result


def _type(value: str):
    try:
        return int(value)
    except ValueError:
        return value


def _map(df: pd.DataFrame) -> dict[int, int]:
    return dict(zip(df['atom_index'].astype(int), df['new_atom_index'].astype(int)))


def modify_types(lines, template_indexes, type_start_index):
    original = _rows(lines, type_start_index, 2, 'Types')
    wanted = set(map(int, template_indexes))
    original_ids = [int(row[0]) for row in original]
    if len(set(original_ids)) != len(original_ids):
        raise ValueError('Types contains duplicate atom IDs')
    available = set(original_ids)
    if not wanted <= available:
        raise ValueError(f'Types missing selected atoms: {sorted(wanted - available)}')
    selected = sorted((int(i), _type(typ)) for i, typ in original if int(i) in wanted)
    if not selected:
        raise ValueError('No atoms selected for reaction template')
    df = pd.DataFrame([
        {'atom_index': i, 'atom_type': t, 'new_atom_index': new}
        for new, (i, t) in enumerate(selected, 1)
    ])
    fmt = ''.join(f'{int(row.new_atom_index)} {row.atom_type}\n' for row in df.itertuples())
    index_changes = _map(df)
    numeric_types = all(isinstance(t, int) for _, t in selected)
    return df, fmt, len(selected), index_changes, numeric_types


def _atom_properties(lines, start, df, cols, section):
    amap = _map(df)
    rows = _rows(lines, start, 1 + len(cols), section)
    records = {}
    for fields in rows:
        atom_id = int(fields[0])
        if atom_id in records:
            raise ValueError(f'{section}: duplicate atom {atom_id}')
        if atom_id in amap:
            vals = [float(value) for value in fields[1:]]
            if not all(math.isfinite(v) for v in vals):
                raise ValueError(f'{section}: invalid numeric value at atom {atom_id}')
            records[amap[atom_id]] = vals
    if set(records) != set(amap.values()):
        raise ValueError(f'{section}: selected atoms without {section} values')
    return ''.join(
        f"{new} " + ' '.join(f'{value:.10g}' for value in records[new]) + '\n'
        for new in sorted(records)
    )


def modify_charges(lines, type_df, charge_start_index):
    return _atom_properties(lines, charge_start_index, type_df, ('q',), 'Charges')


def modify_coords(lines, type_df, coord_start_index):
    return _atom_properties(lines, coord_start_index, type_df, ('x', 'y', 'z'), 'Coords')


def _topology(lines, start, df, section, arity):
    mapping = _map(df)
    result = []
    seen_ids = set()
    for fields in _rows(lines, start, 2 + arity, section):
        old_id = int(fields[0])
        if old_id in seen_ids:
            raise ValueError(f'{section}: duplicate interaction ID {old_id}')
        seen_ids.add(old_id)
        typ = _type(fields[1])
        atom_ids = [int(value) for value in fields[2:]]
        if len(set(atom_ids)) != len(atom_ids):
            raise ValueError(f'{section}: interaction {old_id} contains repeated atoms')
        if all(atom in mapping for atom in atom_ids):
            mapped = [mapping[atom] for atom in atom_ids]
            result.append(f"{len(result) + 1} {typ} {' '.join(map(str, mapped))}\n")
    return ''.join(result), len(result)


def modify_bonds(lines, type_df, bond_start_index, legacy_mode):
    return _topology(lines, bond_start_index, type_df, 'Bonds', 2)


def modify_angles(lines, type_df, angle_start_index, legacy_mode):
    return _topology(lines, angle_start_index, type_df, 'Angles', 3)


def modify_dihedrals(lines, type_df, dihedral_start_index, legacy_mode):
    return _topology(lines, dihedral_start_index, type_df, 'Dihedrals', 4)


def modify_impropers(lines, type_df, improper_start_index, legacy_mode):
    return _topology(lines, improper_start_index, type_df, 'Impropers', 4)
