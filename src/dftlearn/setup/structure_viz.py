"""3D structure helpers for setup-lab visualization and bond handling."""

from __future__ import annotations

from contextlib import suppress
from typing import TYPE_CHECKING

from rdkit import Chem

from dftlearn.io.xyz_structure import element_symbol_from_xyz_label
from dftlearn.setup.distinct_atoms import mol_from_xyz_rows

if TYPE_CHECKING:
    from dftlearn.setup.labeling import SiteGroup
    from dftlearn.setup.session import SetupSession

Row = tuple[str, float, float, float]


def molblock_from_rows(
    rows: list[Row],
    mol: Chem.Mol | None = None,
) -> str:
    """Return an MDL mol block with coordinates taken from ``rows``.

    Uses ``mol`` connectivity when supplied; otherwise infers bonds from distances.
    """
    if mol is None:
        working = mol_from_xyz_rows(list(rows), sanitize=False)
    else:
        working = Chem.Mol(mol)
        if working.GetNumConformers() == 0:
            conf = Chem.Conformer(working.GetNumAtoms())
            working.AddConformer(conf, assignId=True)
        conf = working.GetConformer()
        for idx, (_label, x, y, z) in enumerate(rows):
            conf.SetAtomPosition(idx, (float(x), float(y), float(z)))
    return Chem.MolToMolBlock(working)


def bonds_from_mol(mol: Chem.Mol) -> list[tuple[int, int]]:
    """List unique bonded atom-index pairs from an RDKit molecule."""
    pairs: list[tuple[int, int]] = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if i > j:
            i, j = j, i
        pairs.append((i, j))
    pairs.sort()
    return pairs


def display_bonds(
    rows: list[Row],
    mol: Chem.Mol | None,
) -> list[tuple[int, int]]:
    """Return bond pairs for 3D visualization.

    Prefers RDKit connectivity from the embedded structure. When that graph is
    empty, falls back to valence-limited distance inference.
    """
    if mol is not None and mol.GetNumBonds() > 0:
        return bonds_from_mol(mol)

    inferred = mol_from_xyz_rows(list(rows), sanitize=False)
    return bonds_from_mol(inferred)


def bond_neighbors(mol: Chem.Mol, atom_index: int) -> list[int]:
    """Return sorted neighbor atom indices bonded to ``atom_index``."""
    atom = mol.GetAtomWithIdx(atom_index)
    return sorted(n.GetIdx() for n in atom.GetNeighbors())


def format_bond_label(rows: list[Row], i: int, j: int) -> str:
    """Format a human-readable bond label such as ``C01-C02 (3-7)``."""
    left = rows[i][0]
    right = rows[j][0]
    return f"{left}-{right} ({i}-{j})"


def minimal_generating_indices(site_groups: list[SiteGroup]) -> tuple[int, ...]:
    """Return representative atom indices for enabled site groups."""
    return tuple(
        group.representative_index
        for group in sorted(site_groups, key=lambda g: g.group_id)
        if group.enabled
    )


def session_mol(session: SetupSession) -> Chem.Mol:
    """Rebuild or load the RDKit molecule associated with a setup session."""
    rows = [
        (str(r[0]), float(r[1]), float(r[2]), float(r[3])) for r in session.base_rows
    ]
    mol_block = str(session.meta.get("mol_block", "")).strip()
    if mol_block:
        mol = Chem.MolFromMolBlock(mol_block, removeHs=False, sanitize=False)
        if mol is not None:
            with suppress(Chem.AtomValenceException):
                Chem.SanitizeMol(mol)
            return mol
    return mol_from_xyz_rows(rows, sanitize=False)


def store_mol_block(session: SetupSession, mol: Chem.Mol) -> None:
    """Persist ``mol`` connectivity into ``session.meta['mol_block']``."""
    session.meta["mol_block"] = Chem.MolToMolBlock(mol)


def atom_display_table(
    rows: list[Row],
    site_groups: list[SiteGroup],
) -> list[dict[str, object]]:
    """Build table rows for atom selection in the setup lab."""
    index_to_tag: dict[int, str] = {}
    for group in site_groups:
        for idx in group.atom_indices:
            tag = group.site_tag if group.enabled else f"{group.site_tag}*"
            index_to_tag[idx] = tag

    table: list[dict[str, object]] = []
    for idx, (label, x, y, z) in enumerate(rows):
        sym = element_symbol_from_xyz_label(label)
        table.append(
            {
                "index": idx,
                "label": label,
                "element": sym,
                "x": round(x, 3),
                "y": round(y, 3),
                "z": round(z, 3),
                "site": index_to_tag.get(idx, ""),
            }
        )
    return table
