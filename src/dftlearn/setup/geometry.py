"""StoBe XYZ ordering: generating core sites first, then element blocks."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import TYPE_CHECKING

from dftlearn.io.xyz_structure import element_symbol_from_xyz_label
from dftlearn.setup.types import InitGeometryPlan

if TYPE_CHECKING:
    from dftlearn.setup.types import DistinctAtomGroup


def _element_block_priority(edge_element: str) -> dict[str, int]:
    edge = edge_element.strip().capitalize()
    order = [edge, "H", "N", "O", "S", "P", "F", "Cl", "Br", "I", "Al", "Zn"]
    return {sym: idx for idx, sym in enumerate(order)}


def plan_stobe_geometry(
    rows: list[tuple[str, float, float, float]],
    *,
    edge_element: str,
    distinct_groups: list[DistinctAtomGroup],
) -> InitGeometryPlan:
    """Reorder XYZ rows so generating core sites precede each element block.

    The first ``len(distinct_groups)`` atoms of ``edge_element`` become StoBe
    site folders ``C1``, ``C2``, ... (or ``N1``, etc.). Remaining atoms follow
    in element-block order matching ``molConfig`` ``nElem*`` groups.

    Parameters
    ----------
    rows
        Initial XYZ rows in arbitrary order.
    edge_element
        Core-edge element symbol (``C`` for carbon K-edge).
    distinct_groups
        Distinct equivalence classes; representatives are placed first.

    Returns
    -------
    InitGeometryPlan
        Reordered rows, generating indices, counts, and element-block order.
    """
    edge = edge_element.strip().capitalize()
    n_atoms = len(rows)
    used = [False] * n_atoms

    generating: list[int] = []
    for group in distinct_groups:
        rep = group.representative_index
        if rep < 0 or rep >= n_atoms:
            msg = f"Representative index {rep} out of range"
            raise ValueError(msg)
        generating.append(rep)
        used[rep] = True

    remaining_by_element: dict[str, list[int]] = {}
    for idx, (label, _x, _y, _z) in enumerate(rows):
        if used[idx]:
            continue
        sym = element_symbol_from_xyz_label(label)
        remaining_by_element.setdefault(sym, []).append(idx)

    priority = _element_block_priority(edge)
    element_symbols = {element_symbol_from_xyz_label(r[0]) for r in rows}
    element_group_order = sorted(
        element_symbols,
        key=lambda sym: (priority.get(sym, 999), sym),
    )

    ordered_indices: list[int] = list(generating)
    for sym in element_group_order:
        for idx in remaining_by_element.get(sym, []):
            if idx not in ordered_indices:
                ordered_indices.append(idx)

    if len(ordered_indices) != n_atoms:
        msg = "Geometry reorder dropped atoms"
        raise ValueError(msg)

    counts = Counter(element_symbol_from_xyz_label(label) for label, *_ in rows)
    new_rows: list[tuple[str, float, float, float]] = []
    new_generating: list[int] = []
    elem_counts: Counter[str] = Counter()
    for new_idx, old_idx in enumerate(ordered_indices):
        label, x, y, z = rows[old_idx]
        sym = element_symbol_from_xyz_label(label)
        elem_counts[sym] += 1
        new_label = f"{sym}{elem_counts[sym]:02d}"
        new_rows.append((new_label, x, y, z))
        if old_idx in generating:
            new_generating.append(new_idx)

    return InitGeometryPlan(
        rows=tuple(new_rows),
        edge_element=edge,
        generating_site_indices=tuple(new_generating),
        distinct_groups=tuple(distinct_groups),
        element_counts=dict(counts),
        element_group_order=tuple(element_group_order),
    )


def write_stobe_xyz(
    path: Path,
    plan: InitGeometryPlan,
    *,
    comment: str = "",
) -> None:
    """Write a headerless StoBe-style XYZ file from a geometry plan.

    Parameters
    ----------
    path
        Output ``.xyz`` path.
    plan
        Ordered geometry from :func:`plan_stobe_geometry`.
    comment
        Reserved for future standard XYZ headers; ignored for StoBe output.
    """
    path = Path(path)
    lines = [f"{label} {x:.10f} {y:.10f} {z:.10f}" for label, x, y, z in plan.rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
