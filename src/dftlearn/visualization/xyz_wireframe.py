"""2D XY projections of XYZ geometries with distance-based bonds (no RDKit draw)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from matplotlib.patches import Circle

from dftlearn.io.xyz_structure import element_symbol_from_xyz_label

if TYPE_CHECKING:
    from matplotlib.axes import Axes

_COVALENT_RADII: dict[str, float] = {
    "H": 0.31,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "S": 1.05,
    "P": 1.07,
    "F": 0.57,
    "Cl": 0.99,
    "Br": 1.14,
    "I": 1.33,
    "Zn": 1.22,
    "Cu": 1.32,
    "Fe": 1.25,
    "Ni": 1.24,
    "Mn": 1.39,
    "Co": 1.26,
    "Mg": 1.30,
    "Ca": 1.74,
    "Na": 1.66,
    "K": 2.03,
    "Si": 1.11,
}

_DEFAULT_RADIUS = 0.75

_METALS = frozenset(
    {
        "Zn",
        "Cu",
        "Fe",
        "Ni",
        "Mn",
        "Co",
        "Mg",
        "Ca",
        "Na",
        "K",
        "Li",
        "Pd",
        "Pt",
        "Ag",
        "Au",
        "Cr",
        "Mo",
        "W",
        "Ti",
        "V",
    }
)

_LIGAND_LIGHT = frozenset({"N", "O", "S", "P", "F", "Cl", "Br", "I", "C"})


def _radius(sym: str) -> float:
    return _COVALENT_RADII.get(sym, _DEFAULT_RADIUS)


def infer_bonds_from_xyz_rows(
    rows: list[tuple[str, float, float, float]],
    tolerance: float = 0.32,
    metal_ligand_extra: float = 0.55,
) -> list[tuple[int, int]]:
    """Infer covalent bonds from 3D distances and covalent radii.

    Metal-ligand pairs (metal to N, O, S, P, halogen, or C) use an extra slack
    so Zn-N and similar dative contacts are drawn when the metal sits slightly
    beyond typical covalent sums.

    Parameters
    ----------
    rows : list[tuple[str, float, float, float]]
        ``(label, x, y, z)`` in angstroms, file order defines atom indices.
    tolerance : float, optional
        Extra Å added to every summed covalent radii cutoff.
    metal_ligand_extra : float, optional
        Additional Å for one metal and one ligand/heavy neighbor pair.

    Returns
    -------
    list[tuple[int, int]]
        Unique bonds as ``(i, j)`` with ``i < j``.
    """
    n = len(rows)
    if n < 2:
        return []
    symbols = [element_symbol_from_xyz_label(lab) for lab, _, _, _ in rows]
    pos = np.array([[r[1], r[2], r[3]] for r in rows], dtype=np.float64)
    bonds: list[tuple[int, int]] = []
    for i in range(n):
        for j in range(i + 1, n):
            a, b = symbols[i], symbols[j]
            ra, rb = _radius(a), _radius(b)
            cutoff = ra + rb + tolerance
            if (a in _METALS and b in _LIGAND_LIGHT) or (
                b in _METALS and a in _LIGAND_LIGHT
            ):
                cutoff += metal_ligand_extra
            dist = float(np.linalg.norm(pos[i] - pos[j]))
            if dist <= cutoff:
                bonds.append((i, j))
    return bonds


_ELEMENT_TEXT_COLOR: dict[str, str] = {
    "N": "#1565c0",
    "O": "#c62828",
    "S": "#9e8600",
    "P": "#ef6c00",
    "F": "#2e7d32",
    "Cl": "#2e7d32",
    "Br": "#5d4037",
    "I": "#6a1b9a",
    "Zn": "#424242",
    "Cu": "#5d4037",
    "Fe": "#6d4c41",
    "H": "#212121",
}


def draw_xyz_wireframe_on_ax(
    ax: Axes,
    rows: list[tuple[str, float, float, float]],
    site_atom_colors: dict[int, tuple[float, float, float]],
    *,
    show_hydrogen: bool = False,
    bond_color: str = "#333333",
    bond_lw: float = 1.35,
    bond_alpha: float = 0.9,
    halo_radius_angstrom: float = 0.52,
    halo_alpha: float = 0.5,
    halo_soft_edge: bool = True,
    label_fontsize: int = 9,
    hetero_label_white_pad: float | None = None,
    plot_margins: float = 0.12,
) -> None:
    """Draw a skeletal XY wireframe: bonds only, no atom markers.

    Carbons are implicit (no dots, no ``C`` labels). Heteroatoms are labeled with
    colored element symbols only (no circular markers). Core-excitation sites
    listed in ``site_atom_colors`` get a soft filled disk in data coordinates,
    using the same RGB as the spectrum trace (typically low alpha ~0.5).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes; caller should disable spines/ticks as needed.
    rows : list[tuple[str, float, float, float]]
        Atom rows from :func:`dftlearn.io.xyz_structure.xyz_rows_from_file`.
    site_atom_colors : dict[int, tuple[float, float, float]]
        0-based atom index to RGB in 0-1 for each core site to shade.
    show_hydrogen : bool, optional
        If False (default), hydrogens and their bonds are omitted.
    bond_color : str, optional
        Bond line color.
    bond_lw : float, optional
        Bond line width.
    bond_alpha : float, optional
        Bond line alpha.
    halo_radius_angstrom : float, optional
        Disk radius in Å for each site halo (XY projection).
    halo_alpha : float, optional
        Face alpha for site shading (default ~50%).
    halo_soft_edge : bool, optional
        If True, add a larger, weaker outer disk for a softer falloff.
    label_fontsize : int, optional
        Font size for heteroatom symbols.
    hetero_label_white_pad : float, optional
        If set (default None uses a small pad when drawing Zn/N), ``Zn`` and ``N``
        labels get a tight white box behind the text (no border).
    plot_margins : float, optional
        Matplotlib ``margins`` on x and y after ``set_aspect``; use a small value
        (for example ``0.02``) so the structure fills an inset axes.
    """
    bonds = infer_bonds_from_xyz_rows(rows)
    symbols = [element_symbol_from_xyz_label(lab) for lab, _, _, _ in rows]

    visible: list[int] = []
    for i, sym in enumerate(symbols):
        if sym == "H" and not show_hydrogen:
            continue
        visible.append(i)

    for i, j in bonds:
        if i not in visible or j not in visible:
            continue
        x1, y1 = rows[i][1], rows[i][2]
        x2, y2 = rows[j][1], rows[j][2]
        ax.plot(
            [x1, x2],
            [y1, y2],
            color=bond_color,
            linewidth=bond_lw,
            alpha=bond_alpha,
            solid_capstyle="round",
            zorder=1,
        )

    for atom_idx, rgb in site_atom_colors.items():
        if atom_idx < 0 or atom_idx >= len(rows):
            continue
        if atom_idx not in visible:
            continue
        x, y = rows[atom_idx][1], rows[atom_idx][2]
        if halo_soft_edge:
            outer = Circle(
                (x, y),
                radius=halo_radius_angstrom * 1.45,
                facecolor=rgb,
                edgecolor="none",
                alpha=halo_alpha * 0.35,
                zorder=2,
            )
            ax.add_patch(outer)
        circ = Circle(
            (x, y),
            radius=halo_radius_angstrom,
            facecolor=rgb,
            edgecolor="none",
            alpha=halo_alpha,
            zorder=3,
        )
        ax.add_patch(circ)

    for i in visible:
        sym = symbols[i]
        x, y = rows[i][1], rows[i][2]
        if sym == "C":
            continue
        if sym == "H" and show_hydrogen:
            ax.text(
                x,
                y,
                "H",
                ha="center",
                va="center",
                fontsize=label_fontsize - 1,
                fontweight="bold",
                color=_ELEMENT_TEXT_COLOR["H"],
                zorder=6,
            )
            continue
        if sym == "H":
            continue
        tcol = _ELEMENT_TEXT_COLOR.get(sym, "#424242")
        bbox_kw: dict | None = None
        if sym in ("Zn", "N"):
            pad = 0.08 if hetero_label_white_pad is None else hetero_label_white_pad
            bbox_kw = {
                "boxstyle": f"round,pad={pad}",
                "facecolor": "white",
                "edgecolor": "none",
                "linewidth": 0,
                "alpha": 0.92,
            }
        ax.text(
            x,
            y,
            sym,
            ha="center",
            va="center",
            fontsize=label_fontsize,
            fontweight="bold",
            color=tcol,
            zorder=6,
            bbox=bbox_kw,
        )

    ax.set_aspect("equal")
    ax.margins(plot_margins)
