"""Preview figures for calculation initialization."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.transforms import Bbox

from dftlearn.visualization.xyz_wireframe import draw_xyz_wireframe_on_ax

if TYPE_CHECKING:
    from dftlearn.setup.types import CifStructureMeta, InitGeometryPlan


def write_init_preview_figure(
    path: Path,
    *,
    plan: InitGeometryPlan,
    meta: CifStructureMeta,
    title: str,
) -> None:
    """Write a two-panel init preview: 3D wireframe and distinct-site table.

    Parameters
    ----------
    path
        Output PNG path.
    plan
        Ordered geometry with generating core-site indices.
    meta
        CIF metadata for symmetry annotation.
    title
        Figure suptitle (molecule name).
    """
    path = Path(path)
    rows = list(plan.rows)
    n_sites = len(plan.generating_site_indices)
    cmap = plt.get_cmap("tab10")
    palette = cmap(np.linspace(0, 1, max(n_sites, 1)))
    site_colors: dict[int, tuple[float, float, float]] = {}
    for site_idx, gen_idx in enumerate(plan.generating_site_indices):
        rgb = palette[site_idx % len(palette)]
        site_colors[gen_idx] = (float(rgb[0]), float(rgb[1]), float(rgb[2]))

    fig = plt.figure(figsize=(10.5, 5.2), layout="constrained")
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1.0])
    ax_geom = fig.add_subplot(gs[0, 0])
    ax_tbl = fig.add_subplot(gs[0, 1])
    ax_tbl.axis("off")

    draw_xyz_wireframe_on_ax(
        ax_geom,
        rows,
        site_colors,
        show_hydrogen=False,
        plot_margins=0.06,
    )
    ax_geom.set_title("Generating core sites (XY projection)")
    ax_geom.set_xlabel("x (A)")
    ax_geom.set_ylabel("y (A)")

    sg = meta.space_group or "unknown"
    cell_txt = (
        f"a={meta.cell_a:.3f} b={meta.cell_b:.3f} c={meta.cell_c:.3f} A\n"
        f"alpha={meta.cell_alpha:.2f} beta={meta.cell_beta:.2f} "
        f"gamma={meta.cell_gamma:.2f} deg"
    )
    sym_txt = f"Space group: {sg}\n{cell_txt}\nSymmetry ops: {len(meta.symmetry_ops)}"

    table_rows: list[list[str]] = []
    for site_num, group in enumerate(plan.distinct_groups, start=1):
        table_rows.append(
            [
                f"{plan.edge_element}{site_num}",
                group.xyz_label,
                group.cif_label,
                str(len(group.atom_indices)),
            ]
        )

    ax_tbl.text(
        0.0,
        1.0,
        sym_txt,
        transform=ax_tbl.transAxes,
        va="top",
        fontsize=9,
        family="monospace",
    )
    if table_rows:
        table = ax_tbl.table(
            cellText=table_rows,
            colLabels=["Site", "XYZ", "CIF label", "Class size"],
            loc="upper center",
            cellLoc="left",
            bbox=Bbox.from_bounds(0.0, 0.08, 1.0, 0.55),
        )
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1.0, 1.25)
    else:
        ax_tbl.text(
            0.0,
            0.45,
            "No distinct groups found.",
            transform=ax_tbl.transAxes,
            fontsize=10,
        )

    fig.suptitle(title)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
