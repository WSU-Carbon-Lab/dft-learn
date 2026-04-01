"""Multi-site XAS summary plots with XYZ wireframe structure (trace-matched colors)."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb
from matplotlib.ticker import AutoMinorLocator
from natsort import natsorted

from dftlearn.io.xray_out import (
    collect_site_xray_spectra,
    site_spectra_to_long_frame,
)
from dftlearn.io.xyz_structure import (
    site_label_to_atom_index_from_rows,
    xyz_rows_from_file,
)
from dftlearn.visualization.xyz_wireframe import draw_xyz_wireframe_on_ax


def write_xas_site_report(
    run_root: Path,
    packaged_output_dir: Path,
    xyz_path: Path,
    xray_filename: str = "XrayT001.out",
    dpi: int = 150,
) -> tuple[Path, Path]:
    """Write long-form CSV and a spectral summary figure.

    Spectra use the per-site mean plus individual traces with publication-style
    axes (full spine box, inward major and minor ticks, major and minor grids).
    The molecular XY wireframe sits in a large top-right inset (skeletal carbons,
    heteroatom labels, trace-colored site halos). Site indices match ``rows``
    order (``C1`` = first carbon in the XYZ file).

    Parameters
    ----------
    run_root : pathlib.Path
        StoBe directory containing ``C1/``, ``C2/``, ... site folders.
    packaged_output_dir : pathlib.Path
        Output folder (created if missing); receives ``xray_spectra_long.csv`` and
        ``xas_site_summary.png``.
    xyz_path : pathlib.Path
        Geometry file (with or without XYZ header; StoBe-style labels such as
        ``C01`` supported).
    xray_filename : str, optional
        Spectrum file name inside each site directory.
    dpi : int, optional
        Raster resolution for the PNG figure.

    Returns
    -------
    csv_path : pathlib.Path
        Path to the written long-form table.
    figure_path : pathlib.Path
        Path to the PNG summary figure.

    Raises
    ------
    FileNotFoundError
        If required inputs are missing.
    ValueError
        If spectra cannot be aligned or the XYZ cannot be parsed.
    """
    run_root = Path(run_root).resolve()
    packaged_output_dir = Path(packaged_output_dir).resolve()
    xyz_path = Path(xyz_path).resolve()
    packaged_output_dir.mkdir(parents=True, exist_ok=True)

    energy, spectra = collect_site_xray_spectra(run_root, xray_filename=xray_filename)
    long_df = site_spectra_to_long_frame(energy, spectra)
    csv_path = packaged_output_dir / "xray_spectra_long.csv"
    long_df.to_csv(csv_path, index=False)

    sites = natsorted(spectra.keys())
    mat = np.stack([spectra[s] for s in sites], axis=0)
    mean_curve = np.mean(mat, axis=0)

    rows = xyz_rows_from_file(xyz_path)
    cmap = plt.get_cmap("tab10")
    site_atom_colors: dict[int, tuple[float, float, float]] = {}
    for i, site in enumerate(sites):
        atom_idx = site_label_to_atom_index_from_rows(site, rows)
        site_atom_colors[atom_idx] = to_rgb(cmap(i % 10))

    fig, ax_spec = plt.subplots(figsize=(9.2, 5.4), layout="constrained")

    _pub_fs = 11
    _tick_fs = 10
    ax_spec.plot(
        energy,
        mean_curve,
        color="black",
        linewidth=1.85,
        label="Mean",
        zorder=len(sites) + 2,
    )
    for i, site in enumerate(sites):
        c = cmap(i % 10)
        ax_spec.plot(
            energy,
            spectra[site],
            color=c,
            linewidth=1.05,
            alpha=0.78,
            label=site,
            zorder=i + 1,
        )
    first_non_zero = int(np.where(mean_curve > 0)[0][0])
    lo = max(0, first_non_zero - 10)
    ax_spec.set_xlim(energy[lo], energy.max())
    ax_spec.set_ylim(0, None)
    ax_spec.set_xlabel("Photon energy (eV)", fontsize=_pub_fs)
    ax_spec.set_ylabel("Abs. (arb. units)", fontsize=_pub_fs)
    ax_spec.set_title("Site-resolved X-ray spectra", fontsize=12, pad=8)

    ax_spec.xaxis.set_minor_locator(AutoMinorLocator())
    ax_spec.yaxis.set_minor_locator(AutoMinorLocator())
    ax_spec.tick_params(
        axis="both",
        which="major",
        direction="in",
        top=True,
        right=True,
        length=5,
        width=0.9,
        labelsize=_tick_fs,
    )
    ax_spec.tick_params(
        axis="both",
        which="minor",
        direction="in",
        top=True,
        right=True,
        length=2.5,
        width=0.65,
        labelsize=_tick_fs,
    )
    for spine in ax_spec.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.9)

    ax_spec.set_axisbelow(True)
    ax_spec.grid(
        which="major",
        linestyle="-",
        linewidth=0.55,
        color="0.82",
        alpha=1.0,
        zorder=0,
    )
    ax_spec.grid(
        which="minor",
        linestyle=":",
        linewidth=0.4,
        color="0.88",
        alpha=1.0,
        zorder=0,
    )

    ax_spec.legend(
        loc="lower right",
        fontsize=9,
        framealpha=0.95,
        edgecolor="0.65",
        fancybox=False,
        ncols=len(sites) + 1,
        handlelength=1.35,
        handletextpad=0.5,
        columnspacing=0.9,
    )

    ax_inset = ax_spec.inset_axes(
        (0.55, 0.26, 0.505, 0.715),
        transform=ax_spec.transAxes,
        facecolor="white",
        zorder=10,
    )
    ax_inset.patch.set_edgecolor("0.65")
    ax_inset.patch.set_linewidth(0.9)
    draw_xyz_wireframe_on_ax(
        ax_inset,
        rows,
        site_atom_colors,
        show_hydrogen=False,
        bond_lw=1.45,
        label_fontsize=11,
        plot_margins=0.018,
        halo_radius_angstrom=0.36,
        halo_soft_edge=False,
    )
    ax_inset.axis("off")
    for spine in ax_inset.spines.values():
        spine.set_visible(False)

    figure_path = packaged_output_dir / "xas_site_summary.png"
    fig.savefig(figure_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    return csv_path, figure_path


def load_xray_csv_summary(csv_path: Path) -> pd.DataFrame:
    """Load a long-form X-ray table written by :func:`write_xas_site_report`.

    Parameters
    ----------
    csv_path : pathlib.Path
        Path to ``xray_spectra_long.csv``.

    Returns
    -------
    pandas.DataFrame
        Columns ``energy_ev``, ``abs``, ``site``. Legacy files may use ``value``
        instead of ``abs``; normalize by renaming when present.
    """
    df = pd.read_csv(csv_path)
    if "abs" not in df.columns and "value" in df.columns:
        df = df.rename(columns={"value": "abs"})
    return df
