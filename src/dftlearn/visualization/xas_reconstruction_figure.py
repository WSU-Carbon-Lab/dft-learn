"""Compose a publication-style Delta-KS aligned XAS summary figure."""

from __future__ import annotations

import itertools
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgb, to_rgba
from matplotlib.layout_engine import ConstrainedLayoutEngine
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator
from natsort import natsorted

from dftlearn.io.xyz_structure import (
    site_label_to_atom_index_from_rows,
    xyz_rows_from_file,
)
from dftlearn.visualization.xyz_wireframe import draw_xyz_wireframe_on_ax
from dftlearn.xas.c3_symmetry import write_c3_frame_json
from dftlearn.xas.spectrum import aligned_dipole_tensor_tables, collect_site_tp_xas

if TYPE_CHECKING:
    import pandas as pd
    from matplotlib.axes import Axes

_TENSOR_STYLES: tuple[tuple[str, str, tuple[float, float, float]], ...] = (
    ("abs_xx", r"$I_{xx}$", (0.0, 0.447, 0.698)),
    ("abs_yy", r"$I_{yy}$", (0.902, 0.624, 0.0)),
    ("abs_zz", r"$I_{zz}$", (0.0, 0.62, 0.451)),
)
_TENSOR_LINESTYLES: tuple[str, ...] = ("-", "--", "-.")


def write_xas_reconstruction_report(
    run_root: Path,
    packaged_output_dir: Path,
    *,
    xray_filename: str = "XrayT001.out",
    xyz_path: Path | None = None,
    c3_symmetrize: bool = False,
    dpi: int = 150,
) -> tuple[Path, Path, Path, Path, Path, Path] | None:
    r"""Write TP XAS tables, Delta-KS aligned dipole tensors, and a summary figure.

    Reconstructs xrayspec-equivalent Gaussians from ``{site}.xas`` sticks and
    overlays them on StoBe ``XrayT*.out``. Delta-KS alignment, when SCF totals
    exist, is a rigid x-shift of already-broadened intensities. Cartesian dipole
    components are exported after the same :math:`E^c` translation. The summary
    figure compares StoBe and TP after that shift, with per-site stick skylines
    under each absorption trace and site-mean tensor components.

    When ``c3_symmetrize`` is True, Cartesian OS and tensor spectra use the
    Al-N/O C3-folded frame and ``c3_frame.json`` is written beside the CSVs.

    Parameters
    ----------
    run_root : pathlib.Path
        StoBe run directory with site subfolders.
    packaged_output_dir : pathlib.Path
        Output folder (created if missing).
    xray_filename : str, optional
        StoBe broadened table name inside each site directory.
    xyz_path : pathlib.Path, optional
        Geometry for the summary structure panel and for C3 frame construction.
        When omitted or unreadable for drawing, that panel is dropped and the
        overview spans the top row.
    c3_symmetrize : bool, optional
        Fold Cartesian dipoles under C3 in the Al-N/O molecular frame.
    dpi : int, optional
        PNG resolution.

    Returns
    -------
    tuple of pathlib.Path, or None
        Paths to the long CSV, metrics CSV, sticks CSV, per-site aligned
        dipole-tensor CSV, site-mean tensor CSV, and summary PNG, or ``None``
        if no sticks.

    Raises
    ------
    ValueError
        If stick files exist but cannot be parsed or aligned.
    """
    run_root = Path(run_root).resolve()
    packaged_output_dir = Path(packaged_output_dir).resolve()
    try:
        _energy, spectra, metrics, sticks, c3_frame = collect_site_tp_xas(
            run_root,
            xray_filename=xray_filename,
            xyz_path=xyz_path,
            c3_symmetrize=c3_symmetrize,
        )
    except FileNotFoundError:
        return None
    packaged_output_dir.mkdir(parents=True, exist_ok=True)
    if c3_frame is not None:
        write_c3_frame_json(c3_frame, packaged_output_dir / "c3_frame.json")
    long_path = packaged_output_dir / "xas_tp_reconstructed_long.csv"
    metrics_path = packaged_output_dir / "xas_tp_reconstruction_metrics.csv"
    sticks_path = packaged_output_dir / "xas_tp_sticks_long.csv"
    tensor_long_path = packaged_output_dir / "xas_dipole_tensor_aligned_long.csv"
    tensor_mean_path = packaged_output_dir / "xas_dipole_tensor_aligned_mean.csv"
    spectra.to_csv(long_path, index=False)
    metrics.to_csv(metrics_path, index=False)
    sticks.to_csv(sticks_path, index=False)
    tensor_long, tensor_mean = aligned_dipole_tensor_tables(spectra, metrics)
    tensor_long.to_csv(tensor_long_path, index=False)
    tensor_mean.to_csv(tensor_mean_path, index=False)
    summary_path = packaged_output_dir / "xas_tp_aligned_summary.png"
    resolved_xyz = Path(xyz_path).resolve() if xyz_path is not None else None
    _write_summary_figure(
        spectra,
        metrics,
        sticks,
        summary_path,
        xyz_path=resolved_xyz,
        dpi=dpi,
    )
    return (
        long_path,
        metrics_path,
        sticks_path,
        tensor_long_path,
        tensor_mean_path,
        summary_path,
    )


def _energy_limits(spectra: pd.DataFrame) -> tuple[float, float]:
    """Return the plotted photon-energy window from the spectral table."""
    energy = spectra["energy_ev"].to_numpy(dtype=np.float64)
    return float(np.nanmin(energy)), float(np.nanmax(energy))


def _draw_transition_skyline(
    ax: Axes,
    sticks: pd.DataFrame,
    site_color: dict[str, tuple[float, float, float, float]],
    energy_limits: tuple[float, float],
    *,
    show_legend: bool = True,
) -> Axes:
    """Draw per-transition oscillator-strength sticks colored by site."""
    e_lo, e_hi = energy_limits
    use_aligned = bool(np.any(np.isfinite(sticks["energy_aligned_ev"].to_numpy())))
    energy_col = "energy_aligned_ev" if use_aligned else "energy_ev"
    sites = natsorted(sticks["site"].astype(str).unique())
    handles: list[Line2D] = []
    osc_max = 0.0
    for site in sites:
        chunk = sticks.loc[sticks["site"] == site]
        energy = chunk[energy_col].to_numpy(dtype=np.float64)
        osc = chunk["oscillator_strength"].to_numpy(dtype=np.float64)
        finite = np.isfinite(energy) & np.isfinite(osc)
        in_win = finite & (energy >= e_lo) & (energy <= e_hi)
        energy = energy[in_win]
        osc = osc[in_win]
        color = site_color[site]
        if energy.size:
            osc_max = max(osc_max, float(np.max(osc)))
            segs = np.stack(
                (
                    np.column_stack((energy, np.zeros_like(energy))),
                    np.column_stack((energy, osc)),
                ),
                axis=1,
            )
            ax.add_collection(
                LineCollection(
                    segs.tolist(),
                    colors=(color,),
                    linewidths=0.65,
                    alpha=0.88,
                    zorder=2,
                )
            )
        handles.append(Line2D([0], [0], color=color, linewidth=1.4, label=site))
    ax.set_xlim(e_lo, e_hi)
    ax.set_ylim(0.0, osc_max * 1.08 if osc_max > 0.0 else 1.0)
    if show_legend and len(sites) > 1:
        ax.legend(
            handles=handles,
            loc="upper right",
            fontsize=8,
            framealpha=0.95,
            edgecolor="0.65",
            fancybox=False,
        )
    return ax


def _style_pub_ax(ax: Axes) -> None:
    """Apply inward ticks, full spines, and light grids."""
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(
        axis="both",
        which="major",
        direction="in",
        top=True,
        right=True,
        length=5,
        width=0.9,
        labelsize=10,
    )
    ax.tick_params(
        axis="both",
        which="minor",
        direction="in",
        top=True,
        right=True,
        length=2.5,
        width=0.65,
        labelsize=10,
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.9)
    ax.set_axisbelow(True)
    ax.grid(which="major", linestyle="-", linewidth=0.55, color="0.82")
    ax.grid(which="minor", linestyle=":", linewidth=0.4, color="0.88")


def _write_summary_figure(
    spectra: pd.DataFrame,
    metrics: pd.DataFrame,
    sticks: pd.DataFrame,
    fig_path: Path,
    *,
    xyz_path: Path | None,
    dpi: int,
) -> None:
    """Draw Delta-KS aligned StoBe vs TP, per-site spectra with skylines, and tensor."""
    sites = natsorted(spectra["site"].astype(str).unique())
    cmap = plt.get_cmap("tab10")
    site_color: dict[str, tuple[float, float, float, float]] = {
        site: to_rgba(cmap(i % 10)) for i, site in enumerate(sites)
    }
    shifts = _site_shifts(metrics)
    e_lo, e_hi = _shifted_energy_limits(spectra, shifts)
    n_sites = max(len(sites), 1)
    n_cols = 2
    n_site_rows = max((n_sites + n_cols - 1) // n_cols, 1)
    xyz_rows = _load_xyz_rows(xyz_path)
    header = ["mol", "ov"] if xyz_rows is not None else ["ov", "ov"]
    mosaic = [header]
    pad_keys: list[str] = []
    for r in range(n_site_rows):
        spec_row: list[str] = []
        sky_row: list[str] = []
        for c in range(n_cols):
            idx = r * n_cols + c
            if idx < n_sites:
                spec_row.append(f"s{idx}")
                sky_row.append(f"k{idx}")
            elif n_sites == 1:
                spec_row.append("s0")
                sky_row.append("k0")
            else:
                spec_key = f"pads{idx}"
                sky_key = f"padk{idx}"
                spec_row.append(spec_key)
                sky_row.append(sky_key)
                pad_keys.extend((spec_key, sky_key))
        mosaic.append(spec_row)
        mosaic.append(sky_row)
    mosaic.append(["ten", "ten"])
    height_ratios = [1.22, *([1.05, 0.38] * n_site_rows), 1.12]
    fig_h = 2.45 + 2.55 * n_site_rows + 2.2
    fig, axd = plt.subplot_mosaic(
        mosaic,
        figsize=(7.2, fig_h),
        layout="constrained",
        height_ratios=height_ratios,
        width_ratios=[1.0, 1.12],
    )
    layout_engine = fig.get_layout_engine()
    if isinstance(layout_engine, ConstrainedLayoutEngine):
        layout_engine.set(h_pad=0.03, w_pad=0.04, hspace=0.05, wspace=0.08)
    for key in pad_keys:
        axd[key].set_visible(False)
    letters = _panel_label_letters()
    if xyz_rows is not None:
        _draw_structure_panel(axd["mol"], xyz_rows, sites, site_color)
        _panel_label(axd["mol"], next(letters))
    _draw_overview_panel(axd["ov"], spectra, shifts, sites, site_color, (e_lo, e_hi))
    _panel_label(axd["ov"], next(letters))
    last_row_start = (n_site_rows - 1) * n_cols
    for i, site in enumerate(sites):
        ax_spec = axd[f"s{i}"]
        ax_sky = axd[f"k{i}"]
        ax_sky.sharex(ax_spec)
        _draw_site_panel(
            ax_spec,
            spectra,
            site,
            site_color[site],
            shift_ev=shifts.get(site, float("nan")),
            show_ylabel=i % n_cols == 0,
            show_legend=i == 0,
        )
        _panel_label(ax_spec, next(letters))
        _style_pub_ax(ax_sky)
        _draw_transition_skyline(
            ax_sky,
            sticks.loc[sticks["site"] == site],
            site_color,
            (e_lo, e_hi),
            show_legend=False,
        )
        if i % n_cols == 0:
            ax_sky.set_ylabel("Osc. strength")
        ax_spec.set_xlim(e_lo, e_hi)
        ax_spec.tick_params(axis="x", labelbottom=False)
        if i >= last_row_start:
            ax_sky.set_xlabel("Photon energy (eV)")
        else:
            ax_sky.tick_params(axis="x", labelbottom=False)
    _draw_tensor_panel(axd["ten"], spectra, shifts, (e_lo, e_hi))
    _panel_label(axd["ten"], next(letters))
    axd["ov"].set_xlim(e_lo, e_hi)
    axd["ten"].set_xlim(e_lo, e_hi)
    axd["ten"].set_xlabel("Photon energy (eV)")
    fig.savefig(fig_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _load_xyz_rows(
    xyz_path: Path | None,
) -> list[tuple[str, float, float, float]] | None:
    if xyz_path is None or not xyz_path.is_file():
        return None
    try:
        return xyz_rows_from_file(xyz_path)
    except (OSError, ValueError):
        return None


def _panel_label_letters() -> itertools.chain[str]:
    """Yield panel labels ``a``..``z``, then ``a1``, ``b1``, for many sites."""
    first = (chr(ord("a") + i) for i in range(26))

    def _suffix() -> itertools.chain[str]:
        for n in itertools.count(1):
            for i in range(26):
                yield f"{chr(ord('a') + i)}{n}"

    return itertools.chain(first, _suffix())


def _panel_label(ax: Axes, letter: str) -> None:
    if not letter:
        return
    ax.annotate(
        f"({letter})",
        xy=(0.0, 1.0),
        xycoords="axes fraction",
        xytext=(0, 8),
        textcoords="offset points",
        ha="left",
        va="bottom",
        fontweight="normal",
        fontsize=11,
        annotation_clip=False,
    )


def _site_shifts(metrics: pd.DataFrame) -> dict[str, float]:
    if metrics.empty or "E_c_ev" not in metrics.columns:
        return {}
    sites = metrics["site"].astype(str).to_numpy()
    values = metrics["E_c_ev"].to_numpy(dtype=np.float64)
    return {str(site): float(shift) for site, shift in zip(sites, values, strict=True)}


def _shift_or_zero(shift_ev: float) -> float:
    if np.isfinite(shift_ev):
        return float(shift_ev)
    return 0.0


def _shifted_energy_limits(
    spectra: pd.DataFrame,
    shifts: dict[str, float],
) -> tuple[float, float]:
    lo = float("inf")
    hi = float("-inf")
    for site, chunk in spectra.groupby("site", sort=False):
        energy = chunk["energy_ev"].to_numpy(dtype=np.float64)
        shift = _shift_or_zero(shifts.get(str(site), float("nan")))
        lo = min(lo, float(np.nanmin(energy)) + shift)
        hi = max(hi, float(np.nanmax(energy)) + shift)
    if not np.isfinite(lo):
        return _energy_limits(spectra)
    return lo, hi


def _xshifted_trace(
    spectra: pd.DataFrame,
    site: str,
    column: str,
    shift_ev: float,
) -> tuple[np.ndarray, np.ndarray]:
    chunk = spectra.loc[spectra["site"] == site]
    energy = chunk["energy_ev"].to_numpy(dtype=np.float64)
    values = chunk[column].to_numpy(dtype=np.float64)
    return energy + _shift_or_zero(shift_ev), values


def _mean_xshifted(
    spectra: pd.DataFrame,
    shifts: dict[str, float],
    column: str,
    grid: np.ndarray,
) -> np.ndarray:
    stacked: list[np.ndarray] = []
    sites = natsorted(spectra["site"].astype(str).unique())
    for site in sites:
        x_s, y_s = _xshifted_trace(
            spectra, site, column, shifts.get(site, float("nan"))
        )
        finite = np.isfinite(x_s) & np.isfinite(y_s)
        if not np.any(finite):
            continue
        stacked.append(
            np.interp(grid, x_s[finite], y_s[finite], left=np.nan, right=np.nan)
        )
    if not stacked:
        return np.full(grid.shape, np.nan, dtype=np.float64)
    return np.nanmean(np.vstack(stacked), axis=0)


def _draw_structure_panel(
    ax: Axes,
    rows: list[tuple[str, float, float, float]],
    sites: list[str],
    site_color: dict[str, tuple[float, float, float, float]],
) -> None:
    site_atom_colors: dict[int, tuple[float, float, float]] = {}
    for site in sites:
        try:
            idx = site_label_to_atom_index_from_rows(site, rows)
        except ValueError:
            continue
        rgb = site_color[site][:3]
        site_atom_colors[idx] = (float(rgb[0]), float(rgb[1]), float(rgb[2]))
    draw_xyz_wireframe_on_ax(
        ax,
        rows,
        site_atom_colors,
        show_hydrogen=False,
        bond_lw=1.25,
        label_fontsize=9,
        plot_margins=0.04,
        halo_radius_angstrom=0.40,
        halo_soft_edge=True,
    )
    ax.set_title("Core-excitation sites")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
        spine.set_color("0.65")
    ax.set_aspect("equal")


def _draw_overview_panel(
    ax: Axes,
    spectra: pd.DataFrame,
    shifts: dict[str, float],
    sites: list[str],
    site_color: dict[str, tuple[float, float, float, float]],
    energy_limits: tuple[float, float],
) -> None:
    e_lo, e_hi = energy_limits
    grid = np.linspace(e_lo, e_hi, 2000, dtype=np.float64)
    y_stobe = _mean_xshifted(spectra, shifts, "abs_stobe", grid)
    y_tp = _mean_xshifted(spectra, shifts, "abs_tp", grid)
    _style_pub_ax(ax)
    for site in sites:
        x_s, y_s = _xshifted_trace(
            spectra, site, "abs_tp", shifts.get(site, float("nan"))
        )
        if np.any(np.isfinite(y_s)):
            ax.plot(
                x_s, y_s, color=site_color[site], linewidth=0.85, alpha=0.4, zorder=2
            )
    if np.any(np.isfinite(y_stobe)):
        ax.plot(
            grid,
            y_stobe,
            color="0.15",
            linewidth=1.85,
            label=r"StoBe + $\Delta$KS",
            zorder=4,
        )
    if np.any(np.isfinite(y_tp)):
        ax.plot(
            grid,
            y_tp,
            color=to_rgb("#0072B2"),
            linewidth=1.25,
            linestyle="--",
            label=r"TP + $\Delta$KS",
            zorder=5,
        )
    ax.set_ylabel("Abs. (arb. units)")
    ax.set_title(r"StoBe + $\Delta$KS vs TP + $\Delta$KS")
    ax.set_ylim(bottom=0.0)
    ax.legend(
        loc="upper right",
        fontsize=8,
        framealpha=0.95,
        edgecolor="0.65",
        fancybox=False,
        handlelength=1.8,
    )
    finite_shifts = [
        (site, shifts[site])
        for site in sites
        if np.isfinite(shifts.get(site, float("nan")))
    ]
    if finite_shifts:
        lines = [rf"{site}: $E^c$ = {shift:.2f} eV" for site, shift in finite_shifts]
        ax.text(
            0.03,
            0.97,
            "\n".join(lines),
            transform=ax.transAxes,
            fontsize=8.5,
            va="top",
            ha="left",
            linespacing=1.4,
            zorder=8,
            color="0.1",
            bbox={
                "facecolor": "white",
                "edgecolor": "0.45",
                "linewidth": 0.8,
                "boxstyle": "square,pad=0.4",
                "alpha": 0.96,
            },
        )


def _draw_site_panel(
    ax: Axes,
    spectra: pd.DataFrame,
    site: str,
    color: tuple[float, float, float, float],
    *,
    shift_ev: float,
    show_ylabel: bool,
    show_legend: bool,
) -> None:
    x_stobe, y_stobe = _xshifted_trace(spectra, site, "abs_stobe", shift_ev)
    x_tp, y_tp = _xshifted_trace(spectra, site, "abs_tp", shift_ev)
    _style_pub_ax(ax)
    if np.any(np.isfinite(y_stobe)):
        ax.plot(
            x_stobe,
            y_stobe,
            color="0.35",
            linewidth=1.55,
            label=r"StoBe + $\Delta$KS",
            zorder=3,
        )
    if np.any(np.isfinite(y_tp)):
        ax.plot(
            x_tp,
            y_tp,
            color=color,
            linewidth=1.15,
            linestyle="--",
            label=r"TP + $\Delta$KS",
            zorder=4,
        )
    if np.isfinite(shift_ev):
        ax.set_title(rf"{site}  ($E^c = {shift_ev:.2f}$ eV)", loc="right", fontsize=9)
    else:
        ax.set_title(site, loc="right")
    ax.set_ylim(bottom=0.0)
    if show_ylabel:
        ax.set_ylabel("Abs. (arb. units)")
    if show_legend:
        ax.legend(
            loc="upper right",
            fontsize=7,
            framealpha=0.92,
            edgecolor="0.65",
            fancybox=False,
            handlelength=1.6,
        )


def _draw_tensor_panel(
    ax: Axes,
    spectra: pd.DataFrame,
    shifts: dict[str, float],
    energy_limits: tuple[float, float],
) -> None:
    e_lo, e_hi = energy_limits
    grid = np.linspace(e_lo, e_hi, 2000, dtype=np.float64)
    _style_pub_ax(ax)
    plotted = False
    components: list[np.ndarray] = []
    for (column, label, rgb), ls in zip(
        _TENSOR_STYLES, _TENSOR_LINESTYLES, strict=True
    ):
        values = _mean_xshifted(spectra, shifts, column, grid)
        if not np.any(np.isfinite(values)):
            continue
        ax.plot(
            grid,
            values,
            color=rgb,
            linestyle=ls,
            linewidth=1.45,
            label=label,
            zorder=3,
        )
        components.append(values)
        plotted = True
    if plotted and len(components) == 3:
        iso = (components[0] + components[1] + components[2]) / 3.0
        ax.plot(
            grid,
            iso,
            color="0.15",
            linestyle=":",
            linewidth=1.2,
            label=r"$I_\mathrm{iso}$",
            zorder=4,
        )
    if not plotted:
        y_tp = _mean_xshifted(spectra, shifts, "abs_tp", grid)
        ax.plot(grid, y_tp, color="0.15", linewidth=1.4, label=r"TP + $\Delta$KS")
    ax.set_ylabel("Abs. (arb. units)")
    ax.set_title(r"Dipole tensor (site-mean $I_{ii}$)", loc="right")
    ax.set_ylim(bottom=0.0)
    ax.legend(
        loc="upper right",
        fontsize=8,
        framealpha=0.95,
        edgecolor="0.65",
        fancybox=False,
        ncol=2,
        handlelength=1.8,
    )
