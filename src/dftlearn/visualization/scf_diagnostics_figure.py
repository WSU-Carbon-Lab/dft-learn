"""SCF convergence diagnostics: 3-column (gnd / exc / tp) history + AUC bar charts."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
from matplotlib.ticker import AutoMinorLocator
from natsort import natsorted


def _shared_ylim_with_margin(
    values: np.ndarray,
    *,
    bottom_zero: bool = False,
    margin_frac: float = 0.05,
) -> tuple[float, float]:
    """Return symmetric-ish y limits with a small margin; handles degenerate ranges."""
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return (0.0, 1.0)
    lo = float(finite.min())
    hi = float(finite.max())
    if hi <= lo:
        span = max(abs(lo), 1.0) * 0.05
        return (lo - span, hi + span)
    span = hi - lo
    m = margin_frac * span
    y0 = 0.0 if bottom_zero else lo - m
    y1 = hi + m
    return (y0, y1)


if TYPE_CHECKING:
    import pandas as pd

_CALC_TITLES = {"gnd": "Ground (gnd)", "exc": "Excited (exc)", "tp": "Transition (tp)"}


def _style_axes_pub(ax: plt.Axes, *, with_grid: bool = True) -> None:
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(
        axis="both",
        which="major",
        direction="in",
        top=True,
        right=True,
        length=4.5,
        width=0.85,
        labelsize=9,
    )
    ax.tick_params(
        axis="both",
        which="minor",
        direction="in",
        top=True,
        right=True,
        length=2.2,
        width=0.6,
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.85)
    if with_grid:
        ax.set_axisbelow(True)
        ax.grid(which="major", linestyle="-.", linewidth=0.45, color="0.82", alpha=1.0)
        ax.grid(which="minor", linestyle=":", linewidth=0.35, color="0.9", alpha=1.0)


def write_scf_diagnostics_figure(
    long_df: pd.DataFrame,
    metrics_df: pd.DataFrame,
    packaged_output_dir: Path,
    dpi: int = 150,
) -> Path:
    """Build a 2x3 figure: SCF history (top) and AUC bars (bottom) for gnd/exc/tp.

    Top row: left y-axis total energy (Ha) in scientific notation, right y-axis
    max density with the same y limits across gnd/exc/tp; one color per site
    (solid energy, dashed density). A vertical dotted line marks the median
    first-``diis_active`` iteration across sites in that column.

    Bottom row: grouped bars per site---energy AUC (left axis) and density AUC
    (twin right) in scientific notation, with shared density-AUC limits across
    columns; site-colored fills with thin black edges. No panel letters or
    figure suptitle.

    Parameters
    ----------
    long_df : pandas.DataFrame
        From :func:`dftlearn.io.stobe_scf_convergence.collect_scf_convergence_long`.
    metrics_df : pandas.DataFrame
        From :func:`dftlearn.io.stobe_scf_convergence.scf_convergence_auc_metrics`.
    packaged_output_dir : pathlib.Path
        Directory for ``scf_diagnostics.png``.
    dpi : int, optional
        Raster resolution.

    Returns
    -------
    pathlib.Path
        Path to the written PNG.
    """
    packaged_output_dir = Path(packaged_output_dir).resolve()
    packaged_output_dir.mkdir(parents=True, exist_ok=True)

    sites = natsorted(long_df["site"].unique())
    calcs = ["gnd", "exc", "tp"]
    cmap = plt.get_cmap("tab10")
    site_colors = {s: cmap(i % 10) for i, s in enumerate(sites)}

    dens_ylim = _shared_ylim_with_margin(
        long_df["max_density"].to_numpy(dtype=np.float64),
        bottom_zero=True,
        margin_frac=0.04,
    )
    dens_auc_ylim = _shared_ylim_with_margin(
        metrics_df["density_auc"].to_numpy(dtype=np.float64),
        bottom_zero=True,
        margin_frac=0.05,
    )

    with plt.rc_context(
        {
            "font.family": "serif",
            "font.size": 10,
            "axes.labelsize": 10,
            "mathtext.fontset": "dejavuserif",
        }
    ):
        fig = plt.figure(figsize=(10.2, 6.35), layout="constrained")
        fig.set_constrained_layout_pads(
            w_pad=0.02,
            h_pad=0.03,
            wspace=0.05,
            hspace=0.06,
        )
        gs = gridspec.GridSpec(
            2,
            3,
            figure=fig,
            height_ratios=[1.25, 1.0],
            hspace=0.09,
            wspace=0.11,
        )

        for j, calc in enumerate(calcs):
            ax_t = fig.add_subplot(gs[0, j])
            ax_r = ax_t.twinx()
            diis_starts: list[float] = []
            for site in sites:
                sub = long_df[
                    (long_df["calc_type"] == calc) & (long_df["site"] == site)
                ]
                if sub.empty:
                    continue
                sub = sub.sort_values("iteration")
                col = site_colors[site]
                ax_t.plot(
                    sub["iteration"],
                    sub["total_energy_h"],
                    color=col,
                    linestyle="-",
                    linewidth=1.15,
                    label=site,
                    zorder=3,
                )
                ax_r.plot(
                    sub["iteration"],
                    sub["max_density"],
                    color=col,
                    linestyle="--",
                    linewidth=1.0,
                    alpha=0.88,
                    zorder=2,
                )
                mone = metrics_df[
                    (metrics_df["calc_type"] == calc) & (metrics_df["site"] == site)
                ]
                if not mone.empty:
                    ds = mone["diis_start_iter"].iloc[0]
                    if np.isfinite(ds):
                        diis_starts.append(float(ds))

            if diis_starts:
                dx = float(np.nanmedian(diis_starts))
                ax_t.axvline(
                    dx,
                    color="0.32",
                    linestyle=":",
                    linewidth=1.35,
                    zorder=1,
                    label="DIIS on",
                )
                ax_r.axvline(dx, color="0.32", linestyle=":", linewidth=1.35, zorder=1)

            ax_t.set_title(_CALC_TITLES[calc], fontsize=10.5, pad=3)
            ax_t.set_xlabel("Iteration")
            if j == 0:
                ax_t.set_ylabel("Energy (Ha)")
            ax_t.ticklabel_format(
                axis="y",
                style="sci",
                scilimits=(0, 0),
                useMathText=True,
            )
            ax_r.set_ylim(dens_ylim)
            if j == 2:
                ax_r.set_ylabel("Max density")
            else:
                ax_r.tick_params(axis="y", which="both", labelright=False)
            _style_axes_pub(ax_t)
            _style_axes_pub(ax_r, with_grid=False)
            ax_r.grid(False)
            if j == 2:
                h1, l1 = ax_t.get_legend_handles_labels()
                by_label = dict(zip(l1, h1, strict=False))
                ax_t.legend(
                    by_label.values(),
                    by_label.keys(),
                    loc="upper right",
                    fontsize=7.5,
                    framealpha=0.95,
                    ncol=1,
                )

            ax_b = fig.add_subplot(gs[1, j])
            mcalc = metrics_df[metrics_df["calc_type"] == calc].set_index("site")
            e_auc = [
                float(mcalc.loc[s, "energy_auc"]) if s in mcalc.index else 0.0
                for s in sites
            ]
            d_auc = [
                float(mcalc.loc[s, "density_auc"]) if s in mcalc.index else 0.0
                for s in sites
            ]
            x = np.arange(len(sites), dtype=np.float64)
            w = 0.34
            cols = [site_colors[s] for s in sites]
            ax_b.bar(
                x - w / 2,
                e_auc,
                w,
                color=cols,
                edgecolor="black",
                linewidth=0.55,
                label="Energy AUC",
                zorder=2,
            )
            ax_b2 = ax_b.twinx()
            ax_b2.bar(
                x + w / 2,
                d_auc,
                w,
                color=cols,
                edgecolor="black",
                linewidth=0.55,
                alpha=0.55,
                label="Density AUC",
                zorder=2,
            )
            ax_b.set_xticks(x, sites, rotation=45, ha="right")
            ax_b.set_xlabel("Site")
            if j == 0:
                ax_b.set_ylabel("Energy AUC")
            ax_b.ticklabel_format(
                axis="y",
                style="sci",
                scilimits=(0, 0),
                useMathText=True,
            )
            ax_b2.set_ylim(dens_auc_ylim)
            if j == 2:
                ax_b2.set_ylabel("Density AUC")
            else:
                ax_b2.tick_params(axis="y", which="both", labelright=False)
            ax_b2.ticklabel_format(
                axis="y",
                style="sci",
                scilimits=(0, 0),
                useMathText=True,
            )
            _style_axes_pub(ax_b)
            ax_b2.yaxis.set_minor_locator(AutoMinorLocator())
            ax_b2.tick_params(
                axis="y",
                which="major",
                direction="in",
                length=4.0,
                labelsize=9,
            )
            ax_b2.tick_params(axis="y", which="minor", direction="in", length=2.0)
            for spine in ax_b2.spines.values():
                spine.set_visible(True)
                spine.set_linewidth(0.85)
        out = packaged_output_dir / "scf_diagnostics.png"
        fig.savefig(out, dpi=dpi, bbox_inches="tight", facecolor="white")
        plt.close(fig)

    return out


def write_scf_diagnostics_bundle(
    run_root: Path,
    packaged_output_dir: Path,
    dpi: int = 150,
) -> tuple[Path, Path, Path] | None:
    """Write long + metrics CSVs and the diagnostics figure under ``packaged_output``.

    Parameters
    ----------
    run_root : pathlib.Path
        StoBe run root with per-site ``*gnd.out`` / ``*exc.out`` / ``*tp.out``.
    packaged_output_dir : pathlib.Path
        Output directory (created if needed).
    dpi : int, optional
        Figure DPI.

    Returns
    -------
    tuple[pathlib.Path, pathlib.Path, pathlib.Path] | None
        ``(scf_convergence_long.csv, scf_convergence_metrics.csv, scf_diagnostics.png)``
        or ``None`` if no convergence tables were found.
    """
    from dftlearn.io.stobe_scf_convergence import (
        collect_scf_convergence_long,
        scf_convergence_auc_metrics,
    )

    packaged_output_dir = Path(packaged_output_dir).resolve()
    packaged_output_dir.mkdir(parents=True, exist_ok=True)
    try:
        long_df = collect_scf_convergence_long(run_root)
    except ValueError:
        return None
    metrics_df = scf_convergence_auc_metrics(long_df)
    long_csv = packaged_output_dir / "scf_convergence_long.csv"
    metrics_csv = packaged_output_dir / "scf_convergence_metrics.csv"
    long_df.to_csv(long_csv, index=False)
    metrics_df.to_csv(metrics_csv, index=False)
    fig_path = write_scf_diagnostics_figure(
        long_df,
        metrics_df,
        packaged_output_dir,
        dpi=dpi,
    )
    return (long_csv, metrics_csv, fig_path)
