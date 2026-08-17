"""StoBe diagnostics: per-site four-row figures and optional multi-site SI layout."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator
from natsort import natsorted

from dftlearn.io.stobe_final_energy import collect_delta_ks_site_table
from dftlearn.io.stobe_orbital_table import (
    analyze_exc_orbitals,
    analyze_gnd_orbitals,
    analyze_tp_orbitals,
    element_symbol_from_site,
    format_final_energy_text,
    parse_stobe_orbital_energies_table,
    reference_k_shell_binding_ev,
)

if TYPE_CHECKING:
    import pandas as pd

    from dftlearn.io.stobe_orbital_table import (
        ExcOrbitalReport,
        GndOrbitalReport,
        TpOrbitalReport,
    )


def _style_axes_pub(
    ax: plt.Axes, *, with_grid: bool = True, compact: bool = False
) -> None:
    ls = 6.25 if compact else 8.0
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(
        axis="both",
        which="major",
        direction="in",
        top=True,
        right=True,
        length=3.5 if compact else 4.0,
        width=0.75 if compact else 0.8,
        labelsize=ls,
    )
    ax.tick_params(
        axis="both",
        which="minor",
        direction="in",
        top=True,
        right=True,
        length=2.0,
        width=0.55,
        labelsize=ls,
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
    if with_grid:
        ax.set_axisbelow(True)
        ax.grid(which="major", linestyle="-.", linewidth=0.4, color="0.82", alpha=1.0)
        ax.grid(which="minor", linestyle=":", linewidth=0.35, color="0.9", alpha=1.0)


def _occ_color(occ: float) -> str:
    if occ >= 0.999:
        return "#1f77b4"
    if occ <= 0.001:
        return "#2ca02c"
    if abs(occ - 0.5) < 0.03:
        return "#ff7f0e"
    return "#d62728"


def _plot_alpha_levels(
    ax: plt.Axes,
    df: pd.DataFrame,
    title: str,
    *,
    compact: bool = False,
) -> None:
    """Draw horizontal segments for alpha levels in a valence-focused window."""
    emin = float(df["alpha_ev"].min())
    emax = float(df["alpha_ev"].max())
    lo = max(emin - 2.0, -35.0)
    hi = min(emax + 2.0, 8.0)
    sub = df[(df["alpha_ev"] >= lo) & (df["alpha_ev"] <= hi)].copy()
    if sub.empty:
        sub = df.tail(min(40, len(df)))
    lw = 1.65 if compact else 2.0
    for _, r in sub.iterrows():
        lv = int(r["level"])
        e = float(r["alpha_ev"])
        occ = float(r["alpha_occ"])
        c = _occ_color(occ)
        ax.hlines(e, lv - 0.45, lv + 0.45, colors=c, linewidth=lw, zorder=2)
    fs = 6.5 if compact else 8.0
    ax.set_xlabel("Level index", fontsize=fs)
    ax.set_ylabel(r"$\alpha$ energy (eV)", fontsize=fs)
    ax.set_title(title, fontsize=fs + 0.5)
    _style_axes_pub(ax, compact=compact)


def _save_stobe_diagnostic_png(path: Path, fig: plt.Figure, *, dpi: int) -> None:
    """Write a PNG with a white figure background."""
    fig.savefig(
        path,
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
        edgecolor="none",
        format="png",
    )


def _add_occupation_legend(fig: plt.Figure, *, ncol: int = 4) -> None:
    """Shared legend for alpha occupation coloring (SI-style, figure-level)."""
    handles = [
        Line2D([0], [0], color="#1f77b4", lw=4, solid_capstyle="butt"),
        Line2D([0], [0], color="#2ca02c", lw=4, solid_capstyle="butt"),
        Line2D([0], [0], color="#ff7f0e", lw=4, solid_capstyle="butt"),
        Line2D([0], [0], color="#d62728", lw=4, solid_capstyle="butt"),
    ]
    labels = [
        r"$\alpha$ occ. $\approx$ 1",
        r"$\alpha$ occ. $\approx$ 0",
        r"$\alpha$ occ. $\approx$ 0.5 (TP core)",
        "other partial",
    ]
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=ncol,
        frameon=False,
        fontsize=7.5,
        columnspacing=1.4,
        handletextpad=0.6,
        bbox_to_anchor=(0.5, 0.005),
    )


def _txt_gnd(rep: GndOrbitalReport, bind: float) -> str:
    return (
        f"K-edge ref. binding = {bind:.1f} eV (target -E ~ {-bind:.1f} eV)\n"
        f"Core state: level {rep.core_level},  E = {rep.core_energy_ev:.4f} eV\n"
        f"HOMO: level {rep.homo_level},  E = {rep.homo_energy_ev:.4f} eV\n"
        f"LUMO: level {rep.lumo_level},  E = {rep.lumo_energy_ev:.4f} eV"
    )


def _txt_exc(rep: ExcOrbitalReport, gnd: GndOrbitalReport) -> str:
    core_e = rep.core_energy_ev
    lines = [
        f"Core hole (alpha occ ~0): level {rep.core_level},  E = {core_e:.4f} eV",
    ]
    if rep.homo_level == gnd.lumo_level:
        lines.append(
            f"Excited HOMO (former ground LUMO): level {rep.homo_level},  "
            f"E = {rep.homo_energy_ev:.4f} eV"
        )
    else:
        lines.append(
            f"HOMO: level {rep.homo_level},  E = {rep.homo_energy_ev:.4f} eV"
        )
    lines.extend(
        [
            f"LUMO: level {rep.lumo_level},  E = {rep.lumo_energy_ev:.4f} eV",
            f"Ground ref. HOMO: level {gnd.homo_level} at {gnd.homo_energy_ev:.4f} eV",
            f"Ground ref. LUMO: level {gnd.lumo_level} at {gnd.lumo_energy_ev:.4f} eV",
        ]
    )
    if (
        rep.homo_plus_one_level is not None
        and rep.homo_plus_one_level != rep.homo_level
        and rep.homo_plus_one_energy_ev is not None
    ):
        lines.append(
            f"Prior LUMO now filled: level {rep.homo_plus_one_level},  "
            f"E = {rep.homo_plus_one_energy_ev:.4f} eV"
        )
    return "\n".join(lines)


def _txt_tp(rep: TpOrbitalReport) -> str:
    return (
        f"Core (alpha occ ~{rep.core_alpha_occ:.2f}): level {rep.core_level},  "
        f"E = {rep.core_energy_ev:.4f} eV\n"
        f"HOMO: level {rep.homo_level},  E = {rep.homo_energy_ev:.4f} eV\n"
        f"LUMO: level {rep.lumo_level},  E = {rep.lumo_energy_ev:.4f} eV"
    )


def _populate_site_diagnostic(
    fig: plt.Figure,
    gs: gridspec.GridSpec,
    site: str,
    path_gnd: Path,
    path_exc: Path,
    path_tp: Path,
    row_summary: pd.Series,
    *,
    compact: bool = False,
    panel_label: str | None = None,
) -> None:
    """Draw one site's four-row diagnostic into ``gs`` (4 rows, 2 cols; row 3 spans)."""
    sym = element_symbol_from_site(site)
    fe_gnd = format_final_energy_text(path_gnd)
    fe_exc = format_final_energy_text(path_exc)
    fe_tp = format_final_energy_text(path_tp)

    og = parse_stobe_orbital_energies_table(path_gnd)
    oe = parse_stobe_orbital_energies_table(path_exc)
    ot = parse_stobe_orbital_energies_table(path_tp)

    gnd_rep: GndOrbitalReport | None = None
    exc_rep: ExcOrbitalReport | None = None
    tp_rep: TpOrbitalReport | None = None
    bind_ev = float("nan")
    fail_note = ""
    try:
        bind_ev = reference_k_shell_binding_ev(sym)
        gnd_rep = analyze_gnd_orbitals(og, bind_ev)
        exc_rep = analyze_exc_orbitals(oe, gnd_rep)
        tp_rep = analyze_tp_orbitals(ot)
    except Exception as exc:
        fail_note = f"\n\n[Orbital analysis failed: {exc}]"

    if gnd_rep is not None and np.isfinite(bind_ev):
        m_gnd = _txt_gnd(gnd_rep, bind_ev) + fail_note
    else:
        m_gnd = f"(K-edge / core assignment unavailable){fail_note}"
    if exc_rep is not None and gnd_rep is not None:
        m_exc = _txt_exc(exc_rep, gnd_rep) + fail_note
    else:
        m_exc = "(Excited metrics unavailable)" + fail_note
    if tp_rep is not None:
        m_tp = _txt_tp(tp_rep) + fail_note
    else:
        m_tp = "(TP metrics unavailable)" + fail_note

    mono = 6.2 if compact else 7.5
    ptitle = "Valence alpha levels" if compact else "Alpha levels (valence window)"
    rows = [
        (0, f"Ground ({site}gnd)", fe_gnd, m_gnd, og),
        (1, f"Excited ({site}exc)", fe_exc, m_exc, oe),
        (2, f"Transition potential ({site}tp)", fe_tp, m_tp, ot),
    ]
    for r, title, fe_txt, metrics_txt, odf in rows:
        ax_l = fig.add_subplot(gs[r, 0])
        ax_l.axis("off")
        block = f"{title}\n\n{fe_txt}\n\n--- Orbital metrics ---\n{metrics_txt}"
        ax_l.text(
            0.02,
            0.98,
            block,
            transform=ax_l.transAxes,
            va="top",
            ha="left",
            fontsize=mono,
            family="monospace",
            bbox={
                "boxstyle": "round,pad=0.35",
                "facecolor": "1.0",
                "edgecolor": "0.82",
                "linewidth": 0.45,
            },
        )
        if panel_label is not None and r == 0:
            ax_l.text(
                0.02,
                1.0,
                panel_label,
                transform=ax_l.transAxes,
                va="bottom",
                ha="left",
                fontsize=11 if compact else 12,
                fontweight="bold",
            )
        ax_r = fig.add_subplot(gs[r, 1])
        _plot_alpha_levels(ax_r, odf, ptitle, compact=compact)

    ax_sum = fig.add_subplot(gs[3, :])
    corr = float(row_summary["E_c_deltaKS_ev"])
    dexc = float(row_summary["delta_exc_vs_gnd_ev"])
    dtp = float(row_summary["delta_tp_vs_gnd_ev"])
    xs = np.arange(3)
    labs = [
        r"$\Delta$-KS correction",
        "Excited $-$ ground",
        "TP $-$ ground",
    ]
    vals = [corr, dexc, dtp]
    cols = ["#8c564b", "#1f77b4", "#2ca02c"]
    ax_sum.bar(xs, vals, color=cols, edgecolor="0.25", linewidth=0.45, zorder=3)
    ax_sum.set_xticks(xs)
    rot = 15 if compact else 12
    xfs = 6.5 if compact else 8
    ax_sum.set_xticklabels(labs, rotation=rot, ha="right", fontsize=xfs)
    ax_sum.set_ylabel("Energy (eV)", fontsize=6.8 if compact else 9)
    sum_fs = 7.0 if compact else 9.5
    ax_sum.set_title(f"Summary ({site}): shifts and correction", fontsize=sum_fs)
    ax_sum.axhline(0.0, color="0.45", linestyle="--", linewidth=0.75, zorder=1)
    _style_axes_pub(ax_sum, compact=compact)


def write_stobe_site_diagnostic_figure(
    site: str,
    path_gnd: Path,
    path_exc: Path,
    path_tp: Path,
    row_summary: pd.Series,
    packaged_output_dir: Path,
    dpi: int = 150,
) -> Path:
    """Build a 4-row figure: gnd, exc, tp orbital diagnostics, then summary bars."""
    packaged_output_dir = Path(packaged_output_dir).resolve()
    packaged_output_dir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(11.2, 12.5), layout="constrained")
    gs = gridspec.GridSpec(
        4, 2, figure=fig, height_ratios=[1.0, 1.0, 1.0, 0.72], wspace=0.2
    )
    _populate_site_diagnostic(
        fig,
        gs,
        site,
        path_gnd,
        path_exc,
        path_tp,
        row_summary,
        compact=False,
        panel_label=None,
    )

    out = packaged_output_dir / f"stobe_diagnostics_{site}.png"
    _save_stobe_diagnostic_png(out, fig, dpi=dpi)
    plt.close(fig)
    return out


def write_stobe_combined_diagnostic_figure(
    run_root: Path,
    packaged_output_dir: Path,
    wide: pd.DataFrame,
    dpi: int = 300,
) -> Path | None:
    """Write one multi-panel SI figure (all sites) as a white-background PNG.

    Arranges each site in a 2-by-2 grid with panel labels (a)--(d). Adds a single
    occupation legend. Saves ``stobe_diagnostics_all_sites.png`` under
    ``packaged_output_dir``.

    Parameters
    ----------
    run_root : pathlib.Path
        StoBe run root (site subfolders).
    packaged_output_dir : pathlib.Path
        Output directory (e.g. ``.../packaged_output``).
    wide : pandas.DataFrame
        Wide Delta-KS table from
        :func:`~dftlearn.io.stobe_final_energy.collect_delta_ks_site_table`.
    dpi : int, optional
        Raster resolution for the PNG.

    Returns
    -------
    pathlib.Path | None
        Path to the PNG if at least two sites with gnd/exc/tp outputs exist;
        otherwise ``None``.
    """
    from dftlearn.io.stobe_scf_convergence import discover_site_stobe_out

    packaged_output_dir = Path(packaged_output_dir).resolve()
    run_root = Path(run_root).resolve()
    gnd = dict(discover_site_stobe_out(run_root, "gnd"))
    exc = dict(discover_site_stobe_out(run_root, "exc"))
    tp = dict(discover_site_stobe_out(run_root, "tp"))

    sites: list[str] = []
    for s in natsorted(wide["site"].astype(str).unique()):
        if s in gnd and s in exc and s in tp:
            sites.append(str(s))
    if len(sites) < 2:
        return None

    n = len(sites)
    ncols = 2 if n > 1 else 1
    nrows = int(np.ceil(n / ncols))

    fig_w = 7.4 if n <= 4 else min(7.4 + 0.15 * (nrows - 2), 8.2)
    fig_h = 3.55 * nrows + 0.55
    fig = plt.figure(figsize=(fig_w, fig_h))
    outer = fig.add_gridspec(
        nrows,
        ncols,
        hspace=0.38,
        wspace=0.26,
        left=0.07,
        right=0.98,
        top=0.97,
        bottom=0.11,
    )

    def _panel_letter(idx: int) -> str:
        if idx < 26:
            return f"({chr(ord('a') + idx)})"
        return f"({idx + 1})"

    for k, site in enumerate(sites):
        i, j = divmod(k, ncols)
        inner = outer[i, j].subgridspec(
            4,
            2,
            height_ratios=[1.0, 1.0, 1.0, 0.62],
            hspace=0.28,
            wspace=0.14,
        )
        row = wide.loc[wide["site"] == site].iloc[0]
        _populate_site_diagnostic(
            fig,
            inner,
            site,
            gnd[site],
            exc[site],
            tp[site],
            row,
            compact=True,
            panel_label=_panel_letter(k),
        )

    _add_occupation_legend(fig, ncol=4)
    png_path = packaged_output_dir / "stobe_diagnostics_all_sites.png"
    _save_stobe_diagnostic_png(png_path, fig, dpi=dpi)
    plt.close(fig)
    return png_path


def write_stobe_final_energy_figure(
    df: pd.DataFrame,
    packaged_output_dir: Path,
    dpi: int = 150,
) -> Path:
    """Redirect to the first-site diagnostic figure (deprecated).

    Kept for API compatibility. Prefer :func:`write_stobe_diagnostic_figures_all_sites`.
    """
    packaged_output_dir = Path(packaged_output_dir).resolve()
    if df.empty:
        msg = "Empty dataframe for diagnostic figure"
        raise ValueError(msg)
    site = str(df.iloc[0]["site"])
    row = df.loc[df["site"] == site].iloc[0]
    run_root = packaged_output_dir.parent
    from dftlearn.io.stobe_scf_convergence import discover_site_stobe_out

    gnd = dict(discover_site_stobe_out(run_root, "gnd"))
    exc = dict(discover_site_stobe_out(run_root, "exc"))
    tp = dict(discover_site_stobe_out(run_root, "tp"))
    if site not in gnd or site not in exc or site not in tp:
        msg = f"Missing gnd/exc/tp outputs for site {site}"
        raise ValueError(msg)
    return write_stobe_site_diagnostic_figure(
        site,
        gnd[site],
        exc[site],
        tp[site],
        row,
        packaged_output_dir,
        dpi=dpi,
    )


def write_stobe_diagnostic_figures_all_sites(
    run_root: Path,
    packaged_output_dir: Path,
    wide: pd.DataFrame,
    dpi: int = 150,
) -> list[Path]:
    """Write one ``stobe_diagnostics_<site>.png`` per site row in ``wide``."""
    paths: list[Path] = []
    from dftlearn.io.stobe_scf_convergence import discover_site_stobe_out

    gnd = dict(discover_site_stobe_out(run_root, "gnd"))
    exc = dict(discover_site_stobe_out(run_root, "exc"))
    tp = dict(discover_site_stobe_out(run_root, "tp"))

    for site in natsorted(wide["site"].astype(str).unique()):
        if site not in gnd or site not in exc or site not in tp:
            continue
        row = wide.loc[wide["site"] == site].iloc[0]
        paths.append(
            write_stobe_site_diagnostic_figure(
                site,
                gnd[site],
                exc[site],
                tp[site],
                row,
                packaged_output_dir,
                dpi=dpi,
            )
        )
    return paths


def write_stobe_final_energy_bundle(
    run_root: Path,
    packaged_output_dir: Path,
    dpi: int = 150,
) -> tuple[Path, list[Path]] | None:
    """Write CSV and one orbital/SCF energy summary PNG (all sites)."""
    packaged_output_dir = Path(packaged_output_dir).resolve()
    packaged_output_dir.mkdir(parents=True, exist_ok=True)
    try:
        wide = collect_delta_ks_site_table(run_root)
    except ValueError:
        return None
    if wide.empty:
        return None
    csv_path = packaged_output_dir / "stobe_final_energies.csv"
    wide.to_csv(csv_path, index=False)
    rr = Path(run_root).resolve()
    from dftlearn.visualization.stobe_orbital_energy_summary_figure import (
        write_stobe_orbital_energy_summary_figure,
    )

    summary = write_stobe_orbital_energy_summary_figure(
        rr,
        packaged_output_dir,
        wide,
        dpi=max(dpi, 200),
    )
    if summary is None:
        return None
    return (csv_path, [summary])

