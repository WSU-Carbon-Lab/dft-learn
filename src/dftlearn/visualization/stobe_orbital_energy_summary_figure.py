"""Single-figure StoBe SCF and alpha-level summary for all sites (publication-style)."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd  # noqa: TC002
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator, FuncFormatter
from natsort import natsorted

from dftlearn.io.stobe_orbital_table import (
    analyze_gnd_orbitals,
    element_symbol_from_site,
    parse_stobe_orbital_energies_table,
    reference_k_shell_binding_ev,
)

_LEGEND_FRAME: dict[str, object] = {
    "frameon": True,
    "fancybox": False,
    "facecolor": "white",
    "edgecolor": "0.78",
    "framealpha": 1.0,
    "handlelength": 0.5,
    "markerscale": 0.5,
}


def _apply_tight_negative_energy_ylim(
    ax: plt.Axes,
    arrays: tuple[np.ndarray, ...],
    *,
    pad_frac: float = 0.1,
) -> None:
    """Set y limits from data range with fractional padding (SCF contrast)."""
    parts: list[np.ndarray] = []
    for a in arrays:
        parts.append(np.asarray(a, dtype=float).ravel())
    fin = np.concatenate(parts)
    fin = fin[np.isfinite(fin)]
    if fin.size == 0:
        return
    vmin = float(fin.min())
    vmax = float(fin.max())
    span = vmax - vmin
    if span <= 0.0:
        span = max(abs(vmin), abs(vmax), 1.0) * 0.01
    pad = pad_frac * span
    ax.set_ylim(vmin - pad, vmax + pad)


def _expand_ylim_top(ax: plt.Axes, frac: float = 0.12) -> None:
    """Add empty space at the top of the y-axis (room for in-axes legends)."""
    lo, hi = ax.get_ylim()
    if hi > lo:
        ax.set_ylim(lo, hi + frac * (hi - lo))


def _apply_scf_total_y_tick_format(
    ax: plt.Axes,
    eg: np.ndarray,
    ee: np.ndarray,
    et: np.ndarray,
    *,
    mantissa_decimals: int = 2,
) -> None:
    """Format SCF total-energy y ticks with rounded mantissas and 10^n offset text."""
    parts = [np.asarray(eg).ravel(), np.asarray(ee).ravel(), np.asarray(et).ravel()]
    arr = np.concatenate(parts)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return
    peak = float(np.max(np.abs(arr)))
    if peak < 1e3:
        d = mantissa_decimals

        def _fmt_small(x: float, pos: int) -> str:
            if not np.isfinite(x):
                return ""
            return f"{x:.{d}f}"

        ax.yaxis.set_major_formatter(FuncFormatter(_fmt_small))
        return
    exp = int(np.floor(np.log10(peak)))
    scale = 10.0**exp

    def _fmt(x: float, pos: int) -> str:
        if not np.isfinite(x):
            return ""
        return f"{x / scale:.{mantissa_decimals}f}"

    ax.yaxis.set_major_formatter(FuncFormatter(_fmt))
    ax.yaxis.get_offset_text().set_visible(True)
    ax.yaxis.get_offset_text().set_fontsize(8)
    ax.yaxis.get_offset_text().set_text(rf"$\times 10^{{{exp}}}$")


def _legend_framed(
    ax: plt.Axes,
    *,
    ncol: int,
    loc: str = "upper right",
    bbox_to_anchor: tuple[float, float] | None = None,
    handles: list[object] | None = None,
    **kwargs: object,
) -> None:
    """Add a white boxed legend; pass ``bbox_to_anchor`` for out-of-axes placement."""
    kw: dict[str, object] = {
        "ncol": ncol,
        "loc": loc,
        "fontsize": 8,
        **_LEGEND_FRAME,
    }
    if bbox_to_anchor is not None:
        kw["bbox_to_anchor"] = bbox_to_anchor
    if handles is not None:
        kw["handles"] = handles
    kw.update(kwargs)
    ax.legend(**kw)


def _wide_float_col(wide: pd.DataFrame, sites: list[str], col: str) -> np.ndarray:
    return np.array([float(wide.loc[wide["site"] == s, col].iloc[0]) for s in sites])


def _alpha_ev_at_level(df: pd.DataFrame, level: int) -> float:
    row = df.loc[df["level"] == level]
    if row.empty:
        return float("nan")
    return float(row.iloc[0]["alpha_ev"])


def _style_axes_si(ax: plt.Axes) -> None:
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(
        axis="both",
        which="major",
        direction="in",
        top=True,
        right=True,
        length=3.5,
        width=0.75,
        labelsize=8,
    )
    ax.tick_params(
        axis="both",
        which="minor",
        direction="in",
        top=True,
        right=True,
        length=2.0,
        width=0.5,
        labelsize=8,
    )
    for spine in ax.spines.values():
        spine.set_linewidth(0.75)
    ax.set_axisbelow(True)
    ax.grid(which="major", linestyle="-.", linewidth=0.35, color="0.82", alpha=1.0)


def _save_summary_png(path: Path, fig: plt.Figure, *, dpi: int) -> None:
    fig.savefig(
        path,
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
        edgecolor="none",
        format="png",
    )


def _plot_grouped_totals(
    ax: plt.Axes,
    sites: list[str],
    eg: np.ndarray,
    ee: np.ndarray,
    et: np.ndarray,
) -> None:
    x = np.arange(len(sites), dtype=float)
    w = 0.24
    ec = "0.2"
    ax.bar(x - w, eg, w, label="Ground", color="#4c4c4c", edgecolor=ec, linewidth=0.4)
    ax.bar(x, ee, w, label="Excited", color="#1f77b4", edgecolor=ec, linewidth=0.4)
    ax.bar(x + w, et, w, label="TP", color="#2ca02c", edgecolor=ec, linewidth=0.4)
    ax.set_xticks(x)
    ax.set_xticklabels(sites)
    ax.set_ylabel("SCF total energy (eV)")
    ax.set_title("Total energy (first FINAL ENERGY block)", pad=11)
    _apply_tight_negative_energy_ylim(ax, (eg, ee, et), pad_frac=0.1)
    _expand_ylim_top(ax, frac=0.14)
    _style_axes_si(ax)
    _apply_scf_total_y_tick_format(ax, eg, ee, et, mantissa_decimals=2)
    _legend_framed(ax, ncol=3, loc="upper right", fontsize=7.5)


def _plot_grouped_deltas(
    ax: plt.Axes,
    sites: list[str],
    dks: np.ndarray,
    dexc: np.ndarray,
    dtp: np.ndarray,
    *,
    with_xlabel: bool = False,
) -> None:
    x = np.arange(len(sites), dtype=float)
    w = 0.24
    ec = "0.2"
    kw = {"width": w, "edgecolor": ec, "linewidth": 0.4}
    ax.bar(x - w, dks, label=r"$\Delta$-KS", color="#8c564b", **kw)
    ax.bar(x, dexc, label="Excited $-$ ground", color="#1f77b4", **kw)
    ax.bar(x + w, dtp, label="TP $-$ ground", color="#2ca02c", **kw)
    ax.set_xticks(x)
    ax.set_xticklabels(sites)
    if with_xlabel:
        ax.set_xlabel("Absorption site")
    ax.set_ylabel("Energy (eV)")
    ax.set_title("Spectral / alignment shifts", pad=11)
    ax.axhline(0.0, color="0.45", linestyle="--", linewidth=0.7, zorder=0)
    _apply_tight_negative_energy_ylim(ax, (dks, dexc, dtp), pad_frac=0.08)
    _expand_ylim_top(ax, frac=0.12)
    _style_axes_si(ax)
    _legend_framed(ax, ncol=3, loc="upper right", fontsize=7.5)


def _plot_dumbbells(
    ax: plt.Axes,
    sites: list[str],
    eg: np.ndarray,
    ee: np.ndarray,
    et: np.ndarray,
    *,
    title: str,
    ylabel: str,
    exc_off: float,
    tp_off: float,
    with_xlabel: bool = False,
) -> None:
    """Draw gnd anchor at integer x; exc and TP endpoints horizontally offset."""
    x0 = np.arange(len(sites), dtype=float)
    c_exc = "#1f77b4"
    c_tp = "#2ca02c"
    c_g = "#222222"
    for i in range(len(sites)):
        g, e, t = eg[i], ee[i], et[i]
        xc = float(x0[i])
        if np.isfinite(g):
            ax.scatter(
                [xc],
                [g],
                s=22,
                c=c_g,
                zorder=4,
                edgecolors="0.15",
                linewidths=0.4,
            )
        if np.isfinite(g) and np.isfinite(e):
            ax.plot(
                [xc, xc + exc_off],
                [g, e],
                color=c_exc,
                linewidth=1.75,
                solid_capstyle="round",
                zorder=2,
            )
            ax.scatter(
                [xc + exc_off],
                [e],
                s=20,
                c=c_exc,
                zorder=3,
                edgecolors="0.2",
                linewidths=0.35,
            )
        if np.isfinite(g) and np.isfinite(t):
            ax.plot(
                [xc, xc + tp_off],
                [g, t],
                color=c_tp,
                linewidth=1.75,
                solid_capstyle="round",
                zorder=2,
            )
            ax.scatter(
                [xc + tp_off],
                [t],
                s=20,
                c=c_tp,
                zorder=3,
                edgecolors="0.2",
                linewidths=0.35,
            )
    ax.set_xticks(x0)
    ax.set_xticklabels(sites)
    if with_xlabel:
        ax.set_xlabel("Absorption site")
    ax.set_ylabel(ylabel)
    ax.set_title(title, pad=11)
    _style_axes_si(ax)
    _expand_ylim_top(ax, frac=0.14)
    h_exc = Line2D([0], [0], color=c_exc, lw=2.5, label="vs excited")
    h_tp = Line2D([0], [0], color=c_tp, lw=2.5, label="vs TP")
    h_g = Line2D(
        [0],
        [0],
        marker="o",
        color="w",
        markerfacecolor=c_g,
        markeredgecolor="0.15",
        markersize=7,
        label="Ground anchor",
    )
    _legend_framed(
        ax,
        handles=[h_g, h_exc, h_tp],
        ncol=3,
        loc="upper right",
        fontsize=7.5,
    )


def _plot_frontier_homo_lumo_dumbbells(
    ax: plt.Axes,
    sites: list[str],
    homo_g: np.ndarray,
    homo_e: np.ndarray,
    homo_t: np.ndarray,
    lumo_g: np.ndarray,
    lumo_e: np.ndarray,
    lumo_t: np.ndarray,
    *,
    with_xlabel: bool = False,
) -> None:
    x0 = np.arange(len(sites), dtype=float)
    w = 0.09
    dx = 0.14
    c_g = "#222222"
    for i in range(len(sites)):
        xa, xb = x0[i] - w, x0[i] + w
        hg, he, ht = homo_g[i], homo_e[i], homo_t[i]
        lg, le, lt = lumo_g[i], lumo_e[i], lumo_t[i]
        if np.isfinite(hg):
            ax.scatter(
                [xa],
                [hg],
                s=20,
                c=c_g,
                zorder=4,
                edgecolors="0.15",
                linewidths=0.35,
            )
        if np.isfinite(hg) and np.isfinite(he):
            ax.plot([xa, xa - dx], [hg, he], color="#9467bd", lw=1.65, zorder=2)
            ax.scatter(
                [xa - dx],
                [he],
                s=18,
                c="#9467bd",
                zorder=3,
                edgecolors="0.2",
                linewidths=0.3,
            )
        if np.isfinite(hg) and np.isfinite(ht):
            ax.plot([xa, xa + dx], [hg, ht], color="#bcbd22", lw=1.65, zorder=2)
            ax.scatter(
                [xa + dx],
                [ht],
                s=18,
                c="#bcbd22",
                zorder=3,
                edgecolors="0.2",
                linewidths=0.3,
            )
        if np.isfinite(lg):
            ax.scatter(
                [xb],
                [lg],
                s=20,
                c=c_g,
                zorder=4,
                edgecolors="0.15",
                linewidths=0.35,
            )
        if np.isfinite(lg) and np.isfinite(le):
            ax.plot([xb, xb - dx], [lg, le], color="#1f77b4", lw=1.65, zorder=2)
            ax.scatter(
                [xb - dx],
                [le],
                s=18,
                c="#1f77b4",
                zorder=3,
                edgecolors="0.2",
                linewidths=0.3,
            )
        if np.isfinite(lg) and np.isfinite(lt):
            ax.plot([xb, xb + dx], [lg, lt], color="#2ca02c", lw=1.65, zorder=2)
            ax.scatter(
                [xb + dx],
                [lt],
                s=18,
                c="#2ca02c",
                zorder=3,
                edgecolors="0.2",
                linewidths=0.3,
            )
    ax.set_xticks(x0)
    ax.set_xticklabels(sites)
    if with_xlabel:
        ax.set_xlabel("Absorption site")
    ax.set_ylabel(r"$\alpha$ energy (eV)")
    ax.set_title(
        r"Frontier ($\alpha$ at ground HOMO / LUMO level index in each calculation)",
        pad=11,
    )
    _style_axes_si(ax)
    _expand_ylim_top(ax, frac=0.28)
    h_homo_exc = Line2D([0], [0], color="#9467bd", lw=2.5, label="HOMO vs exc")
    h_homo_tp = Line2D([0], [0], color="#bcbd22", lw=2.5, label="HOMO vs TP")
    h_lumo_exc = Line2D([0], [0], color="#1f77b4", lw=2.5, label="LUMO vs exc")
    h_lumo_tp = Line2D([0], [0], color="#2ca02c", lw=2.5, label="LUMO vs TP")
    h_g = Line2D(
        [0],
        [0],
        marker="o",
        color="w",
        markerfacecolor=c_g,
        markeredgecolor="0.15",
        markersize=6,
        label="Ground anchor",
    )
    _legend_framed(
        ax,
        handles=[h_g, h_homo_exc, h_homo_tp, h_lumo_exc, h_lumo_tp],
        ncol=3,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.0),
        fontsize=7,
    )


def write_stobe_orbital_energy_summary_figure(
    run_root: Path,
    packaged_output_dir: Path,
    wide: pd.DataFrame,
    dpi: int = 200,
) -> Path | None:
    """Build one compact figure: SCF totals, alpha dumbbells, and Delta-KS shifts.

    For each site, reads ``*gnd.out``, ``*exc.out``, ``*tp.out``. Dumbbells use the
    ground-state core / HOMO / LUMO **level indices** and plot ``alpha_ev`` at those
    levels in each calculation (anchor at ground, segments to excited and TP).

    Writes ``stobe_orbital_energy_summary.png`` (white background) under
    ``packaged_output_dir``.

    Parameters
    ----------
    run_root : pathlib.Path
        StoBe run root containing site folders.
    packaged_output_dir : pathlib.Path
        Output directory (e.g. ``packaged_output``).
    wide : pandas.DataFrame
        One row per site from
        :func:`~dftlearn.io.stobe_final_energy.collect_delta_ks_site_table`.
    dpi : int, optional
        PNG resolution.

    Returns
    -------
    pathlib.Path | None
        Path to the PNG if at least one complete site was plotted; otherwise ``None``.
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
    if not sites:
        return None

    eg_ev = _wide_float_col(wide, sites, "E_g_ev")
    ee_ev = _wide_float_col(wide, sites, "E_e_ev")
    et_ev = _wide_float_col(wide, sites, "E_tp_ev")
    dks = _wide_float_col(wide, sites, "E_c_deltaKS_ev")
    dx_exc = _wide_float_col(wide, sites, "delta_exc_vs_gnd_ev")
    dx_tp = _wide_float_col(wide, sites, "delta_tp_vs_gnd_ev")

    core_g: list[float] = []
    core_e: list[float] = []
    core_t: list[float] = []
    homo_g: list[float] = []
    homo_e: list[float] = []
    homo_t: list[float] = []
    lumo_g: list[float] = []
    lumo_e: list[float] = []
    lumo_t: list[float] = []

    for site in sites:
        sym = element_symbol_from_site(site)
        og = parse_stobe_orbital_energies_table(gnd[site])
        oe = parse_stobe_orbital_energies_table(exc[site])
        ot = parse_stobe_orbital_energies_table(tp[site])
        bind = reference_k_shell_binding_ev(sym)
        gnd_r = analyze_gnd_orbitals(og, bind)

        lc, lh, ll = gnd_r.core_level, gnd_r.homo_level, gnd_r.lumo_level
        core_g.append(_alpha_ev_at_level(og, lc))
        core_e.append(_alpha_ev_at_level(oe, lc))
        core_t.append(_alpha_ev_at_level(ot, lc))
        homo_g.append(_alpha_ev_at_level(og, lh))
        homo_e.append(_alpha_ev_at_level(oe, lh))
        homo_t.append(_alpha_ev_at_level(ot, lh))
        lumo_g.append(_alpha_ev_at_level(og, ll))
        lumo_e.append(_alpha_ev_at_level(oe, ll))
        lumo_t.append(_alpha_ev_at_level(ot, ll))

    core_g = np.asarray(core_g, dtype=float)
    core_e = np.asarray(core_e, dtype=float)
    core_t = np.asarray(core_t, dtype=float)
    homo_g = np.asarray(homo_g, dtype=float)
    homo_e = np.asarray(homo_e, dtype=float)
    homo_t = np.asarray(homo_t, dtype=float)
    lumo_g = np.asarray(lumo_g, dtype=float)
    lumo_e = np.asarray(lumo_e, dtype=float)
    lumo_t = np.asarray(lumo_t, dtype=float)

    fig = plt.figure(figsize=(8.2, 9.2))
    gs = fig.add_gridspec(
        4,
        1,
        height_ratios=[1.0, 1.05, 1.05, 0.92],
    )

    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[1, 0], sharex=ax0)
    ax2 = fig.add_subplot(gs[2, 0], sharex=ax0)
    ax3 = fig.add_subplot(gs[3, 0], sharex=ax0)

    _plot_grouped_totals(ax0, sites, eg_ev, ee_ev, et_ev)

    _plot_dumbbells(
        ax1,
        sites,
        core_g,
        core_e,
        core_t,
        title=r"Core $\alpha$ KS level (ground level index)",
        ylabel=r"$\alpha$ energy (eV)",
        exc_off=-0.16,
        tp_off=0.16,
        with_xlabel=False,
    )

    _plot_frontier_homo_lumo_dumbbells(
        ax2,
        sites,
        homo_g,
        homo_e,
        homo_t,
        lumo_g,
        lumo_e,
        lumo_t,
        with_xlabel=False,
    )

    _plot_grouped_deltas(ax3, sites, dks, dx_exc, dx_tp, with_xlabel=True)

    plt.setp(ax0.get_xticklabels(), visible=False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)

    x_right = float(len(sites) - 1) + 0.35
    for ax in (ax0, ax1, ax2, ax3):
        ax.set_xlim(-0.35, x_right)

    fig.subplots_adjust(left=0.11, right=0.92, top=0.95, bottom=0.07, hspace=0.24)

    out = packaged_output_dir / "stobe_orbital_energy_summary.png"
    _save_summary_png(out, fig, dpi=dpi)
    plt.close(fig)
    return out
