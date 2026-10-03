"""OS-elbow and overlap-clustering report tables plus one summary figure."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, ScalarFormatter
from natsort import natsorted

from dftlearn.clustering.iterate import filter_sticks
from dftlearn.clustering.merge import clustering_grid, gaussian_envelope
from dftlearn.clustering.os_elbow import os_percent_elbow
from dftlearn.clustering.selection import select_overlap_threshold
from dftlearn.clustering.types import (
    CLUSTERING_SPEC,
    sticks_from_tp_table,
)
from dftlearn.visualization.pub_style import (
    FIGURE_RC,
    PUB_FS,
    TICK_FS,
    panel_label,
    style_diag_ax,
    style_pub_ax,
)
from dftlearn.xas.spectrum import collect_site_tp_xas

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure
    from matplotlib.transforms import Bbox

    from dftlearn.clustering.types import (
        ClusterResult,
        OverlapSearchResult,
        TransitionSticks,
    )

_INK = (0.13, 0.13, 0.14)
_MUTED = (0.45, 0.45, 0.48)
_ACCENT = (0.90, 0.36, 0.12)
_SITE_RGB: tuple[tuple[float, float, float], ...] = (
    (0.000, 0.447, 0.698),
    (0.835, 0.369, 0.000),
    (0.000, 0.620, 0.451),
    (0.800, 0.475, 0.655),
)
_SQRT_2PI = float(np.sqrt(2.0 * np.pi))
_TENSOR_PANELS: tuple[tuple[str, str, tuple[float, float, float]], ...] = (
    ("iso", r"$I_{\mathrm{iso}}$", (0.15, 0.15, 0.15)),
    ("xx", r"$I_{xx}$", (0.0, 0.447, 0.698)),
    ("yy", r"$I_{yy}$", (0.902, 0.624, 0.0)),
    ("zz", r"$I_{zz}$", (0.0, 0.62, 0.451)),
)
_CARTESIAN_PANELS: tuple[tuple[str, str, tuple[float, float, float]], ...] = (
    _TENSOR_PANELS[1],
    _TENSOR_PANELS[2],
    _TENSOR_PANELS[3],
)


def write_xas_cluster_report(
    run_root: Path,
    packaged_output_dir: Path,
    *,
    xray_filename: str = "XrayT001.out",
    xyz_path: Path | None = None,
    c3_symmetrize: bool = False,
    os_percent: float | None = None,
    n_samples: int = 16,
    lhs_seed: int = 0,
    dpi: int = 150,
) -> tuple[Path, Path, Path, Path, Path] | None:
    """Cluster Delta-KS aligned TP sticks and write tables plus a summary figure.

    Builds clustering sticks from ``collect_site_tp_xas``, chooses an OS%
    cutoff by an elbow test (or ``os_percent``), searches overlap percent
    with a Normal(50%, 10%) sigma ladder and a Gaussian process on BIC,
    and writes membership, effective Gaussians, search diagnostics, and one
    multi-panel PNG. The figure stacks OS and overlap-threshold diagnostics,
    an isotropic envelope, and Cartesian small multiples. Cluster stick
    skylines use OS-weighted site mixture colors keyed on the isotropic panel;
    effective-cluster lines follow the tensor component palette.

    When ``c3_symmetrize`` is True, Cartesian OS components are C3-folded in
    the Al-N/O molecular frame before filtering and overlap clustering.

    Parameters
    ----------
    run_root : pathlib.Path
        StoBe run directory with site folders and ``{site}.xas`` sticks.
    packaged_output_dir : pathlib.Path
        Output folder (created if missing).
    xray_filename : str, optional
        StoBe table used only to assemble the stick catalog.
    xyz_path : pathlib.Path, optional
        Geometry required when ``c3_symmetrize`` is True.
    c3_symmetrize : bool, optional
        Fold Cartesian dipoles under C3 before clustering.
    os_percent : float, optional
        If given, this OS% of windowed max OS is used instead of the elbow.
    n_samples : int, optional
        Unique overlap percents to evaluate in the GP search.
    lhs_seed : int, optional
        Seed reserved for Gaussian-process tie-breaking when the sigma
        ladder is shorter than ``n_samples``.
    dpi : int, optional
        PNG resolution.

    Returns
    -------
    tuple of pathlib.Path, or None
        Members CSV, effective-cluster CSV, search CSV, OS-elbow CSV,
        and summary PNG, or ``None`` if no TP sticks exist.

    Raises
    ------
    ValueError
        If sticks exist but filtering or clustering cannot run.
    """
    run_root = Path(run_root).resolve()
    packaged_output_dir = Path(packaged_output_dir).resolve()
    try:
        _energy, _spectra, _metrics, sticks_df, _c3_frame = collect_site_tp_xas(
            run_root,
            xray_filename=xray_filename,
            xyz_path=xyz_path,
            c3_symmetrize=c3_symmetrize,
        )
    except FileNotFoundError:
        return None
    if sticks_df.empty:
        return None
    sticks = sticks_from_tp_table(sticks_df)
    elbow_pct, percents, n_kept, frac = os_percent_elbow(
        sticks,
        energy_min_ev=CLUSTERING_SPEC.energy_min_ev,
        energy_max_ev=CLUSTERING_SPEC.energy_max_ev,
    )
    chosen_os = float(elbow_pct if os_percent is None else os_percent)
    filtered, _kept = filter_sticks(
        sticks,
        os_percent=chosen_os,
        energy_min_ev=CLUSTERING_SPEC.energy_min_ev,
        energy_max_ev=CLUSTERING_SPEC.energy_max_ev,
    )
    search, result = select_overlap_threshold(
        filtered,
        n_samples=n_samples,
        seed=lhs_seed,
        settings=CLUSTERING_SPEC,
    )
    packaged_output_dir.mkdir(parents=True, exist_ok=True)
    members_path = packaged_output_dir / "xas_cluster_members.csv"
    clusters_path = packaged_output_dir / "xas_effective_clusters.csv"
    search_path = packaged_output_dir / "xas_cluster_search.csv"
    elbow_csv_path = packaged_output_dir / "xas_cluster_os_elbow.csv"
    summary_path = packaged_output_dir / "xas_cluster_summary.png"
    _members_frame(filtered, result).to_csv(members_path, index=False)
    _clusters_frame(result).to_csv(clusters_path, index=False)
    _search_frame(search).to_csv(search_path, index=False)
    _elbow_frame(percents, n_kept, frac, chosen_os).to_csv(elbow_csv_path, index=False)
    _write_summary_figure(
        filtered,
        result,
        search,
        percents,
        n_kept,
        chosen_os,
        summary_path,
        catalog=sticks,
        dpi=dpi,
    )
    return (
        members_path,
        clusters_path,
        search_path,
        elbow_csv_path,
        summary_path,
    )


def _members_frame(sticks: TransitionSticks, result: ClusterResult) -> pd.DataFrame:
    cluster_id = _member_cluster_ids(sticks.energy_ev.shape[0], result.member_indices)
    return pd.DataFrame(
        {
            "site": sticks.site,
            "energy_ev": sticks.energy_ev,
            "oscillator_strength": sticks.oscillator_strength,
            "sigma_ev": sticks.sigma_ev,
            "os_xx": sticks.os_xx,
            "os_yy": sticks.os_yy,
            "os_zz": sticks.os_zz,
            "theta_deg": _theta_deg(sticks.os_xx, sticks.os_yy, sticks.os_zz),
            "cluster_id": cluster_id,
        }
    )


def _clusters_frame(result: ClusterResult) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "cluster_id": np.arange(result.energy_ev.shape[0], dtype=np.int64),
            "energy_ev": result.energy_ev,
            "amplitude": result.amplitude,
            "sigma_ev": result.sigma_ev,
            "os_xx": result.os_xx,
            "os_yy": result.os_yy,
            "os_zz": result.os_zz,
            "theta_deg": result.theta_deg,
            "n_members": result.n_members,
        }
    )


def _search_frame(search: OverlapSearchResult) -> pd.DataFrame:
    selected = np.zeros(search.overlap_percent.shape[0], dtype=np.int64)
    selected[search.selected_index] = 1
    return pd.DataFrame(
        {
            "overlap_percent": search.overlap_percent,
            "n_clusters": search.n_clusters,
            "rss": search.rss,
            "bic": search.bic,
            "selected": selected,
            "lhs_seed": np.full(search.overlap_percent.shape[0], search.lhs_seed),
        }
    )


def _elbow_frame(
    percents: np.ndarray,
    n_kept: np.ndarray,
    frac: np.ndarray,
    chosen_os: float,
) -> pd.DataFrame:
    selected = (np.abs(percents - chosen_os) < 1e-9).astype(np.int64)
    if not np.any(selected):
        selected[int(np.argmin(np.abs(percents - chosen_os)))] = 1
    return pd.DataFrame(
        {
            "os_percent": percents,
            "n_kept": n_kept,
            "retained_os_fraction": frac,
            "selected": selected,
        }
    )


def _member_cluster_ids(
    n: int,
    member_indices: tuple[tuple[int, ...], ...],
) -> np.ndarray:
    cluster_id = np.full(n, -1, dtype=np.int64)
    for cid, members in enumerate(member_indices):
        if not members:
            continue
        cluster_id[np.asarray(members, dtype=np.int64)] = cid
    return cluster_id


def _theta_deg(xx: np.ndarray, yy: np.ndarray, zz: np.ndarray) -> np.ndarray:
    xx = np.asarray(xx, dtype=np.float64)
    yy = np.asarray(yy, dtype=np.float64)
    zz = np.asarray(zz, dtype=np.float64)
    mag = np.sqrt(xx * xx + yy * yy + zz * zz)
    theta = np.full(xx.shape[0], np.nan, dtype=np.float64)
    ok = mag > 0.0
    theta[ok] = np.degrees(np.arccos(np.clip(zz[ok] / mag[ok], -1.0, 1.0)))
    return theta


def _site_names(sites: np.ndarray) -> tuple[str, ...]:
    return tuple(natsorted(np.unique(np.asarray(sites).astype(str))))


def _site_palette(sites: np.ndarray) -> dict[str, tuple[float, float, float]]:
    names = _site_names(sites)
    return {name: _SITE_RGB[i % len(_SITE_RGB)] for i, name in enumerate(names)}


def _cluster_mixture_rgba(
    sticks: TransitionSticks,
    result: ClusterResult,
    palette: dict[str, tuple[float, float, float]],
) -> np.ndarray:
    k = int(result.energy_ev.shape[0])
    out = np.full((k, 4), 0.45, dtype=np.float64)
    out[:, 3] = 1.0
    weights = np.nan_to_num(
        np.asarray(sticks.oscillator_strength, dtype=np.float64),
        nan=0.0,
    )
    labels = np.asarray(sticks.site).astype(str)
    for cid, members in enumerate(result.member_indices):
        if not members:
            continue
        idx = np.asarray(members, dtype=np.int64)
        w = np.maximum(weights[idx], 0.0)
        total = float(np.sum(w))
        rgb = np.zeros(3, dtype=np.float64)
        if total <= 0.0:
            continue
        for name, color in palette.items():
            share = float(np.sum(w[labels[idx] == name]))
            rgb += share * np.asarray(color, dtype=np.float64)
        out[cid, 0:3] = rgb / total
    return out


def _site_legend_handles(
    palette: dict[str, tuple[float, float, float]],
) -> list[Line2D]:
    return [
        Line2D(
            [0],
            [0],
            color=palette[name],
            lw=2.4,
            solid_capstyle="butt",
            label=name,
        )
        for name in palette
    ]


def _draw_site_color_key(
    ax: Axes,
    palette: dict[str, tuple[float, float, float]],
) -> None:
    """Attach a site-color key explaining cluster stick mixture colors."""
    if not palette:
        return
    key = ax.legend(
        handles=_site_legend_handles(palette),
        loc="lower right",
        title="Cluster sticks: OS-weighted site mix",
        fontsize=8,
        title_fontsize=8,
        frameon=True,
        framealpha=0.94,
        edgecolor="0.72",
        fancybox=False,
        handlelength=1.4,
        handletextpad=0.35,
        borderaxespad=0.35,
        labelcolor="0.28",
        ncol=min(len(palette), 4),
    )
    ax.add_artist(key)


def _bbox_overlap(a: Bbox, b: Bbox, *, gap: float = 0.0) -> bool:
    return not (
        a.x1 <= b.x0 + gap
        or b.x1 <= a.x0 + gap
        or a.y1 <= b.y0 + gap
        or b.y1 <= a.y0 + gap
    )


def _overlapping_axes(
    axes: tuple[Axes, ...],
    *,
    glued: frozenset[tuple[int, int]],
) -> list[tuple[int, int]]:
    hits: list[tuple[int, int]] = []
    boxes = [ax.get_position() for ax in axes]
    for i in range(len(axes)):
        for j in range(i + 1, len(axes)):
            if (i, j) in glued or (j, i) in glued:
                continue
            if _bbox_overlap(boxes[i], boxes[j]):
                hits.append((i, j))
    return hits


def _plot_energy_window(
    sticks: TransitionSticks,
    result: ClusterResult,
) -> tuple[float, float]:
    energy = np.concatenate(
        (
            np.asarray(sticks.energy_ev, dtype=np.float64),
            np.asarray(result.energy_ev, dtype=np.float64),
        )
    )
    sigma = np.concatenate(
        (
            np.asarray(sticks.sigma_ev, dtype=np.float64),
            np.asarray(result.sigma_ev, dtype=np.float64),
        )
    )
    amp = np.concatenate(
        (
            np.asarray(sticks.oscillator_strength, dtype=np.float64),
            np.asarray(result.amplitude, dtype=np.float64),
        )
    )
    ok = np.isfinite(energy) & np.isfinite(sigma) & np.isfinite(amp) & (amp > 0.0)
    if not np.any(ok):
        return CLUSTERING_SPEC.energy_min_ev, CLUSTERING_SPEC.energy_max_ev
    energy = energy[ok]
    width = np.maximum(sigma[ok], 0.05)
    e_lo = float(np.min(energy - 2.5 * width))
    e_hi = float(np.max(energy + 3.0 * width))
    e_lo = max(e_lo, CLUSTERING_SPEC.energy_min_ev)
    e_hi = min(e_hi, CLUSTERING_SPEC.energy_max_ev)
    if e_hi - e_lo < 2.0:
        mid = 0.5 * (e_lo + e_hi)
        e_lo = mid - 1.0
        e_hi = mid + 1.0
    return e_lo, e_hi


def _component_amplitude(
    xx: np.ndarray,
    yy: np.ndarray,
    zz: np.ndarray,
    which: str,
) -> np.ndarray:
    xx = np.nan_to_num(np.asarray(xx, dtype=np.float64), nan=0.0)
    yy = np.nan_to_num(np.asarray(yy, dtype=np.float64), nan=0.0)
    zz = np.nan_to_num(np.asarray(zz, dtype=np.float64), nan=0.0)
    match which:
        case "xx":
            return xx
        case "yy":
            return yy
        case "zz":
            return zz
        case "iso":
            return (xx + yy + zz) / 3.0
        case _:
            msg = f"unknown tensor component {which!r}"
            raise ValueError(msg)


def _align_ylabels(fig: Figure, axes: tuple[Axes, ...]) -> None:
    labeled = tuple(ax for ax in axes if ax.get_ylabel())
    if labeled:
        fig.align_ylabels(labeled)


def _write_summary_figure(
    sticks: TransitionSticks,
    result: ClusterResult,
    search: OverlapSearchResult,
    percents: np.ndarray,
    n_kept: np.ndarray,
    chosen_os: float,
    path: Path,
    *,
    catalog: TransitionSticks,
    dpi: int,
) -> None:
    fig, _axes, _glued = _build_summary_figure(
        sticks,
        result,
        search,
        percents,
        n_kept,
        chosen_os,
        catalog=catalog,
    )
    fig.savefig(path, dpi=dpi, facecolor="white")
    plt.close(fig)


def _build_summary_figure(
    sticks: TransitionSticks,
    result: ClusterResult,
    search: OverlapSearchResult,
    percents: np.ndarray,
    n_kept: np.ndarray,
    chosen_os: float,
    *,
    catalog: TransitionSticks,
) -> tuple[Figure, tuple[Axes, ...], frozenset[tuple[int, int]]]:
    _ = catalog
    with plt.rc_context(FIGURE_RC):
        fig = plt.figure(figsize=(7.2, 6.6), facecolor="white")
        gs = GridSpec(
            3,
            1,
            figure=fig,
            height_ratios=[0.58, 1.20, 0.82],
            hspace=0.34,
            left=0.13,
            right=0.97,
            top=0.98,
            bottom=0.10,
        )
        gs_diag = gs[0].subgridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.30)
        ax_os = fig.add_subplot(gs_diag[0, 0])
        ax_ovp = fig.add_subplot(gs_diag[0, 1])
        ax_iso = fig.add_subplot(gs[1])
        gs_cart = gs[2].subgridspec(1, 3, wspace=0.14)
        ax_xx = fig.add_subplot(gs_cart[0, 0])
        ax_yy = fig.add_subplot(gs_cart[0, 1], sharex=ax_xx)
        ax_zz = fig.add_subplot(gs_cart[0, 2], sharex=ax_xx)
        e_lo, e_hi = _plot_energy_window(sticks, result)
        palette = _site_palette(sticks.site)
        mixture = _cluster_mixture_rgba(sticks, result, palette)
        _draw_os_elbow_panel(ax_os, percents, n_kept, chosen_os)
        _draw_ovp_search_panel(ax_ovp, search)
        grid = clustering_grid(CLUSTERING_SPEC)
        iso_spec = _TENSOR_PANELS[0]
        _draw_tensor_compare(
            ax_iso,
            grid,
            sticks,
            result,
            which=iso_spec[0],
            cluster_color=iso_spec[2],
            energy_limits=(e_lo, e_hi),
            mixture=mixture,
            show_legend=True,
            show_xlabel=True,
            show_ylabel=True,
            caption=iso_spec[1],
            normalize=True,
        )
        _draw_site_color_key(ax_iso, palette)
        cart_spec = (ax_xx, ax_yy, ax_zz)
        for i, (key, title, rgb) in enumerate(_CARTESIAN_PANELS):
            _draw_tensor_compare(
                cart_spec[i],
                grid,
                sticks,
                result,
                which=key,
                cluster_color=rgb,
                energy_limits=(e_lo, e_hi),
                mixture=mixture,
                show_legend=False,
                show_xlabel=False,
                show_ylabel=i == 0,
                caption=title,
                caption_ha="right",
                normalize=True,
            )
            cart_spec[i].xaxis.set_major_locator(MaxNLocator(nbins=4, integer=False))
            cart_spec[i].yaxis.set_major_locator(MaxNLocator(nbins=3))
            if i > 0:
                cart_spec[i].tick_params(axis="y", labelleft=False)
        _style_cart_row_xaxis(cart_spec)
        fig.supxlabel("Photon energy (eV)", fontsize=PUB_FS, y=0.03)
        letter_axes: tuple[tuple[Axes, str], ...] = (
            (ax_os, "a"),
            (ax_ovp, "b"),
            (ax_iso, "c"),
            (ax_xx, "d"),
        )
        fig.canvas.draw()
        _align_ylabels(fig, (ax_iso, *cart_spec))
        _align_ylabels(fig, (ax_os, ax_ovp))
        _style_cart_row_xaxis(cart_spec)
        fig.canvas.draw()
        for ax, letter in letter_axes:
            panel_label(ax, letter)
        panel_axes = (ax_os, ax_ovp, ax_iso, *cart_spec)
        return fig, panel_axes, frozenset()


def _choice_marker(ax: Axes, x: float, label: str) -> None:
    ax.axvline(x, color=_ACCENT, ls="-", lw=0.9, zorder=4)
    ha = "left" if x <= 50.0 else "right"
    dx = 1.8 if ha == "left" else -1.8
    ax.text(
        x + dx,
        0.90,
        label,
        transform=ax.get_xaxis_transform(),
        ha=ha,
        va="top",
        fontsize=7,
        color=_ACCENT,
        clip_on=True,
    )


def _draw_os_elbow_panel(
    ax: Axes,
    percents: np.ndarray,
    n_kept: np.ndarray,
    chosen_os: float,
) -> None:
    y = np.asarray(n_kept, dtype=np.float64)
    ax.plot(percents, y, color=_INK, lw=1.35, solid_capstyle="round")
    idx = int(np.argmin(np.abs(percents - chosen_os)))
    _choice_marker(ax, chosen_os, f"{chosen_os:.0f}%  ·  {int(y[idx])} sticks")
    ax.set_xlim(0.0, 100.0)
    ax.set_ylim(0.0, max(float(np.max(y)) * 1.08, 1.0))
    ax.set_xlabel("Oscillator cutoff (%)", fontsize=PUB_FS, labelpad=2)
    ax.set_ylabel("Sticks kept", fontsize=PUB_FS)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
    style_diag_ax(ax)


def _draw_ovp_search_panel(ax: Axes, search: OverlapSearchResult) -> None:
    ovp = search.overlap_percent
    bic = search.bic
    n_cl = np.asarray(search.n_clusters, dtype=np.float64)
    ax.plot(ovp, bic, color=_INK, lw=1.35, label="BIC", solid_capstyle="round")
    ax.scatter(
        ovp,
        bic,
        s=16.0,
        facecolors="white",
        edgecolors=_INK,
        linewidths=0.7,
        zorder=5,
    )
    ax2 = ax.twinx()
    ax2.plot(
        ovp,
        n_cl,
        color=_MUTED,
        lw=1.15,
        label="Clusters",
        solid_capstyle="round",
    )
    ax2.scatter(
        ovp,
        n_cl,
        s=16.0,
        facecolors="white",
        edgecolors=_MUTED,
        linewidths=0.7,
        zorder=5,
    )
    n_sel = int(search.n_clusters[search.selected_index])
    _choice_marker(
        ax,
        search.selected_overlap,
        f"{search.selected_overlap:.0f}%  ·  {n_sel} clusters",
    )
    ax.set_xlim(0.0, 100.0)
    bic_fmt = ScalarFormatter(useMathText=True)
    bic_fmt.set_powerlimits((-3, 4))
    ax.yaxis.set_major_formatter(bic_fmt)
    ax.set_ylabel("BIC", fontsize=PUB_FS)
    ax2.set_ylabel("Clusters", fontsize=PUB_FS)
    ax2.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
    ax.set_xlabel("Overlap threshold (%)", fontsize=PUB_FS, labelpad=2)
    legend = ax.legend(
        handles=[
            Line2D([0], [0], color=_INK, lw=1.35, label="BIC"),
            Line2D([0], [0], color=_MUTED, lw=1.15, label="Clusters"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, 1.12),
        frameon=True,
        facecolor="white",
        framealpha=1.0,
        edgecolor="0.72",
        fancybox=False,
        fontsize=8,
        handlelength=1.4,
        borderaxespad=0.0,
        labelcolor="black",
    )
    legend.set_clip_on(False)
    offset = ax.yaxis.get_offset_text()
    offset.set_color("black")
    offset.set_bbox(
        {"facecolor": "white", "edgecolor": "none", "pad": 1.0, "alpha": 1.0}
    )
    style_diag_ax(ax, right=True)
    ax2.tick_params(axis="y", labelsize=TICK_FS)


def _cluster_bar_heights(
    component_amp: np.ndarray,
    sigma_ev: np.ndarray,
    scale: float,
) -> np.ndarray:
    amp = np.asarray(component_amp, dtype=np.float64)
    sigma = np.asarray(sigma_ev, dtype=np.float64)
    peak = amp / (sigma * _SQRT_2PI)
    if scale > 0.0 and np.isfinite(scale):
        peak = peak / scale
    return peak


def _draw_cluster_os_bars(
    ax: Axes,
    result: ClusterResult,
    component_amp: np.ndarray,
    mixture: np.ndarray,
    scale: float,
) -> None:
    energy = np.asarray(result.energy_ev, dtype=np.float64)
    heights = _cluster_bar_heights(component_amp, result.sigma_ev, scale)
    ok = (
        np.isfinite(energy)
        & np.isfinite(heights)
        & (heights > 0.0)
        & np.isfinite(mixture[:, 0])
    )
    energy = energy[ok]
    heights = heights[ok]
    face = mixture[ok]
    if energy.size == 0:
        return
    segs = np.stack(
        (
            np.column_stack((energy, np.zeros_like(energy))),
            np.column_stack((energy, heights)),
        ),
        axis=1,
    )
    ax.add_collection(
        LineCollection(
            segs.tolist(),
            colors=face.tolist(),
            linewidths=1.0,
            alpha=0.92,
            zorder=2,
        )
    )


def _draw_tensor_compare(
    ax: Axes,
    grid: np.ndarray,
    sticks: TransitionSticks,
    result: ClusterResult,
    *,
    which: str,
    cluster_color: tuple[float, float, float],
    energy_limits: tuple[float, float],
    mixture: np.ndarray,
    show_legend: bool,
    show_xlabel: bool,
    show_ylabel: bool,
    caption: str,
    caption_ha: str = "left",
    normalize: bool,
) -> None:
    dft_amp = _component_amplitude(sticks.os_xx, sticks.os_yy, sticks.os_zz, which)
    cl_amp = _component_amplitude(result.os_xx, result.os_yy, result.os_zz, which)
    dft = gaussian_envelope(grid, sticks.energy_ev, dft_amp, sticks.sigma_ev)
    clustered = gaussian_envelope(grid, result.energy_ev, cl_amp, result.sigma_ev)
    scale = 1.0
    if normalize:
        scale = float(np.nanmax(np.abs(dft))) if dft.size else 1.0
        if not np.isfinite(scale) or scale <= 0.0:
            scale = 1.0
    dft_draw = dft / scale
    cl_draw = clustered / scale
    ax.fill_between(grid, 0.0, dft_draw, color="0.90", linewidth=0.0, zorder=1)
    e_lo, e_hi = energy_limits
    ymax = float(np.nanmax(np.concatenate((dft_draw, cl_draw)))) if dft.size else 1.0
    if not np.isfinite(ymax) or ymax <= 0.0:
        ymax = 1.0
    _draw_cluster_os_bars(ax, result, cl_amp, mixture, scale)
    ax.plot(
        grid,
        dft_draw,
        color=_INK,
        lw=1.45,
        label="OS-filtered DFT",
        zorder=3,
        solid_capstyle="round",
    )
    ax.plot(
        grid,
        cl_draw,
        color=cluster_color,
        lw=1.35,
        ls="--",
        label="Effective clusters",
        zorder=4,
        solid_capstyle="round",
    )
    ax.set_xlim(e_lo, e_hi)
    ax.set_ylim(0.0, ymax * 1.06)
    ax.autoscale(enable=False)
    if normalize and ymax <= 1.15:
        ax.set_yticks([0.0, 0.5, 1.0])
    if show_ylabel:
        ax.set_ylabel("Abs. (arb. units)", fontsize=PUB_FS)
    if show_xlabel:
        ax.set_xlabel("Photon energy (eV)", fontsize=PUB_FS, labelpad=3)
    if caption:
        cap_x = 0.02 if caption_ha == "left" else 0.98
        ax.text(
            cap_x,
            0.86,
            caption,
            transform=ax.transAxes,
            ha=caption_ha,
            va="top",
            fontsize=PUB_FS,
            color="0.28",
        )
    if show_legend:
        spec_legend = ax.legend(
            loc="upper right",
            frameon=True,
            framealpha=0.94,
            edgecolor="0.72",
            fancybox=False,
            fontsize=8,
            handlelength=1.8,
            borderaxespad=0.2,
            labelcolor="0.28",
        )
        ax.add_artist(spec_legend)
    style_pub_ax(ax)


def _style_cart_row_xaxis(axes: tuple[Axes, ...]) -> None:
    """Show shared energy tick labels on every Cartesian small multiple."""
    for ax in axes:
        ax.tick_params(axis="x", labelbottom=True, labeltop=False)
        ax.tick_params(axis="y", labelright=False)
        for tick in ax.get_xticklabels():
            tick.set_visible(True)
