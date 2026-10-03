"""Layout and site-color tests for the XAS cluster summary figure."""

from __future__ import annotations

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.collections import LineCollection
from matplotlib.legend import Legend
from matplotlib.text import Text

from dftlearn.clustering.types import (
    ClusterResult,
    ClusterStage,
    OverlapSearchResult,
    TransitionSticks,
)
from dftlearn.visualization.xas_cluster_figure import (
    _SITE_RGB,
    _build_summary_figure,
    _cluster_mixture_rgba,
    _overlapping_axes,
    _site_palette,
)


def _toy_sticks(n: int = 8) -> TransitionSticks:
    energy = np.linspace(284.5, 291.5, n)
    sites = np.array(["C1", "C2", "C3", "C4"] * (n // 4), dtype=object)
    amp = np.linspace(0.9, 0.25, n)
    return TransitionSticks(
        energy_ev=energy,
        oscillator_strength=amp,
        sigma_ev=np.full(n, 0.22),
        site=sites,
        os_xx=amp * 0.2,
        os_yy=amp * 0.3,
        os_zz=amp * 0.5,
    )


def _toy_result(sticks: TransitionSticks) -> ClusterResult:
    n = int(sticks.energy_ev.shape[0])
    members = (
        tuple(range(n // 2)),
        tuple(range(n // 2, n)),
    )
    k = len(members)
    ov_n = np.eye(n, dtype=np.float64) * 100.0
    ov_k = np.eye(k, dtype=np.float64) * 100.0
    energy = np.array(
        [float(np.mean(sticks.energy_ev[list(m)])) for m in members],
        dtype=np.float64,
    )
    amp = np.array(
        [float(np.sum(sticks.oscillator_strength[list(m)])) for m in members],
        dtype=np.float64,
    )
    xx = np.array([float(np.sum(sticks.os_xx[list(m)])) for m in members])
    yy = np.array([float(np.sum(sticks.os_yy[list(m)])) for m in members])
    zz = np.array([float(np.sum(sticks.os_zz[list(m)])) for m in members])
    mag = np.sqrt(xx * xx + yy * yy + zz * zz)
    theta = np.degrees(np.arccos(np.clip(zz / mag, -1.0, 1.0)))
    return ClusterResult(
        energy_ev=energy,
        amplitude=amp,
        sigma_ev=np.full(k, 0.35),
        os_xx=xx,
        os_yy=yy,
        os_zz=zz,
        theta_deg=theta,
        n_members=np.array([len(m) for m in members], dtype=np.int64),
        member_indices=members,
        n_iterations=1,
        overlap_threshold=50.0,
        final_overlap=ov_k,
        stages=(
            ClusterStage(
                label="OS filtered",
                overlap=ov_n,
                n_peaks=n,
                max_offdiag=0.0,
            ),
            ClusterStage(
                label="Final",
                overlap=ov_k,
                n_peaks=k,
                max_offdiag=0.0,
            ),
        ),
    )


def _toy_search() -> OverlapSearchResult:
    ovp = np.array([0.0, 40.0, 50.0, 90.0])
    return OverlapSearchResult(
        overlap_percent=ovp,
        n_clusters=np.array([8, 5, 4, 2], dtype=np.int64),
        rss=np.array([1.0, 0.6, 0.5, 0.8]),
        bic=np.array([-21000.0, -20500.0, -20000.0, -19800.0]),
        selected_overlap=50.0,
        selected_index=2,
        lhs_seed=0,
    )


def test_cluster_mixture_is_os_weighted() -> None:
    """Mixture RGB is the OS-weighted mean of member site colors."""
    sticks = TransitionSticks(
        energy_ev=np.array([285.0, 285.2, 287.0]),
        oscillator_strength=np.array([1.0, 3.0, 2.0]),
        sigma_ev=np.full(3, 0.2),
        site=np.array(["C1", "C2", "C1"], dtype=object),
        os_xx=np.ones(3),
        os_yy=np.zeros(3),
        os_zz=np.ones(3),
    )
    result = ClusterResult(
        energy_ev=np.array([285.1, 287.0]),
        amplitude=np.array([4.0, 2.0]),
        sigma_ev=np.array([0.3, 0.3]),
        os_xx=np.ones(2),
        os_yy=np.zeros(2),
        os_zz=np.ones(2),
        theta_deg=np.array([45.0, 45.0]),
        n_members=np.array([2, 1], dtype=np.int64),
        member_indices=((0, 1), (2,)),
        n_iterations=1,
        overlap_threshold=40.0,
        final_overlap=np.eye(2) * 100.0,
        stages=(
            ClusterStage(
                label="Final",
                overlap=np.eye(2) * 100.0,
                n_peaks=2,
                max_offdiag=0.0,
            ),
        ),
    )
    palette = _site_palette(sticks.site)
    mix = _cluster_mixture_rgba(sticks, result, palette)
    c1 = np.asarray(palette["C1"], dtype=np.float64)
    c2 = np.asarray(palette["C2"], dtype=np.float64)
    expected = (1.0 * c1 + 3.0 * c2) / 4.0
    np.testing.assert_allclose(mix[0, 0:3], expected)
    np.testing.assert_allclose(mix[1, 0:3], c1)
    assert palette["C1"] == _SITE_RGB[0]
    assert palette["C2"] == _SITE_RGB[1]


def test_summary_figure_layout_and_site_colors() -> None:
    """Summary figure has no panel overlap, cluster OS bars, and no theta panel."""
    sticks = _toy_sticks()
    result = _toy_result(sticks)
    n_clusters = int(result.energy_ev.shape[0])
    percents = np.linspace(0.0, 100.0, 6)
    n_kept = np.array([8, 7, 5, 4, 3, 2], dtype=np.int64)
    fig, axes, glued = _build_summary_figure(
        sticks,
        result,
        _toy_search(),
        percents,
        n_kept,
        20.0,
        catalog=sticks,
    )
    try:
        fig.canvas.draw()
        hits = _overlapping_axes(axes, glued=glued)
        assert hits == []
        ax_os, ax_ovp, ax_iso, ax_xx, ax_yy, ax_zz = axes
        assert len(axes) == 6
        assert glued == frozenset()
        assert ax_os.get_ylabel() == "Sticks kept"
        assert ax_os.get_yscale() == "linear"
        legend_text = ax_ovp.get_legend().get_texts()
        legend_labels = {text.get_text() for text in legend_text}
        assert legend_labels == {"BIC", "Clusters"}
        assert not any(label.startswith("_") for label in legend_labels)
        site_legends = list(
            {
                id(leg): leg
                for leg in ax_iso.get_children()
                if isinstance(leg, Legend)
            }.values()
        )
        assert len(site_legends) == 2
        site_key = next(
            leg
            for leg in site_legends
            if leg.get_title().get_text().startswith("Cluster sticks")
        )
        site_labels = {text.get_text() for text in site_key.get_texts()}
        assert site_labels == {"C1", "C2", "C3", "C4"}
        assert ax_iso.get_ylabel() == "Abs. (arb. units)"
        assert ax_xx.get_ylabel() == "Abs. (arb. units)"
        assert ax_yy.get_ylabel() == ""
        assert ax_zz.get_ylabel() == ""
        abs_labels = [
            ax.get_ylabel() for ax in fig.axes if ax.get_ylabel() == "Abs. (arb. units)"
        ]
        assert len(abs_labels) == 2
        n_theta = sum(1 for ax in fig.axes if r"\theta" in ax.get_ylabel())
        assert n_theta == 0
        n_imshow = sum(
            1
            for ax in axes
            if any(artist.get_label() == "<image>" for artist in ax.get_children())
        )
        assert n_imshow == 0
        letters = {
            child.get_text(): child
            for ax in (ax_os, ax_ovp, ax_iso, ax_xx)
            for child in ax.get_children()
            if isinstance(child, Text)
            and child.get_text().startswith("(")
            and len(child.get_text()) == 3
        }
        assert set(letters) == {f"({c})" for c in "abcd"}
        for child in letters.values():
            assert child.get_color() == "black"
        letter_a = letters["(a)"]
        ax_pos = letter_a.get_position()
        assert ax_pos[0] == pytest.approx(0.02)
        assert ax_pos[1] == pytest.approx(0.98)
        assert ax_iso.get_xlabel() == "Photon energy (eV)"
        assert ax_xx.get_xlabel() == ""
        assert ax_yy.get_xlabel() == ""
        assert ax_zz.get_xlabel() == ""
        fig.canvas.draw()
        for ax in (ax_xx, ax_yy, ax_zz):
            labels = [t for t in ax.get_xticklabels() if t.get_text()]
            assert labels
            assert all(t.get_visible() for t in labels)
        ovp_legend = ax_ovp.get_legend()
        assert ovp_legend.get_frame().get_facecolor()[0] == pytest.approx(1.0, abs=0.05)
        from matplotlib.ticker import ScalarFormatter

        assert isinstance(ax_ovp.yaxis.get_major_formatter(), ScalarFormatter)
        bic_tick_text = "".join(tick.get_text() for tick in ax_ovp.get_yticklabels())
        assert "10" in bic_tick_text
        bars = [
            artist
            for artist in ax_iso.collections
            if isinstance(artist, LineCollection)
            and np.asarray(artist.get_segments()).shape[0] == n_clusters
        ]
        assert len(bars) == 1
        for ax in (ax_xx, ax_yy, ax_zz):
            panel_bars = [
                artist
                for artist in ax.collections
                if isinstance(artist, LineCollection)
                and np.asarray(artist.get_segments()).shape[0] == n_clusters
            ]
            assert len(panel_bars) == 1
        palette = _site_palette(sticks.site)
        mixture = _cluster_mixture_rgba(sticks, result, palette)
        face = np.asarray(bars[0].get_colors())
        assert face.shape[0] == n_clusters
        np.testing.assert_allclose(face[:, 0:3], mixture[:, 0:3], rtol=0.0, atol=1e-6)
    finally:
        plt.close(fig)
