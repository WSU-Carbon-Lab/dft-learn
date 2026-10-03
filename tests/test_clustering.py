"""Tests for Igor-style energetic overlap clustering."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.special import erf

from dftlearn.clustering.grouping import simple_cluster3
from dftlearn.clustering.iterate import cluster_by_overlap, filter_sticks
from dftlearn.clustering.merge import (
    clustering_grid,
    effective_gaussian,
)
from dftlearn.clustering.os_elbow import kneedle_percent, os_percent_elbow
from dftlearn.clustering.overlap import (
    max_offdiag_overlap,
    overlap_matrix,
    percent_overlap,
)
from dftlearn.clustering.selection import (
    DEFAULT_PRIOR_SIGMA,
    _initial_ovps,
    _pick_trace_crossing,
    _sigma_ladder,
    cluster_bic,
    select_overlap_threshold,
)
from dftlearn.clustering.types import (
    TransitionSticks,
    sticks_from_tp_table,
)

_SQRT2 = float(np.sqrt(2.0))
_EXAMPLE = (
    Path(__file__).resolve().parents[1] / "docs" / "stobe" / "example-run" / "znpc-cif"
)


def _independent_unit_overlap(mu1: float, sd1: float, mu2: float, sd2: float) -> float:
    if sd1 == sd2:
        c = 0.5 * (mu1 + mu2)
    elif mu1 < mu2:
        disc = (mu1 - mu2) ** 2 + 2.0 * (sd1**2 - sd2**2) * np.log(sd1 / sd2)
        c = (mu2 * sd1**2 - sd2 * (mu1 * sd2 + sd1 * np.sqrt(disc))) / (sd1**2 - sd2**2)
    else:
        disc = (mu2 - mu1) ** 2 + 2.0 * (sd2**2 - sd1**2) * np.log(sd2 / sd1)
        c = (mu1 * sd2**2 - sd1 * (mu2 * sd1 + sd2 * np.sqrt(disc))) / (sd2**2 - sd1**2)
    if mu1 < mu2:
        a1 = 0.5 * (1.0 + erf((mu1 - c) / (sd1 * _SQRT2)))
        a2 = 0.5 * (1.0 + erf((c - mu2) / (sd2 * _SQRT2)))
    else:
        a1 = 0.5 * (1.0 + erf((c - mu1) / (sd1 * _SQRT2)))
        a2 = 0.5 * (1.0 + erf((mu2 - c) / (sd2 * _SQRT2)))
    return (a1 + a2) * 100.0


def test_percent_overlap_identical_and_far() -> None:
    assert percent_overlap(285.0, 0.2, 1.0, 285.0, 0.2, 1.0) == 100.0
    far = percent_overlap(280.0, 0.2, 1.0, 310.0, 0.2, 1.0)
    assert far == 0.0
    assert percent_overlap(285.0, 0.2, 1.0, 285.0, 0.2, 0.0) == 0.0


def test_percent_overlap_matches_independent_erf() -> None:
    mu1, sd1, mu2, sd2 = 285.0, 0.21, 285.4, 0.21
    got = percent_overlap(mu1, sd1, 1.2, mu2, sd2, 0.8)
    expect = _independent_unit_overlap(mu1, sd1, mu2, sd2)
    np.testing.assert_allclose(got, expect, rtol=0.0, atol=1e-10)
    mu1, sd1, mu2, sd2 = 286.0, 0.15, 287.2, 0.40
    got = percent_overlap(mu1, sd1, 2.0, mu2, sd2, 0.5)
    expect = _independent_unit_overlap(mu1, sd1, mu2, sd2)
    np.testing.assert_allclose(got, expect, rtol=0.0, atol=1e-8)


def test_overlap_matrix_symmetric_diagonal() -> None:
    energy = np.array([284.0, 284.3, 290.0])
    sigma = np.array([0.2, 0.2, 0.3])
    amp = np.array([1.0, 0.8, 0.4])
    ov = overlap_matrix(energy, sigma, amp)
    assert ov.shape == (3, 3)
    assert ov.dtype == np.float64
    np.testing.assert_allclose(np.diag(ov), 100.0)
    np.testing.assert_allclose(ov, ov.T)
    np.testing.assert_allclose(
        ov[0, 1], percent_overlap(284.0, 0.2, 1.0, 284.3, 0.2, 0.8)
    )


def test_simple_cluster3_adjacent_skip_isolated() -> None:
    amp = np.array([1.0, 1.0, 1.0])
    adjacent = np.array([[100.0, 80.0, 0.0], [80.0, 100.0, 0.0], [0.0, 0.0, 100.0]])
    groups = simple_cluster3(adjacent, amp, 50.0)
    assert groups == [[0, 1], [2]]
    skip = np.array([[100.0, 10.0, 80.0], [10.0, 100.0, 10.0], [80.0, 10.0, 100.0]])
    groups = simple_cluster3(skip, amp, 50.0)
    assert groups == [[0, 2], [1]]
    isolated = np.eye(3) * 100.0
    groups = simple_cluster3(isolated, amp, 50.0)
    assert groups == [[0], [1], [2]]


def test_effective_gaussian_single_and_coincident() -> None:
    grid = clustering_grid()
    mu = np.array([290.0])
    amp = np.array([1.5])
    sd = np.array([0.25])
    pos, area, sigma, xx, yy, zz, theta = effective_gaussian(
        mu,
        amp,
        sd,
        grid_ev=grid,
        os_xx=np.array([0.2]),
        os_yy=np.array([0.0]),
        os_zz=np.array([0.4]),
    )
    np.testing.assert_allclose(pos, 290.0, atol=0.03)
    np.testing.assert_allclose(sigma, 0.25, rtol=0.05)
    np.testing.assert_allclose(area, 1.5, rtol=0.05)
    np.testing.assert_allclose(xx, 0.2)
    np.testing.assert_allclose(yy, 0.0)
    np.testing.assert_allclose(zz, 0.4)
    mag = np.hypot(0.2, 0.4)
    np.testing.assert_allclose(theta, np.degrees(np.arccos(0.4 / mag)))
    pos2, area2, _sigma2, xx2, _yy2, _zz2, _th = effective_gaussian(
        np.array([290.0, 290.0]),
        np.array([1.0, 1.0]),
        np.array([0.25, 0.25]),
        grid_ev=grid,
        os_xx=np.array([0.1, 0.2]),
        os_yy=np.array([0.0, 0.0]),
        os_zz=np.array([0.0, 0.0]),
    )
    np.testing.assert_allclose(pos2, 290.0, atol=0.03)
    np.testing.assert_allclose(area2, 2.0, rtol=0.05)
    np.testing.assert_allclose(xx2, 0.3)


def test_cluster_by_overlap_terminates_below_threshold() -> None:
    sticks = TransitionSticks(
        energy_ev=np.array([284.0, 284.15, 295.0]),
        oscillator_strength=np.array([1.0, 0.9, 0.5]),
        sigma_ev=np.full(3, 0.2),
        site=np.array(["C1", "C1", "C2"]),
        os_xx=np.array([0.1, 0.1, 0.0]),
        os_yy=np.zeros(3),
        os_zz=np.array([0.2, 0.2, 0.5]),
    )
    result = cluster_by_overlap(sticks, 50.0)
    assert result.n_iterations >= 1
    assert result.n_iterations <= 20
    assert 1 <= result.energy_ev.shape[0] <= 3
    assert max_offdiag_overlap(result.final_overlap) < 50.0 or result.n_iterations == 20
    flat = [i for g in result.member_indices for i in g]
    assert sorted(flat) == [0, 1, 2]
    assert len(result.stages) == result.n_iterations + 1
    assert result.stages[0].label == "OS filtered"
    assert result.stages[0].n_peaks == 3
    np.testing.assert_allclose(result.stages[-1].overlap, result.final_overlap)
    assert result.stages[-1].n_peaks == result.energy_ev.shape[0]


def test_kneedle_recovers_synthetic_l_curve() -> None:
    x = np.linspace(0.0, 100.0, 201)
    y = np.where(x < 20.0, 100.0 - 0.2 * x, 96.0 - 0.9 * (x - 20.0))
    y = np.clip(y, 1.0, None)
    knee = kneedle_percent(x, y)
    assert 15.0 <= knee <= 30.0


def test_os_percent_elbow_on_bimodal_os() -> None:
    os_vals = np.concatenate([np.full(20, 1.0), np.full(5, 0.02)])
    energy = np.linspace(284.0, 290.0, os_vals.size)
    sticks = TransitionSticks(
        energy_ev=energy,
        oscillator_strength=os_vals,
        sigma_ev=np.full(os_vals.size, 0.2),
        site=np.array(["C1"] * os_vals.size),
        os_xx=np.zeros(os_vals.size),
        os_yy=np.zeros(os_vals.size),
        os_zz=np.zeros(os_vals.size),
    )
    pct, percents, n_kept, _frac = os_percent_elbow(
        sticks,
        energy_min_ev=280.0,
        energy_max_ev=320.0,
    )
    assert percents.shape == n_kept.shape
    assert 0.0 <= pct <= 100.0
    filtered, kept = filter_sticks(
        sticks,
        os_percent=pct,
        energy_min_ev=280.0,
        energy_max_ev=320.0,
    )
    assert kept.size == filtered.energy_ev.shape[0]
    assert filtered.energy_ev.shape[0] >= 1


def test_cluster_bic_rss_and_k_directions() -> None:
    n = 2000
    k = 10
    low = cluster_bic(1.0, n, k)
    high_rss = cluster_bic(4.0, n, k)
    high_k = cluster_bic(1.0, n, 20)
    assert low < high_rss
    assert low < high_k


def test_lhs_seed_reproducible() -> None:
    energy = np.linspace(284.0, 292.0, 8)
    sticks = TransitionSticks(
        energy_ev=energy,
        oscillator_strength=np.linspace(1.0, 0.4, 8),
        sigma_ev=np.full(8, 0.22),
        site=np.array(["C1"] * 8),
        os_xx=np.zeros(8),
        os_yy=np.zeros(8),
        os_zz=np.ones(8) * 0.1,
    )
    a, ra = select_overlap_threshold(sticks, n_samples=4, seed=7)
    b, rb = select_overlap_threshold(sticks, n_samples=4, seed=7)
    np.testing.assert_allclose(a.overlap_percent, b.overlap_percent)
    np.testing.assert_allclose(a.bic, b.bic)
    assert a.selected_overlap == b.selected_overlap
    np.testing.assert_allclose(ra.energy_ev, rb.energy_ev)
    c, _rc = select_overlap_threshold(sticks, n_samples=4, seed=8)
    assert c.lhs_seed == 8
    assert 0.0 in a.overlap_percent
    assert 90.0 in a.overlap_percent
    assert a.selected_overlap in set(a.overlap_percent)


def test_ovp_prior_and_trace_crossing() -> None:
    rng = np.random.default_rng(0)
    init = _initial_ovps(16, 0.0, 90.0, rng)
    assert 0.0 in init
    assert 50.0 in init
    assert 90.0 in init
    ladder = _sigma_ladder(0.0, 90.0)
    first_sigma = [x for x in ladder if abs(x - 50.0) <= DEFAULT_PRIOR_SIGMA]
    second_sigma = [
        x
        for x in ladder
        if DEFAULT_PRIOR_SIGMA < abs(x - 50.0) <= 2.0 * DEFAULT_PRIOR_SIGMA
    ]
    assert len(first_sigma) > len(second_sigma)
    ovp = np.array([0.0, 25.0, 50.0, 75.0, 90.0])
    n_cl = np.array([1, 2, 3, 4, 5])
    bic = np.array([5.0, 4.0, 3.0, 2.0, 1.0])
    assert _pick_trace_crossing(ovp, n_cl, bic) == 50.0


def test_ovp_gp_fills_sample_budget() -> None:
    energy = np.linspace(284.0, 292.0, 8)
    sticks = TransitionSticks(
        energy_ev=energy,
        oscillator_strength=np.linspace(1.0, 0.4, 8),
        sigma_ev=np.full(8, 0.22),
        site=np.array(["C1"] * 8),
        os_xx=np.zeros(8),
        os_yy=np.zeros(8),
        os_zz=np.ones(8) * 0.1,
    )
    search, _result = select_overlap_threshold(sticks, n_samples=8, seed=3)
    assert search.overlap_percent.shape[0] == 8
    assert len(set(search.overlap_percent.tolist())) == 8
    assert 0.0 in search.overlap_percent
    assert 90.0 in search.overlap_percent


def test_sticks_from_tp_table_prefers_aligned_energy() -> None:
    frame = pd.DataFrame(
        {
            "site": ["C2", "C1"],
            "energy_ev": [10.0, 11.0],
            "energy_aligned_ev": [284.2, 285.1],
            "oscillator_strength": [0.2, 0.3],
            "os_xx": [0.01, 0.02],
            "os_yy": [0.0, 0.0],
            "os_zz": [0.1, 0.1],
        }
    )
    sticks = sticks_from_tp_table(frame)
    np.testing.assert_allclose(sticks.energy_ev, [284.2, 285.1])
    assert sticks.sigma_ev.shape == (2,)
    assert np.all(sticks.sigma_ev > 0.0)


@pytest.mark.skipif(not _EXAMPLE.is_dir(), reason="ZnPc example run not present")
def test_znpc_cluster_report_smoke(tmp_path: Path) -> None:
    from dftlearn.visualization.xas_cluster_figure import write_xas_cluster_report

    out = write_xas_cluster_report(_EXAMPLE, tmp_path, n_samples=4)
    if out is None:
        pytest.skip("No TP stick files in ZnPc example")
    members, clusters, search, elbow_csv, summary = out
    for path in (members, clusters, search, elbow_csv, summary):
        assert path.is_file()
    search_df = pd.read_csv(search)
    cluster_df = pd.read_csv(clusters)
    member_df = pd.read_csv(members)
    assert np.isfinite(search_df["bic"]).all()
    n_clusters = len(cluster_df)
    n_members = len(member_df)
    assert 1 <= n_clusters <= n_members
    assert {"site", "cluster_id", "theta_deg"}.issubset(member_df.columns)
    assert int(member_df["cluster_id"].min()) >= 0
    assert int(member_df["cluster_id"].max()) == n_clusters - 1
