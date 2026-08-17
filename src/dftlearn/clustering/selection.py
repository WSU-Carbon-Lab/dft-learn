"""Overlap-percent search by a Normal(50%, 10%) ladder and a Gaussian process.

BIC follows Igor ``calculateBIC`` for the clustered vs OS-filtered DFT
envelope: ``n ln(RSS / n) + k ln(n)`` with ``k`` equal to the number of
effective clusters (amplitude-only). The search always evaluates 0% and
90%, which sit at -5 sigma and +4 sigma of ``Normal(50, 10)``. Remaining
points fill sigma shells from the mean outward, with spacing that grows
as 0.25 k sigma so the first sigma is sampled most densely. If the sample
budget exceeds that ladder, expected-improvement steps of a 1-d Gaussian
process on BIC fill the rest. The selected threshold is the evaluated
overlap where the min-max-scaled BIC and ``n_clusters`` traces are closest.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import numpy as np
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel

from dftlearn.clustering.iterate import cluster_by_overlap
from dftlearn.clustering.merge import clustering_grid, gaussian_envelope
from dftlearn.clustering.types import CLUSTERING_SPEC, OverlapSearchResult

if TYPE_CHECKING:
    from numpy.random import Generator

    from dftlearn.clustering.types import ClusterResult, TransitionSticks
    from dftlearn.io.stobe_xas_sticks import XasSpecSettings

DEFAULT_LHS_SAMPLES = 16
DEFAULT_OVP_MIN = 0.0
DEFAULT_OVP_MAX = 90.0
DEFAULT_LHS_SEED = 0
DEFAULT_PRIOR_MEAN = 50.0
DEFAULT_PRIOR_SIGMA = 10.0
_REFINE_STEP = 0.5
_ACQ_CROSS_WEIGHT = 0.5
_ACQ_CROSS_WIDTH = 10.0
_ANCHORS = (0.0, 90.0)


def cluster_bic(rss: float, n_points: int, n_clusters: int) -> float:
    """Evaluate Igor DFT-vs-cluster BIC.

    Computes ``n * ln(RSS / n) + k * ln(n)`` with ``k = n_clusters``.
    A non-positive RSS is replaced by the smallest positive ``float64``
    so the logarithm stays defined.

    Parameters
    ----------
    rss : float
        Residual sum of squares of the clustered envelope vs the
        OS-filtered DFT envelope. Must be finite.
    n_points : int
        Envelope grid length ``n``. Must be >= 2.
    n_clusters : int
        Number of effective Gaussians ``k``. Must be >= 1.

    Returns
    -------
    float
        BIC value (lower is better).

    Raises
    ------
    ValueError
        If ``rss`` is not finite, ``n_points < 2``, or ``n_clusters < 1``.
    """
    if not np.isfinite(rss):
        msg = "rss must be finite"
        raise ValueError(msg)
    if n_points < 2:
        msg = "n_points must be >= 2"
        raise ValueError(msg)
    if n_clusters < 1:
        msg = "n_clusters must be >= 1"
        raise ValueError(msg)
    rss_safe = max(float(rss), float(np.finfo(np.float64).tiny))
    n = float(n_points)
    k = float(n_clusters)
    return n * float(np.log(rss_safe / n)) + k * float(np.log(n))


def envelope_rss(
    sticks: TransitionSticks,
    result: ClusterResult,
    *,
    settings: XasSpecSettings = CLUSTERING_SPEC,
) -> float:
    """Sum of squared envelope residuals on the clustering grid.

    Parameters
    ----------
    sticks : TransitionSticks
        OS-filtered unclustered peaks (DFT envelope).
    result : ClusterResult
        Effective clusters.
    settings : XasSpecSettings, optional
        Envelope grid.

    Returns
    -------
    float
        ``sum((I_dft - I_clusters)^2)`` on the grid.
    """
    grid = clustering_grid(settings)
    dft = gaussian_envelope(
        grid,
        sticks.energy_ev,
        sticks.oscillator_strength,
        sticks.sigma_ev,
    )
    clustered = gaussian_envelope(
        grid,
        result.energy_ev,
        result.amplitude,
        result.sigma_ev,
    )
    delta = dft - clustered
    return float(np.sum(delta * delta))


def select_overlap_threshold(
    sticks: TransitionSticks,
    *,
    ovp_min: float = DEFAULT_OVP_MIN,
    ovp_max: float = DEFAULT_OVP_MAX,
    n_samples: int = DEFAULT_LHS_SAMPLES,
    seed: int = DEFAULT_LHS_SEED,
    settings: XasSpecSettings = CLUSTERING_SPEC,
) -> tuple[OverlapSearchResult, ClusterResult]:
    """Choose overlap percent by a Normal(50%, 10%) ladder and BIC GP.

    Evaluates 0% (-5 sigma) and 90% (+4 sigma) first, then fills
    ``Normal(50, 10)`` shells with spacing ``0.25 k sigma`` so the first
    sigma is densest. Remaining budget, if any, is spent on
    expected-improvement of a 1-d Gaussian process on BIC, with extra
    weight near the current BIC-vs-``n_clusters`` crossing. The returned
    threshold is the evaluated overlap at that crossing.

    Parameters
    ----------
    sticks : TransitionSticks
        OS- and energy-filtered sticks.
    ovp_min, ovp_max : float
        Inclusive overlap-percent bounds for the search.
    n_samples : int
        Total unique overlap percents to evaluate. Must be >= 2.
    seed : int
        Seed reserved for Gaussian-process tie-breaking when the ladder
        is shorter than ``n_samples``.
    settings : XasSpecSettings, optional
        Envelope grid and clustering widths already on ``sticks``.

    Returns
    -------
    search : OverlapSearchResult
        All evaluated overlap percents with RSS, BIC, and the selection.
    result : ClusterResult
        Clustering at the selected overlap percent.

    Raises
    ------
    ValueError
        If bounds are invalid, ``n_samples < 2``, or clustering fails.
    """
    if ovp_min < 0.0 or ovp_max > 100.0 or ovp_min >= ovp_max:
        msg = "require 0 <= ovp_min < ovp_max <= 100"
        raise ValueError(msg)
    if n_samples < 2:
        msg = "n_samples must be >= 2"
        raise ValueError(msg)
    ladder = _sigma_ladder(ovp_min, ovp_max)
    n_init = min(n_samples, len(ladder))
    rows: list[tuple[float, ClusterResult, float, float]] = []
    seen: set[float] = set()
    for ovp in ladder[:n_init]:
        if ovp in seen:
            continue
        seen.add(ovp)
        rows.append(_evaluate(sticks, ovp, settings=settings))
    while len(rows) < n_samples:
        ovp_arr = np.array([r[0] for r in rows], dtype=np.float64)
        n_arr = np.array([r[1].energy_ev.shape[0] for r in rows], dtype=np.int64)
        bic_arr = np.array([r[3] for r in rows], dtype=np.float64)
        nxt = _gp_next(
            ovp_arr,
            bic_arr,
            n_arr,
            ovp_min=ovp_min,
            ovp_max=ovp_max,
            seen=seen,
        )
        if nxt is None:
            break
        seen.add(nxt)
        rows.append(_evaluate(sticks, nxt, settings=settings))
    rows.sort(key=lambda r: r[0])
    ovp_out = np.array([r[0] for r in rows], dtype=np.float64)
    n_out = np.array([r[1].energy_ev.shape[0] for r in rows], dtype=np.int64)
    rss_out = np.array([r[2] for r in rows], dtype=np.float64)
    bic_out = np.array([r[3] for r in rows], dtype=np.float64)
    selected = _pick_trace_crossing(ovp_out, n_out, bic_out)
    selected_index = int(np.argmin(np.abs(ovp_out - selected)))
    search = OverlapSearchResult(
        overlap_percent=ovp_out,
        n_clusters=n_out,
        rss=rss_out,
        bic=bic_out,
        selected_overlap=float(ovp_out[selected_index]),
        selected_index=selected_index,
        lhs_seed=int(seed),
    )
    return search, rows[selected_index][1]


def _evaluate(
    sticks: TransitionSticks,
    overlap_percent: float,
    *,
    settings: XasSpecSettings,
) -> tuple[float, ClusterResult, float, float]:
    result = cluster_by_overlap(
        sticks,
        overlap_percent,
        settings=settings,
    )
    rss = envelope_rss(sticks, result, settings=settings)
    bic = cluster_bic(rss, settings.n_points, int(result.energy_ev.shape[0]))
    return overlap_percent, result, rss, bic


def _snap_ovp(value: float) -> float:
    return float(np.round(value / _REFINE_STEP) * _REFINE_STEP)


def _sigma_ladder(ovp_min: float, ovp_max: float) -> list[float]:
    mean = DEFAULT_PRIOR_MEAN
    sig = DEFAULT_PRIOR_SIGMA
    lo = _snap_ovp(ovp_min)
    hi = _snap_ovp(ovp_max)
    chosen: list[float] = []

    def add(raw: float) -> None:
        snapped = _snap_ovp(float(np.clip(raw, ovp_min, ovp_max)))
        if lo <= snapped <= hi and snapped not in chosen:
            chosen.append(snapped)

    for anchor in _ANCHORS:
        add(float(anchor))
    add(mean)
    for k in range(1, 6):
        step = _snap_ovp(0.25 * k * sig)
        if step < _REFINE_STEP:
            step = _REFINE_STEP
        max_off = k * sig
        min_off = (k - 1) * sig
        n_step = int(np.floor(max_off / step + 1e-9))
        for j in range(1, n_step + 1):
            off = j * step
            if k > 1 and off <= min_off + 1e-9:
                continue
            add(mean - off)
            add(mean + off)
    return chosen


def _initial_ovps(
    n_init: int,
    ovp_min: float,
    ovp_max: float,
    rng: Generator,
) -> list[float]:
    _ = rng
    return _sigma_ladder(ovp_min, ovp_max)[:n_init]


def _gp_next(
    ovp: np.ndarray,
    bic: np.ndarray,
    n_clusters: np.ndarray,
    *,
    ovp_min: float,
    ovp_max: float,
    seen: set[float],
) -> float | None:
    candidates = _candidate_grid(ovp_min, ovp_max, seen)
    if candidates.size == 0:
        return None
    if ovp.size < 3:
        return _maximin_next(candidates, ovp)
    try:
        gpr = _fit_bic_gp(ovp, bic, ovp_min, ovp_max)
        x_cand = _scale_ovp(candidates, ovp_min, ovp_max)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            mu, std = gpr.predict(x_cand, return_std=True)
    except (np.linalg.LinAlgError, ValueError):
        return _maximin_next(candidates, ovp)
    mu = np.asarray(mu, dtype=np.float64).reshape(-1)
    std = np.asarray(std, dtype=np.float64).reshape(-1)
    if not np.all(np.isfinite(mu)) or not np.all(np.isfinite(std)):
        return _maximin_next(candidates, ovp)
    ei = _expected_improvement(mu, std, float(np.min(bic)))
    x_cross = _pick_trace_crossing(ovp, n_clusters, bic)
    bump = std * np.exp(-0.5 * ((candidates - x_cross) / _ACQ_CROSS_WIDTH) ** 2)
    acq = ei + _ACQ_CROSS_WEIGHT * bump
    if not np.any(np.isfinite(acq)):
        return _maximin_next(candidates, ovp)
    return float(candidates[int(np.nanargmax(acq))])


def _candidate_grid(ovp_min: float, ovp_max: float, seen: set[float]) -> np.ndarray:
    raw = np.arange(
        ovp_min,
        ovp_max + 0.25 * _REFINE_STEP,
        _REFINE_STEP,
        dtype=np.float64,
    )
    snapped = np.array([_snap_ovp(float(x)) for x in raw], dtype=np.float64)
    unique = np.unique(snapped)
    keep = np.array([x not in seen for x in unique], dtype=bool)
    return unique[keep]


def _scale_ovp(ovp: np.ndarray, ovp_min: float, ovp_max: float) -> np.ndarray:
    span = max(ovp_max - ovp_min, 1e-9)
    return ((np.asarray(ovp, dtype=np.float64) - ovp_min) / span).reshape(-1, 1)


def _fit_bic_gp(
    ovp: np.ndarray,
    bic: np.ndarray,
    ovp_min: float,
    ovp_max: float,
) -> GaussianProcessRegressor:
    kernel = ConstantKernel(1.0) * RBF(length_scale=0.3) + WhiteKernel(noise_level=0.05)
    gpr = GaussianProcessRegressor(
        kernel=kernel,
        optimizer=None,
        normalize_y=True,
        alpha=1e-5,
    )
    gpr.fit(_scale_ovp(ovp, ovp_min, ovp_max), bic)
    return gpr


def _expected_improvement(
    mu: np.ndarray,
    sigma: np.ndarray,
    best: float,
) -> np.ndarray:
    sigma = np.maximum(np.asarray(sigma, dtype=np.float64), 1e-12)
    z = (best - mu) / sigma
    return (best - mu) * norm.cdf(z) + sigma * norm.pdf(z)


def _maximin_next(candidates: np.ndarray, ovp: np.ndarray) -> float:
    dist = np.min(np.abs(candidates[:, None] - ovp[None, :]), axis=1)
    return float(candidates[int(np.argmax(dist))])


def _minmax(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=np.float64)
    lo = float(np.min(x))
    hi = float(np.max(x))
    if hi - lo < 1e-12:
        return np.zeros(x.shape, dtype=np.float64)
    return (x - lo) / (hi - lo)


def _pick_trace_crossing(
    ovp: np.ndarray,
    n_clusters: np.ndarray,
    bic: np.ndarray,
) -> float:
    z_bic = _minmax(bic)
    z_n = _minmax(n_clusters)
    return float(ovp[int(np.argmin(np.abs(z_bic - z_n)))])
