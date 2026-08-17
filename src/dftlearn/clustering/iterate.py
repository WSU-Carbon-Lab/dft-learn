"""Iterative overlap clustering until the matrix is below threshold.

Applies Igor ``clusteringTransitions`` with ``simpleCluster3``: first
overlap uses a constant 0.4 eV FWHM; later passes use merged sigmas.
Stops when the maximum off-diagonal overlap is below the threshold or
after 20 grouping passes.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from dftlearn.clustering.grouping import simple_cluster3
from dftlearn.clustering.merge import clustering_grid, merge_groups
from dftlearn.clustering.overlap import (
    first_pass_sigma,
    max_offdiag_overlap,
    overlap_matrix,
)
from dftlearn.clustering.types import (
    CLUSTERING_SPEC,
    ClusterResult,
    ClusterStage,
    TransitionSticks,
)

if TYPE_CHECKING:
    from dftlearn.io.stobe_xas_sticks import XasSpecSettings

MAX_CLUSTER_ITERATIONS = 20


def cluster_by_overlap(
    sticks: TransitionSticks,
    overlap_threshold: float,
    *,
    settings: XasSpecSettings = CLUSTERING_SPEC,
    max_iterations: int = MAX_CLUSTER_ITERATIONS,
) -> ClusterResult:
    """Merge energy-sorted sticks until pairwise overlap is below threshold.

    Parameters
    ----------
    sticks : TransitionSticks
        Energy-sorted, already OS- and energy-filtered peaks.
    overlap_threshold : float
        Minimum percent overlap to join peaks (Igor ``ovpMax``).
    settings : XasSpecSettings, optional
        Envelope grid endpoints and point count.
    max_iterations : int, optional
        Maximum grouping passes (Igor breaks at 20).

    Returns
    -------
    ClusterResult
        Effective Gaussians, membership, merge-stage overlap matrices.

    Raises
    ------
    ValueError
        If there are no sticks, ``overlap_threshold`` is outside
        ``[0, 100]``, or grouping fails to partition the peaks.
    """
    n = sticks.energy_ev.shape[0]
    if n == 0:
        msg = "cluster_by_overlap requires at least one stick"
        raise ValueError(msg)
    if overlap_threshold < 0.0 or overlap_threshold > 100.0:
        msg = "overlap_threshold must be in [0, 100]"
        raise ValueError(msg)
    if max_iterations < 1:
        msg = "max_iterations must be >= 1"
        raise ValueError(msg)
    grid = clustering_grid(settings)
    energy = np.asarray(sticks.energy_ev, dtype=np.float64).copy()
    amp = np.asarray(sticks.oscillator_strength, dtype=np.float64).copy()
    sigma = np.asarray(sticks.sigma_ev, dtype=np.float64).copy()
    xx = np.asarray(sticks.os_xx, dtype=np.float64).copy()
    yy = np.asarray(sticks.os_yy, dtype=np.float64).copy()
    zz = np.asarray(sticks.os_zz, dtype=np.float64).copy()
    members: list[list[int]] = [[i] for i in range(n)]
    ov = overlap_matrix(energy, first_pass_sigma(n), amp)
    stages: list[ClusterStage] = [_overlap_stage("OS filtered", ov)]
    n_iter = 0
    theta = np.full(n, np.nan, dtype=np.float64)
    n_members = np.ones(n, dtype=np.int64)
    member_tuple: tuple[tuple[int, ...], ...] = tuple(tuple(m) for m in members)
    while n_iter < max_iterations:
        groups = simple_cluster3(ov, amp, overlap_threshold)
        n_iter += 1
        (
            energy,
            amp,
            sigma,
            xx,
            yy,
            zz,
            theta,
            n_members,
            member_tuple,
        ) = merge_groups(
            groups,
            energy,
            amp,
            sigma,
            os_xx=xx,
            os_yy=yy,
            os_zz=zz,
            member_indices=members,
            grid_ev=grid,
        )
        members = [list(m) for m in member_tuple]
        ov = overlap_matrix(energy, sigma, amp)
        stages.append(_overlap_stage(f"Merge {n_iter}", ov))
        if max_offdiag_overlap(ov) < overlap_threshold:
            break
    return ClusterResult(
        energy_ev=energy,
        amplitude=amp,
        sigma_ev=sigma,
        os_xx=xx,
        os_yy=yy,
        os_zz=zz,
        theta_deg=theta,
        n_members=n_members,
        member_indices=member_tuple,
        n_iterations=n_iter,
        overlap_threshold=float(overlap_threshold),
        final_overlap=ov,
        stages=tuple(stages),
    )


def _overlap_stage(label: str, overlap: np.ndarray) -> ClusterStage:
    ov = np.asarray(overlap, dtype=np.float64).copy()
    return ClusterStage(
        label=label,
        overlap=ov,
        n_peaks=int(ov.shape[0]),
        max_offdiag=float(max_offdiag_overlap(ov)),
    )


def filter_sticks(
    sticks: TransitionSticks,
    *,
    os_percent: float,
    energy_min_ev: float,
    energy_max_ev: float,
) -> tuple[TransitionSticks, np.ndarray]:
    """Drop sticks outside the energy window or below an OS% of the max.

    Oscillator-strength cutoff is ``(os_percent / 100) * max(OS)`` among
    sticks already inside the energy window (Igor ``tval``).

    Parameters
    ----------
    sticks : TransitionSticks
        Full stick table.
    os_percent : float
        Cutoff as a percent of the windowed maximum OS, in ``[0, 100]``.
    energy_min_ev, energy_max_ev : float
        Inclusive photon-energy window (eV).

    Returns
    -------
    filtered : TransitionSticks
        Surviving sticks in energy order.
    kept_index : numpy.ndarray
        Indices into ``sticks`` that were kept, ``int64``.

    Raises
    ------
    ValueError
        If ``os_percent`` is outside ``[0, 100]``, the window is empty, or
        no stick survives.
    """
    if os_percent < 0.0 or os_percent > 100.0:
        msg = "os_percent must be in [0, 100]"
        raise ValueError(msg)
    energy = np.asarray(sticks.energy_ev, dtype=np.float64)
    os_tot = np.asarray(sticks.oscillator_strength, dtype=np.float64)
    in_window = (energy >= energy_min_ev) & (energy <= energy_max_ev)
    if not np.any(in_window):
        msg = f"No sticks in energy window [{energy_min_ev}, {energy_max_ev}] eV"
        raise ValueError(msg)
    max_os = float(np.max(os_tot[in_window]))
    min_os = (os_percent / 100.0) * max_os
    keep = in_window & (os_tot >= min_os) & np.isfinite(os_tot) & (os_tot > 0.0)
    if not np.any(keep):
        msg = "OS filter removed every stick"
        raise ValueError(msg)
    idx = np.flatnonzero(keep)
    return (
        TransitionSticks(
            energy_ev=sticks.energy_ev[idx],
            oscillator_strength=sticks.oscillator_strength[idx],
            sigma_ev=sticks.sigma_ev[idx],
            site=sticks.site[idx],
            os_xx=sticks.os_xx[idx],
            os_yy=sticks.os_yy[idx],
            os_zz=sticks.os_zz[idx],
        ),
        idx.astype(np.int64, copy=False),
    )
