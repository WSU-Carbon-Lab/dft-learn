"""Labeled stick tables and clustering result records.

This module owns the array layouts passed between overlap, merge, and
selection. It does not run clustering or I/O.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from dftlearn.clustering.overlap import FWHM_TO_SIGMA
from dftlearn.io.stobe_xas_sticks import XasSpecSettings
from dftlearn.xas.spectrum import piecewise_fwhm_ev

CLUSTERING_SPEC = XasSpecSettings(
    energy_min_ev=280.0,
    energy_max_ev=320.0,
    n_points=2000,
    fwhm_low_ev=0.5,
    fwhm_high_ev=12.0,
    e_break_low_ev=288.0,
    e_break_high_ev=320.0,
)


@dataclass(frozen=True, slots=True)
class TransitionSticks:
    """Energy-sorted dipole sticks used as clustering input.

    Attributes
    ----------
    energy_ev : numpy.ndarray
        Photon energies (eV), shape ``(n,)``, ``float64``.
    oscillator_strength : numpy.ndarray
        Dipole oscillator strengths, shape ``(n,)``, ``float64``.
    sigma_ev : numpy.ndarray
        Gaussian standard deviations (eV) from the clustering FWHM
        schedule, shape ``(n,)``.
    site : numpy.ndarray
        Site labels, shape ``(n,)``.
    os_xx, os_yy, os_zz : numpy.ndarray
        Cartesian oscillator-strength components, shape ``(n,)``.
    """

    energy_ev: np.ndarray
    oscillator_strength: np.ndarray
    sigma_ev: np.ndarray
    site: np.ndarray
    os_xx: np.ndarray
    os_yy: np.ndarray
    os_zz: np.ndarray

    def __post_init__(self) -> None:
        """Reject stick tables whose arrays do not share one length."""
        n = int(np.asarray(self.energy_ev).shape[0])
        names = (
            "oscillator_strength",
            "sigma_ev",
            "site",
            "os_xx",
            "os_yy",
            "os_zz",
        )
        for name in names:
            arr = np.asarray(getattr(self, name))
            if arr.shape[0] != n:
                msg = f"{name} length {arr.shape[0]} != energy length {n}"
                raise ValueError(msg)


@dataclass(frozen=True, slots=True)
class ClusterStage:
    """One overlap matrix along the merge sequence.

    The first stage is the OS-filtered stick table with first-pass
    widths. Each later stage is the overlap of effective Gaussians after
    one ``simpleCluster3`` merge pass.

    Attributes
    ----------
    label : str
        Stage name for figure titles (``OS filtered`` or ``Merge n``).
    overlap : numpy.ndarray
        Symmetric percent-overlap matrix, shape ``(n, n)``, ``float64``.
    n_peaks : int
        Number of peaks (matrix order) at this stage.
    max_offdiag : float
        Maximum off-diagonal overlap percent.
    """

    label: str
    overlap: np.ndarray
    n_peaks: int
    max_offdiag: float


@dataclass(frozen=True, slots=True)
class ClusterResult:
    """Effective Gaussians after iterative overlap merging.

    Attributes
    ----------
    energy_ev, amplitude, sigma_ev : numpy.ndarray
        Effective-cluster Gaussian parameters, shape ``(k,)``.
    os_xx, os_yy, os_zz : numpy.ndarray
        Summed Cartesian OS components, shape ``(k,)``.
    theta_deg : numpy.ndarray
        Polar angle of ``(os_xx, os_yy, os_zz)`` from +z, degrees.
    n_members : numpy.ndarray
        Stick counts per cluster, shape ``(k,)``, ``int64``.
    member_indices : tuple of tuple of int
        Original stick indices (into the filtered input table) in each
        cluster, energy-sorted cluster order.
    n_iterations : int
        Number of grouping passes performed (1 to 20).
    overlap_threshold : float
        Overlap percent used for merge decisions.
    final_overlap : numpy.ndarray
        Overlap matrix of the effective clusters, shape ``(k, k)``.
    stages : tuple of ClusterStage
        OS-filtered overlap followed by the overlap after each merge,
        including the terminal matrix (same data as ``final_overlap``).
    """

    energy_ev: np.ndarray
    amplitude: np.ndarray
    sigma_ev: np.ndarray
    os_xx: np.ndarray
    os_yy: np.ndarray
    os_zz: np.ndarray
    theta_deg: np.ndarray
    n_members: np.ndarray
    member_indices: tuple[tuple[int, ...], ...]
    n_iterations: int
    overlap_threshold: float
    final_overlap: np.ndarray
    stages: tuple[ClusterStage, ...]


@dataclass(frozen=True, slots=True)
class OverlapSearchResult:
    """Gaussian-process overlap-percent search with BIC and cluster counts.

    Attributes
    ----------
    overlap_percent : numpy.ndarray
        Evaluated overlap thresholds, shape ``(m,)``.
    n_clusters : numpy.ndarray
        Effective cluster counts, shape ``(m,)``, ``int64``.
    rss : numpy.ndarray
        Envelope residual sum of squares vs the OS-filtered DFT
        envelope, shape ``(m,)``.
    bic : numpy.ndarray
        Igor DFT-vs-cluster BIC values, shape ``(m,)``.
    selected_overlap : float
        Chosen overlap percent (BIC vs ``n_clusters`` trace crossing).
    selected_index : int
        Row of the selected evaluation.
    lhs_seed : int
        RNG seed used for the Normal prior and GP optimizer.
    """

    overlap_percent: np.ndarray
    n_clusters: np.ndarray
    rss: np.ndarray
    bic: np.ndarray
    selected_overlap: float
    selected_index: int
    lhs_seed: int


def sticks_from_tp_table(
    sticks: pd.DataFrame,
    *,
    settings: XasSpecSettings = CLUSTERING_SPEC,
) -> TransitionSticks:
    """Build clustering sticks from a ``collect_site_tp_xas`` stick table.

    Uses ``energy_aligned_ev`` when it is finite for a row, otherwise
    ``energy_ev``. Cartesian columns default to NaN when absent. Widths
    follow ``settings`` (Igor clustering FWHM), not a site ``xas.inp``.

    Parameters
    ----------
    sticks : pandas.DataFrame
        Rows with ``energy_ev``, ``oscillator_strength``, ``site``, and
        optionally ``energy_aligned_ev``, ``os_xx``, ``os_yy``, ``os_zz``.
    settings : XasSpecSettings, optional
        Clustering FWHM schedule.

    Returns
    -------
    TransitionSticks
        Rows sorted by energy with clustering sigmas assigned.

    Raises
    ------
    TypeError
        If ``sticks`` is not a pandas DataFrame.
    ValueError
        If required columns are missing or the table is empty.
    """
    if not isinstance(sticks, pd.DataFrame):
        msg = "sticks must be a pandas DataFrame"
        raise TypeError(msg)
    required = ("energy_ev", "oscillator_strength", "site")
    missing = [c for c in required if c not in sticks.columns]
    if missing:
        msg = f"Stick table missing columns {missing}"
        raise ValueError(msg)
    if sticks.empty:
        msg = "Stick table is empty"
        raise ValueError(msg)
    energy_raw = sticks["energy_ev"].to_numpy(dtype=np.float64)
    if "energy_aligned_ev" in sticks.columns:
        aligned = sticks["energy_aligned_ev"].to_numpy(dtype=np.float64)
        energy = np.where(np.isfinite(aligned), aligned, energy_raw)
    else:
        energy = energy_raw
    os_tot = sticks["oscillator_strength"].to_numpy(dtype=np.float64)
    site = sticks["site"].to_numpy()
    n = energy.shape[0]
    os_xx = _optional_os_column(sticks, "os_xx", n)
    os_yy = _optional_os_column(sticks, "os_yy", n)
    os_zz = _optional_os_column(sticks, "os_zz", n)
    fwhm = piecewise_fwhm_ev(energy, settings)
    sigma = fwhm * FWHM_TO_SIGMA
    order = np.argsort(energy, kind="mergesort")
    return TransitionSticks(
        energy_ev=energy[order],
        oscillator_strength=os_tot[order],
        sigma_ev=np.asarray(sigma, dtype=np.float64)[order],
        site=site[order],
        os_xx=os_xx[order],
        os_yy=os_yy[order],
        os_zz=os_zz[order],
    )


def _optional_os_column(frame: pd.DataFrame, name: str, n: int) -> np.ndarray:
    if name in frame.columns:
        return frame[name].to_numpy(dtype=np.float64)
    return np.full(n, np.nan, dtype=np.float64)
