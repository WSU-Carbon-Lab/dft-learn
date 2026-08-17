"""Envelope effective Gaussians (Igor ``makeMergedPW``).

A cluster is replaced by one Gaussian whose height and FWHM match the
summed member envelope on a uniform energy grid. Cartesian oscillator
strengths are summed. This module does not compute overlap or choose
groups.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from dftlearn.clustering.overlap import FWHM_TO_SIGMA
from dftlearn.clustering.types import CLUSTERING_SPEC

if TYPE_CHECKING:
    from dftlearn.io.stobe_xas_sticks import XasSpecSettings

_SQRT_2PI = float(np.sqrt(2.0 * np.pi))
type _MergedArrays = tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    tuple[tuple[int, ...], ...],
]
type _MergeRec = tuple[
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    int,
    tuple[int, ...],
]


def clustering_grid(settings: XasSpecSettings = CLUSTERING_SPEC) -> np.ndarray:
    """Build the inclusive envelope grid used for effective-cluster fits.

    Parameters
    ----------
    settings : XasSpecSettings, optional
        ``energy_min_ev``, ``energy_max_ev``, and ``n_points``.

    Returns
    -------
    numpy.ndarray
        Shape ``(n_points,)``, ``float64`` energies in eV.

    Raises
    ------
    ValueError
        If ``n_points`` is less than 2.
    """
    if settings.n_points < 2:
        msg = "n_points must be >= 2"
        raise ValueError(msg)
    return np.linspace(
        settings.energy_min_ev,
        settings.energy_max_ev,
        settings.n_points,
        dtype=np.float64,
    )


def gaussian_envelope(
    grid_ev: np.ndarray,
    energy_ev: np.ndarray,
    amplitude: np.ndarray,
    sigma_ev: np.ndarray,
) -> np.ndarray:
    """Sum PDF-normalized Gaussians onto ``grid_ev``.

    Intensity is ``sum_i amp_i * N(E; mu_i, sigma_i)``.

    Parameters
    ----------
    grid_ev : numpy.ndarray
        Evaluation energies (eV), shape ``(m,)``.
    energy_ev, amplitude, sigma_ev : numpy.ndarray
        Peak parameters, each shape ``(k,)``. ``sigma_ev`` must be
        positive.

    Returns
    -------
    numpy.ndarray
        Envelope intensity, shape ``(m,)``, ``float64``.

    Raises
    ------
    ValueError
        If shapes disagree or a sigma is not positive.
    """
    grid = np.asarray(grid_ev, dtype=np.float64)
    mu = np.asarray(energy_ev, dtype=np.float64)
    amp = np.asarray(amplitude, dtype=np.float64)
    sd = np.asarray(sigma_ev, dtype=np.float64)
    if grid.ndim != 1 or mu.ndim != 1:
        msg = "grid_ev and energy_ev must be 1-d"
        raise ValueError(msg)
    if amp.shape != mu.shape or sd.shape != mu.shape:
        msg = "energy_ev, amplitude, and sigma_ev must have the same shape"
        raise ValueError(msg)
    if mu.size == 0:
        return np.zeros(grid.shape[0], dtype=np.float64)
    if np.any(sd <= 0.0):
        msg = "Gaussian standard deviations must be positive"
        raise ValueError(msg)
    z = (grid[:, np.newaxis] - mu[np.newaxis, :]) / sd[np.newaxis, :]
    pdf = np.exp(-0.5 * z * z) / (sd[np.newaxis, :] * _SQRT_2PI)
    return (amp[np.newaxis, :] * pdf).sum(axis=1)


def effective_gaussian(
    energy_ev: np.ndarray,
    amplitude: np.ndarray,
    sigma_ev: np.ndarray,
    *,
    grid_ev: np.ndarray | None = None,
    os_xx: np.ndarray | None = None,
    os_yy: np.ndarray | None = None,
    os_zz: np.ndarray | None = None,
) -> tuple[float, float, float, float, float, float, float]:
    """Fit one envelope Gaussian to the summed member peaks.

    Position is the grid argmax of the envelope. Sigma is FWHM / 2.355
    from the first half-max crossings on each side of the peak (Igor
    ``FindLevel``). Amplitude is ``sqrt(2 pi) * sigma * max_height``.
    Tensor components are sums; theta is the polar angle of the summed
    ``(xx, yy, zz)`` from +z.

    Parameters
    ----------
    energy_ev, amplitude, sigma_ev : numpy.ndarray
        Member peak parameters, each shape ``(p,)``, ``p >= 1``.
    grid_ev : numpy.ndarray, optional
        Envelope grid (eV). Default is :func:`clustering_grid`.
    os_xx, os_yy, os_zz : numpy.ndarray, optional
        Cartesian OS components, each shape ``(p,)``. Missing arrays
        contribute 0.

    Returns
    -------
    energy : float
        Envelope peak position (eV).
    amplitude : float
        Effective Gaussian area.
    sigma : float
        Effective standard deviation (eV).
    os_xx, os_yy, os_zz : float
        Summed Cartesian components.
    theta_deg : float
        Polar angle in degrees, or NaN if the tensor magnitude is 0.

    Raises
    ------
    ValueError
        If there are no members or shapes disagree.
    """
    mu = np.asarray(energy_ev, dtype=np.float64)
    amp = np.asarray(amplitude, dtype=np.float64)
    sd = np.asarray(sigma_ev, dtype=np.float64)
    if mu.size == 0:
        msg = "effective_gaussian requires at least one member"
        raise ValueError(msg)
    if grid_ev is None:
        grid = clustering_grid()
    else:
        grid = np.asarray(grid_ev, dtype=np.float64)
    envelope = gaussian_envelope(grid, mu, amp, sd)
    peak_idx = int(np.argmax(envelope))
    max_height = float(envelope[peak_idx])
    pos = float(grid[peak_idx])
    half = 0.5 * max_height
    en1 = _first_level_crossing(grid[: peak_idx + 1], envelope[: peak_idx + 1], half)
    en2 = _first_level_crossing(grid[peak_idx:], envelope[peak_idx:], half)
    if en1 is None:
        en1 = float(grid[0])
    if en2 is None:
        en2 = float(grid[-1])
    fwhm = max(en2 - en1, 0.0)
    sigma = fwhm * FWHM_TO_SIGMA
    if sigma <= 0.0:
        sigma = float(np.max(sd))
    area = _SQRT_2PI * sigma * max_height
    xx = _sum_or_zero(os_xx, mu.size)
    yy = _sum_or_zero(os_yy, mu.size)
    zz = _sum_or_zero(os_zz, mu.size)
    mag = float(np.sqrt(xx * xx + yy * yy + zz * zz))
    if mag == 0.0:
        theta = float("nan")
    else:
        theta = float(np.degrees(np.arccos(np.clip(zz / mag, -1.0, 1.0))))
    return pos, area, sigma, xx, yy, zz, theta


def merge_groups(
    groups: list[list[int]],
    energy_ev: np.ndarray,
    amplitude: np.ndarray,
    sigma_ev: np.ndarray,
    *,
    os_xx: np.ndarray | None = None,
    os_yy: np.ndarray | None = None,
    os_zz: np.ndarray | None = None,
    member_indices: list[list[int]] | None = None,
    grid_ev: np.ndarray | None = None,
) -> _MergedArrays:
    """Replace each index group with one effective Gaussian.

    Parameters
    ----------
    groups : list of list of int
        Peak indices into the current parameter arrays.
    energy_ev, amplitude, sigma_ev : numpy.ndarray
        Current peak parameters, shape ``(n,)``.
    os_xx, os_yy, os_zz : numpy.ndarray, optional
        Cartesian components, shape ``(n,)``.
    member_indices : list of list of int, optional
        Original stick indices for each current peak. Default is
        ``[[i] for i in range(n)]``.
    grid_ev : numpy.ndarray, optional
        Envelope grid passed to :func:`effective_gaussian`.

    Returns
    -------
    energy, amplitude, sigma, os_xx, os_yy, os_zz, theta, n_members
        Effective parameters, each shape ``(k,)``, energy-sorted.
    members : tuple of tuple of int
        Original stick indices per effective cluster.

    Raises
    ------
    ValueError
        If a group is empty or an index is out of range.
    """
    energy = np.asarray(energy_ev, dtype=np.float64)
    amp = np.asarray(amplitude, dtype=np.float64)
    sd = np.asarray(sigma_ev, dtype=np.float64)
    n = energy.shape[0]
    xx_all = _as_or_zeros(os_xx, n)
    yy_all = _as_or_zeros(os_yy, n)
    zz_all = _as_or_zeros(os_zz, n)
    lineage = [[i] for i in range(n)] if member_indices is None else member_indices
    recs: list[_MergeRec] = []
    for group in groups:
        if not group:
            msg = "merge_groups received an empty group"
            raise ValueError(msg)
        idx = np.asarray(group, dtype=np.int64)
        if np.any(idx < 0) or np.any(idx >= n):
            msg = "group index out of range"
            raise ValueError(msg)
        pos, area, sigma, xx, yy, zz, theta = effective_gaussian(
            energy[idx],
            amp[idx],
            sd[idx],
            grid_ev=grid_ev,
            os_xx=xx_all[idx],
            os_yy=yy_all[idx],
            os_zz=zz_all[idx],
        )
        orig: list[int] = []
        for g in idx.tolist():
            orig.extend(lineage[int(g)])
        recs.append(
            (
                pos,
                area,
                sigma,
                xx,
                yy,
                zz,
                theta,
                int(idx.size),
                tuple(orig),
            )
        )
    recs.sort(key=lambda r: r[0])
    k = len(recs)
    out_e = np.empty(k, dtype=np.float64)
    out_a = np.empty(k, dtype=np.float64)
    out_s = np.empty(k, dtype=np.float64)
    out_x = np.empty(k, dtype=np.float64)
    out_y = np.empty(k, dtype=np.float64)
    out_z = np.empty(k, dtype=np.float64)
    out_t = np.empty(k, dtype=np.float64)
    out_n = np.empty(k, dtype=np.int64)
    members: list[tuple[int, ...]] = []
    for i, rec in enumerate(recs):
        out_e[i] = rec[0]
        out_a[i] = rec[1]
        out_s[i] = rec[2]
        out_x[i] = rec[3]
        out_y[i] = rec[4]
        out_z[i] = rec[5]
        out_t[i] = rec[6]
        out_n[i] = rec[7]
        members.append(rec[8])
    return out_e, out_a, out_s, out_x, out_y, out_z, out_t, out_n, tuple(members)


def _as_or_zeros(values: np.ndarray | None, n: int) -> np.ndarray:
    if values is None:
        return np.zeros(n, dtype=np.float64)
    return np.asarray(values, dtype=np.float64)


def _sum_or_zero(values: np.ndarray | None, n: int) -> float:
    if values is None:
        return 0.0
    arr = np.asarray(values, dtype=np.float64)
    if arr.shape[0] != n:
        msg = "Cartesian OS arrays must match the member count"
        raise ValueError(msg)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return 0.0
    return float(finite.sum())


def _first_level_crossing(
    x: np.ndarray,
    y: np.ndarray,
    level: float,
) -> float | None:
    if x.size == 0:
        return None
    if y[0] == level:
        return float(x[0])
    for i in range(x.size - 1):
        y0 = float(y[i])
        y1 = float(y[i + 1])
        if (y0 - level) * (y1 - level) <= 0.0 and y0 != y1:
            t = (level - y0) / (y1 - y0)
            return float(x[i] + t * (x[i + 1] - x[i]))
    return None
