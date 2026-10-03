"""Percent overlap of Gaussian NEXAFS peaks (Igor ``peakOverlap``).

Pairwise overlap uses the unit-Gaussian error-function construction in
``compGaussErf3`` with ``normal=1``. Amplitudes enter only for identical-peak
and zero-amplitude guards and for the equal-center analytic branch. This
module does not group or merge peaks.
"""

from __future__ import annotations

import numpy as np
from scipy.special import erf

FIRST_OVERLAP_FWHM_EV = 0.4
FWHM_TO_SIGMA = 1.0 / 2.355
_SQRT2 = float(np.sqrt(2.0))
_SQRT_2PI = float(np.sqrt(2.0 * np.pi))
_EQ_ATOL = 1e-15


def percent_overlap(
    mu1: float,
    sd1: float,
    amp1: float,
    mu2: float,
    sd2: float,
    amp2: float,
) -> float:
    """Compute Igor percent overlap of two Gaussians.

    Matches ``peakOverlap`` / ``compGaussErf3(..., normal=1)``: identical
    peaks return 100; a zero amplitude returns 0; equal centers use the
    area ratio ``amp*sd*sqrt(2 pi)``; otherwise unit-Gaussian erf tails
    at the intersection, scaled so the result lies in ``[0, 100]``.

    Parameters
    ----------
    mu1, mu2 : float
        Peak centers (eV).
    sd1, sd2 : float
        Gaussian standard deviations (eV). Must be positive.
    amp1, amp2 : float
        Peak amplitudes (oscillator strength / area).

    Returns
    -------
    float
        Percent overlap in ``[0, 100]``.

    Raises
    ------
    ValueError
        If either standard deviation is not positive.
    """
    if sd1 <= 0.0 or sd2 <= 0.0:
        msg = "Gaussian standard deviations must be positive"
        raise ValueError(msg)
    if amp1 == 0.0 or amp2 == 0.0:
        if mu1 == mu2 and sd1 == sd2 and amp1 == amp2:
            return 100.0
        return 0.0
    if mu1 == mu2 and sd1 == sd2 and amp1 == amp2:
        return 100.0
    if mu1 == mu2:
        area1 = amp1 * sd1 * _SQRT_2PI
        area2 = amp2 * sd2 * _SQRT_2PI
        tot = min(area1, area2)
        if tot <= 0.0:
            return 0.0
        return _clamp_percent((area1 / tot) * 100.0)
    return _clamp_percent(_unit_erf_overlap(mu1, sd1, mu2, sd2))


def overlap_matrix(
    energy_ev: np.ndarray,
    sigma_ev: np.ndarray,
    amplitude: np.ndarray,
) -> np.ndarray:
    """Fill a symmetric percent-overlap matrix for all peak pairs.

    Diagonal entries are 100. The implementation is the broadcast form of
    :func:`percent_overlap` (Igor ``normal=1``).

    Parameters
    ----------
    energy_ev, sigma_ev, amplitude : numpy.ndarray
        Peak parameters, each shape ``(n,)``, ``float64``. ``sigma_ev``
        must be strictly positive.

    Returns
    -------
    numpy.ndarray
        Shape ``(n, n)``, ``float64`` percent overlaps.

    Raises
    ------
    ValueError
        If shapes disagree, ``n`` is 0, or a sigma is not positive.
    """
    mu = np.asarray(energy_ev, dtype=np.float64)
    sd = np.asarray(sigma_ev, dtype=np.float64)
    amp = np.asarray(amplitude, dtype=np.float64)
    if mu.ndim != 1 or sd.shape != mu.shape or amp.shape != mu.shape:
        msg = "energy_ev, sigma_ev, and amplitude must be 1-d arrays of equal length"
        raise ValueError(msg)
    if mu.size == 0:
        msg = "overlap_matrix requires at least one peak"
        raise ValueError(msg)
    if np.any(sd <= 0.0):
        msg = "Gaussian standard deviations must be positive"
        raise ValueError(msg)
    mu_i = mu[:, np.newaxis]
    mu_j = mu[np.newaxis, :]
    sd_i = sd[:, np.newaxis]
    sd_j = sd[np.newaxis, :]
    amp_i = amp[:, np.newaxis]
    amp_j = amp[np.newaxis, :]
    identical = (mu_i == mu_j) & (sd_i == sd_j) & (amp_i == amp_j)
    zero = (amp_i == 0.0) | (amp_j == 0.0)
    same_mu = mu_i == mu_j
    area_i = amp_i * sd_i * _SQRT_2PI
    area_j = amp_j * sd_j * _SQRT_2PI
    tot = np.minimum(area_i, area_j)
    with np.errstate(divide="ignore", invalid="ignore"):
        same_mu_pct = np.where(tot > 0.0, (area_i / tot) * 100.0, 0.0)
    erf_pct = _unit_erf_overlap_broadcast(mu_i, sd_i, mu_j, sd_j)
    out = np.where(
        identical,
        100.0,
        np.where(
            zero,
            0.0,
            np.where(same_mu, same_mu_pct, erf_pct),
        ),
    )
    np.fill_diagonal(out, 100.0)
    return _clamp_percent_array(out)


def first_pass_sigma(n: int) -> np.ndarray:
    """Return the constant first-pass sigma (0.4 eV FWHM) for ``n`` peaks.

    Parameters
    ----------
    n : int
        Peak count. Must be >= 1.

    Returns
    -------
    numpy.ndarray
        Shape ``(n,)``, ``float64``.

    Raises
    ------
    ValueError
        If ``n`` is less than 1.
    """
    if n < 1:
        msg = "n must be >= 1"
        raise ValueError(msg)
    return np.full(n, FIRST_OVERLAP_FWHM_EV * FWHM_TO_SIGMA, dtype=np.float64)


def max_offdiag_overlap(overlap: np.ndarray) -> float:
    """Return the maximum off-diagonal overlap percent.

    Parameters
    ----------
    overlap : numpy.ndarray
        Square overlap matrix.

    Returns
    -------
    float
        Maximum off-diagonal value, or 0 if ``n`` is 1.

    Raises
    ------
    ValueError
        If ``overlap`` is not a square 2-d array.
    """
    ov = np.asarray(overlap, dtype=np.float64)
    if ov.ndim != 2 or ov.shape[0] != ov.shape[1]:
        msg = "overlap must be a square 2-d array"
        raise ValueError(msg)
    n = ov.shape[0]
    if n <= 1:
        return 0.0
    mask = ~np.eye(n, dtype=bool)
    return float(np.max(ov[mask]))


def _unit_erf_overlap(mu1: float, sd1: float, mu2: float, sd2: float) -> float:
    if abs(sd1 - sd2) <= _EQ_ATOL:
        c = 0.5 * (mu1 + mu2)
    elif mu1 < mu2:
        disc = (mu1 - mu2) ** 2 + 2.0 * (sd1**2 - sd2**2) * np.log(sd1 / sd2)
        if disc < 0.0:
            return 0.0
        c = (mu2 * sd1**2 - sd2 * (mu1 * sd2 + sd1 * np.sqrt(disc))) / (sd1**2 - sd2**2)
    else:
        disc = (mu2 - mu1) ** 2 + 2.0 * (sd2**2 - sd1**2) * np.log(sd2 / sd1)
        if disc < 0.0:
            return 0.0
        c = (mu1 * sd2**2 - sd1 * (mu2 * sd1 + sd2 * np.sqrt(disc))) / (sd2**2 - sd1**2)
    if mu1 < mu2:
        a1 = 0.5 * (1.0 + erf((mu1 - c) / (sd1 * _SQRT2)))
        a2 = 0.5 * (1.0 + erf((c - mu2) / (sd2 * _SQRT2)))
    else:
        a1 = 0.5 * (1.0 + erf((c - mu1) / (sd1 * _SQRT2)))
        a2 = 0.5 * (1.0 + erf((mu2 - c) / (sd2 * _SQRT2)))
    return (a1 + a2) * 100.0


def _unit_erf_overlap_broadcast(
    mu_i: np.ndarray,
    sd_i: np.ndarray,
    mu_j: np.ndarray,
    sd_j: np.ndarray,
) -> np.ndarray:
    equal_sd = np.abs(sd_i - sd_j) <= _EQ_ATOL
    c_eq = 0.5 * (mu_i + mu_j)
    lt = mu_i < mu_j
    with np.errstate(divide="ignore", invalid="ignore"):
        disc_lt = (mu_i - mu_j) ** 2 + 2.0 * (sd_i**2 - sd_j**2) * np.log(sd_i / sd_j)
        disc_gt = (mu_j - mu_i) ** 2 + 2.0 * (sd_j**2 - sd_i**2) * np.log(sd_j / sd_i)
        sqrt_lt = np.sqrt(np.maximum(disc_lt, 0.0))
        sqrt_gt = np.sqrt(np.maximum(disc_gt, 0.0))
        c_lt = (mu_j * sd_i**2 - sd_j * (mu_i * sd_j + sd_i * sqrt_lt)) / (
            sd_i**2 - sd_j**2
        )
        c_gt = (mu_i * sd_j**2 - sd_i * (mu_j * sd_i + sd_j * sqrt_gt)) / (
            sd_j**2 - sd_i**2
        )
    c = np.where(equal_sd, c_eq, np.where(lt, c_lt, c_gt))
    a1_lt = 0.5 * (1.0 + erf((mu_i - c) / (sd_i * _SQRT2)))
    a2_lt = 0.5 * (1.0 + erf((c - mu_j) / (sd_j * _SQRT2)))
    a1_gt = 0.5 * (1.0 + erf((c - mu_i) / (sd_i * _SQRT2)))
    a2_gt = 0.5 * (1.0 + erf((mu_j - c) / (sd_j * _SQRT2)))
    ovp = np.where(lt, a1_lt + a2_lt, a1_gt + a2_gt) * 100.0
    bad = (~equal_sd) & np.where(lt, disc_lt < 0.0, disc_gt < 0.0)
    return np.where(np.isfinite(ovp) & ~bad, ovp, 0.0)


def _clamp_percent(value: float) -> float:
    if value < 0.01:
        return 0.0
    if value > 100.0:
        return 100.0
    return float(value)


def _clamp_percent_array(values: np.ndarray) -> np.ndarray:
    out = np.asarray(values, dtype=np.float64)
    out = np.where(out < 0.01, 0.0, out)
    return np.where(out > 100.0, 100.0, out)
