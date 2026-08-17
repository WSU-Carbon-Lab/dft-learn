"""Oscillator-strength cutoff from an elbow (knee) of ``n_kept`` vs OS%.

Scans the Igor ``tval`` axis (percent of the energy-windowed maximum OS)
and picks the knee of the decreasing stick-count curve by maximum
perpendicular distance to the endpoint chord. This module does not cluster.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from dftlearn.clustering.types import TransitionSticks


def os_keep_curve(
    sticks: TransitionSticks,
    *,
    energy_min_ev: float,
    energy_max_ev: float,
    os_percents: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate surviving stick counts along an OS-percent grid.

    Parameters
    ----------
    sticks : TransitionSticks
        Candidate transitions.
    energy_min_ev, energy_max_ev : float
        Inclusive energy window (eV) applied before the OS cutoff.
    os_percents : numpy.ndarray, optional
        Percents in ``[0, 100]``. Default is 0.5% steps from 0 to 100.

    Returns
    -------
    os_percents : numpy.ndarray
        Grid, shape ``(m,)``, ``float64``.
    n_kept : numpy.ndarray
        Surviving counts, shape ``(m,)``, ``int64``.
    retained_os_fraction : numpy.ndarray
        Fraction of windowed total OS retained, shape ``(m,)``.

    Raises
    ------
    ValueError
        If no sticks lie in the energy window or ``os_percents`` is empty.
    """
    energy = np.asarray(sticks.energy_ev, dtype=np.float64)
    os_tot = np.asarray(sticks.oscillator_strength, dtype=np.float64)
    in_window = (energy >= energy_min_ev) & (energy <= energy_max_ev)
    if not np.any(in_window):
        msg = f"No sticks in energy window [{energy_min_ev}, {energy_max_ev}] eV"
        raise ValueError(msg)
    window_os = os_tot[in_window]
    max_os = float(np.max(window_os))
    total_os = float(np.sum(window_os[np.isfinite(window_os) & (window_os > 0.0)]))
    if os_percents is None:
        percents = np.linspace(0.0, 100.0, 201, dtype=np.float64)
    else:
        percents = np.asarray(os_percents, dtype=np.float64)
        if percents.size == 0:
            msg = "os_percents must be non-empty"
            raise ValueError(msg)
    n_kept = np.empty(percents.shape[0], dtype=np.int64)
    frac = np.empty(percents.shape[0], dtype=np.float64)
    for i, pct in enumerate(percents):
        min_os = (float(pct) / 100.0) * max_os
        keep = in_window & np.isfinite(os_tot) & (os_tot > 0.0) & (os_tot >= min_os)
        n_kept[i] = int(np.count_nonzero(keep))
        kept_os = float(np.sum(os_tot[keep])) if total_os > 0.0 else 0.0
        frac[i] = kept_os / total_os if total_os > 0.0 else 0.0
    return percents, n_kept, frac


def kneedle_percent(x: np.ndarray, y: np.ndarray) -> float:
    """Return the ``x`` at maximum chord deviation of a normalized L-curve.

    ``x`` must be strictly increasing. Both axes are scaled to ``[0, 1]``
    and the knee is the interior sample with the largest absolute
    vertical deviation from the endpoint chord. If the curve is flat,
    the first ``x`` is returned.

    Parameters
    ----------
    x, y : numpy.ndarray
        Curve samples, each shape ``(m,)`` with ``m >= 2``.

    Returns
    -------
    float
        ``x`` coordinate of the knee.

    Raises
    ------
    ValueError
        If lengths disagree, ``m < 2``, or ``x`` is not strictly increasing.
    """
    xx = np.asarray(x, dtype=np.float64)
    yy = np.asarray(y, dtype=np.float64)
    if xx.shape != yy.shape or xx.ndim != 1:
        msg = "x and y must be 1-d arrays of equal length"
        raise ValueError(msg)
    if xx.size < 2:
        msg = "kneedle_percent requires at least two samples"
        raise ValueError(msg)
    if np.any(np.diff(xx) <= 0.0):
        msg = "x must be strictly increasing"
        raise ValueError(msg)
    x0 = xx[0]
    x1 = xx[-1]
    y_min = float(np.min(yy))
    y_max = float(np.max(yy))
    dx = x1 - x0
    dy = y_max - y_min
    if dy == 0.0:
        return float(x0)
    xn = (xx - x0) / dx
    yn = (yy - y_min) / dy
    y_chord = yn[0] + xn * (yn[-1] - yn[0])
    delta = yn - y_chord
    if xx.size == 2:
        return float(x0)
    interior = np.abs(delta[1:-1])
    return float(xx[int(np.argmax(interior)) + 1])


def os_percent_elbow(
    sticks: TransitionSticks,
    *,
    energy_min_ev: float,
    energy_max_ev: float,
    os_percents: np.ndarray | None = None,
) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    """Choose an OS% cutoff from the knee of ``n_kept`` vs OS%.

    The knee is clamped so at least one windowed stick survives.

    Parameters
    ----------
    sticks : TransitionSticks
        Candidate transitions.
    energy_min_ev, energy_max_ev : float
        Inclusive energy window (eV).
    os_percents : numpy.ndarray, optional
        OS% grid. Default is 0.5% steps from 0 to 100.

    Returns
    -------
    os_percent : float
        Selected cutoff.
    percents, n_kept, retained_os_fraction : numpy.ndarray
        The scanned curve.

    Raises
    ------
    ValueError
        If the window is empty or every cutoff drops all sticks.
    """
    percents, n_kept, frac = os_keep_curve(
        sticks,
        energy_min_ev=energy_min_ev,
        energy_max_ev=energy_max_ev,
        os_percents=os_percents,
    )
    if np.all(n_kept < 1):
        msg = "OS scan kept zero sticks at every percent"
        raise ValueError(msg)
    viable = n_kept >= 1
    x = percents[viable]
    y = n_kept[viable].astype(np.float64)
    if x.size == 1:
        return float(x[0]), percents, n_kept, frac
    knee = kneedle_percent(x, y)
    # Snap to the scanned grid; never choose a cutoff that drops all sticks.
    idx = int(np.argmin(np.abs(percents - knee)))
    while idx < percents.size and n_kept[idx] < 1:
        idx -= 1
    if idx < 0 or n_kept[idx] < 1:
        idx = int(np.max(np.flatnonzero(n_kept >= 1)))
    return float(percents[idx]), percents, n_kept, frac
