"""Energy-sorted overlap grouping (Igor ``simpleCluster3``).

Peaks must already be ordered by increasing energy. Grouping allows a
single skipped matrix neighbor (degree of adjacency 2) and re-anchors the
overlap row to the strongest oscillator strength in the growing cluster.
This module does not compute overlap or merge envelopes.
"""

from __future__ import annotations

import numpy as np


def simple_cluster3(
    overlap: np.ndarray,
    amplitude: np.ndarray,
    threshold: float,
) -> list[list[int]]:
    """Group energy-sorted peaks by overlap, skipping at most one neighbor.

    Ports ``simpleCluster3``: walk the overlap matrix from the low-energy
    corner, absorb a peak when ``overlap[i, j] >= threshold``, skip one
    failing neighbor, then continue from the highest-amplitude member of
    the current cluster against the last used peak.

    Parameters
    ----------
    overlap : numpy.ndarray
        Square percent-overlap matrix, shape ``(n, n)``, aligned with
        energy-sorted peaks.
    amplitude : numpy.ndarray
        Oscillator strengths, shape ``(n,)``. Used only to re-anchor the
        row index to the strongest member of the growing cluster.
    threshold : float
        Minimum overlap percent to join a peak to the current cluster.

    Returns
    -------
    list of list of int
        Each inner list is the peak indices (into ``overlap``) in one
        cluster, in the order they were absorbed. Every index appears
        once.

    Raises
    ------
    ValueError
        If ``overlap`` is not square, ``amplitude`` length disagrees, or
        ``n`` is 0.
    """
    ov = np.asarray(overlap, dtype=np.float64)
    amp = np.asarray(amplitude, dtype=np.float64)
    if ov.ndim != 2 or ov.shape[0] != ov.shape[1]:
        msg = "overlap must be a square 2-d array"
        raise ValueError(msg)
    n = ov.shape[0]
    if n == 0:
        msg = "simple_cluster3 requires at least one peak"
        raise ValueError(msg)
    if amp.shape != (n,):
        msg = "amplitude must have shape (n,) matching overlap"
        raise ValueError(msg)
    used = np.zeros(n, dtype=np.int8)
    clustered: list[list[int]] = []
    i = 0
    j = 0
    k = 0
    while True:
        _ensure_cluster(clustered, k)
        if i >= n or j >= n:
            i = _first_unused(used)
            j = i
            k += 1
            if i >= n or j >= n:
                break
            continue
        ovp = float(ov[i, j])
        if ovp >= threshold:
            if used[j] == 0:
                used[j] = 1
                clustered[k].append(j)
            j += 1
        else:
            j += 1
            if i >= n or j >= n:
                i = _first_unused(used)
                j = i
                k += 1
                if i >= n or j >= n:
                    break
                continue
            ovp = float(ov[i, j])
            if ovp >= threshold:
                if used[j] == 0:
                    used[j] = 1
                    clustered[k].append(j)
                j += 1
            else:
                i = _max_amplitude_index(amp, clustered[k])
                j = max(_last_used(used), 0)
                while True:
                    _ensure_cluster(clustered, k)
                    if i >= n or j >= n:
                        break
                    ovp = float(ov[i, j])
                    if ovp >= threshold:
                        if used[j] == 0:
                            used[j] = 1
                            clustered[k].append(j)
                        j += 1
                    else:
                        i = _first_unused(used)
                        j = i
                        k += 1
                        break
        if i >= n or j >= n:
            i = _first_unused(used)
            j = i
            k += 1
            if i >= n or j >= n:
                break
    groups = [g for g in clustered if g]
    leftover = [p for p in range(n) if p not in {x for g in groups for x in g}]
    groups.extend([[p] for p in leftover])
    _assert_partition(groups, n)
    return groups


def _ensure_cluster(clustered: list[list[int]], k: int) -> None:
    while len(clustered) <= k:
        clustered.append([])


def _first_unused(used: np.ndarray) -> int:
    unused = np.flatnonzero(used == 0)
    if unused.size == 0:
        return int(used.size)
    return int(unused[0])


def _last_used(used: np.ndarray) -> int:
    unused = np.flatnonzero(used == 0)
    if unused.size == 0:
        return int(used.size) - 1
    return int(unused[0]) - 1


def _max_amplitude_index(amplitude: np.ndarray, members: list[int]) -> int:
    if not members:
        return 0
    idx = np.asarray(members, dtype=np.int64)
    return int(idx[int(np.argmax(amplitude[idx]))])


def _assert_partition(groups: list[list[int]], n: int) -> None:
    seen: list[int] = []
    for group in groups:
        seen.extend(group)
    if sorted(seen) != list(range(n)):
        missing = [p for p in range(n) if p not in seen]
        extra = [p for p in seen if seen.count(p) > 1]
        msg = (
            f"simple_cluster3 did not partition 0..{n - 1}: "
            f"missing={missing} extra={extra}"
        )
        raise ValueError(msg)
