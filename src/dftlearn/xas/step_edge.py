"""Henke ``.nff`` tables and IP-weighted Gaussian step edges.

Pure helpers for bare-atom continuum absorption and DFT ionization-potential
step construction. Experiment anchoring and RMS baseline fits remain in
``python_pipeline.stepEdgeProcs``.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import periodictable
from scipy.interpolate import interp1d
from scipy.special import erf

if TYPE_CHECKING:
    from numpy.typing import NDArray

_NA = 6.0221415e23
_RE_CM = 2.81794e-13
_EV_TO_WAVELENGTH_CM = 1.23984e-4


def parse_henke_nff(path: Path) -> pd.DataFrame:
    """Parse a Henke ``.nff`` scattering-factor table.

    Skips blank lines and ``#`` comments. Sentinel values ``-9999`` become NaN.

    Parameters
    ----------
    path : pathlib.Path
        Path to an element ``.nff`` file (for example ``c.nff``).

    Returns
    -------
    pandas.DataFrame
        Columns ``energy_ev``, ``f1``, ``f2`` as ``float64``.

    Raises
    ------
    FileNotFoundError
        When ``path`` is missing.
    ValueError
        When no numeric rows are parsed.
    """
    path = Path(path)
    if not path.is_file():
        msg = f"Henke .nff file not found: {path}"
        raise FileNotFoundError(msg)
    rows: list[tuple[float, float, float]] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                energy = float(parts[0])
                f1 = float("nan") if parts[1].startswith("-9999") else float(parts[1])
                f2 = float("nan") if parts[2].startswith("-9999") else float(parts[2])
            except ValueError:
                continue
            rows.append((energy, f1, f2))
    if not rows:
        msg = f"No valid Henke rows in {path}"
        raise ValueError(msg)
    return pd.DataFrame.from_records(rows, columns=("energy_ev", "f1", "f2"))


def _atomic_weight(element_symbol: str) -> float:
    """Return elemental mass in g/mol from ``periodictable``."""
    symbol = element_symbol.strip().capitalize()
    element = getattr(periodictable, symbol, None)
    if element is None:
        msg = f"Unknown element symbol for atomic weight: {element_symbol!r}"
        raise ValueError(msg)
    return float(element.mass)


def compound_mu(
    energy_ev: NDArray[np.floating],
    atom_counts: dict[str, int],
    henke_tables: dict[str, pd.DataFrame],
) -> NDArray[np.float64]:
    r"""Compute compound mass absorption coefficient from Henke ``f2`` tables.

    Uses
    :math:`\mu = 2 r_e \lambda N_A \sum_i n_i f_{2,i} / M` with
    :math:`\lambda` in cm from photon energy in eV, matching the staging
    Henke processor.

    Parameters
    ----------
    energy_ev : numpy.ndarray
        Photon energies (eV) for the output grid.
    atom_counts : dict[str, int]
        Element symbol to atom count (for example ``{"C": 6, "H": 6}``).
    henke_tables : dict[str, pandas.DataFrame]
        Map from element symbol to tables from :func:`parse_henke_nff`.

    Returns
    -------
    numpy.ndarray
        Mass absorption coefficient (cm^2/g) on ``energy_ev``.

    Raises
    ------
    ValueError
        When ``atom_counts`` is empty or a required element table is missing.
    """
    if not atom_counts:
        msg = "atom_counts cannot be empty"
        raise ValueError(msg)
    energy_ev = np.asarray(energy_ev, dtype=np.float64)
    f2_sum = np.zeros_like(energy_ev, dtype=np.float64)
    molecular_weight = 0.0
    for element, count in atom_counts.items():
        if count <= 0:
            continue
        if element not in henke_tables:
            msg = f"Henke table missing for element {element!r}"
            raise ValueError(msg)
        table = henke_tables[element]
        valid = table.dropna(subset=["f2"])
        if len(valid) < 2:
            msg = f"Not enough valid f2 points for {element!r}"
            raise ValueError(msg)
        interp = interp1d(
            valid["energy_ev"].to_numpy(dtype=np.float64),
            valid["f2"].to_numpy(dtype=np.float64),
            kind="linear",
            bounds_error=False,
            fill_value="extrapolate",
        )
        f2_sum += float(count) * np.asarray(interp(energy_ev), dtype=np.float64)
        molecular_weight += float(count) * _atomic_weight(element)
    if molecular_weight <= 0.0:
        msg = "molecular weight from atom_counts is non-positive"
        raise ValueError(msg)
    lambda_cm = _EV_TO_WAVELENGTH_CM / energy_ev
    return (2.0 * _RE_CM * lambda_cm * _NA * f2_sum / molecular_weight).astype(
        np.float64
    )


def gaussian_step(
    energy_ev: NDArray[np.floating],
    center_ev: float,
    width_ev: float,
) -> NDArray[np.float64]:
    """Evaluate a unit-height Gaussian cumulative step (erf form).

    Parameters
    ----------
    energy_ev : numpy.ndarray
        Photon energy grid (eV).
    center_ev : float
        Step center (typically an ionization potential).
    width_ev : float
        Positive step width (eV); larger values give a broader edge.

    Returns
    -------
    numpy.ndarray
        Values in ``[0, 1]`` rising through ``center_ev``.

    Raises
    ------
    ValueError
        When ``width_ev`` is not positive.
    """
    if width_ev <= 0.0:
        msg = f"width_ev must be positive, got {width_ev}"
        raise ValueError(msg)
    energy_ev = np.asarray(energy_ev, dtype=np.float64)
    c = 2.0 * np.sqrt(2.0)
    return (0.5 + 0.5 * erf((energy_ev - center_ev) / (width_ev / c))).astype(
        np.float64
    )


def gaussian_step_edge(
    energy_ev: NDArray[np.floating],
    ionization_potentials_ev: NDArray[np.floating],
    *,
    width_ev: float = 0.5,
    total_jump: float = 1.0,
    weight_by_inverse_ip: bool = True,
) -> NDArray[np.float64]:
    """Build an IP-weighted sum of Gaussian steps on ``energy_ev``.

    Each ionization potential contributes a unit step scaled so the heights
    sum to ``total_jump``. When ``weight_by_inverse_ip`` is true, heights are
    proportional to ``1/IP`` (Igor-style), otherwise equal.

    Parameters
    ----------
    energy_ev : numpy.ndarray
        Photon energy grid (eV).
    ionization_potentials_ev : numpy.ndarray
        One IP (eV) per contributing site or atom.
    width_ev : float, optional
        Common step width (eV) for every IP.
    total_jump : float, optional
        Sum of step heights after the edge.
    weight_by_inverse_ip : bool, optional
        Weight step heights by inverse IP when true.

    Returns
    -------
    numpy.ndarray
        Step-edge absorption on ``energy_ev``.

    Raises
    ------
    ValueError
        When no finite IPs remain or ``total_jump`` / ``width_ev`` are invalid.
    """
    if width_ev <= 0.0:
        msg = f"width_ev must be positive, got {width_ev}"
        raise ValueError(msg)
    if total_jump < 0.0:
        msg = f"total_jump must be non-negative, got {total_jump}"
        raise ValueError(msg)
    energy_ev = np.asarray(energy_ev, dtype=np.float64)
    ips = np.asarray(ionization_potentials_ev, dtype=np.float64).ravel()
    ips = ips[np.isfinite(ips)]
    if ips.size == 0:
        msg = "ionization_potentials_ev has no finite values"
        raise ValueError(msg)
    if weight_by_inverse_ip:
        if np.any(ips <= 0.0):
            msg = "ionization potentials must be positive for inverse-IP weights"
            raise ValueError(msg)
        weights = 1.0 / ips
    else:
        weights = np.ones_like(ips)
    heights = total_jump * weights / weights.sum()
    edge = np.zeros_like(energy_ev, dtype=np.float64)
    for ip, height in zip(ips, heights, strict=True):
        edge += float(height) * gaussian_step(energy_ev, float(ip), width_ev)
    return edge
