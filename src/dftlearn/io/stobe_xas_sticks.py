"""StoBe ``fort.11`` / ``C*.xas`` stick spectra and ``*xas.inp`` broadening cards.

This module reads dipole oscillator-strength sticks written by the TP XRAY step
and the RANGE/POINTS/WIDTH cards consumed by ``xrayspec.x``. It does not broaden
spectra, apply Delta-KS shifts, or parse ``XrayT*.out`` tables.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from natsort import natsorted


@dataclass(frozen=True)
class XasSpecSettings:
    """``xrayspec.x`` energy grid and piecewise FWHM schedule.

    Attributes
    ----------
    energy_min_ev, energy_max_ev : float
        Inclusive photon-energy endpoints of ``RANGE`` (eV).
    n_points : int
        Number of samples in ``POINTS`` (must be >= 2).
    fwhm_low_ev, fwhm_high_ev : float
        Gaussian FWHM (eV) below ``e_break_low_ev`` and above
        ``e_break_high_ev``.
    e_break_low_ev, e_break_high_ev : float
        Breakpoints of the linear FWHM ramp (eV).
    """

    energy_min_ev: float
    energy_max_ev: float
    n_points: int
    fwhm_low_ev: float
    fwhm_high_ev: float
    e_break_low_ev: float
    e_break_high_ev: float


CARBON_XAS_SPEC = XasSpecSettings(
    energy_min_ev=280.0,
    energy_max_ev=320.0,
    n_points=2000,
    fwhm_low_ev=0.5,
    fwhm_high_ev=12.0,
    e_break_low_ev=288.0,
    e_break_high_ev=320.0,
)


def _read_stobe_xas_table(path: Path) -> np.ndarray:
    """Read every numeric field from a StoBe ``C*.xas`` stick file."""
    path = Path(path)
    if not path.is_file():
        msg = f"XAS stick file not found: {path}"
        raise FileNotFoundError(msg)
    rows: list[list[float]] = []
    expected: int | None = None
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if expected is None:
                parts = line.split()
                if len(parts) < 2 or parts[0].upper() != "XAS":
                    msg = f"Expected 'XAS <n>' header in {path}"
                    raise ValueError(msg)
                expected = int(parts[1])
                continue
            parts = [_fortran_float_token(tok) for tok in line.split()]
            if len(parts) < 2:
                msg = f"Expected at least two fields per stick row in {path}"
                raise ValueError(msg)
            rows.append([float(p) for p in parts])
    if expected is None:
        msg = f"No XAS header parsed from {path}"
        raise ValueError(msg)
    if not rows:
        msg = f"No stick rows parsed from {path}"
        raise ValueError(msg)
    if len(rows) != expected:
        msg = f"Stick count {len(rows)} != header n={expected} in {path}"
        raise ValueError(msg)
    n_cols = max(len(r) for r in rows)
    table = np.full((len(rows), n_cols), np.nan, dtype=np.float64)
    for i, row in enumerate(rows):
        table[i, : len(row)] = row
    return table


def parse_stobe_xas_sticks(path: Path) -> np.ndarray:
    """Load a StoBe ``C*.xas`` (``fort.11``) dipole stick table as ``float64``.

    The first line is ``XAS <n>``. Each following data line contributes one
    transition: column 0 is excitation energy in Hartree and column 1 is the
    dipole oscillator strength. Remaining fields (quadrupole and tensor
    components) are ignored. Fortran ``D``/``d`` exponents are accepted.

    Parameters
    ----------
    path : pathlib.Path
        Path to a ``C1.xas``-style stick file.

    Returns
    -------
    numpy.ndarray
        Shape ``(n, 2)`` with columns ``(energy_ha, oscillator_strength)`` in
        file order.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    ValueError
        If the header is missing, ``n`` disagrees with the row count, or a
        data line has fewer than two numeric fields.
    """
    return _read_stobe_xas_table(path)[:, :2]


def parse_stobe_xas_dipole_sticks(path: Path) -> np.ndarray:
    r"""Load energy, total OS, and Cartesian dipole matrix elements.

    Column 0 is energy (Hartree), column 1 is the StoBe dipole oscillator
    strength, and columns 3-5 are :math:`\\mu_x, \\mu_y, \\mu_z` in a.u.
    These satisfy :math:`f = (2/3) E_\\mathrm{Ha} |\\mu|^2`.

    Parameters
    ----------
    path : pathlib.Path
        Path to a ``C1.xas``-style stick file with at least six fields.

    Returns
    -------
    numpy.ndarray
        Shape ``(n, 5)``: ``(energy_ha, os, mux, muy, muz)``.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    ValueError
        If a row has fewer than six fields.
    """
    table = _read_stobe_xas_table(path)
    if table.shape[1] < 6:
        msg = f"Expected at least six fields for dipole components in {path}"
        raise ValueError(msg)
    return np.column_stack(
        (table[:, 0], table[:, 1], table[:, 3], table[:, 4], table[:, 5])
    )


def parse_stobe_xas_inp(path: Path) -> XasSpecSettings:
    """Parse ``RANGE``, ``POINTS``, and ``WIDTH`` from a StoBe ``*xas.inp``.

    ``WIDTH`` must be four numbers: FWHM below the first breakpoint, FWHM above
    the second, then the two breakpoints in eV, matching ``xrayspec.x``.

    Parameters
    ----------
    path : pathlib.Path
        Path to ``C1xas.inp`` or equivalent.

    Returns
    -------
    XasSpecSettings
        Grid and FWHM schedule from the file.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    ValueError
        If a required card is missing or malformed.
    """
    path = Path(path)
    if not path.is_file():
        msg = f"XAS input not found: {path}"
        raise FileNotFoundError(msg)
    energy_min: float | None = None
    energy_max: float | None = None
    n_points: int | None = None
    width_vals: list[float] | None = None
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            parts = raw.split()
            if not parts:
                continue
            key = parts[0].upper()
            if key == "RANGE" and len(parts) >= 3:
                energy_min = float(parts[1])
                energy_max = float(parts[2])
            elif key == "POINTS" and len(parts) >= 2:
                n_points = int(parts[1])
            elif key == "WIDTH" and len(parts) >= 5:
                width_vals = [float(p) for p in parts[1:5]]
    if energy_min is None or energy_max is None:
        msg = f"Missing RANGE card in {path}"
        raise ValueError(msg)
    if n_points is None:
        msg = f"Missing POINTS card in {path}"
        raise ValueError(msg)
    if width_vals is None:
        msg = f"Missing WIDTH card with four values in {path}"
        raise ValueError(msg)
    if n_points < 2:
        msg = f"POINTS must be >= 2 in {path}"
        raise ValueError(msg)
    return XasSpecSettings(
        energy_min_ev=energy_min,
        energy_max_ev=energy_max,
        n_points=n_points,
        fwhm_low_ev=width_vals[0],
        fwhm_high_ev=width_vals[1],
        e_break_low_ev=width_vals[2],
        e_break_high_ev=width_vals[3],
    )


def site_xas_stick_paths(run_root: Path) -> list[tuple[str, Path]]:
    """Discover ``site_dir / {site}.xas`` pairs directly under ``run_root``.

    Parameters
    ----------
    run_root : pathlib.Path
        StoBe run directory containing one subdirectory per core-excited site.

    Returns
    -------
    list[tuple[str, pathlib.Path]]
        Naturally sorted ``(site_tag, path)`` pairs.

    Raises
    ------
    FileNotFoundError
        If no matching stick files are found.
    """
    run_root = Path(run_root).resolve()
    found: list[tuple[str, Path]] = []
    for child in run_root.iterdir():
        if not child.is_dir():
            continue
        candidate = child / f"{child.name}.xas"
        if candidate.is_file():
            found.append((child.name, candidate.resolve()))
    if not found:
        msg = f"No site .xas stick files under immediate subdirectories of {run_root}"
        raise FileNotFoundError(msg)
    return natsorted(found, key=lambda pair: pair[0])


def site_xas_inp_path(site_dir: Path, site: str) -> Path | None:
    """Return ``{site}xas.inp`` under ``site_dir`` when that file exists.

    Parameters
    ----------
    site_dir : pathlib.Path
        Site subdirectory (for example ``C1/``).
    site : str
        Site tag used to build ``{site}xas.inp``.

    Returns
    -------
    pathlib.Path or None
        Absolute path when present.
    """
    candidate = Path(site_dir) / f"{site}xas.inp"
    if candidate.is_file():
        return candidate.resolve()
    return None


def _fortran_float_token(token: str) -> str:
    """Normalize a Fortran ``D`` exponent token to a Python float literal."""
    return token.replace("D", "E").replace("d", "e")
