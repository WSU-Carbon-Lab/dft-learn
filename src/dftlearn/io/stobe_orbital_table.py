"""Parse StoBe orbital energy tables and derive level metrics.

The parser targets blocks titled ``ORBITAL ENERGIES (ALL VIRTUALS INCLUDED)``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import periodictable as pt

_ORBITAL_ENERGIES_HEADER = "ORBITAL ENERGIES (ALL VIRTUALS INCLUDED)"

_ORB_TABLE_ROW = re.compile(
    r"^\s*(\d+)\s+([\d.]+)\s+([-+]?\d+\.\d+)\s+\S+\s+\([^)]+\)\s+"
    r"([\d.]+)\s+([-+]?\d+\.\d+)\s+"
)

_OCC_FULL = 1.0 - 1e-5
_OCC_EMPTY = 1e-5
_OCC_HALF = 0.5
_HALF_TOL = 0.02

_K_SHELL_EDGE_BINDING_EV: dict[str, float] = {
    "C": 284.8,
    "N": 409.9,
    "O": 543.1,
    "F": 696.7,
    "Na": 1070.8,
    "Mg": 1305.0,
    "Al": 1559.6,
    "Si": 1839.9,
    "P": 2145.5,
    "S": 2472.0,
    "Cl": 2822.4,
    "K": 3608.4,
    "Ca": 4038.5,
    "Ti": 4966.4,
    "V": 5129.0,
    "Cr": 5989.0,
    "Mn": 6549.0,
    "Fe": 7112.0,
    "Co": 7709.0,
    "Ni": 8333.0,
    "Cu": 8979.0,
    "Zn": 9659.0,
}


def element_symbol_from_site(site: str) -> str:
    """Return the element prefix from a StoBe site tag (e.g. ``C1`` -> ``C``)."""
    m = re.match(r"^([A-Z][a-z]?)(\d+)$", site.strip())
    if not m:
        msg = f"Cannot parse element from site tag {site!r} (expected e.g. C1, Zn12)"
        raise ValueError(msg)
    return m.group(1)


def reference_k_shell_binding_ev(element_symbol: str) -> float:
    """Return reference K-shell absorption-edge binding energy (eV) for the absorber.

    Values are approximate literature/NIST-style edge energies used to locate the
    core KS level nearest ``-binding`` in eV. The element must exist in
    ``periodictable``; the numeric edge energy comes from an internal table
    (not all elements are covered).

    Parameters
    ----------
    element_symbol : str
        Element symbol (e.g. ``C``, ``Zn``).

    Returns
    -------
    float
        Binding energy in eV (positive).

    Raises
    ------
    ValueError
        If the symbol is unknown to ``periodictable`` or has no tabulated K-edge.
    """
    sym = element_symbol.strip()
    key = sym.upper() if len(sym) == 1 else sym[0].upper() + sym[1:].lower()
    try:
        pt.elements.symbol(key)
    except Exception as exc:
        msg = f"Unknown element symbol {element_symbol!r}"
        raise ValueError(msg) from exc
    if key not in _K_SHELL_EDGE_BINDING_EV:
        msg = f"No K-edge binding energy tabulated for {key}"
        raise ValueError(msg)
    return float(_K_SHELL_EDGE_BINDING_EV[key])


def parse_stobe_orbital_energies_table(path: Path) -> pd.DataFrame:
    """Parse the full orbital-energies block from one StoBe ``.out`` file.

    Parameters
    ----------
    path : pathlib.Path
        Path to a StoBe output file containing the orbital table.

    Returns
    -------
    pandas.DataFrame
        Columns ``level``, ``alpha_occ``, ``alpha_ev``, ``beta_occ``, ``beta_ev``.

    Raises
    ------
    FileNotFoundError
        If ``path`` is missing.
    ValueError
        If no data rows were parsed.
    """
    path = Path(path)
    if not path.is_file():
        msg = f"StoBe output not found: {path}"
        raise FileNotFoundError(msg)

    in_table = False
    saw_header = False
    rows: list[dict[str, float | int]] = []

    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if _ORBITAL_ENERGIES_HEADER in line:
                in_table = True
                saw_header = False
                continue
            if not in_table:
                continue
            stripped = line.strip()
            if "Spin alpha" in line and "Spin beta" in line:
                continue
            if "Occup." in line and "Energy(eV)" in line:
                saw_header = True
                continue
            if not saw_header:
                continue
            if stripped.startswith(("MULLIKEN", "ORBITAL ENERGIES")):
                break
            m = _ORB_TABLE_ROW.match(line)
            if not m:
                continue
            rows.append(
                {
                    "level": int(m.group(1)),
                    "alpha_occ": float(m.group(2)),
                    "alpha_ev": float(m.group(3)),
                    "beta_occ": float(m.group(4)),
                    "beta_ev": float(m.group(5)),
                }
            )

    if not rows:
        msg = f"No orbital energy rows parsed from {path}"
        raise ValueError(msg)
    return pd.DataFrame(rows).sort_values("level").reset_index(drop=True)


def _alpha_homo_lumo(
    df: pd.DataFrame, *, excited: bool = False
) -> tuple[pd.Series, pd.Series]:
    """Return alpha HOMO and LUMO.

    For ``excited=True``, the promoted electron can sit above the old HOMO; LUMO is the
    lowest **empty** alpha level strictly above the alpha HOMO energy (excludes the deep
    core hole). For ground-like manifolds (``excited=False``), LUMO is the global lowest
    empty alpha level.
    """
    full = df[df["alpha_occ"] >= _OCC_FULL]
    if full.empty:
        msg = "HOMO requires at least one fully occupied alpha level"
        raise ValueError(msg)
    homo = full.loc[full["alpha_ev"].idxmax()]
    if excited:
        empty = df[
            (df["alpha_occ"] <= _OCC_EMPTY)
            & (df["alpha_ev"] > float(homo["alpha_ev"]) + 1e-9)
        ]
    else:
        empty = df[df["alpha_occ"] <= _OCC_EMPTY]
    if empty.empty:
        msg = "LUMO requires at least one empty alpha level"
        raise ValueError(msg)
    lumo = empty.loc[empty["alpha_ev"].idxmin()]
    return homo, lumo


@dataclass(frozen=True)
class GndOrbitalReport:
    """Ground-state core / HOMO / LUMO characterization (alpha)."""

    binding_ref_ev: float
    core_level: int
    core_energy_ev: float
    homo_level: int
    homo_energy_ev: float
    lumo_level: int
    lumo_energy_ev: float


def analyze_gnd_orbitals(df: pd.DataFrame, binding_ref_ev: float) -> GndOrbitalReport:
    """Locate core level nearest ``-binding_ref_ev`` and alpha HOMO/LUMO."""
    target_ev = -float(binding_ref_ev)
    idx = (df["alpha_ev"] - target_ev).abs().idxmin()
    core = df.loc[idx]
    homo, lumo = _alpha_homo_lumo(df, excited=False)
    return GndOrbitalReport(
        binding_ref_ev=float(binding_ref_ev),
        core_level=int(core["level"]),
        core_energy_ev=float(core["alpha_ev"]),
        homo_level=int(homo["level"]),
        homo_energy_ev=float(homo["alpha_ev"]),
        lumo_level=int(lumo["level"]),
        lumo_energy_ev=float(lumo["alpha_ev"]),
    )


@dataclass(frozen=True)
class ExcOrbitalReport:
    """Excited-state core hole and valence frontier (alpha)."""

    core_level: int
    core_energy_ev: float
    homo_level: int
    homo_energy_ev: float
    lumo_level: int
    lumo_energy_ev: float
    homo_plus_one_level: int | None
    homo_plus_one_energy_ev: float | None


def analyze_exc_orbitals(df: pd.DataFrame, gnd: GndOrbitalReport) -> ExcOrbitalReport:
    """Core hole (first alpha-unoccupied level), HOMO/LUMO, and HOMO+1 (former LUMO)."""
    empty = df[df["alpha_occ"] < _OCC_EMPTY].sort_values("level")
    if empty.empty:
        msg = "Excited state: no empty alpha orbital"
        raise ValueError(msg)
    core = empty.iloc[0]
    homo, lumo = _alpha_homo_lumo(df, excited=True)
    hp1_level: int | None = None
    hp1_e: float | None = None
    row_lumo = df.loc[df["level"] == gnd.lumo_level]
    if not row_lumo.empty and float(row_lumo.iloc[0]["alpha_occ"]) >= _OCC_FULL:
        hp1_level = gnd.lumo_level
        hp1_e = float(row_lumo.iloc[0]["alpha_ev"])
    return ExcOrbitalReport(
        core_level=int(core["level"]),
        core_energy_ev=float(core["alpha_ev"]),
        homo_level=int(homo["level"]),
        homo_energy_ev=float(homo["alpha_ev"]),
        lumo_level=int(lumo["level"]),
        lumo_energy_ev=float(lumo["alpha_ev"]),
        homo_plus_one_level=hp1_level,
        homo_plus_one_energy_ev=hp1_e,
    )


@dataclass(frozen=True)
class TpOrbitalReport:
    """Transition-potential half core hole and valence frontier (alpha)."""

    core_level: int
    core_energy_ev: float
    core_alpha_occ: float
    homo_level: int
    homo_energy_ev: float
    lumo_level: int
    lumo_energy_ev: float


def analyze_tp_orbitals(df: pd.DataFrame) -> TpOrbitalReport:
    """Half-filled core (alpha ~0.5) and alpha HOMO/LUMO."""
    half = df[np.abs(df["alpha_occ"] - _OCC_HALF) < _HALF_TOL].sort_values("level")
    if half.empty:
        msg = "TP: no alpha ~0.5 core orbital"
        raise ValueError(msg)
    core = half.iloc[0]
    homo, lumo = _alpha_homo_lumo(df, excited=False)
    return TpOrbitalReport(
        core_level=int(core["level"]),
        core_energy_ev=float(core["alpha_ev"]),
        core_alpha_occ=float(core["alpha_occ"]),
        homo_level=int(homo["level"]),
        homo_energy_ev=float(homo["alpha_ev"]),
        lumo_level=int(lumo["level"]),
        lumo_energy_ev=float(lumo["alpha_ev"]),
    )


def format_final_energy_text(path: Path) -> str:
    """Return a short multiline summary of the first ``FINAL ENERGY`` block (Ha)."""
    from dftlearn.io.stobe_final_energy import parse_stobe_final_energy_tables

    try:
        tab = parse_stobe_final_energy_tables(path)
    except ValueError:
        return "FINAL ENERGY: (not found)"
    row = tab.loc[tab["block_index"].astype(int) == 0].iloc[0]
    lines = [
        f"Total energy (H) = {row.get('total_energy_h', float('nan')):.8f}",
        f"Nuc-nuc (H)      = {row.get('nuc_nuc_energy_h', float('nan')):.8f}",
        f"El-nuc (H)       = {row.get('el_nuc_energy_h', float('nan')):.8f}",
        f"Kinetic (H)      = {row.get('kinetic_energy_h', float('nan')):.8f}",
        f"Coulomb (H)      = {row.get('coulomb_energy_h', float('nan')):.8f}",
        f"Ex-cor (H)       = {row.get('ex_cor_energy_h', float('nan')):.8f}",
    ]
    return "\n".join(lines)
