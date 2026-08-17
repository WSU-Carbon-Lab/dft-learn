"""Parse StoBe ``FINAL ENERGY / CHARGE / GEOMETRY RESULTS`` blocks from ``.out`` files.

Blocks appear immediately after ``SCF CONVERGED`` in ``*gnd.out``, ``*exc.out``, and
``*tp.out``. Per-site path discovery matches :mod:`dftlearn.io.stobe_scf_convergence`.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from natsort import natsorted

from dftlearn.io.stobe_scf_convergence import discover_site_stobe_out

HA_TO_EV = 27.21138602432245

_FINAL_HEADER = "FINAL ENERGY / CHARGE / GEOMETRY RESULTS"

_ORBITAL_ENERGIES_HEADER = "ORBITAL ENERGIES (ALL VIRTUALS INCLUDED)"

_ORB_TABLE_ROW = re.compile(
    r"^\s*(\d+)\s+([\d.]+)\s+([-+]?\d+\.\d+)\s+\S+\s+\([^)]+\)\s+"
    r"([\d.]+)\s+([-+]?\d+\.\d+)\s+"
)

_UNOCC_EPS = 1e-6

_FLOAT = r"-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?"

_LINE_PATTERNS: tuple[tuple[re.Pattern[str], str], ...] = (
    (re.compile(rf"Total energy\s+\(H\)\s*=\s*({_FLOAT})"), "total_energy_h"),
    (re.compile(rf"Nuc-nuc energy\s+\(H\)\s*=\s*({_FLOAT})"), "nuc_nuc_energy_h"),
    (re.compile(rf"El-nuc energy\s+\(H\)\s*=\s*({_FLOAT})"), "el_nuc_energy_h"),
    (re.compile(rf"Kinetic energy\s+\(H\)\s*=\s*({_FLOAT})"), "kinetic_energy_h"),
    (re.compile(rf"Coulomb energy\s+\(H\)\s*=\s*({_FLOAT})"), "coulomb_energy_h"),
    (re.compile(rf"Ex-cor energy\s+\(H\)\s*=\s*({_FLOAT})"), "ex_cor_energy_h"),
    (
        re.compile(
            rf"<Rho/r12/Rhof>-<Rhof/r12/Rhof>/2\s+\(H\)\s*=\s*({_FLOAT})"
        ),
        "rho_r12_diff_h",
    ),
    (
        re.compile(rf"<Rho/r12/Rhof>/2\s+\(H\)\s*=\s*({_FLOAT})"),
        "rho_r12_half_h",
    ),
    (
        re.compile(rf"Total exchange energy\s+\(H\)\s*=\s*({_FLOAT})"),
        "total_exchange_energy_h",
    ),
    (
        re.compile(rf"Total correlation energy\s+\(H\)\s*=\s*({_FLOAT})"),
        "total_correlation_energy_h",
    ),
)


def _try_parse_line(line: str) -> tuple[str, float] | None:
    """Map one output line to ``(column, value)`` when it matches a known field."""
    for pat, col in _LINE_PATTERNS:
        m = pat.search(line)
        if m:
            return (col, float(m.group(1)))
    return None


def _should_stop_after_line(stripped: str) -> bool:
    """Return whether parsing should end the current energy-component block."""
    if not stripped:
        return False
    return stripped.startswith(
        ("Decomposition of exchange", "GEOMETRY", "Total electron charge")
    )


def parse_stobe_final_energy_tables(path: Path) -> pd.DataFrame:
    """Extract every ``FINAL ENERGY / CHARGE / GEOMETRY RESULTS`` block in one file.

    Reads line-by-line; stops each block at ``Decomposition of exchange``, ``GEOMETRY``,
    or ``Total electron charge``. Multiple blocks (if present) each produce one row with
    increasing ``block_index``.

    Parameters
    ----------
    path : pathlib.Path
        Path to a StoBe ``C1gnd.out``-style output file.

    Returns
    -------
    pandas.DataFrame
        Rows include ``block_index`` and component columns in Hartree (``*_h``), or
        ``NaN`` when a field is absent in that block.

    Raises
    ------
    FileNotFoundError
        If ``path`` is missing.
    ValueError
        If no parsable blocks were found.
    """
    path = Path(path)
    if not path.is_file():
        msg = f"StoBe output not found: {path}"
        raise FileNotFoundError(msg)

    rows: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    block_index = 0

    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if _FINAL_HEADER in line:
                if current is not None and "total_energy_h" in current:
                    rows.append(current)
                current = {"block_index": block_index}
                block_index += 1
                continue
            if current is None:
                continue
            stripped = line.strip()
            if _should_stop_after_line(stripped):
                rows.append(current)
                current = None
                continue
            parsed = _try_parse_line(line)
            if parsed:
                col, val = parsed
                current[col] = val

    if current is not None and "total_energy_h" in current:
        rows.append(current)

    if not rows:
        msg = f"No FINAL ENERGY blocks parsed from {path}"
        raise ValueError(msg)

    cols_order = [
        "block_index",
        "total_energy_h",
        "nuc_nuc_energy_h",
        "el_nuc_energy_h",
        "kinetic_energy_h",
        "coulomb_energy_h",
        "ex_cor_energy_h",
        "rho_r12_diff_h",
        "rho_r12_half_h",
        "total_exchange_energy_h",
        "total_correlation_energy_h",
    ]
    df = pd.DataFrame(rows)
    for c in cols_order:
        if c not in df.columns:
            df[c] = np.nan
    extra = [c for c in df.columns if c not in cols_order]
    return df[cols_order + extra]


def parse_stobe_tp_lumo_alpha_ev(path: Path) -> float:
    """Return alpha-spin LUMO energy (eV) from a StoBe ``*tp.out`` file.

    Reads the ``ORBITAL ENERGIES (ALL VIRTUALS INCLUDED)`` table and returns the
    alpha orbital energy for the **first** level with alpha occupation numerically
    zero (lowest-by-index unoccupied row). StoBe lists levels in energy order, so
    this matches the lowest unoccupied alpha KS level without scanning all virtuals.

    Parameters
    ----------
    path : pathlib.Path
        Path to a ``C1tp.out``-style transition-potential output.

    Returns
    -------
    float
        LUMO energy in eV.

    Raises
    ------
    FileNotFoundError
        If ``path`` is missing.
    ValueError
        If the table or unoccupied rows are not found.
    """
    path = Path(path)
    if not path.is_file():
        msg = f"StoBe output not found: {path}"
        raise FileNotFoundError(msg)

    in_table = False
    saw_header = False

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
            alpha_occ = float(m.group(2))
            alpha_ev = float(m.group(3))
            if alpha_occ <= _UNOCC_EPS:
                return float(alpha_ev)

    msg = f"No alpha LUMO (unoccupied) row parsed from {path}"
    raise ValueError(msg)


def enrich_final_energies_delta_vs_gnd(df: pd.DataFrame) -> pd.DataFrame:
    """Add ``delta_vs_gnd_h`` and ``delta_vs_gnd_ev`` from same-site ground totals.

    For ``calc_type`` ``gnd``, both deltas are ``NaN``. For ``exc`` and ``tp``, subtract
    that site's ``total_energy_h`` at ``block_index`` 0 in the ground-state file when
    present; otherwise ``NaN``.

    Parameters
    ----------
    df : pandas.DataFrame
        Output of :func:`collect_final_energies_long` (must include ``site``,
        ``calc_type``, ``block_index``, ``total_energy_h``).

    Returns
    -------
    pandas.DataFrame
        A copy of ``df`` with ``delta_vs_gnd_h`` and ``delta_vs_gnd_ev`` appended.
    """
    out = df.copy()
    gnd = out[
        (out["calc_type"] == "gnd") & (out["block_index"].astype(int) == 0)
    ][["site", "total_energy_h"]].drop_duplicates(subset=["site"], keep="first")
    gmap = gnd.set_index("site")["total_energy_h"]

    def delta_h(row: pd.Series) -> float:
        if row["calc_type"] == "gnd":
            return float("nan")
        site = row["site"]
        if site not in gmap.index:
            return float("nan")
        base = float(gmap.loc[site])
        return float(row["total_energy_h"]) - base

    out["delta_vs_gnd_h"] = out.apply(delta_h, axis=1)
    out["delta_vs_gnd_ev"] = out["delta_vs_gnd_h"] * HA_TO_EV
    return out


def collect_final_energies_long(
    run_root: Path,
    calc_types: tuple[str, ...] = ("gnd", "exc", "tp"),
) -> pd.DataFrame:
    """Load final-energy component tables for every site and calculation type.

    Skips ``.out`` files that yield no parsable blocks (site absent for that type).

    Parameters
    ----------
    run_root : pathlib.Path
        StoBe run root (same layout as
        :func:`~dftlearn.io.stobe_scf_convergence.discover_site_stobe_out`).
    calc_types : tuple[str, ...], optional
        Suffixes ``gnd``, ``exc``, ``tp``.

    Returns
    -------
    pandas.DataFrame
        Columns include ``site``, ``calc_type``, ``source_path``, block fields from
        :func:`parse_stobe_final_energy_tables`, plus ``delta_vs_gnd_h`` and
        ``delta_vs_gnd_ev`` from :func:`enrich_final_energies_delta_vs_gnd`.

    Raises
    ------
    ValueError
        If no rows were collected from any file.
    """
    pieces: list[pd.DataFrame] = []
    for calc in calc_types:
        for site, out_path in discover_site_stobe_out(run_root, calc):
            try:
                block = parse_stobe_final_energy_tables(out_path)
            except ValueError:
                continue
            block = block.copy()
            block["site"] = site
            block["calc_type"] = calc
            block["source_path"] = str(out_path)
            pieces.append(block)
    if not pieces:
        msg = f"No FINAL ENERGY blocks found under {run_root}"
        raise ValueError(msg)
    merged = pd.concat(pieces, ignore_index=True)
    return enrich_final_energies_delta_vs_gnd(merged)


def final_energy_site_summary(df: pd.DataFrame) -> pd.DataFrame:
    """One row per ``(site, calc_type)`` using ``block_index`` 0 when duplicated.

    Parameters
    ----------
    df : pandas.DataFrame
        Long table from :func:`collect_final_energies_long` or
        :func:`enrich_final_energies_delta_vs_gnd`.

    Returns
    -------
    pandas.DataFrame
        Filtered to ``block_index == 0`` and sorted naturally by ``site``, then
        ``calc_type`` order gnd, exc, tp.
    """
    sub = df.loc[df["block_index"].astype(int) == 0].copy()
    if sub.empty:
        return sub
    sites = natsorted(sub["site"].unique())
    calc_order = ("gnd", "exc", "tp")
    parts: list[pd.DataFrame] = []
    for site in sites:
        for calc in calc_order:
            chunk = sub[(sub["site"] == site) & (sub["calc_type"] == calc)]
            if not chunk.empty:
                parts.append(chunk)
    return pd.concat(parts, ignore_index=True)


def _first_total_energy_h(parsed: pd.DataFrame) -> float:
    """Return ``total_energy_h`` for ``block_index`` 0."""
    sub = parsed.loc[parsed["block_index"].astype(int) == 0, "total_energy_h"]
    return float(sub.iloc[0])


def collect_delta_ks_site_table(run_root: Path) -> pd.DataFrame:
    """Build one row per site: SCF totals and Delta-KS alignment (eV).

    Uses ``FINAL ENERGY`` totals from ``*gnd.out``, ``*exc.out``, ``*tp.out`` (first
    block only) and alpha LUMO (eV) from the TP orbital table. Computes
    ``E_c_deltaKS_ev = E_e_ev - E_g_ev - lumo_tp_alpha_ev`` with SCF totals converted
    from Ha to eV and ``E^l`` the alpha LUMO eigenvalue (eV) from the TP output.

    Parameters
    ----------
    run_root : pathlib.Path
        StoBe run root (see :func:`discover_site_stobe_out`).

    Returns
    -------
    pandas.DataFrame
        Columns ``site``, ``E_g_ev``, ``E_e_ev``, ``E_tp_ev``, ``lumo_tp_alpha_ev``,
        ``delta_exc_vs_gnd_ev``, ``delta_tp_vs_gnd_ev``, ``E_c_deltaKS_ev``, and
        ``source_*`` paths for gnd, exc, tp.

    Raises
    ------
    ValueError
        If no StoBe site outputs are found under ``run_root``.
    """
    run_root = Path(run_root).resolve()
    gnd_paths = dict(discover_site_stobe_out(run_root, "gnd"))
    exc_paths = dict(discover_site_stobe_out(run_root, "exc"))
    tp_paths = dict(discover_site_stobe_out(run_root, "tp"))

    keys_g = set(gnd_paths.keys())
    keys_e = set(exc_paths.keys())
    keys_t = set(tp_paths.keys())
    sites = natsorted(keys_g | keys_e | keys_t)
    if not sites:
        msg = f"No site StoBe outputs found under {run_root}"
        raise ValueError(msg)
    rows: list[dict[str, object]] = []

    for site in sites:
        g_path = gnd_paths.get(site)
        e_path = exc_paths.get(site)
        t_path = tp_paths.get(site)
        row: dict[str, object] = {"site": site}
        if g_path is not None:
            row["source_gnd"] = str(g_path)
            try:
                gdf = parse_stobe_final_energy_tables(g_path)
                eg = _first_total_energy_h(gdf)
                row["E_g_ev"] = eg * HA_TO_EV
            except (ValueError, IndexError):
                row["E_g_ev"] = float("nan")
        else:
            row["source_gnd"] = ""
            row["E_g_ev"] = float("nan")

        if e_path is not None:
            row["source_exc"] = str(e_path)
            try:
                edf = parse_stobe_final_energy_tables(e_path)
                ee = _first_total_energy_h(edf)
                row["E_e_ev"] = ee * HA_TO_EV
            except (ValueError, IndexError):
                row["E_e_ev"] = float("nan")
        else:
            row["source_exc"] = ""
            row["E_e_ev"] = float("nan")

        if t_path is not None:
            row["source_tp"] = str(t_path)
            try:
                tdf = parse_stobe_final_energy_tables(t_path)
                et = _first_total_energy_h(tdf)
                row["E_tp_ev"] = et * HA_TO_EV
            except (ValueError, IndexError):
                row["E_tp_ev"] = float("nan")
            try:
                row["lumo_tp_alpha_ev"] = parse_stobe_tp_lumo_alpha_ev(t_path)
            except ValueError:
                row["lumo_tp_alpha_ev"] = float("nan")
        else:
            row["source_tp"] = ""
            row["E_tp_ev"] = float("nan")
            row["lumo_tp_alpha_ev"] = float("nan")

        eg = row.get("E_g_ev", float("nan"))
        ee = row.get("E_e_ev", float("nan"))
        etp = row.get("E_tp_ev", float("nan"))
        el = row.get("lumo_tp_alpha_ev", float("nan"))

        fin_g = isinstance(eg, float) and np.isfinite(eg)
        fin_e = isinstance(ee, float) and np.isfinite(ee)
        fin_tp = isinstance(etp, float) and np.isfinite(etp)
        if fin_g and fin_e:
            row["delta_exc_vs_gnd_ev"] = float(ee - eg)
        else:
            row["delta_exc_vs_gnd_ev"] = float("nan")
        if fin_g and fin_tp:
            row["delta_tp_vs_gnd_ev"] = float(etp - eg)
        else:
            row["delta_tp_vs_gnd_ev"] = float("nan")
        if (
            isinstance(ee, float)
            and isinstance(eg, float)
            and isinstance(el, float)
            and np.isfinite(ee)
            and np.isfinite(eg)
            and np.isfinite(el)
        ):
            row["E_c_deltaKS_ev"] = float(ee - eg - el)
        else:
            row["E_c_deltaKS_ev"] = float("nan")

        rows.append(row)

    if not rows:
        msg = f"No sites assembled under {run_root}"
        raise ValueError(msg)
    return pd.DataFrame(rows)
