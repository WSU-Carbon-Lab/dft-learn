"""Parse StoBe SCF iteration tables from ``*gnd.out``, ``*exc.out``, and ``*tp.out``.

Output discovery follows the same layouts as ``dftrun run organize``: (1) legacy
per-site directories ``SITE/SITE{gnd|exc|tp}.out``, (2) flat category trees
``GND/``, ``EXC/``, ``TP/`` at the run root, (3) a filtered recursive search for
``*{gnd|exc|tp}.out`` when files live elsewhere.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from natsort import natsorted

_CALC_SUFFIX_TO_CATEGORY: dict[str, str] = {"gnd": "GND", "exc": "EXC", "tp": "TP"}

_LEGACY_SKIP_SITE_DIR_NAMES: frozenset[str] = frozenset(
    {
        "GND",
        "EXC",
        "TP",
        "NEXAFS",
        "packaged_output",
    }
)

_RGLOB_SKIP_DIR_NAMES: frozenset[str] = frozenset(
    {
        "packaged_output",
        ".git",
        "__pycache__",
        ".venv",
        "node_modules",
        "build",
        "dist",
    }
)


def site_tag_from_stobe_out_filename(filename: str, calc_suffix: str) -> str | None:
    """Return the site label embedded in a StoBe ``.out`` name, or ``None``.

    StoBe names are ``{site}{gnd|exc|tp}.out`` (e.g. ``C1gnd.out`` -> ``C1``).

    Parameters
    ----------
    filename : str
        Basename only (no path).
    calc_suffix : str
        One of ``gnd``, ``exc``, ``tp``.

    Returns
    -------
    str | None
        Site tag when ``filename`` ends with ``{calc_suffix}.out`` with a
        non-empty prefix; otherwise ``None``.
    """
    tail = f"{calc_suffix}.out"
    if not filename.endswith(tail):
        return None
    tag = filename[: -len(tail)]
    return tag or None


def _rglob_path_is_skipped(path: Path, run_root: Path) -> bool:
    """Exclude tooling trees; do not skip ``GND``/``EXC``/``TP`` (organized outs)."""
    try:
        rel = path.relative_to(run_root)
    except ValueError:
        return True
    return any(p in _RGLOB_SKIP_DIR_NAMES or p.startswith(".") for p in rel.parts)


def parse_stobe_scf_convergence_table(path: Path) -> pd.DataFrame:
    """Extract the SCF iteration block from one StoBe ``.out`` file.

    Reads from the ``ITER / TOTAL ENERGY / ...`` header until ``SCF CONVERGED``
    or a dashed separator. Lines ``DIIS turned on...`` set ``diis_active`` true
    for all following data rows in that block.

    Parameters
    ----------
    path : pathlib.Path
        Path to a ``C1gnd.out``-style StoBe output file.

    Returns
    -------
    pandas.DataFrame
        Columns ``iteration``, ``total_energy_h``, ``decrease_h``,
        ``aver_density``, ``max_density``, ``diis_error``, ``cpu_s``,
        ``diis_active`` (bool).

    Raises
    ------
    FileNotFoundError
        If ``path`` is missing.
    ValueError
        If no iteration rows were parsed.
    """
    path = Path(path)
    if not path.is_file():
        msg = f"StoBe output not found: {path}"
        raise FileNotFoundError(msg)
    in_header = False
    diis_on = False
    rows: list[dict[str, object]] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if "SCF CONVERGED" in line:
                break
            if "ITER" in line and "TOTAL ENERGY" in line and "DECREASE" in line:
                in_header = True
                continue
            if not in_header:
                continue
            if "DIIS turned on" in line:
                diis_on = True
                continue
            stripped = line.strip()
            if stripped.startswith("---") and rows:
                break
            parts = line.split()
            if len(parts) < 7:
                continue
            if not parts[0].isdigit():
                continue
            try:
                it = int(parts[0])
                te = float(parts[1])
                dec = float(parts[2])
                avd = float(parts[3])
                mxd = float(parts[4])
                diis = float(parts[5])
                cpu = float(parts[6])
            except ValueError:
                continue
            rows.append(
                {
                    "iteration": it,
                    "total_energy_h": te,
                    "decrease_h": dec,
                    "aver_density": avd,
                    "max_density": mxd,
                    "diis_error": diis,
                    "cpu_s": cpu,
                    "diis_active": diis_on,
                }
            )
    if not rows:
        msg = f"No SCF iteration rows parsed from {path}"
        raise ValueError(msg)
    return pd.DataFrame(rows)


def discover_site_stobe_out(
    run_root: Path,
    calc_suffix: str,
) -> list[tuple[str, Path]]:
    """Find ``*`` ``{site}{calc_suffix}.out`` files for one calculation type.

    Resolution order (first path wins per ``site`` tag):

    (1) **Per-site folders** (legacy): ``run_root/SITE/SITE{calc}.out`` for each
    immediate subdirectory ``SITE`` that is not a reserved name (so ``GND`` is
    not treated as a site folder).

    (2) **Category folders** (after ``dftrun run organize``): ``run_root/GND``,
    ``EXC``, or ``TP`` containing flat ``SITE{calc}.out`` files.

    (3) **Recursive fallback**: ``run_root/**/SITE{calc}.out``, skipping only
    non-chemistry trees (``packaged_output``, ``.git``, ``__pycache__``, etc.);
    ``GND``/``EXC``/``TP`` segments are still searched.

    Parameters
    ----------
    run_root : pathlib.Path
        StoBe run directory (the folder passed to ``dftrun postprocess``).
    calc_suffix : str
        One of ``gnd``, ``exc``, ``tp``.

    Returns
    -------
    list[tuple[str, pathlib.Path]]
        Pairs ``(site_tag, path)`` sorted naturally by site tag.

    Raises
    ------
    ValueError
        If ``calc_suffix`` is not ``gnd``, ``exc``, or ``tp``.
    """
    if calc_suffix not in _CALC_SUFFIX_TO_CATEGORY:
        allowed = ", ".join(sorted(_CALC_SUFFIX_TO_CATEGORY))
        msg = f"calc_suffix must be one of {allowed}, got {calc_suffix!r}"
        raise ValueError(msg)

    run_root = Path(run_root).resolve()
    ordered: dict[str, Path] = {}

    def add_if_missing(site: str, path: Path) -> None:
        if site not in ordered:
            ordered[site] = path.resolve()

    for child in sorted(run_root.iterdir(), key=lambda p: p.name):
        if not child.is_dir():
            continue
        if child.name in _LEGACY_SKIP_SITE_DIR_NAMES or child.name.startswith("."):
            continue
        site = child.name
        candidate = child / f"{site}{calc_suffix}.out"
        if candidate.is_file():
            add_if_missing(site, candidate)

    cat_dir = run_root / _CALC_SUFFIX_TO_CATEGORY[calc_suffix]
    if cat_dir.is_dir():
        for f in sorted(cat_dir.glob(f"*{calc_suffix}.out")):
            if not f.is_file():
                continue
            site = site_tag_from_stobe_out_filename(f.name, calc_suffix)
            if site is not None:
                add_if_missing(site, f)

    pattern = f"*{calc_suffix}.out"
    for f in sorted(run_root.rglob(pattern)):
        if not f.is_file():
            continue
        if _rglob_path_is_skipped(f, run_root):
            continue
        site = site_tag_from_stobe_out_filename(f.name, calc_suffix)
        if site is not None:
            add_if_missing(site, f)

    return natsorted(ordered.items(), key=lambda p: p[0])


def collect_scf_convergence_long(
    run_root: Path,
    calc_types: tuple[str, ...] = ("gnd", "exc", "tp"),
) -> pd.DataFrame:
    """Load SCF tables for every site and calculation type into one long frame.

    Skips missing ``.out`` files without error (site simply absent for that
    type).

    Parameters
    ----------
    run_root : pathlib.Path
        StoBe run root: either per-site ``C1/``, ``C2/``, … and/or ``GND/``,
        ``EXC/``, ``TP/`` after ``dftrun run organize`` (see
        :func:`discover_site_stobe_out`).
    calc_types : tuple[str, ...], optional
        StoBe run-type suffixes to scan.

    Returns
    -------
    pandas.DataFrame
        Columns include ``site``, ``calc_type``, plus those from
        :func:`parse_stobe_scf_convergence_table`.

    Raises
    ------
    ValueError
        If no rows were collected from any file.
    """
    pieces: list[pd.DataFrame] = []
    for calc in calc_types:
        for site, out_path in discover_site_stobe_out(run_root, calc):
            try:
                block = parse_stobe_scf_convergence_table(out_path)
            except ValueError:
                continue
            block = block.copy()
            block["site"] = site
            block["calc_type"] = calc
            block["source_path"] = str(out_path)
            pieces.append(block)
    if not pieces:
        msg = f"No SCF convergence data found under {run_root}"
        raise ValueError(msg)
    return pd.concat(pieces, ignore_index=True)


def scf_convergence_auc_metrics(long_df: pd.DataFrame) -> pd.DataFrame:
    """Compute per-(site, calc_type) AUC-style scalars for bar charts.

    Uses the trapezoid rule: energy AUC is the area under ``|decrease_h|`` vs
    ``iteration``; density AUC is the area under ``max_density`` vs
    ``iteration``. ``diis_start_iter`` is the first iteration with
    ``diis_active`` true (NaN if never).

    Parameters
    ----------
    long_df : pandas.DataFrame
        Output of :func:`collect_scf_convergence_long`.

    Returns
    -------
    pandas.DataFrame
        One row per ``(site, calc_type)`` with ``energy_auc``, ``density_auc``,
        ``diis_start_iter``, ``n_iterations``.
    """
    out_rows: list[dict[str, object]] = []
    for (site, calc), grp in long_df.groupby(["site", "calc_type"], sort=False):
        grp = grp.sort_values("iteration")
        it = grp["iteration"].to_numpy(dtype=np.float64)
        dec = np.abs(grp["decrease_h"].to_numpy(dtype=np.float64))
        mxd = grp["max_density"].to_numpy(dtype=np.float64)
        energy_auc = float(np.trapz(dec, it))
        density_auc = float(np.trapz(mxd, it))
        diis_sub = grp.loc[grp["diis_active"], "iteration"]
        diis_start = float(diis_sub.min()) if len(diis_sub) else float("nan")
        out_rows.append(
            {
                "site": site,
                "calc_type": calc,
                "energy_auc": energy_auc,
                "density_auc": density_auc,
                "diis_start_iter": diis_start,
                "n_iterations": int(grp["iteration"].max()),
            }
        )
    return pd.DataFrame(out_rows)
