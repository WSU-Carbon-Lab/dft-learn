"""StoBe run-root discovery and composed multi-calc site loaders.

This module resolves XYZ geometry and per-site GND/EXC/TP / XAS / Xray paths
under a StoBe run directory and composes existing ``dftlearn.io`` parsers into
typed bundles. It does not replace ``python_pipeline.stobeLoader`` side effects
or Streamlit session state.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
from natsort import natsorted

if TYPE_CHECKING:
    import numpy as np

from dftlearn.io.stobe_final_energy import (
    collect_delta_ks_site_table,
    collect_final_energies_long,
)
from dftlearn.io.stobe_scf_convergence import (
    collect_scf_convergence_long,
    discover_site_stobe_out,
)
from dftlearn.io.stobe_xas_sticks import (
    parse_stobe_xas_dipole_sticks,
    parse_stobe_xas_sticks,
    site_xas_stick_paths,
)
from dftlearn.io.xray_out import (
    collect_site_xray_spectra,
    resolve_site_xray_path,
    site_spectra_to_long_frame,
)
from dftlearn.io.xyz_structure import xyz_rows_from_file

_SKIP_SITE_DIR_NAMES: frozenset[str] = frozenset(
    {
        "GND",
        "EXC",
        "TP",
        "NEXAFS",
        "packaged_output",
        "logs",
    }
)


@dataclass(frozen=True)
class StobeSitePaths:
    """Filesystem paths for one core-excited site under a StoBe run root.

    Attributes
    ----------
    site : str
        Site tag such as ``C1``.
    directory : pathlib.Path | None
        Immediate site subdirectory when present.
    gnd_out, exc_out, tp_out : pathlib.Path | None
        StoBe ``.out`` files for ground, excited, and transition-potential calcs.
    xas_sticks : pathlib.Path | None
        TP ``{site}.xas`` stick file when present.
    xray_out : pathlib.Path | None
        Broadened ``XrayT*.out`` (or layout equivalent) when present.
    """

    site: str
    directory: Path | None = None
    gnd_out: Path | None = None
    exc_out: Path | None = None
    tp_out: Path | None = None
    xas_sticks: Path | None = None
    xray_out: Path | None = None


@dataclass
class StobeRunBundle:
    """Composed tables and paths for one StoBe run directory.

    Attributes
    ----------
    run_root : pathlib.Path
        Resolved run root.
    xyz_path : pathlib.Path | None
        Geometry file when resolved.
    xyz_rows : list[tuple[str, float, float, float]]
        Parsed XYZ atom rows (label, x, y, z).
    site_paths : list[StobeSitePaths]
        Per-site path records.
    final_energies_long : pandas.DataFrame
        Long-form FINAL ENERGY rows (may be empty).
    delta_ks_sites : pandas.DataFrame
        Per-site Delta-KS table (may be empty).
    scf_convergence_long : pandas.DataFrame
        Long-form SCF iteration rows (may be empty).
    xray_energy_ev : numpy.ndarray | None
        Shared Xray energy axis when spectra load.
    xray_spectra : dict[str, numpy.ndarray]
        Per-site Xray intensity columns.
    xray_long : pandas.DataFrame
        Long-form Xray table (may be empty).
    xas_sticks : dict[str, numpy.ndarray]
        Per-site oscillator-strength stick tables (energy Ha, OS, ...).
    xas_dipole_sticks : dict[str, numpy.ndarray]
        Per-site dipole stick tables when parseable.
    """

    run_root: Path
    xyz_path: Path | None = None
    xyz_rows: list[tuple[str, float, float, float]] = field(default_factory=list)
    site_paths: list[StobeSitePaths] = field(default_factory=list)
    final_energies_long: pd.DataFrame = field(default_factory=pd.DataFrame)
    delta_ks_sites: pd.DataFrame = field(default_factory=pd.DataFrame)
    scf_convergence_long: pd.DataFrame = field(default_factory=pd.DataFrame)
    xray_energy_ev: np.ndarray | None = None
    xray_spectra: dict[str, np.ndarray] = field(default_factory=dict)
    xray_long: pd.DataFrame = field(default_factory=pd.DataFrame)
    xas_sticks: dict[str, np.ndarray] = field(default_factory=dict)
    xas_dipole_sticks: dict[str, np.ndarray] = field(default_factory=dict)


def find_xyz_files(directory: Path) -> list[Path]:
    """List ``*.xyz`` files directly under ``directory``.

    Parameters
    ----------
    directory : pathlib.Path
        Directory to scan (non-recursive).

    Returns
    -------
    list[pathlib.Path]
        Matching paths in filesystem order from ``Path.glob``.
    """
    return list(Path(directory).glob("*.xyz"))


def find_xyz(run_root: Path, *, xyz: Path | None = None) -> Path:
    """Resolve the XYZ geometry file for a StoBe run root.

    When ``xyz`` is given, that path is resolved and must exist. Otherwise the
    sole ``*.xyz`` under ``run_root`` is used, or the file whose stem matches
    the run directory name when several are present.

    Parameters
    ----------
    run_root : pathlib.Path
        StoBe run directory.
    xyz : pathlib.Path, optional
        Explicit geometry path.

    Returns
    -------
    pathlib.Path
        Resolved XYZ path.

    Raises
    ------
    FileNotFoundError
        When no XYZ file exists or the explicit path is missing.
    ValueError
        When multiple XYZ files exist and none match the run directory name.
    """
    if xyz is not None:
        path = Path(xyz).resolve()
        if not path.is_file():
            msg = f"XYZ geometry not found: {path}"
            raise FileNotFoundError(msg)
        return path

    run_root = Path(run_root).resolve()
    xyz_files = find_xyz_files(run_root)
    if not xyz_files:
        msg = f"No .xyz geometry under {run_root}"
        raise FileNotFoundError(msg)
    if len(xyz_files) == 1:
        return xyz_files[0].resolve()
    dir_name = run_root.name
    for xyz_file in xyz_files:
        if xyz_file.stem == dir_name:
            return xyz_file.resolve()
    names = ", ".join(p.name for p in xyz_files)
    msg = (
        f"Multiple .xyz files under {run_root} ({names}); "
        f"pass xyz= explicitly or name one {dir_name}.xyz"
    )
    raise ValueError(msg)


def discover_site_directories(run_root: Path) -> list[Path]:
    """Return immediate site subdirectories such as ``C1/`` under ``run_root``.

    A site directory name matches uppercase letters followed by digits and is
    not a reserved organize folder (``GND``, ``EXC``, ``TP``, ``NEXAFS``,
    ``packaged_output``, ``logs``).

    Parameters
    ----------
    run_root : pathlib.Path
        StoBe run directory.

    Returns
    -------
    list[pathlib.Path]
        Naturally sorted site directory paths.
    """
    run_root = Path(run_root).resolve()
    found: list[Path] = []
    for child in run_root.iterdir():
        if not child.is_dir():
            continue
        if child.name in _SKIP_SITE_DIR_NAMES or child.name.startswith("."):
            continue
        if re.fullmatch(r"[A-Z]+\d+", child.name):
            found.append(child.resolve())
    return natsorted(found, key=lambda p: p.name)


def discover_stobe_site_paths(
    run_root: Path,
    *,
    xray_filename: str = "XrayT001.out",
) -> list[StobeSitePaths]:
    """Build per-site path records from a StoBe run root.

    Combines site directories with ``discover_site_stobe_out`` for GND/EXC/TP
    and optional XAS / Xray files. Sites that appear only in category folders
    (no ``C1/`` directory) are still included.

    Parameters
    ----------
    run_root : pathlib.Path
        StoBe run directory.
    xray_filename : str, optional
        Preferred spectrum basename inside each site directory.

    Returns
    -------
    list[StobeSitePaths]
        Naturally sorted by site tag.
    """
    run_root = Path(run_root).resolve()
    by_site: dict[str, dict[str, Path | None]] = {}

    def ensure(site: str) -> dict[str, Path | None]:
        if site not in by_site:
            by_site[site] = {
                "directory": None,
                "gnd_out": None,
                "exc_out": None,
                "tp_out": None,
                "xas_sticks": None,
                "xray_out": None,
            }
        return by_site[site]

    for site_dir in discover_site_directories(run_root):
        slot = ensure(site_dir.name)
        slot["directory"] = site_dir

    for suffix, key in (("gnd", "gnd_out"), ("exc", "exc_out"), ("tp", "tp_out")):
        for site, path in discover_site_stobe_out(run_root, suffix):
            slot = ensure(site)
            if slot[key] is None:
                slot[key] = path

    try:
        for site, path in site_xas_stick_paths(run_root):
            slot = ensure(site)
            slot["xas_sticks"] = path
    except FileNotFoundError:
        pass

    for site in list(by_site):
        xray = resolve_site_xray_path(run_root, site, xray_filename=xray_filename)
        if xray is not None:
            by_site[site]["xray_out"] = xray

    records: list[StobeSitePaths] = []
    for site in natsorted(by_site):
        slot = by_site[site]
        records.append(
            StobeSitePaths(
                site=site,
                directory=slot["directory"],
                gnd_out=slot["gnd_out"],
                exc_out=slot["exc_out"],
                tp_out=slot["tp_out"],
                xas_sticks=slot["xas_sticks"],
                xray_out=slot["xray_out"],
            )
        )
    return records


def load_stobe_run(
    run_root: Path,
    *,
    xyz: Path | None = None,
    xray_filename: str = "XrayT001.out",
    require_xyz: bool = False,
) -> StobeRunBundle:
    """Load composed StoBe tables for ``run_root`` using library parsers.

    Missing optional artifacts leave empty frames or dicts rather than failing
    the whole load. XYZ resolution failures raise only when ``require_xyz`` is
    true or when ``xyz`` is passed explicitly.

    Parameters
    ----------
    run_root : pathlib.Path
        StoBe run directory.
    xyz : pathlib.Path, optional
        Explicit geometry path; otherwise auto-detect when possible.
    xray_filename : str, optional
        Preferred ``XrayT*.out`` basename inside site directories.
    require_xyz : bool, optional
        When true, raise if geometry cannot be resolved.

    Returns
    -------
    StobeRunBundle
        Paths plus composed DataFrames and arrays.

    Raises
    ------
    FileNotFoundError, ValueError
        When ``require_xyz`` is true and geometry is missing/ambiguous, or when
        an explicit ``xyz`` path is invalid.
    """
    run_root = Path(run_root).resolve()
    bundle = StobeRunBundle(run_root=run_root)
    bundle.site_paths = discover_stobe_site_paths(run_root, xray_filename=xray_filename)

    try:
        xyz_path = find_xyz(run_root, xyz=xyz)
        bundle.xyz_path = xyz_path
        bundle.xyz_rows = xyz_rows_from_file(xyz_path)
    except (FileNotFoundError, ValueError):
        if require_xyz or xyz is not None:
            raise

    try:
        bundle.final_energies_long = collect_final_energies_long(run_root)
    except (FileNotFoundError, ValueError):
        bundle.final_energies_long = pd.DataFrame()

    try:
        bundle.delta_ks_sites = collect_delta_ks_site_table(run_root)
    except (FileNotFoundError, ValueError):
        bundle.delta_ks_sites = pd.DataFrame()

    try:
        bundle.scf_convergence_long = collect_scf_convergence_long(run_root)
    except (FileNotFoundError, ValueError):
        bundle.scf_convergence_long = pd.DataFrame()

    try:
        energy, spectra = collect_site_xray_spectra(
            run_root, xray_filename=xray_filename
        )
        bundle.xray_energy_ev = energy
        bundle.xray_spectra = spectra
        bundle.xray_long = site_spectra_to_long_frame(energy, spectra)
    except (FileNotFoundError, ValueError):
        pass

    try:
        for site, path in site_xas_stick_paths(run_root):
            bundle.xas_sticks[site] = parse_stobe_xas_sticks(path)
            try:
                bundle.xas_dipole_sticks[site] = parse_stobe_xas_dipole_sticks(path)
            except ValueError:
                continue
    except FileNotFoundError:
        pass

    return bundle
