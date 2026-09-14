"""Assign EXC/TP core-hole occupations from ground-state KS ionization energies.

StoBe ``FSYM scfocc excited`` needs the absorber 1s MO index, not a copied
ZnPc orbital window. After a converged GND calculation this module locates
that orbital as the occupied KS level whose ionization energy
(``-epsilon``) is nearest the element K-edge, then writes
``ALFA 0 1 N 0.0`` (EXC) and ``ALFA 0 1 N 0.5`` (TP) into each site run file.

This module does not run StoBe. Callers must finish GND first.
"""

from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from dftlearn.io.stobe_orbital_table import (
    analyze_gnd_orbitals,
    element_symbol_from_site,
    parse_stobe_orbital_energies_table,
    reference_k_shell_binding_ev,
)
from dftlearn.io.stobe_scf_convergence import discover_site_stobe_out

if TYPE_CHECKING:
    import pandas as pd

    from dftlearn.io.stobe_orbital_table import GndOrbitalReport

ASSIGNMENTS_FILENAME = "core_hole_assignments.json"
DEFAULT_CORE_WINDOW_EV = 80.0

_ALFA_OCC_LINE = re.compile(r"(?m)^ALFA\s+\d+\s+\d+\s+\d+\s+[0-9.]+\s*$")
_MOLCONFIG_ALFA_OCC = re.compile(
    r'^alfaOcc\s*=\s*"[^"]*"',
    re.MULTILINE,
)
_MOLCONFIG_ALFA_OCC_TP = re.compile(
    r'^alfaOccTP\s*=\s*"[^"]*"',
    re.MULTILINE,
)


@dataclass(frozen=True)
class CoreHoleAssignment:
    """Per-site absorber 1s index taken from a GND orbital table."""

    site: str
    core_level: int
    core_energy_ev: float
    ionization_energy_ev: float
    binding_ref_ev: float
    alfa_occ: str
    alfa_occ_tp: str
    gnd_out: str


def format_fsym_alfa_occupation(level: int, occupation: float) -> str:
    """Build a StoBe ``FSYM`` alpha occupation override for one MO.

    Parameters
    ----------
    level : int
        1-based orbital index in the GND KS table (must be >= 1).
    occupation : float
        Alpha occupation to impose (``0.0`` for EXC, ``0.5`` for TP).

    Returns
    -------
    str
        Tokens after ``ALFA``. StoBe reads ``0 nspec orbital occ``, so one
        orbital at ``level`` is ``0 1 {level} {occ}``.

    Raises
    ------
    ValueError
        If ``level`` is less than 1.
    """
    if level < 1:
        msg = f"Core-hole orbital index must be >= 1, got {level}"
        raise ValueError(msg)
    return f"0 1 {level} {occupation:.1f}"


def locate_k_edge_core_orbital(
    df: pd.DataFrame,
    binding_ref_ev: float,
    *,
    max_abs_error_ev: float = DEFAULT_CORE_WINDOW_EV,
) -> GndOrbitalReport:
    """Return the occupied KS 1s nearest the element K-edge ionization energy.

    Parameters
    ----------
    df : pandas.DataFrame
        Table from :func:`parse_stobe_orbital_energies_table`.
    binding_ref_ev : float
        Positive literature K-edge energy (eV); the search target is
        ``-binding_ref_ev``.
    max_abs_error_ev : float, optional
        Maximum allowed ``|epsilon + binding_ref_ev|``. Rejects Al/N/O 1s
        when the absorber is carbon.

    Returns
    -------
    GndOrbitalReport
        Core level plus HOMO/LUMO from :func:`analyze_gnd_orbitals`.

    Raises
    ------
    ValueError
        If no occupied match lies inside ``max_abs_error_ev``.
    """
    report = analyze_gnd_orbitals(df, binding_ref_ev)
    err = abs(report.core_energy_ev + float(binding_ref_ev))
    if err > max_abs_error_ev:
        msg = (
            f"Nearest occupied KS level {report.core_level} at "
            f"{report.core_energy_ev:.3f} eV is {err:.1f} eV from the "
            f"{binding_ref_ev:.1f} eV K-edge; expected the absorber 1s"
        )
        raise ValueError(msg)
    return report


def patch_fsym_alfa_occupation(text: str, occupation: str) -> str:
    """Replace the first four-token ``ALFA`` occupation line in a ``.run`` file.

    Parameters
    ----------
    text : str
        Full EXC or TP run-script contents.
    occupation : str
        Override tokens from :func:`format_fsym_alfa_occupation`.

    Returns
    -------
    str
        Patched script.

    Raises
    ------
    ValueError
        If no four-token ``ALFA`` line is present.
    """
    replacement = f"ALFA {occupation}"
    new, n = _ALFA_OCC_LINE.subn(replacement, text, count=1)
    if n != 1:
        msg = "No FSYM ALFA occupation line (ALFA i j k occ) found in run file"
        raise ValueError(msg)
    return new


def _site_run_file(run_root: Path, site: str, ftype: str) -> Path:
    return run_root / site / f"{site}{ftype}.run"


def _update_molconfig_occupations(
    molconfig: Path, alfa_occ: str, alfa_occ_tp: str
) -> None:
    text = molconfig.read_text(encoding="utf-8")
    text, n1 = _MOLCONFIG_ALFA_OCC.subn(f'alfaOcc = "{alfa_occ}"', text, count=1)
    text, n2 = _MOLCONFIG_ALFA_OCC_TP.subn(
        f'alfaOccTP = "{alfa_occ_tp}"', text, count=1
    )
    if n1 != 1 or n2 != 1:
        msg = f"Could not update alfaOcc lines in {molconfig}"
        raise ValueError(msg)
    molconfig.write_text(text, encoding="utf-8")


def assign_core_holes_from_gnd(
    run_root: Path,
    *,
    sites: list[str] | None = None,
    max_abs_error_ev: float = DEFAULT_CORE_WINDOW_EV,
) -> list[CoreHoleAssignment]:
    """Identify each absorber 1s from GND outputs and patch EXC/TP run files.

    For every site, reads ``GND/{site}gnd.out`` (or the per-site copy), finds
    the occupied KS orbital whose ionization energy is nearest the K-edge of
    the site element, and writes that index into ``{site}exc.run`` and
    ``{site}tp.run``. Writes ``core_hole_assignments.json`` at ``run_root``.
    When every site shares the same core index, also updates ``molConfig.py``.

    Parameters
    ----------
    run_root : pathlib.Path
        StoBe run directory containing site folders and GND outputs.
    sites : list of str, optional
        Restrict assignment to these site tags (e.g. ``C1``). ``None`` uses
        every GND output discovered under ``run_root``.
    max_abs_error_ev : float, optional
        Window around the K-edge used by :func:`locate_k_edge_core_orbital`.

    Returns
    -------
    list of CoreHoleAssignment
        One record per site, natsort order of discovery.

    Raises
    ------
    FileNotFoundError
        If a required GND output or EXC/TP run file is missing.
    ValueError
        If a GND orbital table cannot be parsed or the 1s match fails the
        energy window.
    """
    run_root = Path(run_root).resolve()
    discovered = dict(discover_site_stobe_out(run_root, "gnd"))
    if sites is not None:
        wanted = list(sites)
        missing = [s for s in wanted if s not in discovered]
        if missing:
            msg = (
                "GND output missing for site(s) "
                f"{', '.join(missing)}; run dftrun run gnd first"
            )
            raise FileNotFoundError(msg)
        pairs = [(s, discovered[s]) for s in wanted]
    else:
        pairs = list(discovered.items())
    if not pairs:
        msg = f"No GND outputs under {run_root}; run dftrun run gnd first"
        raise FileNotFoundError(msg)

    assignments: list[CoreHoleAssignment] = []
    for site, discovered_out in pairs:
        gnd_out = run_root / "GND" / f"{site}gnd.out"
        if not gnd_out.is_file():
            gnd_out = discovered_out
        element = element_symbol_from_site(site)
        binding = reference_k_shell_binding_ev(element)
        table = parse_stobe_orbital_energies_table(gnd_out)
        report = locate_k_edge_core_orbital(
            table, binding, max_abs_error_ev=max_abs_error_ev
        )
        alfa = format_fsym_alfa_occupation(report.core_level, 0.0)
        alfa_tp = format_fsym_alfa_occupation(report.core_level, 0.5)
        for ftype, occ in (("exc", alfa), ("tp", alfa_tp)):
            run_file = _site_run_file(run_root, site, ftype)
            if not run_file.is_file():
                msg = f"Missing {ftype} run file {run_file}"
                raise FileNotFoundError(msg)
            patched = patch_fsym_alfa_occupation(
                run_file.read_text(encoding="utf-8"), occ
            )
            run_file.write_text(patched, encoding="utf-8")
        assignments.append(
            CoreHoleAssignment(
                site=site,
                core_level=report.core_level,
                core_energy_ev=report.core_energy_ev,
                ionization_energy_ev=-report.core_energy_ev,
                binding_ref_ev=binding,
                alfa_occ=alfa,
                alfa_occ_tp=alfa_tp,
                gnd_out=str(gnd_out),
            )
        )

    payload = {
        "edge_element": element_symbol_from_site(assignments[0].site),
        "binding_ref_ev": assignments[0].binding_ref_ev,
        "sites": [asdict(row) for row in assignments],
    }
    dest = run_root / ASSIGNMENTS_FILENAME
    dest.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    levels = {row.core_level for row in assignments}
    molconfig = run_root / "molConfig.py"
    if molconfig.is_file() and len(levels) == 1:
        _update_molconfig_occupations(
            molconfig,
            assignments[0].alfa_occ,
            assignments[0].alfa_occ_tp,
        )
    return assignments
