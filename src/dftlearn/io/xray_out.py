"""StoBe ``XrayT*.out`` tables: parse, validate, and merge per-site spectra."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from natsort import natsorted


def parse_xray_out_table(path: Path) -> np.ndarray:
    """Load a two-column StoBe ``XrayT*.out`` spectrum as ``float64`` ``(n, 2)``.

    Each non-empty line must contain photon energy in the first field and a Fortran
    ``D``/``d`` exponent marker in the second; rows preserve file order.

    Parameters
    ----------
    path : pathlib.Path
        Path to an ``XrayT*.out`` file produced by StoBe.

    Returns
    -------
    numpy.ndarray
        Array with columns ``(energy, value)`` in the file's native units.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    ValueError
        If no numeric rows could be parsed or any row has fewer than two fields.
    """
    path = Path(path)
    if not path.is_file():
        msg = f"X-ray table not found: {path}"
        raise FileNotFoundError(msg)
    rows: list[tuple[float, float]] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                msg = f"Expected at least two fields per line in {path}"
                raise ValueError(msg)
            energy = float(parts[0])
            y_str = parts[1].replace("D", "E").replace("d", "e")
            rows.append((energy, float(y_str)))
    if not rows:
        msg = f"No data rows parsed from {path}"
        raise ValueError(msg)
    return np.asarray(rows, dtype=np.float64)


def site_xray_paths(
    run_root: Path,
    xray_filename: str = "XrayT001.out",
) -> list[tuple[str, Path]]:
    """Discover ``site_dir / xray_filename`` pairs directly under ``run_root``.

    ``site_dir`` names (for example ``C1``, ``C10``) are sorted with natural
    ordering. Only immediate child directories are scanned.

    Parameters
    ----------
    run_root : pathlib.Path
        StoBe run directory containing one subdirectory per core-excited site.
    xray_filename : str, optional
        Spectrum file name inside each site directory (default ``XrayT001.out``).

    Returns
    -------
    list[tuple[str, pathlib.Path]]
        Pairs ``(site_tag, path)`` where ``site_tag`` is the subdirectory name.

    Raises
    ------
    FileNotFoundError
        If no matching spectrum files are found.
    """
    run_root = Path(run_root).resolve()
    found: list[tuple[str, Path]] = []
    for child in run_root.iterdir():
        if not child.is_dir():
            continue
        candidate = child / xray_filename
        if candidate.is_file():
            found.append((child.name, candidate.resolve()))
    if not found:
        msg = f"No {xray_filename!r} under immediate subdirectories of {run_root}"
        raise FileNotFoundError(msg)
    return natsorted(found, key=lambda pair: pair[0])


def collect_site_xray_spectra(
    run_root: Path,
    xray_filename: str = "XrayT001.out",
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Load every site spectrum under ``run_root`` and align them to one energy grid.

    All tables must share the same energy axis within floating-point tolerance.

    Parameters
    ----------
    run_root : pathlib.Path
        Directory containing site subfolders (for example ``C1/``, ``C2/``).
    xray_filename : str, optional
        Spectrum file name inside each site directory.

    Returns
    -------
    energy : numpy.ndarray
        Shape ``(n,)``, first column of the reference site table.
    spectra : dict[str, numpy.ndarray]
        Maps each site tag to the second column, shape ``(n,)``.

    Raises
    ------
    FileNotFoundError
        If no spectrum files are found.
    ValueError
        If energy axes differ between sites or a table is empty.
    """
    pairs = site_xray_paths(run_root, xray_filename=xray_filename)
    first_tag, first_path = pairs[0]
    ref = parse_xray_out_table(first_path)
    energy = ref[:, 0].copy()
    spectra: dict[str, np.ndarray] = {first_tag: ref[:, 1].copy()}
    for site, path in pairs[1:]:
        block = parse_xray_out_table(path)
        if block.shape[0] != energy.shape[0]:
            msg = (
                f"Row count mismatch for site {site!r}: "
                f"{block.shape[0]} vs reference {energy.shape[0]} ({first_tag!r})"
            )
            raise ValueError(msg)
        if not np.allclose(block[:, 0], energy, rtol=0.0, atol=1e-6):
            msg = f"Energy axis mismatch for site {site!r} vs {first_tag!r}"
            raise ValueError(msg)
        spectra[site] = block[:, 1].copy()
    return energy, spectra


def site_spectra_to_long_frame(
    energy: np.ndarray,
    spectra: dict[str, np.ndarray],
) -> pd.DataFrame:
    """Stack per-site columns into a long-form table suitable for CSV export.

    Parameters
    ----------
    energy : numpy.ndarray
        Common photon energy axis, shape ``(n,)``.
    spectra : dict[str, numpy.ndarray]
        Per-site intensity (or cross-section) values, each shape ``(n,)``.

    Returns
    -------
    pandas.DataFrame
        Columns ``energy_ev``, ``abs``, ``site`` (site tag such as ``C1``).
        ``abs`` holds absorption strength from the StoBe table (second column).
    """
    pieces: list[pd.DataFrame] = []
    for site in natsorted(spectra.keys()):
        y = spectra[site]
        pieces.append(
            pd.DataFrame(
                {
                    "energy_ev": energy,
                    "abs": y,
                    "site": site,
                }
            )
        )
    return pd.concat(pieces, ignore_index=True)
