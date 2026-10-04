"""Package StoBe run outputs into CSVs and summary figures."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

from dftlearn.io.stobe_run import find_xyz
from dftlearn.visualization.scf_diagnostics_figure import write_scf_diagnostics_bundle
from dftlearn.visualization.stobe_final_energy_figure import (
    write_stobe_final_energy_bundle,
)
from dftlearn.visualization.xas_cluster_figure import write_xas_cluster_report
from dftlearn.visualization.xas_reconstruction_figure import (
    write_xas_reconstruction_report,
)
from dftlearn.visualization.xas_site_figure import write_xas_site_report

if TYPE_CHECKING:
    from rich.console import Console


@dataclass
class PostprocessResult:
    """Paths written by :func:`package_stobe_run`.

    Attributes
    ----------
    packaged_dir : pathlib.Path
        Output directory (``packaged_output`` by default).
    written : list[pathlib.Path]
        Files successfully written, in stage order.
    skipped : list[str]
        Human-readable stage skip reasons (for example missing TP sticks).
    """

    packaged_dir: Path
    written: list[Path] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)


def package_stobe_run(
    run_root: Path,
    *,
    out: Path | None = None,
    xyz: Path | None = None,
    xray_file: str = "XrayT001.out",
    dpi: int = 150,
    os_percent: float | None = None,
    c3_symmetrize: bool = False,
    console: Console | None = None,
) -> PostprocessResult:
    """Write packaged CSVs and figures for a completed StoBe run directory.

    Stages: site Xray long table and figure; optional TP reconstruction and
    overlap clustering; SCF diagnostics; FINAL ENERGY / orbital summary.

    Parameters
    ----------
    run_root : pathlib.Path
        StoBe run root containing site folders (for example ``C1/``) and geometry.
    out : pathlib.Path, optional
        Output folder; defaults to ``run_root/packaged_output``.
    xyz : pathlib.Path, optional
        XYZ geometry file; defaults to auto-detect under ``run_root``.
    xray_file : str, optional
        Spectrum file name inside each site directory.
    dpi : int, optional
        PNG resolution for generated figures.
    os_percent : float, optional
        OS cutoff as percent of windowed max OS; defaults to elbow detection.
    c3_symmetrize : bool, optional
        Fold Cartesian TP dipoles under C3 in the Al-N/O molecular frame
        before reconstruction and clustering.
    console : rich.console.Console, optional
        When provided, prints write/skip messages; otherwise silent.

    Returns
    -------
    PostprocessResult
        Packaged directory, written paths, and skipped stage notes.

    Raises
    ------
    FileNotFoundError, ValueError
        When required inputs such as geometry or X-ray tables are missing.
    """
    run_root = Path(run_root).resolve()
    packaged = Path(out).resolve() if out else (run_root / "packaged_output")
    xyz_path = find_xyz(run_root, xyz=xyz)
    result = PostprocessResult(packaged_dir=packaged)

    def _note(message: str) -> None:
        if console is not None:
            console.print(message)

    csv_p, fig_p = write_xas_site_report(
        run_root,
        packaged,
        xyz_path,
        xray_filename=xray_file,
        dpi=dpi,
    )
    result.written.extend([csv_p, fig_p])
    _note(f"[green]Wrote[/green] {csv_p}")
    _note(f"[green]Wrote[/green] {fig_p}")

    rec = write_xas_reconstruction_report(
        run_root,
        packaged,
        xray_filename=xray_file,
        xyz_path=xyz_path,
        c3_symmetrize=c3_symmetrize,
        dpi=dpi,
    )
    if rec is not None:
        rec_csv, rec_metrics, rec_sticks, rec_tensor, rec_tensor_mean, rec_summary = rec
        for path in (
            rec_csv,
            rec_metrics,
            rec_sticks,
            rec_tensor,
            rec_tensor_mean,
            rec_summary,
        ):
            result.written.append(path)
            _note(f"[green]Wrote[/green] {path}")
        if c3_symmetrize:
            c3_path = packaged / "c3_frame.json"
            result.written.append(c3_path)
            _note(f"[green]Wrote[/green] {c3_path}")
    else:
        skip = "No {site}.xas stick files found (skipped TP XAS reconstruction)."
        result.skipped.append(skip)
        _note(f"[yellow]{skip}[/yellow]")

    clus = write_xas_cluster_report(
        run_root,
        packaged,
        xray_filename=xray_file,
        xyz_path=xyz_path,
        c3_symmetrize=c3_symmetrize,
        os_percent=os_percent,
        dpi=dpi,
    )
    if clus is not None:
        for path in clus:
            result.written.append(path)
            _note(f"[green]Wrote[/green] {path}")
    else:
        skip = "No {site}.xas stick files found (skipped overlap clustering)."
        result.skipped.append(skip)
        _note(f"[yellow]{skip}[/yellow]")

    diag = write_scf_diagnostics_bundle(run_root, packaged, dpi=dpi)
    if diag is not None:
        long_csv, metrics_csv, diag_png = diag
        for path in (long_csv, metrics_csv, diag_png):
            result.written.append(path)
            _note(f"[green]Wrote[/green] {path}")
    else:
        skip = (
            "No SCF convergence tables found "
            "(skipped scf_convergence_*.csv and scf_diagnostics.png)."
        )
        result.skipped.append(skip)
        _note(f"[yellow]{skip}[/yellow]")

    fe = write_stobe_final_energy_bundle(run_root, packaged, dpi=dpi)
    if fe is not None:
        fe_csv, fe_pngs = fe
        result.written.append(fe_csv)
        _note(f"[green]Wrote[/green] {fe_csv}")
        for path in fe_pngs:
            result.written.append(path)
            _note(f"[green]Wrote[/green] {path}")
    else:
        skip = (
            "No FINAL ENERGY blocks found "
            "(skipped stobe_final_energies.csv and "
            "stobe_orbital_energy_summary.png)."
        )
        result.skipped.append(skip)
        _note(f"[yellow]{skip}[/yellow]")

    return result
