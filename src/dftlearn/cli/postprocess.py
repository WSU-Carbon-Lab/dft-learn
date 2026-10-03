"""Post-process a StoBe run: X-ray CSV, TP XAS, clustering, SCF, energies."""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console

from dftlearn.cli.build import _auto_detect_xyz
from dftlearn.visualization.scf_diagnostics_figure import write_scf_diagnostics_bundle
from dftlearn.visualization.stobe_final_energy_figure import (
    write_stobe_final_energy_bundle,
)
from dftlearn.visualization.xas_cluster_figure import write_xas_cluster_report
from dftlearn.visualization.xas_reconstruction_figure import (
    write_xas_reconstruction_report,
)
from dftlearn.visualization.xas_site_figure import write_xas_site_report

_CONSOLE = Console()


def postprocess_run(
    run_root: Path,
    *,
    out: Path | None = None,
    xyz: Path | None = None,
    xray_file: str = "XrayT001.out",
    dpi: int = 150,
    os_percent: float | None = None,
    c3_symmetrize: bool = False,
    console: Console | None = None,
) -> Path:
    """Write packaged CSVs and figures for a completed StoBe run directory.

    Parameters
    ----------
    run_root
        StoBe run root containing site folders (for example ``C1/``) and geometry.
    out
        Output folder; defaults to ``run_root/packaged_output``.
    xyz
        XYZ geometry file; defaults to the sole ``.xyz`` under ``run_root``.
    xray_file
        Spectrum file name inside each site directory.
    dpi
        PNG resolution for generated figures.
    os_percent
        OS cutoff as percent of windowed max OS; defaults to elbow detection.
    c3_symmetrize
        Fold Cartesian TP dipoles under C3 in the Al-N/O molecular frame
        before reconstruction and clustering.
    console
        Rich console for status messages; defaults to a module-level console.

    Returns
    -------
    Path
        Resolved packaged output directory.

    Raises
    ------
    FileNotFoundError, ValueError
        When required inputs such as geometry or X-ray tables are missing.
    """
    out_console = console or _CONSOLE
    run_root = Path(run_root).resolve()
    packaged = Path(out).resolve() if out else (run_root / "packaged_output")
    xyz_path = Path(xyz).resolve() if xyz else _auto_detect_xyz(run_root).resolve()
    csv_p, fig_p = write_xas_site_report(
        run_root,
        packaged,
        xyz_path,
        xray_filename=xray_file,
        dpi=dpi,
    )
    out_console.print(f"[green]Wrote[/green] {csv_p}")
    out_console.print(f"[green]Wrote[/green] {fig_p}")
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
        out_console.print(f"[green]Wrote[/green] {rec_csv}")
        out_console.print(f"[green]Wrote[/green] {rec_metrics}")
        out_console.print(f"[green]Wrote[/green] {rec_sticks}")
        out_console.print(f"[green]Wrote[/green] {rec_tensor}")
        out_console.print(f"[green]Wrote[/green] {rec_tensor_mean}")
        out_console.print(f"[green]Wrote[/green] {rec_summary}")
        if c3_symmetrize:
            out_console.print(f"[green]Wrote[/green] {packaged / 'c3_frame.json'}")
    else:
        out_console.print(
            "[yellow]No {site}.xas stick files found "
            "(skipped TP XAS reconstruction).[/yellow]"
        )
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
            out_console.print(f"[green]Wrote[/green] {path}")
    else:
        out_console.print(
            "[yellow]No {site}.xas stick files found "
            "(skipped overlap clustering).[/yellow]"
        )
    diag = write_scf_diagnostics_bundle(run_root, packaged, dpi=dpi)
    if diag is not None:
        long_csv, metrics_csv, diag_png = diag
        out_console.print(f"[green]Wrote[/green] {long_csv}")
        out_console.print(f"[green]Wrote[/green] {metrics_csv}")
        out_console.print(f"[green]Wrote[/green] {diag_png}")
    else:
        out_console.print(
            "[yellow]No SCF convergence tables found (skipped scf_convergence_*.csv "
            "and scf_diagnostics.png).[/yellow]"
        )
    fe = write_stobe_final_energy_bundle(run_root, packaged, dpi=dpi)
    if fe is not None:
        fe_csv, fe_pngs = fe
        out_console.print(f"[green]Wrote[/green] {fe_csv}")
        for path in fe_pngs:
            out_console.print(f"[green]Wrote[/green] {path}")
    else:
        out_console.print(
            "[yellow]No FINAL ENERGY blocks found (skipped stobe_final_energies.csv "
            "and stobe_orbital_energy_summary.png).[/yellow]"
        )
    return packaged


def postprocess_cmd(
    run_root: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        help="StoBe run root containing site folders (e.g. C1/, C2/) and geometry.",
    ),
    out: Path | None = typer.Option(
        None,
        "--out",
        "-o",
        help="Output folder (default: RUN_ROOT/packaged_output).",
    ),
    xyz: Path | None = typer.Option(
        None,
        "--xyz",
        "-x",
        help="XYZ geometry file; default: sole .xyz under run_root.",
    ),
    xray_file: str = typer.Option(
        "XrayT001.out",
        "--xray-file",
        help="Spectrum file name inside each site directory.",
    ),
    dpi: int = typer.Option(150, "--dpi", min=72, max=600, help="PNG resolution."),
    os_percent: float | None = typer.Option(
        None,
        "--os-percent",
        min=0.0,
        max=100.0,
        help="OS cutoff as percent of windowed max OS (default: elbow).",
    ),
    c3_symmetrize: bool = typer.Option(
        False,
        "--c3-symmetrize",
        help=(
            "Fold Cartesian TP dipoles under C3 in the Al-N/O molecular frame "
            "before reconstruction and clustering (default: off)."
        ),
    ),
) -> None:
    """Extract X-ray tables, CSVs, XAS figure, SCF diagnostics, and final energies."""
    run_root = Path(run_root).resolve()
    try:
        postprocess_run(
            run_root,
            out=out,
            xyz=xyz,
            xray_file=xray_file,
            dpi=dpi,
            os_percent=os_percent,
            c3_symmetrize=c3_symmetrize,
        )
    except (FileNotFoundError, ValueError) as exc:
        _CONSOLE.print(f"[red]{exc}[/red]")
        raise typer.Exit(1) from exc
