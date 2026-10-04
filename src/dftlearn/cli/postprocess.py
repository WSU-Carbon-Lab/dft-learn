"""Post-process a StoBe run: X-ray CSV, TP XAS, clustering, SCF, energies."""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console

from dftlearn.pipeline.postprocess import package_stobe_run

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

    Thin wrapper around :func:`dftlearn.pipeline.package_stobe_run`.

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
    result = package_stobe_run(
        run_root,
        out=out,
        xyz=xyz,
        xray_file=xray_file,
        dpi=dpi,
        os_percent=os_percent,
        c3_symmetrize=c3_symmetrize,
        console=console or _CONSOLE,
    )
    return result.packaged_dir


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
