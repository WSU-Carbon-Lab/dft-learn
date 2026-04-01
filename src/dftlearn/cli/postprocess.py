"""Post-process a finished StoBe run directory: spectra CSV and XAS figures."""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console

from dftlearn.cli.build import _auto_detect_xyz
from dftlearn.visualization.xas_site_figure import write_xas_site_report

_CONSOLE = Console()


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
) -> None:
    """Extract site X-ray tables, write one CSV, and save a summary PNG figure."""
    run_root = Path(run_root).resolve()
    packaged = Path(out).resolve() if out else (run_root / "packaged_output")
    xyz_path = Path(xyz).resolve() if xyz else _auto_detect_xyz(run_root).resolve()
    try:
        csv_p, fig_p = write_xas_site_report(
            run_root,
            packaged,
            xyz_path,
            xray_filename=xray_file,
            dpi=dpi,
        )
    except (FileNotFoundError, ValueError) as exc:
        _CONSOLE.print(f"[red]{exc}[/red]")
        raise typer.Exit(1) from exc
    _CONSOLE.print(f"[green]Wrote[/green] {csv_p}")
    _CONSOLE.print(f"[green]Wrote[/green] {fig_p}")
