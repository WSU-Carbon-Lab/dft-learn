"""dftrun setup: browser setup lab and legacy matplotlib fallback."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

_CONSOLE = Console()


def setup_cmd(
    run_directory: Path | None = typer.Argument(
        None,
        help="Optional run directory to reopen in the setup lab.",
        file_okay=False,
        dir_okay=True,
        readable=True,
    ),
    port: int = typer.Option(8501, "--port", help="Streamlit server port."),
    host: str = typer.Option("localhost", "--host", help="Bind address."),
    matplotlib: bool = typer.Option(
        False,
        "--matplotlib",
        help="Use legacy matplotlib editor instead of the browser lab.",
    ),
    remote: str | None = typer.Option(
        None,
        "--remote",
        help="Rsync finalized outputs to this SSH host alias after save (matplotlib only).",
    ),
    remote_root: str = typer.Option(
        "/home/hduva/projects/dft-runs",
        "--remote-root",
        help="Remote base directory when --remote is set.",
    ),
    no_preview: bool = typer.Option(
        False,
        "--no-preview",
        help="Skip writing init_preview.png on save (matplotlib only).",
    ),
) -> None:
    """Launch the browser setup lab (default) or legacy matplotlib editor.

    Opens a Streamlit workflow: PubChem search, 3D embed, UFF/anneal relaxation,
    and distinct-site labeling. Pass an existing run directory to reopen its session.
    """
    if run_directory is not None and run_directory.name == "serve":
        typer.echo(
            "Use `dftrun setup` to start the browser lab "
            "(the `serve` subcommand is no longer required).",
            err=True,
        )
        raise typer.Exit(1)
    if run_directory is not None and not run_directory.is_dir():
        typer.echo(f"Run directory does not exist: {run_directory}", err=True)
        raise typer.Exit(1)
    if matplotlib:
        _launch_matplotlib(
            run_directory,
            remote=remote,
            remote_root=remote_root,
            no_preview=no_preview,
        )
        return
    _launch_streamlit(run_directory, host=host, port=port)


def _launch_streamlit(
    run_directory: Path | None,
    *,
    host: str,
    port: int,
) -> None:
    lab_file = Path(__file__).resolve().parents[1] / "apps" / "setup_lab.py"
    cmd = [
        sys.executable,
        "-m",
        "streamlit",
        "run",
        str(lab_file),
        "--server.address",
        host,
        "--server.port",
        str(port),
    ]
    env = os.environ.copy()
    if run_directory is not None:
        env["DFTRUN_SETUP_RUN_DIR"] = str(run_directory.resolve())
    _CONSOLE.print(
        f"Starting setup lab at http://{host}:{port}\n"
        "Pipeline: PubChem search -> 3D embed -> UFF/anneal relax -> site labeler"
    )
    try:
        subprocess.run(cmd, check=True, env=env)
    except subprocess.CalledProcessError as exc:
        typer.echo(
            "Could not start Streamlit. Install with: uv sync --extra setup",
            err=True,
        )
        raise typer.Exit(exc.returncode) from exc


def _launch_matplotlib(
    run_directory: Path | None,
    *,
    remote: str | None = None,
    remote_root: str = "/home/hduva/projects/dft-runs",
    no_preview: bool = False,
) -> None:
    if run_directory is None:
        typer.echo("Matplotlib editor requires an existing run directory.", err=True)
        raise typer.Exit(1)

    from dftlearn.setup.interactive_editor import launch_setup_editor
    from dftlearn.setup.session import load_setup_session
    from dftlearn.setup.workflow import finalize_setup_session

    try:
        session = load_setup_session(run_directory)
    except FileNotFoundError as exc:
        typer.echo(f"setup failed: {exc}", err=True)
        raise typer.Exit(1) from exc

    def _save(updated_session) -> None:
        artifacts = finalize_setup_session(
            run_directory,
            updated_session,
            write_preview=not no_preview,
            remote_host=remote,
            remote_root=remote_root,
        )
        summary = json.loads(artifacts.summary_json.read_text(encoding="utf-8"))
        _CONSOLE.print(
            f"Saved {artifacts.geometry_xyz} and {artifacts.molconfig_py}\n"
            f"  core sites: {summary['n_core_sites']}"
        )

    try:
        launch_setup_editor(session, on_save=_save)
    except RuntimeError as exc:
        typer.echo(str(exc), err=True)
        raise typer.Exit(1) from exc

    summary_path = run_directory / "init_summary.json"
    if summary_path.is_file():
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        site_table = Table(title=f"Distinct {summary['edge_element']} sites")
        site_table.add_column("Site")
        site_table.add_column("Name")
        site_table.add_column("CIF")
        site_table.add_column("Class", justify="right")
        for row in summary.get("distinct_sites", []):
            site_table.add_row(
                row["site"],
                row.get("custom_name", ""),
                row["cif_label"],
                str(row["class_size"]),
            )
        _CONSOLE.print(site_table)
