"""dftrun top-level Typer app and subcommands."""

from __future__ import annotations

import typer

from dftlearn.cli import build, postprocess
from dftlearn.cli import run as run_mod

app = typer.Typer(help="dftrun: build and run StoBe DFT workflows.")

app.command(
    "build",
    help="Generate StoBe input files from molConfig + XYZ.",
)(build.build_cmd)
app.command(
    "postprocess",
    help="Collect XrayT*.out spectra, CSV, and XAS summary figure for a run directory.",
)(postprocess.postprocess_cmd)
app.add_typer(run_mod.run_app, name="run")


@app.command("init")
def init_cmd() -> None:
    """Stub for future project init."""
    typer.echo("init: not implemented")
