"""dftrun top-level Typer app and subcommands."""

from __future__ import annotations

import typer

from dftlearn.cli import build, init, postprocess, sync
from dftlearn.cli import remote_cmd as remote_mod
from dftlearn.cli import run as run_mod
from dftlearn.cli import setup as setup_mod

app = typer.Typer(help="dftrun: build and run StoBe DFT workflows.")

app.command(
    "build",
    help="Generate StoBe input files from molConfig + XYZ.",
)(build.build_cmd)
app.command(
    "init",
    help="Initialize a StoBe run from a CIF file (session, optional interactive review).",
)(init.init_cmd)
app.command(
    "setup",
    help="Browser setup lab: PubChem, 3D embed, relax, site labeling.",
)(setup_mod.setup_cmd)
app.command(
    "sync",
    help="Rsync a run directory to a remote StoBe workstation.",
)(sync.sync_cmd)
app.command(
    "postprocess",
    help="Collect XrayT*.out spectra, CSV, and XAS summary figure for a run directory.",
)(postprocess.postprocess_cmd)
app.add_typer(run_mod.run_app, name="run")
app.add_typer(remote_mod.remote_app, name="remote")
