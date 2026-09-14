"""dftrun sync: rsync a run directory to or from a remote StoBe host."""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console

from dftlearn.setup.remote import (
    DEFAULT_REMOTE_ROOT,
    KNOWN_REMOTE_HOSTS,
    pull_run_artifacts,
    pull_run_directory,
    remote_run_directory,
    resolve_ssh_host,
    sync_run_directory,
)

_CONSOLE = Console()


def sync_cmd(
    run_directory: Path = typer.Argument(
        ...,
        help="Local run directory to upload.",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
    ),
    remote: str = typer.Option(
        ...,
        "--remote",
        help=f"SSH host alias ({', '.join(sorted(KNOWN_REMOTE_HOSTS))}).",
    ),
    remote_root: str = typer.Option(
        DEFAULT_REMOTE_ROOT,
        "--remote-root",
        help="Remote base directory; created with mkdir -p when missing.",
    ),
    pull: bool = typer.Option(
        False,
        "--pull",
        help="Download logs and packaged_output from the remote run directory.",
    ),
    full: bool = typer.Option(
        False,
        "--full",
        help="When used with --pull, download the entire remote run directory.",
    ),
) -> None:
    """Rsync a run directory to or from a remote workstation.

    Upload (default) creates ``REMOTE_ROOT/RUN_NAME`` on the remote host before
    transferring files. ``--pull`` downloads ``logs/`` and ``packaged_output/``
    into the local run directory; add ``--full`` to download everything.
    """
    try:
        host = resolve_ssh_host(remote)
    except KeyError as exc:
        typer.echo(str(exc), err=True)
        raise typer.Exit(1) from exc

    remote_dir = remote_run_directory(host, run_directory.name, remote_root)
    local = run_directory.resolve()
    try:
        if pull:
            if full:
                pull_run_directory(local, host=host, remote_dir=remote_dir)
            else:
                pull_run_artifacts(local, host=host, remote_dir=remote_dir)
        else:
            sync_run_directory(local, host=host, remote_dir=remote_dir)
    except (RuntimeError, OSError) as exc:
        typer.echo(f"sync failed: {exc}", err=True)
        raise typer.Exit(1) from exc

    if pull:
        if full:
            _CONSOLE.print(f"Pulled full run {host}:{remote_dir}/\n  -> {local}/")
        else:
            _CONSOLE.print(
                f"Pulled logs/ and packaged_output/ from {host}:{remote_dir}/\n"
                f"  -> {local}/"
            )
    else:
        _CONSOLE.print(f"Synced {local}\n  -> {host}:{remote_dir}/")
