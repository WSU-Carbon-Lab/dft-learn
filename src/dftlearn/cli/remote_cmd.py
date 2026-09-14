"""dftrun remote: install and manage remote StoBe workstations."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import typer
from rich.console import Console

from dftlearn.setup.remote import (
    DEFAULT_REMOTE_ROOT,
    KNOWN_REMOTE_HOSTS,
    inspect_remote_job,
    install_dftrun_on_host,
    pull_run_artifacts,
    remote_run_directory,
    resolve_ssh_host,
    run_remote_diagnose,
    run_remote_reset,
)

if TYPE_CHECKING:
    from dftlearn.setup.remote import RemoteJobStatus

_CONSOLE = Console()
remote_app = typer.Typer(help="Install dftrun and manage remote StoBe hosts.")


def _project_root(source: Path) -> Path:
    resolved = source.resolve()
    if (resolved / "pyproject.toml").is_file() and (
        resolved / "src" / "dftlearn"
    ).is_dir():
        return resolved
    if (resolved.parent / "pyproject.toml").is_file():
        return resolved.parent
    return resolved


@remote_app.command("install")
def remote_install_cmd(
    remote: str = typer.Option(
        ...,
        "--remote",
        help=f"SSH host alias ({', '.join(sorted(KNOWN_REMOTE_HOSTS))}).",
    ),
    source: Path = typer.Option(
        Path.cwd(),
        "--source",
        help="Local dft-learn checkout to rsync (defaults to current directory).",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
    ),
    remote_install_dir: str = typer.Option(
        "/home/hduva/projects/dft-learn",
        "--remote-install-dir",
        help="Directory on the remote host for the checkout and uv tool install.",
    ),
) -> None:
    """Install or refresh ``dftrun`` on a remote workstation over SSH.

    Rsyncs the local repository, runs ``uv tool install --force .`` on the host,
    and verifies ``dftrun --help``. Requires ``uv`` on the remote (installed
    automatically when missing).
    """
    try:
        host = resolve_ssh_host(remote)
    except KeyError as exc:
        typer.echo(str(exc), err=True)
        raise typer.Exit(1) from exc

    root = _project_root(source)
    try:
        install_dftrun_on_host(
            host,
            root,
            remote_install_dir=remote_install_dir,
        )
    except (ValueError, RuntimeError, OSError) as exc:
        typer.echo(f"remote install failed: {exc}", err=True)
        raise typer.Exit(1) from exc

    _CONSOLE.print(
        f"Installed dftrun on {host}\n"
        f"  checkout: {remote_install_dir}\n"
        f"  runs: {DEFAULT_REMOTE_ROOT}/<run-name>/"
    )


@remote_app.command("diagnose")
def remote_diagnose_cmd(
    run_directory: Path = typer.Argument(
        ...,
        help="Local run directory whose remote counterpart will be inspected.",
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
        help="Remote base directory containing the run folder.",
    ),
    pull: bool = typer.Option(
        True,
        "--pull/--no-pull",
        help="Pull logs and diagnostics after inspecting the remote run.",
    ),
) -> None:
    """Run ``dftrun run diagnose`` on a remote run and optionally pull the report."""
    try:
        host = resolve_ssh_host(remote)
    except KeyError as exc:
        typer.echo(str(exc), err=True)
        raise typer.Exit(1) from exc

    remote_dir = remote_run_directory(host, run_directory.name, remote_root)
    code = run_remote_diagnose(host, remote_dir)
    if pull:
        try:
            pull_run_artifacts(
                run_directory.resolve(), host=host, remote_dir=remote_dir
            )
        except (RuntimeError, OSError) as exc:
            typer.echo(f"pull failed: {exc}", err=True)
            raise typer.Exit(1) from exc
        report = run_directory.resolve() / "packaged_output" / "run_diagnostics.txt"
        if report.is_file():
            _CONSOLE.print(report.read_text(encoding="utf-8"))
            _CONSOLE.print(f"[green]Pulled[/green] {report}")
    raise typer.Exit(code)


@remote_app.command("reset")
def remote_reset_cmd(
    run_directory: Path = typer.Argument(
        ...,
        help="Local run directory whose remote counterpart will be reset.",
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
        help="Remote base directory containing the run folder.",
    ),
    yes: bool = typer.Option(
        False,
        "--yes",
        "-y",
        help="Confirm deletion of remote calculation outputs.",
    ),
) -> None:
    """Clear remote StoBe outputs while keeping run inputs."""
    if not yes:
        typer.echo("Refusing to reset without --yes.", err=True)
        raise typer.Exit(1)
    try:
        host = resolve_ssh_host(remote)
    except KeyError as exc:
        typer.echo(str(exc), err=True)
        raise typer.Exit(1) from exc

    remote_dir = remote_run_directory(host, run_directory.name, remote_root)
    code = run_remote_reset(host, remote_dir)
    if code == 0:
        _CONSOLE.print(f"Reset remote run at {host}:{remote_dir}/")
    raise typer.Exit(code)


def _print_remote_job_status(host: str, status: RemoteJobStatus) -> None:
    _CONSOLE.print(
        f"{status.state} on {host}  pid={status.pid}  exit={status.exit_code}\n"
        f"  log: {status.log_path}"
    )
    if status.log_tail:
        _CONSOLE.print(status.log_tail)


@remote_app.command("status")
def remote_status_cmd(
    run_directory: Path = typer.Argument(
        ...,
        help="Local run directory whose remote job will be inspected.",
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
        help="Remote base directory containing the run folder.",
    ),
) -> None:
    """Print whether a detached remote job is still running."""
    try:
        host = resolve_ssh_host(remote)
    except KeyError as exc:
        typer.echo(str(exc), err=True)
        raise typer.Exit(1) from exc

    remote_dir = remote_run_directory(host, run_directory.name, remote_root)
    status = inspect_remote_job(host, remote_dir)
    _print_remote_job_status(host, status)
    if status.state == "running":
        raise typer.Exit(2)
    if status.state == "finished" and status.exit_code not in (0, None):
        raise typer.Exit(int(status.exit_code))
    if status.state == "unknown":
        raise typer.Exit(1)


@remote_app.command("collect")
def remote_collect_cmd(
    run_directory: Path = typer.Argument(
        ...,
        help="Local run directory that receives logs and packaged_output.",
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
        help="Remote base directory containing the run folder.",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Pull artifacts even if the remote job is still running.",
    ),
) -> None:
    """Pull logs and packaged_output after a detached remote run finishes.

    If the job is still running, prints the latest log lines and exits 2.
    Re-run this command later (or after reconnecting) to transfer results.
    """
    try:
        host = resolve_ssh_host(remote)
    except KeyError as exc:
        typer.echo(str(exc), err=True)
        raise typer.Exit(1) from exc

    remote_dir = remote_run_directory(host, run_directory.name, remote_root)
    status = inspect_remote_job(host, remote_dir)
    if status.state == "running" and not force:
        _print_remote_job_status(host, status)
        raise typer.Exit(2)
    if status.state == "unknown" and not force:
        _print_remote_job_status(host, status)
        typer.echo("No finished remote job marker found.", err=True)
        raise typer.Exit(1)

    try:
        pull_run_artifacts(run_directory.resolve(), host=host, remote_dir=remote_dir)
    except (RuntimeError, OSError) as exc:
        typer.echo(f"collect failed: {exc}", err=True)
        raise typer.Exit(1) from exc

    packaged = run_directory.resolve() / "packaged_output"
    _CONSOLE.print(
        f"Pulled logs/ and packaged_output/ from {host}:{remote_dir}/\n  -> {packaged}/"
    )
    if status.state == "finished" and status.exit_code not in (0, None):
        raise typer.Exit(int(status.exit_code))
