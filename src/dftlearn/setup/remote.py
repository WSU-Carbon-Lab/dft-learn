"""SSH helpers for syncing initialized run directories to remote hosts."""

from __future__ import annotations

import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

KNOWN_REMOTE_HOSTS: dict[str, str] = {
    "hduva": "hduva",
    "hduva-workstation": "hduva-workstation",
}

DEFAULT_REMOTE_ROOT = "/home/hduva/projects/dft-runs"

REMOTE_ARTIFACT_DIRS: tuple[str, ...] = ("logs", "packaged_output")

RUN_SYNC_EXCLUDES: tuple[str, ...] = (
    "GND/",
    "EXC/",
    "TP/",
    "NEXAFS/",
    "logs/",
    "packaged_output/",
    "packaged_output.tar.gz",
)


def resolve_ssh_host(alias: str) -> str:
    """Map a user alias to an OpenSSH ``Host`` name from the local config.

    Parameters
    ----------
    alias
        Short name such as ``hduva``.

    Returns
    -------
    str
        Host entry for ``ssh`` and ``rsync``.

    Raises
    ------
    KeyError
        If ``alias`` is not a known remote preset.
    """
    key = alias.strip()
    if key not in KNOWN_REMOTE_HOSTS:
        msg = (
            f"Unknown remote host alias {alias!r}; "
            f"known: {', '.join(sorted(KNOWN_REMOTE_HOSTS))}"
        )
        raise KeyError(msg)
    return KNOWN_REMOTE_HOSTS[key]


def remote_run_directory(host: str, run_name: str, remote_root: str) -> Path:
    """Build the remote POSIX path for a run directory."""
    root = remote_root.rstrip("/")
    return Path(f"{root}/{run_name}")


def ssh_bash(
    host: str,
    command: str,
    *,
    tty: bool = False,
    check: bool = False,
    capture_output: bool = False,
) -> subprocess.CompletedProcess[Any]:
    """Run ``command`` on ``host`` in bash after sourcing ``~/.bashrc``.

    The workstation login shell is fish. OpenSSH runs the remote command
    through that shell, so the ``bash -c`` payload is passed as one quoted
    argument. Otherwise fish re-parses assignments such as ``pid=$!``.
    Invoking bash after ``source ~/.bashrc`` also sets ``LD_LIBRARY_PATH`` so
    ``StoBe.x`` can resolve ``libmkl_intel_lp64.so.2``.

    Parameters
    ----------
    host
        OpenSSH host alias (for example ``hduva``).
    command
        POSIX shell text executed after ``source ~/.bashrc``.
    tty
        When True, allocate a remote TTY (``ssh -t``) for progress output.
    check
        When True, raise ``CalledProcessError`` on a non-zero exit.
    capture_output
        When True, capture stdout and stderr as text.

    Returns
    -------
    subprocess.CompletedProcess
        Result of the ``ssh`` invocation.
    """
    wrapped = f'source "$HOME/.bashrc" >/dev/null 2>&1; {command}'
    argv = ["ssh"]
    if tty:
        argv.append("-t")
    argv.extend([host, f"bash -c {shlex.quote(wrapped)}"])
    return subprocess.run(
        argv,
        check=check,
        capture_output=capture_output,
        text=True,
    )


def sync_run_directory(
    local_dir: Path,
    *,
    host: str,
    remote_dir: Path,
) -> None:
    """Rsync a local run directory to a remote host over SSH.

    Creates ``remote_dir`` on the remote side, then synchronizes inputs with
    ``rsync -az``. Remote calculation trees (``GND/``, ``EXC/``, ``TP/``,
    ``NEXAFS/``, ``logs/``, ``packaged_output/``) are excluded so prior results
    are not deleted or overwritten by an empty local tree.

    Parameters
    ----------
    local_dir
        Initialized run directory on the local machine.
    host
        OpenSSH host alias (for example ``hduva``).
    remote_dir
        Absolute remote destination directory.
    """
    local_dir = Path(local_dir)
    remote = f"{host}:{remote_dir}/"
    ssh_bash(host, f"mkdir -p {shlex.quote(str(remote_dir))}", check=True)
    rsync_cmd = ["rsync", "-az"]
    for pattern in RUN_SYNC_EXCLUDES:
        rsync_cmd.extend(["--exclude", pattern])
    rsync_cmd.extend([f"{local_dir}/", remote])
    proc = subprocess.run(rsync_cmd, check=False, capture_output=True, text=True)
    if proc.returncode != 0:
        msg = f"rsync failed ({proc.returncode}): {proc.stderr.strip()}"
        raise RuntimeError(msg)


def _rsync_from_remote(
    host: str,
    remote_path: Path,
    local_path: Path,
) -> None:
    """Rsync one remote path into a local directory."""
    local_path = Path(local_path)
    local_path.mkdir(parents=True, exist_ok=True)
    remote = f"{host}:{remote_path}/"
    rsync_cmd = [
        "rsync",
        "-az",
        remote,
        f"{local_path}/",
    ]
    proc = subprocess.run(rsync_cmd, check=False, capture_output=True, text=True)
    if proc.returncode != 0:
        msg = f"rsync failed ({proc.returncode}): {proc.stderr.strip()}"
        raise RuntimeError(msg)


def pull_run_artifacts(
    local_dir: Path,
    *,
    host: str,
    remote_dir: Path,
) -> None:
    """Rsync run logs and packaged postprocess outputs from a remote host.

    Pulls ``logs/`` and ``packaged_output/`` only, not full StoBe calculation
    trees (``GND/``, ``EXC/``, ``TP/``, ``NEXAFS/``, site ``*.out`` files).

    Parameters
    ----------
    local_dir
        Local run directory that receives artifact subdirectories.
    host
        OpenSSH host alias (for example ``hduva``).
    remote_dir
        Absolute remote run directory.
    """
    local_dir = Path(local_dir).resolve()
    remote_dir = Path(remote_dir)
    local_dir.mkdir(parents=True, exist_ok=True)
    for name in REMOTE_ARTIFACT_DIRS:
        _rsync_from_remote(host, remote_dir / name, local_dir / name)


def pull_run_directory(
    local_dir: Path,
    *,
    host: str,
    remote_dir: Path,
) -> None:
    """Rsync a remote run directory back to the local machine over SSH.

    Creates ``local_dir`` when missing, then synchronizes remote contents with
    ``rsync -az``. Raises ``RuntimeError`` when ``rsync`` exits non-zero.

    Parameters
    ----------
    local_dir
        Local destination for the run directory.
    host
        OpenSSH host alias (for example ``hduva``).
    remote_dir
        Absolute remote source directory.
    """
    local_dir = Path(local_dir)
    local_dir.mkdir(parents=True, exist_ok=True)
    remote = f"{host}:{remote_dir}/"
    rsync_cmd = [
        "rsync",
        "-az",
        remote,
        f"{local_dir}/",
    ]
    proc = subprocess.run(rsync_cmd, check=False, capture_output=True, text=True)
    if proc.returncode != 0:
        msg = f"rsync failed ({proc.returncode}): {proc.stderr.strip()}"
        raise RuntimeError(msg)


def run_remote_dftrun(
    local_dir: Path,
    *,
    host: str,
    remote_root: str,
    subcommand: str,
    forward_args: list[str] | None = None,
    subcommand_args: list[str] | None = None,
    sync_first: bool = True,
    detach: bool = False,
    postprocess: bool = True,
) -> int:
    """Sync a run directory and execute ``dftrun run SUBCOMMAND`` on a remote host.

    Parameters
    ----------
    local_dir
        Local run directory whose name defines the remote folder under ``remote_root``.
    host
        OpenSSH host alias (for example ``hduva``).
    remote_root
        Remote base directory; created with ``mkdir -p`` when syncing.
    subcommand
        ``dftrun run`` subcommand such as ``all``, ``gnd``, or ``exc``.
    forward_args
        ``dftrun run`` callback flags placed before ``subcommand`` (``--workers``,
        ``--max-workers``, ``--verbose``).
    subcommand_args
        Flags for the run subcommand placed after the directory (``--atom``).
    sync_first
        When True, rsync ``local_dir`` to the remote before launching ``dftrun``.
    detach
        When True, launch the remote job with ``nohup`` and return immediately.
        Results stay on the host until ``dftrun remote collect``.
    postprocess
        When ``detach`` is True, also run ``dftrun postprocess .`` after the DFT
        sequence. Ignored for a blocking remote run (the caller postprocesses).

    Returns
    -------
    int
        Exit code from the remote ``ssh`` session.
    """
    local_dir = Path(local_dir).resolve()
    remote_dir = remote_run_directory(host, local_dir.name, remote_root)
    if sync_first:
        sync_run_directory(local_dir, host=host, remote_dir=remote_dir)

    inner = [
        "dftrun",
        "run",
        *(forward_args or []),
        subcommand,
        ".",
        *(subcommand_args or []),
    ]
    remote_shell = f"cd {shlex.quote(str(remote_dir))} && {shlex.join(inner)}"
    if detach:
        return _start_detached_remote_job(
            host,
            remote_dir,
            remote_shell,
            postprocess=postprocess,
        )
    proc = ssh_bash(host, remote_shell, tty=True)
    return int(proc.returncode)


@dataclass(frozen=True)
class RemoteJobStatus:
    """Reports whether a detached remote ``dftrun`` job is still running.

    Attributes
    ----------
    state
        ``running`` when the recorded PID is alive, ``finished`` when
        ``remote_job.exit`` exists, otherwise ``unknown``.
    pid
        Process id from ``remote_job.pid``, or None when that file is missing.
    exit_code
        Integer from ``remote_job.exit`` after the job ends, otherwise None.
    log_path
        Absolute remote path of ``logs/remote_job.log``.
    log_tail
        Last lines of that log, or an empty string when the file is absent.
    """

    state: str
    pid: int | None
    exit_code: int | None
    log_path: str
    log_tail: str


def _start_detached_remote_job(
    host: str,
    remote_dir: Path,
    run_shell: str,
    *,
    postprocess: bool,
) -> int:
    """Launch ``run_shell`` under nohup; write the PID to ``remote_job.pid``."""
    quoted = shlex.quote(str(remote_dir))
    post_cmd = "dftrun postprocess ." if postprocess else "true"
    job_body = (
        'source "$HOME/.bashrc" >/dev/null 2>&1; '
        'export PATH="$HOME/.local/bin:$PATH"; '
        f"cd {quoted} && mkdir -p logs && rm -f remote_job.exit && "
        f"(echo START && {run_shell} && rc=0 || rc=$?; "
        f'if [ "$rc" -eq 0 ]; then {post_cmd}; rc=$?; '
        f"else dftrun run diagnose . || true; fi; "
        f'echo END exit=$rc; echo "$rc" > remote_job.exit) '
        f"> logs/remote_job.log 2>&1"
    )
    launch = (
        f"nohup bash -c {shlex.quote(job_body)} </dev/null >/dev/null 2>&1 & "
        f'pid=$!; disown "$pid" 2>/dev/null || true; '
        f'echo "$pid" > {quoted}/remote_job.pid; echo "$pid"'
    )
    proc = ssh_bash(host, launch, capture_output=True)
    if proc.returncode != 0:
        detail = (proc.stderr or proc.stdout or "").strip() or "unknown error"
        msg = f"failed to start detached remote job: {detail}"
        raise RuntimeError(msg)
    pid_text = (proc.stdout or "").strip().splitlines()
    if not pid_text or not pid_text[-1].isdigit():
        msg = f"detached remote job did not report a PID: {proc.stdout!r}"
        raise RuntimeError(msg)
    return 0


def inspect_remote_job(host: str, remote_dir: Path) -> RemoteJobStatus:
    """Read pid/exit files and a log tail for a detached remote job.

    Parameters
    ----------
    host
        OpenSSH host alias.
    remote_dir
        Absolute remote run directory.

    Returns
    -------
    RemoteJobStatus
        ``running``, ``finished``, or ``unknown``.
    """
    quoted = shlex.quote(str(remote_dir))
    probe = (
        f"d={quoted}; "
        'pid=\'\'; [ -f "$d/remote_job.pid" ] && pid=$(cat "$d/remote_job.pid"); '
        'ex=\'\'; [ -f "$d/remote_job.exit" ] && ex=$(cat "$d/remote_job.exit"); '
        'if [ -n "$ex" ]; then echo "finished|$pid|$ex"; '
        'elif [ -n "$pid" ] && kill -0 "$pid" 2>/dev/null; '
        'then echo "running|$pid|"; '
        'else echo "unknown|$pid|$ex"; fi; '
        "echo '---LOG---'; "
        'tail -n 20 "$d/logs/remote_job.log" 2>/dev/null || true'
    )
    proc = ssh_bash(host, probe, capture_output=True)
    text = proc.stdout or ""
    header, _, log_tail = text.partition("---LOG---")
    parts = header.strip().split("|")
    state = parts[0] if parts else "unknown"
    pid = int(parts[1]) if len(parts) > 1 and parts[1].strip().isdigit() else None
    raw_exit = parts[2].strip() if len(parts) > 2 else ""
    exit_code = int(raw_exit) if raw_exit.lstrip("-").isdigit() else None
    return RemoteJobStatus(
        state=state,
        pid=pid,
        exit_code=exit_code,
        log_path=f"{remote_dir}/logs/remote_job.log",
        log_tail=log_tail.strip(),
    )


def remote_run_has_postprocess_inputs(host: str, remote_dir: Path) -> bool:
    """Return True when a remote run directory contains StoBe spectrum tables.

    Parameters
    ----------
    host
        OpenSSH host alias (for example ``hduva``).
    remote_dir
        Absolute remote run directory to inspect.

    Returns
    -------
    bool
        True when ``XrayT001.out`` or ``*xas.out`` files exist under the run tree.
    """
    quoted = shlex.quote(str(remote_dir))
    remote_shell = (
        f"d={quoted}; "
        "find \"$d\" -maxdepth 3 \\( -name 'XrayT001.out' -o -name '*xas.out' \\) "
        "-type f -print -quit | grep -q ."
    )
    proc = ssh_bash(host, remote_shell, capture_output=True)
    return proc.returncode == 0


def run_remote_postprocess(
    host: str,
    remote_dir: Path,
) -> int:
    """Run ``dftrun postprocess .`` in a remote run directory.

    Parameters
    ----------
    host
        OpenSSH host alias (for example ``hduva``).
    remote_dir
        Absolute remote run directory containing StoBe outputs.

    Returns
    -------
    int
        Exit code from the remote ``ssh`` session.
    """
    remote_shell = f"cd {shlex.quote(str(remote_dir))} && dftrun postprocess ."
    proc = ssh_bash(host, remote_shell, tty=True)
    return int(proc.returncode)


def run_remote_diagnose(host: str, remote_dir: Path) -> int:
    """Run ``dftrun run diagnose .`` in a remote run directory."""
    remote_shell = f"cd {shlex.quote(str(remote_dir))} && dftrun run diagnose ."
    proc = ssh_bash(host, remote_shell, tty=True)
    return int(proc.returncode)


def run_remote_reset(host: str, remote_dir: Path) -> int:
    """Run ``dftrun run reset . --yes`` in a remote run directory."""
    remote_shell = f"cd {shlex.quote(str(remote_dir))} && dftrun run reset . --yes"
    proc = ssh_bash(host, remote_shell, tty=True)
    return int(proc.returncode)


def install_dftrun_on_host(
    host: str,
    source_dir: Path,
    *,
    remote_install_dir: str = "/home/hduva/projects/dft-learn",
) -> None:
    """Rsync a checkout to ``host`` and install ``dftrun`` with ``uv tool``.

    Parameters
    ----------
    host
        OpenSSH host alias (for example ``hduva``).
    source_dir
        Local ``dft-learn`` repository root containing ``pyproject.toml``.
    remote_install_dir
        Directory on the remote host that receives the checkout.

    Raises
    ------
    ValueError
        If ``source_dir`` is not a project root.
    RuntimeError
        If bootstrap, rsync, or remote install commands fail.
    """
    source_dir = Path(source_dir).resolve()
    pyproject = source_dir / "pyproject.toml"
    if not pyproject.is_file():
        msg = f"Not a project root (missing pyproject.toml): {source_dir}"
        raise ValueError(msg)

    bootstrap = (
        "command -v uv >/dev/null 2>&1 || "
        "curl -LsSf https://astral.sh/uv/install.sh | sh"
    )
    ssh_bash(host, bootstrap, check=True)
    mkdir_remote = f"mkdir -p {shlex.quote(remote_install_dir)}"
    ssh_bash(host, mkdir_remote, check=True)

    remote_target = f"{host}:{remote_install_dir}/"
    rsync_cmd = [
        "rsync",
        "-az",
        "--delete",
        "--exclude",
        ".venv",
        "--exclude",
        "__pycache__",
        "--exclude",
        ".git",
        "--exclude",
        ".cursor",
        "--exclude",
        "streamlit_demo/.venv",
        f"{source_dir}/",
        remote_target,
    ]
    proc = subprocess.run(rsync_cmd, check=False, capture_output=True, text=True)
    if proc.returncode != 0:
        msg = f"rsync failed ({proc.returncode}): {proc.stderr.strip()}"
        raise RuntimeError(msg)

    install_shell = (
        'export PATH="$HOME/.local/bin:$PATH" && '
        f"cd {shlex.quote(remote_install_dir)} && "
        "uv tool install --force . && "
        "dftrun --help >/dev/null"
    )
    proc = ssh_bash(host, install_shell, capture_output=True)
    if proc.returncode != 0:
        detail = proc.stderr.strip() or proc.stdout.strip() or "unknown error"
        msg = f"remote dftrun install failed ({proc.returncode}): {detail}"
        raise RuntimeError(msg)
