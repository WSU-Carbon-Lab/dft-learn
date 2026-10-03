"""Tests for remote sync and run helpers."""

from __future__ import annotations

import shlex
import subprocess
from pathlib import Path

import pytest

from dftlearn.setup.remote import (
    pull_run_artifacts,
    remote_run_directory,
    resolve_ssh_host,
)


def _ssh_inner_command(argv: list[str]) -> str:
    tokens = shlex.split(argv[-1])
    assert tokens[:2] == ["bash", "-c"]
    return tokens[2]


def test_resolve_ssh_host_hduva() -> None:
    assert resolve_ssh_host("hduva") == "hduva"


def test_resolve_ssh_host_unknown() -> None:
    with pytest.raises(KeyError):
        resolve_ssh_host("not-a-host")


def test_remote_run_directory() -> None:
    path = remote_run_directory("hduva", "als3-001", "/home/hduva/projects/dft-runs")
    assert path == Path("/home/hduva/projects/dft-runs/als3-001")


def test_run_remote_dftrun_command_order() -> None:
    """Callback flags precede the subcommand; subcommand flags follow the directory."""
    import shlex

    from dftlearn.setup import remote as remote_mod
    from dftlearn.setup.remote import run_remote_dftrun

    local = Path("/tmp/als3-001")
    captured: list[list[str]] = []

    def fake_run(cmd, check=False, **_kwargs):
        captured.append(cmd)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    original = remote_mod.subprocess.run
    remote_mod.subprocess.run = fake_run  # type: ignore[assignment]
    try:
        run_remote_dftrun(
            local,
            host="hduva",
            remote_root="/home/hduva/projects/dft-runs",
            subcommand="all",
            forward_args=["--max-workers", "4"],
            subcommand_args=["--atom", "C1"],
            sync_first=False,
        )
    finally:
        remote_mod.subprocess.run = original

    argv = captured[0]
    assert argv[:3] == ["ssh", "-t", "hduva"]
    wrapped = _ssh_inner_command(argv)
    assert 'source "$HOME/.bashrc"' in wrapped
    tokens = shlex.split(wrapped.split("&&", maxsplit=1)[1].strip())
    assert tokens[:6] == ["dftrun", "run", "--max-workers", "4", "all", "."]
    assert tokens[6:] == ["--atom", "C1"]


def test_run_remote_dftrun_detach_uses_nohup() -> None:
    """Detached remote jobs launch under nohup without allocating a TTY."""
    from dftlearn.setup import remote as remote_mod
    from dftlearn.setup.remote import run_remote_dftrun

    captured: list[list[str]] = []

    def fake_run(cmd, check=False, **_kwargs):
        captured.append(cmd)
        return subprocess.CompletedProcess(cmd, 0, stdout="4321\n", stderr="")

    original = remote_mod.subprocess.run
    remote_mod.subprocess.run = fake_run  # type: ignore[assignment]
    try:
        code = run_remote_dftrun(
            Path("/tmp/als3-001"),
            host="hduva",
            remote_root="/home/hduva/projects/dft-runs",
            subcommand="all",
            forward_args=["--max-workers", "4"],
            sync_first=False,
            detach=True,
            postprocess=True,
        )
    finally:
        remote_mod.subprocess.run = original

    assert code == 0
    argv = captured[0]
    assert argv[:2] == ["ssh", "hduva"]
    wrapped = _ssh_inner_command(argv)
    assert "nohup bash -c" in wrapped
    assert "dftrun postprocess ." in wrapped
    assert "remote_job.pid" in wrapped
    assert "dftrun run --max-workers 4 all ." in wrapped


def test_inspect_remote_job_parses_running() -> None:
    from dftlearn.setup import remote as remote_mod
    from dftlearn.setup.remote import inspect_remote_job

    def fake_run(cmd, check=False, **_kwargs):
        stdout = "running|4321|\n---LOG---\nSTART\n"
        return subprocess.CompletedProcess(cmd, 0, stdout=stdout, stderr="")

    original = remote_mod.subprocess.run
    remote_mod.subprocess.run = fake_run  # type: ignore[assignment]
    remote_dir = Path("/home/hduva/projects/dft-runs/als3-001")
    try:
        status = inspect_remote_job("hduva", remote_dir)
    finally:
        remote_mod.subprocess.run = original

    assert status.state == "running"
    assert status.pid == 4321
    assert status.exit_code is None
    assert status.log_tail == "START"


def test_inspect_remote_job_parses_finished() -> None:
    from dftlearn.setup import remote as remote_mod
    from dftlearn.setup.remote import inspect_remote_job

    def fake_run(cmd, check=False, **_kwargs):
        stdout = "finished|4321|0\n---LOG---\nEND exit=0\n"
        return subprocess.CompletedProcess(cmd, 0, stdout=stdout, stderr="")

    original = remote_mod.subprocess.run
    remote_mod.subprocess.run = fake_run  # type: ignore[assignment]
    remote_dir = Path("/home/hduva/projects/dft-runs/als3-001")
    try:
        status = inspect_remote_job("hduva", remote_dir)
    finally:
        remote_mod.subprocess.run = original

    assert status.state == "finished"
    assert status.pid == 4321
    assert status.exit_code == 0
    assert "END exit=0" in status.log_tail


def test_pull_run_artifacts_rsync_args() -> None:
    """Artifact pull rsyncs logs/ and packaged_output/ only."""
    from dftlearn.setup import remote as remote_mod

    local = Path("/tmp/als3-001")
    local_resolved = local.resolve()
    captured: list[list[str]] = []

    def fake_run(cmd, check=False, capture_output=False, text=False):
        captured.append(cmd)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    original = remote_mod.subprocess.run
    remote_mod.subprocess.run = fake_run  # type: ignore[assignment]
    try:
        pull_run_artifacts(
            local,
            host="hduva",
            remote_dir=Path("/home/hduva/projects/dft-runs/als3-001"),
        )
    finally:
        remote_mod.subprocess.run = original

    assert len(captured) == 2
    remote_base = "/home/hduva/projects/dft-runs/als3-001"
    assert captured[0][2] == f"hduva:{remote_base}/logs/"
    assert captured[0][3] == f"{local_resolved}/logs/"
    assert captured[1][2] == f"hduva:{remote_base}/packaged_output/"
    assert captured[1][3] == f"{local_resolved}/packaged_output/"


def test_run_remote_postprocess_command() -> None:
    """Remote postprocess runs dftrun postprocess in the run directory."""
    import shlex

    from dftlearn.setup import remote as remote_mod
    from dftlearn.setup.remote import run_remote_postprocess

    captured: list[list[str]] = []

    def fake_run(cmd, check=False, **_kwargs):
        captured.append(cmd)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    original = remote_mod.subprocess.run
    remote_mod.subprocess.run = fake_run  # type: ignore[assignment]
    try:
        run_remote_postprocess(
            "hduva",
            Path("/home/hduva/projects/dft-runs/als3-001"),
        )
    finally:
        remote_mod.subprocess.run = original

    argv = captured[0]
    assert argv[:3] == ["ssh", "-t", "hduva"]
    wrapped = _ssh_inner_command(argv)
    assert "dftrun postprocess ." in wrapped
    tokens = shlex.split(wrapped.split("&&", maxsplit=1)[1].strip())
    assert tokens == ["dftrun", "postprocess", "."]


def test_sync_run_directory_excludes_remote_outputs() -> None:
    """Upload sync excludes calculation trees so remote results are preserved."""
    from dftlearn.setup import remote as remote_mod
    from dftlearn.setup.remote import RUN_SYNC_EXCLUDES, sync_run_directory

    captured: list[list[str]] = []

    def fake_run(cmd, check=False, capture_output=False, text=False):
        captured.append(cmd)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    original = remote_mod.subprocess.run
    remote_mod.subprocess.run = fake_run  # type: ignore[assignment]
    try:
        sync_run_directory(
            Path("/tmp/als3-001"),
            host="hduva",
            remote_dir=Path("/home/hduva/projects/dft-runs/als3-001"),
        )
    finally:
        remote_mod.subprocess.run = original

    rsync_cmd = captured[1]
    assert "--delete" not in rsync_cmd
    for pattern in RUN_SYNC_EXCLUDES:
        assert "--exclude" in rsync_cmd
        assert pattern in rsync_cmd


def test_remote_run_has_postprocess_inputs() -> None:
    from dftlearn.setup import remote as remote_mod
    from dftlearn.setup.remote import remote_run_has_postprocess_inputs

    def fake_run(cmd, check=False, **_kwargs):
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    original = remote_mod.subprocess.run
    remote_mod.subprocess.run = fake_run  # type: ignore[assignment]
    try:
        assert remote_run_has_postprocess_inputs(
            "hduva",
            Path("/home/hduva/projects/dft-runs/als3-001"),
        )
    finally:
        remote_mod.subprocess.run = original


def test_ssh_bash_wraps_bashrc() -> None:
    """Remote commands run under bash after sourcing ~/.bashrc for MKL."""
    from dftlearn.setup import remote as remote_mod
    from dftlearn.setup.remote import ssh_bash

    captured: list[list[str]] = []

    def fake_run(cmd, check=False, capture_output=False, text=False):
        captured.append(cmd)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    original = remote_mod.subprocess.run
    remote_mod.subprocess.run = fake_run  # type: ignore[assignment]
    try:
        ssh_bash("hduva", "dftrun --help", tty=True)
    finally:
        remote_mod.subprocess.run = original

    assert captured[0][:3] == ["ssh", "-t", "hduva"]
    inner = _ssh_inner_command(captured[0])
    assert inner.startswith('source "$HOME/.bashrc"')
    assert inner.endswith("dftrun --help")


def test_project_root_finds_repo() -> None:
    from dftlearn.cli.remote_cmd import _project_root

    root = _project_root(Path(__file__).resolve().parents[1])
    assert (root / "pyproject.toml").is_file()
    assert (root / "src" / "dftlearn").is_dir()
