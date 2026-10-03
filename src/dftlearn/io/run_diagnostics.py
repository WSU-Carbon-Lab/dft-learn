"""Inspect StoBe run directories for toolchain gaps, logs, and missing outputs."""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import TypedDict

from natsort import natsorted

from dftlearn.setup.basis_catalog import validate_molconfig_file


class StoBeOutputSniff(TypedDict):
    """Quick existence and SCF markers for one StoBe ``.out`` file."""

    exists: bool
    size: int
    scf_converged: bool
    has_final_energy: bool
    error_line: str | None

_CALC_TO_DIR = {"gnd": "GND", "exc": "EXC", "tp": "TP"}
_SITE_DIR_PATTERN = re.compile(r"^[A-Z]+\d+$")
_RESET_TREE_DIRS = ("GND", "EXC", "TP", "NEXAFS", "logs", "packaged_output")
_RESET_ROOT_FILES = ("packaged_output.tar.gz",)
_SITE_GLOBS = ("fort.*", "*.out", "*.err", "*.xas", "XrayT001.out", "Molden.molf")
_SUMMARY_PATTERN = re.compile(r"^.+_\d{14}\.txt$")
_ERROR_MARKERS = (
    "error while loading shared libraries",
    "cannot open shared object file",
    "No such file or directory",
    "command not found",
    "STOP",
    "ABORT",
    "Segmentation fault",
    "cannot stat",
)


@dataclass
class ToolchainCheck:
    """One StoBe toolchain probe."""

    name: str
    ok: bool
    detail: str


@dataclass
class JobDiagnostic:
    """Status for one ``*.run`` job under a run directory."""

    site: str
    calc_type: str
    run_file: str
    expected_output: str | None
    output_exists: bool
    output_bytes: int
    bash_returncode: int | None
    scf_converged: bool
    has_final_energy: bool
    log_excerpt: str
    issues: list[str] = field(default_factory=list)


@dataclass
class RunDiagnosticReport:
    """Aggregate inspection results for a StoBe run root."""

    run_root: str
    toolchain: list[ToolchainCheck]
    jobs: list[JobDiagnostic]
    summary: str
    ready_for_postprocess: bool


def expected_stobe_output_path(run_file: Path) -> Path | None:
    """Map a ``SITEcalc.run`` file to its primary StoBe output path.

    Parameters
    ----------
    run_file
        Path such as ``C1/C1gnd.run``.

    Returns
    -------
    pathlib.Path | None
        Expected ``GND/SITEgnd.out``-style path, or ``NEXAFS/SITExas.out`` for XAS.
    """
    run_file = Path(run_file)
    stem = run_file.stem
    if len(stem) <= 3:
        return None
    calc_type = stem[-3:]
    site = stem[:-3]
    run_root = run_file.parent.parent
    if calc_type in _CALC_TO_DIR:
        category = _CALC_TO_DIR[calc_type]
        return run_root / category / f"{site}{calc_type}.out"
    if calc_type == "xas":
        return run_root / "NEXAFS" / f"{site}xas.out"
    return None


def _site_directories(run_root: Path) -> list[Path]:
    return natsorted(
        [
            child
            for child in run_root.iterdir()
            if child.is_dir() and _SITE_DIR_PATTERN.fullmatch(child.name)
        ],
        key=lambda path: path.name,
    )


def _read_log_excerpt(run_file: Path, logs_dir: Path, max_lines: int = 12) -> str:
    site = run_file.parent.name
    calc = run_file.stem[-3:] if len(run_file.stem) > 3 else "run"
    if logs_dir.is_dir():
        matches = sorted(logs_dir.glob(f"{site}_{calc}_*.txt"), reverse=True)
        if matches:
            text = matches[0].read_text(encoding="utf-8", errors="replace")
            lines = text.splitlines()
            return "\n".join(lines[:max_lines])
    err_path = run_file.parent / f"{run_file.name}.err"
    if err_path.is_file():
        lines = err_path.read_text(encoding="utf-8", errors="replace").splitlines()
        return "\n".join(lines[:max_lines])
    out_path = run_file.parent / f"{run_file.name}.out"
    if out_path.is_file():
        lines = out_path.read_text(encoding="utf-8", errors="replace").splitlines()
        return "\n".join(lines[:max_lines])
    return ""


def sniff_stobe_output(path: Path) -> StoBeOutputSniff:
    """Summarize one StoBe ``.out`` file for quick failure detection."""
    path = Path(path)
    if not path.is_file():
        return StoBeOutputSniff(
            exists=False,
            size=0,
            scf_converged=False,
            has_final_energy=False,
            error_line=None,
        )
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    error_line: str | None = None
    for line in lines[:20]:
        lower = line.lower()
        if any(marker.lower() in lower for marker in _ERROR_MARKERS):
            error_line = line
            break
    return StoBeOutputSniff(
        exists=True,
        size=path.stat().st_size,
        scf_converged="SCF CONVERGED" in text,
        has_final_energy="FINAL ENERGY" in text,
        error_line=error_line,
    )


def check_local_toolchain() -> list[ToolchainCheck]:
    """Probe StoBe executables and basis paths on the current machine."""
    checks: list[ToolchainCheck] = []
    stobe = shutil.which("StoBe.x")
    checks.append(
        ToolchainCheck(
            "StoBe.x on PATH",
            stobe is not None,
            stobe or "not found",
        )
    )
    xas = shutil.which("xrayspec.x")
    checks.append(
        ToolchainCheck(
            "xrayspec.x on PATH",
            xas is not None,
            xas or "not found",
        )
    )
    stobe_home = os.environ.get("STOBE", "")
    checks.append(
        ToolchainCheck(
            "STOBE environment variable",
            bool(stobe_home),
            stobe_home or "unset",
        )
    )
    basis = Path("/bin/stobe/Basis/baslib.new7")
    checks.append(
        ToolchainCheck(
            "basis baslib.new7",
            basis.is_file(),
            str(basis),
        )
    )
    if stobe is not None:
        ldd = subprocess.run(
            ["ldd", stobe],
            check=False,
            capture_output=True,
            text=True,
        )
        missing = [
            line.strip()
            for line in (ldd.stdout or "").splitlines()
            if "not found" in line
        ]
        checks.append(
            ToolchainCheck(
                "StoBe.x shared libraries",
                not missing,
                "; ".join(missing) if missing else "resolved",
            )
        )
    return checks


def check_molconfig_basis(run_root: Path) -> ToolchainCheck:
    """Validate ``molConfig.py`` basis names against the bundled baslib catalog."""
    molconfig = Path(run_root) / "molConfig.py"
    if not molconfig.is_file():
        return ToolchainCheck(
            "molConfig basis names",
            True,
            "molConfig.py not present",
        )
    missing = validate_molconfig_file(molconfig)
    if missing:
        return ToolchainCheck(
            "molConfig basis names",
            False,
            "missing from baslib: " + ", ".join(missing),
        )
    live_baslib = Path("/bin/stobe/Basis/baslib.new7")
    if live_baslib.is_file():
        from dftlearn.setup.basis_catalog import (
            basis_strings_from_molconfig_source,
            validate_against_baslib_paths,
        )

        names = basis_strings_from_molconfig_source(
            molconfig.read_text(encoding="utf-8")
        )
        live_missing = validate_against_baslib_paths(names, baslib=live_baslib)
        if live_missing:
            return ToolchainCheck(
                "molConfig basis names (live baslib)",
                False,
                "missing on host: " + ", ".join(live_missing),
            )
    return ToolchainCheck("molConfig basis names", True, "catalog ok")


def collect_job_diagnostics(run_root: Path) -> list[JobDiagnostic]:
    """Inspect every ``*gnd|exc|tp|xas.run`` file under site directories."""
    run_root = Path(run_root).resolve()
    logs_dir = run_root / "logs"
    jobs: list[JobDiagnostic] = []
    for site_dir in _site_directories(run_root):
        for run_file in sorted(site_dir.glob("*.run")):
            calc_type = run_file.stem[-3:] if len(run_file.stem) > 3 else "run"
            if calc_type not in (*_CALC_TO_DIR.keys(), "xas"):
                continue
            expected = expected_stobe_output_path(run_file)
            expected_str = str(expected) if expected is not None else None
            sniff = StoBeOutputSniff(
                exists=False,
                size=0,
                scf_converged=False,
                has_final_energy=False,
                error_line=None,
            )
            if expected is not None:
                sniff = sniff_stobe_output(expected)
            log_excerpt = _read_log_excerpt(run_file, logs_dir)
            issues: list[str] = []
            if expected is not None and not bool(sniff["exists"]):
                rel = expected.relative_to(run_root)
                issues.append(f"missing expected output {rel}")
            elif expected is not None and int(sniff["size"]) < 200:
                issues.append(f"output too small ({sniff['size']} bytes)")
            error_line = sniff.get("error_line")
            if isinstance(error_line, str) and error_line:
                issues.append(error_line)
            scf_expected = calc_type in _CALC_TO_DIR and bool(sniff["exists"])
            if scf_expected and not bool(sniff["scf_converged"]):
                issues.append("SCF did not converge")
            if log_excerpt:
                for marker in _ERROR_MARKERS:
                    if marker.lower() in log_excerpt.lower():
                        issues.append(f"log: {log_excerpt.splitlines()[0]}")
                        break
            jobs.append(
                JobDiagnostic(
                    site=site_dir.name,
                    calc_type=calc_type,
                    run_file=str(run_file.relative_to(run_root)),
                    expected_output=expected_str,
                    output_exists=bool(sniff["exists"]),
                    output_bytes=int(sniff["size"]),
                    bash_returncode=None,
                    scf_converged=bool(sniff["scf_converged"]),
                    has_final_energy=bool(sniff["has_final_energy"]),
                    log_excerpt=log_excerpt,
                    issues=issues,
                )
            )
    return jobs


def build_run_diagnostic_report(run_root: Path) -> RunDiagnosticReport:
    """Build a structured diagnostic report for one run directory."""
    run_root = Path(run_root).resolve()
    toolchain = check_local_toolchain()
    toolchain.append(check_molconfig_basis(run_root))
    jobs = collect_job_diagnostics(run_root)
    bad_jobs = [job for job in jobs if job.issues]
    toolchain_bad = [check for check in toolchain if not check.ok]
    spectrum_ready = any((run_root / "NEXAFS").glob("*xas.out")) or any(
        (site / "XrayT001.out").is_file() for site in _site_directories(run_root)
    )
    lines = [
        f"Run root: {run_root}",
        f"Jobs inspected: {len(jobs)}",
        f"Jobs with issues: {len(bad_jobs)}",
        f"Toolchain checks failed: {len(toolchain_bad)}",
        f"Ready for postprocess: {spectrum_ready}",
    ]
    if toolchain_bad:
        lines.append("")
        lines.append("Toolchain:")
        for check in toolchain_bad:
            lines.append(f"  FAIL {check.name}: {check.detail}")
    if bad_jobs:
        lines.append("")
        lines.append("Sample job issues:")
        for job in bad_jobs[:6]:
            issue = job.issues[0] if job.issues else "unknown"
            lines.append(f"  {job.site} {job.calc_type}: {issue}")
    summary = "\n".join(lines)
    return RunDiagnosticReport(
        run_root=str(run_root),
        toolchain=toolchain,
        jobs=jobs,
        summary=summary,
        ready_for_postprocess=spectrum_ready,
    )


def format_run_diagnostic_report(report: RunDiagnosticReport) -> str:
    """Render a human-readable diagnostic report."""
    lines = [report.summary, ""]
    if report.toolchain:
        lines.append("Toolchain checks:")
        for check in report.toolchain:
            status = "ok" if check.ok else "FAIL"
            lines.append(f"  [{status}] {check.name}: {check.detail}")
        lines.append("")
    if report.jobs:
        lines.append("Jobs:")
        for job in report.jobs:
            flag = "OK" if not job.issues else "FAIL"
            lines.append(
                f"  [{flag}] {job.site} {job.calc_type} "
                f"output={job.output_bytes}B issues={len(job.issues)}"
            )
            if job.issues:
                for issue in job.issues[:2]:
                    lines.append(f"         - {issue}")
            if job.log_excerpt:
                first_log = job.log_excerpt.splitlines()[0]
                lines.append(f"         log: {first_log}")
    return "\n".join(lines)


def write_run_diagnostic_report(
    run_root: Path,
    *,
    out_dir: Path | None = None,
) -> tuple[Path, Path, RunDiagnosticReport]:
    """Write text and JSON diagnostics under ``packaged_output/``.

    Parameters
    ----------
    run_root
        StoBe run directory to inspect.
    out_dir
        Output folder; defaults to ``run_root/packaged_output``.

    Returns
    -------
    text_path : pathlib.Path
        Human-readable report path.
    json_path : pathlib.Path
        Machine-readable report path.
    report : RunDiagnosticReport
        Structured report object.
    """
    run_root = Path(run_root).resolve()
    packaged = Path(out_dir).resolve() if out_dir else (run_root / "packaged_output")
    packaged.mkdir(parents=True, exist_ok=True)
    report = build_run_diagnostic_report(run_root)
    text = format_run_diagnostic_report(report)
    text_path = packaged / "run_diagnostics.txt"
    json_path = packaged / "run_diagnostics.json"
    text_path.write_text(text + "\n", encoding="utf-8")
    json_path.write_text(
        json.dumps(asdict(report), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return text_path, json_path, report


def reset_run_directory(run_root: Path) -> list[str]:
    """Remove calculation outputs and packaged artifacts from a run directory.

    Preserves inputs such as ``molConfig.py``, ``geometry.xyz``, site ``*.run``
    scripts, and setup session files.

    Parameters
    ----------
    run_root
        StoBe run root to reset.

    Returns
    -------
    list[str]
        Short descriptions of removed paths.
    """
    run_root = Path(run_root).resolve()
    removed: list[str] = []
    for name in _RESET_TREE_DIRS:
        target = run_root / name
        if target.exists():
            shutil.rmtree(target)
            removed.append(str(target.relative_to(run_root)) + "/")
    for name in _RESET_ROOT_FILES:
        target = run_root / name
        if target.is_file():
            target.unlink()
            removed.append(name)
    for summary in run_root.glob("*.txt"):
        if _SUMMARY_PATTERN.fullmatch(summary.name):
            summary.unlink()
            removed.append(summary.name)
    for site_dir in _site_directories(run_root):
        for pattern in _SITE_GLOBS:
            for path in site_dir.glob(pattern):
                if path.is_file():
                    path.unlink()
                    removed.append(str(path.relative_to(run_root)))
    for path in run_root.rglob("fort.*"):
        if path.is_file():
            path.unlink()
            removed.append(str(path.relative_to(run_root)))
    return removed
