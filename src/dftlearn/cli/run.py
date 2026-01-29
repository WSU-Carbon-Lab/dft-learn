"""dftrun run: schedule and execute StoBe DFT calculations."""

from __future__ import annotations

import atexit
import multiprocessing
import os
import queue
import shutil
import signal
import subprocess
import sys
import tarfile
import threading
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from threading import Lock, Thread

import typer
from rich.console import Console
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)

from dftlearn.cli.config import load_config

_CONSOLE = Console()


def _resolve_logs_dirs(
    molecule_dir: Path, molecule_name: str, cfg: dict
) -> tuple[Path, Path]:
    log_dir_cfg = cfg.get("log_directory", "directory")
    if log_dir_cfg == "directory":
        logs_dir = molecule_dir / "logs"
        summary_dir = molecule_dir
    else:
        base = Path(log_dir_cfg).expanduser().resolve()
        logs_dir = base / molecule_name
        summary_dir = logs_dir
    return (logs_dir, summary_dir)


def _count_atoms(
    scan_dir: str, calc_types: list[str], atoms: list[str] | None = None
) -> int:
    scan_path = Path(scan_dir)
    directories: set[str] = set()
    for calc_type in calc_types:
        if atoms:
            for a in atoms:
                for run_file in scan_path.glob(f"**/{a}{calc_type}.run"):
                    directories.add(run_file.parent.name)
        else:
            for run_file in scan_path.glob(f"**/*{calc_type}.run"):
                directories.add(run_file.parent.name)
    return len(directories)


def _calculate_optimal_workers(
    scan_dir: str,
    calc_types: list[str],
    atoms: list[str] | None = None,
    max_workers: int | None = None,
) -> int:
    num = _count_atoms(scan_dir, calc_types, atoms)
    cap = (
        max_workers
        if max_workers is not None
        else max(1, multiprocessing.cpu_count() - 1)
    )
    return min(num, cap) if num else 0


@dataclass
class DirectoryJob:
    directory: str
    run_file: Path
    calc_type: str
    status: str = "pending"
    start_time: datetime | None = None
    end_time: datetime | None = None
    duration: float | None = None
    error_msg: str | None = None
    thread: Thread | None = None


class StoBeJobManager:
    def __init__(self, working_dir: Path | None = None) -> None:
        self.working_dir = working_dir or Path.cwd()
        self.jobs: dict = {}
        self.completed_jobs: dict = {}
        self.job_queue: queue.Queue = queue.Queue()
        self.max_concurrent_jobs = 4

    def find_run_files(self, pattern: str = "**/*.run") -> list[Path]:
        out = list(self.working_dir.glob(pattern))
        return sorted(out)

    def execute_run_file(
        self,
        run_file: Path,
        background: bool = True,
        log_path: Path | None = None,
    ) -> dict:
        if not run_file.exists():
            raise FileNotFoundError(f"Run file not found: {run_file}")
        run_dir = run_file.parent
        run_script = run_file.name
        job_info: dict = {
            "run_file": run_file,
            "run_dir": run_dir,
            "start_time": datetime.now(),
            "status": "starting",
            "process": None,
            "output_file": run_dir / f"{run_script}.out",
            "error_file": run_dir / f"{run_script}.err",
        }
        try:
            if background:
                out_f = job_info["output_file"].open("w")
                err_f = job_info["error_file"].open("w")
                proc = subprocess.Popen(
                    ["bash", run_script],
                    cwd=run_dir,
                    stdout=out_f,
                    stderr=err_f,
                    shell=False,
                )
                job_info["process"] = proc
                job_info["status"] = "running"
                job_info["pid"] = proc.pid
            else:
                if log_path is not None:
                    with log_path.open("w") as lf:
                        result = subprocess.run(
                            ["bash", run_script],
                            cwd=run_dir,
                            stdout=lf,
                            stderr=subprocess.STDOUT,
                            text=True,
                        )
                else:
                    result = subprocess.run(
                        ["bash", run_script],
                        cwd=run_dir,
                        capture_output=True,
                        text=True,
                    )
                job_info["returncode"] = result.returncode
                job_info["stdout"] = getattr(result, "stdout", None)
                job_info["stderr"] = getattr(result, "stderr", None)
                job_info["status"] = "completed" if result.returncode == 0 else "failed"
                job_info["end_time"] = datetime.now()
                if result.returncode != 0:
                    job_info["error"] = (
                        getattr(result, "stderr", None)
                        or getattr(result, "stdout", None)
                        or "Unknown error"
                    )
        except Exception as e:
            job_info["status"] = "error"
            job_info["error"] = str(e)
            job_info["end_time"] = datetime.now()
        return job_info


class TyperSchedulerRunner:
    def __init__(self) -> None:
        self.jobs: dict[str, DirectoryJob] = {}
        self.lock = Lock()
        self.completed_count = 0
        self.total_count = 0
        self.start_time: datetime | None = None
        self.console = Console()
        self.logs_dir: Path | None = None
        self.run_timestamp: str = ""
        self.summary_dir: Path | None = None
        self.molecule_name: str = ""
        self.molecule_dir: Path | None = None

    def find_run_files_in_directory(
        self,
        scan_dir: str,
        calc_types: list[str] | None = None,
        atoms: list[str] | None = None,
    ) -> list[Path]:
        scan_path = Path(scan_dir)
        if not scan_path.exists():
            return []
        if calc_types is None:
            calc_types = ["gnd", "exc", "tp", "xas", "seq"]
        patterns: list[str] = []
        if atoms is None:
            for ct in calc_types:
                patterns.append(f"**/*{ct}.run")
        else:
            for a in atoms:
                for ct in calc_types:
                    patterns.append(f"**/{a}{ct}.run")
        run_files: list[Path] = []
        for p in patterns:
            run_files.extend(scan_path.glob(p))
        return sorted(set(run_files))

    def find_directories_with_run_files(
        self,
        scan_dir: str,
        calc_types: list[str],
        atoms: list[str] | None = None,
    ) -> dict[str, list[Path]]:
        scan_path = Path(scan_dir)
        directories: dict[str, list[Path]] = {}
        for calc_type in calc_types:
            if atoms:
                for a in atoms:
                    for run_file in scan_path.glob(f"**/{a}{calc_type}.run"):
                        dn = run_file.parent.name
                        if dn not in directories:
                            directories[dn] = []
                        directories[dn].append(run_file)
            else:
                for run_file in scan_path.glob(f"**/*{calc_type}.run"):
                    dn = run_file.parent.name
                    if dn not in directories:
                        directories[dn] = []
                    directories[dn].append(run_file)
        return directories

    def execute_directory_job(self, job: DirectoryJob) -> DirectoryJob:
        job.status = "running"
        job.start_time = datetime.now()
        done = threading.Event()
        result: dict = {"job_info": None, "exception": None}

        log_path: Path | None = None
        if self.logs_dir is not None and self.run_timestamp:
            log_path = (
                self.logs_dir
                / f"{job.directory}_{job.calc_type}_{self.run_timestamp}.txt"
            )

        def run_calc() -> None:
            try:
                mgr = StoBeJobManager()
                result["job_info"] = mgr.execute_run_file(
                    job.run_file, background=False, log_path=log_path
                )
            except Exception as e:
                result["exception"] = e
            finally:
                done.set()

        t = threading.Thread(target=run_calc, daemon=True)
        t.start()
        while not done.is_set():
            time.sleep(1)
        t.join()
        job.end_time = datetime.now()
        job.duration = (
            (job.end_time - job.start_time).total_seconds() if job.start_time else None
        )
        if result["exception"]:
            job.status = "failed"
            job.error_msg = str(result["exception"])
        elif result["job_info"]:
            ji = result["job_info"]
            job.status = "completed" if ji.get("status") == "completed" else "failed"
            job.error_msg = ji.get("error")
        else:
            job.status = "failed"
            job.error_msg = "No calculation result received"
        return job

    def run_calculations_with_progress(
        self,
        scan_dir: str,
        calc_types: list[str],
        max_workers: int = 4,
        atoms: list[str] | None = None,
        override_logs_dir: Path | None = None,
        override_run_timestamp: str | None = None,
        override_summary_dir: Path | None = None,
        override_molecule_name: str | None = None,
        override_molecule_dir: Path | None = None,
        write_summary: bool = True,
    ) -> list[dict]:
        scan_path = Path(scan_dir).resolve()
        molecule_name = override_molecule_name or scan_path.name
        molecule_dir = override_molecule_dir or scan_path
        if override_logs_dir is not None and override_run_timestamp is not None:
            self.logs_dir = override_logs_dir
            self.run_timestamp = override_run_timestamp
            self.summary_dir = override_summary_dir or override_logs_dir
            self.molecule_name = molecule_name
            self.molecule_dir = molecule_dir
        else:
            cfg = load_config()
            logs_dir, summary_dir = _resolve_logs_dirs(molecule_dir, molecule_name, cfg)
            self.logs_dir = logs_dir
            self.run_timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
            self.summary_dir = summary_dir
            self.molecule_name = molecule_name
            self.molecule_dir = molecule_dir
            logs_dir.mkdir(parents=True, exist_ok=True)

        dirs = self.find_directories_with_run_files(scan_dir, calc_types, atoms)
        if not dirs:
            return []
        for dir_name, run_files in dirs.items():
            for run_file in run_files:
                ct = run_file.stem[-3:] if len(run_file.stem) > 3 else "unknown"
                job = DirectoryJob(directory=dir_name, run_file=run_file, calc_type=ct)
                self.jobs[f"{dir_name}:{run_file.name}"] = job
        self.total_count = len(self.jobs)
        self.completed_count = 0
        self.start_time = datetime.now()

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            TextColumn("-"),
            TimeElapsedColumn(),
            console=self.console,
        ) as progress:
            main_task = progress.add_task(
                f"[cyan]Processing {self.total_count} calculations ({', '.join(calc_types)})",
                total=self.total_count,
            )

            def worker(j: DirectoryJob) -> None:
                res = self.execute_directory_job(j)
                key = f"{j.directory}:{j.run_file.name}"
                with self.lock:
                    self.jobs[key] = res
                    if res.status in ("completed", "failed"):
                        self.completed_count += 1
                        desc = "completed" if res.status == "completed" else "failed"
                        progress.update(
                            main_task,
                            completed=self.completed_count,
                            description=f"[cyan]Processing {self.total_count} calculations | {desc} {j.directory} {j.calc_type} ({res.duration or 0:.1f}s)",
                        )

            active: list[Thread] = []
            items = list(self.jobs.items())
            idx = 0
            while idx < len(items) or active:
                while len(active) < max_workers and idx < len(items):
                    _, j = items[idx]
                    th = Thread(target=worker, args=(j,), daemon=True)
                    j.thread = th
                    th.start()
                    active.append(th)
                    idx += 1
                    time.sleep(0.1)
                active = [t for t in active if t.is_alive()]
                time.sleep(0.5)

        out: list[dict] = []
        for j in self.jobs.values():
            out.append(
                {
                    "directory": j.directory,
                    "calc_type": j.calc_type,
                    "status": j.status,
                    "duration": j.duration,
                    "error": j.error_msg,
                    "start_time": j.start_time,
                    "end_time": j.end_time,
                }
            )
        if (
            write_summary
            and self.summary_dir is not None
            and self.run_timestamp
            and self.molecule_name
        ):
            summary_path = (
                self.summary_dir / f"{self.molecule_name}_{self.run_timestamp}.txt"
            )
            self.summary_dir.mkdir(parents=True, exist_ok=True)
            with summary_path.open("w") as sf:
                sf.write("atom\tcalc_type\tstart_time\tend_time\tduration_s\tstatus\n")
                for j in self.jobs.values():
                    st = j.start_time.isoformat() if j.start_time else ""
                    et = j.end_time.isoformat() if j.end_time else ""
                    dur = f"{j.duration:.2f}" if j.duration is not None else ""
                    sf.write(
                        f"{j.directory}\t{j.calc_type}\t{st}\t{et}\t{dur}\t{j.status}\n"
                    )
        return out


def _cleanup_fort_files(scan_dir: Path) -> None:
    root = Path(scan_dir).resolve()
    for f in root.rglob("fort.*"):
        if f.is_file():
            try:
                f.unlink()
            except OSError:
                pass


def _argv_without_quiet_and_log() -> list[str]:
    skip = {"--quiet", "-q", "--log-file", "-l"}
    out: list[str] = []
    i = 0
    argv = sys.argv[1:]
    while i < len(argv):
        a = argv[i]
        if a in skip:
            if a in ("--log-file", "-l") and i + 1 < len(argv):
                i += 1
            i += 1
            continue
        out.append(a)
        i += 1
    return out


def _create_package(root: Path, output: Path | None = None) -> Path:
    root = Path(root).resolve()
    out = Path(output).resolve() if output else root / "packaged_output.tar.gz"
    subdirs = [root / d for d in ("GND", "EXC", "TP", "NEXAFS") if (root / d).is_dir()]
    if not subdirs:
        raise FileNotFoundError("No GND/, EXC/, TP/, or NEXAFS/ found")
    with tarfile.open(out, "w:gz") as tf:
        for d in subdirs:
            for f in d.rglob("*"):
                if f.is_file():
                    tf.add(f, arcname=f.relative_to(root))
    return out


def _spawn_quiet_background(
    log_file_override: Path | None,
    state: dict,
) -> None:
    cfg = load_config()
    log_dir = cfg.get("log_directory", "directory")
    default_log = Path.cwd() / "output.log"
    if log_dir != "directory" and isinstance(log_dir, str):
        log_path = Path(log_dir).expanduser().resolve()
        if log_path.is_dir():
            default_log = log_path / "output.log"
    log_path = log_file_override or state.get("log_file") or default_log
    log_path = Path(log_path).resolve()
    argv = [sys.executable, "-m", "dftlearn.cli"] + _argv_without_quiet_and_log()
    with log_path.open("w") as lf:
        proc = subprocess.Popen(
            argv,
            stdout=lf,
            stderr=subprocess.STDOUT,
            cwd=os.getcwd(),
            start_new_session=True,
        )
    typer.echo(f"Started in background (PID {proc.pid}), logging to {log_path}")
    raise typer.Exit(0)


_state: dict = {
    "verbose": False,
    "workers": None,
    "auto_workers": True,
    "max_workers": max(1, multiprocessing.cpu_count() - 1),
    "quiet": False,
    "log_file": None,
    "scan_dir": None,
}


def _atexit_cleanup() -> None:
    if _state.get("quiet"):
        return
    scan = _state.get("scan_dir")
    if scan is not None:
        _cleanup_fort_files(Path(scan))


def _signal_cleanup(signum: int | None, frame: object | None) -> None:
    _atexit_cleanup()
    sys.exit(128 + (signum if signum is not None else 0))


atexit.register(_atexit_cleanup)
signal.signal(signal.SIGINT, _signal_cleanup)
try:
    signal.signal(signal.SIGTERM, _signal_cleanup)
except (ValueError, OSError):
    pass

run_app = typer.Typer(help="Schedule and run StoBe DFT calculations.")


@run_app.callback()
def _run_callback(
    workers: int | None = typer.Option(
        None, "--workers", "-w", help="Number of workers"
    ),
    max_workers: int | None = typer.Option(
        None, "--max-workers", "-m", help="Max workers when auto"
    ),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose output"),
    quiet: bool = typer.Option(
        False, "--quiet", "-q", help="Run in background, log to file"
    ),
    log_file: Path | None = typer.Option(
        None, "--log-file", "-l", help="Log file when --quiet"
    ),
) -> None:
    cfg = load_config()
    _state["verbose"] = verbose
    _state["quiet"] = quiet
    _state["log_file"] = log_file.resolve() if log_file else None
    _state["auto_workers"] = workers is None
    _state["workers"] = workers
    _state["max_workers"] = (
        max_workers
        if max_workers is not None
        else int(cfg.get("max_workers", max(1, multiprocessing.cpu_count() - 1)))
    )
    if verbose and not quiet:
        if _state["auto_workers"]:
            typer.echo("Workers: one per atom, max %s" % _state["max_workers"])
        else:
            typer.echo("Using %s workers" % workers)


def _run_calc_type(
    calc_type: str,
    directory: Path,
    atom: list[str] | None,
) -> list[dict]:
    if _state.get("quiet"):
        _spawn_quiet_background(None, _state)
    _state["scan_dir"] = str(directory)
    typer.echo(f"Starting {calc_type} calculations in {directory}")
    if atom:
        typer.echo(f"Target atoms: {atom}")
    workers = _state["workers"]
    if _state["auto_workers"]:
        workers = _calculate_optimal_workers(
            str(directory), [calc_type], atom, _state["max_workers"]
        )
        if _state["verbose"]:
            typer.echo(
                f"Workers: {workers} (one per atom, max %s)" % _state["max_workers"]
            )
    runner = TyperSchedulerRunner()
    return runner.run_calculations_with_progress(
        str(directory), [calc_type], workers or 1, atom
    )


@run_app.command("gnd")
def _gnd(
    directory: Path = typer.Argument(
        ..., exists=True, help="Directory containing run files"
    ),
    atom: list[str] | None = typer.Option(
        None, "--atom", "-a", help="Specific atom(s)"
    ),
    quiet: bool = typer.Option(False, "--quiet", "-q"),
    log_file: Path | None = typer.Option(None, "--log-file", "-l"),
) -> None:
    if quiet:
        _state["quiet"] = True
    if log_file:
        _state["log_file"] = log_file.resolve()
    results = _run_calc_type("gnd", directory, atom)
    c = sum(1 for r in results if r["status"] == "completed")
    f = sum(1 for r in results if r["status"] == "failed")
    typer.echo("\nGround state calculations completed.")
    typer.echo(f"Completed: {c}, Failed: {f}")


@run_app.command("exc")
def _exc(
    directory: Path = typer.Argument(..., exists=True),
    atom: list[str] | None = typer.Option(None, "--atom", "-a"),
    quiet: bool = typer.Option(False, "--quiet", "-q"),
    log_file: Path | None = typer.Option(None, "--log-file", "-l"),
) -> None:
    if quiet:
        _state["quiet"] = True
    if log_file:
        _state["log_file"] = log_file.resolve()
    results = _run_calc_type("exc", directory, atom)
    c = sum(1 for r in results if r["status"] == "completed")
    f = sum(1 for r in results if r["status"] == "failed")
    typer.echo("\nExcited state calculations completed.")
    typer.echo(f"Completed: {c}, Failed: {f}")


@run_app.command("tp")
def _tp(
    directory: Path = typer.Argument(..., exists=True),
    atom: list[str] | None = typer.Option(None, "--atom", "-a"),
    quiet: bool = typer.Option(False, "--quiet", "-q"),
    log_file: Path | None = typer.Option(None, "--log-file", "-l"),
) -> None:
    if quiet:
        _state["quiet"] = True
    if log_file:
        _state["log_file"] = log_file.resolve()
    results = _run_calc_type("tp", directory, atom)
    c = sum(1 for r in results if r["status"] == "completed")
    f = sum(1 for r in results if r["status"] == "failed")
    typer.echo("\nTransition potential calculations completed.")
    typer.echo(f"Completed: {c}, Failed: {f}")


@run_app.command("xas")
def _xas(
    directory: Path = typer.Argument(..., exists=True),
    atom: list[str] | None = typer.Option(None, "--atom", "-a"),
    quiet: bool = typer.Option(False, "--quiet", "-q"),
    log_file: Path | None = typer.Option(None, "--log-file", "-l"),
) -> None:
    if quiet:
        _state["quiet"] = True
    if log_file:
        _state["log_file"] = log_file.resolve()
    results = _run_calc_type("xas", directory, atom)
    c = sum(1 for r in results if r["status"] == "completed")
    f = sum(1 for r in results if r["status"] == "failed")
    typer.echo("\nXAS calculations completed.")
    typer.echo(f"Completed: {c}, Failed: {f}")


@run_app.command("seq")
def _seq(
    directory: Path = typer.Argument(..., exists=True),
    atom: list[str] | None = typer.Option(None, "--atom", "-a"),
    quiet: bool = typer.Option(False, "--quiet", "-q"),
    log_file: Path | None = typer.Option(None, "--log-file", "-l"),
) -> None:
    if quiet:
        _state["quiet"] = True
    if log_file:
        _state["log_file"] = log_file.resolve()
    results = _run_calc_type("seq", directory, atom)
    c = sum(1 for r in results if r["status"] == "completed")
    f = sum(1 for r in results if r["status"] == "failed")
    typer.echo("\nSequential calculations completed.")
    typer.echo(f"Completed: {c}, Failed: {f}")


@run_app.command("all")
def _all(
    directory: Path = typer.Argument(..., exists=True),
    atom: list[str] | None = typer.Option(None, "--atom", "-a"),
    quiet: bool = typer.Option(False, "--quiet", "-q"),
    log_file: Path | None = typer.Option(None, "--log-file", "-l"),
) -> None:
    if quiet:
        _state["quiet"] = True
    if log_file:
        _state["log_file"] = log_file.resolve()
    if _state.get("quiet"):
        _spawn_quiet_background(log_file, _state)
    _state["scan_dir"] = str(directory)
    calc_types = ["gnd", "exc", "tp", "xas"]
    workers = _state["workers"]
    if _state["auto_workers"]:
        workers = _calculate_optimal_workers(
            str(directory), calc_types, atom, _state["max_workers"]
        )
        if _state["verbose"]:
            typer.echo(
                f"Workers: {workers} (one per atom, max %s)" % _state["max_workers"]
            )
    typer.echo(f"Starting all calculations in {directory}")
    typer.echo(f"Sequence: {' -> '.join(calc_types)}")
    scan_path = Path(directory).resolve()
    molecule_name = scan_path.name
    cfg = load_config()
    logs_dir, summary_dir = _resolve_logs_dirs(scan_path, molecule_name, cfg)
    run_ts = datetime.now().strftime("%Y%m%d%H%M%S")
    logs_dir.mkdir(parents=True, exist_ok=True)
    all_results: list[dict] = []
    for ct in calc_types:
        typer.echo(f"\nStarting {ct.upper()} calculations")
        runner = TyperSchedulerRunner()
        res = runner.run_calculations_with_progress(
            str(directory),
            [ct],
            workers or 1,
            atom,
            override_logs_dir=logs_dir,
            override_run_timestamp=run_ts,
            override_summary_dir=summary_dir,
            override_molecule_name=molecule_name,
            override_molecule_dir=scan_path,
            write_summary=False,
        )
        all_results.extend(res)
        if ct != calc_types[-1]:
            typer.echo("Pausing 2 seconds before next calculation type...")
            time.sleep(2)
    summary_path = summary_dir / f"{molecule_name}_{run_ts}.txt"
    summary_dir.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w") as sf:
        sf.write("atom\tcalc_type\tstart_time\tend_time\tduration_s\tstatus\n")
        for r in all_results:
            st = r["start_time"].isoformat() if r.get("start_time") else ""
            et = r["end_time"].isoformat() if r.get("end_time") else ""
            dur = f"{r['duration']:.2f}" if r.get("duration") is not None else ""
            sf.write(
                f"{r['directory']}\t{r['calc_type']}\t{st}\t{et}\t{dur}\t{r['status']}\n"
            )
    c = sum(1 for r in all_results if r["status"] == "completed")
    f = sum(1 for r in all_results if r["status"] == "failed")
    typer.echo("\nAll calculations completed.")
    typer.echo(f"Total processed: {len(all_results)} calculations")
    typer.echo(f"Completed: {c}, Failed: {f}")
    typer.echo(f"Summary: {summary_path}")
    try:
        out = _create_package(Path(directory))
        typer.echo(f"Created {out}")
    except FileNotFoundError:
        pass


@run_app.command("explore")
def _explore(directory: Path = typer.Argument(..., exists=True)) -> None:
    typer.echo(f"Exploring: {directory}")
    typer.echo("-" * 40)
    run_files = list(directory.glob("**/*.run"))
    if not run_files:
        typer.echo("No run files found")
        return
    atoms: set[str] = set()
    calc_types: set[str] = set()
    for rf in run_files:
        stem = rf.stem
        atom_match = stem[:-3] if len(stem) > 3 else stem
        ct = stem[-3:] if len(stem) > 3 else ""
        atoms.add(atom_match)
        calc_types.add(ct)
    typer.echo(f"Available atoms ({len(atoms)}): {sorted(atoms)}")
    typer.echo(f"Available calculation types ({len(calc_types)}): {sorted(calc_types)}")
    typer.echo(f"Total run files: {len(run_files)}")
    typer.echo("\nExample files:")
    for rf in sorted(run_files)[:5]:
        typer.echo(f"  {rf.relative_to(directory)}")
    if len(run_files) > 5:
        typer.echo(f"  ... and {len(run_files) - 5} more")


@run_app.command("organize")
def _organize(directory: Path = typer.Argument(..., exists=True)) -> None:
    root = Path(directory).resolve()
    for sub in ("GND", "EXC", "TP", "NEXAFS"):
        (root / sub).mkdir(parents=True, exist_ok=True)
    to_move: list[tuple[Path, Path]] = []
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        try:
            rel = path.relative_to(root)
        except ValueError:
            continue
        if rel.parts and rel.parts[0] in ("GND", "EXC", "TP", "NEXAFS"):
            continue
        name = path.name
        if name.endswith("gnd.out") or name.endswith("gnd.molden"):
            dest = root / "GND" / name
        elif name.endswith("exc.out") or name.endswith("exc.molden"):
            dest = root / "EXC" / name
        elif name.endswith("tp.out") or name.endswith("tp.molden"):
            dest = root / "TP" / name
        elif name.endswith("xas.out") or name.endswith(".out"):
            dest = root / "NEXAFS" / name
        else:
            continue
        if path.resolve() != dest.resolve():
            to_move.append((path, dest))
    for src, dst in to_move:
        shutil.move(str(src), str(dst))
    typer.echo(f"Organized {len(to_move)} files into GND/, EXC/, TP/, NEXAFS/")


@run_app.command("package")
def _package(
    directory: Path = typer.Argument(..., exists=True),
    output: Path | None = typer.Option(None, "--output", "-o"),
) -> None:
    try:
        out = _create_package(Path(directory), output)
        typer.echo(f"Created {out}")
    except FileNotFoundError as e:
        typer.echo(str(e))
        raise typer.Exit(1) from e
