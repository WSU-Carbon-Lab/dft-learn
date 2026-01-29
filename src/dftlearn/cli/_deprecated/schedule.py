#!/usr/bin/env -S uv run
# /// script
# requires-python = ">=3.8"
# dependencies = [
#     "typer",
#     "rich",
# ]
# ///
"""
StoBe DFT Calculation Scheduler CLI

A command-line interface for scheduling and running StoBe DFT calculations
with parallel execution and real-time progress monitoring.

Usage:
    uv run schedule.py gnd .                     # Run ground states in current directory
    uv run schedule.py gnd /path/to/directory    # Run ground states in specific directory
    uv run schedule.py all .                     # Run all calculation types sequentially
    uv run schedule.py gnd . --atom C1           # Run ground state for specific atom
    uv run schedule.py --workers 8 gnd .         # Use 8 parallel workers
"""

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

# Create console instance
console = Console()


def _count_atoms(scan_dir: str, calc_types: list[str], atoms: list[str] | None = None) -> int:
    scan_path = Path(scan_dir)
    directories: set[str] = set()
    for calc_type in calc_types:
        if atoms:
            for atom in atoms:
                for run_file in scan_path.glob(f"**/{atom}{calc_type}.run"):
                    directories.add(run_file.parent.name)
        else:
            for run_file in scan_path.glob(f"**/*{calc_type}.run"):
                directories.add(run_file.parent.name)
    return len(directories)


def calculate_optimal_workers(
    scan_dir: str,
    calc_types: list[str],
    atoms: list[str] | None = None,
    max_workers: int | None = None,
) -> int:
    num_atoms = _count_atoms(scan_dir, calc_types, atoms)
    cap = max_workers if max_workers is not None else max(1, multiprocessing.cpu_count() - 1)
    return min(num_atoms, cap) if num_atoms else 0


@dataclass
class DirectoryJob:
    """Represents a job for a specific directory (e.g., C1, C2)"""

    directory: str
    run_file: Path
    calc_type: str
    status: str = "pending"  # pending, running, completed, failed
    start_time: datetime | None = None
    end_time: datetime | None = None
    duration: float | None = None
    error_msg: str | None = None
    thread: Thread | None = None


class StoBeJobManager:
    """Manages execution and monitoring of StoBe DFT calculations"""

    def __init__(self, working_dir: Path | None = None):
        self.working_dir = working_dir or Path.cwd()
        self.jobs = {}  # Track running jobs
        self.completed_jobs = {}  # Track completed jobs
        self.job_queue = queue.Queue()
        self.max_concurrent_jobs = 4  # Adjust based on system resources

    def find_run_files(self, pattern: str = "**/*.run") -> list[Path]:
        """Find all run files matching the pattern"""
        run_files = list(self.working_dir.glob(pattern))
        return sorted(run_files)

    def execute_run_file(self, run_file: Path, background: bool = True) -> dict:
        """Execute a single StoBe run file"""
        if not run_file.exists():
            raise FileNotFoundError(f"Run file not found: {run_file}")

        # Change to the run file directory
        run_dir = run_file.parent
        run_script = run_file.name

        # Prepare job info
        job_info = {
            "run_file": run_file,
            "run_dir": run_dir,
            "start_time": datetime.now(),
            "status": "starting",
            "process": None,
            "output_file": run_dir / f"{run_script}.out",
            "error_file": run_dir / f"{run_script}.err",
        }

        try:
            # Execute the run script
            if background:
                process = subprocess.Popen(
                    ["bash", run_script],
                    cwd=run_dir,
                    stdout=open(job_info["output_file"], "w"),
                    stderr=open(job_info["error_file"], "w"),
                    shell=False,
                )
                job_info["process"] = process
                job_info["status"] = "running"
                job_info["pid"] = process.pid
            else:
                # Synchronous execution
                result = subprocess.run(
                    ["bash", run_script], cwd=run_dir, capture_output=True, text=True
                )
                job_info["returncode"] = result.returncode
                job_info["stdout"] = result.stdout
                job_info["stderr"] = result.stderr
                job_info["status"] = "completed" if result.returncode == 0 else "failed"
                job_info["end_time"] = datetime.now()

        except Exception as e:
            job_info["status"] = "error"
            job_info["error"] = str(e)
            job_info["end_time"] = datetime.now()

        return job_info


class TyperSchedulerRunner:
    """Enhanced parallel runner with Typer progress bars"""

    def __init__(self):
        self.jobs: dict[str, DirectoryJob] = {}
        self.lock = Lock()
        self.completed_count = 0
        self.total_count = 0
        self.start_time = None
        self.console = Console()

    def find_run_files_in_directory(
        self,
        scan_dir: str,
        calc_types: list[str] | None = None,
        atoms: list[str] | None = None,
    ):
        """Find run files in a specific scan directory with filtering"""
        scan_path = Path(scan_dir)

        if not scan_path.exists():
            typer.echo(f"❌ Scan directory not found: {scan_dir}")
            return []

        # Build pattern based on filters
        patterns = []

        if calc_types is None:
            calc_types = ["gnd", "exc", "tp", "xas", "seq"]

        if atoms is None:
            # Find all available atoms
            for calc_type in calc_types:
                patterns.append(f"**/*{calc_type}.run")
        else:
            # Specific atoms
            for atom in atoms:
                for calc_type in calc_types:
                    patterns.append(f"**/{atom}{calc_type}.run")

        # Collect all files
        run_files = []
        for pattern in patterns:
            files = list(scan_path.glob(pattern))
            run_files.extend(files)

        return sorted(list(set(run_files)))  # Remove duplicates and sort

    def find_directories_with_run_files(
        self, scan_dir: str, calc_types: list[str], atoms: list[str] | None = None
    ) -> dict[str, list[Path]]:
        """Find all directories containing run files, grouped by directory name"""
        scan_path = Path(scan_dir)
        directories = {}

        for calc_type in calc_types:
            if atoms:
                # Specific atoms
                for atom in atoms:
                    pattern = f"**/{atom}{calc_type}.run"
                    run_files = list(scan_path.glob(pattern))
                    for run_file in run_files:
                        dir_name = run_file.parent.name
                        if dir_name not in directories:
                            directories[dir_name] = []
                        directories[dir_name].append(run_file)
            else:
                # All atoms
                pattern = f"**/*{calc_type}.run"
                run_files = list(scan_path.glob(pattern))
                for run_file in run_files:
                    dir_name = run_file.parent.name
                    if dir_name not in directories:
                        directories[dir_name] = []
                    directories[dir_name].append(run_file)

        return directories

    def execute_directory_job(self, job: DirectoryJob) -> DirectoryJob:
        """Execute a single directory job"""
        job.status = "running"
        job.start_time = datetime.now()

        calculation_done = threading.Event()
        calculation_result: dict = {"job_info": None, "exception": None}

        def run_calculation():
            try:
                temp_manager = StoBeJobManager()
                job_info = temp_manager.execute_run_file(job.run_file, background=False)
                calculation_result["job_info"] = job_info
            except Exception as e:
                calculation_result["exception"] = e
            finally:
                calculation_done.set()

        calc_thread = threading.Thread(target=run_calculation, daemon=True)
        calc_thread.start()

        # Wait for calculation to complete
        while not calculation_done.is_set():
            time.sleep(1)

        calc_thread.join()

        job.end_time = datetime.now()
        job.duration = (job.end_time - job.start_time).total_seconds()

        if calculation_result["exception"]:
            job.status = "failed"
            job.error_msg = str(calculation_result["exception"])
        elif calculation_result["job_info"]:
            job_info = calculation_result["job_info"]
            if job_info["status"] == "completed":
                job.status = "completed"
            else:
                job.status = "failed"
                job.error_msg = job_info.get("error", "Unknown error")
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
    ):
        """Run calculations with Typer progress bars"""
        # Find directories with run files
        directories = self.find_directories_with_run_files(scan_dir, calc_types, atoms)

        if not directories:
            typer.echo(f"❌ No directories with run files found in {scan_dir}")
            return []

        for dir_name, run_files in directories.items():
            for run_file in run_files:
                calc_type = run_file.stem[-3:] if len(run_file.stem) > 3 else "unknown"
                job = DirectoryJob(
                    directory=dir_name, run_file=run_file, calc_type=calc_type
                )
                job_key = f"{dir_name}:{run_file.name}"
                self.jobs[job_key] = job

        self.total_count = len(self.jobs)
        self.completed_count = 0
        self.start_time = datetime.now()

        # Setup progress tracking with only overall progress bar
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            TextColumn("•"),
            TimeElapsedColumn(),
            console=self.console,
        ) as progress:

            # Create main progress task only
            main_task = progress.add_task(
                f"[cyan]Processing {self.total_count} calculations ({', '.join(calc_types)})",
                total=self.total_count
            )

            def directory_worker(job: DirectoryJob):
                """Worker function for each directory thread"""
                result = self.execute_directory_job(job)
                job_key = f"{job.directory}:{job.run_file.name}"

                with self.lock:
                    self.jobs[job_key] = result
                    if result.status in ["completed", "failed"]:
                        self.completed_count += 1
                        status_text = "completed" if result.status == "completed" else "failed"
                        progress.update(
                            main_task,
                            completed=self.completed_count,
                            description=f"[cyan]Processing {self.total_count} calculations | {status_text} {job.directory} {job.calc_type} ({result.duration:.1f}s)"
                        )

            # Execute jobs with threading
            active_threads = []
            job_items = list(self.jobs.items())
            job_index = 0

            while job_index < len(job_items) or active_threads:
                # Start new threads up to max_workers
                while len(active_threads) < max_workers and job_index < len(job_items):
                    dir_name, job = job_items[job_index]
                    thread = Thread(target=directory_worker, args=(job,))
                    thread.daemon = True
                    job.thread = thread
                    thread.start()
                    active_threads.append(thread)
                    job_index += 1
                    time.sleep(0.1)

                # Clean up completed threads
                active_threads = [t for t in active_threads if t.is_alive()]
                time.sleep(0.5)

        # Results summary
        results = []
        for job in self.jobs.values():
            results.append({
                "directory": job.directory,
                "calc_type": job.calc_type,
                "status": job.status,
                "duration": job.duration,
                "error": job.error_msg,
            })

        return results


# Create the main Typer app
app = typer.Typer(help="StoBe DFT Calculation Scheduler - Schedule and run StoBe DFT calculations with parallel execution.")

# Global state for CLI options
state = {
    "verbose": False,
    "workers": None,
    "auto_workers": True,
    "max_workers": max(1, multiprocessing.cpu_count() - 1),
    "quiet": False,
    "log_file": None,
    "scan_dir": None,
}


def _cleanup_fort_files(scan_dir: Path) -> None:
    root = Path(scan_dir).resolve()
    for f in root.rglob("fort.*"):
        if f.is_file():
            try:
                f.unlink()
            except OSError:
                pass


def _atexit_cleanup() -> None:
    if state.get("quiet"):
        return
    scan = state.get("scan_dir")
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


def _spawn_quiet_background(log_file_override: Path | None = None) -> None:
    log_path = log_file_override or state["log_file"] or Path.cwd() / "output.log"
    log_path = Path(log_path).resolve()
    uv = shutil.which("uv") or "uv"
    script = Path(__file__).resolve()
    child_argv = [uv, "run", str(script)] + _argv_without_quiet_and_log()
    with log_path.open("w") as lf:
        proc = subprocess.Popen(
            child_argv,
            stdout=lf,
            stderr=subprocess.STDOUT,
            cwd=os.getcwd(),
            start_new_session=True,
        )
    typer.echo(f"Started in background (PID {proc.pid}), logging to {log_path}")
    raise typer.Exit(0)


@app.callback()
def main(
    workers: int | None = typer.Option(None, "--workers", "-w", help="Number of workers (default: one per atom, capped by --max-workers)"),
    max_workers: int | None = typer.Option(None, "--max-workers", "-m", help="Max workers when auto (default: CPU count - 1)"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose output"),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Run in background and log to output.log"),
    log_file: Path | None = typer.Option(None, "--log-file", "-l", help="Log file when --quiet (default: output.log)"),
):
    """StoBe DFT Calculation Scheduler

    Schedule and run StoBe DFT calculations with parallel execution and progress bars.
    """
    state["quiet"] = quiet
    state["log_file"] = log_file.resolve() if log_file else None
    state["auto_workers"] = workers is None
    state["workers"] = workers
    state["max_workers"] = max_workers if max_workers is not None else max(1, multiprocessing.cpu_count() - 1)
    state["verbose"] = verbose

    if verbose and not quiet:
        if state["auto_workers"]:
            typer.echo("Workers: one per atom, max %s" % state["max_workers"])
        else:
            typer.echo("Using %s workers" % workers)


@app.command()
def gnd(
    directory: Path = typer.Argument(..., exists=True, help="Directory containing run files"),
    atom: list[str] | None = typer.Option(None, "--atom", "-a", help="Specific atom(s) to calculate (e.g., C1, C2)"),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Run in background and log to output.log"),
    log_file: Path | None = typer.Option(None, "--log-file", "-l", help="Log file when --quiet (default: output.log)"),
):
    """Run ground state calculations"""
    if quiet or state.get("quiet"):
        _spawn_quiet_background(log_file)
    state["scan_dir"] = str(directory)
    typer.echo(f"Starting ground state calculations in {directory}")
    if atom:
        typer.echo(f"Target atoms: {atom}")

    workers = state["workers"]
    if state["auto_workers"]:
        workers = calculate_optimal_workers(str(directory), ["gnd"], atom, state["max_workers"])
        if state["verbose"]:
            typer.echo(f"Workers: {workers} (one per atom, max %s)" % state["max_workers"])

    runner = TyperSchedulerRunner()
    results = runner.run_calculations_with_progress(str(directory), ["gnd"], workers, atom)

    completed = sum(1 for r in results if r["status"] == "completed")
    failed = sum(1 for r in results if r["status"] == "failed")

    typer.echo("\n🎉 Ground state calculations completed!")
    typer.echo(f"✅ Completed: {completed}, ❌ Failed: {failed}")


@app.command()
def exc(
    directory: Path = typer.Argument(..., exists=True, help="Directory containing run files"),
    atom: list[str] | None = typer.Option(None, "--atom", "-a", help="Specific atom(s) to calculate"),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Run in background and log to output.log"),
    log_file: Path | None = typer.Option(None, "--log-file", "-l", help="Log file when --quiet (default: output.log)"),
):
    """Run excited state calculations"""
    if quiet or state.get("quiet"):
        _spawn_quiet_background(log_file)
    state["scan_dir"] = str(directory)
    typer.echo(f"Starting excited state calculations in {directory}")

    workers = state["workers"]
    if state["auto_workers"]:
        workers = calculate_optimal_workers(str(directory), ["exc"], atom, state["max_workers"])
        if state["verbose"]:
            typer.echo(f"Workers: {workers} (one per atom, max %s)" % state["max_workers"])

    runner = TyperSchedulerRunner()
    results = runner.run_calculations_with_progress(str(directory), ["exc"], workers, atom)

    completed = sum(1 for r in results if r["status"] == "completed")
    failed = sum(1 for r in results if r["status"] == "failed")

    typer.echo("\n🎉 Excited state calculations completed!")
    typer.echo(f"✅ Completed: {completed}, ❌ Failed: {failed}")


@app.command()
def tp(
    directory: Path = typer.Argument(..., exists=True, help="Directory containing run files"),
    atom: list[str] | None = typer.Option(None, "--atom", "-a", help="Specific atom(s) to calculate"),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Run in background and log to output.log"),
    log_file: Path | None = typer.Option(None, "--log-file", "-l", help="Log file when --quiet (default: output.log)"),
):
    """Run transition potential calculations"""
    if quiet or state.get("quiet"):
        _spawn_quiet_background(log_file)
    state["scan_dir"] = str(directory)
    typer.echo(f"Starting transition potential calculations in {directory}")

    workers = state["workers"]
    if state["auto_workers"]:
        workers = calculate_optimal_workers(str(directory), ["tp"], atom, state["max_workers"])
        if state["verbose"]:
            typer.echo(f"Workers: {workers} (one per atom, max %s)" % state["max_workers"])

    runner = TyperSchedulerRunner()
    results = runner.run_calculations_with_progress(str(directory), ["tp"], workers, atom)

    completed = sum(1 for r in results if r["status"] == "completed")
    failed = sum(1 for r in results if r["status"] == "failed")

    typer.echo("\n🎉 Transition potential calculations completed!")
    typer.echo(f"✅ Completed: {completed}, ❌ Failed: {failed}")


@app.command()
def xas(
    directory: Path = typer.Argument(..., exists=True, help="Directory containing run files"),
    atom: list[str] | None = typer.Option(None, "--atom", "-a", help="Specific atom(s) to calculate"),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Run in background and log to output.log"),
    log_file: Path | None = typer.Option(None, "--log-file", "-l", help="Log file when --quiet (default: output.log)"),
):
    """Run XAS calculations"""
    if quiet or state.get("quiet"):
        _spawn_quiet_background(log_file)
    state["scan_dir"] = str(directory)
    typer.echo(f"Starting XAS calculations in {directory}")

    workers = state["workers"]
    if state["auto_workers"]:
        workers = calculate_optimal_workers(str(directory), ["xas"], atom, state["max_workers"])
        if state["verbose"]:
            typer.echo(f"Workers: {workers} (one per atom, max %s)" % state["max_workers"])

    runner = TyperSchedulerRunner()
    results = runner.run_calculations_with_progress(str(directory), ["xas"], workers, atom)

    completed = sum(1 for r in results if r["status"] == "completed")
    failed = sum(1 for r in results if r["status"] == "failed")

    typer.echo("\n🎉 XAS calculations completed!")
    typer.echo(f"✅ Completed: {completed}, ❌ Failed: {failed}")


@app.command()
def seq(
    directory: Path = typer.Argument(..., exists=True, help="Directory containing run files"),
    atom: list[str] | None = typer.Option(None, "--atom", "-a", help="Specific atom(s) to calculate"),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Run in background and log to output.log"),
    log_file: Path | None = typer.Option(None, "--log-file", "-l", help="Log file when --quiet (default: output.log)"),
):
    """Run sequential calculations"""
    if quiet or state.get("quiet"):
        _spawn_quiet_background(log_file)
    state["scan_dir"] = str(directory)
    typer.echo(f"Starting sequential calculations in {directory}")

    workers = state["workers"]
    if state["auto_workers"]:
        workers = calculate_optimal_workers(str(directory), ["seq"], atom, state["max_workers"])
        if state["verbose"]:
            typer.echo(f"Workers: {workers} (one per atom, max %s)" % state["max_workers"])

    runner = TyperSchedulerRunner()
    results = runner.run_calculations_with_progress(str(directory), ["seq"], workers, atom)

    completed = sum(1 for r in results if r["status"] == "completed")
    failed = sum(1 for r in results if r["status"] == "failed")

    typer.echo("\n🎉 Sequential calculations completed!")
    typer.echo(f"✅ Completed: {completed}, ❌ Failed: {failed}")


@app.command()
def all(
    directory: Path = typer.Argument(..., exists=True, help="Directory containing run files"),
    atom: list[str] | None = typer.Option(None, "--atom", "-a", help="Specific atom(s) to calculate"),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Run in background and log to output.log"),
    log_file: Path | None = typer.Option(None, "--log-file", "-l", help="Log file when --quiet (default: output.log)"),
):
    """Run all calculation types sequentially (gnd -> exc -> tp -> xas)"""
    if quiet or state.get("quiet"):
        _spawn_quiet_background(log_file)
    state["scan_dir"] = str(directory)
    calc_types = ["gnd", "exc", "tp", "xas"]
    all_results = []

    workers = state["workers"]
    if state["auto_workers"]:
        workers = calculate_optimal_workers(str(directory), calc_types, atom, state["max_workers"])
        if state["verbose"]:
            typer.echo(f"Workers: {workers} (one per atom, max %s)" % state["max_workers"])

    typer.echo(f"Starting all calculations in {directory}")
    typer.echo(f"Sequence: {' -> '.join(calc_types)}")

    for calc_type in calc_types:
        typer.echo(f"\nStarting {calc_type.upper()} calculations")

        runner = TyperSchedulerRunner()
        results = runner.run_calculations_with_progress(str(directory), [calc_type], workers, atom)
        all_results.extend(results)

        if calc_type != calc_types[-1]:
            typer.echo("Pausing 2 seconds before next calculation type...")
            time.sleep(2)

    completed = sum(1 for r in all_results if r["status"] == "completed")
    failed = sum(1 for r in all_results if r["status"] == "failed")

    typer.echo("\nAll calculations completed!")
    typer.echo(f"Total processed: {len(all_results)} calculations")
    typer.echo(f"Completed: {completed}, Failed: {failed}")

    try:
        out = _create_package(Path(directory))
        typer.echo(f"Created {out}")
    except FileNotFoundError:
        pass


@app.command()
def explore(
    directory: Path = typer.Argument(..., exists=True, help="Directory to explore"),
):
    """Explore available calculations in a directory"""
    typer.echo(f"🔍 Exploring: {directory}")
    typer.echo("-" * 40)

    # Find all run files
    run_files = list(directory.glob("**/*.run"))

    if not run_files:
        typer.echo("❌ No run files found")
        return

    # Analyze by atom and calculation type
    atoms = set()
    calc_types = set()

    for run_file in run_files:
        filename = run_file.stem
        atom_match = filename[:-3] if len(filename) > 3 else filename
        calc_type = filename[-3:] if len(filename) > 3 else ""
        atoms.add(atom_match)
        calc_types.add(calc_type)

    typer.echo(f"⚛️  Available atoms ({len(atoms)}): {sorted(atoms)}")
    typer.echo(f"🔬 Available calculation types ({len(calc_types)}): {sorted(calc_types)}")
    typer.echo(f"📁 Total run files: {len(run_files)}")

    # Show examples
    typer.echo("\n📄 Example files:")
    for i, run_file in enumerate(sorted(run_files)[:5]):
        rel_path = run_file.relative_to(directory)
        typer.echo(f"  {rel_path}")
    if len(run_files) > 5:
        typer.echo(f"  ... and {len(run_files) - 5} more")


@app.command()
def organize(
    directory: Path = typer.Argument(..., exists=True, help="Directory containing C1/, C2/, ... with outputs"),
):
    """Move *gnd.out, *exc.out, *tp.out, *xas.out, *.out, *.molden from atom dirs into GND/, EXC/, TP/, NEXAFS/"""
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
        if any(rel.parts and rel.parts[0] == d for d in ("GND", "EXC", "TP", "NEXAFS")):
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


@app.command()
def package(
    directory: Path = typer.Argument(..., exists=True, help="Directory containing GND/, EXC/, TP/, NEXAFS/"),
    output: Path | None = typer.Option(None, "--output", "-o", help="Output tarball path (default: packaged_output.tar.gz in directory)"),
):
    """Create packaged_output.tar.gz from GND/, EXC/, TP/, NEXAFS/"""
    try:
        out = _create_package(Path(directory), output)
        typer.echo(f"Created {out}")
    except FileNotFoundError as e:
        typer.echo(str(e))
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
