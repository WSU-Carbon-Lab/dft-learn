# dftrun CLI

Build and run StoBe DFT (transition-potential) workflows. Replaces the legacy `generate.py` / `schedule.py` scripts.

## Installation

```bash
uv tool install --from 'git+https://github.com/<org>/dft-learn' dftrun
```

Or from a clone:

```bash
cd /path/to/dft-learn
uv tool install .
```

Verify:

```bash
uvx dftrun --help
```

## Usage

### Build (generate inputs)

Generate StoBe run files from a directory containing `molConfig.py` and an XYZ geometry:

```bash
uvx dftrun build /path/to/molecule
uvx dftrun build /path/to/molecule /path/to/molecule.xyz
uvx dftrun build /path/to/molecule --xyz /path/to/other.xyz
uvx dftrun build /path/to/molecule -v
```

### Run (execute calculations)

```bash
uvx dftrun run gnd /path/to/molecule
uvx dftrun run exc /path/to/molecule
uvx dftrun run tp /path/to/molecule
uvx dftrun run xas /path/to/molecule
uvx dftrun run all /path/to/molecule
uvx dftrun run all /path/to/molecule --quiet
uvx dftrun run gnd /path/to/molecule --atom C1 --atom C2
uvx dftrun run -v all /path/to/molecule
uvx dftrun run --workers 4 gnd /path/to/molecule
uvx dftrun run explore /path/to/molecule
uvx dftrun run organize /path/to/molecule
uvx dftrun run package /path/to/molecule
```

## Config: `dftrun.toml`

Config is read from the current directory, then parents, then `~/.config/dftrun/dftrun.toml`. Use `dftrun.toml` or `.dftrun.toml`.

Example:

```toml
max_workers = 31
log_directory = "directory"
```

- **max_workers**: Cap when using auto workers (one per atom). Override with `--max-workers` / `-m`.
- **log_directory**: `"directory"` (default) – logs under the molecule directory; or an absolute path – logs under that global directory.

CLI flags override config.

## Logging

- **Per-job logs**: `logs/<atom>_<calctype>_<timestamp>.txt` under the molecule directory (or under `log_directory/<molecule_name>/` when using a global `log_directory`).
- **Summary**: `<molecule_name>_<timestamp>.txt` in the molecule root (or same global subdir), with `atom`, `calc_type`, `start_time`, `end_time`, `duration_s`, `status`.

With `--quiet`, the process runs in the background and writes to a single log file (default `output.log`, or `--log-file`). Per-job and summary logs still apply to the background run.

## Requirements

- Python 3.12+
- [uv](https://docs.astral.sh/uv/)
- StoBe: `STOBE_HOME` (or default `/bin/stobe`), and `StoBe.x`, `xrayspec.x` on `PATH`.

## molConfig

`dftrun build` loads `molConfig.py` from the run directory. See fixture `tests/fixtures/stobe_minimal/molConfig.py` for a minimal example. Legacy generator docs and `molConfig` usage remain in the deprecated scripts under `src/dftlearn/cli/_deprecated/` (reference only).
