"""dftrun.toml config loader."""

from __future__ import annotations

import multiprocessing
import tomllib
from pathlib import Path
from typing import Any

CONFIG_NAMES = ("dftrun.toml", ".dftrun.toml")


def _default_max_workers() -> int:
    return max(1, multiprocessing.cpu_count() - 1)


def find_config(start: Path | None = None) -> Path | None:
    """Locate dftrun.toml from start dir up, then ~/.config/dftrun."""
    start = (start or Path.cwd()).resolve()
    for parent in [start, *start.parents]:
        for name in CONFIG_NAMES:
            p = parent / name
            if p.is_file():
                return p
    global_dir = Path.home() / ".config" / "dftrun" / "dftrun.toml"
    if global_dir.is_file():
        return global_dir
    return None


def load_config(config_path: Path | None = None) -> dict[str, Any]:
    """Load dftrun.toml into a dict (max_workers, log_directory)."""
    path = config_path or find_config()
    out: dict[str, Any] = {
        "max_workers": _default_max_workers(),
        "log_directory": "directory",
    }
    if not path:
        return out
    try:
        with path.open("rb") as f:
            data = tomllib.load(f)
    except OSError:
        return out
    if isinstance(data.get("max_workers"), int) and data["max_workers"] >= 1:
        out["max_workers"] = data["max_workers"]
    if isinstance(data.get("log_directory"), str):
        out["log_directory"] = data["log_directory"]
    return out
