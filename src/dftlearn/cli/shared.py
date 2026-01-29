"""Shared helpers for dftrun CLI (molConfig loading)."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType


def load_mol_config(config_dir: Path) -> ModuleType:
    """Load molConfig.py from config_dir and return the module."""
    config_path = config_dir / "molConfig.py"
    if not config_path.exists():
        msg = f"molConfig.py not found in {config_dir}"
        raise FileNotFoundError(msg)
    spec = importlib.util.spec_from_file_location("molConfig", config_path)
    if spec is None or spec.loader is None:
        msg = f"Could not load molConfig from {config_path}"
        raise ImportError(msg)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod
