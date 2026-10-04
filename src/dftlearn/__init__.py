"""dftlearn: analysis helpers for DFT and StoBe-style spectroscopy workflows."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("dft-learn")
except PackageNotFoundError:
    __version__ = "0.0.0"

__all__: list[str] = ["__version__"]
