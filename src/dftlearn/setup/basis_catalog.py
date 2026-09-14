"""StoBe basis-library catalog and validation against ``baslib*.new7``.

The bundled JSON under ``data/stobe_basis_catalog.json`` lists every
``A-`` / ``O-`` / ``P-`` basis name parsed from the reference StoBe export
(``baslib.new7``, ``baslibD.new7``) plus symmetry groups from ``symbasis.new``.
Callers validate preset strings and ``molConfig.py`` assignments before
scheduling calculations.
"""

from __future__ import annotations

import json
import re
from functools import lru_cache
from importlib import resources
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable

    from dftlearn.setup.types import StoBeBasisPreset

_CATALOG_FILENAME = "stobe_basis_catalog.json"
_BASIS_NAME = re.compile(r"^(A-|O-|P-)[A-Z0-9(+).,\- /;*+|a-z]+")
_MOLCONFIG_BASIS = re.compile(
    r'^\s*elem\d+_(?:A|O|MCP)basis(?:_[ab])?\s*=\s*"([^"]*)"',
    re.MULTILINE,
)


@lru_cache(maxsize=1)
def bundled_basis_catalog() -> dict[str, object]:
    """Load the shipped StoBe basis-name catalog from package data.

    Returns
    -------
    dict
        Keys ``basis_names`` (sorted list of strings) and ``symmetry_groups``.
    """
    data_path = resources.files("dftlearn.setup.data").joinpath(_CATALOG_FILENAME)
    text = data_path.read_text(encoding="utf-8")
    return json.loads(text)


@lru_cache(maxsize=1)
def bundled_basis_names() -> frozenset[str]:
    """Return the set of known StoBe basis names from the bundled catalog."""
    catalog = bundled_basis_catalog()
    names = catalog.get("basis_names", [])
    if not isinstance(names, list):
        msg = f"Invalid catalog: basis_names is not a list in {_CATALOG_FILENAME}"
        raise TypeError(msg)
    return frozenset(str(name) for name in names)


def parse_baslib_basis_names(path: Path) -> frozenset[str]:
    """Parse ``A-`` / ``O-`` / ``P-`` basis names from a StoBe ``baslib*.new7`` file.

    Parameters
    ----------
    path
        Path to ``baslib.new7`` or ``baslibD.new7``.

    Returns
    -------
    frozenset[str]
        Unique basis names in file order (stored sorted in the returned set).

    Raises
    ------
    FileNotFoundError
        When ``path`` is missing.
    """
    path = Path(path)
    if not path.is_file():
        msg = f"StoBe basis library not found: {path}"
        raise FileNotFoundError(msg)
    names: set[str] = set()
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if _BASIS_NAME.match(stripped):
            names.add(stripped)
    return frozenset(names)


def missing_basis_names(
    names: Iterable[str],
    *,
    catalog: Iterable[str] | None = None,
) -> list[str]:
    """List basis strings that are absent from a catalog.

    Parameters
    ----------
    names
        Basis names to check (empty strings are ignored).
    catalog
        Allowed names; defaults to :func:`bundled_basis_names`.

    Returns
    -------
    list[str]
        Missing non-empty names, sorted uniquely.
    """
    allowed = frozenset(catalog) if catalog is not None else bundled_basis_names()
    missing = {name.strip() for name in names if name and name.strip()}
    missing -= allowed
    return sorted(missing)


def validate_basis_preset(
    preset: StoBeBasisPreset,
    *,
    catalog: Iterable[str] | None = None,
) -> list[str]:
    """Validate non-empty strings on a ``StoBeBasisPreset``.

    Parameters
    ----------
    preset
        Element basis preset from ``stobe_basis_for_element``.
    catalog
        Allowed basis names; defaults to the bundled catalog.

    Returns
    -------
    list[str]
        Missing basis names referenced by ``preset``.
    """
    fields = (
        preset.aux_ground,
        preset.aux_excited,
        preset.orbital_ground,
        preset.orbital_excited,
        preset.mcp,
    )
    return missing_basis_names(name for name in fields if name)


def basis_strings_from_molconfig_source(source: str) -> list[str]:
    """Extract quoted basis assignments from ``molConfig.py`` source text."""
    return [
        match.group(1).strip()
        for match in _MOLCONFIG_BASIS.finditer(source)
        if match.group(1).strip()
    ]


def validate_molconfig_basis(
    source: str,
    *,
    catalog: Iterable[str] | None = None,
) -> list[str]:
    """Validate basis strings embedded in a ``molConfig.py`` module.

    Parameters
    ----------
    source
        Full ``molConfig.py`` text.
    catalog
        Allowed basis names; defaults to the bundled catalog.

    Returns
    -------
    list[str]
        Missing basis names assigned in ``source``.
    """
    names = basis_strings_from_molconfig_source(source)
    return missing_basis_names(names, catalog=catalog)


def validate_molconfig_file(
    path: Path,
    *,
    catalog: Iterable[str] | None = None,
) -> list[str]:
    """Validate basis strings in a ``molConfig.py`` file on disk."""
    path = Path(path)
    if not path.is_file():
        msg = f"molConfig not found: {path}"
        raise FileNotFoundError(msg)
    return validate_molconfig_basis(path.read_text(encoding="utf-8"), catalog=catalog)


def validate_against_baslib_paths(
    names: Iterable[str],
    *,
    baslib: Path,
    baslib_d: Path | None = None,
) -> list[str]:
    """Validate names against live ``baslib.new7`` / ``baslibD.new7`` on disk.

    Parameters
    ----------
    names
        Basis names to verify.
    baslib
        Primary basis library path (for example ``/bin/stobe/Basis/baslib.new7``).
    baslib_d
        Optional diffuse/auxiliary library; defaults to ``baslib`` parent /
        ``baslibD.new7``.

    Returns
    -------
    list[str]
        Names not found in either library file.
    """
    baslib = Path(baslib)
    baslib_d_path = (
        Path(baslib_d) if baslib_d is not None else baslib.parent / "baslibD.new7"
    )
    allowed = parse_baslib_basis_names(baslib)
    if baslib_d_path.is_file():
        allowed = allowed | parse_baslib_basis_names(baslib_d_path)
    return missing_basis_names(names, catalog=allowed)
