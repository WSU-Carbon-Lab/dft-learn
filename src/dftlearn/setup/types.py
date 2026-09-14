"""Datatypes for StoBe calculation initialization from crystal structures."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path


@dataclass(frozen=True)
class CifAtomSite:
    """One atom site parsed from a CIF block."""

    label: str
    element: str
    x: float
    y: float
    z: float


@dataclass(frozen=True)
class CifStructureMeta:
    """Crystal metadata and Cartesian sites extracted from a CIF block."""

    block_name: str
    space_group: str | None
    cell_a: float
    cell_b: float
    cell_c: float
    cell_alpha: float
    cell_beta: float
    cell_gamma: float
    sites: tuple[CifAtomSite, ...]
    symmetry_ops: tuple[str, ...] = ()


@dataclass(frozen=True)
class StoBeBasisPreset:
    """StoBe auxiliary, orbital, and optional MCP basis strings for one element."""

    atomic_number: int
    z_eff_ground: int
    z_eff_excited: int
    aux_ground: str
    aux_excited: str
    orbital_ground: str
    orbital_excited: str
    mcp: str | None = None


@dataclass(frozen=True)
class DistinctAtomGroup:
    """Chemically equivalent atoms sharing one canonical RDKit rank."""

    rank: int
    element: str
    atom_indices: tuple[int, ...]
    representative_index: int
    cif_label: str
    xyz_label: str


@dataclass(frozen=True)
class InitGeometryPlan:
    """Ordered StoBe XYZ rows and core-site indices after initialization."""

    rows: tuple[tuple[str, float, float, float], ...]
    edge_element: str
    generating_site_indices: tuple[int, ...]
    distinct_groups: tuple[DistinctAtomGroup, ...]
    element_counts: dict[str, int]
    element_group_order: tuple[str, ...]


@dataclass
class InitArtifacts:
    """Paths and summary data written by calculation initialization."""

    run_directory: Path
    geometry_xyz: Path
    molconfig_py: Path
    source_cif: Path
    preview_png: Path | None
    summary_json: Path
    remote_host: str | None = None
    remote_directory: Path | None = None
    distinct_carbon_groups: list[DistinctAtomGroup] = field(default_factory=list)
