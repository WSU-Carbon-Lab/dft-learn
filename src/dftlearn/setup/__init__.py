"""Calculation initialization: CIF ingest, geometry, basis sets, and molConfig."""

from __future__ import annotations

from dftlearn.setup.basis_sets import stobe_basis_for_element, supported_basis_elements
from dftlearn.setup.cif_io import cif_blocks, load_cif_structure
from dftlearn.setup.distinct_atoms import distinct_groups_for_element
from dftlearn.setup.geometry import plan_stobe_geometry, write_stobe_xyz
from dftlearn.setup.molconfig import (
    UNASSIGNED_ALFA_OCC,
    render_molconfig_py,
    valence_electron_count,
    write_molconfig_py,
)
from dftlearn.setup.preview import write_init_preview_figure
from dftlearn.setup.remote import (
    KNOWN_REMOTE_HOSTS,
    remote_run_directory,
    resolve_ssh_host,
    sync_run_directory,
)
from dftlearn.setup.types import (
    CifAtomSite,
    CifStructureMeta,
    DistinctAtomGroup,
    InitArtifacts,
    InitGeometryPlan,
    StoBeBasisPreset,
)
from dftlearn.setup.workflow import finalize_setup_session, initialize_from_cif

__all__ = [
    "KNOWN_REMOTE_HOSTS",
    "UNASSIGNED_ALFA_OCC",
    "CifAtomSite",
    "CifStructureMeta",
    "DistinctAtomGroup",
    "InitArtifacts",
    "InitGeometryPlan",
    "StoBeBasisPreset",
    "cif_blocks",
    "distinct_groups_for_element",
    "finalize_setup_session",
    "initialize_from_cif",
    "load_cif_structure",
    "plan_stobe_geometry",
    "remote_run_directory",
    "render_molconfig_py",
    "resolve_ssh_host",
    "stobe_basis_for_element",
    "supported_basis_elements",
    "sync_run_directory",
    "valence_electron_count",
    "write_init_preview_figure",
    "write_molconfig_py",
    "write_stobe_xyz",
]
