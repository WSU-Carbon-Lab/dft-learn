"""Parsers and loaders for calculation outputs (StoBe X-ray tables, XYZ geometry)."""

from __future__ import annotations

from dftlearn.io.stobe_scf_convergence import (
    collect_scf_convergence_long,
    discover_site_stobe_out,
    parse_stobe_scf_convergence_table,
    scf_convergence_auc_metrics,
    site_tag_from_stobe_out_filename,
)
from dftlearn.io.xray_out import (
    collect_site_xray_spectra,
    parse_xray_out_table,
    site_spectra_to_long_frame,
    site_xray_paths,
)
from dftlearn.io.xyz_structure import (
    element_symbol_from_xyz_label,
    mol_from_xyz_file,
    site_label_to_atom_index,
    site_label_to_atom_index_from_rows,
    xyz_rows_from_file,
)

__all__ = [
    "collect_scf_convergence_long",
    "collect_site_xray_spectra",
    "discover_site_stobe_out",
    "element_symbol_from_xyz_label",
    "mol_from_xyz_file",
    "parse_stobe_scf_convergence_table",
    "parse_xray_out_table",
    "scf_convergence_auc_metrics",
    "site_label_to_atom_index",
    "site_label_to_atom_index_from_rows",
    "site_spectra_to_long_frame",
    "site_tag_from_stobe_out_filename",
    "site_xray_paths",
    "xyz_rows_from_file",
]
