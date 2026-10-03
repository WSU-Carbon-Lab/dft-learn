"""Parsers and loaders for StoBe X-ray tables, TP sticks, and XYZ geometry."""

from __future__ import annotations

from dftlearn.io.stobe_final_energy import (
    HA_TO_EV,
    collect_delta_ks_site_table,
    collect_final_energies_long,
    enrich_final_energies_delta_vs_gnd,
    final_energy_site_summary,
    ionization_energies_table,
    parse_stobe_final_energy_tables,
    parse_stobe_tp_core_hole_orbital_ev,
    parse_stobe_tp_ionization_potential_ev,
    parse_stobe_tp_lumo_alpha_ev,
)
from dftlearn.io.stobe_scf_convergence import (
    collect_scf_convergence_long,
    discover_site_stobe_out,
    parse_stobe_scf_convergence_table,
    scf_convergence_auc_metrics,
    site_tag_from_stobe_out_filename,
)
from dftlearn.io.stobe_xas_sticks import (
    CARBON_XAS_SPEC,
    XasSpecSettings,
    parse_stobe_xas_dipole_sticks,
    parse_stobe_xas_inp,
    parse_stobe_xas_sticks,
    site_xas_inp_path,
    site_xas_stick_paths,
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
    "CARBON_XAS_SPEC",
    "HA_TO_EV",
    "XasSpecSettings",
    "collect_delta_ks_site_table",
    "collect_final_energies_long",
    "collect_scf_convergence_long",
    "collect_site_xray_spectra",
    "discover_site_stobe_out",
    "element_symbol_from_xyz_label",
    "enrich_final_energies_delta_vs_gnd",
    "final_energy_site_summary",
    "ionization_energies_table",
    "mol_from_xyz_file",
    "parse_stobe_final_energy_tables",
    "parse_stobe_scf_convergence_table",
    "parse_stobe_tp_core_hole_orbital_ev",
    "parse_stobe_tp_ionization_potential_ev",
    "parse_stobe_tp_lumo_alpha_ev",
    "parse_stobe_xas_dipole_sticks",
    "parse_stobe_xas_inp",
    "parse_stobe_xas_sticks",
    "parse_xray_out_table",
    "scf_convergence_auc_metrics",
    "site_label_to_atom_index",
    "site_label_to_atom_index_from_rows",
    "site_spectra_to_long_frame",
    "site_tag_from_stobe_out_filename",
    "site_xas_inp_path",
    "site_xas_stick_paths",
    "site_xray_paths",
    "xyz_rows_from_file",
]
