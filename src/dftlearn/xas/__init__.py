"""Build broadened XAS from TP dipole sticks and Delta-KS energy alignment.

This package converts StoBe ``C*.xas`` sticks to a photon-energy spectrum using
the ``xrayspec.x`` Gaussian FWHM schedule (unshifted, matching ``XrayT*.out``),
then optionally applies the rigid interpolation
:math:`I'(E) = I(E - E^c)` with :math:`E^c = E^e - E^g - E_1`. Optional C3
folding of Cartesian dipoles in the Al-N/O molecular frame is available via
``collect_site_tp_xas(..., c3_symmetrize=True)``. Clustering and experiment
fitting live elsewhere.
"""

from __future__ import annotations

from dftlearn.xas.c3_symmetry import (
    C3Frame,
    build_al_n_o_c3_frame,
    build_c3_frame_from_xyz,
    fold_dipole_c3_tensor,
    fold_dipoles_c3,
    oscillator_strengths_from_tensors,
    write_c3_frame_json,
)
from dftlearn.xas.spectrum import (
    STOBE_XAS_HA_TO_EV,
    XAS_INTENSITY_SCALE,
    aligned_dipole_tensor_tables,
    collect_site_tp_xas,
    dipole_cartesian_oscillator_strengths,
    gaussian_xas_spectrum,
    padded_xas_energy_axis,
    piecewise_fwhm_ev,
    reconstruct_xas_spectrum,
    shift_and_resample_xas_spectrum,
    shift_xas_spectrum,
    xas_delta_ks_shift_ev,
    xas_energy_axis,
)
from dftlearn.xas.step_edge import (
    compound_mu,
    gaussian_step,
    gaussian_step_edge,
    parse_henke_nff,
)

__all__ = [
    "STOBE_XAS_HA_TO_EV",
    "XAS_INTENSITY_SCALE",
    "C3Frame",
    "aligned_dipole_tensor_tables",
    "build_al_n_o_c3_frame",
    "build_c3_frame_from_xyz",
    "collect_site_tp_xas",
    "compound_mu",
    "dipole_cartesian_oscillator_strengths",
    "fold_dipole_c3_tensor",
    "fold_dipoles_c3",
    "gaussian_step",
    "gaussian_step_edge",
    "gaussian_xas_spectrum",
    "oscillator_strengths_from_tensors",
    "padded_xas_energy_axis",
    "parse_henke_nff",
    "piecewise_fwhm_ev",
    "reconstruct_xas_spectrum",
    "shift_and_resample_xas_spectrum",
    "shift_xas_spectrum",
    "write_c3_frame_json",
    "xas_delta_ks_shift_ev",
    "xas_energy_axis",
]
