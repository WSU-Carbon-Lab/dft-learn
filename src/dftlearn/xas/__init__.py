"""Build broadened XAS from TP dipole sticks and Delta-KS energy alignment.

This package converts StoBe ``C*.xas`` sticks to a photon-energy spectrum using
the ``xrayspec.x`` Gaussian FWHM schedule (unshifted, matching ``XrayT*.out``),
then optionally applies the rigid interpolation
:math:`I'(E) = I(E - E^c)` with :math:`E^c = E^e - E^g - E_1`. It does not
cluster peaks or fit experiment.
"""

from __future__ import annotations

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

__all__ = [
    "STOBE_XAS_HA_TO_EV",
    "XAS_INTENSITY_SCALE",
    "aligned_dipole_tensor_tables",
    "collect_site_tp_xas",
    "dipole_cartesian_oscillator_strengths",
    "gaussian_xas_spectrum",
    "padded_xas_energy_axis",
    "piecewise_fwhm_ev",
    "reconstruct_xas_spectrum",
    "shift_and_resample_xas_spectrum",
    "shift_xas_spectrum",
    "xas_delta_ks_shift_ev",
    "xas_energy_axis",
]
