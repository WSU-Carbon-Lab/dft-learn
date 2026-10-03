"""Gaussian XAS reconstruction from TP sticks and Delta-KS alignment shifts.

``xrayspec.x`` (StoBe 2013) evaluates a unit-normalized Gaussian at each dipole
stick, using FWHM as a function of that stick's photon energy, then scales the
sum by 1000. Hartree stick energies are converted with the deMon/StoBe factor
27.2116 eV. Broadening is never given a pre-shift: that spectrum is the
validation target against ``XrayT*.out`` (global shift 0). Delta-KS alignment
is a rigid interpolation of the already-broadened curve,
:math:`I'(E) = I(E - E^c)` with :math:`E^c = E^e - E^g - E_1`. Alignment
reconstructs on a padded energy axis so the shifted intensity remains defined
through the original ``RANGE`` window.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from natsort import natsorted

from dftlearn.io.stobe_final_energy import collect_delta_ks_site_table
from dftlearn.io.stobe_xas_sticks import (
    CARBON_XAS_SPEC,
    parse_stobe_xas_dipole_sticks,
    parse_stobe_xas_inp,
    parse_stobe_xas_sticks,
    site_xas_inp_path,
    site_xas_stick_paths,
)
from dftlearn.io.xray_out import parse_xray_out_table, site_xray_paths
from dftlearn.xas.c3_symmetry import (
    build_c3_frame_from_xyz,
    fold_dipoles_c3,
    oscillator_strengths_from_tensors,
)

if TYPE_CHECKING:
    from dftlearn.io.stobe_xas_sticks import XasSpecSettings
    from dftlearn.xas.c3_symmetry import C3Frame

STOBE_XAS_HA_TO_EV = 27.2116
XAS_INTENSITY_SCALE = 1000.0

_LN2 = float(np.log(2.0))
_FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * _LN2))


def xas_energy_axis(settings: XasSpecSettings) -> np.ndarray:
    """Build the inclusive ``RANGE`` grid used by ``xrayspec.x``.

    Parameters
    ----------
    settings : XasSpecSettings
        ``RANGE`` endpoints and ``POINTS`` count.

    Returns
    -------
    numpy.ndarray
        Shape ``(n_points,)``, ``float64`` energies in eV.
    """
    return np.linspace(
        settings.energy_min_ev,
        settings.energy_max_ev,
        settings.n_points,
        dtype=np.float64,
    )


def piecewise_fwhm_ev(
    energy_ev: np.ndarray,
    settings: XasSpecSettings,
) -> np.ndarray:
    """Evaluate the xrayspec piecewise-linear FWHM schedule at ``energy_ev``.

    FWHM is ``fwhm_low_ev`` at and below ``e_break_low_ev``, ``fwhm_high_ev`` at
    and above ``e_break_high_ev``, and linear between those breakpoints.

    Parameters
    ----------
    energy_ev : numpy.ndarray
        Photon energies (eV), any shape.
    settings : XasSpecSettings
        WIDTH card values.

    Returns
    -------
    numpy.ndarray
        FWHM in eV, same shape as ``energy_ev``, ``float64``.

    Raises
    ------
    ValueError
        If the WIDTH breakpoints are not strictly increasing.
    """
    span = settings.e_break_high_ev - settings.e_break_low_ev
    if span <= 0.0:
        msg = "WIDTH breakpoints must satisfy e_break_high_ev > e_break_low_ev"
        raise ValueError(msg)
    energy = np.asarray(energy_ev, dtype=np.float64)
    ramp = settings.fwhm_low_ev + (
        (settings.fwhm_high_ev - settings.fwhm_low_ev)
        * (energy - settings.e_break_low_ev)
        / span
    )
    return np.where(
        energy <= settings.e_break_low_ev,
        settings.fwhm_low_ev,
        np.where(energy >= settings.e_break_high_ev, settings.fwhm_high_ev, ramp),
    )


def gaussian_xas_spectrum(
    energy_ev: np.ndarray,
    stick_energy_ev: np.ndarray,
    oscillator_strength: np.ndarray,
    fwhm_ev: np.ndarray,
    *,
    intensity_scale: float = XAS_INTENSITY_SCALE,
) -> np.ndarray:
    """Sum unit-normalized Gaussians of dipole sticks onto ``energy_ev``.

    Each stick ``i`` contributes
    ``intensity_scale * OS_i / (sigma_i sqrt(2 pi)) *
    exp(-1/2 ((E - E_i)/sigma_i)^2)``
    with ``sigma_i = FWHM_i / (2 sqrt(2 ln 2))``. This is the xrayspec total
    XAS kernel (dipole only).

    Parameters
    ----------
    energy_ev : numpy.ndarray
        Spectral axis, shape ``(n,)``, eV.
    stick_energy_ev : numpy.ndarray
        Stick centers, shape ``(m,)``, eV.
    oscillator_strength : numpy.ndarray
        Dipole OS, shape ``(m,)``, dimensionless.
    fwhm_ev : numpy.ndarray
        Per-stick Gaussian FWHM, shape ``(m,)``, eV. Must be strictly positive.
    intensity_scale : float, optional
        Global multiplier (xrayspec uses 1000 for total XAS).

    Returns
    -------
    numpy.ndarray
        Shape ``(n,)``, ``float64`` intensity in xrayspec arb. units.

    Raises
    ------
    ValueError
        If stick arrays differ in length, FWHM is not positive, or an axis is
        empty.
    """
    energy = np.asarray(energy_ev, dtype=np.float64)
    centers = np.asarray(stick_energy_ev, dtype=np.float64)
    osc = np.asarray(oscillator_strength, dtype=np.float64)
    fwhm = np.asarray(fwhm_ev, dtype=np.float64)
    if energy.ndim != 1 or energy.size == 0:
        msg = "energy_ev must be a non-empty 1-d array"
        raise ValueError(msg)
    if centers.shape != osc.shape or centers.ndim != 1 or centers.size == 0:
        msg = (
            "stick_energy_ev and oscillator_strength must be non-empty "
            "1-d arrays of equal length"
        )
        raise ValueError(msg)
    if fwhm.shape != centers.shape:
        msg = "fwhm_ev must match stick_energy_ev shape"
        raise ValueError(msg)
    if np.any(fwhm <= 0.0):
        msg = "All FWHM values must be strictly positive"
        raise ValueError(msg)
    sigma = fwhm * _FWHM_TO_SIGMA
    delta = energy[:, np.newaxis] - centers[np.newaxis, :]
    gauss = np.exp(-0.5 * (delta / sigma[np.newaxis, :]) ** 2)
    amp = osc / (sigma * np.sqrt(2.0 * np.pi))
    return np.asarray(
        intensity_scale * (gauss * amp[np.newaxis, :]).sum(axis=1),
        dtype=np.float64,
    )


def xas_delta_ks_shift_ev(
    e_exc_ev: float,
    e_gnd_ev: float,
    e_first_stick_ev: float,
) -> float:
    """Return the rigid XAS alignment :math:`E^c = E^e - E^g - E_1`.

    ``E^e`` and ``E^g`` are SCF total energies in eV. ``E_1`` is the first TP
    dipole stick in eV. Applied after broadening as
    :math:`I'(E) = I(E - E^c)`, this places the first transition at the ΔKS
    1s→LUMO energy without changing the xrayspec FWHM schedule.

    Parameters
    ----------
    e_exc_ev : float
        Excited-state SCF total energy (eV).
    e_gnd_ev : float
        Ground-state SCF total energy (eV).
    e_first_stick_ev : float
        Lowest TP dipole transition energy (eV).

    Returns
    -------
    float
        Alignment shift in eV.
    """
    return float(e_exc_ev) - float(e_gnd_ev) - float(e_first_stick_ev)


def shift_xas_spectrum(
    energy_ev: np.ndarray,
    intensity: np.ndarray,
    shift_ev: float,
) -> np.ndarray:
    """Translate an already-broadened spectrum by ``shift_ev`` on a fixed grid.

    Returns :math:`I(E - E^c)` interpolated onto ``energy_ev``, with zeros
    outside the original support. This is the Igor ``DFTinterp`` rigid shift:
    the xrayspec kernel is left unchanged. For a window-filling tail after a
    negative :math:`E^c`, reconstruct on :func:`padded_xas_energy_axis` first
    (as :func:`collect_site_tp_xas` does) so the interpolant is not clipped.

    Parameters
    ----------
    energy_ev : numpy.ndarray
        Monotonic spectral axis (eV), shape ``(n,)``.
    intensity : numpy.ndarray
        Broadened intensity on ``energy_ev``, shape ``(n,)``.
    shift_ev : float
        Alignment :math:`E^c` in eV. Positive values move features to higher
        photon energy.

    Returns
    -------
    numpy.ndarray
        Shape ``(n,)``, ``float64`` shifted intensity.

    Raises
    ------
    ValueError
        If ``energy_ev`` and ``intensity`` differ in shape, are not 1-d, or
        ``energy_ev`` is not strictly increasing.
    """
    energy = np.asarray(energy_ev, dtype=np.float64)
    y = np.asarray(intensity, dtype=np.float64)
    if energy.shape != y.shape or energy.ndim != 1 or energy.size == 0:
        msg = "energy_ev and intensity must be non-empty 1-d arrays of equal shape"
        raise ValueError(msg)
    if np.any(np.diff(energy) <= 0.0):
        msg = "energy_ev must be strictly increasing"
        raise ValueError(msg)
    sample_at = energy - float(shift_ev)
    return np.interp(sample_at, energy, y, left=0.0, right=0.0).astype(
        np.float64, copy=False
    )


def padded_xas_energy_axis(
    energy_ev: np.ndarray,
    shift_ev: float,
    *,
    extra_ev: float = 2.0,
) -> np.ndarray:
    """Extend a uniform energy grid so :math:`I(E - E^c)` is defined on it.

    The returned axis keeps the original spacing and covers both the native
    window and the window shifted by ``shift_ev``, plus ``extra_ev`` of
    padding on each end.

    Parameters
    ----------
    energy_ev : numpy.ndarray
        Strictly increasing display axis (eV), shape ``(n,)`` with ``n >= 2``.
    shift_ev : float
        Alignment :math:`E^c` in eV.
    extra_ev : float, optional
        Extra padding beyond the required support (eV).

    Returns
    -------
    numpy.ndarray
        Padded ``float64`` axis, strictly increasing.

    Raises
    ------
    ValueError
        If ``energy_ev`` is too short or not strictly increasing.
    """
    energy = np.asarray(energy_ev, dtype=np.float64)
    if energy.ndim != 1 or energy.size < 2:
        msg = "energy_ev must be a 1-d array with at least two points"
        raise ValueError(msg)
    deltas = np.diff(energy)
    if np.any(deltas <= 0.0):
        msg = "energy_ev must be strictly increasing"
        raise ValueError(msg)
    de = float(np.median(deltas))
    e0 = float(energy[0])
    e1 = float(energy[-1])
    shift = float(shift_ev)
    lo = min(e0, e0 - shift) - float(extra_ev)
    hi = max(e1, e1 - shift) + float(extra_ev)
    n = int(np.floor((hi - lo) / de)) + 1
    return lo + de * np.arange(n, dtype=np.float64)


def dipole_cartesian_oscillator_strengths(
    energy_ha: np.ndarray,
    dipole_au: np.ndarray,
) -> np.ndarray:
    r"""Convert Cartesian dipole matrix elements to StoBe oscillator strengths.

    Uses :math:`f_{ii} = (2/3) E_\\mathrm{Ha} \\mu_i^2` so that
    :math:`f_{xx}+f_{yy}+f_{zz}` matches the total dipole OS in ``C*.xas``.

    Parameters
    ----------
    energy_ha : numpy.ndarray
        Transition energies in Hartree, shape ``(m,)``.
    dipole_au : numpy.ndarray
        Dipole matrix elements :math:`(\\mu_x, \\mu_y, \\mu_z)`, shape ``(m, 3)``.

    Returns
    -------
    numpy.ndarray
        Shape ``(m, 3)`` with columns ``(f_xx, f_yy, f_zz)``.

    Raises
    ------
    ValueError
        If shapes are incompatible.
    """
    energy = np.asarray(energy_ha, dtype=np.float64)
    dipole = np.asarray(dipole_au, dtype=np.float64)
    if energy.ndim != 1 or dipole.ndim != 2 or dipole.shape != (energy.size, 3):
        msg = "energy_ha must be (m,) and dipole_au (m, 3)"
        raise ValueError(msg)
    scale = (2.0 / 3.0) * energy[:, np.newaxis]
    return np.asarray(scale * dipole**2, dtype=np.float64)


def shift_and_resample_xas_spectrum(
    energy_src: np.ndarray,
    intensity_src: np.ndarray,
    energy_dst: np.ndarray,
    shift_ev: float,
) -> np.ndarray:
    """Shift a spectrum on ``energy_src`` and interpolate onto ``energy_dst``.

    Parameters
    ----------
    energy_src : numpy.ndarray
        Source axis used for broadening (eV).
    intensity_src : numpy.ndarray
        Intensity on ``energy_src``.
    energy_dst : numpy.ndarray
        Display axis (eV).
    shift_ev : float
        Alignment :math:`E^c` in eV.

    Returns
    -------
    numpy.ndarray
        Intensity on ``energy_dst``.
    """
    shifted = shift_xas_spectrum(energy_src, intensity_src, shift_ev)
    dst = np.asarray(energy_dst, dtype=np.float64)
    src = np.asarray(energy_src, dtype=np.float64)
    return np.interp(dst, src, shifted, left=0.0, right=0.0).astype(
        np.float64, copy=False
    )


def _align_reconstructed_spectrum(
    sticks_ha: np.ndarray,
    oscillator_strength: np.ndarray,
    energy_ev: np.ndarray,
    settings: XasSpecSettings,
    shift_ev: float,
    *,
    ha_to_ev: float,
    intensity_scale: float,
) -> np.ndarray:
    """Broaden sticks on a padded axis and resample :math:`I(E - E^c)`."""
    pair = np.column_stack(
        (
            np.asarray(sticks_ha[:, 0], dtype=np.float64),
            np.asarray(oscillator_strength, dtype=np.float64),
        )
    )
    pad = padded_xas_energy_axis(energy_ev, shift_ev)
    y_pad = reconstruct_xas_spectrum(
        pair,
        pad,
        settings,
        ha_to_ev=ha_to_ev,
        intensity_scale=intensity_scale,
    )
    return shift_and_resample_xas_spectrum(pad, y_pad, energy_ev, shift_ev)


def _site_stick_arrays(
    path: Path,
    *,
    c3_frame: C3Frame | None = None,
) -> tuple[np.ndarray, np.ndarray | None]:
    """Return ``(m, 2)`` total-OS sticks and optional Cartesian OS ``(m, 3)``."""
    try:
        dipole = parse_stobe_xas_dipole_sticks(path)
    except ValueError:
        return parse_stobe_xas_sticks(path), None
    if c3_frame is not None:
        tensors = fold_dipoles_c3(dipole[:, 2:5], c3_frame.rotation)
        cart = oscillator_strengths_from_tensors(dipole[:, 0], tensors)
    else:
        cart = dipole_cartesian_oscillator_strengths(dipole[:, 0], dipole[:, 2:5])
    return dipole[:, :2], cart


def reconstruct_xas_spectrum(
    sticks_ha: np.ndarray,
    energy_ev: np.ndarray,
    settings: XasSpecSettings,
    *,
    ha_to_ev: float = STOBE_XAS_HA_TO_EV,
    intensity_scale: float = XAS_INTENSITY_SCALE,
) -> np.ndarray:
    """Convert Hartree dipole sticks to a broadened spectrum on ``energy_ev``.

    Stick energies are ``energy_ha * ha_to_ev``. FWHM is evaluated at each
    unshifted stick energy so the result matches ``xrayspec.x``.

    Parameters
    ----------
    sticks_ha : numpy.ndarray
        Shape ``(m, 2)`` from :func:`parse_stobe_xas_sticks`.
    energy_ev : numpy.ndarray
        Spectral axis (eV), shape ``(n,)``.
    settings : XasSpecSettings
        xrayspec WIDTH schedule.
    ha_to_ev : float, optional
        Hartree to eV for fort.11 energies (StoBe uses 27.2116).
    intensity_scale : float, optional
        Passed to :func:`gaussian_xas_spectrum`.

    Returns
    -------
    numpy.ndarray
        Shape ``(n,)`` broadened intensity.

    Raises
    ------
    ValueError
        If ``sticks_ha`` is not ``(m, 2)`` with ``m >= 1``.
    """
    sticks = np.asarray(sticks_ha, dtype=np.float64)
    if sticks.ndim != 2 or sticks.shape[1] != 2 or sticks.shape[0] < 1:
        msg = "sticks_ha must have shape (m, 2) with m >= 1"
        raise ValueError(msg)
    stick_ev = sticks[:, 0] * float(ha_to_ev)
    fwhm = piecewise_fwhm_ev(stick_ev, settings)
    return gaussian_xas_spectrum(
        energy_ev,
        stick_ev,
        sticks[:, 1],
        fwhm,
        intensity_scale=intensity_scale,
    )


def collect_site_tp_xas(
    run_root: Path,
    *,
    xray_filename: str = "XrayT001.out",
    ha_to_ev: float = STOBE_XAS_HA_TO_EV,
    intensity_scale: float = XAS_INTENSITY_SCALE,
    xyz_path: Path | None = None,
    c3_symmetrize: bool = False,
) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame, pd.DataFrame, C3Frame | None]:
    r"""Reconstruct per-site XAS from TP sticks and compare to ``XrayT*.out``.

    Broadening matches ``xrayspec.x`` with no energy shift. When ground and
    excited SCF totals exist, :math:`E^c` is applied as a rigid shift of a
    spectrum reconstructed on a padded energy axis so the aligned curve does
    not clip at the high-energy end of ``RANGE``. Cartesian dipole components
    are broadened the same way when ``C*.xas`` provides :math:`\\mu_x,\\mu_y,\\mu_z`.
    Stick tables keep both the unshifted and aligned transition energies.

    When ``c3_symmetrize`` is True, dipoles are rotated into the Al-N/O
    coordination triangle frame and folded under C3 about molecular z before
    Cartesian oscillator strengths are formed.

    Parameters
    ----------
    run_root : pathlib.Path
        StoBe run root with site folders containing ``{site}.xas``.
    xray_filename : str, optional
        StoBe broadened table name used as the validation axis and reference.
    ha_to_ev : float, optional
        fort.11 Hartree conversion.
    intensity_scale : float, optional
        xrayspec intensity scale.
    xyz_path : pathlib.Path, optional
        Geometry used to build the C3 frame when ``c3_symmetrize`` is True.
    c3_symmetrize : bool, optional
        When True, fold Cartesian OS under C3 in the Al-N/O frame.

    Returns
    -------
    energy : numpy.ndarray
        Common photon-energy axis (eV).
    spectra : pandas.DataFrame
        Long table with ``energy_ev``, ``abs_stobe``, ``abs_tp``,
        ``abs_aligned``, ``abs_xx``, ``abs_yy``, ``abs_zz``,
        ``abs_xx_aligned``, ``abs_yy_aligned``, ``abs_zz_aligned``, and
        ``site``. StoBe, aligned, and tensor columns are NaN when the
        corresponding input is missing.
    metrics : pandas.DataFrame
        One row per site: stick count, :math:`E_1`, :math:`E^c`, RMSE and max
        absolute residual vs StoBe, and source paths.
    sticks : pandas.DataFrame
        One row per dipole transition: ``site``, ``energy_ev``,
        ``energy_aligned_ev``, ``oscillator_strength``, and Cartesian
        ``os_xx``, ``os_yy``, ``os_zz`` when dipole components are present.
    c3_frame : C3Frame or None
        Molecular frame used for folding, or ``None`` when not symmetrizing.

    Raises
    ------
    FileNotFoundError
        If no ``{site}.xas`` files are found.
    ValueError
        If a stick file cannot be parsed, site energy axes disagree, or
        ``c3_symmetrize`` is True without a usable ``xyz_path``.
    """
    run_root = Path(run_root).resolve()
    stick_pairs = site_xas_stick_paths(run_root)
    xray_map = dict(_optional_xray_paths(run_root, xray_filename))
    delta_ks = _delta_ks_by_site(run_root)

    c3_frame: C3Frame | None = None
    if c3_symmetrize:
        if xyz_path is None:
            msg = "c3_symmetrize requires xyz_path for the Al-N/O coordination frame"
            raise ValueError(msg)
        c3_frame = build_c3_frame_from_xyz(Path(xyz_path))

    energy: np.ndarray | None = None
    spec_rows: list[pd.DataFrame] = []
    metric_rows: list[dict[str, object]] = []
    stick_rows: list[pd.DataFrame] = []

    for site, stick_path in stick_pairs:
        sticks, cart = _site_stick_arrays(stick_path, c3_frame=c3_frame)
        settings = _settings_for_site(stick_path.parent, site)
        xray_path = xray_map.get(site)
        if xray_path is not None:
            table = parse_xray_out_table(xray_path)
            site_energy = table[:, 0].copy()
            abs_stobe = table[:, 1].copy()
        else:
            site_energy = xas_energy_axis(settings)
            abs_stobe = np.full(site_energy.shape, np.nan, dtype=np.float64)
        if energy is None:
            energy = site_energy
        elif site_energy.shape != energy.shape or not np.allclose(
            site_energy, energy, rtol=0.0, atol=1e-6
        ):
            msg = f"Energy axis mismatch for site {site!r}"
            raise ValueError(msg)
        abs_tp = reconstruct_xas_spectrum(
            sticks,
            site_energy,
            settings,
            ha_to_ev=ha_to_ev,
            intensity_scale=intensity_scale,
        )
        if cart is None:
            abs_xx = np.full(site_energy.shape, np.nan, dtype=np.float64)
            abs_yy = np.full(site_energy.shape, np.nan, dtype=np.float64)
            abs_zz = np.full(site_energy.shape, np.nan, dtype=np.float64)
        else:
            abs_xx = reconstruct_xas_spectrum(
                np.column_stack((sticks[:, 0], cart[:, 0])),
                site_energy,
                settings,
                ha_to_ev=ha_to_ev,
                intensity_scale=intensity_scale,
            )
            abs_yy = reconstruct_xas_spectrum(
                np.column_stack((sticks[:, 0], cart[:, 1])),
                site_energy,
                settings,
                ha_to_ev=ha_to_ev,
                intensity_scale=intensity_scale,
            )
            abs_zz = reconstruct_xas_spectrum(
                np.column_stack((sticks[:, 0], cart[:, 2])),
                site_energy,
                settings,
                ha_to_ev=ha_to_ev,
                intensity_scale=intensity_scale,
            )
        e_first = float(np.min(sticks[:, 0]) * ha_to_ev)
        e_c = _shift_for_site(delta_ks, site, e_first)
        nan_spec = np.full(site_energy.shape, np.nan, dtype=np.float64)
        if np.isfinite(e_c):
            shift = float(e_c)
            abs_aligned = _align_reconstructed_spectrum(
                sticks,
                sticks[:, 1],
                site_energy,
                settings,
                shift,
                ha_to_ev=ha_to_ev,
                intensity_scale=intensity_scale,
            )
            stick_aligned = sticks[:, 0] * ha_to_ev + shift
            if cart is None:
                abs_xx_al = nan_spec.copy()
                abs_yy_al = nan_spec.copy()
                abs_zz_al = nan_spec.copy()
            else:
                abs_xx_al = _align_reconstructed_spectrum(
                    sticks,
                    cart[:, 0],
                    site_energy,
                    settings,
                    shift,
                    ha_to_ev=ha_to_ev,
                    intensity_scale=intensity_scale,
                )
                abs_yy_al = _align_reconstructed_spectrum(
                    sticks,
                    cart[:, 1],
                    site_energy,
                    settings,
                    shift,
                    ha_to_ev=ha_to_ev,
                    intensity_scale=intensity_scale,
                )
                abs_zz_al = _align_reconstructed_spectrum(
                    sticks,
                    cart[:, 2],
                    site_energy,
                    settings,
                    shift,
                    ha_to_ev=ha_to_ev,
                    intensity_scale=intensity_scale,
                )
        else:
            abs_aligned = nan_spec.copy()
            abs_xx_al = nan_spec.copy()
            abs_yy_al = nan_spec.copy()
            abs_zz_al = nan_spec.copy()
            stick_aligned = np.full(sticks.shape[0], np.nan, dtype=np.float64)
        spec_rows.append(
            pd.DataFrame(
                {
                    "energy_ev": site_energy,
                    "abs_stobe": abs_stobe,
                    "abs_tp": abs_tp,
                    "abs_aligned": abs_aligned,
                    "abs_xx": abs_xx,
                    "abs_yy": abs_yy,
                    "abs_zz": abs_zz,
                    "abs_xx_aligned": abs_xx_al,
                    "abs_yy_aligned": abs_yy_al,
                    "abs_zz_aligned": abs_zz_al,
                    "site": site,
                }
            )
        )
        if cart is None:
            os_xx = np.full(sticks.shape[0], np.nan, dtype=np.float64)
            os_yy = np.full(sticks.shape[0], np.nan, dtype=np.float64)
            os_zz = np.full(sticks.shape[0], np.nan, dtype=np.float64)
        else:
            os_xx = cart[:, 0]
            os_yy = cart[:, 1]
            os_zz = cart[:, 2]
        stick_rows.append(
            pd.DataFrame(
                {
                    "site": site,
                    "energy_ev": sticks[:, 0] * ha_to_ev,
                    "energy_aligned_ev": stick_aligned,
                    "oscillator_strength": sticks[:, 1],
                    "os_xx": os_xx,
                    "os_yy": os_yy,
                    "os_zz": os_zz,
                }
            )
        )
        finite_ref = np.isfinite(abs_stobe)
        if np.any(finite_ref):
            delta = abs_tp[finite_ref] - abs_stobe[finite_ref]
            rmse = float(np.sqrt(np.mean(delta**2)))
            max_abs = float(np.max(np.abs(delta)))
            corr = float(np.corrcoef(abs_tp[finite_ref], abs_stobe[finite_ref])[0, 1])
        else:
            rmse = float("nan")
            max_abs = float("nan")
            corr = float("nan")
        inp_path = site_xas_inp_path(stick_path.parent, site)
        metric_rows.append(
            {
                "site": site,
                "n_sticks": int(sticks.shape[0]),
                "e_first_stick_ev": e_first,
                "E_c_ev": e_c,
                "rmse_vs_stobe": rmse,
                "max_abs_vs_stobe": max_abs,
                "corr_vs_stobe": corr,
                "source_xas": str(stick_path),
                "source_xray": str(xray_path) if xray_path is not None else "",
                "source_inp": str(inp_path) if inp_path is not None else "",
            }
        )

    if energy is None or not spec_rows:
        msg = f"No TP stick spectra assembled under {run_root}"
        raise ValueError(msg)
    spectra = pd.concat(spec_rows, ignore_index=True)
    metrics = pd.DataFrame(metric_rows)
    site_order = natsorted(metrics["site"].unique())
    metrics = metrics.set_index("site").loc[site_order].reset_index()
    sticks_df = pd.concat(stick_rows, ignore_index=True)
    sticks_df = sticks_df.set_index("site").loc[site_order].reset_index()
    return energy, spectra, metrics, sticks_df, c3_frame


def aligned_dipole_tensor_tables(
    spectra: pd.DataFrame,
    metrics: pd.DataFrame,
    *,
    n_mean_points: int = 2000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    r"""Build Delta-KS aligned Cartesian dipole XAS tables from reconstructed spectra.

    Photon energy is the native xrayspec axis translated by each site's
    :math:`E^c` (the same x-shift used in the summary figure). Intensities
    :math:`I_{xx}, I_{yy}, I_{zz}` are the unshifted Cartesian reconstructions;
    :math:`I_\mathrm{iso} = (I_{xx}+I_{yy}+I_{zz})/3`. The mean table
    interpolates each site onto a common grid spanning the union of shifted
    windows and averages with NaNs outside each site's support.

    Parameters
    ----------
    spectra : pandas.DataFrame
        Long table from :func:`collect_site_tp_xas` with ``energy_ev``,
        ``abs_xx``, ``abs_yy``, ``abs_zz``, and ``site``.
    metrics : pandas.DataFrame
        One row per site with ``site`` and ``E_c_ev``.
    n_mean_points : int, optional
        Number of samples on the site-mean energy grid. Must be >= 2.

    Returns
    -------
    per_site : pandas.DataFrame
        Columns ``energy_ev``, ``I_xx``, ``I_yy``, ``I_zz``, ``I_iso``,
        ``site``, ``E_c_ev``.
    site_mean : pandas.DataFrame
        Columns ``energy_ev``, ``I_xx``, ``I_yy``, ``I_zz``, ``I_iso``.

    Raises
    ------
    ValueError
        If ``n_mean_points`` is less than 2, required columns are missing, or
        no sites can be assembled.
    """
    required = {"energy_ev", "abs_xx", "abs_yy", "abs_zz", "site"}
    missing = required.difference(spectra.columns)
    if missing:
        msg = f"spectra missing columns {sorted(missing)}"
        raise ValueError(msg)
    if n_mean_points < 2:
        msg = "n_mean_points must be >= 2"
        raise ValueError(msg)
    if "site" not in metrics.columns:
        msg = "metrics must contain a site column"
        raise ValueError(msg)
    e_c_map: dict[str, float] = {}
    if "E_c_ev" in metrics.columns:
        sites_m = metrics["site"].astype(str).to_numpy()
        vals = metrics["E_c_ev"].to_numpy(dtype=np.float64)
        e_c_map = {
            str(site): float(shift) for site, shift in zip(sites_m, vals, strict=True)
        }
    sites = natsorted(spectra["site"].astype(str).unique())
    if not sites:
        msg = "spectra has no site rows"
        raise ValueError(msg)
    per_rows: list[pd.DataFrame] = []
    shifted_lo: list[float] = []
    shifted_hi: list[float] = []
    for site in sites:
        chunk = spectra.loc[spectra["site"] == site]
        energy = chunk["energy_ev"].to_numpy(dtype=np.float64)
        shift = e_c_map.get(site, float("nan"))
        shift_use = float(shift) if np.isfinite(shift) else 0.0
        xx = chunk["abs_xx"].to_numpy(dtype=np.float64)
        yy = chunk["abs_yy"].to_numpy(dtype=np.float64)
        zz = chunk["abs_zz"].to_numpy(dtype=np.float64)
        iso = (xx + yy + zz) / 3.0
        e_shift = energy + shift_use
        shifted_lo.append(float(np.nanmin(e_shift)))
        shifted_hi.append(float(np.nanmax(e_shift)))
        per_rows.append(
            pd.DataFrame(
                {
                    "energy_ev": e_shift,
                    "I_xx": xx,
                    "I_yy": yy,
                    "I_zz": zz,
                    "I_iso": iso,
                    "site": site,
                    "E_c_ev": shift,
                }
            )
        )
    per_site = pd.concat(per_rows, ignore_index=True)
    grid = np.linspace(
        min(shifted_lo),
        max(shifted_hi),
        n_mean_points,
        dtype=np.float64,
    )
    stacked: dict[str, list[np.ndarray]] = {key: [] for key in ("I_xx", "I_yy", "I_zz")}
    for site in sites:
        chunk = per_site.loc[per_site["site"] == site]
        x_s = chunk["energy_ev"].to_numpy(dtype=np.float64)
        for key in stacked:
            y_s = chunk[key].to_numpy(dtype=np.float64)
            finite = np.isfinite(x_s) & np.isfinite(y_s)
            if not np.any(finite):
                stacked[key].append(np.full(grid.shape, np.nan, dtype=np.float64))
                continue
            stacked[key].append(
                np.interp(grid, x_s[finite], y_s[finite], left=np.nan, right=np.nan)
            )
    mean_xx = np.nanmean(np.vstack(stacked["I_xx"]), axis=0)
    mean_yy = np.nanmean(np.vstack(stacked["I_yy"]), axis=0)
    mean_zz = np.nanmean(np.vstack(stacked["I_zz"]), axis=0)
    site_mean = pd.DataFrame(
        {
            "energy_ev": grid,
            "I_xx": mean_xx,
            "I_yy": mean_yy,
            "I_zz": mean_zz,
            "I_iso": (mean_xx + mean_yy + mean_zz) / 3.0,
        }
    )
    return per_site, site_mean


def _optional_xray_paths(
    run_root: Path,
    xray_filename: str,
) -> list[tuple[str, Path]]:
    """Return site XrayT paths, or an empty list when none exist."""
    try:
        return site_xray_paths(run_root, xray_filename=xray_filename)
    except FileNotFoundError:
        return []


def _settings_for_site(site_dir: Path, site: str) -> XasSpecSettings:
    """Load ``{site}xas.inp`` when present; otherwise carbon xrayspec defaults."""
    inp = site_xas_inp_path(site_dir, site)
    if inp is None:
        return CARBON_XAS_SPEC
    return parse_stobe_xas_inp(inp)


def _delta_ks_by_site(run_root: Path) -> pd.DataFrame | None:
    """Load per-site SCF totals, or ``None`` when no FINAL ENERGY blocks exist."""
    try:
        wide = collect_delta_ks_site_table(run_root)
    except ValueError:
        return None
    if wide.empty:
        return None
    return wide.set_index("site")


def _shift_for_site(
    delta_ks: pd.DataFrame | None,
    site: str,
    e_first_stick_ev: float,
) -> float:
    """Return :math:`E^c` for ``site``, or NaN when SCF totals are missing."""
    if delta_ks is None or site not in delta_ks.index:
        return float("nan")
    row = delta_ks.loc[site]
    e_g = float(row["E_g_ev"])
    e_e = float(row["E_e_ev"])
    if not (np.isfinite(e_g) and np.isfinite(e_e)):
        return float("nan")
    return xas_delta_ks_shift_ev(e_e, e_g, e_first_stick_ev)
