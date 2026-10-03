"""Tests for TP stick parsing, xrayspec broadening, and XAS reconstruction."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from dftlearn.io.stobe_final_energy import HA_TO_EV
from dftlearn.io.stobe_xas_sticks import (
    CARBON_XAS_SPEC,
    parse_stobe_xas_dipole_sticks,
    parse_stobe_xas_inp,
    parse_stobe_xas_sticks,
)
from dftlearn.io.xray_out import parse_xray_out_table
from dftlearn.visualization.xas_reconstruction_figure import (
    write_xas_reconstruction_report,
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
)

_EXAMPLE = (
    Path(__file__).resolve().parents[1] / "docs" / "stobe" / "example-run" / "znpc-cif"
)

_FINAL_BLOCK = """
 SCF CONVERGED AFTER   3 ITERATIONS

 FINAL ENERGY / CHARGE / GEOMETRY RESULTS :

 Total energy   (H) =   {energy:.10f}  (incl. numerical value for EXC)
 Nuc-nuc energy (H) =    1.0
 El-nuc energy  (H) =   -2.0
 Kinetic energy (H) =    1.0
 Coulomb energy (H) =    1.0
 Ex-cor energy  (H) =   -0.1

 Decomposition of exchange / correlation :
"""


def test_parse_stobe_xas_sticks_header_and_fortran_d(tmp_path: Path) -> None:
    path = tmp_path / "C1.xas"
    path.write_text(
        "XAS      2\n"
        "         10.5418380964    0.0063083964987634    0.1\n"
        "         10.5485412583    0.130912915578D-03    0.0\n",
        encoding="utf-8",
    )
    arr = parse_stobe_xas_sticks(path)
    assert arr.shape == (2, 2)
    np.testing.assert_allclose(arr[0, 0], 10.5418380964)
    np.testing.assert_allclose(arr[0, 1], 0.0063083964987634)
    np.testing.assert_allclose(arr[1, 1], 0.000130912915578)


def test_parse_stobe_xas_inp(tmp_path: Path) -> None:
    path = tmp_path / "C1xas.inp"
    path.write_text(
        "title\nznpc\nRANGE 280 320\nPOINTS 2000\nWIDTH 0.5 12 288 320\nEND\n",
        encoding="utf-8",
    )
    settings = parse_stobe_xas_inp(path)
    assert settings.energy_min_ev == 280.0
    assert settings.energy_max_ev == 320.0
    assert settings.n_points == 2000
    assert settings.fwhm_low_ev == 0.5
    assert settings.e_break_low_ev == 288.0


def test_piecewise_fwhm_schedule() -> None:
    energy = np.array([280.0, 288.0, 304.0, 320.0, 330.0])
    w = piecewise_fwhm_ev(energy, CARBON_XAS_SPEC)
    np.testing.assert_allclose(w[0], 0.5)
    np.testing.assert_allclose(w[1], 0.5)
    np.testing.assert_allclose(w[2], 6.25)
    np.testing.assert_allclose(w[3], 12.0)
    np.testing.assert_allclose(w[4], 12.0)


def test_gaussian_peak_height_matches_normalized_kernel() -> None:
    energy = np.linspace(280.0, 290.0, 1001)
    e_i = np.array([285.0])
    osc = np.array([0.01])
    fwhm = np.array([0.5])
    y = gaussian_xas_spectrum(energy, e_i, osc, fwhm, intensity_scale=1000.0)
    sigma = 0.5 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    expected = 1000.0 * 0.01 / (sigma * np.sqrt(2.0 * np.pi))
    np.testing.assert_allclose(
        y[np.argmin(np.abs(energy - 285.0))], expected, rtol=1e-12
    )


def test_delta_ks_shift_places_first_stick_on_excitation_energy() -> None:
    e_g = 0.0
    e_e = 290.0
    e_1 = 286.86
    shift = xas_delta_ks_shift_ev(e_e, e_g, e_1)
    np.testing.assert_allclose(e_1 + shift, e_e - e_g)


def test_shift_xas_spectrum_samples_unshifted_curve() -> None:
    energy = np.linspace(0.0, 10.0, 11)
    intensity = energy.copy()
    shifted = shift_xas_spectrum(energy, intensity, 2.0)
    np.testing.assert_allclose(shifted[4], 2.0)
    np.testing.assert_allclose(shifted[0], 0.0)


@pytest.mark.skipif(not _EXAMPLE.is_dir(), reason="example ZnPc run not present")
def test_reconstructed_tp_xas_matches_stobe_xrayt() -> None:
    energy, spectra, metrics, sticks, _frame = collect_site_tp_xas(_EXAMPLE)
    assert energy.shape[0] == 2000
    assert set(metrics["site"]) == {"C1", "C2", "C3", "C4"}
    assert set(sticks["site"].unique()) == {"C1", "C2", "C3", "C4"}
    assert np.all(metrics["corr_vs_stobe"].to_numpy(dtype=np.float64) > 0.9999)
    assert np.all(metrics["max_abs_vs_stobe"].to_numpy(dtype=np.float64) < 0.01)
    c1 = parse_stobe_xas_sticks(_EXAMPLE / "C1" / "C1.xas")
    ref = parse_xray_out_table(_EXAMPLE / "C1" / "XrayT001.out")
    y = reconstruct_xas_spectrum(c1, ref[:, 0], CARBON_XAS_SPEC)
    np.testing.assert_allclose(y, ref[:, 1], atol=0.01, rtol=0.0)


def test_aligned_spectrum_uses_delta_ks_not_ks_lumo(tmp_path: Path) -> None:
    site = tmp_path / "C1"
    site.mkdir()
    e_ha = 10.5
    osc = 0.01
    (site / "C1.xas").write_text(
        f"XAS      1\n         {e_ha:.10f}    {osc:.16f}\n",
        encoding="utf-8",
    )
    (site / "C1xas.inp").write_text(
        "RANGE 280 320\nPOINTS 201\nWIDTH 0.5 12 288 320\n",
        encoding="utf-8",
    )
    gnd_ha = -10.0
    target_shift = 2.0
    exc_ha = gnd_ha + (e_ha * STOBE_XAS_HA_TO_EV + target_shift) / HA_TO_EV
    (site / "C1gnd.out").write_text(
        _FINAL_BLOCK.format(energy=gnd_ha), encoding="utf-8"
    )
    (site / "C1exc.out").write_text(
        _FINAL_BLOCK.format(energy=exc_ha), encoding="utf-8"
    )
    (site / "C1tp.out").write_text(_FINAL_BLOCK.format(energy=-9.5), encoding="utf-8")
    _energy, spectra, metrics, sticks, _frame = collect_site_tp_xas(tmp_path)
    e_c = float(metrics["E_c_ev"].iloc[0])
    np.testing.assert_allclose(e_c, target_shift, atol=1e-6)
    aligned = spectra["abs_aligned"].to_numpy(dtype=np.float64)
    tp = spectra["abs_tp"].to_numpy(dtype=np.float64)
    energy = spectra["energy_ev"].to_numpy(dtype=np.float64)
    np.testing.assert_allclose(
        aligned[20:-20],
        shift_xas_spectrum(energy, tp, e_c)[20:-20],
        atol=1e-8,
        rtol=1e-5,
    )
    e1 = e_ha * STOBE_XAS_HA_TO_EV
    np.testing.assert_allclose(energy[int(np.argmax(aligned))], e1 + e_c, atol=0.25)
    np.testing.assert_allclose(
        float(sticks["energy_aligned_ev"].iloc[0]),
        e1 + e_c,
        atol=1e-9,
    )


def test_parse_missing_xas(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        parse_stobe_xas_sticks(tmp_path / "missing.xas")


def test_intensity_scale_constant() -> None:
    assert XAS_INTENSITY_SCALE == 1000.0
    assert STOBE_XAS_HA_TO_EV == 27.2116


def test_parse_stobe_xas_dipole_sticks_and_cartesian_os(tmp_path: Path) -> None:
    path = tmp_path / "C1.xas"
    energy_ha = 10.5418380964
    mux, muy, muz = 0.00239960212495, -0.00119737773330, 0.0298400912436
    osc = (2.0 / 3.0) * energy_ha * (mux**2 + muy**2 + muz**2)
    path.write_text(
        "XAS      1\n"
        f"         {energy_ha:.10f}    {osc:.16f}    0.0"
        f"  {mux:.16e} {muy:.16e} {muz:.16e}\n",
        encoding="utf-8",
    )
    dipole = parse_stobe_xas_dipole_sticks(path)
    assert dipole.shape == (1, 5)
    cart = dipole_cartesian_oscillator_strengths(dipole[:, 0], dipole[:, 2:5])
    np.testing.assert_allclose(cart.sum(axis=1), dipole[:, 1], rtol=1e-12)


def test_padded_shift_keeps_high_energy_tail() -> None:
    energy = np.linspace(280.0, 320.0, 401)
    intensity = np.linspace(1.0, 8.0, energy.size)
    shift = -1.6
    pad = padded_xas_energy_axis(energy, shift, extra_ev=2.0)
    y_pad = np.interp(pad, energy, intensity, left=intensity[0], right=intensity[-1])
    aligned = shift_and_resample_xas_spectrum(pad, y_pad, energy, shift)
    clipped = shift_xas_spectrum(energy, intensity, shift)
    assert float(clipped[-1]) == 0.0
    expected = float(np.interp(energy[-1] - shift, energy, intensity))
    np.testing.assert_allclose(aligned[-1], expected, rtol=1e-6)
    assert float(aligned[-1]) > 7.0


@pytest.mark.skipif(not _EXAMPLE.is_dir(), reason="example ZnPc run not present")
def test_znpc_aligned_spectrum_has_no_range_cliff() -> None:
    _energy, spectra, _metrics, sticks, _frame = collect_site_tp_xas(_EXAMPLE)
    c1 = spectra.loc[spectra["site"] == "C1"]
    y_al = c1["abs_aligned"].to_numpy(dtype=np.float64)
    y_tp = c1["abs_tp"].to_numpy(dtype=np.float64)
    y_xx = c1["abs_xx"].to_numpy(dtype=np.float64)
    y_yy = c1["abs_yy"].to_numpy(dtype=np.float64)
    y_zz = c1["abs_zz"].to_numpy(dtype=np.float64)
    np.testing.assert_allclose(y_xx + y_yy + y_zz, y_tp, atol=1e-8, rtol=1e-8)
    assert float(y_al[-1]) > 0.05 * float(np.max(y_tp))
    os_sum = (
        sticks.loc[sticks["site"] == "C1", ["os_xx", "os_yy", "os_zz"]]
        .sum(axis=1)
        .to_numpy(dtype=np.float64)
    )
    os_tot = sticks.loc[sticks["site"] == "C1", "oscillator_strength"].to_numpy(
        dtype=np.float64
    )
    np.testing.assert_allclose(os_sum, os_tot, atol=1e-10, rtol=1e-8)
    tensor_long, tensor_mean = aligned_dipole_tensor_tables(spectra, _metrics)
    c1_t = tensor_long.loc[tensor_long["site"] == "C1"]
    np.testing.assert_allclose(
        c1_t["energy_ev"].to_numpy(dtype=np.float64),
        c1["energy_ev"].to_numpy(dtype=np.float64)
        + float(_metrics.loc[_metrics["site"] == "C1", "E_c_ev"].iloc[0]),
        atol=1e-9,
    )
    assert tensor_mean.shape[0] == 2000
    np.testing.assert_allclose(
        tensor_mean["I_iso"].to_numpy(dtype=np.float64),
        (tensor_mean["I_xx"] + tensor_mean["I_yy"] + tensor_mean["I_zz"]).to_numpy(
            dtype=np.float64
        )
        / 3.0,
        atol=1e-10,
        rtol=1e-8,
    )


@pytest.mark.skipif(not _EXAMPLE.is_dir(), reason="example ZnPc run not present")
def test_write_xas_reconstruction_report_summary_png(tmp_path: Path) -> None:
    xyz = _EXAMPLE / "geometry.xyz"
    result = write_xas_reconstruction_report(
        _EXAMPLE,
        tmp_path,
        xyz_path=xyz if xyz.is_file() else None,
    )
    assert result is not None
    _long, _metrics, _sticks, tensor_long, tensor_mean, summary = result
    assert not (tmp_path / "xas_tp_reconstruction_validation.png").exists()
    assert tensor_long.is_file()
    assert tensor_mean.is_file()
    assert summary.is_file()
    assert summary.stat().st_size > 10_000
    long_df = pd.read_csv(tensor_long)
    mean_df = pd.read_csv(tensor_mean)
    assert list(long_df.columns) == [
        "energy_ev",
        "I_xx",
        "I_yy",
        "I_zz",
        "I_iso",
        "site",
        "E_c_ev",
    ]
    assert list(mean_df.columns) == ["energy_ev", "I_xx", "I_yy", "I_zz", "I_iso"]
    np.testing.assert_allclose(
        long_df["I_iso"].to_numpy(dtype=np.float64),
        (long_df["I_xx"] + long_df["I_yy"] + long_df["I_zz"]).to_numpy(dtype=np.float64)
        / 3.0,
        atol=1e-10,
        rtol=1e-8,
    )
