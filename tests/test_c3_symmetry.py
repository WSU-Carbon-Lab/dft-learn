"""Unit tests for Al-N/O C3 frame construction and dipole tensor folding."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from dftlearn.xas.c3_symmetry import (
    build_al_n_o_c3_frame,
    fold_dipole_c3_tensor,
    fold_dipoles_c3,
    oscillator_strengths_from_tensors,
    rotation_matrix_z,
    write_c3_frame_json,
)
from dftlearn.xas.spectrum import (
    collect_site_tp_xas,
    dipole_cartesian_oscillator_strengths,
)

_SQRT3_2 = float(np.sqrt(3.0) * 0.5)


def _ideal_alq3_rows() -> list[tuple[str, float, float, float]]:
    """Al at origin with three N-O bites 120 degrees apart in xy."""
    bites = (
        (1.0, 0.0),
        (-0.5, _SQRT3_2),
        (-0.5, -_SQRT3_2),
    )
    rows: list[tuple[str, float, float, float]] = [
        ("Al01", 0.0, 0.0, 0.0),
    ]
    for i, (cx, cy) in enumerate(bites, start=1):
        rows.append((f"N{i:02d}", cx + 0.2, cy, 0.15))
        rows.append((f"O{i:02d}", cx - 0.2, cy, -0.15))
    return rows


def test_ideal_triangle_frame_is_identity() -> None:
    frame = build_al_n_o_c3_frame(_ideal_alq3_rows())
    np.testing.assert_allclose(frame.rotation, np.eye(3), atol=1e-12)
    assert frame.al_index == 0
    assert len(frame.ligand_pairs) == 3


def test_fold_in_plane_dipole_equalizes_xx_yy() -> None:
    frame = build_al_n_o_c3_frame(_ideal_alq3_rows())
    mu = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    tensor = fold_dipole_c3_tensor(mu, frame.rotation)
    np.testing.assert_allclose(tensor[0, 0], tensor[1, 1], atol=1e-12)
    np.testing.assert_allclose(tensor[2, 2], 0.0, atol=1e-12)
    np.testing.assert_allclose(np.trace(tensor), 1.0, atol=1e-12)


def test_fold_out_of_plane_dipole_leaves_zz() -> None:
    frame = build_al_n_o_c3_frame(_ideal_alq3_rows())
    mu = np.array([0.0, 0.0, 2.0], dtype=np.float64)
    tensor = fold_dipole_c3_tensor(mu, frame.rotation)
    np.testing.assert_allclose(tensor[0, 0], 0.0, atol=1e-12)
    np.testing.assert_allclose(tensor[1, 1], 0.0, atol=1e-12)
    np.testing.assert_allclose(tensor[2, 2], 4.0, atol=1e-12)
    np.testing.assert_allclose(np.trace(tensor), 4.0, atol=1e-12)


def test_fold_normalization_preserves_mu_squared() -> None:
    frame = build_al_n_o_c3_frame(_ideal_alq3_rows())
    mu = np.array([0.3, -0.7, 0.4], dtype=np.float64)
    tensor = fold_dipole_c3_tensor(mu, frame.rotation)
    np.testing.assert_allclose(np.trace(tensor), float(np.dot(mu, mu)), atol=1e-12)


def test_rotation_matrix_z_120() -> None:
    rz = rotation_matrix_z(120.0)
    v = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    out = rz @ v
    np.testing.assert_allclose(out, np.array([-0.5, _SQRT3_2, 0.0]), atol=1e-12)


def test_oscillator_strengths_from_folded_tensor() -> None:
    energy = np.array([10.0], dtype=np.float64)
    tensors = fold_dipoles_c3(
        np.array([[1.0, 0.0, 0.0]], dtype=np.float64),
        np.eye(3, dtype=np.float64),
    )
    cart = oscillator_strengths_from_tensors(energy, tensors)
    np.testing.assert_allclose(cart[0, 0], cart[0, 1], atol=1e-12)
    np.testing.assert_allclose(cart[0, 2], 0.0, atol=1e-12)
    expected_total = (2.0 / 3.0) * 10.0 * 1.0
    np.testing.assert_allclose(cart.sum(), expected_total, atol=1e-12)


def test_write_c3_frame_json(tmp_path: Path) -> None:
    frame = build_al_n_o_c3_frame(_ideal_alq3_rows())
    path = write_c3_frame_json(frame, tmp_path / "c3_frame.json")
    assert path.is_file()
    text = path.read_text(encoding="utf-8")
    assert '"al_index": 0' in text
    assert '"rotation"' in text


def _write_fake_site_with_dipole(
    root: Path,
    *,
    mux: float,
    muy: float,
    muz: float,
) -> None:
    site = root / "C1"
    site.mkdir()
    energy_ha = 10.5
    osc = (2.0 / 3.0) * energy_ha * (mux**2 + muy**2 + muz**2)
    (site / "C1.xas").write_text(
        "XAS      1\n"
        f"         {energy_ha:.10f}    {osc:.16f}    0.0"
        f"  {mux:.16e} {muy:.16e} {muz:.16e}\n",
        encoding="utf-8",
    )
    (site / "C1xas.inp").write_text(
        "RANGE 280 320\nPOINTS 51\nWIDTH 0.5 12 288 320\n",
        encoding="utf-8",
    )


def _write_ideal_xyz(path: Path) -> None:
    rows = _ideal_alq3_rows()
    lines = [f"{len(rows)}", "ideal"]
    for label, x, y, z in rows:
        lines.append(f"{label} {x:.8f} {y:.8f} {z:.8f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_collect_site_tp_xas_c3_equalizes_in_plane(tmp_path: Path) -> None:
    _write_fake_site_with_dipole(tmp_path, mux=1.0, muy=0.0, muz=0.0)
    xyz = tmp_path / "geometry.xyz"
    _write_ideal_xyz(xyz)
    _e, _s, _m, sticks, frame = collect_site_tp_xas(
        tmp_path,
        xyz_path=xyz,
        c3_symmetrize=True,
    )
    assert frame is not None
    row = sticks.iloc[0]
    np.testing.assert_allclose(float(row["os_xx"]), float(row["os_yy"]), atol=1e-12)
    np.testing.assert_allclose(float(row["os_zz"]), 0.0, atol=1e-12)


def test_collect_site_tp_xas_flag_off_matches_raw_cartesian(tmp_path: Path) -> None:
    mux, muy, muz = 0.2, -0.1, 0.5
    _write_fake_site_with_dipole(tmp_path, mux=mux, muy=muy, muz=muz)
    xyz = tmp_path / "geometry.xyz"
    _write_ideal_xyz(xyz)
    _e, _s, _m, sticks, frame = collect_site_tp_xas(
        tmp_path,
        xyz_path=xyz,
        c3_symmetrize=False,
    )
    assert frame is None
    energy_ha = 10.5
    expected = dipole_cartesian_oscillator_strengths(
        np.array([energy_ha]),
        np.array([[mux, muy, muz]]),
    )[0]
    np.testing.assert_allclose(
        sticks.loc[0, ["os_xx", "os_yy", "os_zz"]].to_numpy(dtype=np.float64),
        expected,
        atol=1e-12,
    )


def test_c3_symmetrize_requires_xyz(tmp_path: Path) -> None:
    _write_fake_site_with_dipole(tmp_path, mux=1.0, muy=0.0, muz=0.0)
    with pytest.raises(ValueError, match="xyz_path"):
        collect_site_tp_xas(tmp_path, c3_symmetrize=True)
