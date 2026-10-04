"""Tests for Henke .nff parsing and IP-weighted Gaussian step edges."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from dftlearn.xas.step_edge import (
    compound_mu,
    gaussian_step,
    gaussian_step_edge,
    parse_henke_nff,
)

_HENKE = Path(__file__).resolve().parent / "data" / "henke"


def test_parse_henke_nff_carbon() -> None:
    df = parse_henke_nff(_HENKE / "c.nff")
    assert list(df.columns) == ["energy_ev", "f1", "f2"]
    assert len(df) == 5
    assert np.isclose(df["energy_ev"].iloc[0], 280.0)
    assert np.isclose(df["f2"].iloc[2], 0.80)


def test_parse_henke_nff_missing(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        parse_henke_nff(tmp_path / "missing.nff")


def test_compound_mu_positive_and_finite() -> None:
    tables = {
        "C": parse_henke_nff(_HENKE / "c.nff"),
        "H": parse_henke_nff(_HENKE / "h.nff"),
    }
    energy = np.linspace(282.0, 305.0, 50)
    mu = compound_mu(energy, {"C": 6, "H": 6}, tables)
    assert mu.shape == energy.shape
    assert np.all(np.isfinite(mu))
    assert np.all(mu > 0.0)


def test_compound_mu_missing_element() -> None:
    tables = {"C": parse_henke_nff(_HENKE / "c.nff")}
    with pytest.raises(ValueError, match="Henke table missing"):
        compound_mu(np.array([290.0]), {"C": 1, "O": 1}, tables)


def test_gaussian_step_unit_height() -> None:
    energy = np.linspace(280.0, 300.0, 401)
    step = gaussian_step(energy, center_ev=290.0, width_ev=0.5)
    assert np.isclose(step[0], 0.0, atol=1e-3)
    assert np.isclose(step[-1], 1.0, atol=1e-3)
    assert np.isclose(step[np.argmin(np.abs(energy - 290.0))], 0.5, atol=0.02)


def test_gaussian_step_edge_monotonic_and_jump() -> None:
    energy = np.linspace(280.0, 310.0, 601)
    ips = np.array([288.0, 292.0])
    edge = gaussian_step_edge(
        energy,
        ips,
        width_ev=0.5,
        total_jump=1.0,
        weight_by_inverse_ip=True,
    )
    assert np.all(np.diff(edge) >= -1e-12)
    assert np.isclose(edge[-1], 1.0, atol=1e-3)
    assert edge[0] < 0.05


def test_gaussian_step_edge_equal_weights() -> None:
    energy = np.linspace(280.0, 310.0, 401)
    edge = gaussian_step_edge(
        energy,
        np.array([290.0, 295.0]),
        width_ev=0.4,
        total_jump=2.0,
        weight_by_inverse_ip=False,
    )
    assert np.isclose(edge[-1], 2.0, atol=1e-3)
