"""Tests for StoBe X-ray table parsing."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from dftlearn.io.xray_out import parse_xray_out_table, site_spectra_to_long_frame


def test_parse_xray_out_table_fortran_d(tmp_path: Path) -> None:
    p = tmp_path / "XrayT001.out"
    p.write_text(
        "        280.00000000      0.10000000D+01\n"
        "        280.02001001      0.20000000d-01\n",
        encoding="utf-8",
    )
    arr = parse_xray_out_table(p)
    assert arr.shape == (2, 2)
    np.testing.assert_allclose(arr[0], (280.0, 1.0))
    np.testing.assert_allclose(arr[1], (280.02001001, 0.02))


def test_site_spectra_to_long_frame_orders_sites() -> None:
    energy = np.array([1.0, 2.0])
    spectra = {"C10": np.array([1.0, 2.0]), "C2": np.array([3.0, 4.0])}
    df = site_spectra_to_long_frame(energy, spectra)
    assert list(df.columns) == ["energy_ev", "abs", "site"]
    sites = df["site"].unique().tolist()
    assert sites == ["C2", "C10"]


def test_parse_missing_file(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        parse_xray_out_table(tmp_path / "missing.out")
