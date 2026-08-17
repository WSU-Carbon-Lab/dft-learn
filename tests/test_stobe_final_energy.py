"""Tests for StoBe FINAL ENERGY block parsing and aggregation."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest

from dftlearn.io.stobe_final_energy import (
    HA_TO_EV,
    collect_delta_ks_site_table,
    collect_final_energies_long,
    enrich_final_energies_delta_vs_gnd,
    final_energy_site_summary,
    parse_stobe_final_energy_tables,
    parse_stobe_tp_lumo_alpha_ev,
)

if TYPE_CHECKING:
    from pathlib import Path

MINIMAL_FINAL_TAIL = """
 ------------------------------------------------------------------------------
 SCF CONVERGED AFTER   3 ITERATIONS
 ------------------------------------------------------------------------------


 FINAL ENERGY / CHARGE / GEOMETRY RESULTS :

 Total energy   (H) =   -2438.6024592313  (incl. numerical value for EXC)
 Nuc-nuc energy (H) =    3202.0663146615
 El-nuc energy  (H) =  -12209.0179634253
 Kinetic energy (H) =    2423.5316608697
 Coulomb energy (H) =    4337.4622015007
 Ex-cor energy  (H) =    -192.6446728379

 <Rho/r12/Rhof>-<Rhof/r12/Rhof>/2 (H) =    4337.7936250498
 <Rho/r12/Rhof>/2                 (H) =    4337.5277743277
 Total exchange energy            (H) =    -183.6063174584
 Total correlation energy         (H) =      -9.0383553795

 Decomposition of exchange / correlation :
"""

MINIMAL_ORBITAL_TABLE = """
 ORBITAL ENERGIES (ALL VIRTUALS INCLUDED)

         Spin alpha                              Spin beta
         Occup.    Energy(eV)    Sym  (pos.)     Occup.    Energy(eV)    Sym  (pos.)
    1    1.0000    -10.0000    1A   (   1)     1.0000    -10.0000    1A   (   1)
    2    0.0000     -2.5000    2A   (   2)     0.0000     -2.6000    2A   (   2)
"""


def test_parse_stobe_final_energy_tables(tmp_path: Path) -> None:
    """Parse one block with all standard fields."""
    path = tmp_path / "C1gnd.out"
    path.write_text("head\n" + MINIMAL_FINAL_TAIL, encoding="utf-8")
    df = parse_stobe_final_energy_tables(path)
    assert len(df) == 1
    assert df["block_index"].iloc[0] == 0
    assert np.isclose(df["total_energy_h"].iloc[0], -2438.6024592313)
    assert np.isclose(df["nuc_nuc_energy_h"].iloc[0], 3202.0663146615)
    assert np.isclose(df["total_correlation_energy_h"].iloc[0], -9.0383553795)


def test_collect_and_delta_vs_gnd(tmp_path: Path) -> None:
    """Ground reference drives ``delta_vs_gnd_ev`` for exc and tp."""
    c1 = tmp_path / "C1"
    c1.mkdir()
    gnd_e = -100.0
    exc_e = -90.0
    tp_e = -95.0

    def block(e: float) -> str:
        return MINIMAL_FINAL_TAIL.replace("-2438.6024592313", f"{e:.10f}")

    (c1 / "C1gnd.out").write_text(block(gnd_e), encoding="utf-8")
    (c1 / "C1exc.out").write_text(block(exc_e), encoding="utf-8")
    (c1 / "C1tp.out").write_text(block(tp_e), encoding="utf-8")

    long_df = collect_final_energies_long(tmp_path)
    assert set(long_df["site"]) == {"C1"}
    gnd_row = long_df[(long_df["calc_type"] == "gnd") & (long_df["site"] == "C1")]
    exc_row = long_df[(long_df["calc_type"] == "exc") & (long_df["site"] == "C1")]
    assert np.isnan(gnd_row["delta_vs_gnd_ev"].iloc[0])
    assert np.isclose(exc_row["delta_vs_gnd_h"].iloc[0], exc_e - gnd_e)
    summ = final_energy_site_summary(long_df)
    assert len(summ) == 3


def test_parse_raises_when_missing_block(tmp_path: Path) -> None:
    """No FINAL ENERGY header yields ``ValueError``."""
    path = tmp_path / "empty.out"
    path.write_text("ITER only\n", encoding="utf-8")
    with pytest.raises(ValueError, match="No FINAL ENERGY"):
        parse_stobe_final_energy_tables(path)


def test_enrich_idempotent_gnd_nan(tmp_path: Path) -> None:
    """``enrich_final_energies_delta_vs_gnd`` leaves gnd deltas undefined."""
    import pandas as pd

    df = pd.DataFrame(
        {
            "site": ["C1", "C1"],
            "calc_type": ["gnd", "exc"],
            "block_index": [0, 0],
            "total_energy_h": [-10.0, -9.0],
        }
    )
    out = enrich_final_energies_delta_vs_gnd(df)
    assert np.isnan(out.loc[out["calc_type"] == "gnd", "delta_vs_gnd_h"].iloc[0])
    assert np.isclose(out.loc[out["calc_type"] == "exc", "delta_vs_gnd_h"].iloc[0], 1.0)


def test_parse_stobe_tp_lumo_alpha_ev(tmp_path: Path) -> None:
    """First alpha-unoccupied row yields LUMO in eV."""
    path = tmp_path / "C1tp.out"
    path.write_text(MINIMAL_FINAL_TAIL + MINIMAL_ORBITAL_TABLE, encoding="utf-8")
    assert np.isclose(parse_stobe_tp_lumo_alpha_ev(path), -2.5)


def test_collect_delta_ks_site_table(tmp_path: Path) -> None:
    """Wide table matches :math:`E^c = E^e - E^g - E^l` in eV."""
    c1 = tmp_path / "C1"
    c1.mkdir()
    gnd_ha = -10.0
    exc_ha = -9.0
    tp_ha = -9.5
    lumo_ev = -2.5

    def block(e: float) -> str:
        return MINIMAL_FINAL_TAIL.replace("-2438.6024592313", f"{e:.10f}")

    (c1 / "C1gnd.out").write_text(block(gnd_ha), encoding="utf-8")
    (c1 / "C1exc.out").write_text(block(exc_ha), encoding="utf-8")
    (c1 / "C1tp.out").write_text(
        block(tp_ha) + MINIMAL_ORBITAL_TABLE,
        encoding="utf-8",
    )
    df = collect_delta_ks_site_table(tmp_path)
    assert len(df) == 1
    eg = gnd_ha * HA_TO_EV
    ee = exc_ha * HA_TO_EV
    assert np.isclose(float(df["E_c_deltaKS_ev"].iloc[0]), ee - eg - lumo_ev)
