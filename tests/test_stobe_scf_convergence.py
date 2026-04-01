"""Tests for StoBe SCF convergence table parsing and metrics."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest

if TYPE_CHECKING:
    from pathlib import Path

from dftlearn.io.stobe_scf_convergence import (
    collect_scf_convergence_long,
    parse_stobe_scf_convergence_table,
    scf_convergence_auc_metrics,
)

MINIMAL_SCF_TAIL = """
 ------------------------------------------------------------------------------
 SCF ITERATION STARTS NOW
 ------------------------------------------------------------------------------


ITER      TOTAL ENERGY      DECREASE    AVER-DENSTY   MAX-DENSITY    DIIS     CPU
   1    -1290.57467690     0.00000000    10.000000     0.000000   0.00000    65.18
   2    -1957.77395875   667.19928185     0.324316    33.491079   1.55494     6.67
  DIIS turned on...
   3    -2068.52134157   110.74738282     1.259808   126.224036   0.79796     6.56
 ------------------------------------------------------------------------------
 SCF CONVERGED AFTER   3 ITERATIONS
 ------------------------------------------------------------------------------
"""


def test_parse_stobe_scf_convergence_table_diis_flag(tmp_path: Path) -> None:
    """``diis_active`` flips after the ``DIIS turned on`` line."""
    path = tmp_path / "Xexc.out"
    path.write_text("preamble\n" + MINIMAL_SCF_TAIL, encoding="utf-8")
    df = parse_stobe_scf_convergence_table(path)
    assert len(df) == 3
    assert not bool(df.loc[df["iteration"] == 1, "diis_active"].iloc[0])
    assert not bool(df.loc[df["iteration"] == 2, "diis_active"].iloc[0])
    assert bool(df.loc[df["iteration"] == 3, "diis_active"].iloc[0])
    te3 = df.loc[df["iteration"] == 3, "total_energy_h"].iloc[0]
    assert np.isclose(te3, -2068.52134157)


def test_collect_and_metrics(tmp_path: Path) -> None:
    """Collect one site folder and aggregate AUC metrics."""
    c1 = tmp_path / "C1"
    c1.mkdir()
    (c1 / "C1gnd.out").write_text(MINIMAL_SCF_TAIL, encoding="utf-8")
    long_df = collect_scf_convergence_long(tmp_path, calc_types=("gnd",))
    assert set(long_df["site"]) == {"C1"}
    assert set(long_df["calc_type"]) == {"gnd"}
    m = scf_convergence_auc_metrics(long_df)
    assert len(m) == 1
    assert m["diis_start_iter"].iloc[0] == 3.0
    assert m["energy_auc"].iloc[0] > 0
    assert m["density_auc"].iloc[0] > 0


def test_collect_raises_when_no_data(tmp_path: Path) -> None:
    """Empty run root raises a clear error."""
    with pytest.raises(ValueError, match="No SCF convergence"):
        collect_scf_convergence_long(tmp_path)
