"""Tests for StoBe run discovery and composed load_stobe_run."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from dftlearn.io.stobe_run import (
    discover_site_directories,
    discover_stobe_site_paths,
    find_xyz,
    load_stobe_run,
)

_FINAL_TAIL = """
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

_ORBITAL = """
 ORBITAL ENERGIES (ALL VIRTUALS INCLUDED)

         Spin alpha                              Spin beta
         Occup.    Energy(eV)    Sym  (pos.)     Occup.    Energy(eV)    Sym  (pos.)
    1    1.0000    -10.0000    1A   (   1)     1.0000    -10.0000    1A   (   1)
    2    0.0000     -2.5000    2A   (   2)     0.0000     -2.6000    2A   (   2)
"""

_TP_IP = """
 Orbital energy core hole =    -10.70906 H   (  -291.41063 eV)
 Ionization potential     =    291.41063 eV
"""

_SCF = """
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

_XRAY = (
    "        280.00000000      0.10000000D+01\n"
    "        281.00000000      0.20000000D+01\n"
)


def _final_block(energy_h: float) -> str:
    return _FINAL_TAIL.replace("-2438.6024592313", f"{energy_h:.10f}")


def _write_toy_run(root: Path, *, sites: tuple[str, ...] = ("C1", "C2")) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    (root / "mol.xyz").write_text(
        "2\ncomment\n"
        "C 0.0 0.0 0.0\n"
        "C 1.4 0.0 0.0\n",
        encoding="utf-8",
    )
    for i, site in enumerate(sites):
        site_dir = root / site
        site_dir.mkdir()
        gnd = -10.0 - 0.1 * i
        exc = -9.0 - 0.1 * i
        tp = -9.5 - 0.1 * i
        (site_dir / f"{site}gnd.out").write_text(
            _SCF + _final_block(gnd),
            encoding="utf-8",
        )
        (site_dir / f"{site}exc.out").write_text(
            _SCF + _final_block(exc),
            encoding="utf-8",
        )
        (site_dir / f"{site}tp.out").write_text(
            _SCF + _final_block(tp) + _ORBITAL + _TP_IP,
            encoding="utf-8",
        )
        (site_dir / "XrayT001.out").write_text(_XRAY, encoding="utf-8")
    return root


def test_discover_site_directories_skips_reserved(tmp_path: Path) -> None:
    run = tmp_path / "run"
    (run / "C1").mkdir(parents=True)
    (run / "C2").mkdir()
    (run / "NEXAFS").mkdir()
    (run / "packaged_output").mkdir()
    dirs = discover_site_directories(run)
    assert [p.name for p in dirs] == ["C1", "C2"]


def test_find_xyz_single_and_explicit(tmp_path: Path) -> None:
    run = _write_toy_run(tmp_path / "benzene")
    assert find_xyz(run).name == "mol.xyz"
    other = tmp_path / "alt.xyz"
    other.write_text("1\n\nH 0 0 0\n", encoding="utf-8")
    assert find_xyz(run, xyz=other) == other.resolve()


def test_find_xyz_ambiguous(tmp_path: Path) -> None:
    run = tmp_path / "run"
    run.mkdir()
    (run / "a.xyz").write_text("1\n\nH 0 0 0\n", encoding="utf-8")
    (run / "b.xyz").write_text("1\n\nH 0 0 0\n", encoding="utf-8")
    with pytest.raises(ValueError, match=r"Multiple \.xyz"):
        find_xyz(run)


def test_discover_stobe_site_paths(tmp_path: Path) -> None:
    run = _write_toy_run(tmp_path / "run")
    paths = discover_stobe_site_paths(run)
    assert [p.site for p in paths] == ["C1", "C2"]
    assert paths[0].gnd_out is not None
    assert paths[0].xray_out is not None
    assert paths[0].xray_out.name == "XrayT001.out"


def test_load_stobe_run_bundle(tmp_path: Path) -> None:
    run = _write_toy_run(tmp_path / "run")
    bundle = load_stobe_run(run, require_xyz=True)
    assert bundle.xyz_path is not None
    assert len(bundle.xyz_rows) == 2
    assert set(bundle.final_energies_long["site"]) == {"C1", "C2"}
    assert not bundle.delta_ks_sites.empty
    assert not bundle.scf_convergence_long.empty
    assert bundle.xray_energy_ev is not None
    assert set(bundle.xray_spectra) == {"C1", "C2"}
    assert set(bundle.xray_long["site"]) == {"C1", "C2"}
    np.testing.assert_allclose(bundle.xray_energy_ev[0], 280.0)
