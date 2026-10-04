"""Tests for library package_stobe_run orchestration."""

from __future__ import annotations

from pathlib import Path

from dftlearn.pipeline import package_stobe_run

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


def _write_toy_run(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    (root / "mol.xyz").write_text(
        "2\ncomment\n"
        "C 0.0 0.0 0.0\n"
        "C 1.4 0.0 0.0\n",
        encoding="utf-8",
    )
    for i, site in enumerate(("C1", "C2")):
        site_dir = root / site
        site_dir.mkdir()
        gnd = -10.0 - 0.1 * i
        exc = -9.0 - 0.1 * i
        tp = -9.5 - 0.1 * i
        (site_dir / f"{site}gnd.out").write_text(
            _SCF + _final_block(gnd) + _ORBITAL,
            encoding="utf-8",
        )
        (site_dir / f"{site}exc.out").write_text(
            _SCF + _final_block(exc) + _ORBITAL,
            encoding="utf-8",
        )
        (site_dir / f"{site}tp.out").write_text(
            _SCF + _final_block(tp) + _ORBITAL + _TP_IP,
            encoding="utf-8",
        )
        (site_dir / "XrayT001.out").write_text(_XRAY, encoding="utf-8")
    return root


def test_package_stobe_run_writes_site_scf_final(tmp_path: Path) -> None:
    run = _write_toy_run(tmp_path / "run")
    result = package_stobe_run(run, dpi=72)
    packaged = result.packaged_dir
    assert packaged == (run / "packaged_output").resolve()

    assert (packaged / "xray_spectra_long.csv").is_file()
    assert (packaged / "xas_site_summary.png").is_file()
    assert (packaged / "scf_convergence_long.csv").is_file()
    assert (packaged / "scf_convergence_metrics.csv").is_file()
    assert (packaged / "scf_diagnostics.png").is_file()
    assert (packaged / "stobe_final_energies.csv").is_file()
    assert (packaged / "stobe_orbital_energy_summary.png").is_file()

    written_names = {p.name for p in result.written}
    assert "xray_spectra_long.csv" in written_names
    assert "xas_site_summary.png" in written_names
    assert "scf_diagnostics.png" in written_names
    assert "stobe_final_energies.csv" in written_names
    assert "stobe_orbital_energy_summary.png" in written_names

    assert any("xas" in s.lower() or "stick" in s.lower() for s in result.skipped)
