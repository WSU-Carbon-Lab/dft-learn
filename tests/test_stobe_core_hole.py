"""Tests for MCP valence counting and GND-driven core-hole assignment."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from dftlearn.io.stobe_orbital_table import parse_stobe_orbital_energies_table
from dftlearn.setup.core_hole import (
    assign_core_holes_from_gnd,
    format_fsym_alfa_occupation,
    locate_k_edge_core_orbital,
    patch_fsym_alfa_occupation,
)
from dftlearn.setup.molconfig import UNASSIGNED_ALFA_OCC, valence_electron_count
from dftlearn.setup.types import InitGeometryPlan

if TYPE_CHECKING:
    from pathlib import Path


def _plan(counts: dict[str, int], edge: str = "C") -> InitGeometryPlan:
    order = tuple(counts)
    return InitGeometryPlan(
        rows=tuple(("C01", 0.0, 0.0, 0.0) for _ in range(sum(counts.values()) or 1)),
        edge_element=edge,
        generating_site_indices=(0,),
        distinct_groups=(),
        element_counts=counts,
        element_group_order=order,
    )


def _orb_row(level: int, energy_ev: float, occ: float = 1.0) -> str:
    return (
        f"{level:5d}    {occ:.4f}    {energy_ev:8.4f}    1A   (   {level})     "
        f"{occ:.4f}    {energy_ev:8.4f}    1A   (   {level})"
    )


def _orbital_table_text(rows: list[tuple[int, float, float]]) -> str:
    body = "\n".join(_orb_row(lv, en, occ) for lv, en, occ in rows)
    return (
        " ORBITAL ENERGIES (ALL VIRTUALS INCLUDED)\n"
        "\n"
        "         Spin alpha                              Spin beta\n"
        "         Occup.    Energy(eV)    Sym  (pos.)     Occup.    Energy(eV)"
        "    Sym  (pos.)\n"
        f"{body}\n"
        " MULLIKEN POPULATION ANALYSIS\n"
    )


def test_alq3_valence_electrons_subtract_mcp_carbon_1s() -> None:
    """27 C with MCP spectators yield 186 valence electrons (93/93)."""
    plan = _plan({"C": 27, "H": 18, "N": 3, "O": 3, "Al": 1})
    assert valence_electron_count(plan) == 186


def test_three_carbon_mcp_plus_two_hydrogen() -> None:
    """Three MCP carbons plus two hydrogens are 16 valence electrons."""
    plan = _plan({"C": 3, "H": 2})
    assert valence_electron_count(plan) == 16


def test_format_fsym_alfa_occupation() -> None:
    """EXC zeros one MO; TP half-fills the same index (nspec, orbital, occ)."""
    assert format_fsym_alfa_occupation(8, 0.0) == "0 1 8 0.0"
    assert format_fsym_alfa_occupation(8, 0.5) == "0 1 8 0.5"


def test_locate_carbon_1s_among_deeper_cores(tmp_path: Path) -> None:
    """C 1s is the occupied KS level near -285 eV, not Al/O/N 1s."""
    path = tmp_path / "C1gnd.out"
    path.write_text(
        _orbital_table_text(
            [
                (1, -1559.0, 1.0),
                (2, -543.0, 1.0),
                (3, -410.0, 1.0),
                (8, -278.4, 1.0),
                (9, -12.0, 1.0),
                (10, 1.5, 0.0),
            ]
        ),
        encoding="utf-8",
    )
    df = parse_stobe_orbital_energies_table(path)
    report = locate_k_edge_core_orbital(df, 284.8)
    assert report.core_level == 8
    assert report.core_energy_ev == pytest.approx(-278.4)


def test_locate_rejects_aluminum_1s_as_carbon_edge(tmp_path: Path) -> None:
    """A table without a C 1s fails the K-edge energy window."""
    path = tmp_path / "C1gnd.out"
    path.write_text(
        _orbital_table_text(
            [
                (1, -1559.0, 1.0),
                (2, -12.0, 1.0),
                (3, 1.5, 0.0),
            ]
        ),
        encoding="utf-8",
    )
    df = parse_stobe_orbital_energies_table(path)
    with pytest.raises(ValueError, match="expected the absorber 1s"):
        locate_k_edge_core_orbital(df, 284.8)


def test_patch_fsym_alfa_occupation_replaces_placeholder() -> None:
    """Four-token ALFA override is replaced; total ALFA count is left intact."""
    text = (
        "FSYM scfocc excited\n"
        "ALFA 94\n"
        "BETA 93\n"
        "SYM 1\n"
        f"ALFA {UNASSIGNED_ALFA_OCC}\n"
        "BETA 0 0\n"
        "END\n"
    )
    patched = patch_fsym_alfa_occupation(text, "0 1 8 0.0")
    assert "ALFA 0 1 8 0.0\n" in patched
    assert "ALFA 94\n" in patched


def test_assign_core_holes_from_gnd_patches_exc_and_tp(tmp_path: Path) -> None:
    """GND ionization energy writes per-site EXC 0.0 and TP 0.5 occupations."""
    site = tmp_path / "C1"
    site.mkdir()
    gnd = tmp_path / "GND"
    gnd.mkdir()
    (gnd / "C1gnd.out").write_text(
        _orbital_table_text(
            [
                (1, -1559.0, 1.0),
                (8, -278.4, 1.0),
                (9, -12.0, 1.0),
                (10, 1.5, 0.0),
            ]
        ),
        encoding="utf-8",
    )
    block = (
        "FSYM scfocc excited\n"
        "ALFA 94\n"
        "BETA 93\n"
        "SYM 1\n"
        f"ALFA {UNASSIGNED_ALFA_OCC}\n"
        "BETA 0 0\n"
        "END\n"
    )
    (site / "C1exc.run").write_text(block, encoding="utf-8")
    (site / "C1tp.run").write_text(block, encoding="utf-8")
    (tmp_path / "molConfig.py").write_text(
        f'alfaOcc = "{UNASSIGNED_ALFA_OCC}"\nalfaOccTP = "{UNASSIGNED_ALFA_OCC}"\n',
        encoding="utf-8",
    )

    rows = assign_core_holes_from_gnd(tmp_path)
    assert len(rows) == 1
    assert rows[0].core_level == 8
    assert rows[0].ionization_energy_ev == pytest.approx(278.4)
    assert "ALFA 0 1 8 0.0" in (site / "C1exc.run").read_text(encoding="utf-8")
    assert "ALFA 0 1 8 0.5" in (site / "C1tp.run").read_text(encoding="utf-8")
    mol = (tmp_path / "molConfig.py").read_text(encoding="utf-8")
    assert 'alfaOcc = "0 1 8 0.0"' in mol
    assert 'alfaOccTP = "0 1 8 0.5"' in mol
    assert (tmp_path / "core_hole_assignments.json").is_file()
