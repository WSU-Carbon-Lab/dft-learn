"""Tests for CIF initialization and distinct atom grouping."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from dftlearn.setup import (
    cif_blocks,
    initialize_from_cif,
    load_cif_structure,
    plan_stobe_geometry,
    stobe_basis_for_element,
)
from dftlearn.setup.basis_sets import validate_element_basis
from dftlearn.setup.distinct_atoms import (
    distinct_groups_for_element,
    select_principal_molecule_sites,
)

ALQ3_CIF = (
    Path(__file__).resolve().parents[1]
    / "docs"
    / "stobe"
    / "example-run"
    / "alq3"
    / "cg301291w_si_003.cif"
)


@pytest.fixture
def minimal_cif(tmp_path: Path) -> Path:
    text = """data_test_alq3
_cell_length_a 10.0
_cell_length_b 10.0
_cell_length_c 10.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0.0 0.0 0.0
C2 C 1.4 0.0 0.0
C3 C 2.8 0.0 0.0
H1 H 0.0 1.0 0.0
H2 H 0.0 -1.0 0.0
"""
    path = tmp_path / "minimal.cif"
    path.write_text(text, encoding="utf-8")
    return path


def test_load_cif_structure_minimal(minimal_cif: Path) -> None:
    meta = load_cif_structure(minimal_cif)
    assert meta.block_name == "test_alq3"
    assert len(meta.sites) == 5
    assert meta.sites[0].element == "C"


def test_distinct_carbons_minimal(minimal_cif: Path) -> None:
    meta = load_cif_structure(minimal_cif)
    _rows, _labels, groups = distinct_groups_for_element(meta.sites, element="C")
    assert len(groups) == 1
    assert len(groups[0].atom_indices) == 3


def test_plan_stobe_geometry_generators_first(minimal_cif: Path) -> None:
    meta = load_cif_structure(minimal_cif)
    rows, _labels, groups = distinct_groups_for_element(meta.sites, element="C")
    plan = plan_stobe_geometry(rows, edge_element="C", distinct_groups=groups)
    assert plan.rows[0][0].startswith("C")
    assert len(plan.generating_site_indices) == 1
    assert plan.generating_site_indices == (0,)


def test_initialize_from_cif_writes_artifacts(
    tmp_path: Path, minimal_cif: Path
) -> None:
    run_dir = tmp_path / "run"
    artifacts = initialize_from_cif(
        run_dir,
        minimal_cif,
        name="test-run",
        write_preview=True,
    )
    assert artifacts.geometry_xyz.is_file()
    assert artifacts.molconfig_py.is_file()
    assert artifacts.summary_json.is_file()
    summary = json.loads(artifacts.summary_json.read_text(encoding="utf-8"))
    assert summary["n_core_sites"] == 1
    mol_src = artifacts.molconfig_py.read_text(encoding="utf-8")
    assert "alpha = 3" in mol_src
    assert 'alfaOcc = "0 0 0 0.0"' in mol_src
    assert "basis_sets" in summary
    assert stobe_basis_for_element("C").mcp is not None


@pytest.mark.skipif(not ALQ3_CIF.is_file(), reason="AlQ3 CIF not in repo")
def test_alq3_principal_molecule_drops_acetic_acid() -> None:
    """CIF co-crystallized acetic acid is excluded from the bare AlQ3 molecule."""
    meta = load_cif_structure(ALQ3_CIF)
    assert len(meta.sites) == 60
    molecule = select_principal_molecule_sites(meta.sites)
    counts: dict[str, int] = {}
    for site in molecule:
        counts[site.element] = counts.get(site.element, 0) + 1
    assert len(molecule) == 52
    assert counts == {"C": 27, "H": 18, "N": 3, "O": 3, "Al": 1}
    labels = {site.label for site in molecule}
    assert "C31" not in labels
    assert "O31" not in labels
    assert "O32" not in labels


@pytest.mark.skipif(not ALQ3_CIF.is_file(), reason="AlQ3 CIF not in repo")
def test_alq3_one_ligand_keeps_nine_carbons() -> None:
    """One quinolate wing is nine generating carbons; the molecule stays 52 atoms."""
    from dftlearn.setup.session import create_setup_session

    session = create_setup_session(
        ALQ3_CIF,
        mname="alq3",
        one_ligand=True,
    )
    enabled = [g for g in session.site_groups if g.enabled]
    assert len(enabled) == 9
    assert {g.site_tag for g in enabled} == {f"C{i}" for i in range(1, 10)}
    assert session.meta.get("one_ligand") is True
    assert len(session.meta.get("ligand_carbon_indices", [])) == 9


@pytest.mark.skipif(not ALQ3_CIF.is_file(), reason="AlQ3 CIF not in repo")
def test_alq3_init_uses_bare_molecule(tmp_path: Path) -> None:
    run_dir = tmp_path / "alq3"
    artifacts = initialize_from_cif(run_dir, ALQ3_CIF, name="alq3-003")
    summary = json.loads(artifacts.summary_json.read_text(encoding="utf-8"))
    assert summary["n_core_sites"] >= 9
    xyz = artifacts.geometry_xyz.read_text(encoding="utf-8").splitlines()
    assert len(xyz) == 52
    mol_src = artifacts.molconfig_py.read_text(encoding="utf-8")
    assert "alpha = 93" in mol_src
    assert "beta = 93" in mol_src
    assert 'alfaOcc = "0 0 0 0.0"' in mol_src
    assert stobe_basis_for_element("Al").aux_ground == "A-ALUMINUM (5,4;5,4)"
    assert stobe_basis_for_element("O").orbital_ground == "O-OXYGEN (631/31/1)"
    assert validate_element_basis("O") == []


def test_cif_blocks_empty(tmp_path: Path) -> None:
    path = tmp_path / "empty.cif"
    path.write_text("data_empty\n", encoding="utf-8")
    assert cif_blocks(path) == []
