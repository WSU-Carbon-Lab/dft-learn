"""Tests for bond alignment, structure visualization, and basis overrides."""

from __future__ import annotations

import numpy as np

from dftlearn.setup.alignment import (
    align_bond_to_axis,
    align_frame_from_bond,
    bond_vector,
    rotate_positions,
)
from dftlearn.setup.distinct_atoms import mol_from_xyz_rows
from dftlearn.setup.labeling import site_groups_from_distinct
from dftlearn.setup.session import create_setup_session_from_rows
from dftlearn.setup.structure_viz import (
    bonds_from_mol,
    display_bonds,
    minimal_generating_indices,
    molblock_from_rows,
)
from dftlearn.setup.types import DistinctAtomGroup


def _linear_rows() -> list[tuple[str, float, float, float]]:
    return [
        ("C01", 0.0, 0.0, 0.0),
        ("C02", 1.4, 0.0, 0.0),
        ("H01", 0.0, 1.0, 0.0),
    ]


def test_bond_vector_points_along_x() -> None:
    rows = _linear_rows()
    vec = bond_vector(rows, 0, 1)
    assert np.allclose(vec, [1.4, 0.0, 0.0])


def test_align_bond_to_x_axis() -> None:
    rows = _linear_rows()
    matrix = align_bond_to_axis(rows, 0, 1, "x")
    aligned = rotate_positions(np.array([bond_vector(rows, 0, 1)]), matrix)[0]
    norm = aligned / np.linalg.norm(aligned)
    assert np.allclose(norm, [1.0, 0.0, 0.0], atol=1e-10)


def test_align_frame_from_bond_puts_plane_atom_in_xy() -> None:
    rows = _linear_rows()
    matrix = align_frame_from_bond(rows, 0, 1, 2)
    pos = np.array([[r[1], r[2], r[3]] for r in rows], dtype=np.float64)
    plane = rotate_positions(pos[2:3], matrix)[0]
    assert abs(plane[2]) < 1e-10


def test_mol_from_xyz_rows_tolerates_overbonded_geometry() -> None:
    """Greedy valence pruning prevents sanitize failures on crowded coordinates."""
    rows = [("C01", 0.0, 0.0, 0.0)]
    for idx in range(1, 8):
        angle = 2.0 * np.pi * idx / 7.0
        rows.append(("C02", 1.4 * np.cos(angle), 1.4 * np.sin(angle), 0.0))
    mol = mol_from_xyz_rows(rows, sanitize=True)
    assert mol.GetNumAtoms() == 8
    assert mol.GetAtomWithIdx(0).GetDegree() <= 4


def test_display_bonds_merges_rdkit_and_distance_inference() -> None:
    rows = _linear_rows()
    mol = mol_from_xyz_rows(rows, sanitize=False)
    merged = display_bonds(rows, mol)
    assert len(merged) >= len(bonds_from_mol(mol))


def test_molblock_preserves_bonds() -> None:
    rows = _linear_rows()
    mol = mol_from_xyz_rows(rows)
    block = molblock_from_rows(rows, mol)
    assert "M  END" in block
    rebuilt = mol_from_xyz_rows(rows)
    assert len(bonds_from_mol(rebuilt)) >= 2


def test_basis_override_roundtrip() -> None:
    rows = _linear_rows()
    session = create_setup_session_from_rows(
        rows,
        mname="test",
        edge_element="C",
        atom_labels=["C1", "C2", "H1"],
        mol=mol_from_xyz_rows(rows),
    )
    session.set_basis_override(0, "aux_ground", "CUSTOM-A")
    basis = session.basis_for_atom(0)
    assert basis["aux_ground"] == "CUSTOM-A"
    session.clear_basis_override(0)
    basis = session.basis_for_atom(0)
    assert basis["aux_ground"] != "CUSTOM-A"


def test_minimal_generating_indices() -> None:
    groups = site_groups_from_distinct(
        [
            DistinctAtomGroup(
                rank=1,
                element="C",
                atom_indices=(0, 1),
                representative_index=0,
                cif_label="C1",
                xyz_label="C01",
            ),
            DistinctAtomGroup(
                rank=2,
                element="C",
                atom_indices=(2,),
                representative_index=2,
                cif_label="C3",
                xyz_label="C03",
            ),
        ]
    )
    reps = minimal_generating_indices(groups)
    assert reps == (0, 2)
