"""Tests for setup session, labeling, and alignment helpers."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from dftlearn.setup.alignment import euler_matrix, rotation_align_vector
from dftlearn.setup.labeling import merge_site_groups, site_groups_from_distinct
from dftlearn.setup.session import (
    create_setup_session,
    load_setup_session,
    save_setup_session,
    set_session_alignment_matrix,
    update_session_rotation,
)
from dftlearn.setup.types import DistinctAtomGroup
from dftlearn.setup.workflow import finalize_setup_session


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


def _sample_groups() -> list[DistinctAtomGroup]:
    return [
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


def test_merge_site_groups_combines_classes() -> None:
    groups = site_groups_from_distinct(_sample_groups())
    merged = merge_site_groups(groups, [1, 2])
    enabled = [g for g in merged if g.enabled]
    assert len(enabled) == 1
    assert enabled[0].atom_indices == (0, 1, 2)


def test_session_rotation_composes_with_alignment(minimal_cif: Path) -> None:
    session = create_setup_session(minimal_cif, mname="test")
    base = rotation_align_vector(np.array([1.0, 0.0, 0.0]), np.array([0.0, 0.0, 1.0]))
    set_session_alignment_matrix(session, base, note="test")
    update_session_rotation(session, 10.0, 0.0, 0.0)
    expected = euler_matrix(10.0, 0.0, 0.0) @ base
    actual = np.asarray(session.rotation_matrix)
    assert np.allclose(actual, expected, atol=1e-10)


def test_setup_session_roundtrip(tmp_path: Path, minimal_cif: Path) -> None:
    session = create_setup_session(minimal_cif, mname="test")
    save_setup_session(tmp_path, session)
    loaded = load_setup_session(tmp_path)
    assert loaded.mname == "test"
    assert len(loaded.site_groups) == len(session.site_groups)


def test_finalize_setup_session_writes_files(tmp_path: Path, minimal_cif: Path) -> None:
    session = create_setup_session(minimal_cif, mname="test")
    artifacts = finalize_setup_session(tmp_path, session, write_preview=False)
    assert artifacts.geometry_xyz.is_file()
    assert (tmp_path / "setup_session.json").is_file()
    summary = json.loads(artifacts.summary_json.read_text(encoding="utf-8"))
    assert "alignment" in summary
