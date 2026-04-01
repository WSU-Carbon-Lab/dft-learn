"""Tests for XYZ parsing and site indexing."""

from __future__ import annotations

from pathlib import Path

import pytest

from dftlearn.io.xyz_structure import (
    site_label_to_atom_index,
    site_label_to_atom_index_from_rows,
    xyz_rows_from_file,
)


def test_xyz_headerless_stobe_style(tmp_path: Path) -> None:
    p = tmp_path / "g.xyz"
    p.write_text(
        "C01 0.0 0.0 0.0\n"
        "C02 1.0 0.0 0.0\n"
        "H01 0.5 1.0 0.0\n",
        encoding="utf-8",
    )
    rows = xyz_rows_from_file(p)
    assert len(rows) == 3
    assert rows[0][0] == "C01"


def test_site_label_from_rows_matches_path(tmp_path: Path) -> None:
    p = tmp_path / "g.xyz"
    p.write_text(
        "C01 0 0 0\n"
        "C02 1 0 0\n"
        "N01 2 0 0\n",
        encoding="utf-8",
    )
    rows = xyz_rows_from_file(p)
    assert site_label_to_atom_index_from_rows("C1", rows) == 0
    assert site_label_to_atom_index_from_rows("C2", rows) == 1
    assert site_label_to_atom_index_from_rows("N1", rows) == 2
    assert site_label_to_atom_index("C2", p) == site_label_to_atom_index_from_rows(
        "C2", rows
    )


def test_xyz_standard_two_line_header(tmp_path: Path) -> None:
    p = tmp_path / "g.xyz"
    p.write_text(
        "3\ncomment line\n"
        "C 0.0 0.0 0.0\n"
        "C 1.4 0.0 0.0\n"
        "H 0.7 1.0 0.0\n",
        encoding="utf-8",
    )
    rows = xyz_rows_from_file(p)
    assert len(rows) == 3
    assert site_label_to_atom_index("C1", p) == 0
    assert site_label_to_atom_index("C2", p) == 1


def test_site_label_to_atom_index_requires_match(tmp_path: Path) -> None:
    p = tmp_path / "g.xyz"
    p.write_text("N01 0 0 0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="not found"):
        site_label_to_atom_index("C1", p)
