"""Tests for StoBe basis catalog validation."""

from __future__ import annotations

from pathlib import Path

import pytest

from dftlearn.setup.basis_catalog import (
    bundled_basis_names,
    missing_basis_names,
    validate_basis_preset,
    validate_molconfig_basis,
)
from dftlearn.setup.basis_sets import (
    stobe_basis_for_element,
    validate_all_basis_presets,
    validate_element_basis,
)
from dftlearn.setup.types import StoBeBasisPreset


def test_bundled_catalog_contains_al_and_o_presets() -> None:
    names = bundled_basis_names()
    assert "A-ALUMINUM (5,4;5,4)" in names
    assert "O-OXYGEN (631/31/1)" in names
    assert "O-OXYGEN (32/3)" not in names


def test_validate_all_presets_pass() -> None:
    assert validate_all_basis_presets() == {}


def test_validate_element_basis_o() -> None:
    assert validate_element_basis("O") == []


def test_missing_basis_names_detects_typo() -> None:
    missing = missing_basis_names(["O-OXYGEN (32/3)", "A-CARBON (5,2;5,2)"])
    assert missing == ["O-OXYGEN (32/3)"]


def test_validate_molconfig_basis_from_source() -> None:
    source = (
        'elem4_Abasis = "A-OXYGEN (4,4;4,4)"\n'
        'elem4_Obasis = "O-OXYGEN (32/3)"\n'
    )
    assert validate_molconfig_basis(source) == ["O-OXYGEN (32/3)"]


def test_validate_basis_preset_reports_missing_orbital() -> None:
    preset = StoBeBasisPreset(
        atomic_number=8,
        z_eff_ground=8,
        z_eff_excited=8,
        aux_ground="A-OXYGEN (4,4;4,4)",
        aux_excited="A-OXYGEN (4,4;4,4)",
        orbital_ground="O-OXYGEN (32/3)",
        orbital_excited="O-OXYGEN (32/3)",
        mcp=None,
    )
    assert validate_basis_preset(preset) == ["O-OXYGEN (32/3)"]


def test_als3_molconfig_matches_baslib() -> None:
    molconfig = Path("docs/stobe/example-run/als3-001/molConfig.py")
    if not molconfig.is_file():
        pytest.skip("example run molConfig not present")
    missing = validate_molconfig_basis(molconfig.read_text(encoding="utf-8"))
    assert missing == []


def test_stobe_basis_for_element_o_uses_baslib_name() -> None:
    preset = stobe_basis_for_element("O")
    assert preset.orbital_ground == "O-OXYGEN (631/31/1)"
    assert validate_element_basis("O") == []
