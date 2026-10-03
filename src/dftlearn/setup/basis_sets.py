"""StoBe basis-set presets keyed by element symbol."""

from __future__ import annotations

from dftlearn.setup.basis_catalog import validate_basis_preset
from dftlearn.setup.types import StoBeBasisPreset

_STOBE_BASIS: dict[str, StoBeBasisPreset] = {
    "C": StoBeBasisPreset(
        atomic_number=6,
        z_eff_ground=6,
        z_eff_excited=4,
        aux_ground="A-CARBON (5,2;5,2)",
        aux_excited="A-CARBON(+4) (3,3;3,3)",
        orbital_ground="O-CARBON iii_iglo",
        orbital_excited="O-CARBON(+4) (311/211/1)",
        mcp="P-CARBON(+4) (3,1:8,0)",
    ),
    "H": StoBeBasisPreset(
        atomic_number=1,
        z_eff_ground=1,
        z_eff_excited=1,
        aux_ground="A-HYDROGEN (4,2;4,2)",
        aux_excited="A-HYDROGEN (4,2;4,2)",
        orbital_ground="O-HYDROGEN (311/1) misc",
        orbital_excited="O-HYDROGEN (311/1) misc",
    ),
    "N": StoBeBasisPreset(
        atomic_number=7,
        z_eff_ground=7,
        z_eff_excited=7,
        aux_ground="A-NITROGEN (4,3;4,3)",
        aux_excited="A-NITROGEN (4,3;4,3)",
        orbital_ground="O-NITROGEN (33/3)",
        orbital_excited="O-NITROGEN (33/3)",
    ),
    "O": StoBeBasisPreset(
        atomic_number=8,
        z_eff_ground=8,
        z_eff_excited=8,
        aux_ground="A-OXYGEN (4,4;4,4)",
        aux_excited="A-OXYGEN (4,4;4,4)",
        orbital_ground="O-OXYGEN (631/31/1)",
        orbital_excited="O-OXYGEN (631/31/1)",
    ),
    "Al": StoBeBasisPreset(
        atomic_number=13,
        z_eff_ground=13,
        z_eff_excited=13,
        aux_ground="A-ALUMINUM (5,4;5,4)",
        aux_excited="A-ALUMINUM (5,4;5,4)",
        orbital_ground="O-ALUMINUM (73111/6111/1)",
        orbital_excited="O-ALUMINUM (73111/6111/1)",
    ),
    "Zn": StoBeBasisPreset(
        atomic_number=30,
        z_eff_ground=30,
        z_eff_excited=30,
        aux_ground="A-ZINC (5,5;5,5)",
        aux_excited="A-ZINC (5,5;5,5)",
        orbital_ground="O-ZINC (63321/531/311)",
        orbital_excited="O-ZINC (63321/531/311)",
    ),
}


def stobe_basis_for_element(symbol: str) -> StoBeBasisPreset:
    """Return the StoBe basis preset for ``symbol``.

    Parameters
    ----------
    symbol
        Element symbol such as ``C`` or ``Al``.

    Returns
    -------
    StoBeBasisPreset
        Auxiliary, orbital, and MCP strings for StoBe input generation.

    Raises
    ------
    KeyError
        If no preset is registered for ``symbol``.
    """
    key = symbol.strip().capitalize()
    if key == "Zn":
        return _STOBE_BASIS["Zn"]
    if key not in _STOBE_BASIS:
        msg = f"No StoBe basis preset for element {symbol!r}"
        raise KeyError(msg)
    return _STOBE_BASIS[key]


def supported_basis_elements() -> tuple[str, ...]:
    """Return element symbols with registered StoBe basis presets."""
    return tuple(sorted(_STOBE_BASIS))


def validate_element_basis(symbol: str) -> list[str]:
    """Return basis names in the preset for ``symbol`` missing from the catalog.

    Parameters
    ----------
    symbol
        Element symbol such as ``O`` or ``Al``.

    Returns
    -------
    list[str]
        Missing StoBe basis library names; empty when the preset is valid.
    """
    return validate_basis_preset(stobe_basis_for_element(symbol))


def validate_all_basis_presets() -> dict[str, list[str]]:
    """Validate every registered element preset against the bundled baslib catalog.

    Returns
    -------
    dict[str, list[str]]
        Maps element symbols to missing basis names; omitted when valid.
    """
    invalid: dict[str, list[str]] = {}
    for symbol in _STOBE_BASIS:
        missing = validate_element_basis(symbol)
        if missing:
            invalid[symbol] = missing
    return invalid
