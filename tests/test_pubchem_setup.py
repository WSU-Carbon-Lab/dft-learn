"""Tests for PubChem client and SMILES structure building."""

from __future__ import annotations

import pytest

from dftlearn.setup.pubchem_client import search_pubchem
from dftlearn.setup.session import create_setup_session_from_rows
from dftlearn.setup.structure_build import RelaxMethod, build_3d_from_smiles


@pytest.mark.integration
def test_search_pubchem_alq3_formula() -> None:
    hits = search_pubchem("C27H18AlN3O3", limit=5)
    assert hits
    assert any("Al" in h.molecular_formula for h in hits)


def test_build_3d_from_smiles_alq3_rejects_collapsed_ligands() -> None:
    smiles = (
        "C1=CC2=C(C(=C1)[O-])N=CC=C2."
        "C1=CC2=C(C(=C1)[O-])N=CC=C2."
        "C1=CC2=C(C(=C1)[O-])N=CC=C2.[Al+3]"
    )
    with pytest.raises(ValueError, match="overlapping atoms"):
        build_3d_from_smiles(smiles, relax=RelaxMethod.UFF, relax_steps=200)


def test_create_session_from_pubchem_rows() -> None:
    built = build_3d_from_smiles("c1ccccc1", relax=RelaxMethod.NONE)
    session = create_setup_session_from_rows(
        built.rows,
        mname="benzene-test",
        atom_labels=built.atom_labels,
        source_meta={"pubchem_cid": 241, "smiles": "c1ccccc1"},
        mol=built.mol,
    )
    assert session.meta["source_kind"] == "pubchem"
    assert session.site_groups
