"""PubChem PUG REST client for compound search and structure fetch."""

from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass

_PUG_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_TIMEOUT_S = 30.0


@dataclass(frozen=True)
class PubChemCompound:
    """One PubChem compound hit with identifiers and SMILES."""

    cid: int
    title: str
    molecular_formula: str
    smiles: str
    iupac_name: str


def _fetch_json(url: str) -> dict:
    request = urllib.request.Request(
        url,
        headers={"Accept": "application/json"},
    )
    try:
        with urllib.request.urlopen(request, timeout=_TIMEOUT_S) as response:
            return json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return {}
        raise


def search_pubchem(query: str, *, limit: int = 20) -> list[PubChemCompound]:
    """Search PubChem by name or molecular formula and return compound records.

    Tries a literal name lookup first, then the autocomplete suggestion list,
    then a fast-formula lookup when ``query`` looks like a formula.

    Parameters
    ----------
    query
        Common name, IUPAC fragment, or molecular formula such as ``C27H18AlN3O3``.
    limit
        Maximum number of compounds to return.

    Returns
    -------
    list[PubChemCompound]
        Compounds with CID, title, formula, and SMILES.

    Raises
    ------
    ValueError
        If ``query`` is empty.
    RuntimeError
        If PubChem requests fail unexpectedly.
    """
    text = query.strip()
    if not text:
        msg = "PubChem query must not be empty"
        raise ValueError(msg)

    cids: list[int] = []
    name_url = f"{_PUG_BASE}/compound/name/{urllib.parse.quote(text)}/cids/JSON"
    payload = _fetch_json(name_url)
    cid_list = payload.get("IdentifierList", {}).get("CID", [])
    if cid_list:
        cids.append(int(cid_list[0]))

    if not cids:
        auto_url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/"
            f"{urllib.parse.quote(text)}/json"
        )
        auto = _fetch_json(auto_url)
        suggestions = auto.get("autocomplete", [])[:limit]
        for suggestion in suggestions:
            sug_url = (
                f"{_PUG_BASE}/compound/name/{urllib.parse.quote(suggestion)}/cids/JSON"
            )
            sug_payload = _fetch_json(sug_url)
            for cid in sug_payload.get("IdentifierList", {}).get("CID", []):
                cids.append(int(cid))
                break

    if not cids and _looks_like_formula(text):
        form_url = (
            f"{_PUG_BASE}/compound/fastformula/{urllib.parse.quote(text)}/cids/JSON"
        )
        form_payload = _fetch_json(form_url)
        for cid in form_payload.get("IdentifierList", {}).get("CID", []):
            cids.append(int(cid))

    seen: set[int] = set()
    ordered: list[int] = []
    for cid in cids:
        if cid not in seen:
            seen.add(cid)
            ordered.append(cid)
        if len(ordered) >= limit:
            break

    if not ordered:
        return []

    cid_path = ",".join(str(cid) for cid in ordered)
    prop_url = (
        f"{_PUG_BASE}/compound/cid/{cid_path}/property/"
        "Title,IUPACName,MolecularFormula,CanonicalSMILES,IsomericSMILES/JSON"
    )
    prop_payload = _fetch_json(prop_url)
    compounds: list[PubChemCompound] = []
    for row in prop_payload.get("PropertyTable", {}).get("Properties", []):
        smiles = row.get("SMILES") or row.get("ConnectivitySMILES") or ""
        if not smiles:
            continue
        compounds.append(
            PubChemCompound(
                cid=int(row["CID"]),
                title=str(
                    row.get("Title") or row.get("IUPACName") or f"CID {row['CID']}"
                ),
                molecular_formula=str(row.get("MolecularFormula") or ""),
                smiles=str(smiles),
                iupac_name=str(row.get("IUPACName") or ""),
            )
        )
    return compounds


def fetch_pubchem_smiles(cid: int) -> str:
    """Return the canonical SMILES string for one PubChem CID."""
    url = (
        f"{_PUG_BASE}/compound/cid/{cid}/property/"
        "CanonicalSMILES,IsomericSMILES/JSON"
    )
    payload = _fetch_json(url)
    rows = payload.get("PropertyTable", {}).get("Properties", [])
    if not rows:
        msg = f"No SMILES found for PubChem CID {cid}"
        raise ValueError(msg)
    row = rows[0]
    smiles = row.get("SMILES") or row.get("ConnectivitySMILES")
    if not smiles:
        msg = f"Empty SMILES for PubChem CID {cid}"
        raise ValueError(msg)
    return str(smiles)


def _looks_like_formula(text: str) -> bool:
    if not text:
        return False
    allowed = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789()")
    return all(ch in allowed for ch in text) and any(ch.isdigit() for ch in text)
