"""Persistent setup session for interactive calculation initialization."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from rdkit import Chem

from dftlearn.setup.alignment import (
    apply_rotation_to_rows,
    euler_matrix,
    identity_matrix,
)
from dftlearn.setup.basis_sets import stobe_basis_for_element, supported_basis_elements
from dftlearn.setup.cif_io import load_cif_structure
from dftlearn.setup.distinct_atoms import (
    distinct_atom_groups,
    distinct_groups_for_element,
    ligand_carbon_indices,
    mol_from_xyz_rows,
    select_principal_molecule_sites,
)
from dftlearn.setup.labeling import (
    SiteGroup,
    restrict_site_groups_to_atom_indices,
    site_groups_from_distinct,
)
from dftlearn.setup.types import DistinctAtomGroup

if TYPE_CHECKING:
    from dftlearn.setup.types import CifStructureMeta

SESSION_VERSION = 1
SESSION_FILENAME = "setup_session.json"


@dataclass
class SetupSession:
    """Editable initialization state for interactive StoBe run setup."""

    version: int
    mname: str
    edge_element: str
    cif_path: str
    block_name: str
    cif_labels: list[str]
    base_rows: list[list[float | str]]
    rotation_matrix: list[list[float]]
    alignment_matrix: list[list[float]]
    euler_deg: list[float]
    alignment_axis: list[float]
    alignment_note: str
    view: str
    site_groups: list[SiteGroup]
    meta: dict[str, Any] = field(default_factory=dict)

    @property
    def rotated_rows(self) -> list[tuple[str, float, float, float]]:
        """Return base rows after applying the stored rotation matrix."""
        matrix = np.asarray(self.rotation_matrix, dtype=np.float64)
        raw = [
            (str(r[0]), float(r[1]), float(r[2]), float(r[3])) for r in self.base_rows
        ]
        return apply_rotation_to_rows(raw, matrix)

    def basis_for_atom(self, atom_index: int) -> dict[str, str | None]:
        """Return StoBe basis strings for the element at ``atom_index``."""
        from dftlearn.io.xyz_structure import element_symbol_from_xyz_label

        label = str(self.base_rows[atom_index][0])
        sym = element_symbol_from_xyz_label(label)
        preset = stobe_basis_for_element(sym)
        basis = {
            "element": sym,
            "aux_ground": preset.aux_ground,
            "orbital_ground": preset.orbital_ground,
            "aux_excited": preset.aux_excited,
            "orbital_excited": preset.orbital_excited,
            "mcp": preset.mcp,
        }
        overrides = self.meta.get("basis_overrides", {}).get(str(atom_index), {})
        for key, value in overrides.items():
            if key in basis and value:
                basis[key] = str(value)
        return basis

    def set_basis_override(
        self,
        atom_index: int,
        field: str,
        value: str,
    ) -> None:
        """Store a per-atom basis override in ``meta['basis_overrides']``."""
        allowed = {
            "aux_ground",
            "orbital_ground",
            "aux_excited",
            "orbital_excited",
            "mcp",
        }
        if field not in allowed:
            msg = f"Unsupported basis field {field!r}"
            raise ValueError(msg)
        overrides = dict(self.meta.get("basis_overrides", {}))
        atom_key = str(atom_index)
        atom_overrides = dict(overrides.get(atom_key, {}))
        atom_overrides[field] = value.strip()
        overrides[atom_key] = atom_overrides
        self.meta["basis_overrides"] = overrides

    def clear_basis_override(self, atom_index: int, field: str | None = None) -> None:
        """Remove one or all basis overrides for ``atom_index``."""
        overrides = dict(self.meta.get("basis_overrides", {}))
        atom_key = str(atom_index)
        if atom_key not in overrides:
            return
        if field is None:
            overrides.pop(atom_key, None)
        else:
            atom_overrides = dict(overrides[atom_key])
            atom_overrides.pop(field, None)
            if atom_overrides:
                overrides[atom_key] = atom_overrides
            else:
                overrides.pop(atom_key, None)
        self.meta["basis_overrides"] = overrides

    def group_for_atom(self, atom_index: int) -> SiteGroup | None:
        """Return the enabled site group containing ``atom_index``, if any."""
        for group in self.site_groups:
            if group.enabled and atom_index in group.atom_indices:
                return group
        return None


def create_setup_session(
    cif_path: Path,
    *,
    mname: str,
    edge_element: str = "C",
    block_name: str | None = None,
    keep_all_sites: bool = False,
    one_ligand: bool = False,
) -> SetupSession:
    """Build a new setup session from a CIF file without writing run artifacts.

    By default keeps the covalently bonded principal molecule (metal-containing
    fragment when present) and drops residual solvent or synthesis leftovers.
    Pass ``keep_all_sites=True`` to retain every CIF atom site. Pass
    ``one_ligand=True`` to enable core-edge sites on a single metal-bound
    organic wing (the remaining ligands stay in the geometry).
    """
    cif_path = Path(cif_path)
    meta = load_cif_structure(cif_path, block_name=block_name)
    n_cif_sites = len(meta.sites)
    if not keep_all_sites:
        principal = select_principal_molecule_sites(meta.sites)
        meta = replace(meta, sites=principal)
    unsupported = [
        sym
        for sym in {s.element for s in meta.sites}
        if sym not in supported_basis_elements()
    ]
    if unsupported:
        msg = f"No StoBe basis presets for elements: {', '.join(sorted(unsupported))}"
        raise ValueError(msg)
    rows, cif_labels, groups = distinct_groups_for_element(
        meta.sites,
        element=edge_element,
    )
    site_groups = site_groups_from_distinct(groups)
    base_rows: list[list[float | str]] = [[label, x, y, z] for label, x, y, z in rows]
    mol = mol_from_xyz_rows(rows, sanitize=False)
    meta_dict = _meta_to_dict(meta)
    meta_dict["n_cif_sites"] = n_cif_sites
    meta_dict["n_molecule_sites"] = len(meta.sites)
    meta_dict["kept_all_sites"] = keep_all_sites
    if one_ligand:
        ligand = ligand_carbon_indices(list(rows), ligand_index=0)
        site_groups = restrict_site_groups_to_atom_indices(
            site_groups,
            ligand,
            edge_element,
        )
        meta_dict["one_ligand"] = True
        meta_dict["ligand_carbon_indices"] = list(ligand)
    try:
        meta_dict["mol_block"] = Chem.MolToMolBlock(mol)
    except (RuntimeError, ValueError):
        meta_dict["mol_block"] = None
    return SetupSession(
        version=SESSION_VERSION,
        mname=mname,
        edge_element=edge_element.strip().capitalize(),
        cif_path=str(cif_path.resolve()),
        block_name=meta.block_name,
        cif_labels=cif_labels,
        base_rows=base_rows,
        rotation_matrix=identity_matrix().tolist(),
        alignment_matrix=identity_matrix().tolist(),
        euler_deg=[0.0, 0.0, 0.0],
        alignment_axis=[0.0, 0.0, 1.0],
        alignment_note="identity",
        view="xy",
        site_groups=site_groups,
        meta=meta_dict,
    )


def create_setup_session_from_rows(
    rows: list[tuple[str, float, float, float]],
    *,
    mname: str,
    edge_element: str = "C",
    atom_labels: list[str] | None = None,
    source_meta: dict[str, Any] | None = None,
    cif_path: str = "",
    block_name: str = "",
    mol: Chem.Mol | None = None,
    one_ligand: bool = False,
) -> SetupSession:
    """Build a setup session from explicit XYZ rows (PubChem or relaxed structure).

    Parameters
    ----------
    rows
        ``(label, x, y, z)`` coordinates.
    mname
        Molecule title for StoBe inputs.
    edge_element
        Core-edge element for site grouping.
    atom_labels
        Optional external labels parallel to rows; defaults to XYZ labels.
    source_meta
        Extra metadata persisted under ``meta`` (PubChem CID, SMILES, etc.).
    cif_path
        Optional CIF path when the structure was imported from a crystal file.
    block_name
        Optional CIF block name.

    Returns
    -------
    SetupSession
        Session ready for the browser labeler.

    Raises
    ------
    ValueError
        If basis presets are missing or grouping fails.
    """
    labels = atom_labels or [label for label, *_ in rows]
    if len(labels) != len(rows):
        msg = "atom_labels length must match rows"
        raise ValueError(msg)

    unsupported = {element_symbol_from_row_label(label) for label, *_ in rows}
    missing = [sym for sym in unsupported if sym not in supported_basis_elements()]
    if missing:
        msg = f"No StoBe basis presets for elements: {', '.join(sorted(missing))}"
        raise ValueError(msg)

    grouping_mol = (
        mol
        if mol is not None
        else mol_from_xyz_rows(
            list(rows),
            sanitize=False,
        )
    )
    groups = distinct_atom_groups(
        grouping_mol,
        element=edge_element,
        rows=list(rows),
        cif_labels=labels,
    )
    site_groups = site_groups_from_distinct(groups)
    meta = dict(source_meta or {})
    if one_ligand:
        ligand = ligand_carbon_indices(list(rows), ligand_index=0)
        site_groups = restrict_site_groups_to_atom_indices(
            site_groups,
            ligand,
            edge_element,
        )
        meta["one_ligand"] = True
        meta["ligand_carbon_indices"] = list(ligand)
    base_rows = [[label, x, y, z] for label, x, y, z in rows]
    meta.setdefault("source_kind", "pubchem" if not cif_path else "cif")
    if mol is not None:
        try:
            meta["mol_block"] = Chem.MolToMolBlock(mol)
        except (RuntimeError, ValueError):
            meta["mol_block"] = None
    return SetupSession(
        version=SESSION_VERSION,
        mname=mname,
        edge_element=edge_element.strip().capitalize(),
        cif_path=cif_path,
        block_name=block_name,
        cif_labels=labels,
        base_rows=base_rows,
        rotation_matrix=identity_matrix().tolist(),
        alignment_matrix=identity_matrix().tolist(),
        euler_deg=[0.0, 0.0, 0.0],
        alignment_axis=[0.0, 0.0, 1.0],
        alignment_note="identity",
        view="xy",
        site_groups=site_groups,
        meta=meta,
    )


def element_symbol_from_row_label(label: str) -> str:
    """Return the element symbol encoded in one XYZ row label."""
    from dftlearn.io.xyz_structure import element_symbol_from_xyz_label as _sym

    return _sym(label)


def _meta_to_dict(meta: CifStructureMeta) -> dict[str, Any]:
    return {
        "space_group": meta.space_group,
        "symmetry_ops": list(meta.symmetry_ops),
        "cell": {
            "a": meta.cell_a,
            "b": meta.cell_b,
            "c": meta.cell_c,
            "alpha": meta.cell_alpha,
            "beta": meta.cell_beta,
            "gamma": meta.cell_gamma,
        },
    }


def session_path(run_directory: Path) -> Path:
    """Return the canonical setup session path inside a run directory."""
    return Path(run_directory) / SESSION_FILENAME


def save_setup_session(run_directory: Path, session: SetupSession) -> Path:
    """Serialize ``session`` to ``setup_session.json``."""
    path = session_path(run_directory)
    payload = _session_to_dict(session)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return path


def load_setup_session(run_directory: Path) -> SetupSession:
    """Load ``setup_session.json`` from a run directory."""
    path = session_path(run_directory)
    if not path.is_file():
        msg = f"Setup session not found: {path}"
        raise FileNotFoundError(msg)
    payload = json.loads(path.read_text(encoding="utf-8"))
    return _session_from_dict(payload)


def update_session_rotation(
    session: SetupSession,
    rx_deg: float,
    ry_deg: float,
    rz_deg: float,
) -> SetupSession:
    """Replace manual Euler offsets and recompute the total rotation matrix."""
    manual = euler_matrix(rx_deg, ry_deg, rz_deg)
    base = np.asarray(session.alignment_matrix, dtype=np.float64)
    total = manual @ base
    session.euler_deg = [rx_deg, ry_deg, rz_deg]
    session.rotation_matrix = total.tolist()
    return session


def set_session_alignment_matrix(
    session: SetupSession,
    matrix: np.ndarray,
    *,
    note: str,
) -> SetupSession:
    """Set the base alignment matrix and reset manual Euler offsets to zero."""
    session.alignment_matrix = np.asarray(matrix, dtype=np.float64).tolist()
    session.alignment_note = note
    session.euler_deg = [0.0, 0.0, 0.0]
    session.rotation_matrix = session.alignment_matrix.copy()
    return session


def _site_group_to_dict(group: SiteGroup) -> dict[str, Any]:
    return asdict(group)


def _site_group_from_dict(payload: dict[str, Any]) -> SiteGroup:
    return SiteGroup(**payload)


def _session_to_dict(session: SetupSession) -> dict[str, Any]:
    data = asdict(session)
    data["site_groups"] = [_site_group_to_dict(g) for g in session.site_groups]
    return data


def _session_from_dict(payload: dict[str, Any]) -> SetupSession:
    groups = [_site_group_from_dict(g) for g in payload.get("site_groups", [])]
    return SetupSession(
        version=int(payload.get("version", SESSION_VERSION)),
        mname=str(payload["mname"]),
        edge_element=str(payload["edge_element"]),
        cif_path=str(payload["cif_path"]),
        block_name=str(payload["block_name"]),
        cif_labels=[str(x) for x in payload["cif_labels"]],
        base_rows=payload["base_rows"],
        rotation_matrix=payload["rotation_matrix"],
        alignment_matrix=payload.get(
            "alignment_matrix",
            payload["rotation_matrix"],
        ),
        euler_deg=[float(x) for x in payload.get("euler_deg", [0.0, 0.0, 0.0])],
        alignment_axis=[float(x) for x in payload.get("alignment_axis", [0, 0, 1])],
        alignment_note=str(payload.get("alignment_note", "")),
        view=str(payload.get("view", "xy")),
        site_groups=groups,
        meta=dict(payload.get("meta", {})),
    )


def distinct_groups_from_session(session: SetupSession) -> list[DistinctAtomGroup]:
    """Convert enabled session site groups to distinct groups on rotated rows."""
    rows = session.rotated_rows
    from dftlearn.setup.labeling import active_generating_groups

    active = active_generating_groups(session.site_groups)
    updated: list[DistinctAtomGroup] = []
    for group in active:
        rep = group.representative_index
        updated.append(
            DistinctAtomGroup(
                rank=group.rank,
                element=group.element,
                atom_indices=group.atom_indices,
                representative_index=rep,
                cif_label=group.cif_label,
                xyz_label=rows[rep][0],
            )
        )
    return updated
