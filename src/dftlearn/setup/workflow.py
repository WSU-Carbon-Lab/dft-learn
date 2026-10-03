"""Orchestration for ``dftrun init``: CIF to StoBe run directory."""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import TYPE_CHECKING

from dftlearn.setup.basis_catalog import validate_molconfig_basis
from dftlearn.setup.basis_sets import stobe_basis_for_element, validate_element_basis
from dftlearn.setup.cif_io import cif_blocks
from dftlearn.setup.geometry import plan_stobe_geometry, write_stobe_xyz
from dftlearn.setup.molconfig import render_molconfig_py, write_molconfig_py
from dftlearn.setup.preview import write_init_preview_figure
from dftlearn.setup.remote import (
    remote_run_directory,
    resolve_ssh_host,
    sync_run_directory,
)
from dftlearn.setup.session import (
    create_setup_session,
    distinct_groups_from_session,
    save_setup_session,
)
from dftlearn.setup.types import InitArtifacts

if TYPE_CHECKING:
    from dftlearn.setup.session import (
        SetupSession,
    )
    from dftlearn.setup.types import CifStructureMeta, InitGeometryPlan


def finalize_setup_session(
    run_directory: Path,
    session: SetupSession,
    *,
    write_preview: bool = True,
    remote_host: str | None = None,
    remote_root: str = "/home/hduva/projects/dft-runs",
) -> InitArtifacts:
    """Write StoBe artifacts from an edited setup session.

    Persists ``setup_session.json``, ``geometry.xyz``, ``molConfig.py``,
    ``init_summary.json``, and optional preview figure.

    Parameters
    ----------
    run_directory
        Target run directory.
    session
        Edited setup session state.
    write_preview
        When True, write ``init_preview.png``.
    remote_host
        Optional SSH alias to rsync the run directory after writing.
    remote_root
        Remote base directory for rsync.

    Returns
    -------
    InitArtifacts
        Paths to generated outputs.
    """
    run_directory = Path(run_directory)
    run_directory.mkdir(parents=True, exist_ok=True)
    save_setup_session(run_directory, session)

    rows = session.rotated_rows
    groups = distinct_groups_from_session(session)
    plan = plan_stobe_geometry(
        list(rows),
        edge_element=session.edge_element,
        distinct_groups=groups,
    )

    geometry_path = run_directory / "geometry.xyz"
    write_stobe_xyz(geometry_path, plan, comment=session.mname)

    molconfig_source = render_molconfig_py(
        fname="geometry",
        mname=session.mname,
        plan=plan,
    )
    missing_basis = validate_molconfig_basis(molconfig_source)
    if missing_basis:
        msg = (
            "molConfig references basis names missing from StoBe baslib: "
            + ", ".join(missing_basis)
        )
        raise ValueError(msg)
    for sym in _session_element_symbols(session):
        element_missing = validate_element_basis(sym)
        if element_missing:
            msg = (
                f"Basis preset for {sym} references names missing from baslib: "
                + ", ".join(element_missing)
            )
            raise ValueError(msg)
    molconfig_path = run_directory / "molConfig.py"
    write_molconfig_py(molconfig_path, molconfig_source)

    source_cif = run_directory / "source.cif"
    cif_src = Path(session.cif_path) if session.cif_path else None
    if cif_src is not None and cif_src.is_file():
        shutil.copy2(cif_src, source_cif)
    elif session.meta.get("smiles"):
        source_cif.write_text(
            f"# PubChem structure\n# SMILES: {session.meta.get('smiles')}\n",
            encoding="utf-8",
        )

    meta = _meta_from_session(session)
    preview_path: Path | None = None
    if write_preview:
        preview_path = run_directory / "init_preview.png"
        write_init_preview_figure(
            preview_path,
            plan=plan,
            meta=meta,
            title=session.mname,
        )

    summary_path = run_directory / "init_summary.json"
    summary_cif = cif_src if cif_src is not None else Path(session.cif_path or ".")
    summary = _build_summary_from_session(session, plan, summary_cif)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    remote_directory: Path | None = None
    resolved_remote: str | None = None
    if remote_host:
        resolved_remote = resolve_ssh_host(remote_host)
        remote_directory = remote_run_directory(
            resolved_remote,
            run_directory.name,
            remote_root,
        )
        sync_run_directory(
            run_directory,
            host=resolved_remote,
            remote_dir=remote_directory,
        )

    edge_groups = distinct_groups_from_session(session)
    return InitArtifacts(
        run_directory=run_directory,
        geometry_xyz=geometry_path,
        molconfig_py=molconfig_path,
        source_cif=source_cif,
        preview_png=preview_path,
        summary_json=summary_path,
        remote_host=resolved_remote,
        remote_directory=remote_directory,
        distinct_carbon_groups=edge_groups,
    )


def initialize_from_cif(
    run_directory: Path,
    cif_path: Path,
    *,
    name: str | None = None,
    edge_element: str = "C",
    block_name: str | None = None,
    remote_host: str | None = None,
    remote_root: str = "/home/hduva/projects/dft-runs",
    write_preview: bool = True,
    keep_all_sites: bool = False,
    one_ligand: bool = False,
) -> InitArtifacts:
    """Create a StoBe run directory from a CIF file (batch, non-interactive).

    Parameters
    ----------
    run_directory
        Destination directory for geometry, molConfig, and session files.
    cif_path
        Input CIF.
    name
        Molecule title; defaults to the run directory name.
    edge_element
        Core-edge element for generating sites.
    block_name
        Optional CIF data-block name.
    remote_host
        Optional SSH alias to rsync after writing.
    remote_root
        Remote base directory when ``remote_host`` is set.
    write_preview
        When True, write ``init_preview.png``.
    keep_all_sites
        When True, retain solvent/residual CIF fragments.
    one_ligand
        When True, enable generating core-edge sites on a single metal-bound
        organic wing; other ligands remain in the geometry.
    """
    run_directory = Path(run_directory)
    mname = name or run_directory.name
    session = create_setup_session(
        cif_path,
        mname=mname,
        edge_element=edge_element,
        block_name=block_name,
        keep_all_sites=keep_all_sites,
        one_ligand=one_ligand,
    )
    return finalize_setup_session(
        run_directory,
        session,
        write_preview=write_preview,
        remote_host=remote_host,
        remote_root=remote_root,
    )


def _meta_from_session(session: SetupSession) -> CifStructureMeta:
    from dftlearn.setup.types import CifAtomSite, CifStructureMeta

    cell = session.meta.get("cell", {})
    sites = tuple(
        CifAtomSite(
            label=session.cif_labels[i],
            element=element_symbol_from_row(session.base_rows[i]),
            x=float(session.base_rows[i][1]),
            y=float(session.base_rows[i][2]),
            z=float(session.base_rows[i][3]),
        )
        for i in range(len(session.base_rows))
    )
    return CifStructureMeta(
        block_name=session.block_name,
        space_group=session.meta.get("space_group"),
        cell_a=float(cell.get("a", 1.0)),
        cell_b=float(cell.get("b", 1.0)),
        cell_c=float(cell.get("c", 1.0)),
        cell_alpha=float(cell.get("alpha", 90.0)),
        cell_beta=float(cell.get("beta", 90.0)),
        cell_gamma=float(cell.get("gamma", 90.0)),
        sites=sites,
        symmetry_ops=tuple(session.meta.get("symmetry_ops", [])),
    )


def element_symbol_from_row(row: list[float | str]) -> str:
    """Return the element symbol from one session base row."""
    from dftlearn.io.xyz_structure import element_symbol_from_xyz_label

    return element_symbol_from_xyz_label(str(row[0]))


def _session_element_symbols(session: SetupSession) -> set[str]:
    """Return unique element symbols present in a setup session."""
    return {element_symbol_from_row(row) for row in session.base_rows}


def _build_summary_from_session(
    session: SetupSession,
    plan: InitGeometryPlan,
    cif_path: Path,
) -> dict:
    basis_info = {}
    for sym in plan.element_group_order:
        preset = stobe_basis_for_element(sym)
        basis_info[sym] = {
            "aux_ground": preset.aux_ground,
            "orbital_ground": preset.orbital_ground,
            "aux_excited": preset.aux_excited,
            "orbital_excited": preset.orbital_excited,
            "mcp": preset.mcp,
        }

    enabled_groups = [g for g in session.site_groups if g.enabled]
    source = session.cif_path or session.meta.get("smiles", "")
    blocks = cif_blocks(cif_path) if cif_path.is_file() else []
    return {
        "source_cif": source,
        "source_kind": session.meta.get("source_kind", "unknown"),
        "pubchem_cid": session.meta.get("pubchem_cid"),
        "smiles": session.meta.get("smiles"),
        "setup_session": "setup_session.json",
        "block_name": session.block_name,
        "available_blocks": blocks,
        "space_group": session.meta.get("space_group"),
        "symmetry_ops": session.meta.get("symmetry_ops", []),
        "cell": session.meta.get("cell", {}),
        "edge_element": session.edge_element,
        "n_core_sites": len(plan.generating_site_indices),
        "element_counts": plan.element_counts,
        "element_group_order": list(plan.element_group_order),
        "basis_sets": basis_info,
        "alignment": {
            "note": session.alignment_note,
            "axis": session.alignment_axis,
            "euler_deg": session.euler_deg,
            "alignment_matrix": session.alignment_matrix,
            "rotation_matrix": session.rotation_matrix,
        },
        "distinct_sites": [
            {
                "site": group.site_tag,
                "custom_name": group.custom_name,
                "xyz_label": group.xyz_label,
                "cif_label": group.cif_label,
                "class_size": len(group.atom_indices),
                "rank": group.rank,
                "enabled": group.enabled,
            }
            for group in enabled_groups
        ],
    }


def _build_summary(
    meta: CifStructureMeta,
    plan: InitGeometryPlan,
    cif_path: Path,
    block_name: str | None,
) -> dict:
    basis_info = {}
    for sym in plan.element_group_order:
        preset = stobe_basis_for_element(sym)
        basis_info[sym] = {
            "aux_ground": preset.aux_ground,
            "orbital_ground": preset.orbital_ground,
            "aux_excited": preset.aux_excited,
            "orbital_excited": preset.orbital_excited,
            "mcp": preset.mcp,
        }

    return {
        "source_cif": str(cif_path.resolve()),
        "block_name": meta.block_name,
        "requested_block": block_name,
        "available_blocks": cif_blocks(cif_path),
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
        "edge_element": plan.edge_element,
        "n_core_sites": len(plan.generating_site_indices),
        "element_counts": plan.element_counts,
        "element_group_order": list(plan.element_group_order),
        "basis_sets": basis_info,
        "distinct_sites": [
            {
                "site": f"{plan.edge_element}{i}",
                "xyz_label": g.xyz_label,
                "cif_label": g.cif_label,
                "class_size": len(g.atom_indices),
                "rank": g.rank,
            }
            for i, g in enumerate(plan.distinct_groups, start=1)
        ],
    }
