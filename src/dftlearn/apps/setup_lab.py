"""Streamlit setup lab: PubChem search, 3D build, relaxation, and site labeling."""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import py3Dmol
import streamlit as st
import streamlit.components.v1 as components

from dftlearn.setup.alignment import (
    align_bond_to_axis,
    align_frame_from_bond,
    metal_ligand_axis,
    principal_inertia_axis,
    rotation_align_vector,
)
from dftlearn.setup.labeling import (
    merge_site_groups,
    relabel_site_group,
    renumber_site_tags,
    toggle_site_group,
)
from dftlearn.setup.pubchem_client import search_pubchem
from dftlearn.setup.session import (
    create_setup_session_from_rows,
    load_setup_session,
    save_setup_session,
    set_session_alignment_matrix,
    update_session_rotation,
)
from dftlearn.setup.structure_build import (
    RelaxMethod,
    build_3d_from_smiles,
    rows_to_xyz_text,
)
from dftlearn.setup.structure_viz import (
    atom_display_table,
    bond_neighbors,
    display_bonds,
    format_bond_label,
    minimal_generating_indices,
    session_mol,
    store_mol_block,
)
from dftlearn.setup.workflow import finalize_setup_session

if TYPE_CHECKING:
    from dftlearn.setup.pubchem_client import PubChemCompound
    from dftlearn.setup.session import SetupSession

_PALETTE = [
    "#e6194b",
    "#3cb44b",
    "#4363d8",
    "#f58231",
    "#911eb4",
    "#42d4f4",
    "#f032e6",
    "#bfef45",
    "#fabed4",
    "#469990",
]

_BASIS_FIELDS = (
    ("aux_ground", "A (ground)"),
    ("orbital_ground", "O (ground)"),
    ("aux_excited", "A (excited)"),
    ("orbital_excited", "O (excited)"),
    ("mcp", "MCP"),
)


def _init_state(run_directory: Path | None) -> None:
    defaults = {
        "step": "search",
        "run_directory": str(run_directory) if run_directory else "",
        "query": "C27H18AlN3O3",
        "compounds": [],
        "selected_cid": None,
        "session": None,
        "edge_element": "C",
        "relax_method": RelaxMethod.ANNEAL.value,
        "relax_steps": 500,
        "selected_atom": 0,
        "bond_atom_a": 0,
        "bond_atom_b": 1,
        "plane_atom": 2,
        "show_minimal_generating": False,
        "merge_ids": [],
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value


def _show_structure(
    rows: list[tuple[str, float, float, float]],
    site_groups,
    *,
    mol,
    bonds: list[tuple[int, int]],
    selected_atom: int | None,
    bond_atoms: tuple[int, int] | None,
    minimal_only: bool,
) -> None:
    """Render the structure in py3Dmol with explicit per-bond sticks."""
    xyz = rows_to_xyz_text(list(rows), "setup-lab")
    view = py3Dmol.view(width=780, height=520)
    view.addModel(xyz, "xyz")

    visible = set(range(len(rows)))
    if minimal_only:
        visible = set(minimal_generating_indices(site_groups))

    visible_bonds = [
        (i, j) for i, j in bonds if i in visible and j in visible
    ]
    for i, j in visible_bonds:
        view.addStyle(
            {"bonds": [{"atom1": i, "atom2": j}]},
            {"stick": {"radius": 0.12, "color": "#707070"}},
        )

    view.setStyle({"sphere": {"scale": 0.18}})

    for idx in range(len(rows)):
        if idx not in visible:
            view.setStyle(
                {"serial": idx + 1},
                {"sphere": {"scale": 0.01, "opacity": 0.05}},
            )

    for group in site_groups:
        if not group.enabled:
            continue
        color = _PALETTE[(group.group_id - 1) % len(_PALETTE)]
        for idx in group.atom_indices:
            if idx not in visible:
                continue
            view.setStyle(
                {"serial": idx + 1},
                {"sphere": {"color": color, "scale": 0.28}},
            )

    if bond_atoms is not None:
        i, j = bond_atoms
        if i in visible and j in visible:
            view.addStyle(
                {"serial": i + 1},
                {"sphere": {"color": "#ffd700", "scale": 0.34}},
            )
            view.addStyle(
                {"serial": j + 1},
                {"sphere": {"color": "#ffd700", "scale": 0.34}},
            )
            view.addStyle(
                {"bonds": [{"atom1": i, "atom2": j}]},
                {"stick": {"color": "#ffd700", "radius": 0.22}},
            )

    if selected_atom is not None and selected_atom in visible:
        view.addStyle(
            {"serial": selected_atom + 1},
            {"sphere": {"color": "black", "scale": 0.36}},
        )

    view.zoomTo()
    components.html(view._make_html(), height=540, scrolling=False)


def _groups_dataframe(session: SetupSession) -> pd.DataFrame:
    rows = []
    for group in session.site_groups:
        rows.append(
            {
                "enabled": group.enabled,
                "group_id": group.group_id,
                "site_tag": group.site_tag,
                "custom_name": group.custom_name,
                "cif_label": group.cif_label,
                "class_size": len(group.atom_indices),
                "rank": group.rank,
            }
        )
    return pd.DataFrame(rows)


def _render_search_step() -> None:
    st.subheader("1. PubChem search")
    query = st.text_input("Compound name or formula", key="query")
    st.selectbox("Core-edge element", ["C", "N", "O", "Al"], key="edge_element")
    if st.button("Search PubChem", type="primary"):
        with st.spinner("Querying PubChem..."):
            st.session_state.compounds = search_pubchem(query)
        if not st.session_state.compounds:
            st.warning("No PubChem hits. Try a formula like C27H18AlN3O3.")
            return

    compounds: list[PubChemCompound] = st.session_state.compounds
    if compounds:
        options = {
            f"CID {c.cid} | {c.title} | {c.molecular_formula}": c.cid
            for c in compounds
        }
        label = st.selectbox("Select compound", list(options.keys()))
        st.session_state.selected_cid = options[label]
        selected = next(c for c in compounds if c.cid == st.session_state.selected_cid)
        st.code(selected.smiles, language=None)
        if st.button("Build 3D and continue"):
            st.session_state.step = "relax"
            st.session_state.selected_compound = selected
            st.rerun()


def _render_relax_step() -> None:
    st.subheader("2. Build 3D and relax")
    compound: PubChemCompound = st.session_state.selected_compound
    st.write(f"**{compound.title}** (CID {compound.cid}, {compound.molecular_formula})")
    method = st.selectbox(
        "Relaxation",
        [RelaxMethod.UFF.value, RelaxMethod.ANNEAL.value, RelaxMethod.NONE.value],
        format_func=lambda x: {
            RelaxMethod.UFF.value: "UFF minimization",
            RelaxMethod.ANNEAL.value: "Annealed UFF (thermal sampling proxy)",
            RelaxMethod.NONE.value: "Skip relaxation",
        }[x],
        key="relax_method",
    )
    steps = st.slider(
        "Minimization iterations",
        100,
        2000,
        500,
        step=50,
        key="relax_steps",
    )
    if st.button("Generate structure", type="primary"):
        with st.spinner("Embedding and relaxing..."):
            built = build_3d_from_smiles(
                compound.smiles,
                relax=RelaxMethod(method),
                relax_steps=steps,
            )
            run_name = st.session_state.run_directory or compound.title.replace(
                " ",
                "-",
            )[:40]
            meta = {
                "source_kind": "pubchem",
                "pubchem_cid": compound.cid,
                "pubchem_title": compound.title,
                "smiles": compound.smiles,
                "molecular_formula": compound.molecular_formula,
                "relax_method": method,
                "relax_steps": steps,
                "final_energy": built.final_energy,
            }
            session = create_setup_session_from_rows(
                built.rows,
                mname=run_name,
                edge_element=st.session_state.edge_element,
                atom_labels=built.atom_labels,
                source_meta=meta,
                mol=built.mol,
            )
            st.session_state.session = session
            st.session_state.step = "label"
            st.session_state.selected_atom = 0
            if len(built.rows) > 1:
                st.session_state.bond_atom_b = 1
            st.rerun()

    if st.button("Back to search"):
        st.session_state.step = "search"
        st.rerun()


def _apply_group_edits(session: SetupSession, edited: pd.DataFrame) -> SetupSession:
    by_id = {g.group_id: g for g in session.site_groups}
    for _, row in edited.iterrows():
        gid = int(row["group_id"])
        group = by_id[gid]
        enabled = bool(row["enabled"])
        if enabled != group.enabled:
            session.site_groups = toggle_site_group(
                session.site_groups,
                gid,
                enabled=enabled,
            )
        session.site_groups = relabel_site_group(
            session.site_groups,
            gid,
            site_tag=str(row["site_tag"]),
            custom_name=str(row["custom_name"]),
        )
    session.site_groups = renumber_site_tags(session.site_groups, session.edge_element)
    return session


def _sync_bond_from_atom(session: SetupSession, mol, atom_idx: int) -> None:
    neighbors = bond_neighbors(mol, atom_idx)
    if not neighbors:
        return
    st.session_state.bond_atom_a = atom_idx
    if st.session_state.bond_atom_b not in neighbors:
        st.session_state.bond_atom_b = neighbors[0]


def _render_basis_editor(session: SetupSession, atom_idx: int) -> None:
    basis = session.basis_for_atom(atom_idx)
    st.markdown(f"**Basis set — atom {atom_idx} ({basis['element']})**")
    with st.form(f"basis_form_{atom_idx}"):
        values: dict[str, str] = {}
        for field, label in _BASIS_FIELDS:
            current = basis.get(field) or ""
            values[field] = st.text_input(
                label,
                value=str(current),
                key=f"basis_{field}_{atom_idx}",
            )
        reset = st.form_submit_button("Reset to element default")
        save = st.form_submit_button("Save basis overrides", type="primary")
    if reset:
        session.clear_basis_override(atom_idx)
        st.rerun()
    if save:
        for field, _label in _BASIS_FIELDS:
            session.set_basis_override(atom_idx, field, values[field])
        st.success(f"Saved basis overrides for atom {atom_idx}.")
        st.rerun()


def _render_label_step() -> None:
    st.subheader("3. Label generating sites")
    session: SetupSession = st.session_state.session
    if session is None:
        st.error("No session loaded.")
        return

    rows = session.rotated_rows
    mol = session_mol(session)
    n_atoms = len(rows)
    structure_bonds = display_bonds(list(rows), mol)

    view_cols = st.columns([1.4, 1.0])
    with view_cols[0]:
        st.checkbox(
            "Show minimal generating structure (one representative per site group)",
            key="show_minimal_generating",
        )
        atom_table = pd.DataFrame(atom_display_table(rows, session.site_groups))
        selection = st.dataframe(
            atom_table,
            hide_index=True,
            use_container_width=True,
            on_select="rerun",
            selection_mode="single-row",
            key="atom_table",
        )
        selected_rows = selection.selection.rows  # ty: ignore[unresolved-attribute]
        if selected_rows:
            st.session_state.selected_atom = int(
                atom_table.iloc[selected_rows[0]]["index"]
            )
        atom_idx = int(st.session_state.selected_atom)
        _sync_bond_from_atom(session, mol, atom_idx)

        bond_i = int(st.session_state.bond_atom_a)
        bond_j = int(st.session_state.bond_atom_b)
        all_bonds = structure_bonds
        bond_options = {
            format_bond_label(rows, i, j): (i, j) for i, j in all_bonds
        }
        if bond_options:
            bond_labels = list(bond_options.keys())
            default_label = format_bond_label(rows, bond_i, bond_j)
            if default_label not in bond_options:
                default_label = bond_labels[0]
            picked = st.selectbox(
                "Selected bond",
                bond_labels,
                index=bond_labels.index(default_label),
            )
            bond_i, bond_j = bond_options[picked]
            st.session_state.bond_atom_a = bond_i
            st.session_state.bond_atom_b = bond_j
        else:
            st.warning("No bonds inferred for this structure.")

        _show_structure(
            rows,
            session.site_groups,
            mol=mol,
            bonds=structure_bonds,
            selected_atom=atom_idx,
            bond_atoms=(bond_i, bond_j),
            minimal_only=bool(st.session_state.show_minimal_generating),
        )

    with view_cols[1]:
        group = session.group_for_atom(atom_idx)
        if group is not None:
            st.write(f"Site group: **{group.site_tag}** ({group.custom_name})")
        _render_basis_editor(session, atom_idx)

        st.markdown("**Bond alignment**")
        st.caption("Pick a bond in the table area, then align it to a lab axis.")
        plane_atom = st.number_input(
            "Plane reference atom (for XY frame)",
            min_value=0,
            max_value=max(0, n_atoms - 1),
            value=int(st.session_state.plane_atom),
            step=1,
        )
        st.session_state.plane_atom = int(plane_atom)

        align_row1 = st.columns(3)
        if align_row1[0].button("Bond -> X"):
            matrix = align_bond_to_axis(rows, bond_i, bond_j, "x")
            set_session_alignment_matrix(
                session,
                matrix,
                note=f"bond {bond_i}-{bond_j} -> X",
            )
            st.rerun()
        if align_row1[1].button("Bond -> Y"):
            matrix = align_bond_to_axis(rows, bond_i, bond_j, "y")
            set_session_alignment_matrix(
                session,
                matrix,
                note=f"bond {bond_i}-{bond_j} -> Y",
            )
            st.rerun()
        if align_row1[2].button("Bond -> Z"):
            matrix = align_bond_to_axis(rows, bond_i, bond_j, "z")
            set_session_alignment_matrix(
                session,
                matrix,
                note=f"bond {bond_i}-{bond_j} -> Z",
            )
            st.rerun()

        if st.button("Set XY frame from bond (bond -> X, plane atom in XY)"):
            matrix = align_frame_from_bond(rows, bond_i, bond_j, int(plane_atom))
            set_session_alignment_matrix(
                session,
                matrix,
                note=f"bond {bond_i}-{bond_j} -> X; atom {plane_atom} in XY",
            )
            st.rerun()

        st.markdown("**Manual rotation**")
        with st.form("alignment_form"):
            rx = st.slider("Rot X", -180, 180, int(session.euler_deg[0]))
            ry = st.slider("Rot Y", -180, 180, int(session.euler_deg[1]))
            rz = st.slider("Rot Z", -180, 180, int(session.euler_deg[2]))
            apply_rot = st.form_submit_button("Apply rotation")
        if apply_rot:
            update_session_rotation(session, float(rx), float(ry), float(rz))
            st.rerun()

        preset_cols = st.columns(2)
        if preset_cols[0].button("Align inertia to Z"):
            pos = np.array([[r[1], r[2], r[3]] for r in rows], dtype=np.float64)
            _small, _mid, largest = principal_inertia_axis(pos)
            matrix = rotation_align_vector(largest, np.array([0.0, 0.0, 1.0]))
            set_session_alignment_matrix(session, matrix, note="principal inertia -> Z")
            st.rerun()
        if preset_cols[1].button("Align metal-ligand to Z"):
            matrix = rotation_align_vector(
                metal_ligand_axis(rows),
                np.array([0.0, 0.0, 1.0]),
            )
            set_session_alignment_matrix(session, matrix, note="metal-ligand axis -> Z")
            st.rerun()

        if session.alignment_note:
            st.caption(f"Current alignment: {session.alignment_note}")

    st.markdown("**Distinct site groups**")
    df = _groups_dataframe(session)
    edited = st.data_editor(
        df,
        disabled=["group_id", "cif_label", "class_size", "rank"],
        hide_index=True,
        use_container_width=True,
    )
    session = _apply_group_edits(session, edited)

    merge_ids = st.multiselect(
        "Merge groups (mark as chemically equivalent)",
        options=[g.group_id for g in session.site_groups if g.enabled],
        format_func=lambda gid: next(
            g.site_tag for g in session.site_groups if g.group_id == gid
        ),
    )
    if st.button("Merge selected groups") and len(merge_ids) >= 2:
        session.site_groups = merge_site_groups(session.site_groups, merge_ids)
        session.site_groups = renumber_site_tags(
            session.site_groups,
            session.edge_element,
        )
        st.session_state.session = session
        st.rerun()

    store_mol_block(session, mol)
    st.session_state.session = session

    run_dir_text = st.text_input(
        "Run directory",
        value=st.session_state.run_directory or f"./{session.mname}",
    )
    if st.button("Save run directory", type="primary"):
        run_dir = Path(run_dir_text).expanduser().resolve()
        run_dir.mkdir(parents=True, exist_ok=True)
        save_setup_session(run_dir, session)
        artifacts = finalize_setup_session(run_dir, session, write_preview=True)
        summary = json.loads(artifacts.summary_json.read_text(encoding="utf-8"))
        st.success(
            f"Wrote {artifacts.geometry_xyz.name}, {artifacts.molconfig_py.name}, "
            f"{summary['n_core_sites']} core sites."
        )
        st.code(f"uv run dftrun build {run_dir}", language="bash")


def run_setup_lab(run_directory: Path | None = None) -> None:
    """Launch the Streamlit setup lab application."""
    import os

    st.set_page_config(page_title="dftrun setup lab", layout="wide")
    st.title("dftrun setup lab")
    st.caption("PubChem -> 3D structure -> relaxation -> distinct-site labeling")

    env_run = os.environ.get("DFTRUN_SETUP_RUN_DIR", "").strip()
    if run_directory is None and env_run:
        run_directory = Path(env_run)

    _init_state(run_directory)

    if run_directory and run_directory.is_dir():
        session_file = run_directory / "setup_session.json"
        if session_file.is_file() and st.session_state.session is None:
            st.session_state.session = load_setup_session(run_directory)
            st.session_state.step = "label"
            st.session_state.run_directory = str(run_directory)

    step = st.session_state.step
    if step == "search":
        _render_search_step()
    elif step == "relax":
        _render_relax_step()
    elif step == "label":
        _render_label_step()
    else:
        st.session_state.step = "search"
        st.rerun()


if __name__ == "__main__":
    run_setup_lab()
