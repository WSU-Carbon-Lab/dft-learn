"""Matplotlib interactive editor for StoBe calculation setup sessions."""

from __future__ import annotations

import os
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Button, CheckButtons, RadioButtons, Slider, TextBox

from dftlearn.io.xyz_structure import element_symbol_from_xyz_label
from dftlearn.setup.alignment import (
    metal_ligand_axis,
    principal_inertia_axis,
    projection_indices,
    rotation_align_vector,
)
from dftlearn.setup.labeling import (
    merge_site_groups,
    relabel_site_group,
    renumber_site_tags,
    toggle_site_group,
)
from dftlearn.setup.session import (
    set_session_alignment_matrix,
    update_session_rotation,
)
from dftlearn.visualization.xyz_wireframe import draw_xyz_wireframe_on_ax

if TYPE_CHECKING:
    from collections.abc import Callable

    from dftlearn.setup.session import (
        SetupSession,
    )


class SetupInteractiveEditor:
    """Interactive matplotlib UI for reviewing and editing a setup session."""

    def __init__(
        self,
        session: SetupSession,
        *,
        on_save: Callable[[SetupSession], None] | None = None,
    ) -> None:
        self.session = session
        self.on_save = on_save
        self.selected_group_ids: set[int] = set()
        self.selected_atom_index: int | None = None
        self._group_colors = plt.get_cmap("tab20")

        self.fig = plt.figure(figsize=(13.5, 8.0), layout="constrained")
        gs = self.fig.add_gridspec(4, 3, height_ratios=[4.0, 0.35, 0.35, 0.45])

        self.ax_struct = self.fig.add_subplot(gs[0, :2])
        self.ax_info = self.fig.add_subplot(gs[0, 2])
        self.ax_info.axis("off")

        slider_gs = gs[1, :2].subgridspec(1, 3)
        self.ax_rx = self.fig.add_subplot(slider_gs[0, 0])
        self.ax_ry = self.fig.add_subplot(slider_gs[0, 1])
        self.ax_rz = self.fig.add_subplot(slider_gs[0, 2])

        self.ax_view = self.fig.add_subplot(gs[1, 2])
        self.ax_checks = self.fig.add_subplot(gs[2, :2])
        self.ax_checks.axis("off")

        btn_gs = gs[2, 2].subgridspec(3, 2, hspace=0.45, wspace=0.25)
        self.ax_merge = self.fig.add_subplot(btn_gs[0, 0])
        self.ax_toggle = self.fig.add_subplot(btn_gs[0, 1])
        self.ax_renumber = self.fig.add_subplot(btn_gs[1, 0])
        self.ax_inertia = self.fig.add_subplot(btn_gs[1, 1])
        self.ax_metal = self.fig.add_subplot(btn_gs[2, 0])
        self.ax_save = self.fig.add_subplot(btn_gs[2, 1])

        self.ax_tag = self.fig.add_subplot(gs[3, 0])
        self.ax_name = self.fig.add_subplot(gs[3, 1])
        self.ax_help = self.fig.add_subplot(gs[3, 2])
        self.ax_help.axis("off")

        self._build_widgets()
        self._draw_structure()
        self._refresh_info()
        self.fig.canvas.mpl_connect("button_press_event", self._on_click)

    def _build_widgets(self) -> None:
        rx0, ry0, rz0 = self.session.euler_deg
        self.slider_rx = Slider(
            self.ax_rx,
            "Rot X",
            -180.0,
            180.0,
            valinit=rx0,
            valstep=1.0,
        )
        self.slider_ry = Slider(
            self.ax_ry,
            "Rot Y",
            -180.0,
            180.0,
            valinit=ry0,
            valstep=1.0,
        )
        self.slider_rz = Slider(
            self.ax_rz,
            "Rot Z",
            -180.0,
            180.0,
            valinit=rz0,
            valstep=1.0,
        )
        for slider in (self.slider_rx, self.slider_ry, self.slider_rz):
            slider.on_changed(self._on_rotation)

        self.view_radio = RadioButtons(
            self.ax_view,
            ("xy", "xz", "yz"),
            active=("xy", "xz", "yz").index(self.session.view),
        )
        self.view_radio.on_clicked(self._on_view)

        labels, states = self._check_labels_states()
        self.check_groups = CheckButtons(self.ax_checks, labels, states)
        self.check_groups.on_clicked(self._on_check)

        self.btn_merge = Button(self.ax_merge, "Merge selected")
        self.btn_merge.on_clicked(self._on_merge)
        self.btn_toggle = Button(self.ax_toggle, "Toggle selected")
        self.btn_toggle.on_clicked(self._on_toggle)
        self.btn_renumber = Button(self.ax_renumber, "Renumber sites")
        self.btn_renumber.on_clicked(self._on_renumber)
        self.btn_inertia = Button(self.ax_inertia, "Align inertia Z")
        self.btn_inertia.on_clicked(self._on_align_inertia)
        self.btn_metal = Button(self.ax_metal, "Align metal-lig Z")
        self.btn_metal.on_clicked(self._on_align_metal)
        self.btn_save = Button(self.ax_save, "Save && close")
        self.btn_save.on_clicked(self._on_save)

        init_tag = self._selected_group().site_tag if self._selected_group() else "C1"
        init_name = self._selected_group().custom_name if self._selected_group() else ""
        self.text_tag = TextBox(self.ax_tag, "Site tag ", initial=init_tag)
        self.text_name = TextBox(self.ax_name, "Custom name ", initial=init_name)
        self.text_tag.on_submit(self._on_tag_submit)
        self.text_name.on_submit(self._on_name_submit)

        self.ax_help.text(
            0.0,
            1.0,
            "Click atoms to inspect basis sets.\n"
            "Check groups to select; merge to mark\n"
            "as chemically equivalent.\n"
            "Disabled groups are not generating sites.",
            va="top",
            fontsize=9,
        )

    def _enabled_groups(self) -> list:
        return [g for g in self.session.site_groups if g.enabled]

    def _check_labels_states(self) -> tuple[list[str], list[bool]]:
        labels: list[str] = []
        states: list[bool] = []
        for group in self.session.site_groups:
            if not group.enabled:
                continue
            labels.append(
                f"{group.site_tag} {group.custom_name} ({len(group.atom_indices)})"
            )
            states.append(group.group_id in self.selected_group_ids)
        if not labels:
            labels = ["(no enabled groups)"]
            states = [False]
        return labels, states

    def _selected_group(self):
        if len(self.selected_group_ids) == 1:
            gid = next(iter(self.selected_group_ids))
            for group in self.session.site_groups:
                if group.group_id == gid:
                    return group
        return None

    def _group_color(self, group_id: int) -> tuple[float, float, float]:
        rgb = self._group_colors((group_id - 1) % 20)
        return float(rgb[0]), float(rgb[1]), float(rgb[2])

    def _draw_structure(self) -> None:
        self.ax_struct.clear()
        rows = self.session.rotated_rows
        i_col, j_col = projection_indices(self.session.view)
        site_colors: dict[int, tuple[float, float, float]] = {}
        for group in self.session.site_groups:
            if not group.enabled:
                continue
            color = self._group_color(group.group_id)
            for idx in group.atom_indices:
                site_colors[idx] = color
        draw_rows: list[tuple[str, float, float, float]] = [
            (
                str(rows[i][0]),
                float(rows[i][i_col + 1]),
                float(rows[i][j_col + 1]),
                float(rows[i][3]),
            )
            for i in range(len(rows))
        ]
        draw_xyz_wireframe_on_ax(
            self.ax_struct,
            draw_rows,
            site_colors,
            show_hydrogen=False,
            plot_margins=0.05,
        )
        if self.selected_atom_index is not None:
            x = rows[self.selected_atom_index][i_col + 1]
            y = rows[self.selected_atom_index][j_col + 1]
            self.ax_struct.scatter(
                [x],
                [y],
                s=120,
                facecolors="none",
                edgecolors="black",
                linewidths=2.0,
                zorder=10,
            )
        axis = np.asarray(self.session.alignment_axis, dtype=np.float64)
        if np.linalg.norm(axis) > 1e-8:
            origin = np.mean(
                [[r[1], r[2], r[3]] for r in rows],
                axis=0,
            )
            vec = axis / np.linalg.norm(axis)
            i0, j0 = i_col, j_col
            self.ax_struct.annotate(
                "",
                xy=(origin[i0] + vec[i0], origin[j0] + vec[j0]),
                xytext=(origin[i0], origin[j0]),
                arrowprops={"arrowstyle": "->", "color": "crimson", "lw": 2.0},
            )
        axis_names = ("x", "y", "z")
        self.ax_struct.set_xlabel(f"{axis_names[i_col]} (A)")
        self.ax_struct.set_ylabel(f"{axis_names[j_col]} (A)")
        self.ax_struct.set_title(
            f"{self.session.mname} | {self.session.edge_element} sites: "
            f"{len(self._enabled_groups())} | view={self.session.view.upper()}"
        )
        self.fig.canvas.draw_idle()

    def _refresh_info(self) -> None:
        self.ax_info.clear()
        self.ax_info.axis("off")
        lines = [
            f"Molecule: {self.session.mname}",
            f"Edge element: {self.session.edge_element}",
            f"Alignment: {self.session.alignment_note}",
            "",
        ]
        if self.selected_atom_index is not None:
            idx = self.selected_atom_index
            label = self.session.rotated_rows[idx][0]
            sym = element_symbol_from_xyz_label(str(label))
            cif = self.session.cif_labels[idx]
            basis = self.session.basis_for_atom(idx)
            group = self.session.group_for_atom(idx)
            lines.extend(
                [
                    f"Atom index: {idx}",
                    f"XYZ label: {label}",
                    f"CIF label: {cif}",
                    f"Element: {sym}",
                    "",
                    "Basis (ground):",
                    f"  A: {basis['aux_ground']}",
                    f"  O: {basis['orbital_ground']}",
                    "Basis (excited):",
                    f"  A: {basis['aux_excited']}",
                    f"  O: {basis['orbital_excited']}",
                ]
            )
            if basis["mcp"]:
                lines.append(f"  MCP: {basis['mcp']}")
            if group is not None:
                lines.extend(
                    [
                        "",
                        f"Site group: {group.site_tag}",
                        f"Custom: {group.custom_name}",
                        f"Class size: {len(group.atom_indices)}",
                    ]
                )
        else:
            lines.append("Click an atom to inspect basis sets.")
        self.ax_info.text(
            0.0,
            1.0,
            "\n".join(lines),
            va="top",
            fontsize=9,
            family="monospace",
        )
        self.fig.canvas.draw_idle()

    def _rebuild_checks(self) -> None:
        self.ax_checks.clear()
        self.ax_checks.axis("off")
        labels, states = self._check_labels_states()
        self.check_groups = CheckButtons(self.ax_checks, labels, states)
        self.check_groups.on_clicked(self._on_check)

    def _on_click(self, event) -> None:
        if event.inaxes is not self.ax_struct or event.xdata is None:
            return
        rows = self.session.rotated_rows
        i_col, j_col = projection_indices(self.session.view)
        click = np.array([event.xdata, event.ydata], dtype=np.float64)
        best_idx: int | None = None
        best_dist = 0.45
        for idx, row in enumerate(rows):
            pos = np.array([row[i_col + 1], row[j_col + 1]], dtype=np.float64)
            dist = float(np.linalg.norm(pos - click))
            if dist < best_dist:
                best_dist = dist
                best_idx = idx
        if best_idx is None:
            return
        self.selected_atom_index = best_idx
        group = self.session.group_for_atom(best_idx)
        if group is not None:
            self.selected_group_ids = {group.group_id}
            self.text_tag.set_val(group.site_tag)
            self.text_name.set_val(group.custom_name)
            self._rebuild_checks()
        self._draw_structure()
        self._refresh_info()

    def _on_rotation(self, _val) -> None:
        update_session_rotation(
            self.session,
            self.slider_rx.val,
            self.slider_ry.val,
            self.slider_rz.val,
        )
        self._draw_structure()

    def _on_view(self, label: str | None) -> None:
        if label is None:
            return
        self.session.view = label
        self._draw_structure()

    def _on_check(self, label: str | None) -> None:
        if label is None or label.startswith("("):
            return
        site_tag = label.split()[0]
        for group in self.session.site_groups:
            if group.site_tag == site_tag:
                if group.group_id in self.selected_group_ids:
                    self.selected_group_ids.remove(group.group_id)
                else:
                    self.selected_group_ids.add(group.group_id)
                break
        selected = self._selected_group()
        if selected is not None:
            self.text_tag.set_val(selected.site_tag)
            self.text_name.set_val(selected.custom_name)

    def _on_merge(self, _event) -> None:
        if len(self.selected_group_ids) < 2:
            return
        try:
            self.session.site_groups = merge_site_groups(
                self.session.site_groups,
                sorted(self.selected_group_ids),
            )
        except ValueError:
            return
        self.selected_group_ids = {sorted(self.selected_group_ids)[0]}
        self.session.site_groups = renumber_site_tags(
            self.session.site_groups,
            self.session.edge_element,
        )
        self._rebuild_checks()
        self._draw_structure()
        self._refresh_info()

    def _on_toggle(self, _event) -> None:
        if len(self.selected_group_ids) != 1:
            return
        gid = next(iter(self.selected_group_ids))
        group = next(g for g in self.session.site_groups if g.group_id == gid)
        self.session.site_groups = toggle_site_group(
            self.session.site_groups,
            gid,
            enabled=not group.enabled,
        )
        self.session.site_groups = renumber_site_tags(
            self.session.site_groups,
            self.session.edge_element,
        )
        self._rebuild_checks()
        self._draw_structure()
        self._refresh_info()

    def _on_renumber(self, _event) -> None:
        self.session.site_groups = renumber_site_tags(
            self.session.site_groups,
            self.session.edge_element,
        )
        self._rebuild_checks()
        self._refresh_info()

    def _on_align_inertia(self, _event) -> None:
        raw = [
            (str(r[0]), float(r[1]), float(r[2]), float(r[3]))
            for r in self.session.base_rows
        ]
        pos = np.array([[r[1], r[2], r[3]] for r in raw], dtype=np.float64)
        _small, _mid, largest = principal_inertia_axis(pos)
        matrix = rotation_align_vector(largest, np.array([0.0, 0.0, 1.0]))
        set_session_alignment_matrix(
            self.session,
            matrix,
            note="principal inertia -> Z",
        )
        self.slider_rx.set_val(0.0)
        self.slider_ry.set_val(0.0)
        self.slider_rz.set_val(0.0)
        self._draw_structure()
        self._refresh_info()

    def _on_align_metal(self, _event) -> None:
        raw = [
            (str(r[0]), float(r[1]), float(r[2]), float(r[3]))
            for r in self.session.base_rows
        ]
        axis = metal_ligand_axis(raw)
        matrix = rotation_align_vector(axis, np.array([0.0, 0.0, 1.0]))
        set_session_alignment_matrix(
            self.session,
            matrix,
            note="metal-ligand axis -> Z",
        )
        self.session.alignment_axis = [0.0, 0.0, 1.0]
        self.slider_rx.set_val(0.0)
        self.slider_ry.set_val(0.0)
        self.slider_rz.set_val(0.0)
        self._draw_structure()
        self._refresh_info()

    def _on_tag_submit(self, text: str) -> None:
        if len(self.selected_group_ids) != 1:
            return
        gid = next(iter(self.selected_group_ids))
        self.session.site_groups = relabel_site_group(
            self.session.site_groups,
            gid,
            site_tag=text.strip(),
        )
        self._rebuild_checks()
        self._refresh_info()

    def _on_name_submit(self, text: str) -> None:
        if len(self.selected_group_ids) != 1:
            return
        gid = next(iter(self.selected_group_ids))
        self.session.site_groups = relabel_site_group(
            self.session.site_groups,
            gid,
            custom_name=text.strip(),
        )
        self._rebuild_checks()
        self._refresh_info()

    def _on_save(self, _event) -> None:
        if self.on_save is not None:
            self.on_save(self.session)
        plt.close(self.fig)

    def run(self) -> SetupSession:
        """Block until the user closes the editor window."""
        plt.show()
        return self.session


def interactive_display_available() -> bool:
    """Return True when an interactive Matplotlib backend is likely available."""
    backend = plt.get_backend().lower()
    return not ("agg" in backend and os.environ.get("DISPLAY", "") == "")


def launch_setup_editor(
    session: SetupSession,
    *,
    on_save: Callable[[SetupSession], None],
) -> SetupSession:
    """Open the setup editor, invoking ``on_save`` when the user saves."""
    if not interactive_display_available():
        msg = (
            "Interactive setup requires a display (set DISPLAY or MPLBACKEND). "
            "Use batch init or run setup with X forwarding."
        )
        raise RuntimeError(msg)
    editor = SetupInteractiveEditor(session, on_save=on_save)
    return editor.run()
