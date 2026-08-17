"""Shared Matplotlib defaults for StoBe spectroscopy summary figures."""

from __future__ import annotations

from typing import TYPE_CHECKING

from matplotlib.ticker import AutoMinorLocator

if TYPE_CHECKING:
    from matplotlib.axes import Axes

PUB_FS = 11
TICK_FS = 10
FIGURE_RC: dict[str, object] = {
    "font.size": TICK_FS,
    "axes.labelsize": PUB_FS,
    "axes.titlesize": PUB_FS,
    "xtick.labelsize": TICK_FS,
    "ytick.labelsize": TICK_FS,
    "legend.fontsize": TICK_FS,
    "axes.unicode_minus": False,
}


def style_pub_ax(ax: Axes) -> None:
    """Apply inward ticks, full spines, and light grids."""
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(
        axis="both",
        which="major",
        direction="in",
        top=True,
        right=True,
        length=5,
        width=0.9,
        labelsize=TICK_FS,
    )
    ax.tick_params(
        axis="both",
        which="minor",
        direction="in",
        top=True,
        right=True,
        length=2.5,
        width=0.65,
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.9)
    ax.set_axisbelow(True)
    ax.grid(which="major", linestyle="-", linewidth=0.55, color="0.82")
    ax.grid(which="minor", linestyle=":", linewidth=0.4, color="0.88")


def style_diag_ax(ax: Axes, *, right: bool = False) -> None:
    """Style compact diagnostic axes with outward ticks and no grid."""
    ax.tick_params(
        axis="both",
        which="major",
        direction="out",
        top=False,
        right=right,
        length=4.0,
        width=0.8,
        labelsize=TICK_FS,
        pad=3.0,
    )
    ax.tick_params(axis="both", which="minor", length=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(right)
    for name in ("left", "bottom"):
        ax.spines[name].set_visible(True)
        ax.spines[name].set_linewidth(0.8)
    ax.set_axisbelow(True)
    ax.grid(False)


def panel_label(ax: Axes, letter: str) -> None:
    """Place a panel letter at the top-left corner of the axes."""
    ax.text(
        0.02,
        0.98,
        f"({letter})",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontweight="normal",
        fontsize=PUB_FS,
        color="black",
        clip_on=False,
        zorder=10,
    )
