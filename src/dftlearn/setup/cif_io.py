"""CIF ingestion: block selection, cell metadata, and Cartesian coordinates."""

from __future__ import annotations

import re
from pathlib import Path

import gemmi

from dftlearn.setup.types import CifAtomSite, CifStructureMeta

_FLOAT_RE = re.compile(r"^[\d.+-]+")


def _parse_cif_float(raw: str) -> float:
    token = raw.strip().split("(")[0].strip().strip("'\"")
    if not _FLOAT_RE.match(token):
        msg = f"Could not parse CIF float from {raw!r}"
        raise ValueError(msg)
    return float(token)


def _block_tag(block: gemmi.cif.Block, tag: str) -> str | None:
    value = block.find_value(tag)
    if value == "?":
        return None
    text = str(value).strip().strip("'\"")
    return text or None


def _symmetry_ops_from_block(block: gemmi.cif.Block) -> tuple[str, ...]:
    col = block.find_loop("_symmetry_equiv_pos_as_xyz")
    if col is None:
        return ()
    return tuple(str(col[i]).strip() for i in range(len(col)))


def cif_blocks(cif_path: Path) -> list[tuple[str, int]]:
    """List CIF data blocks with non-empty atom-site tables.

    Parameters
    ----------
    cif_path
        Path to a ``.cif`` file readable by gemmi.

    Returns
    -------
    list[tuple[str, int]]
        ``(block_name, site_count)`` sorted by descending site count.
    """
    cif_path = Path(cif_path)
    doc = gemmi.cif.read(str(cif_path))
    scored: list[tuple[str, int]] = []
    for block in doc:
        try:
            st = gemmi.make_small_structure_from_block(block)
        except (RuntimeError, ValueError, IndexError):
            continue
        n_sites = len(st.sites)
        if n_sites > 0:
            scored.append((block.name, n_sites))
    scored.sort(key=lambda item: item[1], reverse=True)
    return scored


def load_cif_structure(
    cif_path: Path,
    *,
    block_name: str | None = None,
) -> CifStructureMeta:
    """Load one CIF block into Cartesian angstrom coordinates.

    When ``block_name`` is omitted, selects the block with the largest
    atom-site loop so powder multi-block files prefer the main phase.

    Parameters
    ----------
    cif_path
        Path to a CIF file.
    block_name
        Optional explicit block name; must exist and contain atom sites.

    Returns
    -------
    CifStructureMeta
        Cell parameters, optional space group, symmetry operations, and sites.

    Raises
    ------
    FileNotFoundError
        If ``cif_path`` is missing.
    ValueError
        If no suitable block or atom sites are found.
    """
    cif_path = Path(cif_path)
    if not cif_path.is_file():
        msg = f"CIF file not found: {cif_path}"
        raise FileNotFoundError(msg)

    doc = gemmi.cif.read(str(cif_path))
    if block_name is not None:
        block = doc[block_name]
        blocks = [block]
    else:
        ranked = cif_blocks(cif_path)
        if not ranked:
            msg = f"No atom-site blocks found in {cif_path}"
            raise ValueError(msg)
        blocks = [doc[ranked[0][0]]]

    block = blocks[0]
    st = gemmi.make_small_structure_from_block(block)
    if not st.sites:
        msg = f"Block {block.name!r} has no atom sites in {cif_path}"
        raise ValueError(msg)

    cell = st.cell
    sites: list[CifAtomSite] = []
    for site in st.sites:
        pos = cell.orthogonalize(site.fract)
        element = site.element.name.capitalize()
        if element == "D":
            element = "H"
        sites.append(
            CifAtomSite(
                label=site.label,
                element=element,
                x=float(pos.x),
                y=float(pos.y),
                z=float(pos.z),
            )
        )

    return CifStructureMeta(
        block_name=block.name,
        space_group=_block_tag(block, "_symmetry_space_group_name_H-M"),
        cell_a=_parse_cif_float(_block_tag(block, "_cell_length_a") or "1"),
        cell_b=_parse_cif_float(_block_tag(block, "_cell_length_b") or "1"),
        cell_c=_parse_cif_float(_block_tag(block, "_cell_length_c") or "1"),
        cell_alpha=_parse_cif_float(_block_tag(block, "_cell_angle_alpha") or "90"),
        cell_beta=_parse_cif_float(_block_tag(block, "_cell_angle_beta") or "90"),
        cell_gamma=_parse_cif_float(_block_tag(block, "_cell_angle_gamma") or "90"),
        sites=tuple(sites),
        symmetry_ops=_symmetry_ops_from_block(block),
    )
