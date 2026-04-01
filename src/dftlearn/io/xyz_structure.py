"""XYZ geometry: StoBe-style labels, element symbols, and RDKit molecules."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


def element_symbol_from_xyz_label(label: str) -> str:
    """Map an XYZ atom label (for example ``C01``, ``Zn01``) to an RDKit element symbol.

    Leading alphabetic characters define the symbol; two-letter symbols are tried
    before one-letter symbols when both are valid.

    Parameters
    ----------
    label : str
        First field of an XYZ line.

    Returns
    -------
    str
        Canonical element symbol accepted by RDKit (for example ``C``, ``Zn``).

    Raises
    ------
    ValueError
        If no alphabetic prefix exists or no valid symbol can be resolved.
    """
    stripped = label.strip()
    match = re.match(r"^([A-Za-z]+)", stripped)
    if not match:
        msg = f"No element prefix in XYZ label {label!r}"
        raise ValueError(msg)
    alpha = match.group(1)
    table = Chem.GetPeriodicTable()
    for take in (2, 1):
        if len(alpha) >= take:
            sym = alpha[:take].capitalize()
            if table.GetAtomicNumber(sym) > 0:
                return sym
    msg = f"Could not resolve element symbol from XYZ label {label!r}"
    raise ValueError(msg)


def xyz_rows_from_file(path: Path) -> list[tuple[str, float, float, float]]:
    """Parse an XYZ file into atom labels and Cartesian coordinates in angstroms.

    Supports both standard XYZ (first line atom count, second comment, then atoms)
    and headerless StoBe-style files. If the first non-empty line is a single
    integer and at least two more lines exist, the first two lines are skipped;
    otherwise every line is scanned. Atom order in the file defines the index
    sequence used for site mapping (``C01``, ``C``, etc. are treated only as
    labels for element typing).

    Parameters
    ----------
    path : pathlib.Path
        Path to an ``.xyz`` file.

    Returns
    -------
    list[tuple[str, float, float, float]]
        ``(label, x, y, z)`` for each atom row, in file order.

    Raises
    ------
    FileNotFoundError
        If ``path`` is missing.
    ValueError
        If no atom rows were parsed.
    """
    path = Path(path)
    if not path.is_file():
        msg = f"XYZ file not found: {path}"
        raise FileNotFoundError(msg)
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if not lines:
        msg = f"Empty or whitespace-only XYZ file: {path}"
        raise ValueError(msg)
    start = 0
    first_tokens = lines[0].split()
    if len(first_tokens) == 1 and first_tokens[0].isdigit() and len(lines) >= 3:
        start = 2
    rows: list[tuple[str, float, float, float]] = []
    for raw in lines[start:]:
        parts = raw.split()
        if len(parts) < 4:
            continue
        label, xs, ys, zs = parts[0], parts[1], parts[2], parts[3]
        try:
            rows.append((label, float(xs), float(ys), float(zs)))
        except ValueError:
            continue
    if not rows:
        msg = f"No atom rows parsed from {path}"
        raise ValueError(msg)
    return rows


def mol_from_xyz_file(path: Path) -> Chem.Mol:
    """Build a molecule with a single 3D conformer from a StoBe-style XYZ file.

    Parameters
    ----------
    path : pathlib.Path
        XYZ path readable by :func:`xyz_rows_from_file`.

    Returns
    -------
    rdkit.Chem.Mol
        Sanitized molecule; conformer positions match the XYZ file.

    Raises
    ------
    ValueError
        If RDKit sanitization fails.
    """
    rw = Chem.RWMol()
    conf = Chem.Conformer()
    for _idx, (label, x, y, z) in enumerate(xyz_rows_from_file(path)):
        sym = element_symbol_from_xyz_label(label)
        atom_idx = rw.AddAtom(Chem.Atom(sym))
        conf.SetAtomPosition(atom_idx, (x, y, z))
    mol = rw.GetMol()
    mol.AddConformer(conf, assignId=True)
    try:
        Chem.SanitizeMol(mol)
    except Chem.MolSanitizeException as exc:
        msg = f"RDKit could not sanitize molecule from {path}: {exc}"
        raise ValueError(msg) from exc
    return mol


def site_label_to_atom_index_from_rows(
    site_tag: str,
    rows: list[tuple[str, float, float, float]],
) -> int:
    """Map a site tag to a 0-based index using the same row list as the geometry.

    ``C1`` selects the first atom whose element symbol is ``C`` in ``rows`` order,
    ``C2`` the second carbon, and ``N1`` the first nitrogen, independent of XYZ
    labels such as ``C01`` or bare ``C``.

    Parameters
    ----------
    site_tag : str
        Directory-style tag such as ``C1`` or ``N12``.
    rows : list[tuple[str, float, float, float]]
        Atom rows in file order (for example from :func:`xyz_rows_from_file`).

    Returns
    -------
    int
        Index into ``rows`` for that core site.

    Raises
    ------
    ValueError
        If ``site_tag`` is malformed, ``rows`` is empty, or the ordinal is out of
        range for that element.
    """
    if not rows:
        msg = "rows is empty"
        raise ValueError(msg)
    match = re.match(r"^([A-Za-z]+)(\d+)$", site_tag.strip())
    if not match:
        msg = f"Site tag must look like 'C12', got {site_tag!r}"
        raise ValueError(msg)
    elem_raw, num_s = match.group(1), match.group(2)
    ordinal = int(num_s)
    if ordinal < 1:
        msg = f"Site ordinal must be >= 1 in {site_tag!r}"
        raise ValueError(msg)
    table = Chem.GetPeriodicTable()
    target_sym = None
    for take in (2, 1):
        if len(elem_raw) >= take:
            sym = elem_raw[:take].capitalize()
            if table.GetAtomicNumber(sym) > 0:
                target_sym = sym
                break
    if target_sym is None:
        msg = f"Unknown element prefix in site tag {site_tag!r}"
        raise ValueError(msg)
    count = 0
    for idx, (label, _x, _y, _z) in enumerate(rows):
        if element_symbol_from_xyz_label(label) == target_sym:
            count += 1
            if count == ordinal:
                return idx
    msg = f"Ordinal {ordinal} for element {target_sym} not found (site {site_tag!r})"
    raise ValueError(msg)


def site_label_to_atom_index(site_tag: str, xyz_path: Path) -> int:
    """Map a StoBe site folder tag to a 0-based XYZ atom index (file order).

    Equivalent to :func:`site_label_to_atom_index_from_rows` with rows from
    ``xyz_rows_from_file(xyz_path)``.

    Parameters
    ----------
    site_tag : str
        Directory name such as ``C1`` or ``N2``.
    xyz_path : pathlib.Path
        Geometry file paired with the StoBe run.

    Returns
    -------
    int
        Atom index for highlighting or wireframe shading.

    Raises
    ------
    ValueError
        If ``site_tag`` is malformed or the ordinal is out of range.
    FileNotFoundError
        If ``xyz_path`` is missing.
    """
    return site_label_to_atom_index_from_rows(site_tag, xyz_rows_from_file(xyz_path))


def site_highlight_geometry(
    mol: Chem.Mol,
    atom_index: int,
    size: tuple[int, int] = (320, 320),
) -> np.ndarray:
    """Render a 2D depiction with one highlighted atom to an RGBA raster.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Input molecule; coordinates are replaced by RDKit 2D layout for drawing.
    atom_index : int
        0-based atom index to highlight.
    size : tuple[int, int], optional
        Image width and height in pixels.

    Returns
    -------
    numpy.ndarray
        UInt8 RGBA image array suitable for ``matplotlib.axes.Axes.imshow``.

    Raises
    ------
    ValueError
        If ``atom_index`` is out of range.
    """
    from rdkit.Chem import Draw

    n = mol.GetNumAtoms()
    if atom_index < 0 or atom_index >= n:
        msg = f"atom_index {atom_index} out of range for mol with {n} atoms"
        raise ValueError(msg)
    mol = Chem.Mol(mol)
    AllChem.Compute2DCoords(mol)
    img = Draw.MolToImage(
        mol,
        size=size,
        highlightAtoms=[atom_index],
        highlightAtomColors={atom_index: (0.92, 0.20, 0.20)},
        highlightBonds=[],
    )
    return np.asarray(img)
