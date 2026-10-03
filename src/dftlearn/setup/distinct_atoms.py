"""Chemically distinct atom grouping from 3D coordinates and connectivity."""

from __future__ import annotations

from collections import Counter, defaultdict
from typing import TYPE_CHECKING

from rdkit import Chem

from dftlearn.io.xyz_structure import element_symbol_from_xyz_label
from dftlearn.setup.types import DistinctAtomGroup
from dftlearn.visualization.xyz_wireframe import infer_bonds_from_xyz_rows

if TYPE_CHECKING:
    from dftlearn.setup.types import CifAtomSite


def _stobe_xyz_label(element: str, ordinal: int) -> str:
    return f"{element}{ordinal:02d}"


def sites_to_xyz_rows(
    sites: tuple[CifAtomSite, ...],
) -> tuple[list[tuple[str, float, float, float]], list[str]]:
    """Convert CIF sites to StoBe-style XYZ rows preserving input order.

    Parameters
    ----------
    sites
        Cartesian atom sites from :func:`load_cif_structure`.

    Returns
    -------
    rows
        ``(label, x, y, z)`` with labels ``C01``, ``Al01``, etc.
    cif_labels
        Original CIF ``_atom_site_label`` per row index.
    """
    counts: Counter[str] = Counter()
    rows: list[tuple[str, float, float, float]] = []
    cif_labels: list[str] = []
    for site in sites:
        counts[site.element] += 1
        label = _stobe_xyz_label(site.element, counts[site.element])
        rows.append((label, site.x, site.y, site.z))
        cif_labels.append(site.label)
    return rows, cif_labels


def mol_from_xyz_rows(
    rows: list[tuple[str, float, float, float]],
    *,
    sanitize: bool = True,
) -> Chem.Mol:
    """Build an RDKit molecule with inferred bonds from XYZ rows.

    When ``sanitize`` is True, runs :func:`Chem.SanitizeMol` when valence permits.
    On valence failure, returns an unsanitized molecule with greedily pruned bonds.
    """
    rw = Chem.RWMol()
    conf = Chem.Conformer()
    for label, x, y, z in rows:
        sym = element_symbol_from_xyz_label(label)
        idx = rw.AddAtom(Chem.Atom(sym))
        conf.SetAtomPosition(idx, (x, y, z))
    rw.AddConformer(conf, assignId=True)
    candidates = infer_bonds_from_xyz_rows(rows)
    bonds = _greedy_valence_bonds(rows, candidates)
    for i, j in bonds:
        if rw.GetBondBetweenAtoms(i, j) is None:
            rw.AddBond(i, j, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    if not sanitize:
        return mol
    try:
        Chem.SanitizeMol(mol)
    except Chem.AtomValenceException:
        return mol
    return mol


_MAX_VALENCE: dict[str, int] = {
    "H": 1,
    "C": 4,
    "N": 4,
    "O": 2,
    "S": 2,
    "P": 3,
    "F": 1,
    "Cl": 1,
    "Br": 1,
    "I": 1,
    "Al": 6,
    "Zn": 4,
}


def _greedy_valence_bonds(
    rows: list[tuple[str, float, float, float]],
    candidates: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    """Keep shortest candidate bonds that respect per-element valence limits."""
    if not candidates:
        return []
    symbols = [element_symbol_from_xyz_label(label) for label, *_ in rows]
    max_valence = [_MAX_VALENCE.get(sym, 4) for sym in symbols]
    pos = [(r[1], r[2], r[3]) for r in rows]

    def bond_length(ij: tuple[int, int]) -> float:
        i, j = ij
        dx = pos[i][0] - pos[j][0]
        dy = pos[i][1] - pos[j][1]
        dz = pos[i][2] - pos[j][2]
        return (dx * dx + dy * dy + dz * dz) ** 0.5

    valence_used = [0] * len(rows)
    accepted: list[tuple[int, int]] = []
    for i, j in sorted(candidates, key=bond_length):
        if valence_used[i] >= max_valence[i] or valence_used[j] >= max_valence[j]:
            continue
        accepted.append((i, j))
        valence_used[i] += 1
        valence_used[j] += 1
    return accepted


def distinct_atom_groups(
    mol: Chem.Mol,
    *,
    element: str,
    rows: list[tuple[str, float, float, float]],
    cif_labels: list[str],
) -> list[DistinctAtomGroup]:
    """Group atoms of ``element`` by RDKit canonical rank (symmetry class).

    Parameters
    ----------
    mol
        Molecule with 3D conformer and bonds aligned with ``rows``.
    element
        Element symbol filter, for example ``C`` for carbon K-edge sites.
    rows
        XYZ rows parallel to mol atom indices.
    cif_labels
        CIF labels parallel to mol atom indices.

    Returns
    -------
    list[DistinctAtomGroup]
        One group per distinct rank, sorted by descending class size then rank.
    """
    target = element.strip().capitalize()
    try:
        ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    except RuntimeError:
        ranks = list(range(mol.GetNumAtoms()))
    buckets: dict[int, list[int]] = defaultdict(list)
    for idx, atom in enumerate(mol.GetAtoms()):
        if atom.GetSymbol() == target:
            buckets[ranks[idx]].append(idx)

    groups: list[DistinctAtomGroup] = []
    for rank, indices in buckets.items():
        rep = min(indices)
        groups.append(
            DistinctAtomGroup(
                rank=rank,
                element=target,
                atom_indices=tuple(indices),
                representative_index=rep,
                cif_label=cif_labels[rep],
                xyz_label=rows[rep][0],
            )
        )
    groups.sort(key=lambda g: (-len(g.atom_indices), g.rank))
    return groups


_METAL_ELEMENTS = frozenset(
    {
        "Al",
        "Zn",
        "Fe",
        "Ni",
        "Cu",
        "Co",
        "Mn",
        "Ti",
        "V",
        "Cr",
        "Mo",
        "W",
        "Pd",
        "Pt",
        "Au",
        "Ag",
        "Ru",
        "Rh",
        "Ir",
        "Os",
        "Re",
        "Nb",
        "Ta",
        "Sc",
        "Y",
        "La",
        "Ga",
        "In",
        "Sn",
        "Pb",
    }
)


def _connected_components(
    n_atoms: int,
    bonds: list[tuple[int, int]],
) -> list[list[int]]:
    adjacency: list[list[int]] = [[] for _ in range(n_atoms)]
    for i, j in bonds:
        adjacency[i].append(j)
        adjacency[j].append(i)
    seen = [False] * n_atoms
    components: list[list[int]] = []
    for start in range(n_atoms):
        if seen[start]:
            continue
        stack = [start]
        seen[start] = True
        members: list[int] = []
        while stack:
            node = stack.pop()
            members.append(node)
            for neighbor in adjacency[node]:
                if not seen[neighbor]:
                    seen[neighbor] = True
                    stack.append(neighbor)
        components.append(sorted(members))
    return components


def ligand_carbon_indices(
    rows: list[tuple[str, float, float, float]],
    *,
    ligand_index: int = 0,
) -> tuple[int, ...]:
    """Return carbon atom indices for one metal-bound organic ligand.

    Drops metal atoms, finds covalently connected fragments, and keeps fragments
    that contain carbon plus nitrogen or oxygen (quinolate-style wings). Ligands
    are ordered by the smallest carbon index so the first fragment matches the
    earliest XYZ/CIF carbons.

    Parameters
    ----------
    rows
        StoBe-style ``(label, x, y, z)`` rows.
    ligand_index
        Zero-based ligand to keep. ``0`` is the first wing in coordinate order.

    Returns
    -------
    tuple[int, ...]
        Carbon indices of the selected ligand, sorted.

    Raises
    ------
    ValueError
        If no ligand fragments are found or ``ligand_index`` is out of range.
    """
    symbols = [element_symbol_from_xyz_label(label) for label, *_ in rows]
    bonds = infer_bonds_from_xyz_rows(rows)
    metal = {idx for idx, sym in enumerate(symbols) if sym in _METAL_ELEMENTS}
    organic_bonds = [(i, j) for i, j in bonds if i not in metal and j not in metal]
    fragments = _connected_components(len(rows), organic_bonds)
    ligands: list[tuple[int, ...]] = []
    for fragment in fragments:
        elems = {symbols[idx] for idx in fragment}
        carbons = tuple(idx for idx in fragment if symbols[idx] == "C")
        if carbons and (("N" in elems) or ("O" in elems)):
            ligands.append(carbons)
    ligands.sort(key=lambda carbons: (carbons[0], len(carbons)))
    if not ligands:
        msg = "No metal-bound organic ligands with carbon sites were found"
        raise ValueError(msg)
    if ligand_index < 0 or ligand_index >= len(ligands):
        msg = f"ligand_index {ligand_index} is out of range (n_ligands={len(ligands)})"
        raise ValueError(msg)
    return ligands[ligand_index]


def select_principal_molecule_sites(
    sites: tuple[CifAtomSite, ...],
) -> tuple[CifAtomSite, ...]:
    """Keep the covalently bonded molecule, dropping residual solvent fragments.

    Builds a heavy-atom connectivity graph from covalent distances, prefers the
    component that contains a metal, and reattaches hydrogens bonded to that
    component. Isolated synthesis leftovers such as acetic acid are discarded.

    Parameters
    ----------
    sites
        Cartesian CIF sites, including any co-crystallized fragments.

    Returns
    -------
    tuple[CifAtomSite, ...]
        Sites of the principal molecule in original relative order.
    """
    if len(sites) <= 1:
        return sites
    rows, _labels = sites_to_xyz_rows(sites)
    bonds = infer_bonds_from_xyz_rows(rows)
    heavy_bonds = [
        (i, j) for i, j in bonds if sites[i].element != "H" and sites[j].element != "H"
    ]
    components = _connected_components(len(sites), heavy_bonds)
    heavy_components = [
        [idx for idx in comp if sites[idx].element != "H"] for comp in components
    ]
    heavy_components = [comp for comp in heavy_components if comp]
    if not heavy_components:
        return sites

    def score(comp: list[int]) -> tuple[int, int]:
        has_metal = any(sites[idx].element in _METAL_ELEMENTS for idx in comp)
        return (int(has_metal), len(comp))

    chosen_heavy = set(max(heavy_components, key=score))
    hydrogen_links = [
        (i, j)
        for i, j in bonds
        if (sites[i].element == "H") != (sites[j].element == "H")
    ]
    keep: set[int] = set(chosen_heavy)
    for i, j in hydrogen_links:
        if i in chosen_heavy and sites[j].element == "H":
            keep.add(j)
        elif j in chosen_heavy and sites[i].element == "H":
            keep.add(i)
    return tuple(site for idx, site in enumerate(sites) if idx in keep)


def distinct_groups_for_element(
    sites: tuple[CifAtomSite, ...],
    *,
    element: str,
) -> tuple[list[tuple[str, float, float, float]], list[str], list[DistinctAtomGroup]]:
    """Identify symmetry-distinct sites of ``element`` in a CIF structure.

    Parameters
    ----------
    sites
        Cartesian sites from a CIF block.
    element
        Element symbol for core-edge grouping.

    Returns
    -------
    rows
        StoBe-style XYZ rows in CIF order.
    cif_labels
        CIF labels parallel to rows.
    groups
        Distinct equivalence classes for ``element``.
    """
    rows, cif_labels = sites_to_xyz_rows(sites)
    mol = mol_from_xyz_rows(rows)
    groups = distinct_atom_groups(
        mol,
        element=element,
        rows=rows,
        cif_labels=cif_labels,
    )
    return rows, cif_labels, groups
