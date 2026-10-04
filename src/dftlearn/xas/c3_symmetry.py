"""C3 molecular-frame alignment and dipole-tensor folding for TP sticks.

Builds an orthonormal frame from the Al-N / Al-O coordination triangle so the
ligand plane lies in xy (C3 axis along z), then folds each transition dipole
by summing outer-product tensors over rotations of 0, 120, and 240 degrees
about z and normalizing the trace. Cartesian oscillator strengths are taken
from the folded tensor diagonals.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from dftlearn.io.xyz_structure import (
    element_symbol_from_xyz_label,
    xyz_rows_from_file,
)

Mat3 = NDArray[np.float64]
Vec3 = NDArray[np.float64]

_EPS = 1e-12


@dataclass(frozen=True)
class C3Frame:
    """Al-N/O coordination triangle frame for C3 dipole folding.

    Attributes
    ----------
    al_index
        Zero-based XYZ row index of the aluminum atom.
    ligand_pairs
        Three ``(n_index, o_index)`` bite pairs in ligand order.
    ligand_centroids
        Shape ``(3, 3)`` Cartesian centroids of each N-O bite (Angstrom).
    rotation
        Shape ``(3, 3)`` matrix ``R`` such that ``mu_mol = R @ mu_stobe``.
    """

    al_index: int
    ligand_pairs: tuple[tuple[int, int], ...]
    ligand_centroids: Mat3
    rotation: Mat3

    def to_json_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable audit dictionary for ``c3_frame.json``."""
        return {
            "al_index": self.al_index,
            "ligand_pairs": [list(pair) for pair in self.ligand_pairs],
            "ligand_centroids": self.ligand_centroids.tolist(),
            "rotation": self.rotation.tolist(),
        }


def rotation_matrix_z(angle_deg: float) -> Mat3:
    """Build a right-handed rotation matrix about +z by ``angle_deg`` degrees."""
    theta = np.deg2rad(angle_deg)
    c = float(np.cos(theta))
    s = float(np.sin(theta))
    return np.array(
        [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]],
        dtype=np.float64,
    )


def _positions_from_rows(
    rows: list[tuple[str, float, float, float]],
) -> Mat3:
    return np.array([[r[1], r[2], r[3]] for r in rows], dtype=np.float64)


def _indices_by_element(
    rows: list[tuple[str, float, float, float]],
    symbol: str,
) -> list[int]:
    key = symbol.strip().capitalize()
    return [
        i for i, row in enumerate(rows) if element_symbol_from_xyz_label(row[0]) == key
    ]


def pair_n_o_ligand_bites(
    rows: list[tuple[str, float, float, float]],
) -> tuple[int, tuple[tuple[int, int], ...]]:
    """Pair each nitrogen with its nearest oxygen as an Al ligand bite.

    Parameters
    ----------
    rows
        XYZ rows ``(label, x, y, z)``.

    Returns
    -------
    al_index : int
        Index of the single Al atom.
    pairs : tuple of (n_index, o_index)
        Three unique N-O pairs ordered by ascending nitrogen index.

    Raises
    ------
    ValueError
        If the structure does not contain exactly one Al and three N / three O,
        or if pairing cannot assign unique oxygens.
    """
    al = _indices_by_element(rows, "Al")
    nitrogens = _indices_by_element(rows, "N")
    oxygens = _indices_by_element(rows, "O")
    if len(al) != 1:
        msg = f"C3 frame requires exactly one Al atom, found {len(al)}"
        raise ValueError(msg)
    if len(nitrogens) != 3 or len(oxygens) != 3:
        msg = (
            "C3 frame requires three N and three O atoms, "
            f"found N={len(nitrogens)} O={len(oxygens)}"
        )
        raise ValueError(msg)
    pos = _positions_from_rows(rows)
    used_o: set[int] = set()
    pairs: list[tuple[int, int]] = []
    for n_idx in sorted(nitrogens):
        candidates = [
            (float(np.linalg.norm(pos[n_idx] - pos[o_idx])), o_idx)
            for o_idx in oxygens
            if o_idx not in used_o
        ]
        if not candidates:
            msg = "Could not pair each N with a unique O for ligand bites"
            raise ValueError(msg)
        candidates.sort(key=lambda item: (item[0], item[1]))
        o_idx = candidates[0][1]
        used_o.add(o_idx)
        pairs.append((n_idx, o_idx))
    return al[0], tuple(pairs)


def build_al_n_o_c3_frame(
    rows: list[tuple[str, float, float, float]],
) -> C3Frame:
    """Build the molecular C3 frame from the Al-N/O coordination triangle.

    Ligand centroids form a triangle whose plane becomes xy. The triangle
    normal is molecular z. Ligand-0's in-plane Al→centroid direction defines
    +x after projection.

    Parameters
    ----------
    rows
        XYZ geometry rows.

    Returns
    -------
    C3Frame
        Frame with ``rotation`` mapping StoBe vectors into the molecular frame.

    Raises
    ------
    ValueError
        If ligand bites cannot be formed or the triangle is degenerate.
    """
    al_index, pairs = pair_n_o_ligand_bites(rows)
    pos = _positions_from_rows(rows)
    al = pos[al_index]
    centroids = np.array(
        [(pos[n] + pos[o]) * 0.5 for n, o in pairs],
        dtype=np.float64,
    )
    v01 = centroids[1] - centroids[0]
    v02 = centroids[2] - centroids[0]
    normal = np.cross(v01, v02)
    n_norm = float(np.linalg.norm(normal))
    if n_norm < _EPS:
        msg = "Al-N/O ligand triangle is degenerate; cannot define a C3 axis"
        raise ValueError(msg)
    z_hat = normal / n_norm
    lig0 = centroids[0] - al
    x_proj = lig0 - np.dot(lig0, z_hat) * z_hat
    x_norm = float(np.linalg.norm(x_proj))
    if x_norm < _EPS:
        fallback = np.array([1.0, 0.0, 0.0], dtype=np.float64)
        if abs(float(np.dot(fallback, z_hat))) > 0.9:
            fallback = np.array([0.0, 1.0, 0.0], dtype=np.float64)
        x_proj = fallback - np.dot(fallback, z_hat) * z_hat
        x_norm = float(np.linalg.norm(x_proj))
    x_hat = x_proj / x_norm
    y_hat = np.cross(z_hat, x_hat)
    y_hat /= float(np.linalg.norm(y_hat))
    rotation = np.vstack((x_hat, y_hat, z_hat)).astype(np.float64)
    return C3Frame(
        al_index=al_index,
        ligand_pairs=pairs,
        ligand_centroids=centroids,
        rotation=rotation,
    )


def build_c3_frame_from_xyz(path: Path) -> C3Frame:
    """Load XYZ geometry and build :class:`C3Frame`."""
    return build_al_n_o_c3_frame(xyz_rows_from_file(Path(path)))


def fold_dipole_c3_tensor(mu_stobe: Vec3, rotation: Mat3) -> Mat3:
    """Fold one StoBe dipole under C3 about molecular z after frame rotation.

    Computes ``mu' = R @ mu``, then
    ``T = sum_k Rz(k*120) mu' mu'^T Rz^T``, then scales so
    ``Tr(T) = ||mu||^2``.

    Parameters
    ----------
    mu_stobe
        Dipole components in StoBe Cartesian axes, shape ``(3,)``.
    rotation
        Molecular-frame matrix ``R`` from :class:`C3Frame`.

    Returns
    -------
    numpy.ndarray
        Shape ``(3, 3)`` symmetrized absorption tensor.
    """
    mu = np.asarray(mu_stobe, dtype=np.float64).reshape(3)
    r = np.asarray(rotation, dtype=np.float64)
    mu_mol = r @ mu
    tot = np.zeros((3, 3), dtype=np.float64)
    for k in range(3):
        rz = rotation_matrix_z(120.0 * k)
        mu_k = rz @ mu_mol
        tot += np.outer(mu_k, mu_k)
    mu2 = float(np.dot(mu, mu))
    tr = float(np.trace(tot))
    if tr > _EPS and mu2 > _EPS:
        tot *= mu2 / tr
    elif tr <= _EPS:
        tot[:] = 0.0
    return tot


def fold_dipoles_c3(dipoles_stobe: Mat3, rotation: Mat3) -> NDArray[np.float64]:
    """Fold a batch of StoBe dipoles under C3.

    Parameters
    ----------
    dipoles_stobe
        Shape ``(m, 3)`` dipole matrix elements.
    rotation
        Molecular-frame ``R``.

    Returns
    -------
    numpy.ndarray
        Shape ``(m, 3, 3)`` folded tensors.

    Raises
    ------
    ValueError
        If ``dipoles_stobe`` is not ``(m, 3)``.
    """
    mu = np.asarray(dipoles_stobe, dtype=np.float64)
    if mu.ndim != 2 or mu.shape[1] != 3:
        msg = f"dipoles_stobe must have shape (m, 3), got {mu.shape}"
        raise ValueError(msg)
    out = np.empty((mu.shape[0], 3, 3), dtype=np.float64)
    for i in range(mu.shape[0]):
        out[i] = fold_dipole_c3_tensor(mu[i], rotation)
    return out


def oscillator_strengths_from_tensors(
    energy_ha: np.ndarray,
    tensors: NDArray[np.float64],
) -> Mat3:
    r"""Convert folded tensor diagonals to StoBe Cartesian oscillator strengths.

    Uses :math:`f_{ii} = (2/3) E_\mathrm{Ha}\, T_{ii}`.

    Parameters
    ----------
    energy_ha
        Transition energies in Hartree, shape ``(m,)``.
    tensors
        Folded tensors, shape ``(m, 3, 3)``.

    Returns
    -------
    numpy.ndarray
        Shape ``(m, 3)`` with columns ``(f_xx, f_yy, f_zz)``.

    Raises
    ------
    ValueError
        If shapes are incompatible.
    """
    energy = np.asarray(energy_ha, dtype=np.float64)
    tens = np.asarray(tensors, dtype=np.float64)
    if energy.ndim != 1 or tens.shape != (energy.size, 3, 3):
        msg = "energy_ha must be (m,) and tensors (m, 3, 3)"
        raise ValueError(msg)
    diag = np.diagonal(tens, axis1=1, axis2=2)
    scale = (2.0 / 3.0) * energy[:, np.newaxis]
    return np.asarray(scale * diag, dtype=np.float64)


def write_c3_frame_json(frame: C3Frame, path: Path) -> Path:
    """Write ``c3_frame.json`` audit metadata for a packaged postprocess run."""
    path = Path(path)
    path.write_text(json.dumps(frame.to_json_dict(), indent=2) + "\n", encoding="utf-8")
    return path
