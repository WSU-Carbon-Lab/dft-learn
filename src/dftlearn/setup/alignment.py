"""3D alignment: rotation matrices and principal-axis reorientation."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from dftlearn.io.xyz_structure import element_symbol_from_xyz_label

Mat3 = NDArray[np.float64]
Vec3 = NDArray[np.float64]


def identity_matrix() -> Mat3:
    """Return a 3x3 identity rotation matrix."""
    return np.eye(3, dtype=np.float64)


def euler_matrix(rx_deg: float, ry_deg: float, rz_deg: float) -> Mat3:
    """Build a 3x3 rotation matrix from XYZ Euler angles in degrees."""
    rx, ry, rz = np.deg2rad([rx_deg, ry_deg, rz_deg])
    cx, sx = np.cos(rx), np.sin(rx)
    cy, sy = np.cos(ry), np.sin(ry)
    cz, sz = np.cos(rz), np.sin(rz)
    rx_m = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]], dtype=np.float64)
    ry_m = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]], dtype=np.float64)
    rz_m = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]], dtype=np.float64)
    return rz_m @ ry_m @ rx_m


def rotate_positions(
    positions: NDArray[np.float64],
    matrix: Mat3,
) -> NDArray[np.float64]:
    """Apply a 3x3 rotation matrix to ``(n, 3)`` Cartesian coordinates."""
    return positions @ matrix.T


def bond_vector(
    rows: list[tuple[str, float, float, float]],
    atom_i: int,
    atom_j: int,
) -> Vec3:
    """Return the displacement vector from ``atom_i`` to ``atom_j``."""
    pi = np.array(rows[atom_i][1:4], dtype=np.float64)
    pj = np.array(rows[atom_j][1:4], dtype=np.float64)
    return pj - pi


def axis_unit_vector(axis: str) -> Vec3:
    """Return a unit vector for ``axis`` in ``{'x', 'y', 'z'}``."""
    key = axis.strip().lower()
    if key == "x":
        return np.array([1.0, 0.0, 0.0], dtype=np.float64)
    if key == "y":
        return np.array([0.0, 1.0, 0.0], dtype=np.float64)
    if key == "z":
        return np.array([0.0, 0.0, 1.0], dtype=np.float64)
    msg = f"axis must be x, y, or z, got {axis!r}"
    raise ValueError(msg)


def align_bond_to_axis(
    rows: list[tuple[str, float, float, float]],
    atom_i: int,
    atom_j: int,
    axis: str,
) -> Mat3:
    """Build a rotation matrix mapping the bond ``atom_i``->``atom_j`` onto ``axis``."""
    vec = bond_vector(rows, atom_i, atom_j)
    return rotation_align_vector(vec, axis_unit_vector(axis))


def align_frame_from_bond(
    rows: list[tuple[str, float, float, float]],
    bond_i: int,
    bond_j: int,
    plane_atom: int,
) -> Mat3:
    """Align a bond to +X and rotate about X until ``plane_atom`` lies in the XY plane.

    The returned matrix maps the bond vector onto +X, then applies an X rotation so
    the vector from ``bond_j`` to ``plane_atom`` has zero Z component.
    """
    bond_to_x = align_bond_to_axis(rows, bond_i, bond_j, "x")
    pos = np.array([[r[1], r[2], r[3]] for r in rows], dtype=np.float64)
    anchor = rotate_positions(pos[bond_j : bond_j + 1], bond_to_x)[0]
    plane = rotate_positions(pos[plane_atom : plane_atom + 1], bond_to_x)[0]
    reference = plane - anchor
    if np.linalg.norm(reference) < 1e-12:
        return bond_to_x
    angle_deg = -float(np.degrees(np.arctan2(reference[2], reference[1])))
    return euler_matrix(angle_deg, 0.0, 0.0) @ bond_to_x


def rotation_align_vector(source: Vec3, target: Vec3) -> Mat3:
    """Compute the rotation matrix that maps ``source`` onto ``target``.

    Both vectors are normalized before computing the rotation. Returns identity
    when the vectors are parallel or anti-parallel within numerical tolerance.
    """
    src = np.asarray(source, dtype=np.float64)
    tgt = np.asarray(target, dtype=np.float64)
    src_norm = np.linalg.norm(src)
    tgt_norm = np.linalg.norm(tgt)
    if src_norm < 1e-12 or tgt_norm < 1e-12:
        return identity_matrix()
    src_u = src / src_norm
    tgt_u = tgt / tgt_norm
    cross = np.cross(src_u, tgt_u)
    dot = float(np.clip(np.dot(src_u, tgt_u), -1.0, 1.0))
    if np.linalg.norm(cross) < 1e-12:
        if dot > 0:
            return identity_matrix()
        axis = np.array([1.0, 0.0, 0.0], dtype=np.float64)
        if abs(src_u[0]) > 0.9:
            axis = np.array([0.0, 1.0, 0.0], dtype=np.float64)
        cross = np.cross(src_u, axis)
        cross /= np.linalg.norm(cross)
        return euler_matrix(180.0, 0.0, 0.0)
    skew = np.array(
        [
            [0.0, -cross[2], cross[1]],
            [cross[2], 0.0, -cross[0]],
            [-cross[1], cross[0], 0.0],
        ],
        dtype=np.float64,
    )
    return identity_matrix() + skew + skew @ skew * (1.0 / (1.0 + dot))


def principal_inertia_axis(
    positions: NDArray[np.float64],
    masses: NDArray[np.float64] | None = None,
) -> tuple[Vec3, Vec3, Vec3]:
    """Return principal inertia axes sorted by ascending eigenvalue.

    Parameters
    ----------
    positions
        ``(n, 3)`` Cartesian coordinates in angstroms.
    masses
        Optional per-atom masses; defaults to uniform unit mass.

    Returns
    -------
    axes
        Three orthonormal axis vectors (smallest, middle, largest inertia).
    """
    pos = np.asarray(positions, dtype=np.float64)
    if pos.ndim != 2 or pos.shape[1] != 3:
        msg = f"positions must have shape (n, 3), got {pos.shape}"
        raise ValueError(msg)
    if masses is None:
        mass = np.ones(pos.shape[0], dtype=np.float64)
    else:
        mass = np.asarray(masses, dtype=np.float64)
    com = np.average(pos, axis=0, weights=mass)
    centered = pos - com
    inertia = np.zeros((3, 3), dtype=np.float64)
    for vec, m in zip(centered, mass, strict=True):
        r2 = float(np.dot(vec, vec))
        inertia += m * (r2 * np.eye(3) - np.outer(vec, vec))
    evals, evecs = np.linalg.eigh(inertia)
    order = np.argsort(evals)
    axes = tuple(
        np.asarray(evecs[:, order[i]], dtype=np.float64) for i in range(3)
    )
    return axes[0], axes[1], axes[2]


def metal_ligand_axis(
    rows: list[tuple[str, float, float, float]],
    *,
    metal_symbols: tuple[str, ...] = ("Al", "Zn", "Cu", "Fe", "Ni", "Co", "Mg"),
    ligand_symbols: tuple[str, ...] = ("N", "O", "S", "P"),
) -> Vec3:
    """Estimate a metal-to-ligand alignment vector from XYZ rows.

    Uses the centroid of ligand atoms minus the centroid of metal atoms. Falls
    back to the principal inertia axis when no metal is present.
    """
    pos = np.array([[r[1], r[2], r[3]] for r in rows], dtype=np.float64)
    metals: list[NDArray[np.float64]] = []
    ligands: list[NDArray[np.float64]] = []
    for row, p in zip(rows, pos, strict=True):
        sym = element_symbol_from_xyz_label(row[0])
        if sym in metal_symbols:
            metals.append(p)
        elif sym in ligand_symbols:
            ligands.append(p)
    if metals and ligands:
        return np.mean(ligands, axis=0) - np.mean(metals, axis=0)
    _small, _mid, largest = principal_inertia_axis(pos)
    return largest


def apply_rotation_to_rows(
    rows: list[tuple[str, float, float, float]],
    matrix: Mat3,
) -> list[tuple[str, float, float, float]]:
    """Rotate XYZ coordinates while preserving atom labels."""
    pos = np.array([[r[1], r[2], r[3]] for r in rows], dtype=np.float64)
    rotated = rotate_positions(pos, matrix)
    return [
        (label, float(x), float(y), float(z))
        for (label, _x, _y, _z), (x, y, z) in zip(rows, rotated, strict=True)
    ]


def projection_indices(view: str) -> tuple[int, int]:
    """Map view name to coordinate column indices for 2D projection."""
    key = view.strip().lower()
    if key == "xz":
        return 0, 2
    if key == "yz":
        return 1, 2
    return 0, 1
