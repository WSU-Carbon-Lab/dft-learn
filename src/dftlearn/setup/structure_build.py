"""Build 3D coordinates from SMILES and thermally relax structures."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from dftlearn.io.xyz_structure import element_symbol_from_xyz_label


class RelaxMethod(str, Enum):
    """Supported structure relaxation modes."""

    NONE = "none"
    UFF = "uff"
    ANNEAL = "anneal"


@dataclass(frozen=True)
class BuiltStructure:
    """3D structure built from SMILES with optional relaxation metadata."""

    mol: Chem.Mol
    rows: list[tuple[str, float, float, float]]
    atom_labels: list[str]
    relax_method: RelaxMethod
    relax_steps: int
    final_energy: float | None


def mol_to_xyz_rows(
    mol: Chem.Mol,
) -> tuple[list[tuple[str, float, float, float]], list[str]]:
    """Convert an RDKit molecule conformer to StoBe-style XYZ rows.

    Parameters
    ----------
    mol
        Molecule with at least one 3D conformer.

    Returns
    -------
    rows
        ``(label, x, y, z)`` rows with labels ``C01``, ``Al01``, etc.
    atom_labels
        Human-readable labels ``C1``, ``Al1``, ... in file order.
    """
    if mol.GetNumConformers() == 0:
        msg = "mol has no 3D conformer"
        raise ValueError(msg)
    conf = mol.GetConformer()
    counts: dict[str, int] = {}
    rows: list[tuple[str, float, float, float]] = []
    labels: list[str] = []
    for idx in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(idx)
        sym = atom.GetSymbol()
        counts[sym] = counts.get(sym, 0) + 1
        stobe_label = f"{sym}{counts[sym]:02d}"
        pos = conf.GetAtomPosition(idx)
        rows.append((stobe_label, float(pos.x), float(pos.y), float(pos.z)))
        labels.append(f"{sym}{counts[sym]}")
    return rows, labels


def build_3d_from_smiles(
    smiles: str,
    *,
    seed: int = 42,
    relax: RelaxMethod = RelaxMethod.UFF,
    relax_steps: int = 500,
    anneal_cycles: int = 30,
) -> BuiltStructure:
    """Embed a SMILES string in 3D and optionally relax with UFF or annealing.

    Parameters
    ----------
    smiles
        PubChem or user-supplied SMILES.
    seed
        Random seed for ETKDG conformer generation.
    relax
        ``UFF`` minimization, ``ANNEAL`` stochastic annealing plus UFF, or ``NONE``.
    relax_steps
        Maximum UFF minimization iterations per stage.
    anneal_cycles
        Number of anneal cycles when ``relax`` is ``ANNEAL``.

    Returns
    -------
    BuiltStructure
        RDKit molecule, XYZ rows, and relaxation metadata.

    Raises
    ------
    ValueError
        If SMILES parsing or embedding fails.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        msg = f"RDKit could not parse SMILES: {smiles!r}"
        raise ValueError(msg)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()  # ty: ignore[unresolved-attribute]
    params.randomSeed = seed
    code = AllChem.EmbedMolecule(mol, params)  # ty: ignore[unresolved-attribute]
    if code != 0:
        params.useRandomCoords = True
        code = AllChem.EmbedMolecule(mol, params)  # ty: ignore[unresolved-attribute]
    if code != 0:
        msg = "RDKit could not embed a 3D conformer from SMILES"
        raise ValueError(msg)

    final_energy: float | None = None
    if relax is RelaxMethod.UFF:
        final_energy = _uff_minimize(mol, relax_steps)
    elif relax is RelaxMethod.ANNEAL:
        final_energy = _anneal_uff(mol, relax_steps, anneal_cycles, seed=seed)

    rows, labels = mol_to_xyz_rows(mol)
    _assert_nonoverlapping_rows(rows)
    return BuiltStructure(
        mol=mol,
        rows=rows,
        atom_labels=labels,
        relax_method=relax,
        relax_steps=relax_steps,
        final_energy=final_energy,
    )


def rows_to_xyz_text(
    rows: list[tuple[str, float, float, float]],
    comment: str = "",
) -> str:
    """Format XYZ rows as a text block suitable for py3Dmol."""
    lines = [str(len(rows)), comment or "generated"]
    for label, x, y, z in rows:
        sym = element_symbol_from_xyz_label(label)
        lines.append(f"{sym} {x:.6f} {y:.6f} {z:.6f}")
    return "\n".join(lines) + "\n"


def _min_pairwise_distance(
    rows: list[tuple[str, float, float, float]],
) -> tuple[float, str, str]:
    """Return the shortest interatomic distance and the two labels involved."""
    coords = np.array([[x, y, z] for _label, x, y, z in rows], dtype=np.float64)
    labels = [label for label, *_ in rows]
    min_dist = float("inf")
    left = labels[0]
    right = labels[0]
    for i in range(len(rows)):
        deltas = coords[i + 1 :] - coords[i]
        if deltas.size == 0:
            continue
        dists = np.linalg.norm(deltas, axis=1)
        j_rel = int(np.argmin(dists))
        dist = float(dists[j_rel])
        if dist < min_dist:
            min_dist = dist
            left = labels[i]
            right = labels[i + 1 + j_rel]
    return min_dist, left, right


def _assert_nonoverlapping_rows(
    rows: list[tuple[str, float, float, float]],
    *,
    min_angstrom: float = 0.5,
) -> None:
    """Raise when any pair of atoms is closer than ``min_angstrom``."""
    if len(rows) < 2:
        return
    dist, left, right = _min_pairwise_distance(rows)
    if dist < min_angstrom:
        msg = (
            f"Embedded structure has overlapping atoms {left} and {right} "
            f"({dist:.4f} A). Use a crystal CIF for metal complexes such as AlQ3 "
            "instead of SMILES embedding."
        )
        raise ValueError(msg)


def _uff_minimize(mol: Chem.Mol, max_iters: int) -> float:
    ff = AllChem.UFFGetMoleculeForceField(mol)  # ty: ignore[unresolved-attribute]
    if ff is None:
        msg = "UFF force field unavailable for this structure"
        raise ValueError(msg)
    ff.Initialize()
    ff.Minimize(maxIts=max_iters)
    return float(ff.CalcEnergy())


def _anneal_uff(
    mol: Chem.Mol,
    max_iters: int,
    cycles: int,
    *,
    seed: int,
) -> float:
    rng = np.random.default_rng(seed)
    energy = _uff_minimize(mol, max_iters)
    conf = mol.GetConformer()
    n_atoms = mol.GetNumAtoms()
    for cycle in range(cycles):
        temperature = 300.0 * (1.0 - cycle / max(cycles, 1))
        scale = 0.05 * (temperature / 300.0)
        backup = [
            (
                conf.GetAtomPosition(i).x,
                conf.GetAtomPosition(i).y,
                conf.GetAtomPosition(i).z,
            )
            for i in range(n_atoms)
        ]
        for i in range(n_atoms):
            x, y, z = backup[i]
            conf.SetAtomPosition(
                i,
                (
                    x + float(rng.normal(0.0, scale)),
                    y + float(rng.normal(0.0, scale)),
                    z + float(rng.normal(0.0, scale)),
                ),
            )
        try:
            trial = _uff_minimize(mol, max(50, max_iters // 10))
        except ValueError:
            for i, (x, y, z) in enumerate(backup):
                conf.SetAtomPosition(i, (x, y, z))
            continue
        if trial <= energy:
            energy = trial
        else:
            for i, (x, y, z) in enumerate(backup):
                conf.SetAtomPosition(i, (x, y, z))
    energy = _uff_minimize(mol, max_iters)
    return energy
