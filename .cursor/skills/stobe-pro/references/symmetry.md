# Symmetry

> Reference for: StoBe Pro
> Load when: Using symmetry groups, understanding point groups, symmetry operations, reducing computational cost

---

## Overview

Symmetry groups allow StoBe to exploit molecular symmetry to reduce computational cost and simplify orbital analysis. The symmetry library (`symbasis.new`) contains definitions for all supported point groups.

---

## Symmetry Library (symbasis.new)

### File Structure

The symmetry library defines point groups with:
- Group name (e.g., `C1`, `C2`, `D2H`)
- Number of symmetry operations
- Symmetry operation matrices
- Irreducible representations
- Character tables

### Format Example

```
C1
    1    1
         0         0         0         E
    1
A     1
1
STOP
```

**Structure**:
- First line: Group name
- Second line: Number of operations, number of classes
- Operation definitions: Rotation angles, operation type
- Irreducible representations: Labels and dimensions
- Character table: Representation matrices
- `STOP`: End of group definition

### Access

Symmetry library is referenced via `fort.4` in run scripts:
```bash
ln -s ${STOBE}Basis/symbasis.new fort.4
```

---

## Supported Symmetry Groups

### Low Symmetry Groups

- **C1**: No symmetry (most common for large molecules)
- **Ci**: Inversion center only
- **Cs**: Reflection plane only

### Cyclic Groups

- **C2**: Two-fold rotation axis
- **C3**: Three-fold rotation axis
- **C4**: Four-fold rotation axis
- **C5**: Five-fold rotation axis
- **C6**: Six-fold rotation axis

### Dihedral Groups

- **D2**: Two perpendicular two-fold axes
- **D3**: Three-fold axis with perpendicular two-fold axes
- **D4**: Four-fold axis with perpendicular two-fold axes
- **D5**: Five-fold axis with perpendicular two-fold axes
- **D6**: Six-fold axis with perpendicular two-fold axes

### Groups with Vertical Planes

- **C2v**: C2 axis with vertical mirror planes
- **C3v**: C3 axis with vertical mirror planes
- **C4v**: C4 axis with vertical mirror planes
- **C5v**: C5 axis with vertical mirror planes
- **C6v**: C6 axis with vertical mirror planes

### Groups with Horizontal Planes

- **C2h**: C2 axis with horizontal mirror plane
- **C3h**: C3 axis with horizontal mirror plane
- **C4h**: C4 axis with horizontal mirror plane
- **C5h**: C5 axis with horizontal mirror plane
- **C6h**: C6 axis with horizontal mirror plane

- **D2h**: Full octahedral symmetry
- **D3h**: D3 with horizontal plane
- **D4h**: D4 with horizontal plane
- **D5h**: D5 with horizontal plane
- **D6h**: D6 with horizontal plane

### Cubic Groups

- **T**: Tetrahedral rotation group
- **Th**: T with inversion
- **Td**: Full tetrahedral symmetry
- **O**: Octahedral rotation group
- **Oh**: Full octahedral symmetry

### Special Groups

- **S4**: Four-fold improper rotation
- **S6**: Six-fold improper rotation
- **D2d**: D2 with diagonal planes
- **D3d**: D3 with diagonal planes
- **D4d**: D4 with diagonal planes
- **D5d**: D5 with diagonal planes
- **D6d**: D6 with diagonal planes

---

## Using Symmetry

### Specifying Symmetry

In input file:
```stobe
SYMMETRY C1
```

or

```stobe
SYMMETRY D2H
```

### Disabling Symmetry

```stobe
NOSYMM
```

or

```stobe
NOSYMMETRY
```

### Symmetry Requirements

1. **Molecular structure must match symmetry**: Atoms must be positioned according to symmetry operations
2. **Symmetry group must exist**: Check `symbasis.new` for valid groups
3. **Geometry input**: Can use reduced input (only symmetry-inequivalent atoms) or full input

---

## Symmetry Benefits

### Computational Cost Reduction

- Reduces number of integrals to compute
- Blocks Fock matrix by symmetry
- Reduces memory requirements
- Faster SCF convergence

### Orbital Analysis

- Orbitals belong to irreducible representations
- Easier to interpret orbital symmetries
- Simplifies transition analysis
- Clearer visualization

---

## Symmetry Considerations

### When to Use C1 (No Symmetry)

- Large molecules without clear symmetry
- Disordered systems
- Most common choice for complex molecules
- Simplest to set up

### When to Use Higher Symmetry

- Small, symmetric molecules
- Well-defined point groups
- When computational cost is limiting
- When symmetry analysis is needed

### Symmetry Validation

Before using symmetry:
1. Verify molecular structure matches symmetry
2. Check atom positions are symmetry-equivalent
3. Ensure symmetry group exists in `symbasis.new`
4. Test with `NOSYMM` if unsure

---

## Symmetry Operations

### Rotation Operations

- `C2`, `C3`, `C4`, etc.: n-fold rotations
- Defined by rotation axis and angle
- Examples: `C2Z` (rotation about Z axis)

### Reflection Operations

- `SIGX`, `SIGY`, `SIGZ`: Reflection planes
- `SIGH`: Horizontal reflection plane
- `SIGV`: Vertical reflection plane

### Improper Rotations

- `S4`, `S6`: Improper rotations (rotation + reflection)
- Used in S4, S6 groups

### Inversion

- `I`: Inversion center
- Used in Ci, Th, Oh groups

---

## Irreducible Representations

### Notation

Common labels:
- **A, B**: One-dimensional representations
- **E**: Two-dimensional representations
- **T**: Three-dimensional representations
- **g, u**: Even/odd under inversion (for centrosymmetric groups)
- **', "**: Even/odd under horizontal reflection

### Examples

**C2v**:
- A1, A2, B1, B2

**D2h**:
- Ag, Au, B1g, B1u, B2g, B2u, B3g, B3u

**Td**:
- A1, A2, E, T1, T2

---

## Symmetry in Calculations

### Ground State

- Orbitals transform as irreducible representations
- Fock matrix is block-diagonal
- Each block corresponds to one representation

### Excited States

- Excitations between symmetry-adapted orbitals
- Selection rules based on symmetry
- Transition moments depend on symmetry

### X-ray Calculations

- Core-hole breaks symmetry
- Often use C1 for core-hole calculations
- Can use symmetry for non-excited atoms

---

## Symmetry Reduction

### Reduced Input

For symmetric molecules, can specify only symmetry-inequivalent atoms:

**Example for C2v**:
```stobe
SYMMETRY C2V
CARTESIAN ANGSTROMS
N1    0.0000000     0.0000000     1.4302597    7.0
C1    0.0000000     0.0000000    -1.3920631    6.0
H1    0.0000000     0.0000000    -2.4870891    1.0
C2    1.1995239     0.0000000    -0.6738739    6.0
END
```

StoBe generates all symmetry-equivalent atoms automatically.

### Full Input

Can also specify all atoms explicitly:
```stobe
SYMMETRY C2V
CARTESIAN ANGSTROMS
[all atoms specified]
END
```

---

## Finding Symmetry

### Automatic Detection

Use StoBe utility `symfind` or `symgen`:
- Analyzes molecular structure
- Determines point group
- Generates symmetry file

### Manual Determination

1. Identify rotation axes
2. Identify reflection planes
3. Identify inversion center (if present)
4. Match to point group table
5. Verify with `symbasis.new`

---

## Common Patterns

### Pattern 1: Large Molecules

Most large molecules use:
```stobe
SYMMETRY C1
```

No symmetry exploitation, simplest setup.

### Pattern 2: Small Symmetric Molecules

Well-defined symmetry:
```stobe
SYMMETRY C2V
```

or

```stobe
SYMMETRY D2H
```

### Pattern 3: Core-Hole Calculations

Often break symmetry:
```stobe
SYMMETRY C1
```

Even if molecule has symmetry, core-hole breaks it.

---

## Validation

When using symmetry:

- [ ] Symmetry group exists in `symbasis.new`
- [ ] Molecular structure matches symmetry
- [ ] Atom positions are correct
- [ ] Symmetry operations are valid
- [ ] Reduced input (if used) is correct

---

## Troubleshooting

### Invalid Symmetry Group

**Error**: Symmetry group not found

**Fix**: Check `symbasis.new` for valid groups, or use `NOSYMM`

### Structure Doesn't Match Symmetry

**Error**: Symmetry operations don't match atoms

**Fix**: Verify atom positions, or use `NOSYMM`

### Symmetry Breaking

**Error**: Expected symmetry not found in output

**Fix**: Check geometry, verify symmetry operations

---

## Related References

- **Input File Format**: See `input-file-format.md` for symmetry specification
- **Calculation Workflows**: See `calculation-workflows.md` for symmetry in workflows
- **Examples**: See `examples-catalog.md` for symmetry examples
