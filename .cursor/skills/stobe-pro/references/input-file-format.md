# Input File Format

> Reference for: StoBe Pro
> Load when: Parsing or generating .inp files, understanding input structure, validating input files

---

## Overview

StoBe input files (.inp) follow a structured format with specific sections. Understanding this format is essential for generating valid input files and parsing existing ones.

---

## File Structure

StoBe input files consist of six main sections:

1. Header (TITLE, SYMMETRY, CARTESIAN)
2. Geometry (atom coordinates)
3. Calculation Parameters (RUNTYPE, SCFTYPE, etc.)
4. Electronic State (FSYM, ALFA, BETA)
5. Basis Sets (A-, O-, P-, X- prefixes)
6. Termination (END)

---

## Header Section

The header defines basic calculation information:

```
TITLE
<calculation title>
SYMMETRY <symmetry_group>
CARTESIAN ANGSTROMS
```

### TITLE

Single line title for the calculation. Appears immediately after `TITLE` keyword.

```stobe
TITLE
ZnPc Ground State
```

### SYMMETRY

Point group symmetry. Common values:
- `C1`: No symmetry (most common for large molecules)
- `C2`, `C3`, `C4`, `C5`, `C6`: Cyclic groups
- `D2`, `D3`, `D4`, `D5`, `D6`: Dihedral groups
- `C2v`, `C3v`, `C4v`: With vertical mirror planes
- `D2h`, `D3h`, `D4h`: With horizontal mirror planes
- `T`, `Th`, `Td`, `O`, `Oh`: Cubic groups

```stobe
SYMMETRY C1
```

**Note**: Symmetry group must exist in `symbasis.new`. Use `NOSYMM` or `NOSYM` to disable symmetry.

### CARTESIAN

Geometry input format specification:

```stobe
CARTESIAN ANGSTROMS
```

or

```stobe
CARTESIAN BOHR
```

- `ANGSTROMS`: Coordinates in Angstroms (most common)
- `BOHR`: Coordinates in atomic units (Bohr)

---

## Geometry Section

Each atom is defined on a single line with the format:

```
<AtomLabel> <X> <Y> <Z> <EffectiveCharge> <GridPoints>
```

### Atom Label Format

- Element symbol (1-2 characters) + optional number
- Examples: `C01`, `H01`, `N01`, `Zn01`
- Element symbols must be valid (H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, etc.)
- Numbers allow distinguishing multiple atoms of same element

### Coordinate Format

- X, Y, Z: Cartesian coordinates
- Units determined by CARTESIAN keyword (ANGSTROMS or BOHR)
- Can be scientific notation: `1.0965640545`, `-2.7955413241`, `-0.1268843972`

### Effective Nuclear Charge

- `Z_eff`: Effective nuclear charge
- For normal atoms: Usually equals atomic number (6 for C, 1 for H, 7 for N, etc.)
- For core-hole atoms: Modified charge (e.g., 6 instead of 4 for carbon in X-ray calculations)
- For MCP: Reduced charge reflecting model core potential

### Grid Points

- `N_grid`: Number of radial grid points
- Typically `32` for most calculations
- Higher values (e.g., 45) for heavy elements or high accuracy
- Used for exchange-correlation potential evaluation

### Example Geometry Section

```stobe
C01 0.6501129622 -4.1882948611 -0.1553803524     6     32
C02 1.3304511823 -5.3932269664 -0.3028972717     4     32
H01 2.2754431526 -5.4172715017 -0.3985918603     1     32
N01 2.3771305165 -2.429005429 -0.2085415021     7     32
Zn01 0.0 0.0 0.0     30     32
END
```

**Critical**: Geometry section must end with `END` on its own line.

---

## Calculation Parameters Section

Parameters control SCF iteration, convergence, and calculation type:

```stobe
RUNTYPE startup nooptimize
SCFTYPE direct
POTENTIAL nonlocal rpbe pbe
GRID fine
MULTIPLICITY 1
CHARGE 0
MAXCYCLES 300
ECONVERGENCE 0.001
DCONVERGENCE 0.001
DMIXING mdens 0.05
DIIS new 7
ORBI 5d
MULLIKEN on full
VIRT all
SPIN FULL
```

### RUNTYPE

Calculation type:

- `startup nooptimize`: Single point SCF calculation
- `startup optimize`: Geometry optimization
- `restart`: Continue from restart file
- `startup`: Alias for `startup nooptimize`

### SCFTYPE

SCF algorithm type:

- `direct`: Direct SCF (memory intensive, fastest)
- `disk`: Disk-based SCF (slower, less memory)
- `memory`: Memory-based SCF

**Recommendation**: Use `direct` for systems with <1000 basis functions, `disk` for larger systems.

### POTENTIAL

Exchange-correlation functional:

- `nonlocal rpbe pbe`: RPBE functional (recommended for most cases)
- `local`: Local density approximation
- `nonlocal be88 pd86`: B88 exchange, PD86 correlation
- Other functionals available (see StoBeMAN.html)

### GRID

Integration grid quality:

- `fine`: High quality (default for accurate results)
- `medium`: Medium quality
- `coarse`: Low quality (faster, less accurate)

### MULTIPLICITY

Spin multiplicity (2S+1):

- `1`: Singlet (closed shell)
- `2`: Doublet (open shell, one unpaired electron)
- `3`: Triplet (two unpaired electrons)
- Higher multiplicities for high-spin systems

### CHARGE

Total molecular charge:

- `0`: Neutral molecule
- Positive: Cation (e.g., `1` for +1 charge)
- Negative: Anion (e.g., `-1` for -1 charge)

### MAXCYCLES

Maximum number of SCF iterations:

- Default: `300`
- Increase for difficult convergence (e.g., `500`, `1000`)
- Check output if convergence fails

### ECONVERGENCE

Energy convergence threshold:

- Default: `0.001` (Hartree)
- Tighter: `0.000001` for high accuracy
- Looser: `0.01` for quick tests

### DCONVERGENCE

Density convergence threshold:

- Default: `0.001`
- Tighter: `0.000001` for high accuracy
- Usually same as ECONVERGENCE

### DMIXING

Density mixing scheme:

- `mdens 0.05`: Modified density mixing with 0.05 mixing parameter
- `mdens 0.10`: Higher mixing (faster but less stable)
- `mdens 0.01`: Lower mixing (slower but more stable)
- Adjust if convergence problems occur

### DIIS

Direct Inversion in Iterative Subspace:

- `new 7`: Use DIIS with 7 previous iterations
- `on`: Use default DIIS
- `off`: Disable DIIS
- Helps accelerate SCF convergence

### ORBI

Orbital representation:

- `5d`: 5 d-functions (standard)
- `6d`: 6 d-functions (alternative)
- Affects d-orbital representation

### MULLIKEN

Mulliken population analysis:

- `on full`: Full Mulliken analysis
- `on`: Basic Mulliken analysis
- `off`: No Mulliken analysis

### VIRT

Virtual orbital output:

- `all`: Output all virtual orbitals
- `none`: No virtual orbitals
- Number: Output N lowest virtual orbitals

### SPIN

Spin contamination calculation:

- `FULL`: Calculate full spin contamination
- `OFF`: No spin contamination calculation

---

## Electronic State Section

Defines the electronic configuration and occupation.

### Ground State

```stobe
FSYM scfocc
ALFA <n_alpha>
BETA <n_beta>
FILE molden
END
```

- `FSYM scfocc`: Use SCF-occupied orbitals
- `ALFA <n>`: Number of alpha electrons
- `BETA <n>`: Number of beta electrons
- `FILE molden`: Output Molden format file

**Example**:
```stobe
FSYM scfocc
ALFA 147
BETA 147
FILE molden
END
```

### Excited State

```stobe
FSYM scfocc excited
ALFA <n_alpha>
BETA <n_beta>
SYM 1
ALFA <occupation_string>
BETA <occupation_string>
END
```

- `FSYM scfocc excited`: Excited state calculation
- `ALFA <n+1>`: One more alpha electron (core-hole created)
- `BETA <n>`: Same beta electrons
- `SYM 1`: Symmetry block (usually 1)
- Occupation strings define hole state

**Example**:
```stobe
FSYM scfocc excited
ALFA 117
BETA 116
SYM 1
ALFA 0 1 14 0.0
BETA 0 0
END
```

### Transition Potential

```stobe
FSYM scfocc excited
ALFA <n_alpha>
BETA <n_beta>
SYM 1
ALFA <occupation_string>
BETA <occupation_string>
END
MULLIKEN on full
FILE molden
XRAY xas
END
END
```

- Similar to excited state but with `XRAY xas` keyword
- Occupation uses 0.5 (half-electron)
- Produces transition matrix elements

**Example**:
```stobe
FSYM scfocc excited
ALFA 116
BETA 116
SYM 1
ALFA 0 1 14 0.5
BETA 0 0
END
MULLIKEN on full
FILE molden
XRAY xas
END
END
```

### Occupation Strings

Format: `0 <n_occ> <n_virt> <occupation>`

- `0`: Starting orbital index
- `<n_occ>`: Number of occupied orbitals in hole state
- `<n_virt>`: Number of virtual orbitals to include
- `<occupation>`: Occupation number
  - `0.0`: Empty (excited state)
  - `0.5`: Half-filled (transition potential)
  - `1.0`: Filled (ground state)

**Examples**:
- `0 1 14 0.0`: One occupied orbital, 14 virtuals, empty (excited state)
- `0 1 14 0.5`: One occupied orbital, 14 virtuals, half-filled (transition potential)
- `0 1 1 0.0`: One occupied, one virtual, empty

---

## Basis Set Section

Basis sets are specified in order matching the geometry. One basis set per atom.

### Section Order

1. Auxiliary basis sets (A- prefix)
2. Orbital basis sets (O- prefix)
3. Model core potential basis sets (P- prefix, optional)
4. Augmentation basis sets (X- prefix, for X-ray calculations)

### Critical Rules

- **One basis set per atom**: Must match geometry atom count exactly
- **Order matters**: Basis sets must match geometry order
- **Element matching**: Basis set element must match atom element
- **All atoms require**: Auxiliary and orbital basis sets
- **MCP optional**: Only for heavy elements or core-hole atoms
- **Augmentation optional**: Only for X-ray calculations

---

## File Termination

All input files must end with:

```stobe
END
```

**Critical**: Missing END markers cause parsing errors.

---

## Complete Example

```stobe
TITLE
Test Molecule Ground State
SYMMETRY C1
CARTESIAN ANGSTROMS
C01 0.0 0.0 0.0     6     32
H01 1.0 0.0 0.0     1     32
H02 -1.0 0.0 0.0     1     32
END
RUNTYPE startup nooptimize
SCFTYPE direct
POTENTIAL nonlocal rpbe pbe
GRID fine
MULTIPLICITY 1
CHARGE 0
MAXCYCLES 300
ECONVERGENCE 0.001
DCONVERGENCE 0.001
DMIXING mdens 0.05
DIIS new 7
ORBI 5d
MULLIKEN on full
VIRT all
SPIN FULL
FSYM scfocc
ALFA 4
BETA 4
FILE molden
END
A-CARBON (5,2;5,2)
A-HYDROGEN (4,2;4,2)
A-HYDROGEN (4,2;4,2)
O-CARBON iii_iglo
O-HYDROGEN (311/1) misc
O-HYDROGEN (311/1) misc
END
```

---

## Validation Checklist

When generating or parsing input files, verify:

- [ ] TITLE present
- [ ] SYMMETRY specified (or NOSYMM)
- [ ] CARTESIAN with units specified
- [ ] Geometry section has all atoms
- [ ] Geometry section ends with END
- [ ] All required parameters present
- [ ] Electronic state section matches calculation type
- [ ] Electronic state section ends with END
- [ ] Basis set count matches atom count
- [ ] Basis sets match element types
- [ ] File ends with END
- [ ] No syntax errors
- [ ] Effective charges appropriate
- [ ] Occupation numbers match multiplicity

---

## Common Errors

### Mismatched Counts

**Error**: Geometry has N atoms but M basis sets

**Fix**: Ensure exactly one basis set per atom for each type (A-, O-, P-, X-)

### Missing END Markers

**Error**: Sections not terminated

**Fix**: Add END after geometry, after electronic state, and at file end

### Invalid Symmetry

**Error**: Symmetry group not found

**Fix**: Check `symbasis.new` for valid groups, or use `NOSYMM`

### Invalid Basis Format

**Error**: Basis specification syntax error

**Fix**: Check format against examples, verify element name spelling

### Charge Inconsistency

**Error**: Total charge doesn't match electron count

**Fix**: Verify ALFA + BETA = total electrons - CHARGE

---

## Parsing Strategy

```python
def parse_input_file(filename: str) -> dict:
    """Parse StoBe input file into structured format."""
    sections = {
        'title': None,
        'symmetry': None,
        'units': None,
        'geometry': [],
        'parameters': {},
        'electronic': {},
        'basis_sets': {
            'auxiliary': [],
            'orbital': [],
            'mcp': [],
            'augmentation': []
        }
    }
    
    current_section = 'header'
    basis_type = None
    
    with open(filename, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            
            if not line or line.startswith('#'):
                continue
            
            # Header parsing
            if line == 'TITLE':
                sections['title'] = next(f).strip()
            elif line.startswith('SYMMETRY'):
                sections['symmetry'] = line.split()[1] if len(line.split()) > 1 else None
            elif 'CARTESIAN' in line:
                sections['units'] = 'ANGSTROMS' if 'ANGSTROMS' in line else 'BOHR'
                current_section = 'geometry'
            
            # Geometry parsing
            elif current_section == 'geometry':
                if line == 'END':
                    current_section = 'parameters'
                else:
                    atom_data = parse_atom_line(line)
                    sections['geometry'].append(atom_data)
            
            # Parameters parsing
            elif current_section == 'parameters':
                if line.startswith('FSYM'):
                    current_section = 'electronic'
                    sections['parameters']['fsym'] = line.split()[1] if len(line.split()) > 1 else None
                elif '=' not in line and len(line.split()) >= 2:
                    key, *values = line.split()
                    sections['parameters'][key.lower()] = ' '.join(values) if values else None
            
            # Electronic state parsing
            elif current_section == 'electronic':
                if line == 'END':
                    current_section = 'basis'
                    basis_type = 'auxiliary'
                elif line.startswith('ALFA') or line.startswith('BETA'):
                    key, *values = line.split()
                    sections['electronic'][key.lower()] = values
                elif line.startswith('FSYM'):
                    sections['electronic']['fsym'] = line.split()[1] if len(line.split()) > 1 else None
                elif line.startswith('SYM'):
                    sections['electronic']['sym'] = line.split()[1] if len(line.split()) > 1 else None
                elif line.startswith('FILE'):
                    sections['electronic']['file'] = line.split()[1] if len(line.split()) > 1 else None
                elif line.startswith('XRAY'):
                    sections['electronic']['xray'] = line.split()[1] if len(line.split()) > 1 else None
            
            # Basis set parsing
            elif current_section == 'basis':
                if line == 'END':
                    break
                elif line.startswith('A-'):
                    sections['basis_sets']['auxiliary'].append(line)
                elif line.startswith('O-'):
                    sections['basis_sets']['orbital'].append(line)
                elif line.startswith('P-'):
                    sections['basis_sets']['mcp'].append(line)
                elif line.startswith('X-'):
                    sections['basis_sets']['augmentation'].append(line)
    
    return sections

def parse_atom_line(line: str) -> dict:
    """Parse atom geometry line."""
    parts = line.split()
    if len(parts) < 6:
        raise ValueError(f"Invalid atom line: {line}")
    
    label = parts[0]
    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
    z_eff = float(parts[4])
    n_grid = int(parts[5])
    
    return {
        'label': label,
        'coordinates': (x, y, z),
        'effective_charge': z_eff,
        'grid_points': n_grid
    }
```

---

## Related References

- **Basis Sets**: See `basis-sets.md` for detailed basis set specifications
- **Run Scripts**: See `run-scripts.md` for generating .run files
- **Symmetry**: See `symmetry.md` for symmetry group usage
- **Examples**: See `examples-catalog.md` for input file examples
