# X-ray Spectroscopy

> Reference for: StoBe Pro
> Load when: Calculating X-ray absorption spectra, transition potentials, core-hole states, NEXAFS/XANES calculations

---

## Overview

StoBe specializes in X-ray absorption spectroscopy calculations using the transition potential method. This reference covers core-hole calculations, transition potential setup, and spectrum generation for NEXAFS/XANES applications.

---

## X-ray Absorption Theory

### Core-Hole Excitation

X-ray absorption involves:
1. Core electron excited to unoccupied orbitals
2. Core-hole created in 1s orbital
3. Transition matrix elements calculated
4. Spectrum generated with broadening

### Transition Potential Method

- Half-electron in core orbital (0.5 occupation)
- More efficient than full core-hole calculation
- Produces accurate transition energies
- Standard method for X-ray calculations

---

## Core-Hole Calculations

### Effective Nuclear Charge Modification

For core-hole atom:
- **Standard charge**: Use atomic number (e.g., 6 for carbon)
- **Other atoms**: Use modified charge (e.g., 4 for carbon with MCP)

**Example**:
```stobe
C01 0.6501129622 -4.1882948611 -0.1553803524     6     32
C02 1.3304511823 -5.3932269664 -0.3028972717     4     32
```

- C01: Core-hole atom (charge 6)
- C02: Other carbon (charge 4, using MCP)

### Basis Set Modifications

**Core-hole atom**:
- Auxiliary: `A-CARBON (5,2;5,2)` (standard)
- Orbital: `O-CARBON iii_iglo` (standard)
- MCP: `P-CARBON(+4) (3,1:8,0)` (modified charge)
- Augmentation: `X-FIRST`

**Other atoms**:
- Auxiliary: `A-CARBON(+4) (3,3;3,3)` (modified charge)
- Orbital: `O-CARBON(+4) (311/211/1)` (modified charge)
- MCP: `P-CARBON(+4) (3,1:8,0)` (modified charge)
- Augmentation: `X-DUMMY`

---

## Transition Potential Setup

### Input File Structure

```stobe
TITLE
Molecule TP
SYMMETRY C1
CARTESIAN ANGSTROMS
[geometry with modified charges]
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
FSYM scfocc excited
ALFA <n_alpha>
BETA <n_beta>
SYM 1
ALFA 0 1 <n_virt> 0.5
BETA 0 0
END
MULLIKEN on full
FILE molden
XRAY xas
END
END
[basis sets]
END
```

### Key Parameters

**Electronic State**:
```stobe
FSYM scfocc excited
ALFA 116
BETA 116
SYM 1
ALFA 0 1 14 0.5
BETA 0 0
```

- `FSYM scfocc excited`: Excited state calculation
- `ALFA 0 1 14 0.5`: Half-electron occupation
  - `0`: Starting orbital index
  - `1`: One occupied orbital (core orbital)
  - `14`: Number of virtual orbitals
  - `0.5`: Half-electron occupation

**X-ray Keyword**:
```stobe
XRAY xas
```

- Required for transition potential calculations
- Produces transition matrix elements
- Outputs to `fort.11` (moved to `.xas` file)

---

## X-ray Spectrum Generation

### Post-Processing Input

After transition potential calculation, use `xrayspec.x` utility:

```stobe
title
Molecule XAS
PRINT
RANGE 280 320
POINTS 2000
WIDTH 0.5 12 288 320
XRAY xas
TOTAL 1
END
```

### Parameters

**RANGE**: Energy range in eV
- Carbon K-edge: `280 320`
- Nitrogen K-edge: `400 450`
- Oxygen K-edge: `530 580`

**POINTS**: Number of energy points
- Typical: `2000`
- Higher for better resolution

**WIDTH**: Broadening parameters
- Format: `<width> <exponent> <E1> <E2>`
- Example: `0.5 12 288 320`
  - `0.5`: Base width
  - `12`: Exponent
  - `288 320`: Energy range for broadening

**TOTAL**: Total spectrum calculation
- `1`: Calculate total spectrum
- Other options for partial spectra

### Run Script

```bash
#!/bin/csh -f
ln -s ${PWD}/C1.xas fort.1
cat >C1xas.inp<</.
title
ZnPc XAS
PRINT
RANGE 280 320
POINTS 2000
WIDTH 0.5 12 288 320
XRAY xas
TOTAL 1
END
/.
xrayspec.x <C1xas.inp>& C1xas.out
mv XrayT001.out C1.out
rm fort.*
```

---

## Complete NEXAFS Workflow

### Step-by-Step

1. **Ground State**
   ```bash
   ./C1gnd.run
   ```
   - Produces: `C1gnd.out`, `C1gnd.molden`

2. **Excited State** (optional)
   ```bash
   ./C1exc.run
   ```
   - Produces: `C1exc.out`, `C1exc.molden`

3. **Transition Potential**
   ```bash
   ./C1tp.run
   ```
   - Produces: `C1tp.out`, `C1tp.molden`, `C1.xas`

4. **X-ray Spectrum**
   ```bash
   ./C1xas.run
   ```
   - Produces: `C1xas.out`, `C1.out` (spectrum)

### Sequential Script

```bash
#!/bin/csh -f
./C1gnd.run
./C1exc.run
./C1tp.run
./C1xas.run
```

---

## Multi-Atom Calculations

### Pattern: One Atom Per Calculation

For molecules with multiple equivalent atoms:

```bash
# Carbon 1
./C1gnd.run
./C1tp.run
./C1xas.run

# Carbon 2
./C2gnd.run
./C2tp.run
./C2xas.run

# Carbon 3
./C3gnd.run
./C3tp.run
./C3xas.run
```

### Effective Charge Pattern

- First atom (core-hole): Charge 6 (carbon)
- Other atoms: Charge 4 (carbon with MCP)

**Geometry example**:
```stobe
C01 ...     6     32
C02 ...     4     32
C03 ...     4     32
...
```

### Basis Set Pattern

- First atom: Standard basis sets, `X-FIRST`
- Other atoms: Modified charge basis sets, `X-DUMMY`

---

## Energy Ranges

### Common K-edges

**Carbon K-edge**:
- Range: 280-320 eV
- Width: `0.5 12 288 320`

**Nitrogen K-edge**:
- Range: 400-450 eV
- Width: `0.5 12 418 450`

**Oxygen K-edge**:
- Range: 530-580 eV
- Width: `0.5 12 535 560`

### L-edges

For transition metals:
- L3 edge: Lower energy
- L2 edge: Higher energy
- Adjust range based on element

---

## Transition Matrix Elements

### Output Format

Transition potential calculation produces:
- Transition energies
- Transition moments
- Orbital information
- Stored in `.xas` file (fort.11)

### Processing

`xrayspec.x` reads `.xas` file and generates:
- Broadened spectrum
- Energy vs. intensity
- Output in `XrayT001.out`

---

## Advanced Techniques

### Orbital Localization

For symmetric molecules:
```stobe
LOCALIZE
```

- Localizes core orbitals
- Helps identify transitions
- Useful for analysis

### Supersymmetry Freezing

```stobe
SUPSYM
```

- Freezes non-excited core orbitals
- Reduces computational cost
- Maintains accuracy

### Multiple Transitions

Calculate transitions to multiple virtual orbitals:
```stobe
ALFA 0 1 20 0.5
```

- `20`: Number of virtual orbitals
- More virtuals = more transitions
- Higher computational cost

---

## Validation

For X-ray calculations, verify:

- [ ] Effective charges modified correctly
- [ ] Basis sets match core-hole pattern
- [ ] Augmentation basis correct (X-FIRST/X-DUMMY)
- [ ] Occupation uses 0.5 for transition potential
- [ ] XRAY xas keyword present
- [ ] Energy range appropriate for edge
- [ ] Transition potential completed before spectrum

---

## Common Errors

### Wrong Effective Charge

**Error**: Core-hole atom uses modified charge

**Fix**: Core-hole atom should use standard charge (6 for C), others use modified (4 for C)

### Wrong Basis Sets

**Error**: Core-hole atom uses modified charge basis sets

**Fix**: Core-hole atom uses standard basis sets, others use modified

### Missing XRAY Keyword

**Error**: Transition potential doesn't produce .xas file

**Fix**: Add `XRAY xas` keyword in electronic state section

### Wrong Occupation

**Error**: Using 0.0 instead of 0.5

**Fix**: Transition potential requires 0.5 occupation

---

## Related References

- **Calculation Workflows**: See `calculation-workflows.md` for complete workflow patterns
- **Basis Sets**: See `basis-sets.md` for core-hole basis set patterns
- **Input File Format**: See `input-file-format.md` for transition potential input structure
- **Examples**: See `examples-catalog.md` for X-ray calculation examples
