# Quick Reference

> Reference for: StoBe Pro
> Load when: Quick keyword lookup, default parameters, validation checklist, common templates

---

## Overview

Quick reference for common StoBe keywords, parameters, and templates. Use this for fast lookup during input file generation or validation.

---

## File Structure Quick Reference

### Input File (.inp) Sections

```
TITLE
<title>
SYMMETRY <group>
CARTESIAN ANGSTROMS
<geometry>
END
<parameters>
<electronic_state>
<basis_sets>
END
```

### Run Script (.run) Template

```bash
#!/bin/csh -f
ln -s ${STOBE}Basis/baslib.new7 fort.3
ln -s ${STOBE}Basis/symbasis.new fort.4
cat ><name>.inp<</.
[input content]
/.
StoBe.x <<name>.inp>& <name>.out
mv Molden.molf <name>.molden
rm fort.*
```

## Common Keywords

### Geometry
- `CARTESIAN ANGSTROMS`: Cartesian coordinates in Angstroms
- `SYMMETRY <group>`: Point group symmetry

### Calculation Type
- `RUNTYPE startup nooptimize`: Single point
- `RUNTYPE startup optimize`: Geometry optimization
- `SCFTYPE direct`: Direct SCF
- `POTENTIAL nonlocal rpbe pbe`: RPBE functional
- `GRID fine`: Fine integration grid

### Convergence
- `MAXCYCLES 300`: Maximum SCF cycles
- `ECONVERGENCE 0.001`: Energy convergence
- `DCONVERGENCE 0.001`: Density convergence
- `DMIXING mdens 0.05`: Density mixing

### Electronic State
- `MULTIPLICITY 1`: Singlet state
- `CHARGE 0`: Neutral molecule
- `FSYM scfocc`: Ground state
- `FSYM scfocc excited`: Excited state
- `ALFA <n>`: Alpha electrons
- `BETA <n>`: Beta electrons
- `XRAY xas`: X-ray calculation

### Properties
- `MULLIKEN on full`: Mulliken analysis
- `FILE molden`: Molden output
- `VIRT all`: Virtual orbitals

## Basis Set Prefixes

- `A-`: Auxiliary basis (density fitting)
- `O-`: Orbital basis (molecular orbitals)
- `P-`: Model core potential
- `X-`: Augmentation (X-FIRST or X-DUMMY)

## Calculation Types

### Ground State
```
FSYM scfocc
ALFA <n>
BETA <n>
FILE molden
END
```

### Excited State
```
FSYM scfocc excited
ALFA <n+1>
BETA <n>
SYM 1
ALFA 0 1 <n_virt> 0.0
BETA 0 0
END
```

### Transition Potential
```
FSYM scfocc excited
ALFA <n>
BETA <n>
SYM 1
ALFA 0 1 <n_virt> 0.5
BETA 0 0
END
MULLIKEN on full
FILE molden
XRAY xas
END
END
```

## Common Symmetry Groups

- `C1`: No symmetry
- `C2`: Two-fold rotation
- `C3`: Three-fold rotation
- `D2h`: Full octahedral
- `Td`: Tetrahedral
- `Oh`: Octahedral

## Effective Nuclear Charges

Common values:
- Carbon: 6 (normal), 4 (core-hole)
- Nitrogen: 7 (normal), 5 (core-hole)
- Oxygen: 8 (normal), 6 (core-hole)
- Hydrogen: 1
- Transition metals: Atomic number

## Grid Points

Typically: `32` radial grid points per atom

## Default Parameters

```python
DEFAULT_SYM = "C1"
DEFAULT_multiplicity = '1'
DEFAULT_GEOUNITS = 'ANGSTROMS'
DEFAULT_RUNTYPE = 'startup nooptimize'
DEFAULT_SCFTYPE = 'direct'
DEFAULT_POTENTIAL = 'nonlocal rpbe pbe'
DEFAULT_GRID = 'fine'
DEFAULT_CHARGE = '0'
DEFAULT_MAXCYCLES = '300'
DEFAULT_ECONVERGENCE = '0.001'
DEFAULT_DCONVERGENCE = '0.001'
DEFAULT_DMIXING = 'mdens 0.05'
DEFAULT_DIIS = 'new 7'
DEFAULT_ORBI = '5d'
DEFAULT_MULLIKEN = 'on full'
DEFAULT_VIRT = 'all'
DEFAULT_FSYMGND = 'scfocc'
DEFAULT_MOFILE = 'molden'
DEFAULT_FSYMEXC = 'scfocc excited'
DEFAULT_ALPHAOCC = '0 1 1 0.0'
DEFAULT_BETAOCC = '0 0'
DEFAULT_ALPHAOCCTP = '0 1 1 0.5'
DEFAULT_SPIN = "FULL"
DEFAULT_XRAY = "XAS"
```

## File Naming Conventions

- `*gnd.run`: Ground state calculation
- `*exc.run`: Excited state calculation
- `*tp.run`: Transition potential calculation
- `*xas.run`: X-ray spectrum calculation
- `*seq.run`: Sequential batch script

## Output Files

- `*.out`: Main output file
- `*.molden`: Molden format orbitals
- `*.xas`: Transition potential data (fort.11)
- `XrayT001.out`: X-ray spectrum output

## Common Basis Sets

### Carbon
- Auxiliary: `A-CARBON (5,2;5,2)` or `A-CARBON(+4) (3,3;3,3)`
- Orbital: `O-CARBON iii_iglo` or `O-CARBON(+4) (311/211/1)`
- MCP: `P-CARBON(+4) (3,1:8,0)`

### Hydrogen
- Auxiliary: `A-HYDROGEN (4,2;4,2)`
- Orbital: `O-HYDROGEN (311/1) misc`

### Nitrogen
- Auxiliary: `A-NITROGEN (4,3;4,3)`
- Orbital: `O-NITROGEN (33/3)`

### Oxygen
- Auxiliary: `A-OXYGEN (4,3;4,3)`
- Orbital: `O-OXYGEN (33/3)`

## X-ray Energy Ranges

Common ranges for post-processing:
- Carbon K-edge: 280-320 eV
- Nitrogen K-edge: 400-450 eV
- Oxygen K-edge: 530-580 eV

## Validation Checklist

- [ ] Geometry matches basis set count
- [ ] Symmetry group is valid
- [ ] All sections have END markers
- [ ] Effective charges are correct
- [ ] Occupation matches multiplicity
- [ ] Basis sets match elements
- [ ] File paths are correct

## Common Errors

1. **Mismatched counts**: Atom count ≠ basis set count
2. **Invalid symmetry**: Group not in symbasis.new
3. **Missing END**: Sections not terminated
4. **Wrong charge**: Doesn't match electron count
5. **Basis syntax**: Invalid format
6. **Convergence**: SCF didn't converge

## Workflow Order

1. Ground state (gnd)
2. Excited state (exc)
3. Transition potential (tp)
4. X-ray spectrum (xas)

## Python Integration

See `scripts/StoBe_Input_Generator.py` and `scripts/molConfig.py` for examples.
