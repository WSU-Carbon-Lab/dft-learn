# Examples Catalog

> Reference for: StoBe Pro
> Load when: Finding example calculations, learning workflows, understanding usage patterns, locating reference implementations

---

## Overview

This catalog indexes StoBe calculation examples from multiple sources, organized by calculation type, system type, and features. Examples demonstrate best practices and common patterns for various calculation types.

---

## Example Locations

### Official StoBe Examples

**Location**: `$STOBE_HOME/StoBeLinux/StoBeExport/Examples`

The official StoBe distribution includes 32 example calculations demonstrating various features and workflows. These examples are well-documented and include reference output files for validation.

**Access Pattern**:
```bash
# Set STOBE_HOME environment variable or use full path
cd $STOBE_HOME/StoBeLinux/StoBeExport/Examples

# Or use full path
cd /home/hduva/stobe/StoBeLinux/StoBeExport/Examples
```

**File Structure**:
- Each example has a `.run` file (executable shell script)
- Reference output files (`.outref`) for comparison
- Some examples generate additional files (`.res`, `.molf`, `.eps`, etc.)
- `Examples.txt` file contains descriptions of all examples

### Project-Specific Examples

**Location**: `~/projects/dft/`

Real-world calculation examples from active research projects, including:
- Multi-atom NEXAFS calculations
- Complex molecular systems
- Production workflows

**Key Directories**:
- `~/projects/dft/geometry/`: Geometry optimization examples
- `~/projects/dft/tcta_rotated_phyn/`: TCTA molecule calculations
- `~/projects/dft/basis/old/`: Basis set testing examples

**Access Pattern**:
```bash
cd ~/projects/dft/geometry/C1
# Contains: C1gnd.run, C1exc.run, C1tp.run, C1xas.run, C1seq.run
```

### Local Examples

**Location**: `.cursor/skills/stobe-pro/references/examples/`

Copied examples for reference:
- C1, C2, C3 symmetry examples
- Ground, excited, transition potential, X-ray spectrum calculations

---

## Official StoBe Examples Catalog

### Basic SCF Calculations

#### NH3nosym
**File**: `NH3nosym.run`
**Description**: Ammonia molecule, no symmetry, single point SCF
**Key Features**:
- Basic SCF calculation
- No symmetry (C1)
- Demonstrates minimal input requirements

**Use When**: Learning basic StoBe input format

**Location**: `$STOBE_HOME/StoBeLinux/StoBeExport/Examples/NH3nosym.run`

#### NH3ion
**File**: `NH3ion.run`
**Description**: Ammonia molecular ion (NH3+), C3v symmetry, single point SCF
**Key Features**:
- Charged system
- C3v symmetry usage
- Ion calculation

**Use When**: Working with charged systems or exploiting symmetry

**Location**: `$STOBE_HOME/StoBeLinux/StoBeExport/Examples/NH3ion.run`

#### H2ODolg
**File**: `H2ODolg.run`
**Description**: Water molecule using Dolg pseudopotential for oxygen
**Key Features**:
- Pseudopotential usage
- No symmetry
- Experimental feature demonstration

**Use When**: Learning pseudopotential applications

**Location**: `$STOBE_HOME/StoBeLinux/StoBeExport/Examples/H2ODolg.run`

### Geometry Optimization

#### H2Oopt
**File**: `H2Oopt.run`
**Description**: Water molecule geometry optimization, C2v symmetry
**Key Features**:
- Full geometry optimization
- Symmetry exploitation
- Molden format output

**Use When**: Optimizing molecular geometries

**Location**: `$STOBE_HOME/StoBeLinux/StoBeExport/Examples/H2Oopt.run`

#### CO2optvibs
**File**: `CO2optvibs.run`
**Description**: CO2 geometry optimization with vibrational analysis
**Key Features**:
- Multi-step workflow (4 steps)
- Geometry optimization → vibrational analysis
- Symmetry reduction (D4h → C1)
- Vibrational frequency calculation

**Use When**: Calculating vibrational frequencies

**Location**: `$STOBE_HOME/StoBeLinux/StoBeExport/Examples/CO2optvibs.run`

### X-ray Absorption Spectroscopy

#### PyrxasN
**File**: `PyrxasN.run`
**Description**: Pyridine X-ray absorption for N1s excitation
**Key Features**:
- Core-hole calculation
- Transition moment calculation
- Nitrogen K-edge

**Use When**: Calculating NEXAFS spectra for nitrogen

**Location**: `$STOBE_HOME/StoBeLinux/StoBeExport/Examples/PyrxasN.run`

#### PyrxasC
**File**: `PyrxasC.run`
**Description**: Pyridine X-ray absorption for C1s excitation using model potentials
**Key Features**:
- Transition potential method
- Model core potential for non-excited carbons
- Carbon K-edge calculation
- Demonstrates efficient multi-atom approach

**Use When**: Calculating carbon K-edge NEXAFS with MCP

**Location**: `$STOBE_HOME/StoBeLinux/StoBeExport/Examples/PyrxasC.run`

**Key Input Features**:
```stobe
fsym scfocc excited
alfa 17
beta 17
sym 1
alfa 0 1 2 0.5
beta 0 0
xray xas
```

#### PyridGS_TP
**File**: `PyridGS_TP.run`
**Description**: Pyridine ground state and C1s excited states using transition potential approach
**Key Features**:
- Transition potential method
- Multiple symmetry-equivalent atoms
- Orbital localization
- Supersymmetry freezing
- Efficient calculation strategy

**Use When**: Calculating XAS for molecules with symmetry

**Location**: `$STOBE_HOME/StoBeLinux/StoBeExport/Examples/PyridGS_TP.run`

### Property Calculations

#### NH3res
**File**: `NH3res.run`
**Description**: Ammonia molecular ion with property calculations
**Key Features**:
- Mulliken population analysis
- Loewdin population analysis
- Topological atom charges (Bader, Voronoi, Becke)
- Orbital plot output
- Requires NH3ion to be run first

**Use When**: Calculating atomic charges and populations

**Location**: `$STOBE_HOME/StoBeLinux/StoBeExport/Examples/NH3res.run`

#### Pyrpol
**File**: `Pyrpol.run`
**Description**: Pyridine ground state with polarizability calculation
**Key Features**:
- Multi-step workflow (3 steps)
- Symmetry reduction (C2v → C1)
- Polarizability tensor calculation
- Molden format output

**Use When**: Calculating molecular polarizabilities

**Location**: `$STOBE_HOME/StoBeLinux/StoBeExport/Examples/Pyrpol.run`

### Advanced Workflows

#### C2H6NEB
**File**: `C2H6NEB.run`
**Description**: Nudged Elastic Band (NEB) path for ethane rotational barrier
**Key Features**:
- NEB method
- Reaction path calculation
- Multiple images
- Separate restart files per image

**Use When**: Calculating reaction paths and barriers

**Location**: `$STOBE_HOME/StoBeLinux/StoBeExport/Examples/C2H6NEB.run`

#### V2O9H8H
**File**: `V2O9H8H.run`
**Description**: Complex adsorption example with 6 steps
**Key Features**:
- Multi-step workflow
- Restart file combination
- Geometry optimization of adsorbate
- Orbital expansion analysis
- Requires significant computation time

**Use When**: Modeling surface adsorption

**Location**: `$STOBE_HOME/StoBeLinux/StoBeExport/Examples/V2O9H8H.run`

---

## Project Examples

### Geometry Directory

**Location**: `~/projects/dft/geometry/`

#### C1 Examples
**Directory**: `~/projects/dft/geometry/C1/`
**Files**:
- `C1gnd.run`: Ground state calculation
- `C1exc.run`: Excited state calculation
- `C1tp.run`: Transition potential calculation
- `C1xas.run`: X-ray spectrum calculation
- `C1seq.run`: Sequential batch script
- `C1gnd.inp`: Input file example

**System**: Zinc phthalocyanine (ZnPc)
**Features**:
- Large molecule (65 atoms)
- C1 symmetry
- Complete NEXAFS workflow
- Multiple calculation types

**Use When**: 
- Learning complete NEXAFS workflow
- Working with large molecules
- Multi-step calculation sequences

### TCTA Examples

**Location**: `~/projects/dft/tcta_rotated_phyn/`

#### C9, C10, C16 Examples
**Directories**: `C9/`, `C10/`, `C16/`, etc.
**Files per directory**:
- `C<number>gnd.run`: Ground state
- `C<number>exc.run`: Excited state
- `C<number>tp.run`: Transition potential
- `C<number>xas.run`: X-ray spectrum
- `C<number>seq.run`: Sequential script
- `.inp` files for each calculation type

**System**: TCTA (tris(4-carbazoyl-9-ylphenyl)amine) molecule
**Features**:
- Multiple carbon atoms
- Rotated geometry
- Production calculations
- Complete workflow for each atom

**Use When**:
- Calculating NEXAFS for multiple equivalent atoms
- Production workflows
- Batch processing examples

---

## Finding Examples

### By Calculation Type

**Ground State Calculations**:
- Official: `NH3nosym.run`, `NH3ion.run`, `H2ODolg.run`
- Project: `~/projects/dft/geometry/C1/C1gnd.run`
- Project: `~/projects/dft/tcta_rotated_phyn/C*/C*gnd.run`

**Geometry Optimization**:
- Official: `H2Oopt.run`, `CO2optvibs.run`, `Au2.run`

**X-ray Absorption**:
- Official: `PyrxasN.run`, `PyrxasC.run`, `PyridGS_TP.run`
- Project: `~/projects/dft/geometry/C1/C1tp.run`
- Project: `~/projects/dft/tcta_rotated_phyn/C*/C*tp.run`

**Property Calculations**:
- Official: `NH3res.run` (populations)
- Official: `H2Ofukup.run`, `H2Ofukum.run` (Fukui functions)
- Official: `Pyrpol.run` (polarizability)

### By System Type

**Small Molecules**:
- `NH3*.run` (ammonia)
- `H2O*.run` (water)
- `CO2*.run` (carbon dioxide)

**Aromatic Molecules**:
- `Pyr*.run` (pyridine)

**Large Molecules**:
- `~/projects/dft/geometry/C1/*.run` (ZnPc)
- `~/projects/dft/tcta_rotated_phyn/*/*.run` (TCTA)

### By Feature

**Symmetry Usage**:
- `NH3ion.run` (C3v)
- `H2Oopt.run` (C2v)
- `CO2optvibs.run` (D4h)

**Model Core Potentials**:
- `PyrxasC.run` (carbon MCP)
- `Auatom.run`, `Au2.run` (gold MCP)

**Multi-Step Workflows**:
- `CO2optvibs.run` (4 steps)
- `Pyrpol.run` (3 steps)
- `V2O9H8H.run` (6 steps)

---

## Running Examples

### Prerequisites

1. StoBe installation accessible
2. Basis set files available (`baslib.new7`, `symbasis.new`)
3. Proper environment variables set (if using `$STOBE_HOME`)

### Execution

```bash
# Make executable
chmod +x example.run

# Run example
./example.run

# Check output
cat example.out
```

### Validation

Compare output with reference:
```bash
diff example.out example.outref
```

---

## Related References

- **Input File Format**: See `input-file-format.md` for input structure
- **Run Scripts**: See `run-scripts.md` for script creation
- **Calculation Workflows**: See `calculation-workflows.md` for workflow patterns
- **X-ray Spectroscopy**: See `xray-spectroscopy.md` for X-ray examples
