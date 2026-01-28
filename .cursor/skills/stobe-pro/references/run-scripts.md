# Run Scripts

> Reference for: StoBe Pro
> Load when: Creating or modifying .run shell scripts, setting up calculations, automating execution

---

## Overview

Run scripts (.run) are C shell scripts that automate StoBe calculation execution. They handle basis set linking, input file generation, program execution, and output processing.

---

## Standard Structure

Every run script follows this pattern:

```bash
#!/bin/csh -f
ln -s <basis_path>/baslib.new7 fort.3
ln -s <basis_path>/symbasis.new fort.4
cat ><filename>.inp<</.
[input file content]
/.
<stobe_executable> <<filename>.inp>& <filename>.out
[output processing]
rm fort.*
```

---

## Key Components

### 1. Shebang

```bash
#!/bin/csh -f
```

- Uses C shell (`csh`)
- `-f` flag: Fast startup (skip .cshrc)

### 2. Basis Set Library Links

```bash
ln -s ${STOBE}Basis/baslib.new7 fort.3
ln -s ${STOBE}Basis/symbasis.new fort.4
```

**Alternative paths**:
```bash
# Using environment variable
ln -s ${STOBE}Basis/baslib.new7 fort.3

# Using absolute path
ln -s /home/hduva/stobe/StoBeLinux/Basis/baslib.new7 fort.3

# Using relative path
ln -s ../Basis/baslib.new7 fort.3
```

**File mappings**:
- `fort.3` → `baslib.new7` (basis set library)
- `fort.4` → `symbasis.new` (symmetry definitions)

### 3. Input File Generation

```bash
cat ><filename>.inp<</.
[input content here]
/.
```

- Uses heredoc syntax (`<</.` ... `/.`)
- Generates input file from script content
- Input content follows standard .inp format

### 4. Execution

```bash
StoBe.x <<filename>.inp>& <filename>.out
```

**Variations**:
```bash
# Standard execution
StoBe.x <input.inp>& output.out

# With full path
~/STOBE/Source/StoBe.x <input.inp>& output.out

# With environment variable
${STOBE}Source/StoBe.x <input.inp>& output.out
```

**Redirection**:
- `<input.inp`: Input file
- `>& output.out`: Redirect stdout and stderr

### 5. Output Processing

```bash
mv Molden.molf <filename>.molden
mv fort.11 <filename>.xas
rm fort.*
```

Common output files:
- `Molden.molf` → Rename to `.molden`
- `fort.11` → X-ray transition data (for TP calculations)
- `fort.2` → Restart file (optional, usually not moved)
- `fort.*` → Clean up temporary files

---

## Calculation Types

### Ground State (gnd.run)

**Purpose**: Standard SCF calculation

**Outputs**:
- `.out`: Main output file
- `.molden`: Molden format orbitals

**Example**:
```bash
#!/bin/csh -f
ln -s ${STOBE}Basis/baslib.new7 fort.3
ln -s ${STOBE}Basis/symbasis.new fort.4
cat >C1gnd.inp<</.
TITLE
ZnPc GND
SYMMETRY C1
CARTESIAN ANGSTROMS
[geometry]
END
[parameters]
FSYM scfocc
ALFA 147
BETA 147
FILE molden
END
[basis sets]
END
/.
StoBe.x <C1gnd.inp>& C1gnd.out
mv Molden.molf C1gnd.molden
rm fort.*
```

### Excited State (exc.run)

**Purpose**: Core-hole excited state calculation

**Outputs**:
- `.out`: Main output file
- `.molden`: Molden format orbitals

**Key differences**:
- Uses `FSYM scfocc excited`
- ALFA increased by 1
- Includes occupation strings

**Example**:
```bash
#!/bin/csh -f
ln -s ${STOBE}Basis/baslib.new7 fort.3
ln -s ${STOBE}Basis/symbasis.new fort.4
cat >C1exc.inp<</.
TITLE
ZnPc EXC
SYMMETRY C1
CARTESIAN ANGSTROMS
[geometry]
END
[parameters]
FSYM scfocc excited
ALFA 117
BETA 116
SYM 1
ALFA 0 1 14 0.0
BETA 0 0
END
MULLIKEN on full
FILE molden
END
[basis sets]
END
/.
StoBe.x <C1exc.inp>& C1exc.out
mv Molden.molf C1exc.molden
rm fort.*
```

### Transition Potential (tp.run)

**Purpose**: Transition potential for X-ray calculations

**Outputs**:
- `.out`: Main output file
- `.molden`: Molden format orbitals
- `.xas`: Transition potential data (from fort.11)

**Key differences**:
- Uses `XRAY xas` keyword
- Occupation uses 0.5 (half-electron)
- Moves `fort.11` to `.xas` file

**Example**:
```bash
#!/bin/csh -f
ln -s ${STOBE}Basis/baslib.new7 fort.3
ln -s ${STOBE}Basis/symbasis.new fort.4
cat >C1tp.inp<</.
TITLE
ZnPc TP
SYMMETRY C1
CARTESIAN ANGSTROMS
[geometry]
END
[parameters]
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
[basis sets]
END
/.
StoBe.x <C1tp.inp>& C1tp.out
mv Molden.molf C1tp.molden
mv fort.11 C1.xas
rm fort.*
```

### X-ray Spectrum (xas.run)

**Purpose**: Post-processing of transition potential results

**Input**: `.xas` file from transition potential calculation

**Outputs**:
- `.out`: Spectrum output
- `XrayT001.out`: Formatted spectrum file

**Uses**: `xrayspec.x` utility (not `StoBe.x`)

**Example**:
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

**X-ray input parameters**:
- `RANGE <E_min> <E_max>`: Energy range in eV
- `POINTS <n>`: Number of energy points
- `WIDTH <params>`: Broadening parameters

---

## File Naming Conventions

### Standard Patterns

- `<atom><number>gnd.run`: Ground state
- `<atom><number>exc.run`: Excited state
- `<atom><number>tp.run`: Transition potential
- `<atom><number>xas.run`: X-ray spectrum
- `<atom><number>seq.run`: Sequential batch script

**Examples**:
- `C1gnd.run`, `C1exc.run`, `C1tp.run`, `C1xas.run`
- `C2gnd.run`, `C2exc.run`, `C2tp.run`, `C2xas.run`

### Sequential Scripts

Sequential scripts run multiple calculations in order:

```bash
#!/bin/csh -f
chmod +x ./C1gnd.run
./C1gnd.run

chmod +x ./C1exc.run
./C1exc.run

chmod +x ./C1tp.run
./C1tp.run

chmod +x ./C1xas.run
./C1xas.run
```

---

## Advanced Patterns

### Directory Organization

```bash
#!/bin/csh -f
mkdir -p ../gnd
mkdir -p ../molden
ln -s ${STOBE}Basis/baslib.new7 fort.3
ln -s ${STOBE}Basis/symbasis.new fort.4
cat >C1gnd.inp<</.
[input content]
/.
StoBe.x <C1gnd.inp>& ../gnd/C1gnd.out
mv Molden.molf ../molden/C1gnd.molden
rm fort.*
```

### Conditional Execution

```bash
#!/bin/csh -f
if (! -f C1gnd.res) then
    echo "Running ground state..."
    ./C1gnd.run
else
    echo "Ground state already exists"
endif
```

### Error Handling

```bash
#!/bin/csh -f
set status = 0
ln -s ${STOBE}Basis/baslib.new7 fort.3
ln -s ${STOBE}Basis/symbasis.new fort.4
cat >C1gnd.inp<</.
[input content]
/.
StoBe.x <C1gnd.inp>& C1gnd.out
set status = $status
if ($status != 0) then
    echo "Error: Calculation failed"
    exit 1
endif
mv Molden.molf C1gnd.molden
rm fort.*
```

---

## Environment Variables

Common variables used:

- `${STOBE}`: StoBe installation directory
- `${PWD}`: Current working directory
- `${HOME}`: User home directory

**Setting STOBE**:
```bash
setenv STOBE /home/hduva/stobe/StoBeLinux/
```

**Using in script**:
```bash
ln -s ${STOBE}Basis/baslib.new7 fort.3
```

---

## Path Handling

### Absolute Paths

```bash
ln -s /home/hduva/stobe/StoBeLinux/Basis/baslib.new7 fort.3
```

### Relative Paths

```bash
# From Examples directory
ln -s ../Basis/baslib.new7 fort.3

# From calculation directory
ln -s ../../Basis/baslib.new7 fort.3
```

### Environment Variables

```bash
# Recommended approach
ln -s ${STOBE}Basis/baslib.new7 fort.3
```

---

## Output File Management

### Standard Outputs

```bash
# Main output
StoBe.x <input.inp>& output.out

# Molden file
mv Molden.molf output.molden

# X-ray data (TP calculations)
mv fort.11 output.xas

# Restart file (optional)
# mv fort.2 output.res
```

### Organized Outputs

```bash
# Create directories
mkdir -p ../gnd
mkdir -p ../exc
mkdir -p ../tp
mkdir -p ../nexafs
mkdir -p ../molden

# Move outputs
mv C1gnd.out ../gnd/
mv C1exc.out ../exc/
mv C1tp.out ../tp/
mv C1.out ../nexafs/
mv *.molden ../molden/
```

---

## Common Patterns

### Pattern 1: Single Calculation

```bash
#!/bin/csh -f
ln -s ${STOBE}Basis/baslib.new7 fort.3
ln -s ${STOBE}Basis/symbasis.new fort.4
cat >input.inp<</.
[input]
/.
StoBe.x <input.inp>& input.out
mv Molden.molf input.molden
rm fort.*
```

### Pattern 2: Multi-Step Workflow

```bash
#!/bin/csh -f
# Ground state
./C1gnd.run

# Excited state
./C1exc.run

# Transition potential
./C1tp.run

# X-ray spectrum
./C1xas.run
```

### Pattern 3: Batch Processing

```bash
#!/bin/csh -f
foreach atom (C1 C2 C3 C4 C5)
    cd $atom
    ./${atom}gnd.run
    ./${atom}exc.run
    ./${atom}tp.run
    ./${atom}xas.run
    cd ..
end
```

---

## Validation

Before running scripts, verify:

- [ ] Script is executable (`chmod +x script.run`)
- [ ] Basis set paths are correct
- [ ] StoBe executable path is correct
- [ ] Input file content is valid
- [ ] Output directory exists (if creating)
- [ ] Environment variables are set (if using)

---

## Troubleshooting

### Permission Denied

```bash
chmod +x script.run
```

### Basis Files Not Found

Check paths:
```bash
ls -l ${STOBE}Basis/baslib.new7
ls -l ${STOBE}Basis/symbasis.new
```

### StoBe Executable Not Found

Check path:
```bash
which StoBe.x
ls -l ${STOBE}Source/StoBe.x
```

### Input File Not Generated

Check heredoc syntax:
- Must use `<</.` to start
- Must use `/.` to end
- No spaces around delimiters

---

## Related References

- **Input File Format**: See `input-file-format.md` for .inp file structure
- **Calculation Workflows**: See `calculation-workflows.md` for multi-step patterns
- **Examples**: See `examples-catalog.md` for example run scripts
