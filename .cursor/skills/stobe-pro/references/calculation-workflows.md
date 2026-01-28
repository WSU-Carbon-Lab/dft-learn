# Calculation Workflows

> Reference for: StoBe Pro
> Load when: Setting up multi-step calculations, NEXAFS workflows, batch processing, understanding calculation dependencies

---

## Overview

StoBe calculations often require multiple steps executed in sequence. Understanding workflow patterns is essential for setting up complete calculation sequences, especially for X-ray spectroscopy applications.

---

## Standard NEXAFS Workflow

The most common workflow for X-ray absorption spectroscopy:

### Step 1: Ground State (gnd)

**Purpose**: Initial SCF calculation to obtain reference molecular orbitals

**Input**: Molecular geometry

**Output**: 
- Ground state orbitals
- Total energy
- `.out` file
- `.molden` file

**Key Parameters**:
```stobe
RUNTYPE startup nooptimize
FSYM scfocc
ALFA <n>
BETA <n>
```

**Run Script**: `*gnd.run`

### Step 2: Excited State (exc)

**Purpose**: Core-hole excited state calculation

**Input**: Ground state results (can use restart file)

**Output**:
- Core-hole state orbitals
- Excited state energy
- `.out` file
- `.molden` file

**Key Parameters**:
```stobe
RUNTYPE startup nooptimize
FSYM scfocc excited
ALFA <n+1>
BETA <n>
SYM 1
ALFA 0 1 <n_virt> 0.0
BETA 0 0
```

**Run Script**: `*exc.run`

**Dependencies**: Requires ground state calculation first

### Step 3: Transition Potential (tp)

**Purpose**: Calculate transition matrix elements for X-ray absorption

**Input**: Ground state geometry with half-core-hole

**Output**:
- Transition matrix elements
- `.out` file
- `.molden` file
- `.xas` file (fort.11)

**Key Parameters**:
```stobe
RUNTYPE startup nooptimize
FSYM scfocc excited
ALFA <n>
BETA <n>
SYM 1
ALFA 0 1 <n_virt> 0.5
BETA 0 0
XRAY xas
```

**Run Script**: `*tp.run`

**Dependencies**: Can use ground state geometry, but independent calculation

### Step 4: X-ray Spectrum (xas)

**Purpose**: Post-process transition potential results to generate spectrum

**Input**: `.xas` file from transition potential calculation

**Output**:
- Absorption spectrum
- `.out` file
- `XrayT001.out` file

**Key Parameters**:
```stobe
RANGE <E_min> <E_max>
POINTS <n_points>
WIDTH <width_params>
XRAY xas
```

**Run Script**: `*xas.run`

**Dependencies**: Requires transition potential calculation first

**Utility**: Uses `xrayspec.x` (not `StoBe.x`)

---

## Sequential Execution

### Sequential Script Pattern

Create `*seq.run` script:

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

### Execution Order

1. Ground state must run first
2. Excited state can use ground state results
3. Transition potential is independent but typically run after ground state
4. X-ray spectrum requires transition potential results

---

## Multi-Atom Workflows

For systems with multiple equivalent atoms:

### Pattern: One Calculation Per Atom

```bash
# Atom C1
./C1gnd.run
./C1exc.run
./C1tp.run
./C1xas.run

# Atom C2
./C2gnd.run
./C2exc.run
./C2tp.run
./C2xas.run

# Atom C3
./C3gnd.run
./C3exc.run
./C3tp.run
./C3xas.run
```

### Batch Processing

```bash
#!/bin/csh -f
foreach atom (C1 C2 C3 C4 C5 C6 C7 C8 C9 C10)
    cd $atom
    ./${atom}gnd.run
    ./${atom}exc.run
    ./${atom}tp.run
    ./${atom}xas.run
    cd ..
end
```

### Parallel Execution

Use job scheduler or parallel execution:
```bash
# Run all ground states in parallel
for atom in C1 C2 C3 C4 C5; do
    (cd $atom && ./${atom}gnd.run) &
done
wait
```

---

## Workflow Variations

### Variation 1: Geometry Optimization First

```bash
# Step 0: Optimize geometry
./molecule_opt.run

# Step 1: Ground state with optimized geometry
./molecule_gnd.run

# Step 2-4: Standard NEXAFS workflow
./molecule_exc.run
./molecule_tp.run
./molecule_xas.run
```

### Variation 2: Restart from Previous

```bash
# Step 1: Ground state
./C1gnd.run

# Step 2: Excited state using restart
# Modify run script to use:
# RUNTYPE restart
# RESINPUT C1gnd.res
./C1exc_restart.run
```

### Variation 3: Multiple Core Holes

```bash
# Calculate for different atoms
./C1tp.run  # Carbon 1
./C2tp.run  # Carbon 2
./N1tp.run  # Nitrogen 1
```

---

## Workflow Dependencies

### Dependency Graph

```
Ground State (gnd)
    │
    ├──→ Excited State (exc) [optional]
    │
    └──→ Transition Potential (tp)
            │
            └──→ X-ray Spectrum (xas)
```

### Required Dependencies

- **Excited state** → Requires ground state (can use restart)
- **X-ray spectrum** → Requires transition potential (needs .xas file)

### Optional Dependencies

- **Transition potential** → Can use ground state geometry but independent
- **Excited state** → Can run independently but typically uses ground state

---

## Workflow Organization

### Directory Structure

```
molecule/
├── C1/
│   ├── C1gnd.run
│   ├── C1exc.run
│   ├── C1tp.run
│   ├── C1xas.run
│   ├── C1seq.run
│   └── outputs/
├── C2/
│   └── ...
└── shared/
    └── geometry.xyz
```

### Output Organization

```
calculations/
├── gnd/
│   ├── C1gnd.out
│   ├── C2gnd.out
│   └── ...
├── exc/
│   ├── C1exc.out
│   └── ...
├── tp/
│   ├── C1tp.out
│   └── ...
├── nexafs/
│   ├── C1.out
│   └── ...
└── molden/
    ├── C1gnd.molden
    └── ...
```

---

## Workflow Automation

### Python Automation

```python
import subprocess
import os

def run_nexafs_workflow(atom_name: str, base_dir: str) -> None:
    """Run complete NEXAFS workflow for one atom."""
    atom_dir = os.path.join(base_dir, atom_name)
    
    # Ground state
    subprocess.run(['./', f'{atom_name}gnd.run'], cwd=atom_dir, check=True)
    
    # Excited state
    subprocess.run(['./', f'{atom_name}exc.run'], cwd=atom_dir, check=True)
    
    # Transition potential
    subprocess.run(['./', f'{atom_name}tp.run'], cwd=atom_dir, check=True)
    
    # X-ray spectrum
    subprocess.run(['./', f'{atom_name}xas.run'], cwd=atom_dir, check=True)

# Batch processing
for atom in ['C1', 'C2', 'C3']:
    run_nexafs_workflow(atom, './calculations')
```

### Shell Script Automation

```bash
#!/bin/csh -f
set atoms = (C1 C2 C3 C4 C5)

foreach atom ($atoms)
    echo "Processing $atom"
    cd $atom
    
    ./${atom}gnd.run
    if ($status != 0) then
        echo "Error in ground state for $atom"
        exit 1
    endif
    
    ./${atom}exc.run
    ./${atom}tp.run
    ./${atom}xas.run
    
    cd ..
end
```

---

## Error Handling in Workflows

### Checkpoint Pattern

```bash
#!/bin/csh -f
# Check if ground state already completed
if (-f C1gnd.out && -f C1gnd.molden) then
    echo "Ground state already exists, skipping"
else
    ./C1gnd.run
endif

# Check if excited state can proceed
if (-f C1gnd.molden) then
    ./C1exc.run
else
    echo "Error: Ground state not found"
    exit 1
endif
```

### Validation Checks

```bash
#!/bin/csh -f
# Run ground state
./C1gnd.run

# Validate output
grep -q "CONVERGED" C1gnd.out
if ($status != 0) then
    echo "Warning: Ground state did not converge"
    exit 1
endif

# Proceed with next step
./C1exc.run
```

---

## Workflow Best Practices

### 1. Validate Each Step

- Check convergence in output files
- Verify output files are created
- Validate file sizes are reasonable

### 2. Organize Outputs

- Use consistent naming
- Organize by calculation type
- Keep input files for reproducibility

### 3. Document Parameters

- Record calculation parameters
- Note any modifications
- Document dependencies

### 4. Handle Failures

- Check for errors at each step
- Don't proceed if previous step failed
- Log errors for debugging

### 5. Optimize Execution

- Run independent calculations in parallel
- Use job schedulers for large batches
- Monitor resource usage

---

## Common Workflow Patterns

### Pattern 1: Single Atom NEXAFS

```bash
./C1gnd.run
./C1exc.run
./C1tp.run
./C1xas.run
```

### Pattern 2: Multi-Atom Sequential

```bash
foreach atom (C1 C2 C3)
    ./${atom}gnd.run
    ./${atom}exc.run
    ./${atom}tp.run
    ./${atom}xas.run
end
```

### Pattern 3: Batch Ground States First

```bash
# Run all ground states
foreach atom (C1 C2 C3 C4 C5)
    ./${atom}gnd.run &
end
wait

# Then run excited states
foreach atom (C1 C2 C3 C4 C5)
    ./${atom}exc.run &
end
wait
```

### Pattern 4: Restart-Based

```bash
# Ground state
./C1gnd.run

# Excited state from restart
# Modify to use: RUNTYPE restart
./C1exc_restart.run
```

---

## Workflow Validation

Before running workflows, verify:

- [ ] All run scripts exist
- [ ] Dependencies are clear
- [ ] Output directories exist
- [ ] Basis set files accessible
- [ ] StoBe executable available
- [ ] Sufficient disk space
- [ ] Calculation order is correct

---

## Troubleshooting Workflows

### Step Failed

**Issue**: One step in workflow failed

**Fix**: 
1. Check error messages in output file
2. Verify input file is valid
3. Check dependencies are met
4. Fix issue and re-run failed step

### Missing Dependencies

**Issue**: Step requires output from previous step

**Fix**:
1. Verify previous step completed successfully
2. Check output files exist
3. Verify file names match
4. Re-run previous step if needed

### Out of Order Execution

**Issue**: Steps run in wrong order

**Fix**:
1. Use sequential script
2. Add dependency checks
3. Verify execution order

---

## Related References

- **X-ray Spectroscopy**: See `xray-spectroscopy.md` for X-ray workflow details
- **Run Scripts**: See `run-scripts.md` for script creation
- **Examples**: See `examples-catalog.md` for workflow examples
