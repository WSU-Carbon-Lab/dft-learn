# Output Files

> Reference for: StoBe Pro
> Load when: Parsing StoBe output files, extracting results, validating calculations, checking convergence

---

## Overview

StoBe output files (.out) contain extensive information about calculations including geometry, SCF iterations, final energies, orbital information, and property calculations. Understanding output file structure is essential for result extraction and validation.

---

## Output File Structure

### Key Sections

1. **Header Information**
   - Program version and date
   - Input file reference
   - Calculation type

2. **Geometry Information**
   - Input geometry
   - Final geometry (if optimized)
   - Symmetry information

3. **SCF Iteration Information**
   - Iteration-by-iteration energy
   - Convergence criteria
   - Density changes

4. **Final Results**
   - Total energy
   - Orbital energies
   - Convergence status

5. **Property Output**
   - Mulliken populations (if requested)
   - Orbital information
   - Transition moments (for X-ray calculations)

---

## Parsing Output Files

### Key Information to Extract

#### Convergence Status

Look for:
```
CONVERGED
```

or

```
NOT CONVERGED
```

#### Final Energy

Pattern:
```
TOTAL ENERGY = <value>
```

or

```
FINAL ENERGY = <value>
```

#### SCF Iterations

Pattern:
```
ITERATION    ENERGY         DELTA E        DELTA D
    1    -1234.567890    0.00000000    0.12345678
    2    -1234.678901   -0.11101100    0.01234567
    ...
```

#### Orbital Energies

Pattern:
```
ORBITAL ENERGIES
   1    -10.1234
   2     -9.8765
   ...
```

#### Mulliken Populations

Pattern:
```
MULLIKEN POPULATION ANALYSIS
ATOM    CHARGE    POPULATION
  C1     0.123     4.877
  C2     0.456     5.544
  ...
```

---

## Parsing Strategy

### Python Parsing Example

```python
from typing import Dict, List, Optional
import re

def parse_stobe_output(filename: str) -> Dict:
    """Parse StoBe output file and extract key information."""
    results = {
        'converged': False,
        'energy': None,
        'iterations': [],
        'orbitals': [],
        'populations': {},
        'calculation_type': None
    }
    
    with open(filename, 'r') as f:
        content = f.read()
        lines = content.split('\n')
        
        # Check convergence
        if 'CONVERGED' in content:
            results['converged'] = True
        
        # Extract final energy
        energy_match = re.search(r'TOTAL ENERGY\s*=\s*([-\d.]+)', content)
        if energy_match:
            results['energy'] = float(energy_match.group(1))
        
        # Extract SCF iterations
        in_iterations = False
        for i, line in enumerate(lines):
            if 'ITERATION' in line and 'ENERGY' in line:
                in_iterations = True
                continue
            if in_iterations:
                if line.strip() and not line.startswith(' '):
                    in_iterations = False
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        iteration = int(parts[0])
                        energy = float(parts[1])
                        results['iterations'].append({
                            'iteration': iteration,
                            'energy': energy
                        })
                    except (ValueError, IndexError):
                        continue
        
        # Extract orbital energies
        in_orbitals = False
        for line in lines:
            if 'ORBITAL ENERGIES' in line:
                in_orbitals = True
                continue
            if in_orbitals:
                if line.strip() and not line[0].isdigit():
                    if 'MULLIKEN' in line or 'END' in line:
                        in_orbitals = False
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        orbital_num = int(parts[0])
                        energy = float(parts[1])
                        results['orbitals'].append({
                            'orbital': orbital_num,
                            'energy': energy
                        })
                    except (ValueError, IndexError):
                        continue
        
        # Extract Mulliken populations
        in_mulliken = False
        for line in lines:
            if 'MULLIKEN' in line and 'POPULATION' in line:
                in_mulliken = True
                continue
            if in_mulliken:
                if 'ATOM' in line and 'CHARGE' in line:
                    continue
                if not line.strip() or line.startswith('---'):
                    continue
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        atom = parts[0]
                        charge = float(parts[1])
                        population = float(parts[2])
                        results['populations'][atom] = {
                            'charge': charge,
                            'population': population
                        }
                    except (ValueError, IndexError):
                        continue
    
    return results
```

---

## Restart Files

### Format

Restart files (`.res`) are binary files containing:
- Molecular geometry
- Molecular orbitals
- Density matrices
- SCF convergence information

### Usage

Restart files enable:
- Continuing interrupted calculations
- Starting new calculations from previous results
- Transferring information between calculations

### Conversion Utilities

**Binary to ASCII**:
```bash
transf input.res output.asc
```

**Symmetry reduction**:
```bash
desymm input.res output.res C2V C1
```

**Analysis**:
```bash
anlyz input.res
```

### Restart File Usage

In input file:
```stobe
RUNTYPE restart
RESINPUT previous.res
```

---

## Molden Files

### Format

Molden format (`.molden`) contains:
- Molecular geometry
- Molecular orbitals
- Orbital energies
- Basis set information

### Usage

Molden files can be:
- Viewed with Molden, VMD, etc.
- Used for visualization
- Analyzed with other tools
- Converted to other formats

### Generation

Include in input:
```stobe
FILE molden
```

Output file: `Molden.molf` (renamed to `.molden` in run script)

---

## X-ray Output Files

### Transition Potential Output

File: `fort.11` (moved to `.xas`)

Contains:
- Transition energies
- Transition moments
- Orbital information

### Spectrum Output

File: `XrayT001.out`

Contains:
- Energy vs. intensity
- Broadened spectrum
- Can be plotted directly

Format:
```
ENERGY    INTENSITY
280.0     0.001234
280.1     0.001456
...
```

---

## Validation

### Convergence Check

```python
def check_convergence(output_file: str) -> bool:
    """Check if calculation converged."""
    with open(output_file, 'r') as f:
        content = f.read()
        return 'CONVERGED' in content
```

### Energy Extraction

```python
def extract_energy(output_file: str) -> Optional[float]:
    """Extract final energy from output."""
    import re
    with open(output_file, 'r') as f:
        content = f.read()
        match = re.search(r'TOTAL ENERGY\s*=\s*([-\d.]+)', content)
        if match:
            return float(match.group(1))
    return None
```

### Iteration Analysis

```python
def analyze_iterations(output_file: str) -> Dict:
    """Analyze SCF iteration convergence."""
    iterations = parse_stobe_output(output_file)['iterations']
    
    if not iterations:
        return {'error': 'No iterations found'}
    
    energies = [it['energy'] for it in iterations]
    energy_changes = [energies[i+1] - energies[i] 
                     for i in range(len(energies)-1)]
    
    return {
        'total_iterations': len(iterations),
        'final_energy': energies[-1],
        'energy_change': energy_changes[-1] if energy_changes else 0,
        'converged': abs(energy_changes[-1]) < 0.001 if energy_changes else False
    }
```

---

## Common Output Patterns

### Successful Calculation

```
*******************************************************************************
 StoBe-deMon SOFTWARE
 Stockholm-Berlin version 3.3
*******************************************************************************

TITLE
Test Calculation

...

ITERATION    ENERGY         DELTA E        DELTA D
    1    -1234.567890    0.00000000    0.12345678
    2    -1234.678901   -0.11101100    0.01234567
    ...
   10    -1234.789012   -0.00000123    0.00000012

CONVERGED

TOTAL ENERGY = -1234.78901234
```

### Failed Convergence

```
...

ITERATION    ENERGY         DELTA E        DELTA D
    1    -1234.567890    0.00000000    0.12345678
    ...
  300    -1234.789012   -0.00123456    0.01234567

NOT CONVERGED AFTER 300 ITERATIONS
```

---

## Error Detection

### Common Error Patterns

**Convergence failure**:
```
NOT CONVERGED
```

**SCF error**:
```
ERROR IN SCF
```

**Basis set error**:
```
BASIS SET ERROR
```

**Symmetry error**:
```
SYMMETRY ERROR
```

**Memory error**:
```
INSUFFICIENT MEMORY
```

---

## Output File Validation

When parsing output files, check:

- [ ] File exists and is readable
- [ ] Contains expected sections
- [ ] Convergence status present
- [ ] Final energy extracted
- [ ] Iteration count reasonable
- [ ] No error messages
- [ ] File size is reasonable

---

## Related References

- **Input File Format**: See `input-file-format.md` for input structure
- **Calculation Workflows**: See `calculation-workflows.md` for workflow validation
- **Examples**: See `examples-catalog.md` for output file examples
