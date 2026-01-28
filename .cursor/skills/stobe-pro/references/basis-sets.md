# Basis Sets

> Reference for: StoBe Pro
> Load when: Specifying basis sets, working with auxiliary/orbital/MCP/augmentation bases, understanding basis set formats

---

## Overview

StoBe uses four types of basis sets, each serving a specific purpose in the calculation. Understanding basis set selection and format is critical for accurate calculations.

---

## Basis Set Types

### 1. Auxiliary Basis Sets (A- prefix)

**Purpose**: Density fitting for exchange-correlation potential

**Format**: `A-<ELEMENT>(<modifier>) <specification>`

**Examples**:
```stobe
A-CARBON (5,2;5,2)
A-CARBON(+4) (3,3;3,3)
A-HYDROGEN (4,2;4,2)
A-NITROGEN (4,3;4,3)
A-ZINC (5,5;5,5)
```

**Specification Format**: `(n1,n2;n3,n4)`
- Numbers indicate contraction scheme
- Semicolon separates different angular momentum components

**Modifiers**:
- `(+4)`: Modified effective nuclear charge (for core-hole atoms)
- No modifier: Standard charge

### 2. Orbital Basis Sets (O- prefix)

**Purpose**: Molecular orbital expansion

**Format**: `O-<ELEMENT>(<modifier>) <specification>`

**Examples**:
```stobe
O-CARBON iii_iglo
O-CARBON(+4) (311/211/1)
O-HYDROGEN (311/1) misc
O-NITROGEN (33/3)
O-ZINC (63321/531/311)
```

**Specification Types**:

1. **Predefined names**:
   - `iii_iglo`: IGLO-III basis set
   - Other predefined names available

2. **Explicit notation**: `(n1/n2/n3)`
   - Numbers indicate contraction
   - Slashes separate different shells
   - Example: `(311/211/1)` means 3s1p1 / 2s1p1 / 1s

**Modifiers**:
- `(+4)`: Modified charge for core-hole atoms
- No modifier: Standard charge

### 3. Model Core Potential (P- prefix)

**Purpose**: Replace core electrons for heavy elements or core-hole calculations

**Format**: `P-<ELEMENT>(<modifier>) <specification>`

**Examples**:
```stobe
P-CARBON(+4) (3,1:8,0)
```

**Specification Format**: `(n1,n2:n3,n4)`
- Colon (`:`) separates different components
- Used for heavy elements or core-hole calculations

**When to Use**:
- Heavy elements (transition metals, etc.)
- Core-hole calculations (reduces computational cost)
- Large systems where core electrons can be approximated

### 4. Augmentation Basis Sets (X- prefix)

**Purpose**: Additional basis functions for X-ray calculations

**Format**: `X-<option>`

**Options**:
- `X-FIRST`: Apply augmentation to first atom (core-hole site)
- `X-DUMMY`: No augmentation for this atom

**Usage**:
- Only for X-ray absorption calculations
- One `X-FIRST` per calculation (for core-hole atom)
- All other atoms use `X-DUMMY`

**Example**:
```stobe
X-DUMMY
X-DUMMY
X-FIRST
X-DUMMY
...
```

---

## Basis Set Library (baslib.new7)

### File Structure

- Binary/compressed format
- Contains predefined basis sets for all elements
- Referenced via `fort.3` symbolic link in run scripts
- Basis sets accessed by element name and specification

### Accessing Basis Sets

Basis sets are referenced by:
1. Element name (e.g., `CARBON`, `HYDROGEN`, `NITROGEN`)
2. Optional modifiers (e.g., `(+4)` for modified charge)
3. Specification string (e.g., `(5,2;5,2)`, `iii_iglo`)

### Common Basis Sets

#### Carbon

**Auxiliary**:
- `A-CARBON (5,2;5,2)`: Standard auxiliary basis
- `A-CARBON(+4) (3,3;3,3)`: Modified charge auxiliary basis

**Orbital**:
- `O-CARBON iii_iglo`: IGLO-III orbital basis
- `O-CARBON(+4) (311/211/1)`: Modified charge orbital basis

**MCP**:
- `P-CARBON(+4) (3,1:8,0)`: Model core potential for core-hole

#### Hydrogen

**Auxiliary**:
- `A-HYDROGEN (4,2;4,2)`: Standard auxiliary basis
- `A-HYDROGEN (3,1;3,1)`: Smaller auxiliary basis

**Orbital**:
- `O-HYDROGEN (311/1) misc`: Standard orbital basis

#### Nitrogen

**Auxiliary**:
- `A-NITROGEN (4,3;4,3)`: Standard auxiliary basis
- `A-NITROGEN (5,2;5,2)`: Larger auxiliary basis

**Orbital**:
- `O-NITROGEN (33/3)`: Standard orbital basis
- `O-NITROGEN (7111/411/1)`: Larger orbital basis

#### Oxygen

**Auxiliary**:
- `A-OXYGEN (4,3;4,3)`: Standard auxiliary basis

**Orbital**:
- `O-OXYGEN (33/3)`: Standard orbital basis

#### Zinc

**Auxiliary**:
- `A-ZINC (5,5;5,5)`: Standard auxiliary basis

**Orbital**:
- `O-ZINC (63321/531/311)`: Standard orbital basis

---

## Basis Set Selection Rules

### Rule 1: One Basis Set Per Atom

- Each atom in geometry must have exactly one basis set of each type
- Order must match geometry order exactly
- Missing or extra basis sets cause errors

### Rule 2: Element Matching

- Basis set element must match atom element
- `A-CARBON` for carbon atoms
- `A-HYDROGEN` for hydrogen atoms
- Cannot mix elements

### Rule 3: Modified Charge Basis Sets

For core-hole calculations:
- Core-hole atom: Use modified charge basis sets `(+4)`
- Other atoms of same element: Use standard basis sets
- Example: First carbon uses `A-CARBON (5,2;5,2)`, others use `A-CARBON(+4) (3,3;3,3)`

### Rule 4: MCP Usage

- Use MCP for heavy elements (Zn, transition metals, etc.)
- Use MCP for core-hole atoms to reduce cost
- One MCP basis set per atom that needs it
- Not required for light elements (H, C, N, O)

### Rule 5: Augmentation Basis

- Only for X-ray calculations
- One `X-FIRST` for core-hole atom
- `X-DUMMY` for all other atoms
- Must match atom count

---

## Basis Set Order

Basis sets appear in this order in input files:

1. **Auxiliary basis sets** (A- prefix)
   - One per atom, matching geometry order
   - Required for all atoms

2. **Orbital basis sets** (O- prefix)
   - One per atom, matching geometry order
   - Required for all atoms

3. **Model core potential** (P- prefix)
   - One per atom that needs MCP
   - Optional, only for certain atoms

4. **Augmentation basis** (X- prefix)
   - One per atom
   - Optional, only for X-ray calculations

**Example**:
```stobe
A-CARBON (5,2;5,2)
A-CARBON(+4) (3,3;3,3)
A-HYDROGEN (4,2;4,2)
O-CARBON iii_iglo
O-CARBON(+4) (311/211/1)
O-HYDROGEN (311/1) misc
P-CARBON(+4) (3,1:8,0)
P-CARBON(+4) (3,1:8,0)
X-FIRST
X-DUMMY
X-DUMMY
```

---

## Basis Set Specifications

### Auxiliary Basis Format

Format: `(n1,n2;n3,n4)`

- `n1, n2`: First angular momentum component
- `n3, n4`: Second angular momentum component
- Semicolon separates components

**Examples**:
- `(5,2;5,2)`: Standard carbon auxiliary
- `(3,3;3,3)`: Modified charge carbon auxiliary
- `(4,2;4,2)`: Hydrogen auxiliary
- `(4,3;4,3)`: Nitrogen/oxygen auxiliary

### Orbital Basis Format

Two formats:

1. **Predefined names**:
   - `iii_iglo`: IGLO-III basis
   - Other predefined names available

2. **Explicit notation**: `(n1/n2/n3)`
   - `n1`: Contraction for first shell
   - `n2`: Contraction for second shell
   - `n3`: Contraction for third shell
   - Slashes separate shells

**Examples**:
- `(311/211/1)`: 3s1p1 / 2s1p1 / 1s
- `(311/1)`: 3s1p1 / 1s (hydrogen)
- `(33/3)`: 3s3p / 3s (nitrogen/oxygen)
- `(63321/531/311)`: Large basis (zinc)

### MCP Basis Format

Format: `(n1,n2:n3,n4)`

- Colon (`:`) separates components
- Used for model core potential

**Example**:
- `(3,1:8,0)`: Carbon MCP specification

---

## Core-Hole Basis Set Pattern

For X-ray absorption calculations:

### Core-Hole Atom

- **Auxiliary**: `A-CARBON (5,2;5,2)` (standard, not modified)
- **Orbital**: `O-CARBON iii_iglo` (standard, not modified)
- **MCP**: `P-CARBON(+4) (3,1:8,0)` (modified charge)
- **Augmentation**: `X-FIRST`

### Other Atoms of Same Element

- **Auxiliary**: `A-CARBON(+4) (3,3;3,3)` (modified charge)
- **Orbital**: `O-CARBON(+4) (311/211/1)` (modified charge)
- **MCP**: `P-CARBON(+4) (3,1:8,0)` (modified charge)
- **Augmentation**: `X-DUMMY`

**Key Point**: Core-hole atom uses standard basis sets, other atoms use modified charge basis sets.

---

## Validation

When specifying basis sets, verify:

- [ ] Count matches atom count exactly
- [ ] Order matches geometry order
- [ ] Elements match atom elements
- [ ] Modified charge basis sets used correctly
- [ ] MCP specified for atoms that need it
- [ ] Augmentation basis matches calculation type
- [ ] Syntax is correct (parentheses, semicolons, slashes)

---

## Common Errors

### Mismatched Count

**Error**: N atoms but M basis sets

**Fix**: Ensure exactly one basis set per atom for each type

### Wrong Element

**Error**: `A-CARBON` for hydrogen atom

**Fix**: Match basis set element to atom element

### Invalid Format

**Error**: Syntax error in basis specification

**Fix**: Check parentheses, semicolons, slashes match format

### Missing MCP

**Error**: Heavy element without MCP

**Fix**: Add MCP basis set for heavy elements

### Wrong Augmentation

**Error**: Multiple `X-FIRST` or missing `X-DUMMY`

**Fix**: One `X-FIRST` for core-hole atom, `X-DUMMY` for others

---

## Basis Set Selection Guide

### For Ground State Calculations

- Use standard basis sets for all atoms
- No modified charge basis sets
- No MCP unless heavy elements
- No augmentation basis

### For Core-Hole Calculations

- Core-hole atom: Standard basis sets
- Other atoms: Modified charge basis sets `(+4)`
- Use MCP for core-hole atom
- Use augmentation: `X-FIRST` for core-hole, `X-DUMMY` for others

### For Heavy Elements

- Use MCP to reduce computational cost
- Standard auxiliary and orbital basis sets
- Consider larger basis sets for accuracy

---

## Related References

- **Input File Format**: See `input-file-format.md` for basis set section structure
- **X-ray Spectroscopy**: See `xray-spectroscopy.md` for core-hole basis set patterns
- **Examples**: See `examples-catalog.md` for basis set examples
