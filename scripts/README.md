# StoBe DFT Calculation Suite

A complete suite of modern command-line tools for StoBe DFT calculations, featuring input file generation and parallel calculation scheduling with beautiful terminal interfaces.

## 🎯 Overview

This suite provides two complementary CLI tools:

1. **📝 Input Generator (`generate.py`)** - Generate StoBe DFT calculation input files from molecular configurations
2. **🚀 Calculation Scheduler (`schedule.py`)** - Execute StoBe calculations in parallel with real-time progress tracking

Both tools are built with modern Python practices, featuring automatic dependency management, rich terminal output, and robust error handling.

---

# 📝 StoBe Input Generator CLI

Generate StoBe DFT calculation input files from molecular configurations. This tool replaces the original StoBe Input Generator with a flexible, user-friendly CLI that supports auto-detection, dynamic configuration loading, and proper file organization.

## 🚀 Features

### Input Generator Features

- **UV Script Format**: Self-contained script with automatic dependency management
- **Auto-Detection**: Automatically detects XYZ files in the target directory
- **Dynamic Configuration**: Loads `molConfig.py` from target directories without hard-coded paths
- **Robust File Organization**: Fixed bug where C10+ files were incorrectly placed in C1 directory
- **Rich Terminal Output**: Beautiful progress display with colors and status updates
- **Flexible Input Methods**: Multiple ways to specify XYZ files and directories
- **Error Handling**: Comprehensive error checking with clear, colored messages

## 📦 Installation

No installation required! The script uses UV's inline dependencies.

### Prerequisites

- Python 3.10+
- [UV package manager](https://docs.astral.sh/uv/)
- StoBe environment setup

Install UV if you haven't already:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

## 🔧 Usage

### Input Generator Usage

Generate input files for your StoBe calculations:

```bash
# Auto-detect XYZ file in the directory
uv run generate.py tcta_flat

# Specify XYZ file explicitly
uv run generate.py tcta_flat tcta_flat.xyz

# Use --xyz flag for external files
uv run generate.py tcta_flat --xyz /path/to/geometry.xyz

# Verbose output for debugging
uv run generate.py tcta_flat --verbose
```

### Directory Structure

Your calculation directory should contain:

```
my_calculation/
├── molConfig.py          # Configuration file
└── geometry.xyz          # Molecular geometry file (auto-detected)
```

**Note**: The `molConfig.py` file no longer needs `fname` or `mname` variables - these are auto-detected!

## ⚙️ Configuration

### Updated molConfig.py Format

The CLI loads configuration dynamically, eliminating the need for hard-coded paths:

```python
# *****StoBe Input Generator Configuration****
aname = "C"              # Element name for calculations
nFiles = 18              # Number of atoms for calculations
difElems = 3             # Number of different elements

# Effective Nuclear Charges
element1a = 6            # Non-screened effective nuclear charge
element1b = 4            # Excitation state screened charge
element2 = 1             # Hydrogen
element3 = 7             # Nitrogen

# Element counts
nElem1 = 54              # Carbon atoms
nElem2 = 36              # Hydrogen atoms
nElem3 = 4               # Nitrogen atoms

# Basis sets
elem1_Abasis_a = "A-CARBON (5,2;5,2)"
elem1_Abasis_b = "A-CARBON(+4) (3,3;3,3)"
elem2_Abasis = "A-HYDROGEN (2;2)"
elem3_Abasis = "A-NITROGEN (5,2;5,2)"

elem1_Obasis_a = "O-CARBON (4,2,1;4,2,1)"
elem1_Obasis_b = "O-CARBON(+4) (4,2,1;4,2,1)"
elem2_Obasis = "O-HYDROGEN (2,1;2,1)"
elem3_Obasis = "O-NITROGEN (4,2,1;4,2,1)"

# Electron configuration
alpha = 194              # Spin up electrons
beta = 194               # Spin down electrons
alfaOcc = "0 1 5 0.0"    # Alpha occupation for excited state
betaOcc = "0 0"          # Beta occupation
alfaOccTP = "0 1 5 0.5"  # Alpha occupation for transition potential

# Optional parameters (defaults will be used if not specified)
title = "My Calculation"
multiplicity = "1"
charge = "0"
# ... (other StoBe parameters)
```

## 📁 Output Structure

The generator creates a well-organized directory structure:

```
tcta_flat/
├── molConfig.py
├── tcta_flat.xyz
├── C1/
│   ├── C1gnd.run       # Ground state calculation
│   ├── C1exc.run       # Excited state calculation
│   ├── C1tp.run        # Transition potential calculation
│   ├── C1xas.run       # XAS calculation
│   └── C1seq.run       # Sequential batch script
├── C2/
│   ├── C2gnd.run
│   ├── C2exc.run
│   ├── C2tp.run
│   ├── C2xas.run
│   └── C2seq.run
...
└── C18/
    ├── C18gnd.run
    ├── C18exc.run
    ├── C18tp.run
    ├── C18xas.run
    └── C18seq.run
```

### Generated File Types

- **`*.gnd.run`** - Ground state calculations
- **`*.exc.run`** - Excited state calculations
- **`*.tp.run`** - Transition potential calculations
- **`*.xas.run`** - X-ray absorption spectroscopy calculations
- **`*.seq.run`** - Sequential batch calculations (runs all types in order)

## 🐛 Bug Fixes

### Fixed: C10+ File Organization Bug

**Issue**: Run files for atoms with indices > 9 (C10, C11, C12, etc.) were incorrectly placed in the C1 directory due to string prefix matching.

**Root Cause**: `"C10gnd.run".startswith("C1")` returned `True`, causing misplacement.

**Solution**: Implemented precise regex-based matching:

```python
# Before (buggy):
if run_file.startswith(atom_dir):

# After (fixed):
pattern = rf"^{re.escape(aname)}{i}(?:\D|$)"
if re.match(pattern, run_file):
```

**Result**: All files now correctly organized in their respective directories:

- ✅ C1 directory contains only C1 files
- ✅ C10 directory contains only C10 files
- ✅ C15 directory contains only C15 files
- ✅ Works correctly for any number of atoms

## 📋 CLI Reference

```
Usage: generate.py [OPTIONS] RUN_DIRECTORY [XYZ_FILE]

Generate StoBe DFT calculation input files.

Arguments:
  RUN_DIRECTORY  Directory containing molConfig.py [required]
  XYZ_FILE       Path to XYZ geometry file [optional]

Options:
  --xyz TEXT     Alternative way to specify XYZ file
  --verbose -v   Verbose output with detailed information
  --help         Show this message and exit
```

### Examples

#### Example 1: Simple Auto-Detection

```bash
uv run generate.py tcta_flat
```

**Output**:

```
🧪 StoBe DFT Input Generator
Auto-detected XYZ file: /path/to/tcta_flat/tcta_flat.xyz
Loading configuration from /path/to/tcta_flat/molConfig.py
✅ Successfully generated StoBe input files for 18 C atoms
Output directory: /path/to/tcta_flat
```

#### Example 2: Explicit XYZ File

```bash
uv run generate.py my_calc my_geometry.xyz
```

#### Example 3: External XYZ File

```bash
uv run generate.py my_calc --xyz /home/user/geometries/molecule.xyz
```

#### Example 4: Verbose Mode

```bash
uv run generate.py tcta_flat --verbose
```

**Additional Output**:

```
XYZ file: /path/to/tcta_flat/tcta_flat.xyz
Run directory: /path/to/tcta_flat
Using existing tcta_flat.xyz file
```

## 🔄 Migration from Original

### Quick Migration Steps

1. **Remove hard-coded paths** from your `molConfig.py`:

   ```python
   # Remove these lines:
   # fname = "path/to/file"
   # mname = "molecule_name"
   ```

2. **Place files in calculation directory**:

   ```bash
   my_calculation/
   ├── molConfig.py    # Updated configuration
   └── geometry.xyz    # Your molecular geometry
   ```

3. **Run the new CLI**:
   ```bash
   uv run generate.py my_calculation
   ```

### Comparison: Original vs CLI

| Feature               | Original Script           | New CLI                      |
| --------------------- | ------------------------- | ---------------------------- |
| **Dependencies**      | Import molConfig directly | Dynamic loading              |
| **File paths**        | Hard-coded fname/mname    | Auto-detected                |
| **Interface**         | Manual script execution   | Modern CLI with help         |
| **Progress**          | Basic text output         | Rich progress with colors    |
| **Error handling**    | Minimal                   | Comprehensive with colors    |
| **Flexibility**       | Single configuration      | Multiple directories         |
| **File organization** | String prefix (buggy)     | Regex-based (fixed)          |
| **XYZ input**         | Hard-coded path           | Auto-detect + flexible input |

## 🔬 StoBe Integration

### Environment Requirements

Ensure these StoBe environment variables are set:

```bash
export STOBE_HOME="/path/to/stobe"
export PATH="$STOBE_HOME/Source:$PATH"
```

### Supported Elements

The CLI includes energy ranges for XAS calculations:

- **Carbon (C)**: 280-320 eV
- **Nitrogen (N)**: 400-450 eV
- **Oxygen (O)**: 530-580 eV
- **Nickel (Ni)**: 830-880 eV

Additional elements can be added by extending the `energy_ranges` dictionary in the source code.

## 🚀 Advanced Usage

### Working with Multiple Molecules

```bash
# Generate files for different molecules
uv run generate.py tcta_flat
uv run generate.py basis/znpc_default
uv run generate.py PEO/bare
uv run generate.py nio
```

### Batch Processing

```bash
#!/bin/bash
# Generate input files for multiple calculations
for dir in tcta_flat basis/znpc_* PEO/*; do
    if [ -f "$dir/molConfig.py" ]; then
        echo "Processing $dir..."
        uv run generate.py "$dir"
    fi
done
```

### Custom XYZ Files

```bash
# Use different geometry for same configuration
uv run generate.py tcta_flat --xyz tcta_optimized.xyz
uv run generate.py tcta_flat --xyz tcta_distorted.xyz
```

---

# 🛠️ Development & Technical Details

## Requirements

- Python 3.10+
- [UV package manager](https://docs.astral.sh/uv/)
- StoBe environment setup

### Auto-Managed Dependencies

Both tools automatically manage their dependencies via UV:

- **Input Generator**: `typer`, `rich`, `natsort`
- **Scheduler**: `typer`, `rich`

## Code Structure

### Input Generator (`generate.py`)

```python
# Key components:
- load_mol_config()     # Dynamic configuration loading
- auto_detect_xyz()     # XYZ file auto-detection
- makeRunFile()         # Generate calculation input files
- listFiles()           # Organize files with fixed regex matching
- makeXASrun()          # Generate XAS-specific runs
- makeSEQrun()          # Generate sequential batch scripts
```

### Scheduler (`schedule.py`)

```python
# Key components:
- explore_command()     # Directory exploration and analysis
- run_calculation()     # Parallel calculation execution
- progress_tracking()   # Real-time progress display
- worker_management()   # Configurable parallel workers
```

## 🚀 Advanced Workflows

### Batch Processing Multiple Molecules

```bash
#!/bin/bash
# Complete workflow for multiple molecules
for dir in tcta_flat basis/znpc_* PEO/*; do
    if [ -f "$dir/molConfig.py" ]; then
        echo "🔧 Generating inputs for $dir..."
        uv run generate.py "$dir"

        echo "🚀 Running calculations for $dir..."
        uv run schedule.py all "$dir"
    fi
done
```

### Custom Calculation Sequences

```bash
# Generate inputs
uv run generate.py my_molecule

# Run ground states first (parallel)
uv run schedule.py --workers 8 gnd my_molecule

# Wait for ground states, then run excited states
uv run schedule.py --workers 6 exc my_molecule

# Finally run XAS calculations
uv run schedule.py --workers 4 xas my_molecule
```

### Testing and Validation

```bash
# Quick test on single atom
uv run generate.py test_molecule
uv run schedule.py gnd test_molecule --atom C1 --verbose

# Systematic testing
uv run schedule.py explore my_molecule  # Check what's available
uv run schedule.py --workers 2 gnd my_molecule --atom C1 --atom C2
```

## 📝 License

This tool is part of the StoBe DFT calculation suite. Please refer to your StoBe license for usage terms.

## 🤝 Contributing

Improvements and bug reports are welcome! Key areas for enhancement:

- Additional element support for XAS calculations
- Extended basis set libraries
- Enhanced error diagnostics
- Additional output formats

## 📚 Related Tools

- **StoBe Scheduler CLI**: For parallel execution of generated calculations
- **StoBe Analysis Tools**: For processing calculation results
- **Molecular Visualization**: For geometry verification

---

# 🚀 StoBe Calculation Scheduler CLI

A command-line interface for scheduling and running StoBe DFT calculations with parallel execution and **real-time progress bars**.

## 🎯 Scheduler Features

- **Parallel Execution**: Run multiple calculations simultaneously with configurable worker limits
- **Rich Progress Bars**: Live updating progress bars with animated spinners using Typer/Rich
- **Flexible Targeting**: Run calculations on entire directories or specific atoms
- **Multiple Calculation Types**: Support for ground state (gnd), excited state (exc), transition potential (tp), XAS, and sequential calculations
- **Beautiful Terminal Output**: Professional-grade progress display with colors and animations
- **UV Script**: Self-contained script with automatic dependency management

## 🔧 Scheduler Usage

### Basic Scheduler Commands

All commands follow the pattern: `uv run schedule.py [command] [directory] [options]`

#### 1. Explore a Directory

See what calculations are available in a directory:

```bash
uv run schedule.py explore tcta_flat
uv run schedule.py explore .
```

This shows available atoms, calculation types, and total run files.

#### 2. Run Ground State Calculations

Run all ground state calculations in a directory:

```bash
# Current directory
uv run schedule.py gnd .

# Specific directory
uv run schedule.py gnd /path/to/directory

# Specific atom only
uv run schedule.py gnd tcta_flat --atom C1

# Multiple specific atoms
uv run schedule.py gnd tcta_flat --atom C1 --atom C2 --atom C3
```

#### 3. Run Other Calculation Types

```bash
# Excited state calculations
uv run schedule.py exc tcta_flat

# Transition potential calculations
uv run schedule.py tp tcta_flat

# XAS calculations
uv run schedule.py xas tcta_flat

# Sequential calculations
uv run schedule.py seq tcta_flat
```

#### 4. Run All Calculation Types

Run ground state → excited state → transition potential → XAS sequentially:

```bash
uv run schedule.py all tcta_flat
uv run schedule.py all .
```

Each calculation type runs in parallel within itself, but the types run sequentially.

### Advanced Scheduler Options

#### Parallel Workers

Control the number of parallel workers (default: 4):

```bash
# Use 8 parallel workers
uv run schedule.py --workers 8 gnd tcta_flat

# Use 2 parallel workers for slower systems
uv run schedule.py --workers 2 all .
```

#### Verbose Output

Enable verbose output for debugging:

```bash
uv run schedule.py --verbose gnd tcta_flat
```

### Real-time Progress Display

During execution, you'll see beautiful progress bars like this:

```
🎯 Starting ground state calculations in tcta_flat
⚛️  Target atoms: ['C1', 'C2', 'C3']

⠹ Overall Progress (1 types) ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  60% 0:02:15
⠙ C1: ✅ Completed (12.3s)    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:12
⠹ C2: Running... 8s          ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  75% 0:00:08
⠋ C3: Pending...             ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━   0% 0:00:00
```

Key features of the progress display:

- **Animated Spinners**: Each calculation gets a rotating spinner (⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏)
- **Individual Progress Bars**: Real-time progress for each directory/atom
- **Overall Progress**: Main progress bar showing total completion
- **Live Timing**: Elapsed time for each calculation and overall progress
- **Status Updates**: Clear indication of completed (✅), running, pending, and failed calculations

## 📋 Scheduler Examples

### Quick Testing

```bash
# Explore what's available
uv run schedule.py explore tcta_flat

# Run just one atom for testing
uv run schedule.py gnd tcta_flat --atom C1

# Run a few atoms with more workers
uv run schedule.py --workers 6 gnd tcta_flat --atom C1 --atom C2
```

### Production Runs

```bash
# Run all ground states with maximum parallelism
uv run schedule.py --workers 16 gnd tcta_flat

# Run complete calculation sequence
uv run schedule.py --workers 8 all tcta_flat

# Run specific calculation type on all atoms
uv run schedule.py --workers 12 exc basis/znpc_default
```

### Working with Different Directories

```bash
# Current directory
uv run schedule.py gnd .

# Relative path
uv run schedule.py gnd ../other_project

# Absolute path
uv run schedule.py gnd /home/user/calculations/my_molecule

# Basis directory
uv run schedule.py gnd basis/znpc_321
```

## 📊 Scheduler Output and Results

- **Live Progress**: Real-time updates every 2 seconds during execution
- **Completion Summary**: Final statistics showing completed, failed, and timing information
- **Error Reporting**: Clear indication of any failed calculations with error messages
- **File Generation**: Standard StoBe output files (.out, .err) in respective directories

## 🔧 Troubleshooting

### Common Issues

1. **Permission Errors**: Make sure the script is executable and you have write permissions in the target directories

2. **Missing Run Files**: Use `explore` command to verify run files exist in the expected locations

3. **Memory Issues**: Reduce the number of workers if you encounter memory problems:

   ```bash
   uv run schedule.py --workers 2 gnd tcta_flat
   ```

4. **Stuck Calculations**: The script shows live progress - if a calculation appears stuck, you can cancel with Ctrl+C

### Performance Tuning

- **Workers**: Start with 4-8 workers and adjust based on your system resources
- **Memory**: Each worker runs a separate StoBe calculation, so monitor memory usage
- **I/O**: For systems with slow I/O, fewer workers may be more efficient

## 🔗 Integration

This CLI can be easily integrated into batch processing scripts, job schedulers, or other workflow automation tools. The exit codes and structured output make it suitable for automated environments.

---

# 🌟 Complete Workflow Example

Here's how to use both tools together for a complete StoBe DFT workflow:

## 1. Generate Input Files

```bash
# Generate input files for your molecule
uv run generate.py tcta_flat

# Verify the structure
ls tcta_flat/
# Output: molConfig.py tcta_flat.xyz C1/ C2/ C3/ ... C18/
```

## 2. Run Calculations

```bash
# Start with ground state calculations
uv run schedule.py gnd tcta_flat

# Run excited state calculations
uv run schedule.py exc tcta_flat

# Or run everything in sequence
uv run schedule.py all tcta_flat
```

## 3. Monitor Progress

Watch the beautiful real-time progress bars as your calculations run in parallel!

## 🔄 Complete Migration Example

### From Original Scripts

```bash
# Old way (manual script editing and execution)
# 1. Edit StoBe_Input_Generator_v4.py with hard-coded paths
# 2. python StoBe_Input_Generator_v4.py
# 3. Manually run each calculation
# 4. Deal with file organization bugs

# New way (modern CLI workflow)
uv run generate.py my_molecule     # Generate files
uv run schedule.py all my_molecule # Run calculations
```

### Updated molConfig.py

```python
# Remove these old variables:
# fname = "hardcoded/path"  ❌
# mname = "molecule_name"   ❌

# Keep these essential variables:
aname = "C"                # ✅
nFiles = 18               # ✅
# ... rest of configuration
```

## 📝 License

This tool suite is part of the StoBe DFT calculation ecosystem. Please refer to your StoBe license for usage terms.

## 🤝 Contributing

Improvements and bug reports are welcome for both tools! Key areas for enhancement:

### Input Generator

- Additional element support for XAS calculations
- Extended basis set libraries
- Enhanced error diagnostics
- Additional output formats

### Scheduler

- Advanced job queuing systems
- Integration with cluster schedulers (SLURM, PBS)
- Results analysis and visualization
- Checkpoint and resume functionality

## 📚 Related Tools & Ecosystem

- **StoBe Analysis Tools**: For processing calculation results
- **Molecular Visualization**: For geometry verification and results display
- **Data Processing**: Scripts for parsing and analyzing output files
- **Cluster Integration**: HPC and cloud computing adapters

---

**🎯 Ready to run StoBe DFT calculations with confidence!**

_The modern CLI suite provides robust file organization, parallel execution, flexible input handling, and beautiful user interfaces for your DFT calculations._
