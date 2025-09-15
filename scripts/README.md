# StoBe DFT Calculation Scheduler CLI

A command-line interface for scheduling and running StoBe DFT calculations with parallel execution and **real-time Typer progress bars**.

## Features

- **Parallel Execution**: Run multiple calculations simultaneously with configurable worker limits
- **Rich Progress Bars**: Live updating progress bars with animated spinners using Typer/Rich
- **Flexible Targeting**: Run calculations on entire directories or specific atoms
- **Multiple Calculation Types**: Support for ground state (gnd), excited state (exc), transition potential (tp), XAS, and sequential calculations
- **Beautiful Terminal Output**: Professional-grade progress display with colors and animations
- **UV Script**: Self-contained script with automatic dependency management

## Installation

This script uses [UV](https://docs.astral.sh/uv/) for dependency management. Make sure you have UV installed:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

No additional installation is needed - UV will automatically handle the dependencies (`typer` and `rich`) when you run the script.

## Usage

### Basic Commands

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

### Advanced Options

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

During execution, you'll see beautiful Typer progress bars like this:

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
- **Status Updates**: Clear indication of completed (✅), running, pending, and failed calculations## Examples

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

## Output and Results

- **Live Progress**: Real-time updates every 2 seconds during execution
- **Completion Summary**: Final statistics showing completed, failed, and timing information
- **Error Reporting**: Clear indication of any failed calculations with error messages
- **File Generation**: Standard StoBe output files (.out, .err) in respective directories

## Troubleshooting

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

## Technical Details

- **Dependencies**: Automatically managed by UV (`typer`, `rich`)
- **Python Version**: Requires Python >=3.8
- **Execution Model**: Thread-per-directory parallel execution
- **Progress System**: Typer/Rich progress bars with real-time updates and animations
- **Output**: Beautiful terminal display with colored status indicators and spinners

## Integration

This CLI can be easily integrated into batch processing scripts, job schedulers, or other workflow automation tools. The exit codes and structured output make it suitable for automated environments.

**🚀 The Typer-based system provides beautiful, efficient parallel StoBe calculations with stunning real-time progress bars and professional terminal formatting!**
