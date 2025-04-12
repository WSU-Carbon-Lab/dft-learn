"""
Demonstration script for StoBe data processing with xarray DataTree.
"""

import sys
from pathlib import Path

# Add the project directory to the Python path
project_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_dir))

# Import the modules
from src.dftlearn.io.datatree import create_datatree, save_datatree
from src.dftlearn.io.stobe import (
    extract_lumo_energy,
    extract_total_energy,
    process_excitation_center,
    save_to_csv,
)


def main():
    """Process StoBe data and create xarray DataTree."""
    print("StoBe Data Processing Demonstration")
    print("===================================")

    # Set directory and prefix
    # Set directory and prefix
    directory = "output"
    prefix = "C1"
    print(f"Processing files with prefix '{prefix}' in directory '{directory}'")

    # 1. Extract basic data
    print("\n1. Basic Data Extraction")
    print("----------------------")

    gnd_file = Path(directory) / f"{prefix}gnd.out"
    exc_file = Path(directory) / f"{prefix}exc.out"
    tp_file = Path(directory) / f"{prefix}tp.out"
    gnd_energy = extract_total_energy(str(gnd_file))
    exc_energy = extract_total_energy(str(exc_file))
    tp_energy = extract_total_energy(str(tp_file))
    lumo_energy = extract_lumo_energy(str(tp_file))

    print(f"Ground state energy: {gnd_energy} Hartree")
    print(f"Excited state energy: {exc_energy} Hartree")
    print(f"TP energy: {tp_energy} Hartree")
    print(f"LUMO energy: {lumo_energy} eV")

    # 2. Process excitation center
    print("\n2. Process Excitation Center")
    print("--------------------------")

    data = process_excitation_center(directory, prefix)

    print(f"Energy correction: {data['energy_correction']} eV")

    if data.get("oscillator_strengths"):
        print(
            f"Found {len(data['oscillator_strengths']['energies'])} transitions"
        )

    if data.get("geometry"):
        print(f"Found {len(data['geometry'])} atoms in geometry")

    # 3. Save spectrum to CSV
    # 3. Save spectrum to CSV
    print("\n3. Save Spectrum to CSV")
    print("--------------------")

    csv_file = Path(directory) / f"{prefix}_spectrum.csv"
    save_to_csv(data, csv_file)
    print(f"Saved spectrum to {csv_file}")
    # 4. Create xarray DataTree
    print("\n4. Create xarray DataTree")
    print("----------------------")

    try:
        import numpy as np
        import xarray as xr

        # Create DataTree
        tree = create_datatree(directory, prefix)
        print(f"Created DataTree with {len(tree.groups)} groups:")
        for group in tree.groups:
            print(f"  - {group}")

        # Save DataTree to NetCDF
        # Save DataTree to NetCDF
        netcdf_file = Path(directory) / f"{prefix}_data.nc"
        save_datatree(tree, netcdf_file)
    except ImportError as e:
        print(f"Could not create DataTree: {e}")
        print("Install xarray and numpy to use DataTree functionality.")


if __name__ == "__main__":
    main()
