"""
Test script for StoBe output processing
"""

import sys
from pathlib import Path

# Add the project directory to the Python path
project_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_dir))

from src.dftlearn.io.stobe import (
    classify_stobe_file,
    extract_ionization_potential,
    extract_lumo_energy,
    extract_nexafs_params,
    extract_oscillator_strengths,
    extract_total_energy,
    process_excitation_center,
    read_output,
)


def test_ground_state():
    """Test reading a ground state file."""
    file_path = "output/C1gnd.out"
    print(f"Testing ground state file processing: {file_path}")

    # Test file classification
    file_type = classify_stobe_file(file_path)
    print(f"File type: {file_type}")

    # Test energy extraction
    energy = extract_total_energy(file_path)
    print(f"Total energy: {energy} Hartree")

    # Test full file reading
    data = read_output(file_path)
    print(f"Keys in output: {list(data.keys())}")
    print(f"Geometry entries: {len(data['geometry'])}")
    if data["geometry"]:
        print(f"First atom: {data['geometry'][0]}")


def test_transition_potential():
    """Test reading a transition potential file."""
    file_path = "output/C1tp.out"
    print(f"\nTesting transition potential file processing: {file_path}")

    # Test energy extraction
    energy = extract_total_energy(file_path)
    print(f"Total energy: {energy} Hartree")

    # Test LUMO energy extraction
    lumo = extract_lumo_energy(file_path)
    print(f"LUMO energy: {lumo} eV")

    # Test IP extraction
    ip = extract_ionization_potential(file_path)
    print(f"Ionization potential: {ip} eV")

    # Test oscillator strengths extraction
    os_data = extract_oscillator_strengths(file_path)
    print(f"Found {len(os_data['energies'])} transitions")
    if os_data["energies"]:
        print(
            f"First transition: {os_data['energies'][0]} eV, strength: {os_data['strengths'][0]}"
        )


def test_excitation_center():
    """Test processing a complete excitation center."""
    directory = "output"
    prefix = "C1"
    print(f"\nTesting complete excitation center processing for {prefix}")

    # Process the excitation center
    params = extract_nexafs_params(directory, prefix)
    print(f"Ground state energy: {params['ground_state_energy']} Hartree")
    print(f"Excited state energy: {params['excited_state_energy']} Hartree")
    print(f"TP energy: {params['tp_energy']} Hartree")
    print(f"Energy correction: {params['energy_correction']} eV")

    # Process the full center
    data = process_excitation_center(directory, prefix)
    if "interpolated_spectrum" in data:
        print(
            f"Generated interpolated spectrum with {len(data['interpolated_spectrum']['energy'])} points"
        )
    else:
        print("No interpolated spectrum generated")


def test_create_xarray_tree():
    """Test creating an xarray DataTree from StoBe output."""
    try:
        import xarray as xr
    except ImportError:
        print("\nxarray not installed, skipping DataTree test")
        return

    directory = "output"
    prefix = "C1"
    print(f"\nTesting creation of xarray DataTree for {prefix}")

    # Process the excitation center
    data = process_excitation_center(directory, prefix)

    # Create dataset for ground state
    coords = {}
    data_vars = {}

    # Add ground state energy
    if data["ground_state_energy"] is not None:
        data_vars["energy"] = data["ground_state_energy"]

    # Add geometry if available
    if "geometry" in data:
        # Create atom coordinates
        atoms = []
        x = []
        y = []
        z = []
        element = []

        for atom in data["geometry"]:
            atoms.append(atom["atom"])
            x.append(atom["x"])
            y.append(atom["y"])
            z.append(atom["z"])
            element.append(
                atom["atom"][0]
            )  # First character is usually the element symbol

        coords["atom"] = atoms
        data_vars["x"] = ("atom", x)
        data_vars["y"] = ("atom", y)
        data_vars["z"] = ("atom", z)
        data_vars["element"] = ("atom", element)

    # Create the dataset
    ds_ground = xr.Dataset(
        data_vars=data_vars, coords=coords, attrs={"type": "ground_state"}
    )
    print(
        f"Created ground state dataset with {len(ds_ground['atom']) if 'atom' in ds_ground.coords else 0} atoms"
    )

    # Create dataset for spectrum
    if "interpolated_spectrum" in data:
        spectrum_data = {
            "energy": data["interpolated_spectrum"]["energy"],
            "intensity": ("energy", data["interpolated_spectrum"]["intensity"]),
        }
        ds_spectrum = xr.Dataset(
            data_vars={"intensity": spectrum_data["intensity"]},
            coords={"energy": spectrum_data["energy"]},
            attrs={"type": "spectrum"},
        )
        print(
            f"Created spectrum dataset with {len(ds_spectrum['energy'])} points"
        )

        # Create a DataTree
        tree = xr.DataTree()
        tree["ground_state"] = ds_ground
        tree["spectrum"] = ds_spectrum

        print("Successfully created xarray DataTree")
        print(f"Tree structure: {list(tree.groups.keys())}")
    else:
        print("No spectrum data available for DataTree")


if __name__ == "__main__":
    test_ground_state()
    test_transition_potential()
    test_excitation_center()
    test_create_xarray_tree()
