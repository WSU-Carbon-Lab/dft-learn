"""
Module for creating xarray DataTree objects from StoBe calculations
"""

import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from .stobe import (
    process_excitation_center,
    process_multiple_centers,
    combine_spectra,
    HARTREE_TO_EV
)

def create_datatree(directory: str, prefix: str) -> Any:
    """
    Create an xarray DataTree from StoBe calculations for a single excitation center.

    Args:
        directory: Directory containing StoBe output files
        prefix: Prefix for the excitation center (e.g., "C1")

    Returns:
        xarray.DataTree: DataTree containing all data for the excitation center
    """
    try:
        import xarray as xr
        import numpy as np
    except ImportError:
        raise ImportError("xarray and numpy must be installed to create DataTree objects")

    # Process the excitation center
    data = process_excitation_center(directory, prefix)

    # Create a DataTree
    tree = xr.DataTree()

    # Create ground state dataset
    if data['ground_state_energy'] is not None:
        # Create dataset with coordinates
        coords = {}
        data_vars = {
            'energy': data['ground_state_energy'],
            'energy_ev': data['ground_state_energy'] * HARTREE_TO_EV
        }

        # Add geometry if available
        if data['geometry']:
            # Create atom coordinates
            atom_ids = []
            x, y, z = [], [], []
            elements = []

            for atom in data['geometry']:
                atom_ids.append(atom['atom'])
                x.append(atom['x'])
                y.append(atom['y'])
                z.append(atom['z'])
                elements.append(atom['atom'][0])  # First character is usually the element

            coords['atom'] = atom_ids
            data_vars['x'] = ('atom', x)
            data_vars['y'] = ('atom', y)
            data_vars['z'] = ('atom', z)
            data_vars['element'] = ('atom', elements)

        ds_ground = xr.Dataset(
            data_vars=data_vars,
            coords=coords,
            attrs={'type': 'ground_state', 'prefix': prefix}
        )
        tree['ground_state'] = ds_ground

    # Create excited state dataset
    if data['excited_state_energy'] is not None:
        ds_excited = xr.Dataset(
            data_vars={
                'energy': data['excited_state_energy'],
                'energy_ev': data['excited_state_energy'] * HARTREE_TO_EV
            },
            attrs={'type': 'excited_state', 'prefix': prefix}
        )
        tree['excited_state'] = ds_excited

    # Create transition potential dataset
    if data['tp_energy'] is not None:
        tp_vars = {
            'energy': data['tp_energy'],
            'energy_ev': data['tp_energy'] * HARTREE_TO_EV
        }

        # Add LUMO energy if available
        if data['lumo_energy'] is not None:
            tp_vars['lumo_energy'] = data['lumo_energy']

        # Add IP if available
        if data['ip'] is not None:
            tp_vars['ionization_potential'] = data['ip']

        ds_tp = xr.Dataset(
            data_vars=tp_vars,
            attrs={'type': 'transition_potential', 'prefix': prefix}
        )
        tree['transition_potential'] = ds_tp

    # Add oscillator strengths if available
    if data['oscillator_strengths'] and data['oscillator_strengths']['energies']:
        os_data = data['oscillator_strengths']

        # Create coordinates for transitions
        n_transitions = len(os_data['energies'])

        # Create dataset for oscillator strengths
        ds_os = xr.Dataset(
            data_vars={
                'energy': ('transition', os_data['energies']),
                'strength': ('transition', os_data['strengths']),
                'tdm_x': ('transition', os_data['tdm_x']),
                'tdm_y': ('transition', os_data['tdm_y']),
                'tdm_z': ('transition', os_data['tdm_z'])
            },
            coords={
                'transition': np.arange(n_transitions)
            },
            attrs={'type': 'oscillator_strengths', 'prefix': prefix}
        )

        # Add corrected energies if available
        if 'corrected_energies' in os_data:
            ds_os['corrected_energy'] = ('transition', os_data['corrected_energies'])

        tree['oscillator_strengths'] = ds_os

    # Add spectrum if available
    if 'interpolated_spectrum' in data and data['interpolated_spectrum']:
        spectrum = data['interpolated_spectrum']

        ds_spectrum = xr.Dataset(
            data_vars={
                'intensity': ('energy', spectrum['intensity'])
            },
            coords={
                'energy': spectrum['energy']
            },
            attrs={
                'type': 'spectrum',
                'prefix': prefix,
                'energy_correction': data['energy_correction'] if data['energy_correction'] is not None else 0.0
            }
        )
        tree['spectrum'] = ds_spectrum

    # Add XAS data if available
    if data['xas_spectrum'] and data['xas_spectrum']['energies']:
        xas = data['xas_spectrum']

        ds_xas = xr.Dataset(
            data_vars={
                'intensity': ('energy', xas['intensities'])
            },
            coords={
                'energy': xas['energies']
            },
            attrs={'type': 'xas', 'prefix': prefix}
        )
        tree['xas'] = ds_xas

    return tree

def create_multi_center_datatree(directory: str, prefixes: List[str]) -> Any:
    """
    Create an xarray DataTree from StoBe calculations for multiple excitation centers.

    Args:
        directory: Directory containing StoBe output files
        prefixes: List of excitation center prefixes (e.g., ["C1", "C2", "N1"])

    Returns:
        xarray.DataTree: DataTree containing data for all excitation centers
    """
    try:
        import xarray as xr
    except ImportError:
        raise ImportError("xarray must be installed to create DataTree objects")

    # Create the main DataTree
    tree = xr.DataTree()

    # Process each center and add to the tree
    for prefix in prefixes:
        center_tree = create_datatree(directory, prefix)
        if center_tree is not None:
            tree[prefix] = center_tree

    # Add combined spectrum if possible
    try:
        # Process all centers to get their data
        centers_data = process_multiple_centers(directory, prefixes)

        # Combine the spectra
        combined = combine_spectra(centers_data)

        if combined['energy'].size > 0:
            # Create a dataset for the combined spectrum
            ds_combined = xr.Dataset(
                data_vars={
                    'intensity': ('energy', combined['intensity'])
                },
                coords={
                    'energy': combined['energy']
                },
                attrs={'type': 'combined_spectrum'}
            )
            tree['combined_spectrum'] = ds_combined
    except Exception as e:
        print(f"Could not create combined spectrum: {str(e)}")

    return tree

def save_datatree(tree: Any, output_file: str) -> None:
    """
    Save an xarray DataTree to a NetCDF file.

    Args:
        tree: xarray DataTree to save
        output_file: Path to the output NetCDF file
    """
    try:
        tree.to_netcdf(output_file)
        print(f"DataTree saved to {output_file}")
    except Exception as e:
        print(f"Error saving DataTree: {str(e)}")

def load_datatree(file_path: str) -> Any:
    """
    Load an xarray DataTree from a NetCDF file.

    Args:
        file_path: Path to the NetCDF file

    Returns:
        xarray.DataTree: Loaded DataTree
    """
    try:
        import xarray as xr
        return xr.open_datatree(file_path)
    except Exception as e:
        print(f"Error loading DataTree: {str(e)}")
        return None
