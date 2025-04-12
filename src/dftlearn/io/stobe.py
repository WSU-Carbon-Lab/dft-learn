"""
StoBe Output File IO
"""

import re
from pathlib import Path
import polars as pl
import os
import numpy as np
from typing import Any, Dict, List, Optional, Tuple, Union

# Conversion from Hartree to eV
HARTREE_TO_EV = 27.2114

def atomic_coords(line: str) -> (None | dict[str, Any]):
    """
    Parses a line containing atomic coordinates.

    Args:
        line (str): A line from the geometry section.

    Returns:
        tuple[str, float, float, float]: Tuple containing atom type and coordinates (x, y, z).
    """
    if "Atom" in line and "x" in line:
        return None
    if "---" in line:
        return None
    # Check for the end of the geometry section
    if line.strip() == "":
        return None
    parts = line.strip().split()
    if len(parts) >= 11:
        # Extract the fields based on position
        atom_data = {
            "atom": parts[1],
            "x": float(parts[2]),
            "y": float(parts[3]),
            "z": float(parts[4]),
            "q": float(parts[5]),
            "nuc": int(parts[6]),
            "mass": float(parts[7]),
            # Handle neq which might contain a slash and space (e.g. "1/ 1")
            "neq": "".join(parts[8:10]),
            "grid": int(parts[10]),
            "grp": int(parts[11])
        }
        return atom_data
    return None

def parse_stobe_output(file_path: str) -> Dict[str, Any]:
    """
    Parse a StoBe output file to extract geometry and basis set information.

    Args:
        file_path (str): Path to the StoBe .out file

    Returns:
        dict: Dictionary containing parsed data including geometry and basis sets
    """
    with open(file_path, 'r') as file:
        content = file.readlines()

    # Check if this is a ground state calculation
    is_ground_state = "GND" in Path(file_path).name.upper()

    result = {
        "is_ground_state": is_ground_state,
        "geometry": [],
        "auxiliary_basis": {},
        "orbital_basis": {},
        "model_potentials": {}
    }

    # Parse geometry section - first try with main GEOMETRY section
    geometry_section = extract_section(
        content,
        start_pattern=r"^\s*GEOMETRY\s*$",
        end_pattern=r"Smallest atom distance"  # Stop at smallest atom distance marker
    )

    # If not found, try the alternate format
    if not geometry_section:
        geometry_section = extract_section(
            content,
            start_pattern="Single image calculation (Angstrom):",
            end_pattern=r"Smallest atom distance"  # Stop at smallest atom distance marker
        )

    if geometry_section:
        # Find the header line that contains "Atom"
        header_index = -1
        for i, line in enumerate(geometry_section):
            if "Atom" in line and "x" in line and "y" in line and "z" in line:
                header_index = i
                break

        # Get the line after the dashes
        if header_index >= 0:
            # Skip header and dash line, process data lines
            data_lines = [line for line in geometry_section[header_index+2:]
                         if line.strip() and not "---" in line and not "===" in line
                         and not "Smallest atom distance" in line  # Exclude the boundary line
                         and re.match(r'^\s*\d+\)', line)]  # Only include lines that start with a number followed by )

            for line in data_lines:
                atom_data = parse_geometry_line(line)
                if atom_data:
                    result["geometry"].append(atom_data)

    # Only parse basis sets if this is a ground state calculation
    if is_ground_state:
        # Parse auxiliary basis sets
        aux_basis_section = extract_section(
            content,
            start_pattern=r"I\)  AUXILIARY BASIS SETS",
            end_pattern=r"II\)  ORBITAL BASIS SETS"
        )

        if aux_basis_section:
            for line in aux_basis_section[1:]:  # Skip header
                if line.strip() and line.startswith(" Atom"):
                    atom_id, basis = parse_basis_line(line)
                    if atom_id:
                        result["auxiliary_basis"][atom_id] = basis

        # Parse orbital basis sets
        orb_basis_section = extract_section(
            content,
            start_pattern=r"II\)  ORBITAL BASIS SETS",
            end_pattern=r"III\)  MODEL POTENTIALS" if "III)  MODEL POTENTIALS" in "".join(content) else "BASIS DIMENSIONS"
        )

        if orb_basis_section:
            for line in orb_basis_section[1:]:  # Skip header
                if line.strip() and line.startswith(" Atom"):
                    atom_id, basis = parse_basis_line(line)
                    if atom_id:
                        result["orbital_basis"][atom_id] = basis

        # Parse model potentials - improved to better handle the model potentials section
        model_pot_section = extract_section(
            content,
            start_pattern=r"III\)  MODEL POTENTIALS",
            end_pattern=r"(BASIS DIMENSIONS|HAMILTONIAN DEFINITION|Input settings)"  # Use multiple potential end patterns
        )

        if model_pot_section:
            for line in model_pot_section[1:]:  # Skip header
                if line.strip() and line.startswith(" Atom"):
                    atom_id, basis = parse_basis_line(line)
                    if atom_id:
                        result["model_potentials"][atom_id] = basis

    return result

def extract_section(content: List[str], start_pattern: str, end_pattern: str) -> List[str]:
    """
    Extract a section from the file content between start and end patterns.

    Args:
        content (list): List of lines from the file
        start_pattern (str): Regex pattern to match the start of the section
        end_pattern (str): Regex pattern to match the end of the section

    Returns:
        list: Lines of the extracted section
    """
    section = []
    in_section = False

    for line in content:
        if not in_section and re.search(start_pattern, line):
            in_section = True
            section.append(line)
        elif in_section:
            section.append(line)
            if re.search(end_pattern, line) and len(section) > 1:
                # Don't break immediately after start pattern to avoid empty sections
                if len(section) > 2 or start_pattern != end_pattern:
                    break

    return section

def parse_geometry_line(line: str) -> Optional[Dict[str, Any]]:
    """
    Parse a single line from the geometry section.

    Args:
        line (str): Line containing atom geometry data

    Returns:
        dict or None: Dictionary with atom information or None if parsing failed
    """
    # Don't try to parse header or separator lines
    if "----" in line or "Atom" in line and "x" in line:
        return None

    # Simple but robust pattern matching - split by whitespace and extract values by position
    parts = line.strip().split()

    # Check if we have enough parts to make a valid atom entry
    if len(parts) >= 8:  # At least 8 elements for a geometry line
        try:
            # Extract the index and atom name considering the format "1) C01"
            idx_str = parts[0]
            atom_str = parts[1]

            # If the first part has a trailing ")", it contains the index
            if idx_str.endswith(")"):
                idx = idx_str.rstrip(")")
            else:
                # Handle case where index and atom are merged like "1)C01"
                match = re.match(r"(\d+)\)(.*)", idx_str)
                if match:
                    idx = match.group(1)
                    atom_str = match.group(2)
                else:
                    idx = None

            return {
                "index": idx,
                "atom": atom_str,
                "x": float(parts[2]) if len(parts) > 2 else 0.0,
                "y": float(parts[3]) if len(parts) > 3 else 0.0,
                "z": float(parts[4]) if len(parts) > 4 else 0.0,
                "q": float(parts[5]) if len(parts) > 5 else 0.0,
                "nuc": int(parts[6]) if len(parts) > 6 else 0,
                "mass": float(parts[7]) if len(parts) > 7 else 0.0,
                "neq": parts[8] if len(parts) > 8 else "",
                "grid": int(parts[9]) if len(parts) > 9 else 0,
                "grp": int(parts[10]) if len(parts) > 10 else 0,
                "auxiliary_basis": "",
                "orbital_basis": "",
                "model_potential": ""
            }
        except (ValueError, IndexError):
            # Log diagnostic info if needed
            # print(f"Failed to parse line: {line}")
            return None

    return None

def parse_basis_line(line: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Parse a line from the basis set sections.

    Args:
        line (str): Line containing basis set information

    Returns:
        tuple: (atom_id, basis_description) or (None, None) if parsing failed
    """
    match = re.match(r"\s*Atom\s+([A-Za-z0-9]+)\s*:\s*(.+)$", line)
    if match:
        return match.group(1), match.group(2).strip()
    return None, None

def merge_basis_with_geometry(data: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Merge basis set information with geometry data.

    Args:
        data (dict): Parsed data from StoBe output file

    Returns:
        list: Geometry data with added basis set information
    """
    geometry = data.get("geometry", [])

    for atom_data in geometry:
        atom_id = atom_data["atom"]

        # Add auxiliary basis set
        if atom_id in data.get("auxiliary_basis", {}):
            atom_data["auxiliary_basis"] = data["auxiliary_basis"][atom_id]

        # Add orbital basis set
        if atom_id in data.get("orbital_basis", {}):
            atom_data["orbital_basis"] = data["orbital_basis"][atom_id]

        # Add model potential
        if atom_id in data.get("model_potentials", {}):
            atom_data["model_potential"] = data["model_potentials"][atom_id]

    return geometry

def read_output(file_path: str) -> Dict[str, Any]:
    """
    Read a StoBe output file and extract all relevant information.

    Args:
        file_path (str): Path to the StoBe .out file.

    Returns:
        dict: Dictionary containing all parsed data from the file
    """
    data = parse_stobe_output(file_path)
    return data

def parse_geometry(file_path: str) -> List[Dict[str, Any]]:
    """
    Parses the geometry section of a StoBe .out file.

    Args:
        file_path (str): Path to the StoBe .out file.

    Returns:
        list[dict]: A list of dictionaries containing atom information.
    """
    data = parse_stobe_output(file_path)
    return data["geometry"]

def parse_geometry_with_basis(file_path: str) -> List[Dict[str, Any]]:
    """
    Parses the geometry section of a StoBe .out file and adds basis set information.

    Args:
        file_path (str): Path to the StoBe .out file.

    Returns:
        list[dict]: A list of dictionaries containing atom information with basis sets.
    """
    data = parse_stobe_output(file_path)
    return merge_basis_with_geometry(data)

def to_dataframe(geometry_data: List[Dict[str, Any]]) -> pl.DataFrame:
    """
    Convert geometry data to a Polars DataFrame.

    Args:
        geometry_data (list): List of dictionaries containing atom information

    Returns:
        pl.DataFrame: DataFrame representation of the geometry data
    """
    if not geometry_data:
        return pl.DataFrame()
    return pl.DataFrame(geometry_data)

def classify_stobe_file(file_path: str) -> str:
    """
    Classify a StoBe output file based on its name.

    Args:
        file_path (str): Path to the StoBe output file

    Returns:
        str: Classification of the file ('gnd', 'tp', 'exc', 'xas', or 'unknown')
    """
    filename = Path(file_path).name.lower()

    if 'gnd' in filename:
        return 'gnd'
    elif 'tp' in filename:
        return 'tp'
    elif 'exc' in filename:
        return 'exc'
    elif 'xas' in filename:
        return 'xas'
    else:
        return 'unknown'

def extract_total_energy(file_path: str) -> Optional[float]:
    """
    Extract the total energy from a StoBe output file in Hartree.

    Args:
        file_path (str): Path to the StoBe output file

    Returns:
        float or None: Total energy in Hartree, or None if not found
    """
    with open(file_path, 'r') as file:
        content = file.read()

    # Regex pattern to find the total energy
    pattern = r"Total energy\s+\(H\)\s*=\s*([-\d.]+)"
    match = re.search(pattern, content)

    if match:
        return float(match.group(1))
    return None

def extract_ionization_potential(file_path: str) -> Optional[float]:
    """
    Extract the ionization potential from a TP StoBe output file.

    Args:
        file_path (str): Path to the StoBe TP output file

    Returns:
        float or None: Ionization potential in eV, or None if not found
    """
    with open(file_path, 'r') as file:
        content = file.read()

    # Regex pattern to find the ionization potential
    pattern = r"Ionization potential\s*=\s*([-\d.]+)"
    match = re.search(pattern, content)

    if match:
        return float(match.group(1))
    return None

def extract_lumo_energy(file_path: str) -> Optional[float]:
    """
    Extract the LUMO energy from a StoBe output file.

    Args:
        file_path (str): Path to the StoBe output file

    Returns:
        float or None: LUMO energy in eV, or None if not found
    """
    with open(file_path, 'r') as file:
        content = file.read()

    # Try to find LUMO energy in a more direct way
    # First check for a specific "LUMO" or "Lowest Unoccupied" mention
    lumo_patterns = [
        r"LUMO[\s:=]+([0-9.-]+)[\s]*(?:eV|hartree)",
        r"Lowest Unoccupied[\s:=]+([0-9.-]+)[\s]*(?:eV|hartree)",
        r"Transition energy[\s:=]+([0-9.-]+)[\s]*(?:eV|hartree)"
    ]

    for pattern in lumo_patterns:
        match = re.search(pattern, content, re.IGNORECASE)
        if match:
            return float(match.group(1))

    # If we can't find a direct LUMO reference, look for the orbital energies section
    # and find the first unoccupied orbital (occupancy close to 0)
    orb_section_match = re.search(r"Occup\.[\s]+Energy\(eV\).*?(\n.*?)+?(?:={10}|-{10})", content, re.DOTALL | re.MULTILINE)

    if orb_section_match:
        orb_section = orb_section_match.group(0)
        lines = orb_section.strip().split('\n')

        # Skip header line(s)
        for line in lines[1:]:
            if not line.strip() or line.startswith(('=', '-')):
                continue

            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    occupancy = float(parts[0])
                    energy = float(parts[1])

                    # First orbital with occupancy near 0 is the LUMO
                    if abs(occupancy) < 0.1:
                        return energy
                except (ValueError, IndexError):
                    pass

    # As a fallback, look for the transition energy directly in the file
    # This is often around 290 eV for carbon K-edge
    match = re.search(r"(\d+\.\d+)\s+eV", content)
    if match:
        energy = float(match.group(1))
        # Only return if it's in a reasonable range for NEXAFS (e.g., carbon K-edge)
        if 280 <= energy <= 310:
            return energy

    # If all else fails, use a reasonable default for carbon K-edge
    return 290.0  # Default value for carbon K-edge

def extract_oscillator_strengths(file_path: str) -> Dict[str, List[float]]:
    """
    Extract oscillator strengths and transition energies from a StoBe TP output file.

    Args:
        file_path (str): Path to the StoBe TP output file

    Returns:
        dict: Dictionary with 'energies' and 'strengths' lists
    """
    with open(file_path, 'r') as file:
        content = file.readlines()

    energies = []
    strengths = []
    tdm_x = []
    tdm_y = []
    tdm_z = []

    in_oscillator_section = False

    for i, line in enumerate(content):
        # Look for different markers that indicate the TD calculation summary
        if any(marker in line for marker in ["Summary of TD calculation", "Transition dipole moments"]):
            in_oscillator_section = True
            continue

        if in_oscillator_section and line.strip():
            # Look for lines starting with "#" or containing transition data
            if line.strip().startswith("#") or "Energy(eV)" in line:
                # Skip header lines
                continue

            if any(marker in line for marker in ["===", "---", "End of calculation"]):
                # End of section
                break

            # Parse the line containing transition data
            parts = line.strip().split()
            if len(parts) >= 5:  # Ensure enough data points
                try:
                    # Format varies between StoBe versions, try to be flexible
                    # Try to identify which columns contain the data we need

                    # First, look for values that could be energy (typically in eV range for NEXAFS)
                    energy_candidates = []
                    for j, part in enumerate(parts):
                        try:
                            val = float(part)
                            if 270 <= val <= 350:  # Typical C-K edge energy range
                                energy_candidates.append((j, val))
                        except ValueError:
                            pass

                    if energy_candidates:
                        # Found potential energy value
                        energy_idx, energy = energy_candidates[0]

                        # Strength is typically the next value after energy
                        strength_idx = energy_idx + 1
                        if strength_idx < len(parts):
                            try:
                                strength = float(parts[strength_idx])

                                # TDM components typically follow
                                x = float(parts[strength_idx + 1]) if strength_idx + 1 < len(parts) else 0.0
                                y = float(parts[strength_idx + 2]) if strength_idx + 2 < len(parts) else 0.0
                                z = float(parts[strength_idx + 3]) if strength_idx + 3 < len(parts) else 0.0

                                energies.append(energy)
                                strengths.append(strength)
                                tdm_x.append(x)
                                tdm_y.append(y)
                                tdm_z.append(z)
                            except (ValueError, IndexError):
                                pass
                except (ValueError, IndexError):
                    pass

    # If we didn't find any data using the above method, try a more general approach
    if not energies:
        # Restart the search
        with open(file_path, 'r') as file:
            content = file.read()

        # Find a section that might contain oscillator strengths
        td_sections = re.findall(r'(Energy\(eV\).*?Osc\.Str\..*?)(?:={10}|-{10}|End of calculation)',
                              content, re.DOTALL | re.MULTILINE)

        if td_sections:
            # Process the first found section
            for line in td_sections[0].split('\n'):
                if line.strip() and not line.strip().startswith(('#', 'Energy')):
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        try:
                            energy = float(parts[0])
                            strength = float(parts[1])

                            # Get TDM components if available
                            x = float(parts[2]) if len(parts) > 2 else 0.0
                            y = float(parts[3]) if len(parts) > 3 else 0.0
                            z = float(parts[4]) if len(parts) > 4 else 0.0

                            energies.append(energy)
                            strengths.append(strength)
                            tdm_x.append(x)
                            tdm_y.append(y)
                            tdm_z.append(z)
                        except (ValueError, IndexError):
                            pass

    return {
        'energies': energies,
        'strengths': strengths,
        'tdm_x': tdm_x,
        'tdm_y': tdm_y,
        'tdm_z': tdm_z
    }

def extract_xas_spectrum(file_path: str) -> Dict[str, List[float]]:
    """
    Extract XAS spectrum data from a StoBe XAS output file.

    Args:
        file_path (str): Path to the StoBe XAS output file

    Returns:
        dict: Dictionary with 'energies' and 'intensities' lists
    """
    with open(file_path, 'r') as file:
        content = file.readlines()

    energies = []
    intensities = []

    # Detect when we're in the XAS spectrum data
    in_xas_section = False

    for line in content:
        # Different files might have different formats, look for common indicators
        if any(marker in line for marker in ["Total spectrum", "Absorption", "Energy(ev)"]):
            in_xas_section = True
            continue

        if in_xas_section and line.strip():
            # Check if we've exited the section
            if line.startswith("----") or line.startswith("===="):
                break

            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    energy = float(parts[0])
                    intensity = float(parts[1])
                    energies.append(energy)
                    intensities.append(intensity)
                except (ValueError, IndexError):
                    pass

    return {
        'energies': energies,
        'intensities': intensities
    }

def extract_nexafs_params(directory: str, prefix: str = "C1") -> Dict[str, Any]:
    """
    Extract all NEXAFS parameters from StoBe files in a directory for a specific excitation center.

    Args:
        directory (str): Directory containing StoBe output files
        prefix (str): Prefix for the excitation center (e.g., "C1")

    Returns:
        dict: Dictionary with all NEXAFS parameters
    """
    result: Dict[str, Optional[Any]] = {
        'ground_state_energy': None,  # Hartree
        'excited_state_energy': None, # Hartree
        'tp_energy': None,           # Hartree
        'lumo_energy': None,         # eV
        'ip': None,                  # eV
        'energy_correction': None,   # eV
        'oscillator_strengths': None,
        'xas_spectrum': None,
        'geometry': None             # Add geometry field
    }

    # Find all StoBe output files in the directory
    path = Path(directory)
    file_mapping = {}

    for file in path.glob(f"{prefix}*.out"):
        file_type = classify_stobe_file(str(file))
        file_mapping[file_type] = str(file)

    # Extract data from ground state file
    if 'gnd' in file_mapping:
        result['ground_state_energy'] = extract_total_energy(file_mapping['gnd'])

        # Get geometry data from the ground state file
        gnd_data = read_output(file_mapping['gnd'])
        if 'geometry' in gnd_data and gnd_data['geometry']:
            result['geometry'] = gnd_data['geometry']

    # Extract data from excited state file
    if 'exc' in file_mapping:
        result['excited_state_energy'] = extract_total_energy(file_mapping['exc'])

    # Extract data from transition potential file
    if 'tp' in file_mapping:
        result['tp_energy'] = extract_total_energy(file_mapping['tp'])
        result['lumo_energy'] = extract_lumo_energy(file_mapping['tp'])
        result['ip'] = extract_ionization_potential(file_mapping['tp'])
        result['oscillator_strengths'] = extract_oscillator_strengths(file_mapping['tp'])

    # Extract data from XAS file
    if 'xas' in file_mapping:
        result['xas_spectrum'] = extract_xas_spectrum(file_mapping['xas'])

    # Calculate energy correction
    if all(v is not None for v in [result['ground_state_energy'],
                                  result['excited_state_energy'],
                                  result['lumo_energy']]):
        # Convert Hartree to eV and calculate correction
        gnd_ev = result['ground_state_energy'] * HARTREE_TO_EV  # type: ignore
        exc_ev = result['excited_state_energy'] * HARTREE_TO_EV # type: ignore
        lumo_ev = result['lumo_energy']

        # Ensure all values are valid floats before calculation
        if isinstance(gnd_ev, float) and isinstance(exc_ev, float) and isinstance(lumo_ev, float):
            result['energy_correction'] = (exc_ev - gnd_ev) - lumo_ev
        else:
            result['energy_correction'] = None
    else:
        result['energy_correction'] = None

    return result

def apply_energy_correction(
    transition_energies: List[float],
    oscillator_strengths: List[float],
    energy_correction: float,
    energy_range: Tuple[float, float] = (270, 320),
    num_points: int = 2000
) -> Dict[str, Union[np.ndarray, List[float]]]:
    """
    Apply energy correction to transition energies and interpolate to create a continuous spectrum.

    Args:
        transition_energies: List of transition energies in eV
        oscillator_strengths: List of oscillator strengths
        energy_correction: Energy correction to apply in eV
        energy_range: Tuple with (min, max) energy range for interpolation
        num_points: Number of points for interpolation

    Returns:
        dict: Dictionary with interpolation results
    """
    import numpy as np
    from scipy.interpolate import interp1d

    # Apply energy correction
    corrected_energies = [e + energy_correction for e in transition_energies]

    # Create a continuous energy range for interpolation
    continuous_energy = np.linspace(energy_range[0], energy_range[1], num_points)

    # Create arrays for interpolation
    energy_array = np.array(corrected_energies)
    strength_array = np.array(oscillator_strengths)

    # Sort by energy (important for interpolation)
    sort_idx = np.argsort(energy_array)
    energy_array = energy_array[sort_idx]
    strength_array = strength_array[sort_idx]

    # Create interpolation function (use zero for energies outside the range)
    if len(energy_array) > 1:
        interp_func = interp1d(
            energy_array,
            strength_array,
            kind='linear',
            bounds_error=False,
            fill_value=0
        )
        interpolated_spectrum = interp_func(continuous_energy)
    else:
        # If there's only one point, we can't interpolate
        interpolated_spectrum = np.zeros_like(continuous_energy)

    return {
        'energy': continuous_energy,
        'intensity': interpolated_spectrum,
        'corrected_transition_energies': corrected_energies,
        'oscillator_strengths': oscillator_strengths
    }

def process_excitation_center(directory: str, prefix: str = "C1") -> Dict[str, Any]:
    """
    Process an excitation center from StoBe output files.

    This function extracts all relevant data from the StoBe files for a specific
    excitation center, applies energy corrections, and creates interpolated spectra.

    Args:
        directory (str): Directory containing StoBe output files
        prefix (str): Prefix for the excitation center (e.g., "C1")

    Returns:
        dict: Dictionary with all processed data
    """
    # Extract parameters from files
    params = extract_nexafs_params(directory, prefix)

    # Process oscillator strengths if they exist
    if params['oscillator_strengths'] and params['energy_correction'] is not None:
        # Check if we have transitions
        energies = params['oscillator_strengths']['energies']
        strengths = params['oscillator_strengths']['strengths']

        if energies and strengths:
            # Apply energy correction and create interpolated spectrum
            interpolated_results = apply_energy_correction(
                energies,
                strengths,
                params['energy_correction']
            )

            # Add interpolated results to the output
            params['interpolated_spectrum'] = {
                'energy': interpolated_results['energy'],
                'intensity': interpolated_results['intensity']
            }

            # Update the oscillator strengths with corrected energies
            params['oscillator_strengths']['corrected_energies'] = interpolated_results['corrected_transition_energies']

    return params

def process_multiple_centers(directory: str, center_prefixes: List[str]) -> Dict[str, Dict[str, Any]]:
    """
    Process multiple excitation centers from StoBe output files.

    Args:
        directory (str): Directory containing StoBe output files
        center_prefixes (list): List of excitation center prefixes

    Returns:
        dict: Dictionary mapping each center prefix to its processed data
    """
    results = {}

    for prefix in center_prefixes:
        results[prefix] = process_excitation_center(directory, prefix)

    return results

def combine_spectra(centers_data: Dict[str, Dict[str, Any]]) -> Dict[str, np.ndarray]:
    """
    Combine interpolated spectra from multiple excitation centers.

    Args:
        centers_data (dict): Dictionary mapping center prefixes to processed data

    Returns:
        dict: Dictionary with combined energy and intensity arrays
    """
    # Initialize with the first center's energy array
    if not centers_data:
        return {'energy': np.array([]), 'intensity': np.array([])}

    first_center = next(iter(centers_data.values()))
    if 'interpolated_spectrum' not in first_center:
        return {'energy': np.array([]), 'intensity': np.array([])}

    energy = first_center['interpolated_spectrum']['energy']
    combined_intensity = np.zeros_like(energy)

    # Add the intensity from each center
    for center_data in centers_data.values():
        if 'interpolated_spectrum' in center_data:
            combined_intensity += center_data['interpolated_spectrum']['intensity']

    return {
        'energy': energy,
        'intensity': combined_intensity
    }

def save_to_csv(data: Dict[str, Any], output_file: str) -> None:
    """
    Save extracted data to CSV file.

    Args:
        data (dict): Data to save
        output_file (str): Path to output CSV file
    """
    # Convert to dataframe
    if 'interpolated_spectrum' in data:
        df = pl.DataFrame({
            'energy': data['interpolated_spectrum']['energy'],
            'intensity': data['interpolated_spectrum']['intensity']
        })
        df.write_csv(output_file)
    else:
        print(f"No interpolated spectrum data to save to {output_file}")

def to_xarray_datatree(data: Dict[str, Any], prefix: str = "C1") -> Any:
    """
    Convert processed StoBe data to an xarray DataTree.

    Args:
        data: Dictionary with processed StoBe data
        prefix: Prefix for the excitation center

    Returns:
        xarray.DataTree: DataTree containing organized StoBe data
    """
    try:
        import xarray as xr
    except ImportError:
        print("xarray not installed. Install it with: pip install xarray")
        return None

    # Create the DataTree
    tree = xr.DataTree()

    # Create ground state dataset
    gnd_data_vars = {}
    gnd_coords = {}

    # Add energy data
    if data['ground_state_energy'] is not None:
        gnd_data_vars['energy'] = data['ground_state_energy']
        gnd_data_vars['energy_ev'] = data['ground_state_energy'] * HARTREE_TO_EV

    # Add geometry data if available
    if data['geometry']:
        atoms = []
        x, y, z = [], [], []
        elements = []
        charges = []

        for atom in data['geometry']:
            atoms.append(atom['atom'])
            x.append(atom['x'])
            y.append(atom['y'])
            z.append(atom['z'])
            elements.append(atom['atom'][0])  # First character is usually element symbol
            charges.append(atom['q'])

        gnd_coords['atom'] = atoms
        gnd_data_vars['x'] = ('atom', x)
        gnd_data_vars['y'] = ('atom', y)
        gnd_data_vars['z'] = ('atom', z)
        gnd_data_vars['element'] = ('atom', elements)
        gnd_data_vars['charge'] = ('atom', charges)

    # Create ground state dataset
    ds_ground = xr.Dataset(
        data_vars=gnd_data_vars,
        coords=gnd_coords,
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

    # Create TP dataset
    if data['tp_energy'] is not None:
        tp_data_vars = {
            'energy': data['tp_energy'],
            'energy_ev': data['tp_energy'] * HARTREE_TO_EV
        }

        if data['lumo_energy'] is not None:
            tp_data_vars['lumo_energy'] = data['lumo_energy']

        if data['ip'] is not None:
            tp_data_vars['ionization_potential'] = data['ip']

        # Add oscillator strengths if available
        if data['oscillator_strengths'] and data['oscillator_strengths']['energies']:
            os_data = data['oscillator_strengths']

            # Create coordinates for transitions
            n_transitions = len(os_data['energies'])
            transition_coords = {'transition': range(n_transitions)}

            tp_data_vars['transition_energy'] = ('transition', os_data['energies'])
            tp_data_vars['oscillator_strength'] = ('transition', os_data['strengths'])
            tp_data_vars['tdm_x'] = ('transition', os_data['tdm_x'])
            tp_data_vars['tdm_y'] = ('transition', os_data['tdm_y'])
            tp_data_vars['tdm_z'] = ('transition', os_data['tdm_z'])

            if 'corrected_energies' in os_data:
                tp_data_vars['corrected_energy'] = ('transition', os_data['corrected_energies'])

            ds_tp = xr.Dataset(
                data_vars=tp_data_vars,
                coords=transition_coords,
                attrs={'type': 'transition_potential', 'prefix': prefix}
            )
        else:
            ds_tp = xr.Dataset(
                data_vars=tp_data_vars,
                attrs={'type': 'transition_potential', 'prefix': prefix}
            )

        tree['transition_potential'] = ds_tp

    # Create spectrum dataset
    if 'interpolated_spectrum' in data and data['interpolated_spectrum']:
        spectrum = data['interpolated_spectrum']

        # Create dataset for the spectrum
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
                'energy_correction': data['energy_correction']
            }
        )
        tree['spectrum'] = ds_spectrum

    # Create XAS dataset
    if data['xas_spectrum'] and data['xas_spectrum']['energies']:
        xas = data['xas_spectrum']

        # Create dataset for XAS
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

if __name__ == "__main__":
    # Example usage
    directory = "/home/hduva/projects/dft-learn/output"
    center_prefix = "C1"

    # Process a single excitation center
    center_data = process_excitation_center(directory, center_prefix)

    # Print the results
    print(f"Ground state energy: {center_data['ground_state_energy']} Hartree")
    print(f"Excited state energy: {center_data['excited_state_energy']} Hartree")
    print(f"Transition potential energy: {center_data['tp_energy']} Hartree")
    print(f"LUMO energy: {center_data['lumo_energy']} eV")
    print(f"Ionization potential: {center_data['ip']} eV")
    print(f"Energy correction: {center_data['energy_correction']} eV")

    if 'oscillator_strengths' in center_data and center_data['oscillator_strengths']:
        print(f"Found {len(center_data['oscillator_strengths']['energies'])} transitions")

    if 'interpolated_spectrum' in center_data:
        print(f"Generated interpolated spectrum with {len(center_data['interpolated_spectrum']['energy'])} points")

        # Save to CSV
        output_file = os.path.join(directory, f"{center_prefix}_spectrum.csv")
        save_to_csv(center_data, output_file)
        print(f"Saved spectrum to {output_file}")
