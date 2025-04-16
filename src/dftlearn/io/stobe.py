import re
from datetime import datetime
from pathlib import Path
from typing import Any

import polars as pl
import pandas as pd
import xarray as xr

def parse_geometry_line(line: str) -> dict[str, Any] | None:
    """
    Parse a single line from the geometry section.

    Args:
        line (str): Line containing atom geometry data

    Returns
    -------
        dict or None: Dictionary with atom information or None if parsing failed
    """
    # Don't try to parse header or separator lines
    if "----" in line or ("Atom" in line and "x" in line):
        return None

    # Simple but robust pattern matching - split by whitespace
    parts = line.strip().split()

    # Check if we have enough parts to make a valid atom entry
    if len(parts) >= 8:  # At least 8 elements for a geometry line
        try:
            # Extract the index and atom name considering the format "1) C01"
            idx_str = parts[0]
            atom_str = parts[1]

            # If the first part has a trailing ")", it contains the index
            if idx_str.endswith(")"):
                _ = idx_str.rstrip(")")
            else:
                # Handle case where index and atom are merged like "1)C01"
                match = re.match(r"(\d+)\)(.*)", idx_str)
                if match:
                    atom_str = match.group(2)

            return {
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
            }
        except (ValueError, IndexError):
            # Log diagnostic info if needed
            # print(f"Failed to parse line: {line}")
            return None

    return None


def extract_section(
    content: list[str], start_pattern: str, end_pattern: str
) -> list[str]:
    """
    Extract a section from the file content between start and end patterns.

    Args:
        content (list): List of lines from the file
        start_pattern (str): Regex pattern to match the start of the section
        end_pattern (str): Regex pattern to match the end of the section

    Returns
    -------
        list: Lines of the extracted section
    """
    section = []
    in_section = False
    for line in content:
        # Check if we're at the start of the section
        if not in_section and re.search(start_pattern, line):
            in_section = True
            section.append(line)
        # If we're already in the section
        elif in_section:
            section.append(line)
            # Check if we've reached the end of the section
            if re.search(end_pattern, line):
                break

    return section

def extract_geometry(lines: list[str]) -> pl.DataFrame:
    """
    Extract the geometry section from the GND file.
    """
    geometry_section = extract_section(
        lines,
        start_pattern=r"INPUT GEOMETRY \(input file\)",
        end_pattern=r"Smallest atom distance.*?=.*?\d+\.\d+",
    )
    geometry = pl.DataFrame(
        [
            parse_geometry_line(line)
            for line in geometry_section
            if line.strip() and parse_geometry_line(line) is not None
        ]
    )
    return geometry

def extract_basis_set(line):
    """
    Extract the basis set information from a line.

    Parameters
    ----------
    line : str
        A line from the basis set section of the GND file.

    Returns
    -------
    dict or None
        A dictionary with the basis set information or None if parsing failed.
    """
    # ensure the line is not the header or separator line
    if "I" in line or "II" in line or "III" in line:
        return None
    # parse the line by atom, and basis set
    parts = line.split()
    if len(parts) < 2:
        return None
    return {
        "atom": parts[1],
        "basis": " ".join(parts[3:]),
    }


def extract_basis_sets(lines: list[str]) -> pl.DataFrame:
    """
    Extract the basis sets from the GND file.

    Parameters
    ----------
    lines : list[str]
        The lines from the GND file.

    Returns
    -------
    pl.DataFrame
        A DataFrame containing the basis set information.
    """
    # Extract different basis set sections
    section_configs = [
        {
            "name": "aux",
            "start": r"I\)  AUXILIARY BASIS SETS",
            "end": r"II\)  ORBITAL BASIS SETS",
            "suffix": "",
        },
        {
            "name": "orbital",
            "start": r"II\)  ORBITAL BASIS SETS",
            "end": r" BASIS DIMENSIONS",
            "suffix": "_orbital",
        },
        {
            "name": "model_core",
            "start": r"III\)  MODEL POTENTIALS",
            "end": r" \(NEW\) SYMMETRIZATION INFORMATION",
            "suffix": "_model_core",
        },
    ]

    # Extract and process each section
    dfs = {}
    for config in section_configs:
        section = extract_section(
            lines,
            start_pattern=config["start"],
            end_pattern=config["end"],
        )

        # Filter and parse data
        data = [
            entry
            for line in section
            if line.strip() and (entry := extract_basis_set(line)) is not None
        ]

        # Create DataFrame or empty placeholder
        dfs[config["name"]] = (
            pl.DataFrame(data) if data else pl.DataFrame({"atom": [], "basis": []})
        )

    # Combine DataFrames using left joins to preserve all rows from orbital basis
    basis_sets = (
        dfs["orbital"]
        .join(dfs["aux"], on="atom", how="left", suffix=" auxiliary")
        .join(
            dfs["model_core"],
            on="atom",
            how="left",
            suffix=" model core potential",
        )
    )

    return basis_sets

def extract_energy(line: str) -> dict[str, Any] | None:
    """
    Extract energy information from a line.

    Args:
        line (str): Line containing energy data
        pattern (str): Regex pattern to match the energy data

    Returns
    -------
        dict or None: Dictionary with energy information or None if parsing failed
    """
    if "<Rho/" in line:
        return None
    # Energies have the pattern
    # \sKind of energy (H) = value
    pattern = r"([\w\s-]+?)\s*\((\w+)\)\s*=\s*([+-]?\d+\.\d+)"
    match = re.search(pattern, line)
    if match:
        energy_type = match.group(1).strip()  # Strip whitespace from the type
        energy_value = float(match.group(3)) * 27.2114  # convert from Hartree to eV
        return {
            "type": energy_type,
            "value": energy_value,
        }
    return None


def extract_energy_section(lines: list[str]) -> list[dict[str, Any]]:
    """
    Extract energy information from the energy section.

    Args:
        lines (list): Lines from the energy section

    Returns
    -------
        list: List of dictionaries with energy information
    """
    energy_section = extract_section(
        lines,
        start_pattern=r" FINAL ENERGY",
        end_pattern=r" Decomposition of",
    )
    return [
        entry
        for line in energy_section
        if line.strip() and (entry := extract_energy(line)) is not None
    ]

def extract_orbital_energies(lines: list[str]) -> pl.DataFrame:
    """
    Extract orbital energy information from the GND file.

    Parameters
    ----------
    lines : list[str]
        The lines from the GND file.

    Returns
    -------
    dict
        A dictionary containing the orbital energy information for alpha and beta spins.
    """
    orbital_section = extract_section(
        lines,
        start_pattern=r"ORBITAL ENERGIES \(ALL VIRTUALS INCLUDED\)",
        end_pattern=r" LISTING OF SPIN ALPHA ORBITALS ONLY",
    )

    # Skip header lines (2 lines for header, 1 for column names)
    data_lines = orbital_section[3:]

    # Initialize lists to store the parsed data
    orbitals = []
    alpha_occupations = []
    alpha_energies = []
    alpha_symmetries = []
    alpha_positions = []
    beta_occupations = []
    beta_energies = []
    beta_symmetries = []
    beta_positions = []

    # Parse each data line
    for line in data_lines:
        if not line.strip():  # Skip empty lines
            continue

        parts = line.split()
        if len(parts) >= 11:  # Check if we have enough parts for a complete line
            try:
                orbital_num = int(parts[0])
                orbitals.append(orbital_num)

                # Alpha spin data
                alpha_occupations.append(float(parts[1]))
                alpha_energies.append(float(parts[2]))
                alpha_symmetries.append(parts[3])
                alpha_positions.append(int(parts[5].strip("()")))

                # Beta spin data
                beta_occupations.append(float(parts[6]))
                beta_energies.append(float(parts[7]))
                beta_symmetries.append(parts[8])
                beta_positions.append(int(parts[10].strip("()")))
            except (ValueError, IndexError):
                # Skip lines that don't match the expected format
                continue

    return pl.DataFrame(
        {
            "orbital": orbitals,
            "alpha_occupation": alpha_occupations,
            "alpha_energy": alpha_energies,
            "alpha_symmetry": alpha_symmetries,
            "alpha_position": alpha_positions,
            "beta_occupation": beta_occupations,
            "beta_energy": beta_energies,
            "beta_symmetry": beta_symmetries,
            "beta_position": beta_positions,
        }
    )

def read_stobe_ground(
    file: str | Path,
    author: str | None = None,
    comment: str | None = None,
    save_checkpoint=True,
) -> xr.DataTree:
    """
    Read the GND file and extract relevant data into an xarray DataTree.

    This function extracts geometry, basis sets, energy, and spin information
    from the GND file and organizes it into a structured xarray DataTree.


    Parameters
    ----------
    file : str or Path
        The path to the GND file.

    author : str
        The name of the author.

    comment : str
        A comment regarding the data extraction.

    save_checkpoint : bool, optional
        Whether to save a checkpoint of the extraction process (default is True).

    Returns
    -------
    xr.Dataset
        A dataset containing the extracted data.

    Details
    -------
    This loads the GND file, extracts the geometry, basis sets, energy, and occupation.
    This creates three groups of data in the DataTree:
    - GEOMETRY: Contains the geometry information of the atoms.
    - BASIS: Contains the basis set information for the atoms.
    - GND: Contains the ground state energy and spin information.
        - energy: Contains the ground state energy information.
        - spin: Contains the spin occupation information.

    Example
    -------
    >>> dt = read_stobe_ground(
    ...     "path/to/gnd/file.gnd", author="John Doe", comment="Test extraction"
    ... )
    >>> print(dt)
    <xarray.DataTree>
    <GEOMETRY>: ...
    <BASIS>: ...
    <GND>: <energy>: ... <spin>: ...
    """
    file = Path(file)
    with file.open("r") as f:
        lines = f.readlines()
    geometry_section = extract_geometry(lines)
    basis_section = extract_basis_sets(lines)
    ground_energy_section = extract_energy_section(lines)
    gnd_spin_section = extract_orbital_energies(lines)
    #  create a dataset with each energy type as a separate variable
    geometry_ds = xr.Dataset(
        {
            "x": (("atom"), geometry_section["x"].to_numpy()),
            "y": (("atom"), geometry_section["y"].to_numpy()),
            "z": (("atom"), geometry_section["z"].to_numpy()),
            "q": (("atom"), geometry_section["q"].to_numpy()),
            "nuc": (("atom"), geometry_section["nuc"].to_numpy()),
            "mass": (("atom"), geometry_section["mass"].to_numpy()),
            "neq": (("atom"), geometry_section["neq"].to_numpy()),
            "grid": (("atom"), geometry_section["grid"].to_numpy()),
            "grp": (("atom"), geometry_section["grp"].to_numpy()),
        },
        coords={
            "atom": geometry_section["atom"].to_numpy(),
        },
    )
    basis_ds = xr.Dataset(
        {
            "orbital_basis": (("atom"), basis_section["basis"].to_numpy()),
            "auxiliary_basis": (("atom"), basis_section["basis auxiliary"].to_numpy()),
            "model_core_basis": (
                ("atom"),
                basis_section["basis model core potential"].to_numpy(),
            ),
        },
        coords={
            "atom": basis_section["atom"].to_numpy(),
        },
    )
    energy_ds = xr.Dataset(
        {
            entry["type"].replace(" ", "_").lower(): ([], entry["value"])
            for entry in ground_energy_section
        },
    )
    # extract molecular orbitals
    orbital_section = extract_section(
        lines,
        start_pattern=r" Alpha occupation:",
        end_pattern=r" Beta occupation:",
    )
    # get the alpha occupation numbers for the HOMO

    homo_orbital = int(orbital_section[0].split(":")[1].strip())
    lumo_orbital = homo_orbital + 1
    homo_energy = (
        gnd_spin_section["alpha_energy"].to_numpy()[homo_orbital]
        if homo_orbital is not None
        else None
    )
    lumo_energy = (
        gnd_spin_section["alpha_energy"].to_numpy()[lumo_orbital]
        if lumo_orbital is not None
        else None
    )
    energy_ds["lumo-energy"] = lumo_energy
    energy_ds["homo-energy"] = homo_energy
    # search for a core hole with .5 occupation, if it dow not exit then don't add it
    # do the dextract_basis_set
    core_hole_orbital = gnd_spin_section.filter(pl.col("alpha_occupation").eq(0.5))[
        "orbital"
    ]

    spin_ds = xr.Dataset(
        {
            "homo": homo_orbital,
            "lumo": lumo_orbital,
            "alpha_occupation": (
                ("orbital"),
                gnd_spin_section["alpha_occupation"].to_numpy(),
            ),
            "alpha_energy": (("orbital"), gnd_spin_section["alpha_energy"].to_numpy()),
            "alpha_symmetry": (
                ("orbital"),
                gnd_spin_section["alpha_symmetry"].to_numpy(),
            ),
            "alpha_position": (
                ("orbital"),
                gnd_spin_section["alpha_position"].to_numpy(),
            ),
            "beta_occupation": (
                ("orbital"),
                gnd_spin_section["beta_occupation"].to_numpy(),
            ),
            "beta_energy": (("orbital"), gnd_spin_section["beta_energy"].to_numpy()),
            "beta_symmetry": (
                ("orbital"),
                gnd_spin_section["beta_symmetry"].to_numpy(),
            ),
            "beta_position": (
                ("orbital"),
                gnd_spin_section["beta_position"].to_numpy(),
            ),
        },
        coords={
            "orbital": gnd_spin_section["orbital"].to_numpy(),
        },
    )
    if not core_hole_orbital.is_empty():
        spin_ds["core-hole"] = core_hole_orbital.to_numpy()[0]
        energy_ds["core-hole-energy"] = gnd_spin_section["alpha_energy"].to_numpy()[
            core_hole_orbital.to_numpy()[0]
        ]
    # combine to a single datatree
    # combine to a single datatree
    dt = xr.DataTree()
    gnd = xr.DataTree()
    dt["GEOMETRY"] = geometry_ds
    dt["BASIS"] = basis_ds
    gnd["energy"] = energy_ds
    gnd["spin"] = spin_ds
    dt["GND"] = gnd

    # add metadata
    dt.attrs["date"] = datetime.now().strftime("%Y-%m-%d")
    dt.attrs["source"] = str(file)
    dt.attrs["description"] = "GND file extracted data"
    dt.attrs["version"] = "1.0"
    dt.attrs["additional_info"] = "This dataset contains extracted data from GND files."

    if comment:
        dt.attrs["comment"] = comment
    if author:
        dt.attrs["author"] = author

    if save_checkpoint:
        dt.to_netcdf(
            Path(file).parent / "gnd.nc",
            mode="w",
            engine="h5netcdf",
        )
    return dt

def read_stobe_excited(
    file: str | Path,
    author: str | None = None,
    comment: str | None = None,
    *,
    save_checkpoint: bool = True,
) -> xr.DataTree:
    """
    Read the EXC file and extract relevant data into an xarray DataTree.

    This function extracts geometry, basis sets, energy, and spin information
    from the EXC file and organizes it into a structured xarray DataTree.


    Parameters
    ----------
    file : str or Path
        The path to the EXC file.

    author : str
        The name of the author.

    comment : str
        A comment regarding the data extraction.

    save_checkpoint : bool, optional
        Whether to save a checkpoint of the extraction process (default is True).

    Returns
    -------
    xr.Dataset
        A dataset containing the extracted data.

    Details
    -------
    This loads the EXC file, extracts the geometry, basis sets, energy, and occupation.
    This creates three groups of data in the DataTree:
    - GEOMETRY: Contains the geometry information of the atoms.
    - BASIS: Contains the basis set information for the atoms.
    - GND: Contains the ground state energy and spin information.
        - energy: Contains the ground state energy information.
        - spin: Contains the spin occupation information.

    Example
    -------
    >>> dt = read_stobe_ground(
    ...     "path/to/gnd/file.gnd", author="John Doe", comment="Test extraction"
    ... )
    >>> print(dt)
    <xarray.DataTree>
    <GEOMETRY>: ...
    <BASIS>: ...
    <EXC>: <energy>: ... <spin>: ...
    """
    bad_name_data_tree = read_stobe_ground(file, author, comment, save_checkpoint=False)
    # rename for excited state
    exc = bad_name_data_tree["GND"].copy()
    exc.name = "EXC"
    del bad_name_data_tree["GND"]
    bad_name_data_tree["EXC"] = exc

    # up

    bad_name_data_tree.attrs["description"] = "EXC file extracted data"

    if save_checkpoint:
        bad_name_data_tree.to_netcdf(
            Path(file).parent / "exc.nc",
            mode="w",
            engine="h5netcdf",
        )
    return bad_name_data_tree

def process_transition_line(line):
    """
    Process a line from the transition section of the EXC file.
    """
    parts = line.split()
    if len(parts) < 9:
        return None
    return {
        "energy": float(parts[2]),
        "oscillator_strength": float(parts[3]),
        "oslx": float(parts[4]),
        "osly": float(parts[5]),
        "oslz": float(parts[6]),
        "osc_r2": float(parts[7]),
        "<r2>": float(parts[8]),
    }


def extract_transitions(lines: list[str]) -> pl.DataFrame:
    """
    Extract the transitions from the EXC file.

    Parameters
    ----------
    lines : list[str]
        The lines from the EXC file.

    Returns
    -------
    pl.DataFrame
        A DataFrame containing the transition information.
    """
    transition_section = extract_section(
        lines,
        start_pattern=r"\s+E \(eV\)\s+OSCL\s+oslx\s+osly\s+oslz\s+osc\(r2\)\s+<r2> ",
        end_pattern=r"\s+SUM\s+[\d.]+\s+[\d.]+",
    )
    transitions = [
        process_transition_line(line)
        for line in transition_section
        if line.strip() and process_transition_line(line) is not None
    ]
    return pl.DataFrame(transitions)


def read_stobe_tp(
    file: str | Path,
    author: str | None = None,
    comment: str | None = None,
    *,
    save_checkpoint: bool = True,
) -> xr.DataTree:
    """
    Read the TP file and extract relevant data into an xarray DataTree.

    This function extracts geometry, basis sets, energy, spin information,
    and transitions from the TP file and organizes it into a structured xarray DataTree.

    Parameters
    ----------
    file : str or Path
        The path to the TP file.

    author : str
        The name of the author.

    comment : str
        A comment regarding the data extraction.

    save_checkpoint : bool, optional
        Whether to save a checkpoint of the extraction process (default is True).

    Returns
    -------
    xr.Dataset
        A dataset containing the extracted data.

    Details
    -------
    This loads the TP file, extracts the geometry, basis sets, energy, and occupation.
    This creates three groups of data in the DataTree:
    - GEOMETRY: Contains the geometry information of the atoms.
    - BASIS: Contains the basis set information for the atoms.
    - TP: Contains the excited state energy and spin information.
        - energy: Contains the excited state energy information.
        - spin: Contains the spin occupation information.
        - transitions: Contains the transition information.

    Example
    -------
    >>> dt = read_stobe_tp(
    ...     "path/to/tp/file.tp", author="John Doe", comment="Test extraction"
    ... )
    >>> print(dt)
    <xarray.DataTree>
    <GEOMETRY>: ...
    <BASIS>: ...
    <TP>: <energy>: ... <spin>: ... <transitions>: ...
    """
    bad_name_data_tree = read_stobe_excited(
        file, author, comment, save_checkpoint=False
    )
    # rename for excited state
    tp = bad_name_data_tree["EXC"].copy()
    tp.name = "TP"
    del bad_name_data_tree["EXC"]
    # add transitions
    # First, read the file content
    with Path(file).open("r") as f:
        lines = f.readlines()

    # extract the core hole information
    core_hole_section = extract_section(
        lines,
        start_pattern=r"Orbital energy core hole",
        end_pattern=r"Ionization potential",
    )
    core_hole_energy = (
        core_hole_section[0].split("=")[1].split("(")[1].split("e")[0].strip()
    )
    ridgid_shift = core_hole_section[1].split("=")[1].split("e")[0].strip()
    ionization_potential = core_hole_section[2].split("=")[1].split("e")[0].strip()
    tp["energy"]["core_hole_energy"] = float(core_hole_energy)
    tp["energy"]["ridgid_shift"] = float(ridgid_shift)
    tp["energy"]["ionization_potential"] = float(ionization_potential)

    # Extract transitions data
    transitions_df = extract_transitions(lines)

    # Create a proper xarray Dataset for transitions
    transitions_ds = xr.Dataset(
        {
            "oscillator_strength": (
                ("energy"),
                transitions_df["oscillator_strength"].to_numpy(),
            ),
            "oslx": (("energy"), transitions_df["oslx"].to_numpy()),
            "osly": (("energy"), transitions_df["osly"].to_numpy()),
            "oslz": (("energy"), transitions_df["oslz"].to_numpy()),
            "osc_r2": (("energy"), transitions_df["osc_r2"].to_numpy()),
            "r2": (("energy"), transitions_df["<r2>"].to_numpy()),
        },
        coords={
            "energy": transitions_df["energy"].to_numpy(),
        },
    )

    # Add the transitions as a separate group in the TP DataTree
    tp["transitions"] = transitions_ds
    bad_name_data_tree["TP"] = tp

    # Update description
    bad_name_data_tree.attrs["description"] = "TP file extracted data"

    if save_checkpoint:
        bad_name_data_tree.to_netcdf(
            Path(file).parent / "tp.nc",
            mode="w",
            engine="h5netcdf",
        )
    return bad_name_data_tree

def read_stobe(
    output_directory: str | Path,
    author: str | None = None,
    comment: str | None = None,
    *,
    save_checkpoint: bool = True,
) -> xr.DataTree:
    """
    Read primary output files and extract relevant data into an xarray DataTree.

    This function extracts geometry, basis sets, energy, and spin information
    from the GND file and organizes it into a structured xarray DataTree.

    Parameters
    ----------
    output_directory : str or Path
        The path to the output directory containing the GND and EXC files.

    author : str
        The name of the author.

    comment : str
        A comment regarding the data extraction.

    save_checkpoint : bool, optional
        Whether to save a checkpoint of the extraction process (default is True).

    Returns
    -------
    xr.Dataset
        A dataset containing the extracted data.
    """
    output_directory = Path(output_directory)
    gnd_path = output_directory / f"{output_directory.name}gnd.out"
    exc_path = output_directory / f"{output_directory.name}exc.out"
    tp_path = output_directory / f"{output_directory.name}tp.out"

    # Read the GND file
    gnd_section = read_stobe_ground(
        gnd_path, author=author, comment=comment, save_checkpoint=False
    )
    # Read the EXC file
    exc_section = read_stobe_excited(
        exc_path, author=author, comment=comment, save_checkpoint=False
    )
    # Read the TP file
    tp_section = read_stobe_tp(
        tp_path, author=author, comment=comment, save_checkpoint=False
    )

    # Combine the sections into a single DataTree
    combined_data_tree = xr.DataTree()
    combined_data_tree["GEOMETRY"] = gnd_section["GEOMETRY"]
    combined_data_tree["BASIS"] = gnd_section["BASIS"]
    combined_data_tree["GND"] = gnd_section["GND"]
    combined_data_tree["EXC"] = exc_section["EXC"]
    combined_data_tree["TP"] = tp_section["TP"]

    # Add metadata to the DataTree
    combined_data_tree.attrs["date"] = datetime.now().strftime("%Y-%m-%d")
    combined_data_tree.attrs["description"] = (
        "Combined data from GND and EXC files extracted data"
    )
    if comment:
        combined_data_tree.attrs["comment"] = comment
    if author:
        combined_data_tree.attrs["author"] = author

    if save_checkpoint:
        combined_data_tree.to_netcdf(
            output_directory / "out.nc",
            mode="w",
            engine="h5netcdf",
        )
    return combined_data_tree


def read_calculations(
    directory: str | Path,
    calculation_backend="stobe",
    author: str | None = None,
    comment: str | None = None,
    *,
    save_checkpoint: bool = True,
) -> xr.DataTree:
    """
    Extract all the data from several calculations and combine them into a single DataTree.

    Data from different calculations (identified by subdirectories) are concatenated
    along a new 'excitation_atom' dimension.

    Parameters
    ----------
    directory : str or Path
        The parent directory containing subdirectories for each calculation (e.g., C1, C2).
    calculation_backend : str, optional
        The backend used for the calculations (default is "stobe").
    author : str, optional
        The name of the author.
    comment : str, optional
        A comment regarding the data extraction.
    save_checkpoint : bool, optional
        Whether to save intermediate checkpoints for each calculation (default is True).

    Returns
    -------
    xr.DataTree
        A DataTree containing the combined data from all calculations.
    """
    directory = Path(directory)
    excitation_atoms = [e.name for e in directory.iterdir() if e.is_dir()]
    try:
        excitation_atoms.sort(key=lambda x: int(re.split(r"([A-Za-z]+)(\d+)", x)[2]))
    except (IndexError, ValueError):
        excitation_atoms.sort()

    if not excitation_atoms:
        error_message = f"No calculation subdirectories found in {directory}"
        raise FileNotFoundError(error_message)

    # Lists to store datasets from each calculation
    geometry_list, basis_list = [], []
    gnd_energy_list, gnd_spin_list = [], []
    exc_energy_list, exc_spin_list = [], []
    tp_energy_list, tp_spin_list, tp_transitions_list = [], [], []

    for excitation_atom in excitation_atoms:
        atom_directory = directory / excitation_atom
        # check if the directory contains the files
        match calculation_backend:
            case "stobe":
                # Pass save_checkpoint to read_stobe for individual file saving
                data_tree = read_stobe(
                    atom_directory,
                    author=author,
                    comment=comment,
                    save_checkpoint=save_checkpoint,  # Control individual saving here
                )
            case _:
                error_message = f"Unknown calculation backend: {calculation_backend}"
                raise ValueError(error_message)

        # Append datasets to lists
        geometry_list.append(data_tree["GEOMETRY"].ds)
        basis_list.append(data_tree["BASIS"].ds)
        gnd_energy_list.append(data_tree["GND/energy"].ds)
        gnd_spin_list.append(data_tree["GND/spin"].ds)
        exc_energy_list.append(data_tree["EXC/energy"].ds)
        exc_spin_list.append(data_tree["EXC/spin"].ds)
        tp_energy_list.append(data_tree["TP/energy"].ds)
        tp_spin_list.append(data_tree["TP/spin"].ds)
        # Handle case where TP/transitions might be missing
        if "transitions" in data_tree["TP"]:
            tp_transitions_list.append(data_tree["TP/transitions"].ds)
        else:
            pass

    # Create the coordinate for the new dimension
    excitation_coord = pd.Index(excitation_atoms, name="excitation_atom")

    # Concatenate datasets along the new 'excitation_atom' dimension
    # Use DataTree.from_dict for cleaner construction
    combined_data = {
        "GEOMETRY": xr.concat(geometry_list, dim=excitation_coord),
        "BASIS": xr.concat(basis_list, dim=excitation_coord),
        "GND/energy": xr.concat(gnd_energy_list, dim=excitation_coord),
        "GND/spin": xr.concat(gnd_spin_list, dim=excitation_coord),
        "EXC/energy": xr.concat(exc_energy_list, dim=excitation_coord),
        "EXC/spin": xr.concat(exc_spin_list, dim=excitation_coord),
        "TP/energy": xr.concat(tp_energy_list, dim=excitation_coord),
        "TP/spin": xr.concat(tp_spin_list, dim=excitation_coord),
    }
    if tp_transitions_list:  # Only add transitions if they were found
        combined_data["TP/transitions"] = xr.concat(
            tp_transitions_list, dim=excitation_coord
        )

    combined_dt = xr.DataTree.from_dict(combined_data)

    # Add metadata
    combined_dt.attrs["date"] = datetime.now().strftime("%Y-%m-%d")
    combined_dt.attrs["description"] = (
        "Combined data from several calculations extracted data"
    )
    combined_dt.attrs["source_directory"] = str(directory)
    combined_dt.attrs["version"] = "1.1"  # Increment version due to structure change
    if comment:
        combined_dt.attrs["comment"] = comment
    if author:
        combined_dt.attrs["author"] = author

    if save_checkpoint:
        final_save_path = directory / "combined_calculations.nc"
        combined_dt.to_netcdf(final_save_path, mode="w", engine="h5netcdf")
        print(f"Saved combined DataTree to {final_save_path}")

    return combined_dt
