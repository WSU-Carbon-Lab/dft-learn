"""I/O functions for StoBe output files."""

import os
import re
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd

CPU_COUNT = os.cpu_count() or 1


class EnergyPattern(Enum):
    """
    Enum to hold patterns for parsing energy values from StoBe output.
    """

    TOTAL_ENERGY_H = (
        "Total energy   (H)",
        r"Total energy\s+\(H\)\s*=\s*([-+]?[0-9]*\.?[0-9]+)",
        ["Total energy (H)"],
    )
    NUC_NUC_ENERGY_H = (
        "Nuc-nuc energy (H)",
        r"Nuc-nuc energy\s+\(H\)\s*=\s*([-+]?[0-9]*\.?[0-9]+)",
        ["Nuc-nuc energy (H)"],
    )
    EL_NUC_ENERGY_H = (
        "El-nuc energy  (H)",
        r"El-nuc energy\s+\(H\)\s*=\s*([-+]?[0-9]*\.?[0-9]+)",
        ["El-nuc energy (H)"],
    )
    KINETIC_ENERGY_H = (
        "Kinetic energy (H)",
        r"Kinetic energy\s+\(H\)\s*=\s*([-+]?[0-9]*\.?[0-9]+)",
        ["Kinetic energy (H)"],
    )
    COULOMB_ENERGY_H = (
        "Coulomb energy (H)",
        r"Coulomb energy\s+\(H\)\s*=\s*([-+]?[0-9]*\.?[0-9]+)",
        ["Coulomb energy (H)"],
    )
    EX_COR_ENERGY_H = (
        "Ex-cor energy  (H)",
        r"Ex-cor energy\s+\(H\)\s*=\s*([-+]?[0-9]*\.?[0-9]+)",
        ["Ex-cor energy (H)"],
    )
    ORBITAL_ENERGY_CORE_HOLE = (
        "Orbital energy core hole",
        r"Orbital energy core hole\s*=\s*([-+]?[0-9]*\.?[0-9]+)\s*H\s*\(\s*([-+]?[0-9]*\.?[0-9]+)\s*eV\s*\)",  # noqa: E501
        ["Orbital energy core hole (H)", "Orbital energy core hole (eV)"],
    )
    RIGID_SPECTRAL_SHIFT_EV = (
        "Rigid spectral shift",
        r"Rigid spectral shift\s*=\s*([-+]?[0-9]*\.?[0-9]+)\s*eV",
        ["Rigid spectral shift (eV)"],
    )
    IONIZATION_POTENTIAL_EV = (
        "Ionization potential",
        r"Ionization potential\s*=\s*([-+]?[0-9]*\.?[0-9]+)\s*eV",
        ["Ionization potential (eV)"],
    )

    def __init__(self, search_str, regex, keys):
        self.search_str = search_str
        self.regex = re.compile(regex)
        self.keys = keys


def _energy_parser(lines: list[str]) -> dict:
    """Extract energy information from lines of a StoBe output file."""
    energies: dict[str, float] = {
        "Total energy (H)": 0.0,
        "Nuc-nuc energy (H)": 0.0,
        "El-nuc energy (H)": 0.0,
        "Kinetic energy (H)": 0.0,
        "Coulomb energy (H)": 0.0,
        "Ex-cor energy (H)": 0.0,
        "Orbital energy core hole (H)": 0.0,
        "Orbital energy core hole (eV)": 0.0,
        "Rigid spectral shift (eV)": 0.0,
        "Ionization potential (eV)": 0.0,
    }

    for line in lines:
        for pattern in EnergyPattern:
            if pattern.search_str in line:
                match = pattern.regex.search(line)
                if match:
                    for i, key in enumerate(pattern.keys):
                        energies[key] = float(match.group(i + 1))
                    break  # Move to the next line once a pattern is matched

    return energies


def _basis_parser(line: str) -> tuple[str | None, str | None]:
    """Extract basis set information from a line of a StoBe output file."""
    line = line.strip()

    # Only process lines that start with "Atom " and contain a colon
    if line.startswith("Atom ") and ":" in line:
        try:
            parts = line.split(":", 1)
            if len(parts) == 2:
                # Extract atom identifier (everything after "Atom " and before ":")
                atom_part = parts[0].replace("Atom ", "").strip()
                basis_part = parts[1].strip()

                # Key fix: Only accept atom identifiers that contain letters
                # (element symbols)
                # This excludes pure numbers like "1", "2", "3" from the
                # exchange/correlation section
                # and only accepts identifiers like "C1", "Cu1", "N1", etc.
            if atom_part and re.match(r"^[A-Za-z]+\d*$", atom_part):
                return atom_part, basis_part
        except ValueError as e:
            # This might occur if string operations fail unexpectedly.
            print(
                f"Warning: Could not parse basis line due to ValueError: {line.strip()} -> {e}"  # noqa: E501
            )
        except IndexError as e:
            # This might occur if splitting the line doesn't produce expected parts.
            print(
                f"Warning: Could not parse basis line due to IndexError: {line.strip()} -> {e}"  # noqa: E501
            )

    return None, None


@dataclass
class OutputData:
    """
    Dataclass to hold extracted DataFrames from StoBe output.
    """

    df_alpha: pd.DataFrame
    df_beta: pd.DataFrame
    df_auxiliary: pd.DataFrame
    df_orbital: pd.DataFrame
    df_model: pd.DataFrame
    df_energies: pd.DataFrame
    df_xray_transitions: pd.DataFrame
    df_atomic_coordinates: pd.DataFrame

    @property
    def combined(self) -> pd.DataFrame:
        """Combine auxiliary, orbital, and model DataFrames."""
        return self.df_auxiliary.merge(
            self.df_orbital, on=["Atom", "Originating File"], how="outer"
        ).merge(self.df_model, on=["Atom", "Originating File"], how="outer")


def _concatenate_outputs(outputs: list[OutputData]) -> OutputData:
    """Concatenate multiple OutputData instances into one."""
    return OutputData(
        pd.concat([o.df_alpha for o in outputs], ignore_index=True),
        pd.concat([o.df_beta for o in outputs], ignore_index=True),
        pd.concat([o.df_auxiliary for o in outputs], ignore_index=True),
        pd.concat([o.df_orbital for o in outputs], ignore_index=True),
        pd.concat([o.df_model for o in outputs], ignore_index=True),
        pd.concat([o.df_energies for o in outputs], ignore_index=True),
        pd.concat([o.df_xray_transitions for o in outputs], ignore_index=True),
        pd.concat([o.df_atomic_coordinates for o in outputs], ignore_index=True),
    )


def _add_filename(output: OutputData, file_name: str) -> OutputData:
    for df in [
        output.df_auxiliary,
        output.df_orbital,
        output.df_model,
        output.df_alpha,
        output.df_beta,
        output.df_energies,
        output.df_xray_transitions,
        output.df_atomic_coordinates,
    ]:
        df["Originating File"] = file_name
    return output


def _output(file_path: str | Path, originating_atom: str) -> OutputData:
    """
    Read a StoBe output file and extract various pieces of information into DataFrames.

    Parameters
    ----------
    file_path : str | Path
        Path to the StoBe output file.
    originating_atom : str
        Identifier for the originating atom (e.g., "C1", "N1").
    """
    data_alpha = []
    data_beta = []
    auxiliary_basis = []
    orbital_basis = []
    model_potential = []
    xray_transitions = []
    atomic_coordinates = []
    first_xray_energy: float | None = None
    file_path = Path(file_path)
    with file_path.open() as file:
        lines = file.readlines()

    energies = _energy_parser(lines)
    energies["Atom"] = originating_atom
    match file_path.name:
        case name if "tp.out" in name:
            energies["Calculation Type"] = "TP"
        case name if "gnd.out" in name:
            energies["Calculation Type"] = "GND"
        case name if "exc.out" in name:
            energies["Calculation Type"] = "EXC"
        case _:
            energies["Calculation Type"] = None

    start_index: int | None = None
    end_index: int | None = None
    xray_start: bool = False
    atomic_start_index: int | None = None
    atomic_end_index: int | None = None
    current_section: str | None = None

    for i, line in enumerate(lines):
        if "         Spin alpha                              Spin beta" in line:
            start_index = i + 2
        elif " IV)" in line:
            end_index = i
        elif "I)  AUXILIARY BASIS SETS" in line:
            current_section = "auxiliary"
        elif "II)  ORBITAL BASIS SETS" in line:
            current_section = "orbital"
        elif "III)  MODEL POTENTIALS" in line:
            current_section = "model"
        elif current_section in ["auxiliary", "orbital", "model"]:
            # Check for section end markers or next section start
            if (
                "BASIS DIMENSIONS" in line
                or "IV)" in line
                or "WARNING!" in line
                or line.strip() == ""
                or (current_section == "auxiliary" and "II)" in line)
                or (current_section == "orbital" and "III)" in line)
            ):
                # Reset section if we hit a boundary
                if "II)" in line:
                    current_section = "orbital"
                elif "III)" in line:
                    current_section = "model"
                elif "BASIS DIMENSIONS" in line or "IV)" in line or "WARNING!" in line:
                    current_section = None
            else:
                # Only parse lines that look like atom definitions
                atom, basis = _basis_parser(line)
                if atom and basis:
                    if current_section == "auxiliary":
                        auxiliary_basis.append([atom, basis])
                    elif current_section == "orbital":
                        orbital_basis.append([atom, basis])
                    elif current_section == "model":
                        model_potential.append([atom, basis])
        elif (
            "           E (eV)   OSCL       oslx       osly       oslz         osc(r2)       <r2>"  # noqa: E501
            in line
        ):
            xray_start = True
        elif xray_start and "-----" in line:
            continue
        elif xray_start and line.strip() and line.startswith(" #"):
            try:
                index = int(line[2:6].strip())
                e_ev = float(line[7:17].strip())
                oscl = float(line[18:28].strip())
                oslx = float(line[29:39].strip())
                osly = float(line[40:50].strip())
                oslz = float(line[51:61].strip())
                osc_r2 = float(line[62:72].strip())
                r2 = float(line[73:].strip())

                # Capture the first X-ray transition energy for LUMO_En
                if index == 1 and first_xray_energy is None:
                    first_xray_energy = e_ev

                xray_transitions.append(
                    {
                        "Index": index,
                        "E": e_ev,
                        "OS": oscl,
                        "osx": oslx,
                        "osy": osly,
                        "osz": oslz,
                        "osc(r2)": osc_r2,
                        "<r2>": r2,
                    }
                )
            except ValueError as e:
                print(f"Error parsing line: {line}\nError: {e}")
        elif "Single image calculation (Angstrom):" in line:
            atomic_start_index = i + 3
        elif "Smallest atom distance" in line:
            atomic_end_index = i

    if start_index is not None and end_index is not None:
        for line in lines[start_index:end_index]:
            if line.strip() == "":
                continue
            components = [x for x in line.split() if x.strip()]
            if len(components) >= 9:
                mo_index_alpha, occup_alpha, energy_alpha, sym_alpha = components[:4]
                mo_index_beta, occup_beta, energy_beta, sym_beta = components[5:9]
                mo_index_alpha = mo_index_alpha.strip(")")
                mo_index_beta = mo_index_beta.strip(")")
                data_alpha.append(
                    {
                        "MO_Index": mo_index_alpha,
                        "Occup.": occup_alpha,
                        "Energy(eV)": energy_alpha,
                        "Sym.": sym_alpha,
                    }
                )
                data_beta.append(
                    {
                        "MO_Index": mo_index_beta,
                        "Occup.": occup_beta,
                        "Energy(eV)": energy_beta,
                        "Sym.": sym_beta,
                    }
                )

    if atomic_start_index is not None and atomic_end_index is not None:
        atomic_coordinates_lines = lines[atomic_start_index:atomic_end_index]
        for line in atomic_coordinates_lines:
            if line.strip() and not any(
                col in line
                for col in [
                    "Atom",
                    "x",
                    "y",
                    "z",
                    "q",
                    "nuc",
                    "mass",
                    "neq",
                    "grid",
                    "grp",
                ]
            ):  # Skip empty lines and the header
                split_line = line.split()
                if len(split_line) >= 11:
                    atom_info = split_line[1]  # Use the atom type and number
                    atomic_coordinates.append([atom_info, *split_line[2:11]])

    # Add the first X-ray transition energy to the energies dictionary
    energies["LUMO_En"] = first_xray_energy

    df_alpha = pd.DataFrame(data_alpha)
    df_beta = pd.DataFrame(data_beta)
    df_auxiliary = pd.DataFrame(
        auxiliary_basis, columns=pd.Index(["Atom", "Auxiliary Basis"])
    )
    df_orbital = pd.DataFrame(
        orbital_basis, columns=pd.Index(["Atom", "Orbital Basis"])
    )
    df_model = pd.DataFrame(
        model_potential, columns=pd.Index(["Atom", "Model Potential"])
    )
    df_energies = pd.DataFrame([energies])
    df_xray_transitions = pd.DataFrame(xray_transitions)
    df_atomic_coordinates = pd.DataFrame(
        atomic_coordinates,
        columns=pd.Index(
            ["Atom", "x", "y", "z", "q", "nuc", "mass", "neq", "grid", "grp"]
        ),
    )

    numeric_columns = ["x", "y", "z", "q", "nuc", "mass", "neq", "grid", "grp"]
    df_atomic_coordinates[numeric_columns] = df_atomic_coordinates[
        numeric_columns
    ].apply(pd.to_numeric, errors="coerce")

    return OutputData(
        df_alpha,
        df_beta,
        df_auxiliary,
        df_orbital,
        df_model,
        df_energies,
        df_xray_transitions,
        df_atomic_coordinates,
    )


def process_file(file_path: str | Path) -> OutputData | None:
    """
    Extract data from a single StoBe output file.
    """
    file_path = Path(file_path)
    try:
        originating_atom = extract_atom_id(file_path)
        output = _output(file_path, originating_atom)

        file_name = file_path.name
        output = _add_filename(output, file_name)
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        return None

    return output


def extract_atom_id(file_name: Path) -> str:
    """Extract atom identifier from filename."""
    base = Path(file_name).name.lower()
    match = re.match(r"([a-z]+[0-9]+)", base)
    return match.group(1).upper() if match else base.split(".")[0].upper()


def sort_dataframe_naturally(df, column):
    """Sorts a DataFrame naturally by the specified column."""
    try:
        from natsort import natsorted

        df[column] = df[column].astype(str)
        sorted_index = natsorted(df[column].tolist())
        df = df.set_index(column).loc[sorted_index].reset_index()
    except ImportError:
        print("Warning: natsort not available, using standard sort")
        df = df.sort_values(column)
    return df


def process_calculation_type(
    directory: str | Path,
    calc_type: Literal["GND", "EXC", "TP", "NEXAFS"],
    verbose=True,
    n_jobs=CPU_COUNT,
):
    """
    Process StoBe output files for a specific calculation type.

    Parameters
    ----------
    directory : str | Path
        Path to the base directory containing the calculation type folder.
    calc_type : Literal["GND", "EXC", "TP", "NEXAFS"]
        The calculation type to process (e.g., 'GND', 'EXC', 'TP', 'NEXAFS').
    verbose : bool
        Whether to print progress information.
    n_jobs : int, optional
        Number of parallel processes to use.

    Returns
    -------
    OutputData
        A tuple of DataFrames: (combined, energies, alpha, beta, xray, coords)
    """
    folder_path = Path(directory) / calc_type
    if not folder_path.is_dir():
        if verbose:
            print(
                f"Directory not found for calculation type '{calc_type}': {folder_path}"
            )
        return OutputData(
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
        )

    suffix = f"{calc_type.lower()}.out"
    if calc_type == "NEXAFS":  # NEXAFS is a special case
        suffix = "nexafs.out"

    file_paths = [str(p) for p in folder_path.glob(f"*{suffix}")]

    if not file_paths:
        if verbose:
            print(f"No '{suffix}' files found in {folder_path}")
        return OutputData(
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
        )

    if verbose:
        print(f"Processing {len(file_paths)} files for calculation type '{calc_type}'.")

    outputs: list[OutputData] = []

    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        results = executor.map(process_file, file_paths)
        for result in results:
            if result:
                outputs.append(result)

    return _concatenate_outputs(outputs)


def process_directory(
    directory,
    width1,
    width2,
    ewid1,
    ewid2,
    calc_types: list[Literal["GND", "EXC", "TP", "NEXAFS"]],
    verbose=True,
    n_jobs=None,
):
    """
    Process all .out files in specified subdirectories (e.g., GND, EXC, TP).

    Parameters
    ----------
    directory : str
        Path to directory containing StoBe output folders.
    width1, width2 : float
        Gaussian widths for spectral broadening (sigma values).
    ewid1, ewid2 : float
        Energy range for broadening transition.
    calc_types : list of str, optional
        A list of calculation types to process (e.g., ['GND', 'EXC', 'TP']).
        If None, defaults to ['GND', 'EXC', 'TP', 'NEXAFS'].
    verbose : bool
        Whether to print progress information.
    n_jobs : int, optional
        Number of parallel processes to use.

    Returns
    -------
    tuple : (basis_sets, energy_results, orbital_alpha, orbital_beta, xray_transitions,
    atomic_coordinates)
    """
    import numpy as np

    if calc_types is None:
        calc_types = ["GND", "EXC", "TP", "NEXAFS"]

    outputs: list[OutputData] = []

    for calc_type in calc_types:
        output = process_calculation_type(
            directory, calc_type, verbose=verbose, n_jobs=n_jobs
        )
        outputs.append(output)

    # Combine all results
    output = _concatenate_outputs(outputs)

    # Sort and clean energy results
    if not output.df_energies.empty:
        output.df_energies = sort_dataframe_naturally(output.df_energies, "Atom")
        output.df_energies["Atom"] = output.df_energies["Atom"].str.upper()
        output.df_energies = output.df_energies.drop_duplicates().reset_index(drop=True)

    # Add atom identifiers to dataframes
    orbital_alpha = output.df_alpha.copy()
    orbital_beta = output.df_beta.copy()
    xray_transitions = output.df_xray_transitions.copy()

    # Process orbital_alpha
    if not orbital_alpha.empty:
        orbital_alpha["Atom"] = orbital_alpha["Originating File"].apply(extract_atom_id)
        orbital_alpha = orbital_alpha[
            ["Atom"] + [col for col in orbital_alpha.columns if col != "Atom"]
        ]

    # Process orbital_beta
    if not orbital_beta.empty:
        orbital_beta["Atom"] = orbital_beta["Originating File"].apply(extract_atom_id)
        orbital_beta = orbital_beta[
            ["Atom"] + [col for col in orbital_beta.columns if col != "Atom"]
        ]

    # Process xray_transitions
    if not xray_transitions.empty:
        xray_transitions["Atom"] = xray_transitions["Originating File"].apply(
            extract_atom_id
        )
        xray_transitions = xray_transitions[
            ["Atom"] + [col for col in xray_transitions.columns if col != "Atom"]
        ]

    # Apply spectral broadening if X-ray transitions exist
    if not xray_transitions.empty and "E" in xray_transitions.columns:

        def broad(E):
            if ewid1 > E:
                return width1
            elif ewid2 < E:
                return width2
            else:
                return width1 + (width2 - width1) * (E - ewid1) / (ewid2 - ewid1)

        xray_transitions["width"] = xray_transitions.apply(
            lambda row: broad(row["E"]), axis=1
        )

        # Calculate normalized oscillator strength components and angles
        mag = np.sqrt(
            xray_transitions["osx"] ** 2
            + xray_transitions["osy"] ** 2
            + xray_transitions["osz"] ** 2
        )
        # Avoid division by zero
        mag = np.where(mag == 0, 1, mag)

        xray_transitions["normalized_osx"] = xray_transitions["osx"] / mag
        xray_transitions["normalized_osy"] = xray_transitions["osy"] / mag
        xray_transitions["normalized_osz"] = xray_transitions["osz"] / mag
        xray_transitions["normalized_magnitude"] = np.sqrt(
            xray_transitions["normalized_osx"] ** 2
            + xray_transitions["normalized_osy"] ** 2
            + xray_transitions["normalized_osz"] ** 2
        )

        # Calculate angles with respect to z-axis
        dot = xray_transitions["normalized_osz"]
        theta = np.arccos(np.clip(dot, -1, 1))  # Clip to avoid numerical errors
        theta_deg = np.degrees(theta)
        theta_deg = np.where(theta_deg > 90, 180 - theta_deg, theta_deg)
        xray_transitions["theta"] = theta_deg

    # Restructure energy results and apply corrections
    if not output.df_energies.empty:
        restructured_energy_results = restructure_energies_dataframe(output.df_energies)

        if verbose:
            print(
                "Available energy columns:",
                restructured_energy_results.columns.tolist(),
            )

        # Apply energy corrections to xray_transitions if both dataframes
        # exist and have data
        if not xray_transitions.empty and not restructured_energy_results.empty:
            if verbose:
                print("Applying Energy Corrections")
            xray_transitions = apply_energy_corrections(
                xray_transitions, restructured_energy_results, verbose=verbose
            )

            # Recalculate broadening with corrected energies
            if "E" in xray_transitions.columns:
                xray_transitions["width"] = xray_transitions["E"].apply(broad)
    else:
        restructured_energy_results = pd.DataFrame()

    if verbose:
        print(
            f"Processing completed. Found {len(restructured_energy_results)} "
            f"unique atoms with {len(xray_transitions)} total transitions."
        )

    return (
        output.combined,
        restructured_energy_results,
        orbital_alpha,
        orbital_beta,
        xray_transitions,
        output.df_atomic_coordinates,
    )


def restructure_energies_dataframe(df_energies):
    """
    Restructures the energies dataframe.

    From having 3 rows per atom (one for each calculation type)
    to having 1 row per atom with prefixed columns for each calculation type.

    Parameters
    ----------
    df_energies : pd.DataFrame
        Original dataframe with columns including 'Atom', 'Calculation Type',
        and various energy columns

    Returns
    -------
    pd.DataFrame : Restructured dataframe with one row per atom and prefixed columns
    """
    if df_energies.empty:
        return df_energies

    # Get the list of columns to restructure (exclude 'Atom' and 'Calculation Type')
    columns_to_restructure = [
        col for col in df_energies.columns if col not in ["Atom", "Calculation Type"]
    ]

    # Create an empty list to store the restructured data
    restructured_data = []

    # Get unique atoms
    unique_atoms = df_energies["Atom"].unique()

    for atom in unique_atoms:
        atom_data = df_energies[df_energies["Atom"] == atom]

        # Initialize the row dictionary with the atom name
        new_row = {"Atom": atom}

        # For each calculation type, add prefixed columns
        for _, row in atom_data.iterrows():
            calc_type = row["Calculation Type"]
            if calc_type:  # Only process if calculation type is not None
                for col in columns_to_restructure:
                    prefixed_col = f"{calc_type}_{col}"
                    new_row[prefixed_col] = row[col]

        restructured_data.append(new_row)

    # Create the new dataframe
    df_restructured = pd.DataFrame(restructured_data)

    # Clean up by removing redundant columns - keep only TP versions of these specific
    # columns
    columns_to_remove = []
    cleanup_columns = [
        "Ionization potential (eV)",
        "Orbital energy core hole (H)",
        "Orbital energy core hole (eV)",
        "LUMO_En",
        "Rigid spectral shift (eV)",
    ]

    for col in cleanup_columns:
        # Remove GND and EXC versions, keep only TP versions
        gnd_col = f"GND_{col}"
        exc_col = f"EXC_{col}"

        if gnd_col in df_restructured.columns:
            columns_to_remove.append(gnd_col)
        if exc_col in df_restructured.columns:
            columns_to_remove.append(exc_col)

    # Drop the redundant columns
    df_restructured = df_restructured.drop(columns=columns_to_remove, errors="ignore")
    df_restructured_final = calculate_energy_correction_restructured(df_restructured)

    return df_restructured_final


def calculate_energy_correction_restructured(df_restructured):
    """
    Calculate Energy_Correction for the restructured dataframe.

    Formula: EXC_Total energy - GND_Total energy - TP_LUMO_En
    Converts Hartree columns to eV and performs calculation.

    Parameters
    ----------
    df_restructured : pd.DataFrame
        Restructured dataframe with prefixed columns

    Returns
    -------
    pd.DataFrame : Dataframe with added 'Energy_Correction (eV)' column
    """
    df_result = df_restructured.copy()

    # Define required columns
    gnd_col = "GND_Total energy (H)"
    exc_col = "EXC_Total energy (H)"
    lumo_col = "TP_LUMO_En"

    # Check if required columns exist
    required_columns = [gnd_col, exc_col, lumo_col]
    missing_columns = [col for col in required_columns if col not in df_result.columns]

    if missing_columns:
        print(
            f"Warning: Missing columns for energy correction calculation: {missing_columns}"  # noqa: E501
        )
        df_result["Energy_Correction (eV)"] = None
        return df_result

    # Conversion factor: 1 Hartree = 27.2114 eV
    hartree_to_ev = 27.2114

    # Calculate energy correction
    def calculate_correction(row):
        exc_energy_h = row.get(exc_col)
        gnd_energy_h = row.get(gnd_col)
        lumo_en = row.get(lumo_col)

        # Check if all values are available and not None/NaN
        if (
            exc_energy_h is not None
            and gnd_energy_h is not None
            and lumo_en is not None
            and pd.notna(exc_energy_h)
            and pd.notna(gnd_energy_h)
            and pd.notna(lumo_en)
        ):
            # Convert Hartree to eV
            exc_energy_ev = exc_energy_h * hartree_to_ev
            gnd_energy_ev = gnd_energy_h * hartree_to_ev

            # Calculate energy correction
            energy_correction = exc_energy_ev - gnd_energy_ev - lumo_en

            return energy_correction
        else:
            return None

    df_result["Energy_Correction (eV)"] = df_result.apply(calculate_correction, axis=1)

    return df_result


def apply_energy_corrections(xray_transitions_df, energy_results_df, verbose=True):
    """
    Apply energy corrections from the restructured energy dataframe.

    Parameters
    ----------
    xray_transitions_df : pd.DataFrame
        DataFrame containing X-ray transitions with 'E' and 'Atom' columns
    energy_results_df : pd.DataFrame
        Restructured energy DataFrame with 'Atom' and 'Energy_Correction (eV)' columns
    verbose : bool
        Whether to print progress information

    Returns
    -------
    pd.DataFrame : Modified xray_transitions dataframe with corrected energies and
    original energies preserved
    """
    if xray_transitions_df.empty or energy_results_df.empty:
        if verbose:
            print(
                "Warning: Cannot apply energy corrections: one or both dataframes are empty"  # noqa: E501
            )
        return xray_transitions_df

    # Check if Energy_Correction column exists
    if "Energy_Correction (eV)" not in energy_results_df.columns:
        if verbose:
            print(
                "Warning: Energy_Correction (eV) column not found in energy results. Skipping corrections."  # noqa: E501
            )
        return xray_transitions_df

    # Create a copy to avoid modifying the original
    corrected_df = xray_transitions_df.copy()

    # Store original energies before correction
    corrected_df["E_original"] = corrected_df["E"].copy()

    # Create a mapping of atom to energy correction
    energy_corrections = energy_results_df.set_index("Atom")[
        "Energy_Correction (eV)"
    ].to_dict()

    # Track corrections applied
    corrections_applied = 0
    atoms_without_corrections = []

    # Apply corrections for each atom
    for atom in corrected_df["Atom"].unique():
        if atom in energy_corrections:
            correction = energy_corrections[atom]
            if pd.notna(correction):  # Only apply if correction is not NaN
                # Apply correction: corrected_energy = original_energy + correction
                mask = corrected_df["Atom"] == atom
                corrected_df.loc[mask, "E"] = corrected_df.loc[mask, "E"] + correction
                corrections_applied += 1
                if verbose:
                    num_transitions = mask.sum()
                    print(
                        f"Applied correction of {correction:.4f} eV to {num_transitions} transitions for atom {atom}"  # noqa: E501
                    )
            else:
                atoms_without_corrections.append(atom)
        else:
            atoms_without_corrections.append(atom)

    # Add a column to indicate which energies were corrected
    corrected_df["Energy_Corrected"] = corrected_df["Atom"].map(
        lambda x: x in energy_corrections and pd.notna(energy_corrections.get(x, None))
    )

    # Add the correction amount for reference
    corrected_df["Applied_Correction"] = corrected_df["Atom"].map(
        lambda x: energy_corrections.get(x, 0.0)
        if x in energy_corrections and pd.notna(energy_corrections.get(x, None))
        else 0.0
    )

    # Report summary
    if verbose:
        print("Energy corrections applied successfully!")
        print(f"- Corrections applied to {corrections_applied} unique atoms")
        print(
            f"- Total transitions corrected: {corrected_df['Energy_Corrected'].sum()}"
        )

        if atoms_without_corrections:
            unique_missing = set(atoms_without_corrections)
            print(f"- No corrections available for atoms: {', '.join(unique_missing)}")

    return corrected_df


def tdm_tensor(j, i, maxRot, normalized_osx, normalized_osy, normalized_osz):
    """
    Python implementation of the IGOR makeTDMTensor function.

    Parameters
    ----------
    j : int
        Row index parameter (though not used in the same way as IGOR version)
    i : int
        Transition index
    maxRot : int
        Maximum rotation symmetry (principal symmetry axis)
    normalized_osx, normalized_osy, normalized_osz : float
        Normalized cartesian components of transition dipole moment

    Returns
    -------
    dict : Dictionary containing the 9 tensor elements
    (xx, xy, xz, yx, yy, yz, zx, zy, zz)
    """
    # Calculate rotation degrees
    # rot_degrees = 360.0 / maxRot

    # Initialize total tensor as 3x3 matrix
    tot_tensor = np.zeros((3, 3))

    # TDM components
    tdm_components = np.array([normalized_osx, normalized_osy, normalized_osz])

    # Add up the TDM tensors for each symmetry equivalent center
    for k in range(maxRot):  # noqa: B007
        # c_rot = rot_degrees * k

        # Create individual TDM tensor (3x3 matrix)
        # Based on the IGOR code structure: outer product of TDM vector with itself
        tdm_tensor = np.outer(tdm_components, tdm_components)

        # Add to total tensor
        tot_tensor += tdm_tensor

    # Extract the 9 tensor elements and return as dictionary
    tensor_elements = {
        "xx": tot_tensor[0, 0],
        "xy": tot_tensor[0, 1],
        "xz": tot_tensor[0, 2],
        "yx": tot_tensor[1, 0],
        "yy": tot_tensor[1, 1],
        "yz": tot_tensor[1, 2],
        "zx": tot_tensor[2, 0],
        "zy": tot_tensor[2, 1],
        "zz": tot_tensor[2, 2],
    }

    return tensor_elements


def load_data(
    directory: str | Path,
    width1: float,
    width2: float,
    max_energy: float,
    maxRot=4,
    verbose=False,
    n_jobs=CPU_COUNT,
):
    """
    Load data from a directory, processes it, and returns the results as a dictionary.

    Non-Streamlit version of the data loading function.

    Parameters
    ----------
    directory : str
        Path to directory containing StoBe output folders
    width1, width2 : float
        Gaussian widths for spectral broadening (sigma values)
    max_energy : float
        Maximum energy for processing
    maxRot : int, optional
        Principal symmetry axis for TDM tensor calculation (default=4)
    verbose : bool
        Whether to print progress information
    n_jobs : int, optional
        Number of parallel processes to use

    Returns
    -------
    dict : Dictionary containing all processed dataframes
    """
    import time

    if not Path(directory).is_dir():
        print(f"Invalid directory: {directory}. Please enter a valid directory path.")
        return None

    if verbose:
        print(f"Processing directory: {directory}")
        print(f"Using maxRot = {maxRot} for TDM tensor calculation")

    start_time = time.time()

    # Process the directory
    (
        basis_sets,
        energy_results,
        orbital_alpha,
        orbital_beta,
        xray_transitions,
        atomic_coordinates,
    ) = process_directory(
        directory,
        width1,
        width2,
        290,
        max_energy,
        ["GND", "EXC", "TP", "NEXAFS"],
        verbose=verbose,
        n_jobs=n_jobs,
    )

    # Calculate TDM tensor elements for xray_transitions
    if not xray_transitions.empty and all(
        col in xray_transitions.columns
        for col in ["normalized_osx", "normalized_osy", "normalized_osz"]
    ):
        if verbose:
            print("Calculating TDM tensor elements...")

        # Initialize new columns
        tensor_columns = ["xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"]
        for col in tensor_columns:
            xray_transitions[col] = 0.0

        # Calculate tensor elements for each transition
        for idx, row in xray_transitions.iterrows():
            tensor_elements = tdm_tensor(
                j=0,  # Not used in the same way as IGOR version
                i=idx,
                maxRot=maxRot,
                normalized_osx=row["normalized_osx"],
                normalized_osy=row["normalized_osy"],
                normalized_osz=row["normalized_osz"],
            )

            # Assign tensor elements to dataframe
            for col in tensor_columns:
                xray_transitions.at[idx, col] = tensor_elements[col]

        if verbose:
            print(
                f"TDM tensor calculation completed for {len(xray_transitions)} transitions"  # noqa: E501
            )

    elif not xray_transitions.empty:
        print(
            "Warning: Required columns (normalized_osx, normalized_osy, normalized_osz) not found in xray_transitions"  # noqa: E501
        )

    end_time = time.time()

    if verbose:
        print(f"Processing completed in {end_time - start_time:.2f} seconds.")

    # Create the data dictionary
    data = {
        "basis_sets": basis_sets,
        "energy_results": energy_results,
        "orbital_alpha": orbital_alpha,
        "orbital_beta": orbital_beta,
        "xray_transitions": xray_transitions,
        "atomic_coordinates": atomic_coordinates,
    }

    # Print summary if verbose
    if verbose:
        print("\n" + "=" * 50)
        print("PROCESSING COMPLETE - DATA OVERVIEW")
        print("=" * 50)

        # Summary metrics
        unique_atoms = (
            energy_results["Atom"].nunique() if not energy_results.empty else 0
        )
        total_transitions = len(xray_transitions) if not xray_transitions.empty else 0
        alpha_orbitals = len(orbital_alpha) if not orbital_alpha.empty else 0
        atomic_coords = len(atomic_coordinates) if not atomic_coordinates.empty else 0

        print(f"Unique Atoms: {unique_atoms}")
        print(f"X-ray Transitions: {total_transitions}")
        print(f"Alpha Orbitals: {alpha_orbitals}")
        print(f"Atomic Coordinates: {atomic_coords}")

        # Show energy corrections info if available
        if (
            not xray_transitions.empty
            and "Energy_Corrected" in xray_transitions.columns
        ):
            corrected_count = xray_transitions["Energy_Corrected"].sum()
            total_count = len(xray_transitions)
            print(
                f"Energy corrections applied to {corrected_count}/{total_count} transitions"  # noqa: E501
            )

        print("=" * 50)

    return data
