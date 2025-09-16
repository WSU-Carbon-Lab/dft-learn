#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "typer",
#     "rich",
#     "natsort",
# ]
# ///

"""
StoBe Input Generator CLI

Generate StoBe DFT calculation input files from molecular configurations.

Usage:
    uv run generate.py [run_directory] [xyz_file]
    uv run generate.py tcta_flat                              # Auto-detect XYZ file
    uv run generate.py tcta_flat tcta_flat/tcta_flat.xyz       # Explicit XYZ file
    uv run generate.py tcta_flat --xyz another_file.xyz       # Use --xyz flag
"""

import glob
import importlib.util
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Optional

import typer
from rich import print
from rich.console import Console
from rich.panel import Panel
from rich.progress import track

console = Console()

# Import the functions from the original generator
# We'll define them here directly to avoid import complexity
from natsort import natsort_keygen, ns

# StoBe defaults
DEFAULT_SYM = "C1"
DEFAULT_MULT = "1"
DEFAULT_GEOUNITS = "ANGSTROMS"
DEFAULT_RUNTYPE = "startup nooptimize"
DEFAULT_SCFTYPE = "direct"
DEFAULT_POTENTIAL = "nonlocal rpbe pbe"
DEFAULT_GRID = "fine"
DEFAULT_CHARGE = "0"
DEFAULT_MAXCYCLES = "300"
DEFAULT_ECONVERGENCE = "0.001"
DEFAULT_DCONVERGENCE = "0.001"
DEFAULT_DMIXING = "mdens 0.05"
DEFAULT_DIIS = "new 7"
DEFAULT_ORBI = "5d"
DEFAULT_MULLIKEN = "on full"
DEFAULT_VIRT = "all"
DEFAULT_FSYMGND = "scfocc"
DEFAULT_MOFILE = "molden"
DEFAULT_FSYMEXC = "scfocc excited"
DEFAULT_ALPHAOCC = "0 1 1 0.0"
DEFAULT_BETAOCC = "0 0"
DEFAULT_ALPHAOCCTP = "0 1 1 0.5"
DEFAULT_SPIN = "FULL"
DEFAULT_XRAY = "XAS"

# StoBe paths
STOBE_HOME: Path = Path(os.environ.get("STOBE_HOME", "/bin/stobe"))
stobe_exe_path = shutil.which("StoBe.x")
if stobe_exe_path:
    STOBE: Path = Path(stobe_exe_path)
else:
    STOBE: Path = STOBE_HOME / "Source/StoBe.x"

xasrun_exe_path = shutil.which("xrayspec.x")
if xasrun_exe_path:
    XASRUN: Path = Path(xasrun_exe_path)
else:
    XASRUN: Path = STOBE_HOME / "Source/xrayspec.x"
SYMBASIS: Path = STOBE_HOME / "Basis/symbasis.new"
BASISLIB: Path = STOBE_HOME / "Basis/baslib.new7"

app = typer.Typer(
    help="StoBe DFT Input Generator - Generate calculation input files from molecular configurations",
    rich_markup_mode="rich"
)


def load_mol_config(config_dir: Path):
    """Load molConfig.py from the specified directory"""
    config_path = config_dir / "molConfig.py"

    if not config_path.exists():
        console.print(f"[red]Error: molConfig.py not found in {config_dir}[/red]")
        raise typer.Exit(1)

    # Load the module dynamically
    spec = importlib.util.spec_from_file_location("molConfig", config_path)
    if spec is None or spec.loader is None:
        console.print(f"[red]Error: Could not load molConfig.py from {config_path}[/red]")
        raise typer.Exit(1)

    mol_config = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mol_config)

    return mol_config


def find_xyz_files(directory: Path) -> list[Path]:
    """Find all .xyz files in the directory"""
    xyz_files = list(directory.glob("*.xyz"))
    return xyz_files


def auto_detect_xyz(run_directory: Path) -> Path:
    """Auto-detect the XYZ file in the run directory"""
    xyz_files = find_xyz_files(run_directory)

    if not xyz_files:
        console.print(f"[red]Error: No .xyz files found in {run_directory}[/red]")
        raise typer.Exit(1)

    if len(xyz_files) == 1:
        console.print(f"[green]Auto-detected XYZ file: {xyz_files[0]}[/green]")
        return xyz_files[0]

    # Multiple files found - try to find one matching the directory name
    dir_name = run_directory.name
    for xyz_file in xyz_files:
        if xyz_file.stem == dir_name:
            console.print(f"[green]Auto-detected XYZ file: {xyz_file} (matches directory name)[/green]")
            return xyz_file

    # If no exact match, prompt user
    console.print(f"[yellow]Multiple .xyz files found in {run_directory}:[/yellow]")
    for i, xyz_file in enumerate(xyz_files, 1):
        console.print(f"  {i}. {xyz_file.name}")

    choice = typer.prompt(f"Select XYZ file (1-{len(xyz_files)})", type=int)
    if 1 <= choice <= len(xyz_files):
        return xyz_files[choice - 1]
    else:
        console.print("[red]Invalid selection[/red]")
        raise typer.Exit(1)


def makeRunFile(
    fname,
    mname,
    aname,
    nAtoms,
    alphaGND,
    betaGND,
    element1a,
    element1b,
    element2,
    element3,
    element4,
    nElem1,
    nElem2,
    nElem3,
    nElem4,
    nElem5,
    nElem6,
    difElements,
    elem1_Abasis_a,
    elem1_Abasis_b,
    elem2_Abasis,
    elem3_Abasis,
    elem4_Abasis,
    elem5_Abasis,
    elem6_Abasis,
    elem1_Obasis_a,
    elem1_Obasis_b,
    elem2_Obasis,
    elem3_Obasis,
    elem4_Obasis,
    elem5_Obasis,
    elem6_Obasis,
    elem1_MCPbasis,
    title,
    sym,
    geom,
    runtype,
    scftype,
    potential,
    grid,
    multiplicity,
    charge,
    maxcycles,
    econvergence,
    dconvergence,
    dmixing,
    diis,
    orbi,
    mulliken,
    virt,
    spin,
    mofile,
    fsymGND,
    fsymEXC,
    fsymTP,
    alfaOcc,
    alfaOccTP,
    betaOcc,
    ftype,
):
    """Generate run files with proper geometry parsing"""
    # Apply defaults
    sym = sym or DEFAULT_SYM
    geom = geom or DEFAULT_GEOUNITS
    runtype = runtype or DEFAULT_RUNTYPE
    scftype = scftype or DEFAULT_SCFTYPE
    potential = potential or DEFAULT_POTENTIAL
    grid = grid or DEFAULT_GRID
    multiplicity = multiplicity or DEFAULT_MULT
    charge = charge or DEFAULT_CHARGE
    maxcycles = maxcycles or DEFAULT_MAXCYCLES
    econvergence = econvergence or DEFAULT_ECONVERGENCE
    dconvergence = dconvergence or DEFAULT_DCONVERGENCE
    dmixing = dmixing or DEFAULT_DMIXING
    diis = diis or DEFAULT_DIIS
    orbi = orbi or DEFAULT_ORBI
    mulliken = mulliken or DEFAULT_MULLIKEN
    virt = virt or DEFAULT_VIRT
    fsymGND = fsymGND or DEFAULT_FSYMGND
    mofile = mofile or DEFAULT_MOFILE
    alfaOcc = alfaOcc or DEFAULT_ALPHAOCC
    betaOcc = betaOcc or DEFAULT_BETAOCC
    fsymEXC = fsymEXC or DEFAULT_FSYMEXC
    fsymTP = fsymTP or DEFAULT_FSYMEXC
    spin = spin or DEFAULT_SPIN

    # Create run files for each atom
    for i in range(1, nAtoms + 1):
        run_file = f"{aname}{i}{ftype}.run"

        with open(run_file, "w+", newline="\n") as f:
            f.write("#!/bin/csh -f\n")
            f.write(f"ln -s {BASISLIB} fort.3\n")
            f.write(f"ln -s {SYMBASIS} fort.4\n")
            f.write(f"cat >{aname}{i}{ftype}.inp<</.\n")
            f.write("TITLE\n")
            f.write(f"{title} {ftype.upper()}\n")
            f.write(f"SYMMETRY {sym}\n")
            f.write(f"CARTESIAN {geom}\n")

        # Process geometry from XYZ file with proper element handling
        with open(f"{fname}.xyz") as xyz_file:
            xyz_lines = xyz_file.readlines()

        with open(run_file, "a+", newline="\n") as f:
            n1 = n2 = n3 = n4 = n5 = n6 = 1

            for line in xyz_lines:
                line = line.strip()
                if not line or len(line.split()) < 4:
                    continue

                # Parse atom info from XYZ line
                parts = line.split()
                if len(parts) >= 4:
                    # Handle different numbers of elements
                    if difElements == 3:
                        if n1 <= nElem1:
                            # Element 1 (target element for calculations)
                            if n1 == i:
                                # This is the excited atom
                                new_line = f"{line}     {element1a}     32\n"
                            else:
                                # Regular atom of same element
                                new_line = f"{line}     {element1b}     32\n"
                            f.write(new_line)
                            n1 += 1
                        elif n2 <= nElem2:
                            # Element 2
                            new_line = f"{line}     {element2}     32\n"
                            f.write(new_line)
                            n2 += 1
                        elif n3 <= nElem3:
                            # Element 3
                            new_line = f"{line}     {element3}     32\n"
                            f.write(new_line)
                            n3 += 1

                    elif difElements == 2:
                        if n1 <= nElem1:
                            if n1 == i:
                                new_line = f"{line}     {element1a}     32\n"
                            else:
                                new_line = f"{line}     {element1b}     32\n"
                            f.write(new_line)
                            n1 += 1
                        elif n2 <= nElem2:
                            new_line = f"{line}     {element2}     32\n"
                            f.write(new_line)
                            n2 += 1

                    elif difElements == 1:
                        if n1 <= nElem1:
                            if n1 == i:
                                new_line = f"{line}     {element1a}     32\n"
                            else:
                                new_line = f"{line}     {element1b}     32\n"
                            f.write(new_line)
                            n1 += 1

                    # Add support for more elements if needed
                    elif difElements >= 4:
                        if n1 <= nElem1:
                            if n1 == i:
                                new_line = f"{line}     {element1a}     32\n"
                            else:
                                new_line = f"{line}     {element1b}     32\n"
                            f.write(new_line)
                            n1 += 1
                        elif n2 <= nElem2:
                            new_line = f"{line}     {element2}     32\n"
                            f.write(new_line)
                            n2 += 1
                        elif n3 <= nElem3:
                            new_line = f"{line}     {element3}     32\n"
                            f.write(new_line)
                            n3 += 1
                        elif n4 <= nElem4:
                            new_line = f"{line}     {element4}     32\n"
                            f.write(new_line)
                            n4 += 1

        # Continue with run file
        with open(run_file, "a", newline="\n") as f:
            f.write("\nEND\n")

            # Write calculation parameters
            f.write(f"RUNTYPE {runtype}\n")
            f.write(f"SCFTYPE {scftype}\n")
            f.write(f"POTENTIAL {potential}\n")
            f.write(f"GRID {grid}\n")
            f.write(f"MULTIPLICITY {multiplicity}\n")
            f.write(f"CHARGE {charge}\n")
            f.write(f"MAXCYCLES {maxcycles}\n")
            f.write(f"ECONVERGENCE {econvergence}\n")
            f.write(f"DCONVERGENCE {dconvergence}\n")
            f.write(f"DMIXING {dmixing}\n")
            f.write(f"DIIS {diis}\n")
            f.write(f"ORBI {orbi}\n")
            f.write(f"MULLIKEN {mulliken}\n")
            f.write(f"VIRT {virt}\n")

            # File type specific sections
            if ftype == "gnd":
                f.write(f"FSYM {fsymGND}\n")
                f.write(f"ALFA {alphaGND}\n")
                f.write(f"BETA {betaGND}\n")
                f.write(f"FILE {mofile}\n")
                f.write("END\n")
            elif ftype == "exc":
                f.write(f"FSYM {fsymEXC}\n")
                f.write(f"ALFA {alphaGND + 1}\n")
                f.write(f"BETA {betaGND}\n")
                f.write("SYM 1\n")
                f.write(f"ALFA {alfaOcc}\n")
                f.write(f"BETA {betaOcc}\n")
                f.write("END\n")
                f.write(f"MULLIKEN {mulliken}\n")
                f.write(f"FILE {mofile}\n")
                f.write("END\n")
            elif ftype == "tp":
                f.write(f"FSYM {fsymTP}\n")
                f.write(f"ALFA {alphaGND}\n")
                f.write(f"BETA {betaGND}\n")
                f.write("SYM 1\n")
                f.write(f"ALFA {alfaOccTP}\n")
                f.write(f"BETA {betaOcc}\n")
                f.write("END\n")
                f.write("MULLIKEN on full\n")
                f.write(f"FILE {mofile}\n")
                f.write("XRAY xas\n")
                f.write("END\n")
                f.write("END\n")

            # Add basis set sections (simplified - could be expanded)
            # Auxiliary basis sets
            f.write("END\n")

            # Finalize run file
            f.write("/.\n")
            f.write(f"StoBe.x <{aname}{i}{ftype}.inp>& {aname}{i}{ftype}.out\n")
            f.write(f"mv Molden.molf {aname}{i}{ftype}.molden\n")
            if ftype == "tp":
                f.write(f"mv fort.11 {aname}{i}.xas\n")
            f.write("rm fort.*\n")


def makeXASrun(mname, aname, nAtoms, title):
    """Generate XAS run files"""
    # Element-specific energy ranges
    energy_ranges = {
        "C": ("280 320", "0.6 12 288 320"),
        "N": ("400 450", "0.6 12 418 450"),
        "O": ("530 580", "0.6 12 535 560"),
        "Ni": ("830 880", "0.6 12 830 880"),
    }

    element = aname.upper()
    if element not in energy_ranges:
        console.print(f"[yellow]Warning: Energy range not defined for element {element}[/yellow]")
        return

    erange, ebroad = energy_ranges[element]

    for i in range(1, nAtoms + 1):
        run_file = f"{aname}{i}xas.run"
        with open(run_file, "w+", newline="\n") as f:
            f.write("#!/bin/csh -f\n")
            f.write(f"ln -s ~/{mname}/{aname}{i}/{aname}{i}.xas fort.1\n")
            f.write(f"cat >{aname}{i}xas.inp<</.\n")
            f.write("title\n")
            f.write(f"{title} XAS\n")
            f.write("PRINT\n")
            f.write(f"RANGE {erange}\n")
            f.write("POINTS 2000\n")
            f.write(f"WIDTH {ebroad}\n")
            f.write("XRAY xas\n")
            f.write("TOTAL 1\n")
            f.write("END\n")
            f.write("/.\n")
            f.write(f"{XASRUN} <{aname}{i}xas.inp>& {aname}{i}xas.out\n")
            f.write(f"cp XrayT001.out {aname}{i}.out\n")
            f.write("rm fort.*\n")


def makeSEQrun(nAtoms, aname):
    """Generate sequential run files"""
    for i in range(1, nAtoms + 1):
        with open(f"{aname}{i}seq.run", "w+", newline="\n") as f:
            f.write(f"chmod +x ./{aname}{i}gnd.run\n")
            f.write(f"./{aname}{i}gnd.run\n\n")
            f.write(f"chmod +x ./{aname}{i}exc.run\n")
            f.write(f"./{aname}{i}exc.run\n\n")
            f.write(f"chmod +x ./{aname}{i}tp.run\n")
            f.write(f"./{aname}{i}tp.run\n\n")
            f.write(f"chmod +x ./{aname}{i}xas.run\n")
            f.write(f"./{aname}{i}xas.run\n")


def makeFolders(aName, nAtoms, mname):
    """Create directories for calculations"""
    for i in range(1, nAtoms + 1):
        dir_name = f"{aName}{i}"
        os.makedirs(dir_name, exist_ok=True)

    os.makedirs(mname, exist_ok=True)


def listFiles(aname, nAtoms):
    """Organize run files into folders"""
    run_files = glob.glob("*.run")
    natSortKey = natsort_keygen(key=lambda y: y.lower(), alg=ns.IGNORECASE)
    run_files.sort(key=natSortKey)

    for i in range(1, nAtoms + 1):
        atom_dir = f"{aname}{i}"
        # Create regex pattern that matches exactly the atom number followed by a non-digit
        # Pattern: aname + number + (non-digit or end of string)
        pattern = rf"^{re.escape(aname)}{i}(?:\D|$)"

        for run_file in run_files:
            if re.match(pattern, run_file):
                dest_path = os.path.join(atom_dir, run_file)
                if os.path.exists(run_file):
                    shutil.move(run_file, dest_path)


@app.command()
def main(
    run_directory: str = typer.Argument(..., help="Directory containing molConfig.py"),
    xyz_file: Optional[str] = typer.Argument(None, help="Path to XYZ geometry file"),
    xyz: Optional[str] = typer.Option(None, "--xyz", help="Alternative way to specify XYZ file"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose output"),
):
    """
    Generate StoBe DFT calculation input files.

    Examples:
        uv run generate.py tcta_flat                              # Auto-detect XYZ
        uv run generate.py tcta_flat tcta_flat.xyz                # Explicit XYZ
        uv run generate.py tcta_flat --xyz another_file.xyz       # Using --xyz flag
    """

    console.print(Panel(
        "[bold blue]StoBe Input Generator[/bold blue]\n"
        "Generating DFT calculation input files...",
        title="🧪 StoBe DFT",
        border_style="blue"
    ))

    # Convert to Path objects and get absolute path early
    run_dir = Path(run_directory).resolve()

    if not run_dir.exists():
        console.print(f"[red]Error: Directory {run_directory} does not exist[/red]")
        raise typer.Exit(1)    # Determine XYZ file - convert to absolute path immediately
    xyz_path = None
    if xyz:
        xyz_path = Path(xyz).resolve()
    elif xyz_file:
        xyz_path = Path(xyz_file).resolve()
    else:
        xyz_path = auto_detect_xyz(run_dir).resolve()  # Make absolute

    if not xyz_path.exists():
        console.print(f"[red]Error: XYZ file {xyz_path} does not exist[/red]")
        raise typer.Exit(1)

    # Load configuration
    console.print(f"[blue]Loading configuration from {run_dir}/molConfig.py[/blue]")
    mol_config = load_mol_config(run_dir)

    if verbose:
        console.print(f"[dim]XYZ file: {xyz_path}[/dim]")
        console.print(f"[dim]Run directory: {run_dir}[/dim]")

    # Extract configuration values
    mname = run_directory  # Use directory name as molecule name
    fname = xyz_path.stem  # Use XYZ file stem as base name

    # Get values from mol_config with error handling
    try:
        aname = mol_config.aname
        nFiles = mol_config.nFiles
        difElems = mol_config.difElems

        # Element data
        alpha = mol_config.alpha
        beta = mol_config.beta
        element1a = mol_config.element1a
        element1b = mol_config.element1b
        element2 = getattr(mol_config, 'element2', 0)
        element3 = getattr(mol_config, 'element3', 0)
        element4 = getattr(mol_config, 'element4', 0)

        # Element counts
        nElem1 = mol_config.nElem1
        nElem2 = getattr(mol_config, 'nElem2', 0)
        nElem3 = getattr(mol_config, 'nElem3', 0)
        nElem4 = getattr(mol_config, 'nElem4', 0)
        nElem5 = getattr(mol_config, 'nElem5', 0)
        nElem6 = getattr(mol_config, 'nElem6', 0)

        # Basis sets
        elem1_Abasis_a = mol_config.elem1_Abasis_a
        elem1_Abasis_b = mol_config.elem1_Abasis_b
        elem2_Abasis = getattr(mol_config, 'elem2_Abasis', "")
        elem3_Abasis = getattr(mol_config, 'elem3_Abasis', "")
        elem4_Abasis = getattr(mol_config, 'elem4_Abasis', "")
        elem5_Abasis = getattr(mol_config, 'elem5_Abasis', "")
        elem6_Abasis = getattr(mol_config, 'elem6_Abasis', "")

        elem1_Obasis_a = mol_config.elem1_Obasis_a
        elem1_Obasis_b = mol_config.elem1_Obasis_b
        elem2_Obasis = getattr(mol_config, 'elem2_Obasis', "")
        elem3_Obasis = getattr(mol_config, 'elem3_Obasis', "")
        elem4_Obasis = getattr(mol_config, 'elem4_Obasis', "")
        elem5_Obasis = getattr(mol_config, 'elem5_Obasis', "")
        elem6_Obasis = getattr(mol_config, 'elem6_Obasis', "")

        elem1_MCPbasis = getattr(mol_config, 'elem1_MCPbasis', "")

        # Electron configuration
        alfaOcc = mol_config.alfaOcc
        betaOcc = mol_config.betaOcc
        alfaOccTP = mol_config.alfaOccTP

        # Other parameters
        title = getattr(mol_config, 'title', mname)
        multiplicity = getattr(mol_config, 'multiplicity', None)

    except AttributeError as e:
        console.print(f"[red]Error: Missing required parameter in molConfig.py: {e}[/red]")
        raise typer.Exit(1)

    # Change to run directory for file generation
    original_cwd = os.getcwd()
    os.chdir(run_dir)

    try:
        # Check if XYZ file is already in working directory with correct name
        local_xyz = Path(f"{fname}.xyz")

        if not local_xyz.exists():
            # Copy XYZ file to working directory
            shutil.copy2(xyz_path, local_xyz)
            if verbose:
                console.print(f"[dim]Copied {xyz_path} to {local_xyz}[/dim]")
        elif verbose:
            console.print(f"[dim]Using existing {local_xyz} file[/dim]")

        # Generate files with progress tracking
        file_types = ["gnd", "exc", "tp"]

        for ftype in file_types:
            makeRunFile(
                fname, mname, aname, nFiles, alpha, beta,
                element1a, element1b, element2, element3, element4,
                nElem1, nElem2, nElem3, nElem4, nElem5, nElem6, difElems,
                elem1_Abasis_a, elem1_Abasis_b, elem2_Abasis, elem3_Abasis,
                elem4_Abasis, elem5_Abasis, elem6_Abasis,
                elem1_Obasis_a, elem1_Obasis_b, elem2_Obasis, elem3_Obasis,
                elem4_Obasis, elem5_Obasis, elem6_Obasis, elem1_MCPbasis,
                title, getattr(mol_config, 'sym', None),
                getattr(mol_config, 'geom', None),
                getattr(mol_config, 'runtype', None),
                getattr(mol_config, 'scftype', None),
                getattr(mol_config, 'potential', None),
                getattr(mol_config, 'grid', None),
                multiplicity,
                getattr(mol_config, 'charge', None),
                getattr(mol_config, 'maxcycles', None),
                getattr(mol_config, 'econvergence', None),
                getattr(mol_config, 'dconvergence', None),
                getattr(mol_config, 'dmixing', None),
                getattr(mol_config, 'diis', None),
                getattr(mol_config, 'orbi', None),
                getattr(mol_config, 'mulliken', None),
                getattr(mol_config, 'virt', None),
                getattr(mol_config, 'spin', None),
                getattr(mol_config, 'mofile', None),
                getattr(mol_config, 'fsymGND', None),
                getattr(mol_config, 'fsymEXC', None),
                getattr(mol_config, 'fsymEXC', None),  # fsymTP
                alfaOcc, alfaOccTP, betaOcc, ftype
            )

        makeXASrun(mname, aname, nFiles, title)
        makeSEQrun(nFiles, aname)

        # Create atom directories directly in current directory (run_dir)
        for i in range(1, nFiles + 1):
            dir_name = f"{aname}{i}"
            os.makedirs(dir_name, exist_ok=True)

        listFiles(aname, nFiles)

        console.print(f"[bold green]✅ Successfully generated StoBe input files for {nFiles} {aname} atoms[/bold green]")
        console.print(f"[cyan]Output directory: {run_dir}[/cyan]")

    finally:
        os.chdir(original_cwd)


if __name__ == "__main__":
    app()
