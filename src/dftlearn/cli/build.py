"""dftrun build: generate StoBe input files from molConfig + XYZ."""

# ruff: noqa: B023

from __future__ import annotations

import os
import re
import shutil
from pathlib import Path

import typer
from natsort import natsort_keygen, ns
from rich.console import Console
from rich.panel import Panel

from dftlearn.cli.shared import load_mol_config

_CONSOLE = Console()

_DEFAULT_SYM = "C1"
_DEFAULT_MULT = "1"
_DEFAULT_GEOUNITS = "ANGSTROMS"
_DEFAULT_RUNTYPE = "startup nooptimize"
_DEFAULT_SCFTYPE = "direct"
_DEFAULT_POTENTIAL = "nonlocal rpbe pbe"
_DEFAULT_GRID = "fine"
_DEFAULT_CHARGE = "0"
_DEFAULT_MAXCYCLES = "300"
_DEFAULT_ECONVERGENCE = "0.001"
_DEFAULT_DCONVERGENCE = "0.001"
_DEFAULT_DMIXING = "mdens 0.05"
_DEFAULT_DIIS = "new 7"
_DEFAULT_ORBI = "5d"
_DEFAULT_MULLIKEN = "on full"
_DEFAULT_VIRT = "all"
_DEFAULT_FSYMGND = "scfocc"
_DEFAULT_MOFILE = "molden"
_DEFAULT_FSYMEXC = "scfocc excited"
_DEFAULT_ALPHAOCC = "0 1 1 0.0"
_DEFAULT_BETAOCC = "0 0"
_DEFAULT_ALPHAOCCTP = "0 1 1 0.5"
_DEFAULT_SPIN = "FULL"

_stobe_home = Path(os.environ.get("STOBE_HOME", "/bin/stobe"))
_stobe_exe = shutil.which("StoBe.x")
STOBE = Path(_stobe_exe) if _stobe_exe else _stobe_home / "Source" / "StoBe.x"
_xas_exe = shutil.which("xrayspec.x")
XASRUN = Path(_xas_exe) if _xas_exe else _stobe_home / "Source" / "xrayspec.x"
SYMBASIS = _stobe_home / "Basis" / "symbasis.new"
BASISLIB = _stobe_home / "Basis" / "baslib.new7"


def _find_xyz_files(directory: Path) -> list[Path]:
    return list(directory.glob("*.xyz"))


def _auto_detect_xyz(run_directory: Path) -> Path:
    xyz_files = _find_xyz_files(run_directory)
    if not xyz_files:
        raise typer.Exit(1)
    if len(xyz_files) == 1:
        return xyz_files[0]
    dir_name = run_directory.name
    for xyz_file in xyz_files:
        if xyz_file.stem == dir_name:
            return xyz_file
    choice = typer.prompt(f"Select XYZ file (1-{len(xyz_files)})", type=int)
    if 1 <= choice <= len(xyz_files):
        return xyz_files[choice - 1]
    raise typer.Exit(1)


def _basis_line(
    group: int,
    is_excited: bool,
    e1a: str,
    e1b: str,
    e2: str,
    e3: str,
    e4: str,
    e5: str,
    e6: str,
) -> str:
    if group == 1:
        return e1a if is_excited else e1b
    if group == 2:
        return e2
    if group == 3:
        return e3
    if group == 4:
        return e4
    if group == 5:
        return e5
    if group == 6:
        return e6
    return ""


def _orbital_line(
    group: int,
    is_excited: bool,
    e1a: str,
    e1b: str,
    e2: str,
    e3: str,
    e4: str,
    e5: str,
    e6: str,
) -> str:
    if group == 1:
        return e1a if is_excited else e1b
    if group == 2:
        return e2
    if group == 3:
        return e3
    if group == 4:
        return e4
    if group == 5:
        return e5
    if group == 6:
        return e6
    return ""


def _make_run_file(
    fname: str,
    mname: str,
    aname: str,
    n_atoms: int,
    alpha_gnd: int,
    beta_gnd: int,
    element1a: int,
    element1b: int,
    element2: int,
    element3: int,
    element4: int,
    element5: int,
    element6: int,
    n_elem1: int,
    n_elem2: int,
    n_elem3: int,
    n_elem4: int,
    n_elem5: int,
    n_elem6: int,
    dif_elements: int,
    ab1a: str,
    ab1b: str,
    ab2: str,
    ab3: str,
    ab4: str,
    ab5: str,
    ab6: str,
    ob1a: str,
    ob1b: str,
    ob2: str,
    ob3: str,
    ob4: str,
    ob5: str,
    ob6: str,
    mcp1: str,
    title: str,
    sym: str | None,
    geom: str | None,
    runtype: str | None,
    scftype: str | None,
    potential: str | None,
    grid: str | None,
    multiplicity: str | None,
    charge: str | None,
    maxcycles: str | None,
    econvergence: str | None,
    dconvergence: str | None,
    dmixing: str | None,
    diis: str | None,
    orbi: str | None,
    mulliken: str | None,
    virt: str | None,
    spin: str | None,
    mofile: str | None,
    fsym_gnd: str | None,
    fsym_exc: str | None,
    fsym_tp: str | None,
    alfa_occ: str,
    alfa_occ_tp: str,
    beta_occ: str,
    ftype: str,
) -> None:
    sym = sym or _DEFAULT_SYM
    geom = geom or _DEFAULT_GEOUNITS
    runtype = runtype or _DEFAULT_RUNTYPE
    scftype = scftype or _DEFAULT_SCFTYPE
    potential = potential or _DEFAULT_POTENTIAL
    grid = grid or _DEFAULT_GRID
    multiplicity = multiplicity or _DEFAULT_MULT
    charge = charge or _DEFAULT_CHARGE
    maxcycles = maxcycles or _DEFAULT_MAXCYCLES
    econvergence = econvergence or _DEFAULT_ECONVERGENCE
    dconvergence = dconvergence or _DEFAULT_DCONVERGENCE
    dmixing = dmixing or _DEFAULT_DMIXING
    diis = diis or _DEFAULT_DIIS
    orbi = orbi or _DEFAULT_ORBI
    mulliken = mulliken or _DEFAULT_MULLIKEN
    virt = virt or _DEFAULT_VIRT
    fsym_gnd = fsym_gnd or _DEFAULT_FSYMGND
    mofile = mofile or _DEFAULT_MOFILE
    alfa_occ = alfa_occ or _DEFAULT_ALPHAOCC
    beta_occ = beta_occ or _DEFAULT_BETAOCC
    fsym_exc = fsym_exc or _DEFAULT_FSYMEXC
    fsym_tp = fsym_tp or _DEFAULT_FSYMEXC
    spin = spin or _DEFAULT_SPIN

    for i in range(1, n_atoms + 1):
        run_file = f"{aname}{i}{ftype}.run"
        with Path(run_file).open("w+", newline="\n") as f:
            f.write("#!/bin/bash\n")
            f.write(f"ln -sf {BASISLIB} fort.3\n")
            f.write(f"ln -sf {SYMBASIS} fort.4\n")
            f.write(f"cat >{aname}{i}{ftype}.inp<</.\n")
            f.write("TITLE\n")
            f.write(f"{title} {ftype.upper()}\n")
            f.write(f"SYMMETRY {sym}\n")
            f.write(f"CARTESIAN {geom}\n")

        with Path(f"{fname}.xyz").open() as xyz_file:
            xyz_lines = xyz_file.readlines()

        atoms_for_basis: list[tuple[int, bool]] = []
        n1 = n2 = n3 = n4 = n5 = n6 = 1

        with Path(run_file).open("a+", newline="\n") as f:

            def wr(
                g: int,
                ex: bool,
                el_a: int | str,
                el_b: int | str,
                el_oth: int | str,
                nl: bool = True,
            ) -> None:
                ch = el_a if (g == 1 and ex) else (el_b if g == 1 else el_oth)
                x = round(float(parts[1]), 6)
                y = round(float(parts[2]), 6)
                z = round(float(parts[3]), 6)
                row = f"{label}     {x:10.6f}    {y:10.6f}    {z:10.6f}     {ch}     32"
                f.write(row + ("\n" if nl else ""))
                atoms_for_basis.append((g, ex))

            def last_line() -> bool:
                if dif_elements == 1:
                    return n1 == n_elem1
                if dif_elements == 2:
                    return n2 == n_elem2
                if dif_elements == 3:
                    return n3 == n_elem3
                if dif_elements == 4:
                    return n4 == n_elem4
                if dif_elements == 5:
                    return n5 == n_elem5
                return dif_elements >= 6 and n6 == n_elem6

            def nl() -> bool:
                return not last_line()

            for line in xyz_lines:
                line = line.strip()
                if not line or len(line.split()) < 4:
                    continue
                parts = line.split()
                if len(parts) < 4:
                    continue
                label = parts[0]

                if dif_elements == 1:
                    if n1 <= n_elem1:
                        wr(1, n1 == i, element1a, element1b, element1b, nl())
                        n1 += 1
                elif dif_elements == 2:
                    if n1 <= n_elem1:
                        wr(1, n1 == i, element1a, element1b, element1b, True)
                        n1 += 1
                    elif n2 <= n_elem2:
                        wr(2, False, element1a, element1b, element2, nl())
                        n2 += 1
                elif dif_elements == 3:
                    if n1 <= n_elem1:
                        wr(1, n1 == i, element1a, element1b, element1b, True)
                        n1 += 1
                    elif n2 <= n_elem2:
                        wr(2, False, element1a, element1b, element2, True)
                        n2 += 1
                    elif n3 <= n_elem3:
                        wr(3, False, element1a, element1b, element3, nl())
                        n3 += 1
                elif dif_elements == 4:
                    if n1 <= n_elem1:
                        wr(1, n1 == i, element1a, element1b, element1b, True)
                        n1 += 1
                    elif n2 <= n_elem2:
                        wr(2, False, element1a, element1b, element2, True)
                        n2 += 1
                    elif n3 <= n_elem3:
                        wr(3, False, element1a, element1b, element3, True)
                        n3 += 1
                    elif n4 <= n_elem4:
                        wr(4, False, element1a, element1b, element4, nl())
                        n4 += 1
                elif dif_elements == 5:
                    if n1 <= n_elem1:
                        wr(1, n1 == i, element1a, element1b, element1b, True)
                        n1 += 1
                    elif n2 <= n_elem2:
                        wr(2, False, element1a, element1b, element2, True)
                        n2 += 1
                    elif n3 <= n_elem3:
                        wr(3, False, element1a, element1b, element3, True)
                        n3 += 1
                    elif n4 <= n_elem4:
                        wr(4, False, element1a, element1b, element4, True)
                        n4 += 1
                    elif n5 <= n_elem5:
                        wr(5, False, element1a, element1b, element5, nl())
                        n5 += 1
                elif dif_elements >= 6:
                    if n1 <= n_elem1:
                        wr(1, n1 == i, element1a, element1b, element1b, True)
                        n1 += 1
                    elif n2 <= n_elem2:
                        wr(2, False, element1a, element1b, element2, True)
                        n2 += 1
                    elif n3 <= n_elem3:
                        wr(3, False, element1a, element1b, element3, True)
                        n3 += 1
                    elif n4 <= n_elem4:
                        wr(4, False, element1a, element1b, element4, True)
                        n4 += 1
                    elif n5 <= n_elem5:
                        wr(5, False, element1a, element1b, element5, True)
                        n5 += 1
                    elif n6 <= n_elem6:
                        wr(6, False, element1a, element1b, element6, nl())
                        n6 += 1

        with Path(run_file).open("a", newline="\n") as f:
            f.write("\nEND\n")
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
            f.write(f"SPIN {spin}\n")

            if ftype == "gnd":
                f.write(f"FSYM {fsym_gnd}\n")
                f.write(f"ALFA {alpha_gnd}\n")
                f.write(f"BETA {beta_gnd}\n")
                f.write(f"FILE {mofile}\n")
                f.write("END\n")
            elif ftype == "exc":
                f.write(f"FSYM {fsym_exc}\n")
                f.write(f"ALFA {alpha_gnd + 1}\n")
                f.write(f"BETA {beta_gnd}\n")
                f.write("SYM 1\n")
                f.write(f"ALFA {alfa_occ}\n")
                f.write(f"BETA {beta_occ}\n")
                f.write("END\n")
                f.write(f"MULLIKEN {mulliken}\n")
                f.write(f"FILE {mofile}\n")
                f.write("END\n")
            elif ftype == "tp":
                f.write(f"FSYM {fsym_tp}\n")
                f.write(f"ALFA {alpha_gnd}\n")
                f.write(f"BETA {beta_gnd}\n")
                f.write("SYM 1\n")
                f.write(f"ALFA {alfa_occ_tp}\n")
                f.write(f"BETA {beta_occ}\n")
                f.write("END\n")
                f.write("MULLIKEN on full\n")
                f.write(f"FILE {mofile}\n")
                f.write("XRAY xas\n")
                f.write("END\n")
                f.write("END\n")

            for g, ex in atoms_for_basis:
                s = _basis_line(g, ex, ab1a, ab1b, ab2, ab3, ab4, ab5, ab6)
                if s:
                    f.write(s + "\n")
            for g, ex in atoms_for_basis:
                s = _orbital_line(g, ex, ob1a, ob1b, ob2, ob3, ob4, ob5, ob6)
                if s:
                    f.write(s + "\n")
            n_e1 = sum(1 for g, _ in atoms_for_basis if g == 1)
            if mcp1 and n_e1:
                for _ in range(n_e1):
                    f.write(mcp1 + "\n")
            for g, ex in atoms_for_basis:
                if g == 1 and ex:
                    f.write("X-FIRST\n")
                    break
            f.write("END\n")
            f.write("/.\n")
            out_dir = {"gnd": "GND", "exc": "EXC", "tp": "TP"}[ftype]
            f.write(f"mkdir -p ../{out_dir}\n")
            out_line = (
                f"StoBe.x <{aname}{i}{ftype}.inp>& ../{out_dir}/{aname}{i}{ftype}.out\n"
            )
            f.write(out_line)
            f.write(f"mv Molden.molf ../{out_dir}/{aname}{i}{ftype}.molden\n")
            if ftype == "tp":
                f.write(f"mv fort.11 {aname}{i}.xas\n")
            f.write("rm fort.*\n")


def _make_xas_run(mname: str, aname: str, n_atoms: int, title: str) -> None:
    energy_ranges: dict[str, tuple[str, str]] = {
        "C": ("280 320", "0.6 12 288 320"),
        "N": ("400 450", "0.6 12 418 450"),
        "O": ("530 580", "0.6 12 535 560"),
        "Ni": ("830 880", "0.6 12 830 880"),
    }
    elem = aname.upper()
    if elem not in energy_ranges:
        return
    erange, ebroad = energy_ranges[elem]
    for i in range(1, n_atoms + 1):
        run_file = f"{aname}{i}xas.run"
        with Path(run_file).open("w+", newline="\n") as f:
            f.write("#!/bin/bash\n")
            f.write(f'ln -sf "${{PWD}}/{aname}{i}.xas" fort.1\n')
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
            f.write("mkdir -p ../NEXAFS\n")
            f.write(f"{XASRUN} <{aname}{i}xas.inp>& ../NEXAFS/{aname}{i}xas.out\n")
            f.write(f"cp XrayT001.out ../NEXAFS/{aname}{i}.out\n")
            f.write("rm fort.*\n")


def _make_seq_run(n_atoms: int, aname: str) -> None:
    for i in range(1, n_atoms + 1):
        with Path(f"{aname}{i}seq.run").open("w+", newline="\n") as f:
            f.write(f"chmod +x ./{aname}{i}gnd.run\n")
            f.write(f"./{aname}{i}gnd.run\n\n")
            f.write(f"chmod +x ./{aname}{i}exc.run\n")
            f.write(f"./{aname}{i}exc.run\n\n")
            f.write(f"chmod +x ./{aname}{i}tp.run\n")
            f.write(f"./{aname}{i}tp.run\n\n")
            f.write(f"chmod +x ./{aname}{i}xas.run\n")
            f.write(f"./{aname}{i}xas.run\n")


def _list_files(aname: str, n_atoms: int) -> None:
    run_files = list(Path().glob("*.run"))
    nat_sort = natsort_keygen(key=lambda y: y.lower(), alg=ns.IGNORECASE)
    run_files.sort(key=lambda p: nat_sort(p.name))
    for i in range(1, n_atoms + 1):
        atom_dir = Path(f"{aname}{i}")
        pattern = rf"^{re.escape(aname)}{i}(?:\D|$)"
        for rf in run_files:
            if re.match(pattern, rf.name):
                dest = atom_dir / rf.name
                if rf.exists():
                    shutil.move(str(rf), str(dest))


def build_cmd(
    run_directory: str = typer.Argument(..., help="Directory containing molConfig.py"),
    xyz_file: str | None = typer.Argument(None, help="Path to XYZ geometry file"),
    xyz: str | None = typer.Option(None, "--xyz", "-x", help="Alternative XYZ file"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose output"),
) -> None:
    """Generate StoBe run files from molConfig + XYZ."""
    run_dir = Path(run_directory).resolve()
    if not run_dir.exists():
        _CONSOLE.print(f"[red]Error: Directory {run_directory} does not exist[/red]")
        raise typer.Exit(1)

    xyz_path: Path | None = None
    if xyz:
        xyz_path = Path(xyz).resolve()
    elif xyz_file:
        xyz_path = Path(xyz_file).resolve()
    else:
        try:
            xyz_path = _auto_detect_xyz(run_dir).resolve()
        except typer.Exit:
            _CONSOLE.print(f"[red]Error: No .xyz files found in {run_dir}[/red]")
            raise
        if len(_find_xyz_files(run_dir)) > 1:
            _CONSOLE.print(f"[dim]Auto-detected XYZ: {xyz_path}[/dim]")
    if xyz_path is None or not xyz_path.exists():
        _CONSOLE.print("[red]Error: XYZ file does not exist[/red]")
        raise typer.Exit(1)

    try:
        mol_config = load_mol_config(run_dir)
    except (FileNotFoundError, ImportError) as e:
        _CONSOLE.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1) from e

    mname = run_dir.name
    fname = xyz_path.stem

    def g(obj: object, key: str, default: object = None) -> object:
        return getattr(obj, key, default)

    def _int(v: object) -> int:
        if v is None:
            return 0
        return int(v)  # type: ignore[arg-type]

    def _opt_str(v: object) -> str | None:
        return str(v) if v is not None else None

    try:
        aname = mol_config.aname
        n_files = int(mol_config.nFiles)
        dif_elems = int(mol_config.difElems)
        alpha = _int(mol_config.alpha)
        beta = _int(mol_config.beta)
        e1a = _int(mol_config.element1a)
        e1b = _int(mol_config.element1b)
        e2 = _int(g(mol_config, "element2"))
        e3 = _int(g(mol_config, "element3"))
        e4 = _int(g(mol_config, "element4"))
        e5 = _int(g(mol_config, "element5"))
        e6 = _int(g(mol_config, "element6"))
        ne1 = _int(mol_config.nElem1)
        ne2 = _int(g(mol_config, "nElem2"))
        ne3 = _int(g(mol_config, "nElem3"))
        ne4 = _int(g(mol_config, "nElem4"))
        ne5 = _int(g(mol_config, "nElem5"))
        ne6 = _int(g(mol_config, "nElem6"))
        ab1a = mol_config.elem1_Abasis_a
        ab1b = mol_config.elem1_Abasis_b
        ab2 = str(g(mol_config, "elem2_Abasis") or "")
        ab3 = str(g(mol_config, "elem3_Abasis") or "")
        ab4 = str(g(mol_config, "elem4_Abasis") or "")
        ab5 = str(g(mol_config, "elem5_Abasis") or "")
        ab6 = str(g(mol_config, "elem6_Abasis") or "")
        ob1a = mol_config.elem1_Obasis_a
        ob1b = mol_config.elem1_Obasis_b
        ob2 = str(g(mol_config, "elem2_Obasis") or "")
        ob3 = str(g(mol_config, "elem3_Obasis") or "")
        ob4 = str(g(mol_config, "elem4_Obasis") or "")
        ob5 = str(g(mol_config, "elem5_Obasis") or "")
        ob6 = str(g(mol_config, "elem6_Obasis") or "")
        mcp1 = str(g(mol_config, "elem1_MCPbasis") or "")
        alfa_occ = mol_config.alfaOcc
        beta_occ = mol_config.betaOcc
        alfa_occ_tp = mol_config.alfaOccTP
        title = str(g(mol_config, "title") or mname)
        mult = _opt_str(g(mol_config, "multiplicity"))
    except (AttributeError, TypeError, ValueError) as e:
        _CONSOLE.print(f"[red]Error: Missing parameter in molConfig.py: {e}[/red]")
        raise typer.Exit(1) from e

    orig = Path.cwd()
    os.chdir(run_dir)
    try:
        local_xyz = Path(f"{fname}.xyz")
        if not local_xyz.exists():
            shutil.copy2(xyz_path, local_xyz)
        if verbose:
            _CONSOLE.print(f"[dim]XYZ: {local_xyz}, run_dir: {run_dir}[/dim]")

        for ftype in ("gnd", "exc", "tp"):
            _make_run_file(
                fname,
                mname,
                aname,
                n_files,
                alpha,
                beta,
                e1a,
                e1b,
                e2,
                e3,
                e4,
                e5,
                e6,
                ne1,
                ne2,
                ne3,
                ne4,
                ne5,
                ne6,
                dif_elems,
                ab1a,
                ab1b,
                ab2,
                ab3,
                ab4,
                ab5,
                ab6,
                ob1a,
                ob1b,
                ob2,
                ob3,
                ob4,
                ob5,
                ob6,
                mcp1,
                title,
                _opt_str(g(mol_config, "sym")),
                _opt_str(g(mol_config, "geom")),
                _opt_str(g(mol_config, "runtype")),
                _opt_str(g(mol_config, "scftype")),
                _opt_str(g(mol_config, "potential")),
                _opt_str(g(mol_config, "grid")),
                mult,
                _opt_str(g(mol_config, "charge")),
                _opt_str(g(mol_config, "maxcycles")),
                _opt_str(g(mol_config, "econvergence")),
                _opt_str(g(mol_config, "dconvergence")),
                _opt_str(g(mol_config, "dmixing")),
                _opt_str(g(mol_config, "diis")),
                _opt_str(g(mol_config, "orbi")),
                _opt_str(g(mol_config, "mulliken")),
                _opt_str(g(mol_config, "virt")),
                _opt_str(g(mol_config, "spin")),
                _opt_str(g(mol_config, "mofile")),
                _opt_str(g(mol_config, "fsymGND")),
                _opt_str(g(mol_config, "fsymEXC")),
                _opt_str(g(mol_config, "fsymEXC")),
                alfa_occ,
                alfa_occ_tp,
                beta_occ,
                ftype,
            )
        _make_xas_run(mname, aname, n_files, title)
        _make_seq_run(n_files, aname)
        for i in range(1, n_files + 1):
            Path(f"{aname}{i}").mkdir(parents=True, exist_ok=True)
        _list_files(aname, n_files)

        msg = (
            f"Generated StoBe input files for {n_files} {aname} atoms.\n"
            f"Output: {run_dir}"
        )
        _CONSOLE.print(Panel(msg, title="dftrun build", border_style="blue"))
    finally:
        os.chdir(orig)
