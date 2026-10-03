"""dftrun init: scaffold StoBe runs from CIF structures."""

from __future__ import annotations

import json
from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from dftlearn.setup import cif_blocks
from dftlearn.setup.remote import KNOWN_REMOTE_HOSTS
from dftlearn.setup.session import (
    create_setup_session,
    create_setup_session_from_rows,
    save_setup_session,
)
from dftlearn.setup.workflow import finalize_setup_session

_CONSOLE = Console()


def init_cmd(
    run_directory: Path = typer.Argument(
        ...,
        help="Directory to create or populate with setup_session.json.",
        exists=False,
        file_okay=False,
        dir_okay=True,
        writable=True,
    ),
    cif: Path | None = typer.Option(
        None,
        "--cif",
        help="Input CIF file (single-molecule or asymmetric unit).",
        exists=True,
        dir_okay=False,
        readable=True,
    ),
    pubchem: str | None = typer.Option(
        None,
        "--pubchem",
        help="PubChem name or formula; uses the first hit (e.g. C27H18AlN3O3).",
    ),
    name: str | None = typer.Option(
        None,
        "--name",
        help="Molecule title for molConfig and figures; defaults to run directory name.",
    ),
    edge_element: str = typer.Option(
        "C",
        "--edge-element",
        help="Core-edge element for StoBe site folders (e.g. C for K-edge).",
    ),
    block: str | None = typer.Option(
        None,
        "--block",
        help="CIF data block name; auto-selects the largest atom-site block when omitted.",
    ),
    remote: str | None = typer.Option(
        None,
        "--remote",
        help=f"Sync run directory to SSH host alias ({', '.join(sorted(KNOWN_REMOTE_HOSTS))}).",
    ),
    remote_root: str = typer.Option(
        "/home/hduva/projects/dft-runs",
        "--remote-root",
        help="Remote base directory when --remote is set.",
    ),
    no_preview: bool = typer.Option(
        False,
        "--no-preview",
        help="Skip writing init_preview.png.",
    ),
    list_blocks: bool = typer.Option(
        False,
        "--list-blocks",
        help="List CIF blocks and exit without writing a run directory.",
    ),
    interactive: bool = typer.Option(
        False,
        "--interactive",
        "-i",
        help="Open the browser setup lab after creating the session.",
    ),
    batch: bool = typer.Option(
        False,
        "--batch",
        help="Write geometry and molConfig immediately without interactive review.",
    ),
    keep_all_sites: bool = typer.Option(
        False,
        "--keep-all-sites",
        help="Keep residual solvent/synthesis fragments from the CIF (default: drop them).",
    ),
    one_ligand: bool = typer.Option(
        False,
        "--one-ligand",
        help="Enable core-edge sites on a single metal-bound organic wing only.",
    ),
) -> None:
    """Initialize a StoBe run from a CIF file or PubChem compound.

    By default creates ``setup_session.json`` for interactive review with
    ``dftrun setup``. Use ``--batch`` to write ``geometry.xyz`` and
    ``molConfig.py`` immediately, or ``--interactive`` to open the browser lab.
    """
    if list_blocks:
        if cif is None:
            typer.echo("init failed: --list-blocks requires --cif", err=True)
            raise typer.Exit(1)
        blocks = cif_blocks(cif)
        if not blocks:
            typer.echo(f"No atom-site blocks in {cif}")
            raise typer.Exit(1)
        table = Table(title=f"CIF blocks in {cif.name}")
        table.add_column("Block")
        table.add_column("Sites", justify="right")
        for block_name, count in blocks:
            table.add_row(block_name, str(count))
        _CONSOLE.print(table)
        raise typer.Exit(0)

    if batch and interactive:
        typer.echo(
            "init failed: use either --batch or --interactive, not both", err=True
        )
        raise typer.Exit(1)

    if cif is None and pubchem is None:
        typer.echo("init failed: provide --cif or --pubchem", err=True)
        raise typer.Exit(1)
    if cif is not None and pubchem is not None:
        typer.echo("init failed: use either --cif or --pubchem, not both", err=True)
        raise typer.Exit(1)

    mname = name or run_directory.name
    source_label = cif.name if cif is not None else f"PubChem:{pubchem}"
    try:
        if pubchem is not None:
            session = _session_from_pubchem(
                pubchem,
                mname=mname,
                edge_element=edge_element,
                one_ligand=one_ligand,
            )
        else:
            assert cif is not None
            session = create_setup_session(
                cif,
                mname=mname,
                edge_element=edge_element,
                block_name=block,
                keep_all_sites=keep_all_sites,
                one_ligand=one_ligand,
            )
    except (ValueError, FileNotFoundError) as exc:
        typer.echo(f"init failed: {exc}", err=True)
        raise typer.Exit(1) from exc

    run_directory.mkdir(parents=True, exist_ok=True)
    save_setup_session(run_directory, session)

    if interactive:
        from dftlearn.cli.setup import _launch_streamlit

        save_setup_session(run_directory, session)
        _launch_streamlit(run_directory, host="localhost", port=8501)
        return

    if batch:
        try:
            artifacts = finalize_setup_session(
                run_directory,
                session,
                write_preview=not no_preview,
                remote_host=remote,
                remote_root=remote_root,
            )
        except (ValueError, KeyError, RuntimeError) as exc:
            typer.echo(f"init failed: {exc}", err=True)
            raise typer.Exit(1) from exc
        summary = json.loads(artifacts.summary_json.read_text(encoding="utf-8"))
        _print_init_summary(
            run_directory, source_label, edge_element, summary, artifacts
        )
        return

    _CONSOLE.print(
        f"Draft session written to {run_directory / 'setup_session.json'}\n"
        f"  {edge_element} site groups: {len(session.site_groups)}\n"
        "Run `dftrun setup` (or `dftrun setup RUN_DIR`) to open the browser lab:\n"
        "  PubChem search -> 3D embed -> relaxation -> site labeler"
    )


def _session_from_pubchem(
    query: str,
    *,
    mname: str,
    edge_element: str,
    one_ligand: bool = False,
):
    """Build a setup session from the first PubChem hit for ``query``."""
    from dftlearn.setup.pubchem_client import search_pubchem
    from dftlearn.setup.structure_build import RelaxMethod, build_3d_from_smiles

    hits = search_pubchem(query, limit=1)
    if not hits:
        msg = f"No PubChem hits for {query!r}"
        raise ValueError(msg)
    compound = hits[0]
    built = build_3d_from_smiles(
        compound.smiles,
        relax=RelaxMethod.ANNEAL,
        relax_steps=500,
    )
    meta = {
        "source_kind": "pubchem",
        "pubchem_cid": compound.cid,
        "pubchem_title": compound.title,
        "smiles": compound.smiles,
        "molecular_formula": compound.molecular_formula,
        "relax_method": RelaxMethod.ANNEAL.value,
        "relax_steps": 500,
        "final_energy": built.final_energy,
    }
    return create_setup_session_from_rows(
        built.rows,
        mname=mname,
        edge_element=edge_element,
        atom_labels=built.atom_labels,
        source_meta=meta,
        mol=built.mol,
        one_ligand=one_ligand,
    )


def _print_init_summary(
    run_directory: Path,
    source_label: str,
    edge_element: str,
    summary: dict,
    artifacts,
) -> None:
    n_sites = summary["n_core_sites"]
    sg = summary.get("space_group", "?")
    _CONSOLE.print(
        f"Initialized {run_directory} from {source_label}\n"
        f"  block: {summary['block_name']}\n"
        f"  space group: {sg}\n"
        f"  {edge_element} core sites: {n_sites}\n"
        f"  geometry: {run_directory / 'geometry.xyz'}\n"
        f"  molConfig: {run_directory / 'molConfig.py'}"
    )
    if artifacts is not None and artifacts.preview_png:
        _CONSOLE.print(f"  preview: {artifacts.preview_png}")
    if artifacts is not None and artifacts.remote_directory and artifacts.remote_host:
        _CONSOLE.print(
            f"  remote: {artifacts.remote_host}:{artifacts.remote_directory}"
        )

    site_table = Table(title=f"Distinct {edge_element} generating sites")
    site_table.add_column("Site")
    site_table.add_column("Name")
    site_table.add_column("XYZ")
    site_table.add_column("CIF")
    site_table.add_column("Class", justify="right")
    for row in summary.get("distinct_sites", []):
        site_table.add_row(
            str(row.get("site", "")),
            str(row.get("custom_name", "")),
            str(row.get("xyz_label", "")),
            str(row.get("cif_label", "")),
            str(row.get("class_size", "")),
        )
    _CONSOLE.print(site_table)
