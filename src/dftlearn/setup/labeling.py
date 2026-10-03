"""Site-group merging, relabeling, and conversion to geometry plans."""

from __future__ import annotations

from dataclasses import dataclass, replace

from dftlearn.setup.types import DistinctAtomGroup


@dataclass
class SiteGroup:
    """User-editable equivalence class for one core-edge generating site."""

    group_id: int
    rank: int
    element: str
    atom_indices: tuple[int, ...]
    representative_index: int
    cif_label: str
    xyz_label: str
    site_tag: str
    custom_name: str
    enabled: bool = True


def site_groups_from_distinct(groups: list[DistinctAtomGroup]) -> list[SiteGroup]:
    """Build editable site groups from auto-detected distinct classes."""
    site_groups: list[SiteGroup] = []
    for idx, group in enumerate(groups, start=1):
        tag = f"{group.element}{idx}"
        site_groups.append(
            SiteGroup(
                group_id=idx,
                rank=group.rank,
                element=group.element,
                atom_indices=group.atom_indices,
                representative_index=group.representative_index,
                cif_label=group.cif_label,
                xyz_label=group.xyz_label,
                site_tag=tag,
                custom_name=group.cif_label,
                enabled=True,
            )
        )
    return site_groups


def active_generating_groups(site_groups: list[SiteGroup]) -> list[DistinctAtomGroup]:
    """Convert enabled site groups into distinct groups for geometry planning."""
    active = [g for g in site_groups if g.enabled]
    active.sort(key=lambda g: g.group_id)
    result: list[DistinctAtomGroup] = []
    for group in active:
        result.append(
            DistinctAtomGroup(
                rank=group.rank,
                element=group.element,
                atom_indices=group.atom_indices,
                representative_index=group.representative_index,
                cif_label=group.cif_label,
                xyz_label=group.xyz_label,
            )
        )
    return result


def merge_site_groups(
    site_groups: list[SiteGroup],
    group_ids: list[int],
) -> list[SiteGroup]:
    """Merge multiple site groups into the first selected group.

    Parameters
    ----------
    site_groups
        Mutable list of site groups to update in place logically.
    group_ids
        Two or more ``group_id`` values to merge.

    Returns
    -------
    list[SiteGroup]
        New list with merged groups disabled except the primary.

    Raises
    ------
    ValueError
        If fewer than two groups are supplied or IDs are missing.
    """
    if len(group_ids) < 2:
        msg = "Select at least two groups to merge"
        raise ValueError(msg)
    by_id = {g.group_id: g for g in site_groups}
    missing = [gid for gid in group_ids if gid not in by_id]
    if missing:
        msg = f"Unknown group ids: {missing}"
        raise ValueError(msg)
    primary_id = group_ids[0]
    primary = by_id[primary_id]
    merged_indices: list[int] = list(primary.atom_indices)
    for gid in group_ids[1:]:
        other = by_id[gid]
        if other.element != primary.element:
            msg = "Cannot merge groups with different elements"
            raise ValueError(msg)
        merged_indices.extend(other.atom_indices)
    merged_tuple = tuple(sorted(set(merged_indices)))
    rep = min(merged_tuple)
    updated: list[SiteGroup] = []
    for group in site_groups:
        if group.group_id == primary_id:
            updated.append(
                replace(
                    group,
                    atom_indices=merged_tuple,
                    representative_index=rep,
                    enabled=True,
                )
            )
        elif group.group_id in group_ids[1:]:
            updated.append(replace(group, enabled=False))
        else:
            updated.append(group)
    return updated


def relabel_site_group(
    site_groups: list[SiteGroup],
    group_id: int,
    *,
    site_tag: str | None = None,
    custom_name: str | None = None,
) -> list[SiteGroup]:
    """Update the StoBe site tag or custom name for one site group."""
    updated: list[SiteGroup] = []
    found = False
    for group in site_groups:
        if group.group_id != group_id:
            updated.append(group)
            continue
        found = True
        updated.append(
            replace(
                group,
                site_tag=site_tag if site_tag is not None else group.site_tag,
                custom_name=custom_name
                if custom_name is not None
                else group.custom_name,
            )
        )
    if not found:
        msg = f"Unknown group id {group_id}"
        raise ValueError(msg)
    return updated


def restrict_site_groups_to_atom_indices(
    site_groups: list[SiteGroup],
    atom_indices: tuple[int, ...] | set[int],
    edge_element: str,
) -> list[SiteGroup]:
    """Enable only site groups whose atoms intersect ``atom_indices``.

    Parameters
    ----------
    site_groups
        Current generating-site groups.
    atom_indices
        Atom indices of the ligand (or other subset) to keep as core sites.
    edge_element
        Core-edge element used to renumber enabled tags (``C1``, ``C2``, ...).

    Returns
    -------
    list[SiteGroup]
        Copy of ``site_groups`` with other groups disabled and tags renumbered.
    """
    keep = set(atom_indices)
    updated = [
        replace(
            group,
            enabled=bool(keep.intersection(group.atom_indices)),
        )
        for group in site_groups
    ]
    return renumber_site_tags(updated, edge_element)


def toggle_site_group(
    site_groups: list[SiteGroup],
    group_id: int,
    *,
    enabled: bool,
) -> list[SiteGroup]:
    """Enable or disable one site group as a generating core site."""
    return [
        replace(group, enabled=enabled if group.group_id == group_id else group.enabled)
        for group in site_groups
    ]


def renumber_site_tags(
    site_groups: list[SiteGroup],
    edge_element: str,
) -> list[SiteGroup]:
    """Renumber enabled site tags as ``C1``, ``C2``, ... in group-id order."""
    edge = edge_element.strip().capitalize()
    enabled = sorted(
        (g for g in site_groups if g.enabled),
        key=lambda g: g.group_id,
    )
    tag_map = {g.group_id: f"{edge}{idx}" for idx, g in enumerate(enabled, start=1)}
    return [
        replace(group, site_tag=tag_map.get(group.group_id, group.site_tag))
        for group in site_groups
    ]
