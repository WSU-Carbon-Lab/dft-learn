"""Energetic-overlap clustering of NEXAFS transitions.

This package implements the Igor overlap-merge core: percent overlap of
Gaussians, ``simpleCluster3`` grouping, envelope effective clusters, an OS%
elbow cutoff, and overlap-threshold selection by a Normal(50%, 10%) sigma
ladder plus a Gaussian process on BIC. It does not classify transition
symmetry, refit amplitudes to experiment, or build film tensors.
"""

from __future__ import annotations

from dftlearn.clustering.grouping import simple_cluster3
from dftlearn.clustering.iterate import cluster_by_overlap
from dftlearn.clustering.merge import effective_gaussian, merge_groups
from dftlearn.clustering.os_elbow import os_percent_elbow
from dftlearn.clustering.overlap import (
    FIRST_OVERLAP_FWHM_EV,
    FWHM_TO_SIGMA,
    overlap_matrix,
    percent_overlap,
)
from dftlearn.clustering.selection import (
    cluster_bic,
    select_overlap_threshold,
)
from dftlearn.clustering.types import (
    CLUSTERING_SPEC,
    ClusterResult,
    ClusterStage,
    OverlapSearchResult,
    TransitionSticks,
    sticks_from_tp_table,
)

__all__ = [
    "CLUSTERING_SPEC",
    "FIRST_OVERLAP_FWHM_EV",
    "FWHM_TO_SIGMA",
    "ClusterResult",
    "ClusterStage",
    "OverlapSearchResult",
    "TransitionSticks",
    "cluster_bic",
    "cluster_by_overlap",
    "effective_gaussian",
    "merge_groups",
    "os_percent_elbow",
    "overlap_matrix",
    "percent_overlap",
    "select_overlap_threshold",
    "simple_cluster3",
    "sticks_from_tp_table",
]
