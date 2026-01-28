"""Sklearn-compatible clustering algorithms for spectral peaks."""

import numpy as np
import pandas as pd
from sklearn.utils.validation import check_is_fitted

from dftlearn._backends import get_backend
from dftlearn.analysis.overlap import OverlapCalculator
from dftlearn.pipeline.estimators import SpectrumClusterer


class OverlapClusterer(SpectrumClusterer):
    """
    Cluster spectral peaks by Gaussian overlap.

    This sklearn-compatible clusterer groups peaks based on their overlap
    percentage. It first calculates an overlap matrix, then applies either
    sequential or skip-tolerant clustering to identify groups of overlapping peaks.

    Parameters
    ----------
    overlap_threshold : float, default=50.0
        Minimum overlap percentage (0-100) to merge peaks into same cluster.
    method : {'sequential', 'skip_tolerant'}, default='sequential'
        Clustering algorithm:
        - 'sequential': Group peaks by overlap with first peak in cluster
        - 'skip_tolerant': Allow gaps in overlap sequence
    n_skipped : int, default=1
        For skip_tolerant method: maximum number of peaks below threshold
        to skip before terminating cluster.
    backend : {'auto', 'cpu', 'gpu'}, default='auto'
        Computation backend for overlap matrix calculation.

    Attributes
    ----------
    labels_ : ndarray, shape (n_samples,)
        Cluster label for each peak (set during fit).
    n_clusters_ : int
        Number of clusters found (set during fit).
    overlap_matrix_ : ndarray, shape (n_samples, n_samples)
        Computed overlap matrix (set during fit).
    cluster_info_ : list of lists
        Peak indices grouped by cluster (set during fit).

    Examples
    --------
    >>> import pandas as pd
    >>> import numpy as np
    >>> df = pd.DataFrame({
    ...     'E': [285.0, 285.5, 290.0],
    ...     'width': [0.5, 0.5, 0.5],
    ...     'OS': [1.0, 0.8, 0.5]
    ... })
    >>> clusterer = OverlapClusterer(overlap_threshold=50.0)
    >>> clusterer.fit(df)
    OverlapClusterer(overlap_threshold=50.0)
    >>> clusterer.n_clusters_
    2
    >>> clusterer.labels_
    array([0, 0, 1])

    Notes
    -----
    This clusterer is designed for spectroscopy data where peaks are represented
    as Gaussian functions with energy, width, and amplitude parameters.

    The sequential method groups peaks in energy order, adding each peak to the
    current cluster if it overlaps with the first peak in that cluster.

    The skip_tolerant method allows skipping up to n_skipped peaks below the
    threshold before terminating a cluster, enabling capture of peaks that are
    adjacent but not directly overlapping.
    """

    def __init__(
        self,
        overlap_threshold: float = 50.0,
        method: str = "sequential",
        n_skipped: int = 1,
        backend: str = "auto",
    ):
        """
        Initialize overlap-based clusterer.

        Parameters
        ----------
        overlap_threshold : float, default=50.0
            Minimum overlap percentage to merge peaks.
        method : {'sequential', 'skip_tolerant'}, default='sequential'
            Clustering algorithm.
        n_skipped : int, default=1
            Number of peaks to skip (skip_tolerant only).
        backend : {'auto', 'cpu', 'gpu'}, default='auto'
            Computation backend.
        """
        self.overlap_threshold = overlap_threshold
        self.method = method
        self.n_skipped = n_skipped
        self.backend = backend

    def fit(self, X, y=None):
        """
        Compute overlap matrix and perform clustering.

        Parameters
        ----------
        X : DataFrame
            Peak data with columns 'E' (energy), 'width', and 'OS' (oscillator strength).
        y : array-like, optional
            Ignored. Present for sklearn API compatibility.

        Returns
        -------
        self : object
            Fitted clusterer instance.

        Raises
        ------
        ValueError
            If required columns are missing or method is invalid.
        """
        # Validate input
        if not isinstance(X, pd.DataFrame):
            raise TypeError(
                f"X must be a pandas DataFrame, got {type(X).__name__}. "
                f"Required columns: 'E', 'width', 'OS'"
            )

        required_cols = ["E", "width", "OS"]
        missing_cols = [col for col in required_cols if col not in X.columns]
        if missing_cols:
            raise ValueError(
                f"Missing required columns: {missing_cols}. "
                f"Available columns: {list(X.columns)}"
            )

        if len(X) == 0:
            raise ValueError("Cannot cluster empty DataFrame")

        # Calculate overlap matrix
        overlap_calc = OverlapCalculator(backend=self.backend)
        result = overlap_calc.calculate_from_dataframe(X)
        self.overlap_matrix_ = result.overlap_matrix

        # Sort peaks by energy
        sorted_indices = np.argsort(X["E"].values)

        # Get backend for clustering
        backend_instance = get_backend(self.backend)

        # Perform clustering
        if self.method == "sequential":
            clusters = backend_instance.sequential_clustering(
                self.overlap_matrix_, sorted_indices, self.overlap_threshold
            )
        elif self.method == "skip_tolerant":
            clusters = backend_instance.skip_tolerant_clustering(
                self.overlap_matrix_,
                sorted_indices,
                self.overlap_threshold,
                self.n_skipped,
            )
        else:
            raise ValueError(
                f"Unknown clustering method: {self.method!r}. "
                f"Expected 'sequential' or 'skip_tolerant'"
            )

        # Store cluster information
        self.cluster_info_ = clusters
        self.n_clusters_ = len(clusters)

        # Convert clusters to labels array
        self.labels_ = np.zeros(len(X), dtype=np.int32)
        for cluster_id, peak_indices in enumerate(clusters):
            for peak_idx in peak_indices:
                self.labels_[peak_idx] = cluster_id

        return self

    def predict(self, X):
        """
        Return cluster labels for fitted data.

        Parameters
        ----------
        X : DataFrame
            Peak data (must be same as fitted data).

        Returns
        -------
        labels : ndarray, shape (n_samples,)
            Cluster label for each peak.

        Notes
        -----
        This method returns the labels computed during fit(). It does not
        predict labels for new data, as overlap-based clustering requires
        recalculating the full overlap matrix.

        To cluster new data, call fit() again.
        """
        check_is_fitted(self, ["labels_"])
        return self.labels_

    def fit_predict(self, X, y=None):
        """
        Fit clusterer and return labels.

        Parameters
        ----------
        X : DataFrame
            Peak data with columns 'E', 'width', 'OS'.
        y : array-like, optional
            Ignored.

        Returns
        -------
        labels : ndarray, shape (n_samples,)
            Cluster label for each peak.
        """
        return self.fit(X, y).labels_

    def get_cluster_peaks(self, X):
        """
        Get peaks grouped by cluster.

        Parameters
        ----------
        X : DataFrame
            Peak data (must be same as fitted data).

        Returns
        -------
        clusters : list of DataFrames
            Each DataFrame contains peaks from one cluster.

        Examples
        --------
        >>> clusterer = OverlapClusterer()
        >>> clusterer.fit(df)  # doctest: +SKIP
        >>> cluster_peaks = clusterer.get_cluster_peaks(df)  # doctest: +SKIP
        >>> len(cluster_peaks)  # Number of clusters  # doctest: +SKIP
        2
        """
        check_is_fitted(self, ["cluster_info_"])

        clusters = []
        for peak_indices in self.cluster_info_:
            cluster_df = X.iloc[peak_indices].copy()
            clusters.append(cluster_df)

        return clusters

    def get_cluster_centers(self):
        """
        Get representative energy for each cluster.

        Returns
        -------
        centers : ndarray, shape (n_clusters,)
            Mean energy of peaks in each cluster.

        Examples
        --------
        >>> clusterer = OverlapClusterer()
        >>> clusterer.fit(df)  # doctest: +SKIP
        >>> centers = clusterer.get_cluster_centers()  # doctest: +SKIP
        """
        check_is_fitted(self, ["cluster_info_", "overlap_matrix_"])

        # This requires the original data, which we don't store
        # Raise informative error
        raise NotImplementedError(
            "get_cluster_centers() requires storing original data. "
            "Use get_cluster_peaks(X) instead to retrieve cluster information."
        )

    def __repr__(self) -> str:
        """Return string representation."""
        return (
            f"{self.__class__.__name__}("
            f"overlap_threshold={self.overlap_threshold}, "
            f"method={self.method!r}, "
            f"n_skipped={self.n_skipped}, "
            f"backend={self.backend!r})"
        )
