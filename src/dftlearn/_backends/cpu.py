"""CPU backend using NumPy and Numba for accelerated computation."""

import numpy as np
from numba import jit, prange


class CPUBackend:
    """
    CPU-accelerated backend using Numba JIT compilation.

    This backend provides optimized implementations of core algorithms using
    Numba's just-in-time compilation and parallel execution capabilities.

    Methods
    -------
    calculate_overlap_matrix(energies, widths, amplitudes)
        Calculate overlap matrix between Gaussian peaks.
    sequential_clustering(overlap_matrix, sorted_indices, threshold)
        Perform sequential clustering based on overlap with first peak.
    skip_tolerant_clustering(overlap_matrix, sorted_indices, threshold, n_skipped)
        Perform skip-tolerant clustering allowing gaps.

    Notes
    -----
    All methods use Numba's @jit decorator with nopython=True for maximum
    performance. The overlap matrix calculation uses parallel execution
    via prange for additional speedup on multi-core CPUs.
    """

    @staticmethod
    @jit(nopython=True, parallel=True, fastmath=True)
    def _gaussian_numba(x, amp, center, sigma):
        """
        Fast Gaussian function optimized with Numba.

        Parameters
        ----------
        x : ndarray
            Energy points.
        amp : float
            Amplitude (oscillator strength).
        center : float
            Peak center (energy).
        sigma : float
            Peak width (standard deviation).

        Returns
        -------
        ndarray
            Gaussian values at x.
        """
        normalization = 1.0 / (sigma * np.sqrt(2.0 * np.pi))
        return amp * normalization * np.exp(-0.5 * ((x - center) / sigma) ** 2)

    @staticmethod
    @jit(nopython=True, parallel=True, fastmath=True)
    def _calculate_overlap_matrix_numba(energies, widths, amplitudes):
        """
        Calculate overlap matrix between Gaussian peaks using Numba.

        Optimized implementation using:
        - Parallel execution (prange)
        - Fast math operations
        - Early termination for distant peaks
        - Vectorized Gaussian calculations
        - Trapezoidal integration

        Parameters
        ----------
        energies : ndarray, shape (n,)
            Peak centers.
        widths : ndarray, shape (n,)
            Peak widths (standard deviations).
        amplitudes : ndarray, shape (n,)
            Peak amplitudes (oscillator strengths).

        Returns
        -------
        overlap_matrix : ndarray, shape (n, n)
            Symmetric matrix where overlap_matrix[i, j] is the percentage
            overlap between peaks i and j.
        """
        n = len(energies)
        overlap_matrix = np.zeros((n, n), dtype=np.float64)

        # Set diagonal to 100%
        for i in range(n):
            overlap_matrix[i, i] = 100.0

        # Calculate upper triangle (exploit symmetry)
        for i in prange(n):
            for j in range(i + 1, n):
                area1 = amplitudes[i]
                area2 = amplitudes[j]
                min_area = min(area1, area2)

                if min_area > 0:
                    center1, sigma1, amp1 = energies[i], widths[i], amplitudes[i]
                    center2, sigma2, amp2 = energies[j], widths[j], amplitudes[j]

                    separation = abs(center1 - center2)
                    combined_width = (sigma1 + sigma2) * 0.5

                    # Early termination for distant peaks (>8 sigma apart)
                    if separation > 8.0 * combined_width:
                        percent_overlap = 0.0
                    else:
                        # Optimized integration bounds
                        sigma_min = min(sigma1, sigma2)
                        x_min = min(center1 - 5.0 * sigma1, center2 - 5.0 * sigma2)
                        x_max = max(center1 + 5.0 * sigma1, center2 + 5.0 * sigma2)

                        # Adaptive number of integration points
                        n_points = min(
                            10000, max(1000, int((x_max - x_min) / (sigma_min * 0.1)))
                        )
                        dx = (x_max - x_min) / (n_points - 1)

                        # Vectorized energy grid
                        x_vals = x_min + np.arange(n_points) * dx

                        # Pre-compute normalization constants
                        norm1 = amp1 / (sigma1 * np.sqrt(2.0 * np.pi))
                        norm2 = amp2 / (sigma2 * np.sqrt(2.0 * np.pi))
                        inv_2sigma1_sq = -0.5 / (sigma1 * sigma1)
                        inv_2sigma2_sq = -0.5 / (sigma2 * sigma2)

                        # Vectorized Gaussian calculations
                        diff1 = x_vals - center1
                        diff2 = x_vals - center2
                        g1_vals = norm1 * np.exp(inv_2sigma1_sq * diff1 * diff1)
                        g2_vals = norm2 * np.exp(inv_2sigma2_sq * diff2 * diff2)

                        # Minimum envelope (overlap region)
                        overlap_vals = np.minimum(g1_vals, g2_vals)

                        # Trapezoidal integration
                        overlap_area = np.trapz(overlap_vals, dx=dx)

                        # Percentage of smaller peak
                        percent_overlap = max(
                            0.0, min(100.0, (overlap_area / min_area) * 100.0)
                        )

                    overlap_matrix[i, j] = percent_overlap
                    overlap_matrix[j, i] = percent_overlap  # Symmetric
                else:
                    overlap_matrix[i, j] = 0.0
                    overlap_matrix[j, i] = 0.0

        return overlap_matrix

    def calculate_overlap_matrix(self, energies, widths, amplitudes):
        """
        Calculate overlap matrix between Gaussian peaks.

        Parameters
        ----------
        energies : array-like, shape (n,)
            Peak centers (eV).
        widths : array-like, shape (n,)
            Peak widths/standard deviations (eV).
        amplitudes : array-like, shape (n,)
            Peak amplitudes/oscillator strengths.

        Returns
        -------
        overlap_matrix : ndarray, shape (n, n)
            Symmetric matrix of percentage overlaps.

        Examples
        --------
        >>> backend = CPUBackend()
        >>> energies = np.array([285.0, 286.0, 290.0])
        >>> widths = np.array([0.5, 0.6, 0.5])
        >>> amplitudes = np.array([1.0, 0.8, 0.5])
        >>> overlap = backend.calculate_overlap_matrix(energies, widths, amplitudes)
        >>> overlap[0, 0]
        100.0
        >>> overlap[0, 2] < 1.0  # Distant peaks have low overlap
        True
        """
        return self._calculate_overlap_matrix_numba(
            np.asarray(energies, dtype=np.float64),
            np.asarray(widths, dtype=np.float64),
            np.asarray(amplitudes, dtype=np.float64),
        )

    @staticmethod
    @jit(nopython=True)
    def _sequential_clustering_numba(overlap_matrix, sorted_indices, threshold):
        """
        Sequential clustering checking overlap with FIRST peak in each cluster.

        Parameters
        ----------
        overlap_matrix : ndarray, shape (n, n)
            Overlap percentages between peaks.
        sorted_indices : ndarray, shape (n,)
            Peak indices sorted by energy.
        threshold : float
            Minimum overlap percentage to merge peaks.

        Returns
        -------
        clusters : list of lists
            Each inner list contains indices of peaks in that cluster.
        """
        clusters = [[sorted_indices[0]]]  # Start with first peak
        current_cluster_idx = 0

        for i in range(1, len(sorted_indices)):
            current_peak_idx = sorted_indices[i]
            first_peak_idx = clusters[current_cluster_idx][0]  # Check with FIRST peak

            # Check overlap with first peak in current cluster
            overlap = overlap_matrix[current_peak_idx, first_peak_idx]

            if overlap >= threshold:
                # Add to current cluster
                clusters[current_cluster_idx].append(current_peak_idx)
            else:
                # Start new cluster
                clusters.append([current_peak_idx])
                current_cluster_idx += 1

        return clusters

    def sequential_clustering(self, overlap_matrix, sorted_indices, threshold):
        """
        Perform sequential clustering based on overlap with first peak.

        This method groups peaks sequentially by energy, adding each peak to
        the current cluster if it overlaps with the first peak in that cluster.

        Parameters
        ----------
        overlap_matrix : ndarray, shape (n, n)
            Overlap percentages between peaks.
        sorted_indices : ndarray, shape (n,)
            Peak indices sorted by energy.
        threshold : float
            Minimum overlap percentage to merge peaks (0-100).

        Returns
        -------
        clusters : list of lists
            Each inner list contains indices of peaks in that cluster.

        Examples
        --------
        >>> backend = CPUBackend()
        >>> overlap = np.array([[100, 60, 10], [60, 100, 15], [10, 15, 100]])
        >>> indices = np.array([0, 1, 2])
        >>> clusters = backend.sequential_clustering(overlap, indices, 50.0)
        >>> len(clusters)
        2
        >>> clusters[0]  # First two peaks cluster together
        [0, 1]
        >>> clusters[1]  # Third peak separate
        [2]
        """
        return self._sequential_clustering_numba(overlap_matrix, sorted_indices, threshold)

    @staticmethod
    @jit(nopython=True)
    def _skip_tolerant_clustering_numba(
        overlap_matrix, sorted_indices, threshold, n_skipped
    ):
        """
        Skip-tolerant clustering allowing up to n_skipped peaks below threshold.

        Parameters
        ----------
        overlap_matrix : ndarray, shape (n, n)
            Overlap percentages between peaks.
        sorted_indices : ndarray, shape (n,)
            Peak indices sorted by energy.
        threshold : float
            Minimum overlap percentage to merge peaks.
        n_skipped : int
            Number of peaks below threshold to skip before terminating cluster.

        Returns
        -------
        clusters : list of lists
            Each inner list contains indices of peaks in that cluster.
        """
        n_peaks = len(sorted_indices)
        clustered = [False] * n_peaks  # Track which peaks are clustered
        clusters = []

        i = 0
        while i < n_peaks:
            # Find next unclustered peak
            while i < n_peaks and clustered[i]:
                i += 1

            if i >= n_peaks:
                break

            # Start new cluster with this peak
            current_cluster = [sorted_indices[i]]
            first_peak_idx = sorted_indices[i]
            clustered[i] = True

            # Check subsequent peaks for inclusion
            skip_count = 0
            j = i + 1

            while j < n_peaks and skip_count <= n_skipped:
                if clustered[j]:
                    j += 1
                    continue

                current_peak_idx = sorted_indices[j]
                overlap = overlap_matrix[current_peak_idx, first_peak_idx]

                if overlap >= threshold:
                    # Add to cluster and reset skip count
                    current_cluster.append(current_peak_idx)
                    clustered[j] = True
                    skip_count = 0
                else:
                    # Skip this peak
                    skip_count += 1
                    if skip_count > n_skipped:
                        break

                j += 1

            clusters.append(current_cluster)
            i += 1

        return clusters

    def skip_tolerant_clustering(
        self, overlap_matrix, sorted_indices, threshold, n_skipped=1
    ):
        """
        Perform skip-tolerant clustering allowing gaps in overlap.

        This method allows skipping up to n_skipped peaks below the overlap
        threshold before terminating a cluster, enabling capture of peaks that
        are adjacent but not directly overlapping.

        Parameters
        ----------
        overlap_matrix : ndarray, shape (n, n)
            Overlap percentages between peaks.
        sorted_indices : ndarray, shape (n,)
            Peak indices sorted by energy.
        threshold : float
            Minimum overlap percentage to merge peaks (0-100).
        n_skipped : int, default=1
            Number of peaks below threshold to skip before terminating cluster.

        Returns
        -------
        clusters : list of lists
            Each inner list contains indices of peaks in that cluster.

        Examples
        --------
        >>> backend = CPUBackend()
        >>> # Peaks: 0-1 overlap, 1-2 don't, 2-3 overlap with 0
        >>> overlap = np.array([
        ...     [100, 60, 10, 55],
        ...     [60, 100, 5, 8],
        ...     [10, 5, 100, 12],
        ...     [55, 8, 12, 100]
        ... ])
        >>> indices = np.array([0, 1, 2, 3])
        >>> clusters = backend.skip_tolerant_clustering(overlap, indices, 50.0, n_skipped=2)
        >>> # Peak 2 is skipped, but peak 3 is included due to overlap with peak 0
        >>> 3 in clusters[0]
        True
        """
        return self._skip_tolerant_clustering_numba(
            overlap_matrix, sorted_indices, threshold, n_skipped
        )
