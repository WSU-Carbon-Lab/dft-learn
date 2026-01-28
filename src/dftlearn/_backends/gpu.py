"""GPU backend using CuPy for CUDA-accelerated computation."""

import cupy as cp
import numpy as np


class GPUBackend:
    """
    GPU-accelerated backend using CuPy for CUDA computation.

    This backend provides GPU-optimized implementations using CuPy, which offers
    a NumPy-compatible interface with CUDA acceleration. For algorithms that
    don't benefit from GPU parallelization (e.g., sequential clustering), this
    backend falls back to CPU implementations.

    Methods
    -------
    calculate_overlap_matrix(energies, widths, amplitudes)
        Calculate overlap matrix using GPU acceleration.
    sequential_clustering(overlap_matrix, sorted_indices, threshold)
        Sequential clustering (falls back to CPU).
    skip_tolerant_clustering(overlap_matrix, sorted_indices, threshold, n_skipped)
        Skip-tolerant clustering (falls back to CPU).

    Notes
    -----
    Requires CuPy and a CUDA-capable GPU. Install with:
        pip install cupy-cuda12x

    GPU acceleration is most beneficial for:
    - Large overlap matrix calculations (>1000 peaks)
    - Spectrum generation on dense energy grids
    - Vectorized Gaussian operations

    GPU is typically NOT beneficial for:
    - Sequential algorithms (clustering)
    - Small datasets (<1000 peaks)
    - Algorithms with high CPU-GPU memory transfer overhead
    """

    # CuPy RawKernel for Gaussian overlap calculation
    # This CUDA kernel computes overlaps in parallel on the GPU
    _overlap_kernel_code = r'''
    extern "C" __global__
    void calculate_gaussian_overlap(
        const double* energies,
        const double* widths,
        const double* amplitudes,
        double* overlap_matrix,
        int n
    ) {
        // Thread indices
        int i = blockIdx.x * blockDim.x + threadIdx.x;
        int j = blockIdx.y * blockDim.y + threadIdx.y;

        if (i >= n || j >= n) return;

        // Diagonal elements are 100%
        if (i == j) {
            overlap_matrix[i * n + j] = 100.0;
            return;
        }

        // Only compute upper triangle (exploit symmetry)
        if (i > j) return;

        double center1 = energies[i];
        double center2 = energies[j];
        double sigma1 = widths[i];
        double sigma2 = widths[j];
        double amp1 = amplitudes[i];
        double amp2 = amplitudes[j];

        double min_area = fmin(amp1, amp2);

        if (min_area <= 0.0) {
            overlap_matrix[i * n + j] = 0.0;
            overlap_matrix[j * n + i] = 0.0;
            return;
        }

        // Early termination for distant peaks
        double separation = fabs(center1 - center2);
        double combined_width = (sigma1 + sigma2) * 0.5;

        if (separation > 8.0 * combined_width) {
            overlap_matrix[i * n + j] = 0.0;
            overlap_matrix[j * n + i] = 0.0;
            return;
        }

        // Integration bounds
        double sigma_min = fmin(sigma1, sigma2);
        double x_min = fmin(center1 - 5.0 * sigma1, center2 - 5.0 * sigma2);
        double x_max = fmax(center1 + 5.0 * sigma1, center2 + 5.0 * sigma2);

        // Number of integration points (adaptive)
        int n_points = (int)fmin(10000.0, fmax(1000.0, (x_max - x_min) / (sigma_min * 0.1)));
        double dx = (x_max - x_min) / (n_points - 1);

        // Pre-compute normalization constants
        const double SQRT_2PI = 2.50662827463100050241;
        double norm1 = amp1 / (sigma1 * SQRT_2PI);
        double norm2 = amp2 / (sigma2 * SQRT_2PI);
        double inv_2sigma1_sq = -0.5 / (sigma1 * sigma1);
        double inv_2sigma2_sq = -0.5 / (sigma2 * sigma2);

        // Trapezoidal integration
        double overlap_area = 0.0;
        for (int k = 0; k < n_points; k++) {
            double x = x_min + k * dx;

            // Gaussian 1
            double diff1 = x - center1;
            double g1 = norm1 * exp(inv_2sigma1_sq * diff1 * diff1);

            // Gaussian 2
            double diff2 = x - center2;
            double g2 = norm2 * exp(inv_2sigma2_sq * diff2 * diff2);

            // Overlap (minimum)
            double overlap_val = fmin(g1, g2);

            // Trapezoidal rule
            if (k == 0 || k == n_points - 1) {
                overlap_area += 0.5 * overlap_val;
            } else {
                overlap_area += overlap_val;
            }
        }
        overlap_area *= dx;

        // Percentage of smaller peak
        double percent_overlap = fmax(0.0, fmin(100.0, (overlap_area / min_area) * 100.0));

        // Store symmetric result
        overlap_matrix[i * n + j] = percent_overlap;
        overlap_matrix[j * n + i] = percent_overlap;
    }
    '''

    def __init__(self):
        """Initialize GPU backend and compile CUDA kernels."""
        self._overlap_kernel = cp.RawKernel(
            self._overlap_kernel_code, "calculate_gaussian_overlap"
        )

    def calculate_overlap_matrix(self, energies, widths, amplitudes):
        """
        Calculate overlap matrix using GPU acceleration.

        This method transfers data to GPU, runs a CUDA kernel for parallel
        computation, and returns the result to CPU memory.

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
        >>> backend = GPUBackend()  # doctest: +SKIP
        >>> energies = np.array([285.0, 286.0, 290.0])  # doctest: +SKIP
        >>> widths = np.array([0.5, 0.6, 0.5])  # doctest: +SKIP
        >>> amplitudes = np.array([1.0, 0.8, 0.5])  # doctest: +SKIP
        >>> overlap = backend.calculate_overlap_matrix(energies, widths, amplitudes)  # doctest: +SKIP
        >>> overlap[0, 0]  # doctest: +SKIP
        100.0

        Notes
        -----
        GPU acceleration is beneficial for n > ~1000 peaks. For smaller
        datasets, CPU backend may be faster due to memory transfer overhead.
        """
        n = len(energies)

        # Transfer data to GPU
        energies_gpu = cp.asarray(energies, dtype=cp.float64)
        widths_gpu = cp.asarray(widths, dtype=cp.float64)
        amplitudes_gpu = cp.asarray(amplitudes, dtype=cp.float64)
        overlap_gpu = cp.zeros((n, n), dtype=cp.float64)

        # Launch CUDA kernel
        # Use 16x16 thread blocks for good occupancy
        block_size = (16, 16)
        grid_size = ((n + 15) // 16, (n + 15) // 16)

        self._overlap_kernel(
            grid_size,
            block_size,
            (energies_gpu, widths_gpu, amplitudes_gpu, overlap_gpu, n),
        )

        # Transfer result back to CPU
        return cp.asnumpy(overlap_gpu)

    def sequential_clustering(self, overlap_matrix, sorted_indices, threshold):
        """
        Sequential clustering (CPU fallback).

        Clustering algorithms are inherently sequential and don't benefit from
        GPU parallelization. This method falls back to the CPU implementation.

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

        Notes
        -----
        This method uses the CPU backend because sequential clustering
        algorithms have data dependencies that prevent efficient parallelization.
        """
        # Import here to avoid circular dependency
        from dftlearn._backends.cpu import CPUBackend

        cpu_backend = CPUBackend()
        return cpu_backend.sequential_clustering(
            overlap_matrix, sorted_indices, threshold
        )

    def skip_tolerant_clustering(
        self, overlap_matrix, sorted_indices, threshold, n_skipped=1
    ):
        """
        Skip-tolerant clustering (CPU fallback).

        Clustering algorithms are inherently sequential and don't benefit from
        GPU parallelization. This method falls back to the CPU implementation.

        Parameters
        ----------
        overlap_matrix : ndarray, shape (n, n)
            Overlap percentages between peaks.
        sorted_indices : ndarray, shape (n,)
            Peak indices sorted by energy.
        threshold : float
            Minimum overlap percentage to merge peaks (0-100).
        n_skipped : int, default=1
            Number of peaks below threshold to skip.

        Returns
        -------
        clusters : list of lists
            Each inner list contains indices of peaks in that cluster.

        Notes
        -----
        This method uses the CPU backend because skip-tolerant clustering
        has data dependencies that prevent efficient parallelization.
        """
        # Import here to avoid circular dependency
        from dftlearn._backends.cpu import CPUBackend

        cpu_backend = CPUBackend()
        return cpu_backend.skip_tolerant_clustering(
            overlap_matrix, sorted_indices, threshold, n_skipped
        )
