"""Unit tests for backend dispatch and computation."""

import numpy as np
import pytest

from dftlearn._backends import get_backend
from dftlearn._backends.cpu import CPUBackend


class TestBackendDispatch:
    """Test backend selection and dispatch."""

    def test_get_cpu_backend(self):
        """Test CPU backend selection."""
        backend = get_backend("cpu")
        assert isinstance(backend, CPUBackend)

    def test_get_auto_backend(self):
        """Test auto backend selection."""
        backend = get_backend("auto")
        # Should return either CPU or GPU backend
        assert backend is not None
        assert hasattr(backend, "calculate_overlap_matrix")

    def test_invalid_backend(self):
        """Test invalid backend raises ValueError."""
        with pytest.raises(ValueError, match="Unknown backend"):
            get_backend("invalid")


class TestCPUBackend:
    """Test CPU backend computation."""

    def test_overlap_matrix_diagonal(self, cpu_backend):
        """Test that diagonal elements are 100%."""
        energies = np.array([285.0, 286.0, 287.0])
        widths = np.array([0.5, 0.5, 0.5])
        amplitudes = np.array([1.0, 1.0, 1.0])

        overlap = cpu_backend.calculate_overlap_matrix(energies, widths, amplitudes)

        assert overlap.shape == (3, 3)
        np.testing.assert_array_almost_equal(np.diag(overlap), [100.0, 100.0, 100.0])

    def test_overlap_matrix_symmetric(self, cpu_backend):
        """Test that overlap matrix is symmetric."""
        energies = np.array([285.0, 286.0, 290.0])
        widths = np.array([0.5, 0.6, 0.5])
        amplitudes = np.array([1.0, 0.8, 0.5])

        overlap = cpu_backend.calculate_overlap_matrix(energies, widths, amplitudes)

        np.testing.assert_array_almost_equal(overlap, overlap.T)

    def test_sequential_clustering(self, cpu_backend, simple_overlap_matrix):
        """Test sequential clustering algorithm."""
        sorted_indices = np.array([0, 1, 2, 3, 4])
        threshold = 50.0

        clusters = cpu_backend.sequential_clustering(
            simple_overlap_matrix, sorted_indices, threshold
        )

        # Expected: peaks 0-1 cluster, peak 2 alone, peaks 3-4 cluster
        assert len(clusters) >= 2
        assert 0 in clusters[0] and 1 in clusters[0]

    def test_skip_tolerant_clustering(self, cpu_backend, simple_overlap_matrix):
        """Test skip-tolerant clustering algorithm."""
        sorted_indices = np.array([0, 1, 2, 3, 4])
        threshold = 50.0

        clusters = cpu_backend.skip_tolerant_clustering(
            simple_overlap_matrix, sorted_indices, threshold, n_skipped=1
        )

        # Should produce clusters
        assert len(clusters) >= 1
        assert isinstance(clusters, list)
