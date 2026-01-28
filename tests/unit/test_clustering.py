"""Unit tests for clustering algorithms."""

import numpy as np
import pandas as pd
import pytest
from sklearn.utils.estimator_checks import parametrize_with_checks

from dftlearn.analysis.clustering import OverlapClusterer


class TestOverlapClusterer:
    """Test OverlapClusterer estimator."""

    def test_initialization(self):
        """Test clusterer initialization."""
        clusterer = OverlapClusterer(
            overlap_threshold=60.0, method="sequential", backend="cpu"
        )
        assert clusterer.overlap_threshold == 60.0
        assert clusterer.method == "sequential"
        assert clusterer.backend == "cpu"

    def test_fit_basic(self, sample_peaks):
        """Test basic fitting."""
        clusterer = OverlapClusterer(overlap_threshold=50.0)
        clusterer.fit(sample_peaks)

        assert hasattr(clusterer, "labels_")
        assert hasattr(clusterer, "n_clusters_")
        assert hasattr(clusterer, "overlap_matrix_")
        assert len(clusterer.labels_) == len(sample_peaks)

    def test_fit_overlapping_peaks(self, sample_overlapping_peaks):
        """Test clustering with known overlap structure."""
        # Peaks 0-2 should cluster, peak 3 alone, peaks 4-5 cluster
        clusterer = OverlapClusterer(overlap_threshold=50.0, method="sequential")
        clusterer.fit(sample_overlapping_peaks)

        # Should identify at least 2 clusters (maybe 3)
        assert clusterer.n_clusters_ >= 2

        # Peaks 0 and 1 should be in same cluster (0.5 eV apart)
        assert clusterer.labels_[0] == clusterer.labels_[1]

        # Peak 3 should be separate (5 eV from others)
        assert clusterer.labels_[3] != clusterer.labels_[0]

    def test_fit_predict(self, sample_peaks):
        """Test fit_predict method."""
        clusterer = OverlapClusterer(overlap_threshold=50.0)
        labels = clusterer.fit_predict(sample_peaks)

        assert len(labels) == len(sample_peaks)
        assert hasattr(clusterer, "n_clusters_")

    def test_predict_after_fit(self, sample_peaks):
        """Test predict returns fitted labels."""
        clusterer = OverlapClusterer(overlap_threshold=50.0)
        clusterer.fit(sample_peaks)
        labels = clusterer.predict(sample_peaks)

        np.testing.assert_array_equal(labels, clusterer.labels_)

    def test_predict_before_fit_raises(self, sample_peaks):
        """Test that predict before fit raises error."""
        clusterer = OverlapClusterer(overlap_threshold=50.0)

        with pytest.raises(Exception):  # sklearn raises NotFittedError
            clusterer.predict(sample_peaks)

    def test_sequential_method(self, sample_peaks):
        """Test sequential clustering method."""
        clusterer = OverlapClusterer(
            overlap_threshold=50.0, method="sequential", backend="cpu"
        )
        clusterer.fit(sample_peaks)

        assert clusterer.n_clusters_ >= 1
        assert len(clusterer.cluster_info_) == clusterer.n_clusters_

    def test_skip_tolerant_method(self, sample_peaks):
        """Test skip-tolerant clustering method."""
        clusterer = OverlapClusterer(
            overlap_threshold=50.0, method="skip_tolerant", n_skipped=2, backend="cpu"
        )
        clusterer.fit(sample_peaks)

        assert clusterer.n_clusters_ >= 1
        assert len(clusterer.cluster_info_) == clusterer.n_clusters_

    def test_invalid_method_raises(self, sample_peaks):
        """Test that invalid method raises ValueError."""
        clusterer = OverlapClusterer(
            overlap_threshold=50.0, method="invalid_method"
        )

        with pytest.raises(ValueError, match="Unknown clustering method"):
            clusterer.fit(sample_peaks)

    def test_fit_non_dataframe_raises(self):
        """Test that non-DataFrame input raises TypeError."""
        clusterer = OverlapClusterer()
        X = np.array([[1, 2, 3], [4, 5, 6]])

        with pytest.raises(TypeError, match="DataFrame"):
            clusterer.fit(X)

    def test_fit_missing_columns_raises(self):
        """Test that missing required columns raises ValueError."""
        clusterer = OverlapClusterer()
        df = pd.DataFrame({"E": [285.0, 286.0]})  # Missing width and OS

        with pytest.raises(ValueError, match="Missing required columns"):
            clusterer.fit(df)

    def test_fit_empty_dataframe_raises(self):
        """Test that empty DataFrame raises ValueError."""
        clusterer = OverlapClusterer()
        df = pd.DataFrame({"E": [], "width": [], "OS": []})

        with pytest.raises(ValueError, match="empty"):
            clusterer.fit(df)

    def test_get_cluster_peaks(self, sample_overlapping_peaks):
        """Test getting peaks grouped by cluster."""
        clusterer = OverlapClusterer(overlap_threshold=50.0)
        clusterer.fit(sample_overlapping_peaks)

        cluster_peaks = clusterer.get_cluster_peaks(sample_overlapping_peaks)

        assert len(cluster_peaks) == clusterer.n_clusters_
        assert all(isinstance(df, pd.DataFrame) for df in cluster_peaks)

        # Total peaks across clusters should equal original
        total_peaks = sum(len(df) for df in cluster_peaks)
        assert total_peaks == len(sample_overlapping_peaks)

    def test_different_thresholds(self, sample_peaks):
        """Test that different thresholds give different clustering."""
        clusterer_low = OverlapClusterer(overlap_threshold=30.0)
        clusterer_high = OverlapClusterer(overlap_threshold=70.0)

        clusterer_low.fit(sample_peaks)
        clusterer_high.fit(sample_peaks)

        # Higher threshold should generally give more clusters
        assert clusterer_high.n_clusters_ >= clusterer_low.n_clusters_

    def test_cluster_info_attribute(self, sample_peaks):
        """Test cluster_info_ attribute structure."""
        clusterer = OverlapClusterer(overlap_threshold=50.0)
        clusterer.fit(sample_peaks)

        assert hasattr(clusterer, "cluster_info_")
        assert isinstance(clusterer.cluster_info_, list)
        assert len(clusterer.cluster_info_) == clusterer.n_clusters_

        # Each cluster should be a list of indices
        for cluster in clusterer.cluster_info_:
            assert isinstance(cluster, list)
            assert all(isinstance(idx, (int, np.integer)) for idx in cluster)

    def test_overlap_matrix_stored(self, sample_peaks):
        """Test that overlap matrix is stored."""
        clusterer = OverlapClusterer(overlap_threshold=50.0)
        clusterer.fit(sample_peaks)

        assert hasattr(clusterer, "overlap_matrix_")
        assert clusterer.overlap_matrix_.shape == (len(sample_peaks), len(sample_peaks))

        # Should be symmetric
        np.testing.assert_array_almost_equal(
            clusterer.overlap_matrix_, clusterer.overlap_matrix_.T
        )

        # Diagonal should be 100
        np.testing.assert_array_almost_equal(
            np.diag(clusterer.overlap_matrix_), np.full(len(sample_peaks), 100.0)
        )

    def test_repr(self):
        """Test string representation."""
        clusterer = OverlapClusterer(
            overlap_threshold=60.0, method="skip_tolerant", n_skipped=2
        )
        repr_str = repr(clusterer)

        assert "OverlapClusterer" in repr_str
        assert "60.0" in repr_str
        assert "skip_tolerant" in repr_str

    def test_get_params(self):
        """Test get_params for sklearn compatibility."""
        clusterer = OverlapClusterer(
            overlap_threshold=55.0, method="sequential", backend="cpu"
        )
        params = clusterer.get_params()

        assert params["overlap_threshold"] == 55.0
        assert params["method"] == "sequential"
        assert params["backend"] == "cpu"

    def test_set_params(self):
        """Test set_params for sklearn compatibility."""
        clusterer = OverlapClusterer()
        clusterer.set_params(overlap_threshold=65.0, method="skip_tolerant")

        assert clusterer.overlap_threshold == 65.0
        assert clusterer.method == "skip_tolerant"


@pytest.mark.gpu
class TestOverlapClustererGPU:
    """Test GPU backend for clustering."""

    def test_gpu_backend_fit(self, sample_peaks):
        """Test fitting with GPU backend."""
        pytest.importorskip("cupy")

        clusterer = OverlapClusterer(overlap_threshold=50.0, backend="gpu")
        clusterer.fit(sample_peaks)

        assert hasattr(clusterer, "labels_")
        assert hasattr(clusterer, "n_clusters_")

    def test_cpu_gpu_equivalence(self, sample_peaks):
        """Test that CPU and GPU give same clustering."""
        pytest.importorskip("cupy")

        clusterer_cpu = OverlapClusterer(overlap_threshold=50.0, backend="cpu")
        clusterer_gpu = OverlapClusterer(overlap_threshold=50.0, backend="gpu")

        clusterer_cpu.fit(sample_peaks)
        clusterer_gpu.fit(sample_peaks)

        # Should produce same number of clusters
        assert clusterer_cpu.n_clusters_ == clusterer_gpu.n_clusters_

        # Should produce same labels
        np.testing.assert_array_equal(clusterer_cpu.labels_, clusterer_gpu.labels_)
