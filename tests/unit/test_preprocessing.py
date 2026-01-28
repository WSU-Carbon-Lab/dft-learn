"""Unit tests for preprocessing transformers."""

import numpy as np
import pandas as pd
import pytest
from sklearn.pipeline import Pipeline

from dftlearn.preprocessing.transformers import (
    EnergyFilter,
    NormalizedOSTransformer,
    OscillatorStrengthFilter,
)


class TestEnergyFilter:
    """Test EnergyFilter transformer."""

    def test_initialization(self):
        """Test filter initialization."""
        filter_obj = EnergyFilter(energy_threshold=300.0)
        assert filter_obj.energy_threshold == 300.0
        assert filter_obj.verbose is False

    def test_fit(self):
        """Test fitting."""
        df = pd.DataFrame({"E": [285, 290, 325], "OS": [1.0, 0.8, 0.5]})
        filter_obj = EnergyFilter(energy_threshold=300.0)

        filter_obj.fit(df)

        assert hasattr(filter_obj, "n_peaks_before_")
        assert hasattr(filter_obj, "n_peaks_after_")
        assert filter_obj.n_peaks_before_ == 3
        assert filter_obj.n_peaks_after_ == 2

    def test_transform(self):
        """Test transformation."""
        df = pd.DataFrame({"E": [285, 290, 325], "OS": [1.0, 0.8, 0.5]})
        filter_obj = EnergyFilter(energy_threshold=300.0)

        result = filter_obj.fit_transform(df)

        assert len(result) == 2
        assert result["E"].max() <= 300.0
        assert 325 not in result["E"].values

    def test_transform_preserves_other_columns(self):
        """Test that other columns are preserved."""
        df = pd.DataFrame(
            {"E": [285, 290, 325], "OS": [1.0, 0.8, 0.5], "width": [0.5, 0.6, 0.5]}
        )
        filter_obj = EnergyFilter(energy_threshold=300.0)

        result = filter_obj.fit_transform(df)

        assert "OS" in result.columns
        assert "width" in result.columns

    def test_fit_non_dataframe_raises(self):
        """Test that non-DataFrame raises TypeError."""
        filter_obj = EnergyFilter()
        X = np.array([[1, 2], [3, 4]])

        with pytest.raises(TypeError, match="DataFrame"):
            filter_obj.fit(X)

    def test_fit_missing_column_raises(self):
        """Test that missing E column raises ValueError."""
        filter_obj = EnergyFilter()
        df = pd.DataFrame({"OS": [1.0, 0.8]})

        with pytest.raises(ValueError, match="'E' column"):
            filter_obj.fit(df)

    def test_verbose_mode(self, capsys):
        """Test verbose output."""
        df = pd.DataFrame({"E": [285, 290, 325], "OS": [1.0, 0.8, 0.5]})
        filter_obj = EnergyFilter(energy_threshold=300.0, verbose=True)

        filter_obj.fit(df)
        captured = capsys.readouterr()

        assert "Original data: 3 peaks" in captured.out
        assert "After energy filter" in captured.out

    def test_get_params(self):
        """Test get_params for sklearn compatibility."""
        filter_obj = EnergyFilter(energy_threshold=310.0, verbose=True)
        params = filter_obj.get_params()

        assert params["energy_threshold"] == 310.0
        assert params["verbose"] is True

    def test_set_params(self):
        """Test set_params for sklearn compatibility."""
        filter_obj = EnergyFilter()
        filter_obj.set_params(energy_threshold=315.0, verbose=False)

        assert filter_obj.energy_threshold == 315.0
        assert filter_obj.verbose is False


class TestNormalizedOSTransformer:
    """Test NormalizedOSTransformer."""

    def test_initialization(self):
        """Test transformer initialization."""
        transformer = NormalizedOSTransformer()
        assert transformer.os_column == "OS"
        assert transformer.normalized_column == "normalized_OS"

    def test_fit(self):
        """Test fitting."""
        df = pd.DataFrame({"OS": [1.0, 0.5, 0.25]})
        transformer = NormalizedOSTransformer()

        transformer.fit(df)

        assert hasattr(transformer, "max_os_")
        assert transformer.max_os_ == 1.0

    def test_transform(self):
        """Test transformation."""
        df = pd.DataFrame({"OS": [1.0, 0.5, 0.25]})
        transformer = NormalizedOSTransformer()

        result = transformer.fit_transform(df)

        assert "normalized_OS" in result.columns
        np.testing.assert_array_almost_equal(
            result["normalized_OS"].values, [1.0, 0.5, 0.25]
        )

    def test_transform_with_different_max(self):
        """Test normalization with different max value."""
        df = pd.DataFrame({"OS": [2.0, 1.0, 0.5]})
        transformer = NormalizedOSTransformer()

        result = transformer.fit_transform(df)

        assert transformer.max_os_ == 2.0
        np.testing.assert_array_almost_equal(
            result["normalized_OS"].values, [1.0, 0.5, 0.25]
        )

    def test_transform_preserves_original(self):
        """Test that original OS column is preserved."""
        df = pd.DataFrame({"OS": [1.0, 0.5, 0.25]})
        transformer = NormalizedOSTransformer()

        result = transformer.fit_transform(df)

        assert "OS" in result.columns
        np.testing.assert_array_equal(result["OS"].values, df["OS"].values)

    def test_custom_column_names(self):
        """Test with custom column names."""
        df = pd.DataFrame({"amplitude": [2.0, 1.0]})
        transformer = NormalizedOSTransformer(
            os_column="amplitude", normalized_column="norm_amp"
        )

        result = transformer.fit_transform(df)

        assert "norm_amp" in result.columns
        assert result["norm_amp"].max() == 1.0

    def test_fit_non_dataframe_raises(self):
        """Test that non-DataFrame raises TypeError."""
        transformer = NormalizedOSTransformer()
        X = np.array([[1.0], [0.5]])

        with pytest.raises(TypeError, match="DataFrame"):
            transformer.fit(X)

    def test_fit_missing_column_raises(self):
        """Test that missing column raises ValueError."""
        transformer = NormalizedOSTransformer()
        df = pd.DataFrame({"E": [285, 286]})

        with pytest.raises(ValueError, match="'OS' column"):
            transformer.fit(df)

    def test_transform_before_fit_raises(self):
        """Test that transform before fit raises error."""
        transformer = NormalizedOSTransformer()
        df = pd.DataFrame({"OS": [1.0, 0.5]})

        with pytest.raises(Exception):  # sklearn raises NotFittedError
            transformer.transform(df)


class TestOscillatorStrengthFilter:
    """Test OscillatorStrengthFilter."""

    def test_initialization(self):
        """Test filter initialization."""
        filter_obj = OscillatorStrengthFilter(os_threshold_percent=5.0)
        assert filter_obj.os_threshold_percent == 5.0

    def test_fit(self):
        """Test fitting."""
        df = pd.DataFrame({"OS": [1.0, 0.5, 0.01]})
        filter_obj = OscillatorStrengthFilter(os_threshold_percent=10.0)

        filter_obj.fit(df)

        assert hasattr(filter_obj, "max_os_")
        assert hasattr(filter_obj, "os_threshold_")
        assert filter_obj.max_os_ == 1.0
        assert filter_obj.os_threshold_ == 0.1  # 10% of 1.0

    def test_transform(self):
        """Test transformation."""
        df = pd.DataFrame({"OS": [1.0, 0.5, 0.01]})
        filter_obj = OscillatorStrengthFilter(os_threshold_percent=10.0)

        result = filter_obj.fit_transform(df)

        assert len(result) == 2  # Only first two peaks >= 10%
        assert 0.01 not in result["OS"].values

    def test_transform_2_percent_threshold(self):
        """Test with 2% threshold (default)."""
        df = pd.DataFrame({"OS": [1.0, 0.5, 0.03, 0.01]})
        filter_obj = OscillatorStrengthFilter(os_threshold_percent=2.0)

        result = filter_obj.fit_transform(df)

        assert len(result) == 3  # First three >= 2%
        assert 0.01 not in result["OS"].values

    def test_custom_column_name(self):
        """Test with custom column name."""
        df = pd.DataFrame({"amplitude": [1.0, 0.5, 0.01]})
        filter_obj = OscillatorStrengthFilter(
            os_threshold_percent=10.0, os_column="amplitude"
        )

        result = filter_obj.fit_transform(df)

        assert len(result) == 2

    def test_fit_stores_statistics(self):
        """Test that fit stores filtering statistics."""
        df = pd.DataFrame({"OS": [1.0, 0.5, 0.01, 0.005]})
        filter_obj = OscillatorStrengthFilter(os_threshold_percent=2.0)

        filter_obj.fit(df)

        assert filter_obj.n_peaks_before_ == 4
        assert filter_obj.n_peaks_after_ == 2  # 2 peaks >= 2% (1.0 and 0.5)

    def test_verbose_mode(self, capsys):
        """Test verbose output."""
        df = pd.DataFrame({"OS": [1.0, 0.5, 0.01]})
        filter_obj = OscillatorStrengthFilter(
            os_threshold_percent=10.0, verbose=True
        )

        filter_obj.fit(df)
        captured = capsys.readouterr()

        assert "Maximum OS" in captured.out
        assert "Threshold" in captured.out
        assert "10.0%" in captured.out


class TestPreprocessingPipeline:
    """Test combining preprocessing transformers in a pipeline."""

    def test_energy_and_os_pipeline(self):
        """Test pipeline with energy and OS filtering."""
        df = pd.DataFrame(
            {
                "E": [285, 290, 295, 300, 330],
                "OS": [1.0, 0.8, 0.5, 0.03, 0.02],
                "width": [0.5, 0.5, 0.6, 0.5, 0.5],
            }
        )

        pipeline = Pipeline(
            [
                ("energy_filter", EnergyFilter(energy_threshold=320.0)),
                ("os_filter", OscillatorStrengthFilter(os_threshold_percent=5.0)),
            ]
        )

        result = pipeline.fit_transform(df)

        # Should remove E=330 (energy filter) and E=300 (OS filter)
        assert len(result) == 3
        assert result["E"].max() <= 320
        assert result["OS"].min() >= 0.05

    def test_full_preprocessing_pipeline(self):
        """Test complete preprocessing pipeline."""
        df = pd.DataFrame(
            {"E": [285, 290, 295, 330], "OS": [1.0, 0.8, 0.5, 0.02], "width": [0.5] * 4}
        )

        pipeline = Pipeline(
            [
                ("energy_filter", EnergyFilter(energy_threshold=320.0)),
                ("normalize", NormalizedOSTransformer()),
                ("os_filter", OscillatorStrengthFilter(os_threshold_percent=10.0)),
            ]
        )

        result = pipeline.fit_transform(df)

        # Should have normalized_OS column
        assert "normalized_OS" in result.columns

        # Should filter by both energy and OS
        # Energy filter removes 330, OS filter keeps all >= 10% (1.0, 0.8, 0.5)
        assert len(result) == 3  # 285, 290, 295 remain
        assert result["normalized_OS"].max() == 1.0

    def test_pipeline_with_get_params(self):
        """Test that pipeline parameters are accessible."""
        pipeline = Pipeline(
            [
                ("energy_filter", EnergyFilter(energy_threshold=310.0)),
                ("os_filter", OscillatorStrengthFilter(os_threshold_percent=3.0)),
            ]
        )

        params = pipeline.get_params()

        assert params["energy_filter__energy_threshold"] == 310.0
        assert params["os_filter__os_threshold_percent"] == 3.0

    def test_pipeline_with_set_params(self):
        """Test setting pipeline parameters."""
        pipeline = Pipeline(
            [
                ("energy_filter", EnergyFilter()),
                ("os_filter", OscillatorStrengthFilter()),
            ]
        )

        pipeline.set_params(
            energy_filter__energy_threshold=315.0,
            os_filter__os_threshold_percent=5.0,
        )

        assert pipeline.named_steps["energy_filter"].energy_threshold == 315.0
        assert pipeline.named_steps["os_filter"].os_threshold_percent == 5.0
