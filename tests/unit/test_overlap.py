"""Unit tests for overlap calculations."""

import numpy as np
import pandas as pd
import pytest

from dftlearn.analysis.overlap import OverlapCalculator, OverlapResult


class TestOverlapResult:
    """Test OverlapResult dataclass."""

    def test_creation(self):
        """Test creating OverlapResult."""
        overlap_matrix = np.array([[100, 50], [50, 100]])
        energies = np.array([285.0, 286.0])
        widths = np.array([0.5, 0.5])
        amplitudes = np.array([1.0, 0.8])

        result = OverlapResult(
            overlap_matrix=overlap_matrix,
            energies=energies,
            widths=widths,
            amplitudes=amplitudes,
        )

        assert result.overlap_matrix.shape == (2, 2)
        assert len(result.energies) == 2
        np.testing.assert_array_equal(result.energies, energies)

    def test_validation_shape_mismatch(self):
        """Test that shape mismatch raises ValueError."""
        with pytest.raises(ValueError, match="shape.*does not match"):
            OverlapResult(
                overlap_matrix=np.array([[100, 50], [50, 100]]),
                energies=np.array([285.0]),  # Wrong length
                widths=np.array([0.5]),
                amplitudes=np.array([1.0]),
            )

    def test_validation_widths_mismatch(self):
        """Test that widths length mismatch raises ValueError."""
        with pytest.raises(ValueError, match="widths.*does not match"):
            OverlapResult(
                overlap_matrix=np.array([[100, 50], [50, 100]]),
                energies=np.array([285.0, 286.0]),
                widths=np.array([0.5]),  # Wrong length
                amplitudes=np.array([1.0, 0.8]),
            )


class TestOverlapCalculator:
    """Test OverlapCalculator class."""

    def test_initialization(self):
        """Test calculator initialization."""
        calc = OverlapCalculator(backend="cpu")
        assert calc.backend == "cpu"

    def test_initialization_auto(self):
        """Test auto backend selection."""
        calc = OverlapCalculator(backend="auto")
        assert calc.backend == "auto"
        assert calc._backend_instance is not None

    def test_calculate_diagonal(self):
        """Test that diagonal elements are 100%."""
        calc = OverlapCalculator(backend="cpu")
        energies = np.array([285.0, 286.0, 287.0])
        widths = np.array([0.5, 0.5, 0.5])
        amplitudes = np.array([1.0, 1.0, 1.0])

        result = calc.calculate(energies, widths, amplitudes)

        np.testing.assert_array_almost_equal(
            np.diag(result.overlap_matrix), [100.0, 100.0, 100.0]
        )

    def test_calculate_symmetric(self):
        """Test that overlap matrix is symmetric."""
        calc = OverlapCalculator(backend="cpu")
        energies = np.array([285.0, 286.0, 290.0])
        widths = np.array([0.5, 0.6, 0.5])
        amplitudes = np.array([1.0, 0.8, 0.5])

        result = calc.calculate(energies, widths, amplitudes)

        np.testing.assert_array_almost_equal(
            result.overlap_matrix, result.overlap_matrix.T
        )

    def test_calculate_distant_peaks(self):
        """Test that distant peaks have low overlap."""
        calc = OverlapCalculator(backend="cpu")
        energies = np.array([285.0, 295.0])  # 10 eV apart
        widths = np.array([0.5, 0.5])
        amplitudes = np.array([1.0, 1.0])

        result = calc.calculate(energies, widths, amplitudes)

        # Peaks 10 eV apart with 0.5 eV width should have ~0% overlap
        assert result.overlap_matrix[0, 1] < 0.1

    def test_calculate_close_peaks(self):
        """Test that close peaks have high overlap."""
        calc = OverlapCalculator(backend="cpu")
        energies = np.array([285.0, 285.3])  # 0.3 eV apart
        widths = np.array([0.5, 0.5])  # 1 sigma = 0.5 eV
        amplitudes = np.array([1.0, 1.0])

        result = calc.calculate(energies, widths, amplitudes)

        # Peaks 0.3 eV apart with 0.5 eV width should have significant overlap
        assert result.overlap_matrix[0, 1] > 50.0

    def test_calculate_list_inputs(self):
        """Test that list inputs work."""
        calc = OverlapCalculator(backend="cpu")
        energies = [285.0, 286.0, 287.0]
        widths = [0.5, 0.5, 0.5]
        amplitudes = [1.0, 0.8, 0.6]

        result = calc.calculate(energies, widths, amplitudes)

        assert result.overlap_matrix.shape == (3, 3)
        assert isinstance(result.energies, np.ndarray)

    def test_calculate_length_mismatch(self):
        """Test that length mismatch raises ValueError."""
        calc = OverlapCalculator(backend="cpu")

        with pytest.raises(ValueError, match="same length"):
            calc.calculate(
                energies=[285.0, 286.0],
                widths=[0.5],  # Wrong length
                amplitudes=[1.0, 0.8],
            )

    def test_calculate_empty_input(self):
        """Test that empty input raises ValueError."""
        calc = OverlapCalculator(backend="cpu")

        with pytest.raises(ValueError, match="empty"):
            calc.calculate(energies=[], widths=[], amplitudes=[])

    def test_calculate_from_dataframe(self, sample_peaks):
        """Test calculating from DataFrame."""
        calc = OverlapCalculator(backend="cpu")

        result = calc.calculate_from_dataframe(sample_peaks)

        assert result.overlap_matrix.shape == (len(sample_peaks), len(sample_peaks))
        assert len(result.energies) == len(sample_peaks)

    def test_calculate_from_dataframe_custom_columns(self):
        """Test DataFrame with custom column names."""
        calc = OverlapCalculator(backend="cpu")

        df = pd.DataFrame(
            {
                "energy": [285.0, 286.0],
                "sigma": [0.5, 0.5],
                "amplitude": [1.0, 0.8],
            }
        )

        result = calc.calculate_from_dataframe(
            df, energy_col="energy", width_col="sigma", amplitude_col="amplitude"
        )

        assert result.overlap_matrix.shape == (2, 2)

    def test_calculate_from_dataframe_missing_column(self):
        """Test that missing column raises KeyError."""
        calc = OverlapCalculator(backend="cpu")

        df = pd.DataFrame({"E": [285.0, 286.0], "width": [0.5, 0.5]})  # Missing OS

        with pytest.raises(KeyError, match="not found"):
            calc.calculate_from_dataframe(df)

    def test_repr(self):
        """Test string representation."""
        calc = OverlapCalculator(backend="cpu")
        repr_str = repr(calc)

        assert "OverlapCalculator" in repr_str
        assert "cpu" in repr_str


@pytest.mark.gpu
class TestOverlapCalculatorGPU:
    """Test GPU backend for overlap calculation."""

    def test_gpu_backend(self):
        """Test GPU backend initialization."""
        pytest.importorskip("cupy")
        calc = OverlapCalculator(backend="gpu")
        assert calc.backend == "gpu"

    def test_cpu_gpu_equivalence(self, sample_peaks):
        """Test that CPU and GPU give same results."""
        pytest.importorskip("cupy")

        calc_cpu = OverlapCalculator(backend="cpu")
        calc_gpu = OverlapCalculator(backend="gpu")

        result_cpu = calc_cpu.calculate_from_dataframe(sample_peaks)
        result_gpu = calc_gpu.calculate_from_dataframe(sample_peaks)

        np.testing.assert_allclose(
            result_cpu.overlap_matrix, result_gpu.overlap_matrix, rtol=1e-5
        )
