"""Shared pytest fixtures for dftlearn tests."""

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def sample_peaks():
    """
    Generate synthetic peak data for testing.

    Returns
    -------
    DataFrame
        Synthetic peak data with columns:
        - E: Energy (eV)
        - width: Peak width/standard deviation (eV)
        - OS: Oscillator strength
        - normalized_OS: Normalized oscillator strength
        - theta: Orientation angle (degrees)
        - xx, yy, zz: Diagonal tensor components
        - xy, xz, yz: Off-diagonal tensor components

    Examples
    --------
    >>> def test_clustering(sample_peaks):
    ...     assert len(sample_peaks) == 10
    ...     assert 'E' in sample_peaks.columns
    """
    np.random.seed(42)  # Reproducible test data

    n_peaks = 10
    return pd.DataFrame(
        {
            "E": np.linspace(280, 290, n_peaks),
            "width": np.random.uniform(0.5, 1.5, n_peaks),
            "OS": np.random.uniform(0.1, 1.0, n_peaks),
            "normalized_OS": np.random.uniform(0.0, 1.0, n_peaks),
            "theta": np.random.uniform(0, 90, n_peaks),
            "xx": np.random.uniform(0, 1, n_peaks),
            "yy": np.random.uniform(0, 1, n_peaks),
            "zz": np.random.uniform(0, 1, n_peaks),
            "xy": np.zeros(n_peaks),
            "xz": np.zeros(n_peaks),
            "yz": np.zeros(n_peaks),
        }
    )


@pytest.fixture
def sample_overlapping_peaks():
    """
    Generate synthetic peaks with known overlaps for testing clustering.

    Returns
    -------
    DataFrame
        Peak data with:
        - Peaks 0-2 close together (should cluster)
        - Peak 3 isolated (separate cluster)
        - Peaks 4-5 close together (should cluster)

    Examples
    --------
    >>> def test_clustering_known_structure(sample_overlapping_peaks):
    ...     # Test that clusterer identifies the 3 expected clusters
    ...     pass
    """
    return pd.DataFrame(
        {
            "E": np.array([285.0, 285.5, 286.0, 290.0, 292.0, 292.5]),
            "width": np.array([0.5, 0.5, 0.5, 0.4, 0.6, 0.6]),
            "OS": np.array([1.0, 0.8, 0.6, 0.5, 0.7, 0.9]),
            "normalized_OS": np.array([1.0, 0.8, 0.6, 0.5, 0.7, 0.9]),
            "theta": np.array([45, 50, 48, 60, 55, 52]),
            "xx": np.array([0.3, 0.3, 0.3, 0.4, 0.2, 0.2]),
            "yy": np.array([0.4, 0.4, 0.4, 0.3, 0.5, 0.5]),
            "zz": np.array([0.3, 0.3, 0.3, 0.3, 0.3, 0.3]),
            "xy": np.zeros(6),
            "xz": np.zeros(6),
            "yz": np.zeros(6),
        }
    )


@pytest.fixture
def sample_experimental_spectra():
    """
    Generate synthetic experimental NEXAFS spectra at multiple angles.

    Returns
    -------
    dict
        Dictionary mapping angle (degrees) to (energy, intensity) tuples:
        - angles: [20, 45, 70, 90]
        - energy: 500-point grid from 280 to 295 eV
        - intensity: Noisy Gaussian-like spectrum

    Examples
    --------
    >>> def test_fitting(sample_experimental_spectra):
    ...     spectrum_20deg = sample_experimental_spectra[20]
    ...     energy, intensity = spectrum_20deg
    ...     assert len(energy) == len(intensity)
    """
    np.random.seed(42)

    angles = [20, 45, 70, 90]
    energy = np.linspace(280, 295, 500)
    spectra = {}

    for angle in angles:
        # Generate synthetic spectrum with angular dependence
        # Simulate cos²θ dichroism
        angle_factor = 0.5 + 0.5 * np.cos(np.radians(angle)) ** 2

        # Main peak + noise
        intensity = (
            angle_factor * np.exp(-((energy - 285) / 2) ** 2)
            + 0.3 * np.exp(-((energy - 287.5) / 1.5) ** 2)
            + np.random.normal(0, 0.01, len(energy))
        )

        # Ensure non-negative
        intensity = np.maximum(intensity, 0)

        spectra[angle] = (energy, intensity)

    return spectra


@pytest.fixture(params=["cpu", "gpu"])
def backend(request):
    """
    Parametrize tests across CPU and GPU backends.

    This fixture allows the same test to run with both backends,
    ensuring numerical equivalence and catching backend-specific bugs.

    Parameters
    ----------
    request : FixtureRequest
        Pytest fixture request with param attribute.

    Returns
    -------
    str
        Backend name: 'cpu' or 'gpu'

    Examples
    --------
    >>> def test_overlap_matrix(backend):
    ...     from dftlearn._backends import get_backend
    ...     b = get_backend(backend)
    ...     # Test runs twice: once with CPU, once with GPU
    """
    if request.param == "gpu":
        pytest.importorskip("cupy")  # Skip GPU tests if CuPy not available
    return request.param


@pytest.fixture
def cpu_backend():
    """
    Get CPU backend instance.

    Returns
    -------
    CPUBackend
        CPU backend for testing.

    Examples
    --------
    >>> def test_cpu_specific(cpu_backend):
    ...     result = cpu_backend.calculate_overlap_matrix(...)
    """
    from dftlearn._backends.cpu import CPUBackend

    return CPUBackend()


@pytest.fixture
def gpu_backend():
    """
    Get GPU backend instance (requires CuPy).

    Returns
    -------
    GPUBackend
        GPU backend for testing.

    Examples
    --------
    >>> @pytest.mark.gpu
    ... def test_gpu_specific(gpu_backend):
    ...     result = gpu_backend.calculate_overlap_matrix(...)
    """
    pytest.importorskip("cupy")
    from dftlearn._backends.gpu import GPUBackend

    return GPUBackend()


@pytest.fixture
def simple_overlap_matrix():
    """
    Create a simple overlap matrix for clustering tests.

    Returns
    -------
    ndarray, shape (5, 5)
        Overlap matrix with known structure:
        - Peaks 0-1 overlap highly (should cluster)
        - Peak 2 isolated
        - Peaks 3-4 overlap highly (should cluster)

    Examples
    --------
    >>> def test_clustering_simple(simple_overlap_matrix):
    ...     # Test clustering on known overlap structure
    ...     pass
    """
    return np.array(
        [
            [100.0, 75.0, 10.0, 5.0, 2.0],  # Peak 0
            [75.0, 100.0, 15.0, 8.0, 3.0],  # Peak 1 (overlaps with 0)
            [10.0, 15.0, 100.0, 12.0, 8.0],  # Peak 2 (isolated)
            [5.0, 8.0, 12.0, 100.0, 80.0],  # Peak 3
            [2.0, 3.0, 8.0, 80.0, 100.0],  # Peak 4 (overlaps with 3)
        ]
    )


@pytest.fixture
def energy_grid():
    """
    Standard energy grid for NEXAFS spectra.

    Returns
    -------
    ndarray
        Energy points from 280 to 320 eV with 0.1 eV spacing.

    Examples
    --------
    >>> def test_spectrum_generation(energy_grid):
    ...     assert len(energy_grid) == 401
    ...     assert energy_grid[0] == 280.0
    """
    return np.linspace(280, 320, 401)


@pytest.fixture(scope="session")
def test_data_dir(tmp_path_factory):
    """
    Create temporary directory for test data files.

    Parameters
    ----------
    tmp_path_factory : TempPathFactory
        Pytest temporary path factory.

    Returns
    -------
    Path
        Path to temporary test data directory.

    Examples
    --------
    >>> def test_file_io(test_data_dir):
    ...     test_file = test_data_dir / "test.csv"
    ...     # Write and read test file
    """
    return tmp_path_factory.mktemp("test_data")
