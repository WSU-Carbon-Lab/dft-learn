"""Overlap calculation between Gaussian spectral peaks."""

from dataclasses import dataclass

import numpy as np
import pandas as pd

from dftlearn._backends import get_backend


@dataclass
class OverlapResult:
    """
    Container for overlap calculation results.

    Attributes
    ----------
    overlap_matrix : ndarray, shape (n, n)
        Symmetric matrix of percentage overlaps between peaks.
    energies : ndarray, shape (n,)
        Peak center energies (eV).
    widths : ndarray, shape (n,)
        Peak widths/standard deviations (eV).
    amplitudes : ndarray, shape (n,)
        Peak amplitudes/oscillator strengths.

    Examples
    --------
    >>> result = OverlapResult(
    ...     overlap_matrix=np.array([[100, 50], [50, 100]]),
    ...     energies=np.array([285.0, 286.0]),
    ...     widths=np.array([0.5, 0.5]),
    ...     amplitudes=np.array([1.0, 0.8])
    ... )
    >>> result.overlap_matrix[0, 1]
    50.0
    """

    overlap_matrix: np.ndarray
    energies: np.ndarray
    widths: np.ndarray
    amplitudes: np.ndarray

    def __post_init__(self):
        """Validate data after initialization."""
        n = len(self.energies)
        if self.overlap_matrix.shape != (n, n):
            raise ValueError(
                f"Overlap matrix shape {self.overlap_matrix.shape} "
                f"does not match number of peaks {n}"
            )
        if len(self.widths) != n:
            raise ValueError(
                f"Number of widths {len(self.widths)} "
                f"does not match number of peaks {n}"
            )
        if len(self.amplitudes) != n:
            raise ValueError(
                f"Number of amplitudes {len(self.amplitudes)} "
                f"does not match number of peaks {n}"
            )


class OverlapCalculator:
    """
    Calculate overlap between Gaussian spectral peaks.

    This class provides a high-level interface for calculating overlap matrices
    between Gaussian peaks using either CPU (Numba) or GPU (CuPy) backends.

    The overlap between two Gaussian peaks is calculated as the percentage of
    the smaller peak's area that overlaps with the larger peak, using the
    minimum envelope method.

    Parameters
    ----------
    backend : {'auto', 'cpu', 'gpu'}, default='auto'
        Computation backend:
        - 'auto': Automatically select GPU if available, otherwise CPU
        - 'cpu': Force CPU backend (NumPy/Numba)
        - 'gpu': Force GPU backend (CuPy) - requires CUDA

    Attributes
    ----------
    backend : str
        Selected backend name.

    Methods
    -------
    calculate(energies, widths, amplitudes)
        Calculate overlap matrix from peak parameters.
    calculate_from_dataframe(df, energy_col='E', width_col='width', amplitude_col='OS')
        Calculate overlap matrix from DataFrame.

    Examples
    --------
    >>> import numpy as np
    >>> calc = OverlapCalculator(backend='cpu')
    >>> energies = np.array([285.0, 286.0, 290.0])
    >>> widths = np.array([0.5, 0.6, 0.5])
    >>> amplitudes = np.array([1.0, 0.8, 0.5])
    >>> result = calc.calculate(energies, widths, amplitudes)
    >>> result.overlap_matrix[0, 0]
    100.0
    >>> result.overlap_matrix[0, 2] < 5.0  # Distant peaks have low overlap
    True

    Notes
    -----
    The overlap calculation uses the minimum envelope method, where the overlap
    area is the integral of min(g1(x), g2(x)) over all x, divided by the area
    of the smaller peak.

    For GPU acceleration (backend='gpu'), CuPy must be installed. GPU is most
    beneficial for large datasets (>1000 peaks).
    """

    def __init__(self, backend: str = "auto"):
        """
        Initialize overlap calculator with specified backend.

        Parameters
        ----------
        backend : {'auto', 'cpu', 'gpu'}, default='auto'
            Computation backend to use.
        """
        self.backend = backend
        self._backend_instance = get_backend(backend)

    def calculate(
        self,
        energies: np.ndarray | list,
        widths: np.ndarray | list,
        amplitudes: np.ndarray | list,
    ) -> OverlapResult:
        """
        Calculate overlap matrix from peak parameters.

        Parameters
        ----------
        energies : array-like, shape (n,)
            Peak center energies (eV).
        widths : array-like, shape (n,)
            Peak widths/standard deviations (eV).
        amplitudes : array-like, shape (n,)
            Peak amplitudes/oscillator strengths.

        Returns
        -------
        OverlapResult
            Container with overlap matrix and peak parameters.

        Raises
        ------
        ValueError
            If input arrays have different lengths.

        Examples
        --------
        >>> calc = OverlapCalculator()
        >>> energies = [285.0, 286.0, 287.0]
        >>> widths = [0.5, 0.5, 0.5]
        >>> amplitudes = [1.0, 0.8, 0.6]
        >>> result = calc.calculate(energies, widths, amplitudes)
        >>> result.overlap_matrix.shape
        (3, 3)
        """
        # Convert to numpy arrays
        energies = np.asarray(energies, dtype=np.float64)
        widths = np.asarray(widths, dtype=np.float64)
        amplitudes = np.asarray(amplitudes, dtype=np.float64)

        # Validate inputs
        if not (len(energies) == len(widths) == len(amplitudes)):
            raise ValueError(
                f"Input arrays must have same length: "
                f"energies={len(energies)}, widths={len(widths)}, "
                f"amplitudes={len(amplitudes)}"
            )

        if len(energies) == 0:
            raise ValueError("Cannot calculate overlap for empty peak list")

        # Calculate overlap matrix using backend
        overlap_matrix = self._backend_instance.calculate_overlap_matrix(
            energies, widths, amplitudes
        )

        return OverlapResult(
            overlap_matrix=overlap_matrix,
            energies=energies,
            widths=widths,
            amplitudes=amplitudes,
        )

    def calculate_from_dataframe(
        self,
        df: pd.DataFrame,
        energy_col: str = "E",
        width_col: str = "width",
        amplitude_col: str = "OS",
    ) -> OverlapResult:
        """
        Calculate overlap matrix from DataFrame.

        Convenience method for working with pandas DataFrames containing peak data.

        Parameters
        ----------
        df : DataFrame
            Peak data with columns for energy, width, and amplitude.
        energy_col : str, default='E'
            Column name for peak energies.
        width_col : str, default='width'
            Column name for peak widths.
        amplitude_col : str, default='OS'
            Column name for peak amplitudes/oscillator strengths.

        Returns
        -------
        OverlapResult
            Container with overlap matrix and peak parameters.

        Raises
        ------
        KeyError
            If specified columns are not found in DataFrame.

        Examples
        --------
        >>> import pandas as pd
        >>> df = pd.DataFrame({
        ...     'E': [285.0, 286.0, 287.0],
        ...     'width': [0.5, 0.5, 0.5],
        ...     'OS': [1.0, 0.8, 0.6]
        ... })
        >>> calc = OverlapCalculator()
        >>> result = calc.calculate_from_dataframe(df)
        >>> result.overlap_matrix.shape
        (3, 3)
        """
        # Validate columns exist
        missing_cols = []
        for col in [energy_col, width_col, amplitude_col]:
            if col not in df.columns:
                missing_cols.append(col)

        if missing_cols:
            raise KeyError(
                f"Columns not found in DataFrame: {missing_cols}. "
                f"Available columns: {list(df.columns)}"
            )

        return self.calculate(
            energies=df[energy_col].values,
            widths=df[width_col].values,
            amplitudes=df[amplitude_col].values,
        )

    def __repr__(self) -> str:
        """Return string representation."""
        return f"OverlapCalculator(backend={self.backend!r})"
