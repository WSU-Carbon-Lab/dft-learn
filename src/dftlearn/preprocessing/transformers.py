"""Sklearn-compatible transformers for preprocessing NEXAFS data."""

import pandas as pd
from sklearn.utils.validation import check_is_fitted

from dftlearn.pipeline.estimators import SpectrumTransformer


class EnergyFilter(SpectrumTransformer):
    """
    Filter spectral peaks by energy threshold.

    This transformer removes peaks above a specified energy threshold,
    useful for focusing analysis on a specific energy range.

    Parameters
    ----------
    energy_threshold : float, default=320.0
        Maximum energy (eV) to include. Peaks above this are removed.
    verbose : bool, default=False
        If True, print filtering statistics.

    Attributes
    ----------
    n_peaks_before_ : int
        Number of peaks before filtering (set during fit).
    n_peaks_after_ : int
        Number of peaks after filtering (set during fit).

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({'E': [285, 290, 325], 'OS': [1.0, 0.8, 0.5]})
    >>> filter_obj = EnergyFilter(energy_threshold=300.0)
    >>> result = filter_obj.fit_transform(df)
    >>> len(result)
    2

    Notes
    -----
    This transformer is stateless - fit() does not learn anything from the data.
    It is included for sklearn pipeline compatibility.
    """

    def __init__(self, energy_threshold: float = 320.0, verbose: bool = False):
        """
        Initialize energy filter.

        Parameters
        ----------
        energy_threshold : float, default=320.0
            Maximum energy to include.
        verbose : bool, default=False
            Whether to print filtering statistics.
        """
        self.energy_threshold = energy_threshold
        self.verbose = verbose

    def fit(self, X, y=None):
        """
        Fit the transformer (no-op, included for sklearn compatibility).

        Parameters
        ----------
        X : DataFrame
            Peak data with 'E' column.
        y : array-like, optional
            Ignored.

        Returns
        -------
        self : object
            Fitted transformer.
        """
        if not isinstance(X, pd.DataFrame):
            raise TypeError(f"X must be a pandas DataFrame, got {type(X).__name__}")

        if "E" not in X.columns:
            raise ValueError(
                f"DataFrame must have 'E' column. Available columns: {list(X.columns)}"
            )

        self.n_peaks_before_ = len(X)
        self.n_peaks_after_ = len(X[X["E"] <= self.energy_threshold])

        if self.verbose:
            print(f"Original data: {self.n_peaks_before_} peaks")
            print(
                f"After energy filter (E <= {self.energy_threshold}): "
                f"{self.n_peaks_after_} peaks"
            )

        return self

    def transform(self, X):
        """
        Filter peaks by energy threshold.

        Parameters
        ----------
        X : DataFrame
            Peak data with 'E' column.

        Returns
        -------
        X_filtered : DataFrame
            Peaks with E <= energy_threshold.
        """
        if not isinstance(X, pd.DataFrame):
            raise TypeError(f"X must be a pandas DataFrame, got {type(X).__name__}")

        if "E" not in X.columns:
            raise ValueError(
                f"DataFrame must have 'E' column. Available columns: {list(X.columns)}"
            )

        return X[X["E"] <= self.energy_threshold].copy()


class NormalizedOSTransformer(SpectrumTransformer):
    """
    Add normalized oscillator strength column.

    This transformer adds a 'normalized_OS' column to the DataFrame,
    normalizing oscillator strengths by the maximum value.

    Parameters
    ----------
    os_column : str, default='OS'
        Column name for oscillator strength.
    normalized_column : str, default='normalized_OS'
        Name for the new normalized column.
    verbose : bool, default=False
        If True, print normalization statistics.

    Attributes
    ----------
    max_os_ : float
        Maximum oscillator strength in fitted data (set during fit).

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({'OS': [1.0, 0.5, 0.25]})
    >>> transformer = NormalizedOSTransformer()
    >>> result = transformer.fit_transform(df)
    >>> result['normalized_OS'].tolist()
    [1.0, 0.5, 0.25]

    Notes
    -----
    The normalization is: normalized_OS = OS / max(OS)

    This ensures the maximum normalized value is 1.0, making it easy to
    filter by percentage thresholds.
    """

    def __init__(
        self,
        os_column: str = "OS",
        normalized_column: str = "normalized_OS",
        verbose: bool = False,
    ):
        """
        Initialize normalized OS transformer.

        Parameters
        ----------
        os_column : str, default='OS'
            Column name for oscillator strength.
        normalized_column : str, default='normalized_OS'
            Name for normalized column.
        verbose : bool, default=False
            Whether to print statistics.
        """
        self.os_column = os_column
        self.normalized_column = normalized_column
        self.verbose = verbose

    def fit(self, X, y=None):
        """
        Compute maximum oscillator strength.

        Parameters
        ----------
        X : DataFrame
            Peak data with oscillator strength column.
        y : array-like, optional
            Ignored.

        Returns
        -------
        self : object
            Fitted transformer.
        """
        if not isinstance(X, pd.DataFrame):
            raise TypeError(f"X must be a pandas DataFrame, got {type(X).__name__}")

        if self.os_column not in X.columns:
            raise ValueError(
                f"DataFrame must have '{self.os_column}' column. "
                f"Available columns: {list(X.columns)}"
            )

        self.max_os_ = X[self.os_column].max()

        if self.verbose:
            print(f"Maximum {self.os_column}: {self.max_os_:.6f}")

        return self

    def transform(self, X):
        """
        Add normalized oscillator strength column.

        Parameters
        ----------
        X : DataFrame
            Peak data.

        Returns
        -------
        X_transformed : DataFrame
            Input data with normalized_OS column added.
        """
        check_is_fitted(self, ["max_os_"])

        if not isinstance(X, pd.DataFrame):
            raise TypeError(f"X must be a pandas DataFrame, got {type(X).__name__}")

        if self.os_column not in X.columns:
            raise ValueError(
                f"DataFrame must have '{self.os_column}' column. "
                f"Available columns: {list(X.columns)}"
            )

        X_transformed = X.copy()
        X_transformed[self.normalized_column] = X[self.os_column] / self.max_os_

        return X_transformed


class OscillatorStrengthFilter(SpectrumTransformer):
    """
    Filter peaks by normalized oscillator strength threshold.

    This transformer filters peaks based on a percentage threshold of the
    maximum oscillator strength. It automatically normalizes the data during
    fit() and filters during transform().

    Parameters
    ----------
    os_threshold_percent : float, default=2.0
        Minimum oscillator strength as percentage of maximum (0-100).
    os_column : str, default='OS'
        Column name for oscillator strength.
    verbose : bool, default=False
        If True, print filtering statistics.

    Attributes
    ----------
    max_os_ : float
        Maximum oscillator strength in fitted data (set during fit).
    os_threshold_ : float
        Absolute threshold value (set during fit).
    n_peaks_before_ : int
        Number of peaks before filtering (set during fit).
    n_peaks_after_ : int
        Number of peaks after filtering (set during fit).

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({'OS': [1.0, 0.5, 0.01]})
    >>> filter_obj = OscillatorStrengthFilter(os_threshold_percent=10.0)
    >>> result = filter_obj.fit_transform(df)
    >>> len(result)
    2

    Notes
    -----
    This filter is useful for removing weak peaks that don't significantly
    contribute to the spectrum. A threshold of 2% means peaks with oscillator
    strength less than 2% of the maximum are removed.
    """

    def __init__(
        self,
        os_threshold_percent: float = 2.0,
        os_column: str = "OS",
        verbose: bool = False,
    ):
        """
        Initialize oscillator strength filter.

        Parameters
        ----------
        os_threshold_percent : float, default=2.0
            Minimum OS as percentage of maximum.
        os_column : str, default='OS'
            Column name for oscillator strength.
        verbose : bool, default=False
            Whether to print statistics.
        """
        self.os_threshold_percent = os_threshold_percent
        self.os_column = os_column
        self.verbose = verbose

    def fit(self, X, y=None):
        """
        Compute maximum OS and threshold.

        Parameters
        ----------
        X : DataFrame
            Peak data with oscillator strength column.
        y : array-like, optional
            Ignored.

        Returns
        -------
        self : object
            Fitted transformer.
        """
        if not isinstance(X, pd.DataFrame):
            raise TypeError(f"X must be a pandas DataFrame, got {type(X).__name__}")

        if self.os_column not in X.columns:
            raise ValueError(
                f"DataFrame must have '{self.os_column}' column. "
                f"Available columns: {list(X.columns)}"
            )

        self.max_os_ = X[self.os_column].max()
        self.os_threshold_ = (self.os_threshold_percent / 100.0) * self.max_os_

        self.n_peaks_before_ = len(X)
        self.n_peaks_after_ = len(X[X[self.os_column] >= self.os_threshold_])

        if self.verbose:
            print(f"Maximum {self.os_column}: {self.max_os_:.6f}")
            print(
                f"Threshold ({self.os_threshold_percent}%): {self.os_threshold_:.6f}"
            )
            print(
                f"After {self.os_threshold_percent}% OS threshold: "
                f"{self.n_peaks_after_} peaks"
            )

        return self

    def transform(self, X):
        """
        Filter peaks by oscillator strength threshold.

        Parameters
        ----------
        X : DataFrame
            Peak data.

        Returns
        -------
        X_filtered : DataFrame
            Peaks with OS >= threshold.
        """
        check_is_fitted(self, ["os_threshold_"])

        if not isinstance(X, pd.DataFrame):
            raise TypeError(f"X must be a pandas DataFrame, got {type(X).__name__}")

        if self.os_column not in X.columns:
            raise ValueError(
                f"DataFrame must have '{self.os_column}' column. "
                f"Available columns: {list(X.columns)}"
            )

        return X[X[self.os_column] >= self.os_threshold_].copy()
