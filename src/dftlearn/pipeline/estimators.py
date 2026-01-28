"""Base estimator classes for sklearn-compatible spectroscopy analysis."""

from abc import ABC, abstractmethod

from sklearn.base import BaseEstimator, ClusterMixin, TransformerMixin


class SpectrumTransformer(BaseEstimator, TransformerMixin, ABC):
    """
    Base class for spectrum data transformations.

    This abstract base class provides the foundation for sklearn-compatible
    transformers that process spectroscopy data. All transformers inherit
    from this class and must implement fit() and transform() methods.

    Inheriting Classes
    ------------------
    - EnergyFilter: Filter peaks by energy threshold
    - OscillatorStrengthFilter: Filter peaks by oscillator strength
    - NormalizedOSTransformer: Add normalized oscillator strength column
    - ClusterGaussianFitter: Fit Gaussians to clustered peaks

    Methods
    -------
    fit(X, y=None)
        Fit the transformer to data.
    transform(X)
        Transform the data.
    fit_transform(X, y=None)
        Fit and transform in one step (provided by TransformerMixin).

    Notes
    -----
    This class automatically provides get_params() and set_params() methods
    through sklearn's BaseEstimator, enabling:
    - Hyperparameter tuning with GridSearchCV/RandomizedSearchCV
    - Pipeline composition
    - Model persistence (pickling)

    Examples
    --------
    >>> from sklearn.base import BaseEstimator, TransformerMixin
    >>> class MyTransformer(SpectrumTransformer):
    ...     def __init__(self, threshold=1.0):
    ...         self.threshold = threshold
    ...
    ...     def fit(self, X, y=None):
    ...         return self
    ...
    ...     def transform(self, X):
    ...         return X[X > self.threshold]
    """

    def get_params(self, deep=True):
        """
        Get parameters for this estimator.

        Parameters
        ----------
        deep : bool, default=True
            If True, will return the parameters for this estimator and
            contained subobjects that are estimators.

        Returns
        -------
        params : dict
            Parameter names mapped to their values.
        """
        return super().get_params(deep)

    def set_params(self, **params):
        """
        Set the parameters of this estimator.

        Parameters
        ----------
        **params : dict
            Estimator parameters.

        Returns
        -------
        self : object
            Estimator instance.
        """
        return super().set_params(**params)

    @abstractmethod
    def fit(self, X, y=None):
        """
        Fit the transformer to data.

        Parameters
        ----------
        X : array-like or DataFrame
            Training data.
        y : array-like, optional
            Target values (ignored, exists for sklearn compatibility).

        Returns
        -------
        self : object
            Fitted transformer instance.
        """
        pass

    @abstractmethod
    def transform(self, X):
        """
        Transform the data.

        Parameters
        ----------
        X : array-like or DataFrame
            Data to transform.

        Returns
        -------
        X_transformed : array-like or DataFrame
            Transformed data.
        """
        pass


class SpectrumClusterer(BaseEstimator, ClusterMixin, ABC):
    """
    Base class for spectrum clustering algorithms.

    This abstract base class provides the foundation for sklearn-compatible
    clustering algorithms for spectroscopy data. All clusterers inherit from
    this class and must implement fit() and predict() methods.

    Inheriting Classes
    ------------------
    - OverlapClusterer: Cluster peaks by Gaussian overlap

    Attributes
    ----------
    labels_ : ndarray, shape (n_samples,)
        Cluster labels for each sample (set during fit).
    n_clusters_ : int
        Number of clusters found (set during fit).

    Methods
    -------
    fit(X, y=None)
        Fit the clusterer to data.
    predict(X)
        Predict cluster labels.
    fit_predict(X, y=None)
        Fit and predict in one step (provided by ClusterMixin).

    Notes
    -----
    Unlike sklearn's standard clusterers, spectrum clusterers typically
    work with DataFrames containing peak information (energy, width, amplitude)
    rather than raw feature matrices.

    The fit() method computes overlap matrices and performs clustering, storing
    the results in labels_ and n_clusters_ attributes.

    Examples
    --------
    >>> import pandas as pd
    >>> from sklearn.base import BaseEstimator, ClusterMixin
    >>> class MyClusterer(SpectrumClusterer):
    ...     def __init__(self, threshold=50.0):
    ...         self.threshold = threshold
    ...
    ...     def fit(self, X, y=None):
    ...         # Perform clustering
    ...         self.labels_ = [0, 0, 1, 1]  # Example labels
    ...         self.n_clusters_ = 2
    ...         return self
    ...
    ...     def predict(self, X):
    ...         return self.labels_
    """

    @abstractmethod
    def fit(self, X, y=None):
        """
        Fit the clusterer to data.

        Parameters
        ----------
        X : DataFrame
            Peak data with columns like 'E' (energy), 'width', 'OS' (oscillator strength).
        y : array-like, optional
            Target values (ignored, exists for sklearn compatibility).

        Returns
        -------
        self : object
            Fitted clusterer instance with labels_ and n_clusters_ attributes set.
        """
        pass

    @abstractmethod
    def predict(self, X):
        """
        Predict cluster labels for data.

        Parameters
        ----------
        X : DataFrame
            Peak data to predict clusters for.

        Returns
        -------
        labels : ndarray, shape (n_samples,)
            Cluster labels for each peak.
        """
        pass


class SpectrumEstimator(BaseEstimator, ABC):
    """
    Base class for spectrum estimation algorithms.

    This abstract base class provides the foundation for sklearn-compatible
    estimators that don't fit into the transformer or clusterer categories.
    Examples include complete analysis pipelines and fitting algorithms.

    Inheriting Classes
    ------------------
    - NEXAFSPipeline: Complete NEXAFS analysis workflow
    - MultiSpectrumFitter: Multi-spectrum Gaussian fitting
    - TensorFitter: Molecular absorption tensor fitting
    - StepEdgeBuilder: Theoretical absorption edge construction

    Methods
    -------
    fit(X, y=None)
        Fit the estimator to data.
    predict(X)
        Make predictions using the fitted estimator (optional).
    score(X, y)
        Return a score for the predictions (optional).

    Notes
    -----
    This is a more flexible base class for algorithms that don't follow
    the standard transformer or clusterer patterns. Subclasses should
    implement fit() at minimum, and optionally predict() and score()
    depending on the algorithm's purpose.

    Examples
    --------
    >>> from sklearn.base import BaseEstimator
    >>> class MyEstimator(SpectrumEstimator):
    ...     def __init__(self, parameter=1.0):
    ...         self.parameter = parameter
    ...
    ...     def fit(self, X, y=None):
    ...         # Perform fitting
    ...         self.is_fitted_ = True
    ...         return self
    ...
    ...     def predict(self, X):
    ...         # Make predictions
    ...         return X * self.parameter
    """

    @abstractmethod
    def fit(self, X, y=None):
        """
        Fit the estimator to data.

        Parameters
        ----------
        X : array-like or DataFrame
            Training data.
        y : array-like, optional
            Target values.

        Returns
        -------
        self : object
            Fitted estimator instance.
        """
        pass
