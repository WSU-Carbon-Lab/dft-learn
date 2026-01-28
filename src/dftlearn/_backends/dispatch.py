"""Backend selection and dispatch logic for CPU/GPU computation."""

from typing import Literal

Backend = Literal["auto", "cpu", "gpu"]


def get_backend(backend: Backend = "auto"):
    """
    Get computation backend for analysis algorithms.

    This function provides automatic backend selection with graceful fallback
    from GPU to CPU when CuPy is not available.

    Parameters
    ----------
    backend : {'auto', 'cpu', 'gpu'}, default='auto'
        Backend selection strategy:

        - 'auto': Try GPU (CuPy) first, fallback to CPU (Numba) if unavailable
        - 'cpu': Force CPU backend (NumPy/Numba)
        - 'gpu': Force GPU backend (CuPy) - raises error if unavailable

    Returns
    -------
    backend : CPUBackend or GPUBackend
        Backend instance with methods for overlap calculation, clustering, etc.

    Raises
    ------
    RuntimeError
        If backend='gpu' is requested but CuPy is not available.
    ValueError
        If an unknown backend name is provided.

    Examples
    --------
    >>> # Automatic backend selection
    >>> backend = get_backend('auto')
    >>> overlap_matrix = backend.calculate_overlap_matrix(energies, widths, amplitudes)

    >>> # Force CPU backend
    >>> cpu_backend = get_backend('cpu')

    >>> # Force GPU backend (requires CuPy)
    >>> gpu_backend = get_backend('gpu')  # doctest: +SKIP

    Notes
    -----
    The GPU backend requires CuPy and a CUDA-capable GPU. For installation
    instructions, see: https://docs.cupy.dev/en/stable/install.html

    Performance considerations:
    - GPU is typically faster for large arrays (>1000 peaks)
    - CPU has lower overhead for small arrays (<1000 peaks)
    - Use 'auto' for best performance across different problem sizes
    """
    if backend == "cpu":
        from dftlearn._backends.cpu import CPUBackend

        return CPUBackend()

    elif backend == "gpu":
        try:
            from dftlearn._backends.gpu import GPUBackend

            return GPUBackend()
        except ImportError as e:
            raise RuntimeError(
                "GPU backend requested but CuPy is not available. "
                "Install with: pip install cupy-cuda12x\n"
                "See: https://docs.cupy.dev/en/stable/install.html"
            ) from e

    elif backend == "auto":
        try:
            import cupy as cp  # noqa: F401

            from dftlearn._backends.gpu import GPUBackend

            return GPUBackend()
        except ImportError:
            from dftlearn._backends.cpu import CPUBackend

            return CPUBackend()

    else:
        raise ValueError(
            f"Unknown backend: {backend!r}. "
            f"Expected one of: 'auto', 'cpu', 'gpu'"
        )


def should_use_gpu(array_size: int, threshold: int = 1000) -> bool:
    """
    Determine if GPU should be used based on array size.

    GPU acceleration is beneficial for large arrays where parallelization
    overhead is amortized. For small arrays, CPU execution is typically faster
    due to lower memory transfer overhead.

    Parameters
    ----------
    array_size : int
        Number of elements in the array to process.
    threshold : int, default=1000
        Minimum array size for GPU to be beneficial.

    Returns
    -------
    bool
        True if GPU should be used, False otherwise.

    Examples
    --------
    >>> should_use_gpu(500)
    False
    >>> should_use_gpu(5000)
    True
    >>> should_use_gpu(800, threshold=500)
    True

    Notes
    -----
    The default threshold of 1000 is a rough heuristic. Optimal threshold
    depends on:
    - GPU model and memory bandwidth
    - CPU model and number of cores
    - Algorithm complexity
    - Memory access patterns

    For your specific hardware, benchmark with different array sizes to
    determine the optimal threshold.
    """
    return array_size >= threshold
