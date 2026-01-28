"""Backend dispatch system for CPU (Numba) and GPU (CuPy) acceleration."""

from dftlearn._backends.dispatch import get_backend

__all__ = ["get_backend"]
