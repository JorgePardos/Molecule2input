"""Optical chemical structure recognition: image -> molecule.

This is the only statistically uncertain layer of m2i. Everything downstream is
deterministic RDKit, which is why the human verification step sits exactly at
this boundary.
"""

from __future__ import annotations

from .base import BackendError, OCSRBackend
from .registry import available_backends, describe_backends, recognize

__all__ = [
    "BackendError",
    "OCSRBackend",
    "available_backends",
    "describe_backends",
    "recognize",
]
