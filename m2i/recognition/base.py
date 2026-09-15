"""The contract every recognition backend must satisfy.

Vision backends (DECIMER) live in their own virtual environments
and are driven over a subprocess, but they return exactly this object, so
nothing downstream needs to change when they are plugged in.
"""

from __future__ import annotations

from pathlib import Path
from typing import Protocol, runtime_checkable

from ..types import RecognitionResult


class BackendError(RuntimeError):
    """The backend could not produce a structure."""


@runtime_checkable
class OCSRBackend(Protocol):
    name: str
    #: Human-readable note on what this backend is good at.
    description: str
    #: True when this backend returns a molblock (2D layout + wedge bonds),
    #: which is what makes geometric stereochemistry possible.
    returns_molblock: bool

    def available(self) -> tuple[bool, str]:
        """(usable, reason). The reason explains how to install it if not."""

    def recognize(self, image_path: Path) -> RecognitionResult:
        """Read the image and return a structure."""
