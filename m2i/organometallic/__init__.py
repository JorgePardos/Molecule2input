"""Metal complexes from drawings: read the page, pick the isomer, build in 3D.

The route for complexes drawn in ChemDraw or as molfiles. For complexes with a
crystal structure, the CIF route (``m2i.crystal``) is better still: there the
geometry is measured, not built.
"""

from __future__ import annotations

from .build import BuiltComplex, build
from .drawing import OrganometallicDrawing, has_metal, read, read_raw
from .geometry import Arrangement, arrangements, is_chiral, mirror_partners

__all__ = [
    "Arrangement",
    "BuiltComplex",
    "OrganometallicDrawing",
    "arrangements",
    "build",
    "has_metal",
    "is_chiral",
    "mirror_partners",
    "read",
    "read_raw",
]
