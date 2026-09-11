"""The backend where the chemist is the recogniser.

Always available, always the fallback, and the reference implementation of the
backend contract. It is also what makes m2i useful before any vision model is
installed: draw the molecule, read it yourself, type the SMILES.
"""

from __future__ import annotations

from pathlib import Path

from ..types import RecognitionResult
from .base import BackendError


class ManualBackend:
    name = "manual"
    description = "SMILES (or molblock) supplied by the user"
    returns_molblock = True  # if the user hands one over

    def __init__(self, smiles: str | None = None, molblock: str | None = None) -> None:
        self.smiles = smiles
        self.molblock = molblock

    def available(self) -> tuple[bool, str]:
        return True, "always available"

    def recognize(self, image_path: Path | None = None, **options) -> RecognitionResult:
        if self.molblock:
            return RecognitionResult(
                smiles=self.smiles or "",
                molblock=self.molblock,
                confidence=1.0,
                backend=self.name,
            )
        if not self.smiles:
            raise BackendError(
                "the manual backend needs a --smiles value (or a molfile via "
                "--molfile); no vision model is installed yet"
            )
        return RecognitionResult(
            smiles=self.smiles, molblock=None, confidence=1.0, backend=self.name
        )


def from_molfile(path: Path, smiles: str | None = None) -> ManualBackend:
    """Build a manual backend from an existing .mol/.sdf drawn in ChemDraw."""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    # An SDF may hold several records; only the first is used.
    molblock = text.split("$$$$")[0]
    return ManualBackend(smiles=smiles, molblock=molblock)
