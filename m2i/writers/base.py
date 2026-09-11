"""Writer interface plus the geometry formatting and basis-set sanity checks
shared by every program.
"""

from __future__ import annotations

from pathlib import Path
from typing import Protocol, Sequence

from ..types import IssueLog, JobSpec

ATOMIC_NUMBERS = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9,
    "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17,
    "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25,
    "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32,
    "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39,
    "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46,
    "Ag": 47, "Cd": 48, "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I": 53,
    "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Hf": 72, "Ta": 73, "W": 74,
    "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81,
    "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86,
}

# Highest atomic number a basis-set family covers with an all-electron or
# built-in-ECP definition. Deliberately conservative: this drives a warning,
# never a hard failure.
BASIS_MAX_Z = {
    "sto-3g": 53,
    "3-21g": 54,
    "6-31g": 36,
    "6-311g": 36,
    "cc-pv": 36,
    "aug-cc-pv": 36,
    "def2": 86,
    "lanl2": 86,
    "sdd": 86,
    "midix": 36,
}


class WriterError(ValueError):
    """The job cannot be written as requested."""


class InputWriter(Protocol):
    """What every program writer must provide."""

    program: str
    extension: str

    def render(self, job: JobSpec, log: IssueLog) -> str:
        """Return the complete file contents."""

    def write(self, job: JobSpec, path: Path, log: IssueLog) -> Path:
        """Write the file and return the path actually used."""


class BaseWriter:
    """Shared implementation: rendering to text, then to disk."""

    program = "base"
    extension = ".txt"

    def render(self, job: JobSpec, log: IssueLog) -> str:  # pragma: no cover
        raise NotImplementedError

    def write(self, job: JobSpec, path: Path, log: IssueLog) -> Path:
        path = Path(path)
        if path.suffix.lower() != self.extension:
            path = path.with_suffix(self.extension)
        path.parent.mkdir(parents=True, exist_ok=True)
        # Quantum-chemistry programs are line-oriented and several are picky
        # about stray carriage returns, so newlines are written explicitly.
        path.write_text(self.render(job, log), encoding="utf-8", newline="\n")
        return path


def format_geometry(
    elements: Sequence[str],
    coords: Sequence[Sequence[float]],
    *,
    width: int = 16,
    decimals: int = 8,
) -> str:
    """Fixed-width Cartesian block, one atom per line."""
    lines = []
    for symbol, (x, y, z) in zip(elements, coords):
        lines.append(
            f" {symbol:<3s}"
            f"{x:>{width}.{decimals}f}"
            f"{y:>{width}.{decimals}f}"
            f"{z:>{width}.{decimals}f}"
        )
    return "\n".join(lines)


def check_basis_coverage(job: JobSpec, log: IssueLog) -> None:
    """Warn when the basis-set family is unlikely to cover a heavy element."""
    basis = (job.profile.basis or "").lower().replace(" ", "")
    if not basis or basis in ("gen", "genecp"):
        return

    limit = None
    for family, max_z in BASIS_MAX_Z.items():
        if basis.startswith(family) or family in basis:
            limit = max_z
            break
    if limit is None:
        return

    heavy = sorted(
        {
            symbol
            for symbol in set(job.elements)
            if ATOMIC_NUMBERS.get(symbol, 0) > limit
        },
        key=lambda s: ATOMIC_NUMBERS.get(s, 0),
    )
    if heavy:
        log.warn(
            "basis.coverage",
            f"Basis {job.profile.basis} probably has no definition for "
            f"{', '.join(heavy)}. Switch to a def2 basis, or use gen/genecp with "
            "an explicit basis block in the profile's extra_sections.",
        )


def check_size(job: JobSpec, log: IssueLog, *, soft_limit: int = 150) -> None:
    if job.n_atoms > soft_limit:
        log.info(
            "job.size",
            f"{job.n_atoms} atoms: expect a long run at {job.profile.method}/"
            f"{job.profile.basis}.",
        )


# -- registry ------------------------------------------------------------


def known_programs() -> tuple[str, ...]:
    return ("gaussian", "orca", "xyz", "sdf")


def get_writer(program: str) -> InputWriter:
    program = (program or "").lower()
    if program == "gaussian":
        from .gaussian import GaussianWriter

        return GaussianWriter()
    if program == "orca":
        from .orca import OrcaWriter

        return OrcaWriter()
    if program == "xyz":
        from .geometry import XyzWriter

        return XyzWriter()
    if program == "sdf":
        from .geometry import SdfWriter

        return SdfWriter()
    raise WriterError(
        f"no writer for program {program!r}; known: {', '.join(known_programs())}"
    )
