"""Data passed between pipeline stages.

Nothing here knows about RDKit internals except ``MoleculeSpec.mol``; the
writers only ever see plain elements/coordinates so a new backend or a new
quantum-chemistry program can be added without touching the middle.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable, Sequence

INFO = "info"
WARNING = "warning"
ERROR = "error"


@dataclass
class Issue:
    """Something the user should know before trusting the generated input."""

    level: str
    code: str
    message: str

    def __str__(self) -> str:  # pragma: no cover - trivial
        return f"[{self.level}] {self.message}"


class IssueLog:
    """Ordered collection of issues, shared by every stage of the pipeline."""

    def __init__(self, issues: Iterable[Issue] = ()) -> None:
        self.issues: list[Issue] = list(issues)

    def add(self, level: str, code: str, message: str) -> Issue:
        issue = Issue(level, code, message)
        self.issues.append(issue)
        return issue

    def info(self, code: str, message: str) -> Issue:
        return self.add(INFO, code, message)

    def warn(self, code: str, message: str) -> Issue:
        return self.add(WARNING, code, message)

    def error(self, code: str, message: str) -> Issue:
        return self.add(ERROR, code, message)

    def extend(self, other: "IssueLog | Iterable[Issue]") -> None:
        self.issues.extend(other.issues if isinstance(other, IssueLog) else other)

    @property
    def warnings(self) -> list[Issue]:
        return [i for i in self.issues if i.level == WARNING]

    @property
    def errors(self) -> list[Issue]:
        return [i for i in self.issues if i.level == ERROR]

    def has_errors(self) -> bool:
        return bool(self.errors)

    def to_list(self) -> list[dict]:
        return [{"level": i.level, "code": i.code, "message": i.message} for i in self.issues]

    def __iter__(self):
        return iter(self.issues)

    def __len__(self) -> int:
        return len(self.issues)


@dataclass
class RecognitionResult:
    """What an OCSR backend returns. The contract Phase 2 backends must meet."""

    smiles: str
    molblock: str | None = None
    confidence: float | None = None
    backend: str = "unknown"
    raw: dict | None = None

    def to_dict(self) -> dict:
        return {
            "smiles": self.smiles,
            "backend": self.backend,
            "confidence": self.confidence,
            "has_molblock": self.molblock is not None,
        }


@dataclass
class StereoCenter:
    """A tetrahedral (or otherwise atom-centred) stereo element."""

    atom_index: int
    symbol: str
    label: str | None  # CIP descriptor: 'R', 'S', 'r', 's', or None if unspecified
    specified: bool

    def describe(self) -> str:
        where = f"{self.symbol}{self.atom_index + 1}"
        return f"{where}: {self.label}" if self.specified else f"{where}: unspecified"


@dataclass
class StereoBond:
    """A double bond whose configuration is (or should be) defined."""

    begin_index: int
    end_index: int
    label: str | None  # 'E', 'Z', or None if unspecified
    specified: bool

    def describe(self) -> str:
        where = f"{self.begin_index + 1}={self.end_index + 1}"
        return f"{where}: {self.label}" if self.specified else f"{where}: unspecified"


@dataclass
class StereoSummary:
    centers: list[StereoCenter] = field(default_factory=list)
    bonds: list[StereoBond] = field(default_factory=list)
    #: Centres RDKit considered potentially stereogenic that chemistry rules
    #: out -- see stereo.summarize(). Kept for transparency, never warned about.
    ignored_centers: list[StereoCenter] = field(default_factory=list)

    @property
    def unspecified_centers(self) -> list[StereoCenter]:
        return [c for c in self.centers if not c.specified]

    @property
    def unspecified_bonds(self) -> list[StereoBond]:
        return [b for b in self.bonds if not b.specified]

    @property
    def is_complete(self) -> bool:
        return not self.unspecified_centers and not self.unspecified_bonds

    def fingerprint(self) -> dict[str, str]:
        """Assigned descriptors keyed by position, for the 2D vs 3D comparison."""
        out: dict[str, str] = {}
        for c in self.centers:
            if c.specified and c.label:
                out[f"atom:{c.atom_index}"] = c.label
        for b in self.bonds:
            if b.specified and b.label:
                out[f"bond:{b.begin_index}-{b.end_index}"] = b.label
        return out

    def describe(self) -> str:
        parts = [c.describe() for c in self.centers] + [b.describe() for b in self.bonds]
        return "; ".join(parts) if parts else "no stereocentres"

    def to_dict(self) -> dict:
        return {
            "centers": [
                {
                    "atom_index": c.atom_index,
                    "symbol": c.symbol,
                    "label": c.label,
                    "specified": c.specified,
                }
                for c in self.centers
            ],
            "bonds": [
                {
                    "begin_index": b.begin_index,
                    "end_index": b.end_index,
                    "label": b.label,
                    "specified": b.specified,
                }
                for b in self.bonds
            ],
            "ignored_centers": [
                {"atom_index": c.atom_index, "symbol": c.symbol}
                for c in self.ignored_centers
            ],
        }


@dataclass
class MoleculeSpec:
    """A validated molecule, ready to be embedded in 3D."""

    mol: Any  # rdkit.Chem.Mol, without explicit hydrogens
    smiles: str  # canonical, isomeric
    inchikey: str
    formula: str
    charge: int
    multiplicity: int
    stereo: StereoSummary
    name: str = "molecule"
    source: RecognitionResult | None = None

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "smiles": self.smiles,
            "inchikey": self.inchikey,
            "formula": self.formula,
            "charge": self.charge,
            "multiplicity": self.multiplicity,
            "stereo": self.stereo.to_dict(),
            "source": self.source.to_dict() if self.source else None,
        }


@dataclass
class Conformer:
    """One 3D geometry with explicit hydrogens."""

    index: int
    elements: Sequence[str]
    coords: Sequence[Sequence[float]]  # (n_atoms, 3) in angstrom
    energy: float | None = None  # force-field energy
    energy_unit: str = "kcal/mol"
    force_field: str = "MMFF94s"
    relative_energy: float | None = None

    @property
    def n_atoms(self) -> int:
        return len(self.elements)

    def to_dict(self) -> dict:
        return {
            "index": self.index,
            "n_atoms": self.n_atoms,
            "energy": self.energy,
            "relative_energy": self.relative_energy,
            "energy_unit": self.energy_unit,
            "force_field": self.force_field,
        }


@dataclass
class JobSpec:
    """Everything a writer needs: geometry + electronic state + recipe."""

    name: str
    title: str
    elements: Sequence[str]
    coords: Sequence[Sequence[float]]
    charge: int
    multiplicity: int
    profile: Any  # m2i.config.JobProfile
    metadata: dict = field(default_factory=dict)
    # Optional RDKit molecule carrying this exact geometry, for the writers
    # that need connectivity (SDF). Never required by the QM writers.
    mol: Any = None

    @property
    def n_atoms(self) -> int:
        return len(self.elements)


@dataclass
class PipelineResult:
    """Outcome of a full run, including everything needed for provenance."""

    molecule: MoleculeSpec | None
    conformers: list[Conformer] = field(default_factory=list)
    written_files: list[str] = field(default_factory=list)
    issues: IssueLog = field(default_factory=IssueLog)
    profile: Any = None

    @property
    def ok(self) -> bool:
        return self.molecule is not None and not self.issues.has_errors()
