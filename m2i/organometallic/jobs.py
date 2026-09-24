"""From a built complex to input files, and the words to describe its isomers."""

from __future__ import annotations

import platform
import re
from datetime import datetime, timezone
from pathlib import Path

from .. import __version__
from ..chem import metals
from ..crystal.coordination import analyse
from ..crystal.jobs import centred
from ..recognition.base import BackendError
from ..report import validate
from ..types import IssueLog, JobSpec
from ..writers import get_writer
from .build import BuiltComplex
from .drawing import OrganometallicDrawing
from .geometry import AMBIGUOUS, Arrangement, is_chiral, mirror_partners

SAFE = re.compile(r"[^A-Za-z0-9._-]+")


def describe(arrangement: Arrangement, drawing: OrganometallicDrawing) -> str:
    """'octahedral; trans: P2/P3, N1/C4...' in the drawing's own atom names."""
    pairs = sorted(
        " / ".join(sorted(_short(drawing, d) for d in pair)) for pair in arrangement.trans_pairs()
    )
    text = arrangement.polyhedron
    if pairs:
        text += "; trans: " + ", ".join(pairs)
    if is_chiral(arrangement, drawing.labels, drawing.chelated):
        text += " (chiral)"
    return text


def mirror_only(options: list[Arrangement], drawing: OrganometallicDrawing) -> bool:
    """The two best differ only as mirror images (a drawing without wedges)."""
    return mirror_partners(options[:2], drawing.labels, drawing.chelated).get(0) == 1


def ambiguous(options: list[Arrangement]) -> bool:
    """The drawing does not decide between the two best arrangements."""
    return len(options) > 1 and options[1].misfit - options[0].misfit < AMBIGUOUS


def default_name(species, source: Path) -> str:
    """What the files are called when no name is asked for."""
    return SAFE.sub("_", f"{species.formula}_{Path(source).stem}")


def write(
    built: BuiltComplex,
    drawing: OrganometallicDrawing,
    *,
    source: Path,
    charge: int,
    multiplicity: int,
    profile,
    output_dir: Path,
    log: IssueLog,
    oxidation: dict[str, int] | None = None,
    name: str | None = None,
) -> list[str]:
    species = built.species
    problem = metals.parity_problem(species.elements, charge, multiplicity)
    if problem:
        raise BackendError(f"impossible electronic state: {problem}")

    stem = SAFE.sub("_", name) if name else default_name(species, source)
    job = JobSpec(
        name=stem,
        title=(
            f"{species.formula} | built from {Path(source).name}, "
            f"{describe(built.arrangement, drawing)} | charge {charge} mult {multiplicity} "
            f"| m2i {__version__}"
        ),
        elements=species.elements,
        coords=centred(species).tolist(),
        charge=charge,
        multiplicity=multiplicity,
        profile=profile,
        metadata={"formula": species.formula, "source": Path(source).name},
    )
    output_dir = Path(output_dir)
    written = [str(get_writer(profile.program).write(job, output_dir / stem, log))]
    record = {
        "m2i_version": __version__,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python_version": platform.python_version(),
        "source": {"file": Path(source).name, "kind": "organometallic drawing"},
        "species": {
            "formula": species.formula,
            "atoms": species.n_atoms,
            "arrangement": describe(built.arrangement, drawing),
            "fit_to_drawing": round(built.arrangement.misfit, 3),
            "coordination": [line for c in analyse(species) for line in c.describe()],
        },
        "drawing": {
            "substituents_added": drawing.completed,
            "hydrogens_assumed": drawing.assumed,
        },
        "geometry": "distance geometry + UFF pre-optimisation: a starting point for a "
                    "quantum-chemical optimisation, not a final structure",
        "electronic_state": {
            "charge": charge,
            "multiplicity": multiplicity,
            "oxidation_states": oxidation or {},
            "decided_by": "user",
        },
        "notes": [message for _, _, message in drawing.notes + built.notes],
        "profile": profile.to_dict(),
        "issues": log.to_list(),
        "written_files": list(written),
    }
    written.append(str(validate.write_provenance(output_dir / f"{stem}.m2i.json", record)))
    return written


def _short(drawing: OrganometallicDrawing, donor: int) -> str:
    return drawing.names.get(donor, str(donor)).split(" of ")[0]
