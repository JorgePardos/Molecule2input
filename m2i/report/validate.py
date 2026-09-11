"""Provenance records and human-readable issue reporting.

Every generated input gets a sidecar JSON so that months later it is still
possible to answer "where did this geometry come from, and did anything warn
me at the time?".
"""

from __future__ import annotations

import json
import platform
from datetime import datetime, timezone
from pathlib import Path

from rdkit import rdBase

from .. import __version__
from ..types import ERROR, INFO, WARNING, IssueLog, PipelineResult

LEVEL_PREFIX = {INFO: "  -", WARNING: "  !", ERROR: "  X"}


def provenance(result: PipelineResult, extra: dict | None = None) -> dict:
    """Everything needed to reproduce or audit this run."""
    molecule = result.molecule
    record = {
        "m2i_version": __version__,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "rdkit_version": rdBase.rdkitVersion,
        "python_version": platform.python_version(),
        "molecule": molecule.to_dict() if molecule else None,
        "profile": result.profile.to_dict() if result.profile else None,
        "conformers": [c.to_dict() for c in result.conformers],
        "written_files": result.written_files,
        "issues": result.issues.to_list(),
    }
    if extra:
        record.update(extra)
    return record


def write_provenance(path: Path, record: dict) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return path


def format_issues(log: IssueLog, *, verbose: bool = False) -> str:
    """Render the issue log for the terminal, warnings and errors first."""
    if not len(log):
        return ""

    shown = [i for i in log if verbose or i.level != INFO]
    if not shown:
        return ""

    order = {ERROR: 0, WARNING: 1, INFO: 2}
    shown.sort(key=lambda i: order.get(i.level, 3))

    lines = []
    for issue in shown:
        prefix = LEVEL_PREFIX.get(issue.level, "  -")
        lines.append(f"{prefix} {issue.message}")
    return "\n".join(lines)


def summarize(result: PipelineResult) -> str:
    """One-paragraph summary shown before the confirmation prompt."""
    molecule = result.molecule
    if molecule is None:
        return "No molecule was produced."

    lines = [
        f"  SMILES        {molecule.smiles}",
        f"  Formula       {molecule.formula}",
        f"  InChIKey      {molecule.inchikey or '(unavailable)'}",
        f"  Charge/mult   {molecule.charge} / {molecule.multiplicity}",
        f"  Stereo        {molecule.stereo.describe()}",
    ]
    if result.conformers:
        energies = ", ".join(
            f"{c.relative_energy:+.2f}" if c.relative_energy is not None else "?"
            for c in result.conformers
        )
        lines.append(
            f"  Conformers    {len(result.conformers)} "
            f"(relative {result.conformers[0].force_field} energies, kcal/mol: {energies})"
        )
    return "\n".join(lines)
