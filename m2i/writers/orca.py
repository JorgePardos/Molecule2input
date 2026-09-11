"""ORCA input files (.inp)."""

from __future__ import annotations

from ..types import IssueLog, JobSpec
from .base import BaseWriter, check_basis_coverage, check_size, format_geometry

JOB_KEYWORDS = {
    "opt": "Opt",
    "freq": "Freq",
    "numfreq": "NumFreq",
    "sp": "SP",
    "optts": "OptTS",
    "irc": "IRC",
    "nmr": "NMR",
    "scan": "Opt",  # the scan itself lives in a %geom block
}


class OrcaWriter(BaseWriter):
    program = "orca"
    extension = ".inp"

    def render(self, job: JobSpec, log: IssueLog) -> str:
        check_basis_coverage(job, log)
        check_size(job, log)
        profile = job.profile
        resources = profile.resources

        lines = [f"# {job.title or job.name}"]
        lines.append(self._keyword_line(job))

        if resources.nproc > 1:
            lines.append(f"%pal nprocs {resources.nproc} end")
        lines.append(f"%maxcore {resources.orca_maxcore()}")

        solvent_block = self._solvent_block(job, log)
        if solvent_block:
            lines.append(solvent_block)

        for section in profile.extra_sections:
            lines.append(section.rstrip("\n"))

        lines.append("")
        lines.append(f"* xyz {job.charge} {job.multiplicity}")
        lines.append(format_geometry(job.elements, job.coords, width=16, decimals=8))
        lines.append("*")
        lines.append("")
        return "\n".join(lines)

    def _keyword_line(self, job: JobSpec) -> str:
        profile = job.profile
        parts = ["!"]
        if profile.method:
            parts.append(profile.method)
        if profile.basis:
            parts.append(profile.basis)
        if profile.dispersion:
            parts.append(profile.dispersion.upper())
        parts += [JOB_KEYWORDS.get(j, j) for j in profile.jobs if j]
        parts += [k for k in profile.extra_keywords if k]
        # CPCM is a keyword; SMD needs its own block (handled separately).
        if profile.solvent and profile.solvent.model != "smd":
            parts.append(f"CPCM({profile.solvent.name})")
        return " ".join(parts)

    def _solvent_block(self, job: JobSpec, log: IssueLog) -> str:
        solvent = job.profile.solvent
        if not solvent or solvent.model != "smd":
            return ""
        return (
            "%cpcm\n"
            "  smd true\n"
            f'  SMDsolvent "{solvent.name}"\n'
            "end"
        )
