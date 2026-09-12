"""Gaussian input files (.gjf).

Gaussian's parser is positional: Link0, route, blank, title, blank, charge and
multiplicity, geometry, blank, optional sections each followed by a blank line.
A missing blank line is the single most common reason a hand-made input dies on
line 1, so the layout is built explicitly here rather than by string concat.
"""

from __future__ import annotations

from ..types import IssueLog, JobSpec
from .base import BaseWriter, check_size, format_geometry, report_ecp, uncovered_elements

# Job keywords that need no argument and map straight through.
JOB_KEYWORDS = {
    "opt": "opt",
    "freq": "freq",
    "sp": "sp",
    "nmr": "nmr",
    "td": "td",
    "irc": "irc",
    "scan": "scan",
    "stable": "stable=opt",
    "polar": "polar",
}


class GaussianWriter(BaseWriter):
    program = "gaussian"
    extension = ".gjf"

    def render(self, job: JobSpec, log: IssueLog) -> str:
        check_size(job, log)
        self._check_gen_basis(job, log)
        self._ecp = self._ecp_block(job, log)

        blocks = [self._main_step(job, log)]
        for extra in job.profile.link_jobs:
            blocks.append(self._link_step(job, extra, log))
        return "--Link1--\n".join(blocks)

    # -- steps -----------------------------------------------------------

    def _main_step(self, job: JobSpec, log: IssueLog) -> str:
        lines: list[str] = []
        lines += self._link0(job)
        lines.append(self._route(job, job.profile.jobs, job.profile.extra_keywords))
        lines.append("")
        lines.append(job.title or job.name)
        lines.append("")
        lines.append(f"{job.charge} {job.multiplicity}")
        lines.append(format_geometry(job.elements, job.coords))
        lines.append("")  # terminates the geometry
        if self._ecp:
            lines.append(self._ecp)
            lines.append("")
        for section in job.profile.extra_sections:
            lines.append(section.rstrip("\n"))
            lines.append("")
        return "\n".join(lines) + "\n"

    def _link_step(self, job: JobSpec, extra: dict, log: IssueLog) -> str:
        """A follow-up job reading the geometry and wavefunction from the .chk."""
        profile = job.profile
        jobs = extra.get("jobs", ())
        if isinstance(jobs, str):
            jobs = jobs.split()
        keywords = list(extra.get("extra_keywords", ()))
        # Without these the second step would re-read the (absent) geometry.
        keywords = ["geom=check", "guess=read"] + keywords

        method = extra.get("method", profile.method)
        basis = extra.get("basis", profile.basis)

        lines = self._link0(job)
        lines.append(
            self._route(
                job,
                tuple(jobs),
                tuple(keywords),
                method=method,
                basis=basis,
                dispersion=extra.get("dispersion", profile.dispersion),
            )
        )
        lines.append("")
        lines.append(f"{job.title or job.name} - {' '.join(jobs) or 'follow-up'}")
        lines.append("")
        lines.append(f"{job.charge} {job.multiplicity}")
        lines.append("")
        if self._ecp and basis == profile.basis:
            lines.append(self._ecp)
            lines.append("")
        for section in extra.get("extra_sections", ()):
            lines.append(str(section).rstrip("\n"))
            lines.append("")
        return "\n".join(lines) + "\n"

    # -- pieces ----------------------------------------------------------

    def _link0(self, job: JobSpec) -> list[str]:
        resources = job.profile.resources
        return [
            f"%chk={job.name}.chk",
            f"%mem={resources.mem}",
            f"%nprocshared={resources.nproc}",
        ]

    def _route(
        self,
        job: JobSpec,
        jobs: tuple[str, ...],
        extra_keywords: tuple[str, ...],
        *,
        method: str | None = None,
        basis: str | None = None,
        dispersion: str | None = "",
    ) -> str:
        profile = job.profile
        method = method if method is not None else profile.method
        basis = basis if basis is not None else profile.basis
        dispersion = profile.dispersion if dispersion == "" else dispersion

        if self._ecp and basis == profile.basis:
            basis = "genecp"

        parts = ["#p"]
        parts += [JOB_KEYWORDS.get(j, j) for j in jobs if j]
        if method and basis:
            parts.append(f"{method}/{basis}")
        elif method:
            parts.append(method)
        if dispersion:
            parts.append(f"empiricaldispersion={dispersion}")
        if profile.solvent:
            model = profile.solvent.model or "iefpcm"
            parts.append(f"scrf=({model},solvent={profile.solvent.name})")
        parts += [k for k in extra_keywords if k]
        return " ".join(parts)

    _ecp: str = ""

    def _ecp_block(self, job: JobSpec, log: IssueLog) -> str:
        """A genecp basis + ECP section when the basis stops short of an element.

        Gaussian stops at the first atom its basis does not define. A Pople
        basis ends at Kr, so any 4d/5d metal needs its own basis and an
        effective core potential: the usual recipe is LANL2DZ on the metal and
        the chosen basis everywhere else.
        """
        heavy = uncovered_elements(job)
        if not heavy:
            return ""
        ecp = job.profile.ecp_basis or "LANL2DZ"
        light = [e for e in dict.fromkeys(job.elements) if e not in heavy]
        report_ecp(job, heavy, ecp, job.profile.basis, log)
        lines = []
        if light:
            lines += [f"{' '.join(light)} 0", job.profile.basis, "****"]
        lines += [f"{' '.join(heavy)} 0", ecp, "****", "", f"{' '.join(heavy)} 0", ecp]
        return "\n".join(lines)

    def _check_gen_basis(self, job: JobSpec, log: IssueLog) -> None:
        basis = (job.profile.basis or "").lower()
        if basis in ("gen", "genecp") and not job.profile.extra_sections:
            log.error(
                "gaussian.gen_without_block",
                f"basis is '{basis}' but the profile defines no extra_sections, so "
                "the file has no basis-set definition and Gaussian will stop at "
                "the route section.",
            )
