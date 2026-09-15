"""Calculation recipes as the output step offers them, and electronic-state choices.

Ports of ``format_settings``, ``_profiles_for`` and ``crystal_page.multiplicity``
from the Streamlit app: which settings a format needs, and building a
validated profile from what the visitor typed.
"""

from __future__ import annotations

from dataclasses import replace

from pydantic import BaseModel

from .. import config
from ..chem import metals
from ..crystal import jobs as crystal_jobs

QM_PROGRAMS = ("gaussian", "orca")
FORMATS = {"gaussian": "Gaussian (.gjf)", "orca": "ORCA (.inp)", "xyz": "XYZ (.xyz)", "sdf": "SDF (.sdf)"}
DISPERSION = {"gaussian": config.GAUSSIAN_DISPERSION, "orca": config.ORCA_DISPERSION}
NO_DISPERSION = "none"


class Settings(BaseModel):
    """The output step's form."""

    format: str = "gaussian"
    recipe: str | None = None
    method: str | None = None
    basis: str | None = None
    dispersion: str | None = None
    solvent: str | None = None
    nproc: int | None = None
    mem: str | None = None


def recipes(program: str) -> list[dict]:
    found = []
    for name in sorted(config.available_profiles()):
        try:
            profile = config.load_profile(name)
        except Exception:  # noqa: BLE001 - a broken local recipe is skipped, not fatal
            continue
        if profile.program != program:
            continue
        found.append({
            "name": name,
            "description": profile.description,
            "method": profile.method,
            "basis": profile.basis,
            "dispersion": (profile.dispersion or NO_DISPERSION).lower(),
            "solvent": profile.solvent.name if profile.solvent else "",
            "jobs": list(profile.jobs),
            "nproc": profile.resources.nproc,
            "mem": profile.resources.mem,
        })
    return found


def dispersion_choices(program: str) -> list[str]:
    return [NO_DISPERSION, *DISPERSION.get(program, ())]


def build(settings: Settings) -> config.JobProfile:
    """A validated profile, or ``config.ProfileError`` with the reason."""
    fmt = settings.format
    if fmt not in FORMATS:
        raise config.ProfileError(f"unknown format {fmt!r}")
    if fmt not in QM_PROGRAMS:
        return config.JobProfile(name=f"{fmt}_geometry", program=fmt, method="", basis="", jobs=())

    names = [r["name"] for r in recipes(fmt)]
    if settings.recipe and settings.recipe in names:
        base = config.load_profile(settings.recipe)
    elif names:
        base = config.load_profile(names[0])
    else:
        base = config.JobProfile(program=fmt)

    def text(value, default):
        return default if value is None else value.strip()

    dispersion = text(settings.dispersion, (base.dispersion or NO_DISPERSION)).lower()
    solvent = text(settings.solvent, base.solvent.name if base.solvent else "")
    profile = replace(
        base,
        method=text(settings.method, base.method),
        basis=text(settings.basis, base.basis),
        dispersion=None if dispersion in ("", NO_DISPERSION) else dispersion,
        solvent=config.Solvent(model=base.solvent.model if base.solvent else "", name=solvent)
        if solvent else None,
        resources=config.Resources(
            mem=text(settings.mem, base.resources.mem) or base.resources.mem,
            nproc=int(settings.nproc or base.resources.nproc),
            maxcore_mb=base.resources.maxcore_mb,
        ),
    )
    profile.validate()
    return profile


def spin_choices(species, charge: int | None, oxidation: dict[str, int]) -> dict:
    """The spin states to offer, as the Streamlit pages did.

    With one metal of known oxidation state, only the states the electron count
    allows; otherwise a free multiplicity, starting from the lowest allowed.
    """
    options = []
    if charge is not None and len(species.metals) == 1 and species.metals[0] in oxidation:
        element = species.metals[0]
        options = [
            {"multiplicity": o.multiplicity, "unpaired": o.unpaired, "label": o.label}
            for o in metals.spin_options(element, oxidation[element])
            if metals.parity_problem(species.elements, charge, o.multiplicity) is None
        ]
    lowest = (
        crystal_jobs.default_multiplicity(species, charge, oxidation) if charge is not None else None
    )
    return {
        "hints": crystal_jobs.spin_hint(species, oxidation),
        "options": options,
        "lowest": lowest,
    }


def parity(species, charge: int, multiplicity: int) -> str | None:
    return metals.parity_problem(species.elements, charge, multiplicity)
