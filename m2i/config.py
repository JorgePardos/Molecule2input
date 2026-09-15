"""Job profiles: method/basis/solvent/resources, loaded from YAML.

A profile is everything the writers need that does not depend on the molecule.
Built-in profiles live in ``m2i/profiles``; the user can add their own in
``./profiles`` or ``~/.m2i/profiles`` and they take precedence by name.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

import yaml

PROGRAMS = ("gaussian", "orca", "xyz", "sdf")

# Solvent-model keywords accepted per program. Used for validation only; the
# solvent *name* is passed through verbatim because each program keeps its own
# (long, frequently updated) table of solvents.
SOLVENT_MODELS = {
    "gaussian": ("pcm", "cpcm", "iefpcm", "smd", "dipole"),
    "orca": ("cpcm", "smd"),
}

# Gaussian dispersion keywords for EmpiricalDispersion=
GAUSSIAN_DISPERSION = ("gd3", "gd3bj", "gd2")

# ORCA dispersion keywords appended to the ! line
ORCA_DISPERSION = ("d3", "d3bj", "d3zero", "d4")


class ProfileError(ValueError):
    """Raised when a profile is malformed or names an unknown option."""


@dataclass(frozen=True)
class Solvent:
    model: str
    name: str

    @classmethod
    def from_obj(cls, obj: Any) -> "Solvent | None":
        if obj is None:
            return None
        if isinstance(obj, str):  # bare solvent name, model defaults per program
            return cls(model="", name=obj)
        if isinstance(obj, dict):
            name = obj.get("name")
            if not name:
                raise ProfileError("solvent needs a 'name'")
            return cls(model=str(obj.get("model", "")).lower(), name=str(name))
        raise ProfileError(f"cannot read solvent from {obj!r}")


@dataclass(frozen=True)
class Resources:
    mem: str = "4GB"
    nproc: int = 4
    # ORCA needs per-core memory; derived from mem/nproc when not given.
    maxcore_mb: int | None = None

    def orca_maxcore(self) -> int:
        if self.maxcore_mb:
            return self.maxcore_mb
        return max(256, int(_mem_to_mb(self.mem) * 0.75) // max(1, self.nproc))


@dataclass(frozen=True)
class JobProfile:
    """A reusable calculation recipe."""

    name: str = "custom"
    program: str = "gaussian"
    method: str = "b3lyp"
    basis: str = "6-31G(d)"
    dispersion: str | None = None
    solvent: Solvent | None = None
    jobs: tuple[str, ...] = ("opt", "freq")
    resources: Resources = field(default_factory=Resources)
    # Free-form extras merged verbatim into the route/keyword line.
    extra_keywords: tuple[str, ...] = ()
    # Raw text blocks appended after the coordinates (gen basis, ModRedundant...).
    extra_sections: tuple[str, ...] = ()
    # Follow-up jobs chained with --Link1-- (Gaussian) or written as extra files.
    link_jobs: tuple[dict, ...] = ()
    description: str = ""
    # Basis (with its ECP) for elements the main basis does not cover, such as
    # a metal beyond Kr with a Pople basis. None = the writer's default:
    # LANL2DZ in Gaussian, def2-TZVP with def2-ECP in ORCA.
    ecp_basis: str | None = None

    def validate(self) -> None:
        if self.program not in PROGRAMS:
            raise ProfileError(
                f"unknown program {self.program!r}; expected one of {', '.join(PROGRAMS)}"
            )
        if self.program in ("gaussian", "orca"):
            if not self.method:
                raise ProfileError("profile needs a 'method'")
            if not self.basis:
                raise ProfileError("profile needs a 'basis'")
        if self.dispersion:
            known = (
                GAUSSIAN_DISPERSION if self.program == "gaussian" else ORCA_DISPERSION
            )
            if self.program in ("gaussian", "orca") and self.dispersion.lower() not in known:
                raise ProfileError(
                    f"dispersion {self.dispersion!r} not recognised for {self.program}; "
                    f"expected one of {', '.join(known)}"
                )
        if self.solvent and self.solvent.model:
            known_models = SOLVENT_MODELS.get(self.program)
            if known_models and self.solvent.model not in known_models:
                raise ProfileError(
                    f"solvent model {self.solvent.model!r} not recognised for "
                    f"{self.program}; expected one of {', '.join(known_models)}"
                )
        if self.resources.nproc < 1:
            raise ProfileError("resources.nproc must be >= 1")

    # -- construction ----------------------------------------------------

    @classmethod
    def from_dict(cls, data: dict, *, name: str | None = None) -> "JobProfile":
        data = dict(data)
        unknown = set(data) - {
            "name",
            "program",
            "method",
            "basis",
            "dispersion",
            "solvent",
            "jobs",
            "resources",
            "extra_keywords",
            "extra_sections",
            "link_jobs",
            "description",
            "ecp_basis",
        }
        if unknown:
            raise ProfileError(f"unknown profile keys: {', '.join(sorted(unknown))}")

        res = data.get("resources") or {}
        if not isinstance(res, dict):
            raise ProfileError("'resources' must be a mapping")
        resources = Resources(
            mem=str(res.get("mem", "4GB")),
            nproc=int(res.get("nproc", 4)),
            maxcore_mb=int(res["maxcore_mb"]) if res.get("maxcore_mb") else None,
        )

        jobs = data.get("jobs", ("opt", "freq"))
        if isinstance(jobs, str):
            jobs = jobs.split()

        profile = cls(
            name=str(data.get("name", name or "custom")),
            program=str(data.get("program", "gaussian")).lower(),
            method=str(data.get("method", "b3lyp")),
            basis=str(data.get("basis", "6-31G(d)")),
            dispersion=(str(data["dispersion"]) if data.get("dispersion") else None),
            solvent=Solvent.from_obj(data.get("solvent")),
            jobs=tuple(str(j).lower() for j in jobs),
            resources=resources,
            extra_keywords=tuple(_as_list(data.get("extra_keywords"))),
            extra_sections=tuple(_as_list(data.get("extra_sections"))),
            link_jobs=tuple(data.get("link_jobs") or ()),
            description=str(data.get("description", "")),
            ecp_basis=(str(data["ecp_basis"]) if data.get("ecp_basis") else None),
        )
        profile.validate()
        return profile

    def merged_with(self, **overrides: Any) -> "JobProfile":
        """Return a copy with the non-None overrides applied (CLI flags win)."""
        clean = {k: v for k, v in overrides.items() if v is not None}
        if not clean:
            return self
        if "solvent" in clean and not isinstance(clean["solvent"], (Solvent, type(None))):
            clean["solvent"] = Solvent.from_obj(clean["solvent"])

        new_program = clean.get("program", self.program)
        if new_program != self.program and "dispersion" not in clean and self.dispersion:
            # Dispersion keywords are program-specific: --program orca on a
            # Gaussian profile would otherwise fail validation on 'gd3bj'.
            clean["dispersion"] = translate_dispersion(
                self.dispersion, self.program, new_program
            )

        profile = replace(self, **clean)
        profile.validate()
        return profile

    def to_dict(self) -> dict:
        d = {
            "name": self.name,
            "program": self.program,
            "method": self.method,
            "basis": self.basis,
            "jobs": list(self.jobs),
            "resources": {"mem": self.resources.mem, "nproc": self.resources.nproc},
        }
        if self.dispersion:
            d["dispersion"] = self.dispersion
        if self.solvent:
            d["solvent"] = {"model": self.solvent.model, "name": self.solvent.name}
        if self.extra_keywords:
            d["extra_keywords"] = list(self.extra_keywords)
        if self.description:
            d["description"] = self.description
        if self.ecp_basis:
            d["ecp_basis"] = self.ecp_basis
        return d


# -- dispersion keywords across programs ---------------------------------

# Same physics, different spelling. Only pairs that mean the same correction.
DISPERSION_EQUIVALENTS = {
    ("gaussian", "orca"): {"gd3bj": "d3bj", "gd3": "d3zero"},
    ("orca", "gaussian"): {"d3bj": "gd3bj", "d3zero": "gd3", "d3": "gd3"},
}


def translate_dispersion(value: str, from_program: str, to_program: str) -> str | None:
    """Rewrite a dispersion keyword when the target program is changed.

    Returns None when the target program needs no keyword (xyz/sdf output) and
    raises when the correction has no equivalent, because silently dropping a
    dispersion correction changes the chemistry.
    """
    if to_program not in ("gaussian", "orca"):
        return None
    if from_program == to_program:
        return value

    table = DISPERSION_EQUIVALENTS.get((from_program, to_program), {})
    translated = table.get(value.lower())
    if translated:
        return translated
    raise ProfileError(
        f"dispersion {value!r} has no {to_program} equivalent; pick a profile for "
        f"{to_program} or set the dispersion explicitly in a profile"
    )


# -- profile discovery ---------------------------------------------------

BUILTIN_PROFILE_DIR = Path(__file__).parent / "profiles"


def profile_search_path() -> list[Path]:
    """Directories searched for profiles, most specific first."""
    paths: list[Path] = []
    env = os.environ.get("M2I_PROFILE_PATH")
    if env:
        paths += [Path(p) for p in env.split(os.pathsep) if p]
    paths.append(Path.cwd() / "profiles")
    paths.append(Path.home() / ".m2i" / "profiles")
    paths.append(BUILTIN_PROFILE_DIR)
    return paths


def available_profiles() -> dict[str, Path]:
    """Map profile name -> file, with earlier search paths shadowing later ones."""
    found: dict[str, Path] = {}
    for directory in profile_search_path():
        if not directory.is_dir():
            continue
        for path in sorted(directory.glob("*.yaml")) + sorted(directory.glob("*.yml")):
            found.setdefault(path.stem, path)
    return found


def load_profile(name_or_path: str) -> JobProfile:
    """Load a profile by name (searched) or by explicit file path."""
    candidate = Path(name_or_path)
    if candidate.suffix in (".yaml", ".yml") and candidate.is_file():
        return _load_profile_file(candidate)

    profiles = available_profiles()
    if name_or_path in profiles:
        return _load_profile_file(profiles[name_or_path])

    known = ", ".join(sorted(profiles)) or "(none found)"
    raise ProfileError(f"profile {name_or_path!r} not found. Available: {known}")


def _load_profile_file(path: Path) -> JobProfile:
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except yaml.YAMLError as exc:
        raise ProfileError(f"{path}: invalid YAML: {exc}") from exc
    if not isinstance(data, dict):
        raise ProfileError(f"{path}: expected a mapping at the top level")
    try:
        return JobProfile.from_dict(data, name=path.stem)
    except ProfileError as exc:
        raise ProfileError(f"{path}: {exc}") from exc


# -- helpers -------------------------------------------------------------


def _as_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    return [str(v) for v in value]


def _mem_to_mb(mem: str) -> int:
    """Parse '8GB' / '8000MB' / '8' (GB assumed) into megabytes."""
    text = str(mem).strip().lower().replace(" ", "")
    for suffix, factor in (("tb", 1024 * 1024), ("gb", 1024), ("mb", 1), ("kb", 0)):
        if text.endswith(suffix):
            number = text[: -len(suffix)]
            try:
                return max(1, int(float(number) * factor))
            except ValueError:
                break
    try:
        return int(float(text) * 1024)  # bare number: gigabytes
    except ValueError:
        return 4096
