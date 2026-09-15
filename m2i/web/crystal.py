"""Crystal structures (.cif): screen 2f. A port of ``gui/crystal_page.py``.

The page's state -- species chosen, hydrogen counts, ring priorities, mirror,
oxidation state, charge -- is sent whole to ``/state`` on every change, and the
answer is everything the page shows. One route to keep consistent instead of
one per widget.
"""

from __future__ import annotations

from pathlib import Path

from fastapi import APIRouter, Depends, UploadFile
from pydantic import BaseModel, Field

from .. import config
from ..crystal import CrystalError, read_cif
from ..crystal import hydrogens as hydrogen_rules
from ..crystal import jobs as crystal_jobs
from ..crystal.coordination import analyse
from ..types import IssueLog
from . import recipes, shapes
from .organic import CIF_SUFFIXES, error, save_upload
from .session import Session, current

router = APIRouter(prefix="/api")


@router.post("/source/cif")
def source_cif(file: UploadFile, normalise: bool = True, session: Session = Depends(current)):
    with session.lock:
        path, digest = save_upload(session, file, CIF_SUFFIXES)
        key = f"cif:{digest}:{int(normalise)}"
        if key not in session.crystals:
            log = IssueLog()
            try:
                reading = read_cif(path, log, normalise_hydrogens=normalise)
            except CrystalError as exc:
                error(str(exc))
            # Stored under a hash: show and name the files by the original.
            reading.path = Path(Path(file.filename or "structure.cif").name)
            session.crystals[key] = (reading, list(log))
        reading, issues = session.crystals[key]
        a, b, c, alpha, beta, gamma = reading.cell
        return {
            "key": key,
            "shown_name": reading.path.name,
            "source": f"COD {reading.cod_id}" if reading.cod_id else reading.path.name,
            "spacegroup": reading.spacegroup,
            "z": reading.z,
            "cell": [a, b, c, alpha, beta, gamma],
            "normalise": normalise,
            "species": [
                {"index": i, "formula": s.formula, "copies": s.copies, "atoms": s.n_atoms,
                 "metals": sorted(set(s.metals))}
                for i, s in enumerate(reading.species)
            ],
            "extended": list(reading.extended),
            "issues": shapes.issues(issues),
        }


class CrystalState(BaseModel):
    species: int = 0
    counts: dict[str, int] = Field(default_factory=dict)
    add_hydrogens: bool | None = None  # None: add them when they match the formula
    allow_missing: bool = False
    priorities: dict[str, str] = Field(default_factory=dict)
    mirror: bool = False
    oxidation: dict[str, int] = Field(default_factory=dict)
    charge: int | None = None


def _reading(session: Session, key: str):
    entry = session.crystals.get(key)
    if entry is None:
        error("That crystal is no longer in this session; upload it again.", 404)
    return entry


def _resolve(session: Session, key: str, state: CrystalState):
    """The species that will be written, and everything decided on the way."""
    reading, issues = _reading(session, key)
    if not 0 <= state.species < len(reading.species):
        error("There is no such species in this crystal.")
    index = state.species
    original = reading.species[index]
    hydrogens = {"needed": False, "elsewhere": False}
    species, added, allow_missing = original, None, False
    if reading.missing_hydrogens:
        if hydrogen_rules.lacking(reading, index):
            proposal = hydrogen_rules.plan(reading, index)
            defaults = {s.label: s.count for s in proposal.sites}
            changed = {k: v for k, v in state.counts.items() if k in defaults and v != defaults[k]}
            final = proposal.with_counts(changed)
            add = final.matches if state.add_hydrogens is None else state.add_hydrogens
            hydrogens = {
                "needed": True,
                "sites": [
                    {"label": s.label, "count": s.count, "reason": s.reason, "certain": s.certain,
                     "proposed": defaults[s.label]}
                    for s in final.sites
                ],
                "total": final.total, "target": final.target, "matches": final.matches, "add": add,
            }
            if add:
                cache_key = ("completed", key, index, tuple(sorted((s.label, s.count) for s in final.sites)))
                if cache_key not in session.crystals:
                    session.crystals[cache_key] = hydrogen_rules.apply(original, final)
                species, added = session.crystals[cache_key]
            else:
                allow_missing = state.allow_missing
        else:
            hydrogens = {"needed": False, "elsewhere": True}
            allow_missing = True
    return reading, issues, index, original, species, added, allow_missing, hydrogens


@router.post("/crystal/{key}/state")
def crystal_state(key: str, state: CrystalState, session: Session = Depends(current)):
    with session.lock:
        reading, issues, index, original, species, added, allow_missing, hydrogens = _resolve(session, key, state)
        centres = analyse(original)
        questions = []
        for centre in centres:
            for ring in centre.planar:
                if not ring.certain:
                    questions.append({"ring": ring.ring_label, "reason": ring.reason,
                                      "options": [ring.first, ring.second],
                                      "chosen": state.priorities.get(ring.ring_label, ring.first)})
        priorities = {q["ring"]: q["chosen"] for q in questions}
        notes, chiral = crystal_jobs.chirality_notes(reading, centres, priorities)
        asked = crystal_jobs.questions(reading, index)
        charge = state.charge if state.charge is not None else asked.suggested_charge
        blockers = []
        if hydrogens.get("needed") and not hydrogens.get("add") and not allow_missing:
            blockers.append("Hydrogens are missing from the CIF (see above); nothing will be written.")
        return {
            "species": {"index": index, "formula": species.formula, "atoms": species.n_atoms,
                        "metals": sorted(set(species.metals)), "completed": added is not None},
            "hydrogens": hydrogens,
            "viewer": shapes.viewer(species),
            "coordination": shapes.coordination(original),
            "chirality": {"chiral": chiral, "notes": notes, "questions": questions},
            "electronic": {
                "metals": asked.metals,
                "needs_answers": asked.needs_answers,
                "suggested_charge": asked.suggested_charge,
                "charge_reason": asked.charge_reason,
                "charge": charge,
                **recipes.spin_choices(species, charge, state.oxidation),
            },
            "blockers": blockers,
            "issues": shapes.issues(issues),
        }


class CrystalGenerate(CrystalState):
    multiplicity: int
    settings: recipes.Settings


@router.post("/crystal/{key}/generate")
def crystal_generate(key: str, body: CrystalGenerate, session: Session = Depends(current)):
    with session.lock:
        reading, issues, index, original, species, added, allow_missing, hydrogens = _resolve(session, key, body)
        if body.settings.format == "sdf":
            error("A CIF has no bond orders, so SDF is not offered: choose Gaussian, ORCA or XYZ.")
        if body.charge is None:
            error("Set the charge to continue.")
        if hydrogens.get("needed") and not hydrogens.get("add") and not allow_missing:
            error("Hydrogens are missing from the CIF (see above); nothing will be written.")
        problem = recipes.parity(species, body.charge, body.multiplicity)
        if problem:
            error(f"Impossible electronic state: {problem}.")
        try:
            profile = recipes.build(body.settings)
        except config.ProfileError as exc:
            error(str(exc))
        log = IssueLog(issues)
        try:
            written = crystal_jobs.write(
                reading, index, charge=body.charge, multiplicity=body.multiplicity,
                profile=profile, output_dir=session.output, log=log, oxidation=body.oxidation,
                allow_missing_hydrogens=allow_missing,
                species=species if added is not None else None, hydrogens_added=added,
                mirror=body.mirror, priorities=body.priorities,
            )
        except CrystalError as exc:
            error(str(exc))
        return _written(session, written, body.settings, profile, log, issues)


def _written(session, written, settings, profile, log, earlier) -> dict:
    """The result screen's payload for a crystal or a complex: one input, one record."""
    files, extras = [], []
    for path in map(Path, written):
        entry = session.offer([path])[0]
        if path.name.endswith(".m2i.json"):
            entry["kind"] = "provenance"
            extras.append(entry)
        else:
            entry.update({"kind": "input", "preview": path.read_text(encoding="utf-8", errors="replace"),
                          "conformer": None})
            files.append(entry)
    return {
        "format": settings.format,
        "format_label": recipes.FORMATS[settings.format],
        "jobs": " ".join(profile.jobs) or "single point",
        "files": files,
        "extras": extras,
        "issues": shapes.new_issues(log, earlier),
    }
