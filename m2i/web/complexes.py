"""Metal complexes drawn in ChemDraw: screen 2e. A port of ``gui/organometallic_page.py``.

As for crystals, the whole page state goes to ``/state`` and the answer is the
whole page: what the drawing leaves out, the arrangements, the built complex
and the electronic-state choices.
"""

from __future__ import annotations

from pathlib import Path

from fastapi import APIRouter, Depends
from pydantic import BaseModel, Field

from .. import config
from ..organometallic import arrangements, build, mirror_partners, read
from ..organometallic import jobs as om_jobs
from ..recognition.base import BackendError
from ..types import IssueLog
from . import recipes, shapes
from .crystal import _written
from .organic import chosen_name, error
from .session import Session, current

router = APIRouter(prefix="/api")

#: Arrangements shown; beyond these the drawing clearly does not fit.
SHOWN = 12


class ComplexState(BaseModel):
    substituents: dict[str, str] = Field(default_factory=dict)
    arrangement: int = 0
    seed: int = 0xF00D
    oxidation: dict[str, int] = Field(default_factory=dict)
    charge: int | None = None


def _drawing(session: Session, key: str, substituents: dict[str, str]):
    source = session.sources.get(key)
    if source is None or source.kind != "om":
        error("That drawing is no longer in this session; upload it again.", 404)
    cache_key = (key, tuple(sorted(substituents.items())))
    if cache_key not in session.om_drawings:
        try:
            session.om_drawings[cache_key] = read(source.path, substituents)
        except BackendError as exc:
            error(str(exc))
    return source, session.om_drawings[cache_key]


def _gaps(session: Session, key: str) -> dict:
    """What the drawing leaves out, from the drawing as it is, without additions."""
    _, base = _drawing(session, key, {})
    return {**base.assumed, **base.completed}


def _resolve(session: Session, key: str, state: ComplexState):
    substituents = {k: v.strip() for k, v in state.substituents.items() if v and v.strip()}
    source, drawing = _drawing(session, key, substituents)
    options = arrangements(drawing.donors, drawing.drawn, drawing.chelated, labels=drawing.labels)
    if not options:
        error(
            f"{len(drawing.donors)} atoms are bound to {drawing.metal_symbol}; complexes with 2 to 6 "
            "donors can be built from a drawing. Is one of those bonds drawn by mistake? A bare line "
            "ending at the metal is read as a methyl."
        )
    choice = state.arrangement if 0 <= state.arrangement < min(len(options), SHOWN) else 0
    build_key = (key, tuple(sorted(substituents.items())), choice, state.seed)
    if build_key not in session.om_built:
        try:
            session.om_built[build_key] = build(drawing, options[choice], seed=state.seed)
        except BackendError as exc:
            error(str(exc))
    return source, drawing, options, choice, session.om_built[build_key]


@router.post("/om/{key}/state")
def complex_state(key: str, state: ComplexState, session: Session = Depends(current)):
    with session.lock:
        gaps = _gaps(session, key)
        source, drawing, options, choice, built = _resolve(session, key, state)
        mirrors = mirror_partners(options, drawing.labels, drawing.chelated)
        species = built.species
        charge = state.charge if state.charge is not None else drawing.drawn_charge
        return {
            "shown_name": source.shown_name,
            "drawing": {
                "metal": drawing.metal_symbol,
                "donors": [drawing.names[d] for d in drawing.donors],
                "drawn_charge": drawing.drawn_charge,
                "notes": [{"level": level, "code": code, "message": message}
                          for level, code, message in drawing.notes],
                "gaps": [
                    {"atom": atom, "bonds_missing": len(added), "groups": state.substituents.get(atom, "")}
                    for atom, added in gaps.items()
                ],
            },
            "arrangements": [
                {"index": i, "label": om_jobs.describe(option, drawing), "fit": round(option.misfit, 2),
                 "mirror_of": mirrors.get(i)}
                for i, option in enumerate(options[:SHOWN])
            ],
            "ambiguous": om_jobs.ambiguous(options),
            "mirror_only": om_jobs.mirror_only(options, drawing),
            "arrangement": choice,
            "built": {
                "formula": species.formula,
                "name": om_jobs.default_name(species, Path(source.shown_name or source.path.name)),
                "viewer": shapes.viewer(species),
                "coordination": shapes.coordination(species),
                "notes": [message for _, _, message in built.notes],
            },
            "electronic": {
                "metals": [drawing.metal_symbol],
                "charge": charge,
                "charge_reason": f"the formal charges drawn add up to {drawing.drawn_charge:+d}; a "
                                 "drawing often leaves charges out, so check it",
                **recipes.spin_choices(species, charge, state.oxidation),
            },
        }


class ComplexGenerate(ComplexState):
    multiplicity: int
    settings: recipes.Settings
    name: str | None = None  # blank: named after the formula and the drawing


@router.post("/om/{key}/generate")
def complex_generate(key: str, body: ComplexGenerate, session: Session = Depends(current)):
    with session.lock:
        source, drawing, options, choice, built = _resolve(session, key, body)
        if body.charge is None:
            error("Set the charge to continue.")
        problem = recipes.parity(built.species, body.charge, body.multiplicity)
        if problem:
            error(f"Impossible electronic state: {problem}.")
        try:
            profile = recipes.build(body.settings)
        except config.ProfileError as exc:
            error(str(exc))
        log = IssueLog()
        try:
            written = om_jobs.write(
                built, drawing, source=Path(source.shown_name or source.path.name), charge=body.charge,
                multiplicity=body.multiplicity, profile=profile, output_dir=session.output, log=log,
                oxidation=body.oxidation, name=chosen_name(body.name),
            )
        except BackendError as exc:
            error(str(exc))
        return _written(session, written, body.settings, profile, log, [])
