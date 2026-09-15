"""Whether a structure needs a human look before any input is written.

Only a picture can be misread. A SMILES, a ChemDraw file or a molfile says
exactly what it contains, so there is nothing to confirm -- whatever m2i finds
wrong with it (an undefined stereocentre, a valence it had to repair) is still
reported, and errors still stop it, but nobody is asked "is this what you
meant?" about something they typed or drew in ChemDraw themselves.

A reading of a picture is accepted without a look only when nothing points at
a mistake:

- the model's own confidence below ``CONFIDENT``;
- any stereochemistry at all. DECIMER writes stereocentres and double-bond
  geometry as tokens instead of measuring them off the drawing;
- anything m2i itself warned about while reading it: a repaired valence, a
  dropped fragment, an undefined stereocentre.

Why all three, and not the confidence alone: 60 molecules rendered by RDKit,
clean and degraded like a photo (rotation, blur, noise, JPEG), read by DECIMER
through m2i. 113 of 120 readings were exact. Of the 7 wrong ones, 4 came back
at 0.93-0.99 confidence -- a penicillin with its stereocentres wrong, a
resveratrol that lost its E, an amphetamine with a stereocentre invented --
so no threshold catches them; the stereochemistry rule does. With all three
rules, none of the 7 would have been accepted unseen, and 70 of the 113 right
readings would have gone through without a question. Rendered depictions are
not hand drawings, and real photos will be misread more often; the rules are
meant to stay on the safe side of that.

Accepting is recorded in the provenance with its reasons, like confirming.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..types import IssueLog, MoleculeSpec

#: A picture reading at or above this confidence, with nothing else against it,
#: is accepted without asking.
CONFIDENT = 0.9


@dataclass
class Review:
    needed: bool
    #: why a look is needed -- or, when it is not, why not
    reasons: list[str] = field(default_factory=list)

    def summary(self) -> str:
        return "; ".join(self.reasons)


def assess(molecule: MoleculeSpec, log: IssueLog) -> Review:
    """Decide from the structure, where it came from, and what was noticed."""
    source = molecule.source
    if source is None or not is_picture_reading(source.backend):
        how = "typed in" if source is None or source.backend == "manual" else "read from the file"
        return Review(False, [f"the structure was {how}, not recognised from a picture"])

    reasons = []
    if source.confidence is None:
        reasons.append(f"{source.backend} gave no confidence for this reading")
    elif source.confidence < CONFIDENT:
        reasons.append(f"{source.backend} is only {_percent(source.confidence)} confident")
    stereo = molecule.stereo
    if stereo.centers or stereo.bonds:
        reasons.append(
            f"its stereochemistry ({stereo.describe()}) was written by the model, "
            "not measured off the drawing"
        )
    for issue in log.warnings:
        if issue.code != "recognition.low_confidence":  # already said above
            reasons.append(issue.message)

    if reasons:
        return Review(True, reasons)
    return Review(False, [
        f"{source.backend} is {_percent(source.confidence)} confident, the molecule has no "
        "stereochemistry, and nothing looked wrong"
    ])


def _percent(value: float) -> str:
    """Truncated, never rounded up: 0.996 is not 100% confident."""
    return f"{int(value * 100)}%"


def is_picture_reading(backend: str) -> bool:
    from ..backends import SPECS

    return backend in SPECS


def record(log: IssueLog, review: Review, *, confirmed: bool | None) -> None:
    """Keep the decision with the other notes, so the provenance says it."""
    if not review.needed:
        log.info("review.not_needed", f"Accepted without review: {review.summary()}.")
    elif confirmed:
        log.info("review.confirmed", f"Confirmed by the user after review: {review.summary()}.")
    else:
        log.info("review.waived", f"Review waived (--yes or batch): {review.summary()}.")
