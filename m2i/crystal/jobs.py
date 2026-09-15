"""From a species of the crystal to input files, shared by the CLI and the GUI.

The geometry is written exactly as rebuilt from the crystal. Everything that
could be wrong with it is checked before a file exists: hydrogens that were
never located, and a charge/multiplicity pair the electron count forbids.
"""

from __future__ import annotations

import platform
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from .. import __version__
from ..chem import metals
from ..report import validate
from ..types import IssueLog, JobSpec
from ..writers import get_writer
from .chirality import mirrored
from .coordination import MetalCentre, analyse
from .structure import CrystalError, CrystalReading, Species

SAFE = re.compile(r"[^A-Za-z0-9._-]+")


@dataclass
class ElectronicQuestions:
    """What to ask the user about a species, with m2i's best suggestions."""

    metals: list[str]
    suggested_charge: int | None
    charge_reason: str
    centres: list[MetalCentre] = field(default_factory=list)

    @property
    def needs_answers(self) -> bool:
        """A metal means the charge and spin cannot be taken on trust."""
        return bool(self.metals)


def questions(reading: CrystalReading, index: int) -> ElectronicQuestions:
    species = reading.species[index]
    others = [(s.formula, s.copies) for i, s in enumerate(reading.species) if i != index]
    charge, reason = metals.suggest_charge(species.formula, species.copies, others)
    return ElectronicQuestions(
        metals=sorted(set(species.metals)),
        suggested_charge=charge,
        charge_reason=reason,
        centres=analyse(species),
    )


def spin_hint(species: Species, oxidation: dict[str, int]) -> list[str]:
    """Human-readable spin options for each metal whose oxidation state is known."""
    lines = []
    for element in sorted(set(species.metals)):
        if element not in oxidation:
            continue
        state = oxidation[element]
        shell = metals.valence_electron_count(element, state)
        options = metals.spin_options(element, state)
        if shell is None or not options:
            lines.append(f"{element}({state:+d}): no d/f-shell guidance for this element")
            continue
        kind, n = shell
        choices = ", ".join(f"{o.multiplicity} ({o.label})" for o in options)
        lines.append(f"{element}({state:+d}) is {kind}{n}: multiplicity {choices}")
    if len(species.metals) > 1:
        lines.append(
            "Several metal centres: the total multiplicity depends on how their "
            "spins couple (ferro- or antiferromagnetically), which m2i cannot know."
        )
    return lines


def default_multiplicity(species: Species, charge: int, oxidation: dict[str, int]) -> int:
    """The lowest multiplicity the electron count allows -- offered, never imposed."""
    if len(set(species.metals)) == 1 and len(species.metals) == 1:
        element = species.metals[0]
        if element in oxidation:
            for option in metals.spin_options(element, oxidation[element]):
                if metals.parity_problem(species.elements, charge, option.multiplicity) is None:
                    return option.multiplicity
    electrons = metals.electron_count(species.elements, charge)
    return 1 if electrons % 2 == 0 else 2


def chirality_notes(
    reading: CrystalReading,
    centres: list[MetalCentre],
    priorities: dict[str, str] | None = None,
) -> tuple[list[str], bool]:
    """What there is to say about chirality, and whether the species is chiral.

    ``priorities`` maps a ring label to the substituent label the user ranked
    first, for rings whose CIP ranking m2i could not settle.
    """
    notes, chiral = [], False
    for centre in centres:
        if centre.helicity in ("Delta", "Lambda"):
            chiral = True
            notes.append(f"{centre.label}: {centre.helicity} ({centre.helicity_detail})")
        elif centre.helicity == "mixed":
            notes.append(f"{centre.label}: chelate helicities disagree ({centre.helicity_detail})")
        for ring in resolve_planar(centre, priorities):
            chiral = True
            notes.append(f"{centre.label}: {ring.describe()}")
    if not chiral:
        return notes, False

    if reading.racemic:
        notes.append(
            "The crystal is racemic (its space group has an inversion, mirror or "
            "glide): both enantiomers are present and this is one of them. Energies "
            "are the same for both; for chiroptical properties choose the one you "
            "need with the mirror-image option."
        )
    elif reading.flack:
        notes.append(
            f"Enantiopure crystal; absolute configuration from the refinement, Flack "
            f"parameter {reading.flack} (reliable when close to 0 with a small "
            "uncertainty; near 1 means the model is inverted)."
        )
    else:
        notes.append(
            "Enantiopure crystal, but the CIF gives no Flack parameter: the absolute "
            "configuration may not have been determined."
        )
    return notes, True


def resolve_planar(centre: MetalCentre, priorities: dict[str, str] | None):
    """Planar chirality of a centre, with any ranking the user decided applied."""
    resolved = []
    for ring in centre.planar:
        chosen = (priorities or {}).get(ring.ring_label)
        if chosen and chosen == ring.second:
            ring = ring.flipped()
        elif chosen == ring.first and not ring.certain:
            ring = type(ring)(**{**ring.__dict__, "certain": True, "reason": "priority set by the user"})
        resolved.append(ring)
    return resolved


def write(
    reading: CrystalReading,
    index: int,
    *,
    charge: int,
    multiplicity: int,
    profile,
    output_dir: Path,
    log: IssueLog,
    oxidation: dict[str, int] | None = None,
    name: str | None = None,
    allow_missing_hydrogens: bool = False,
    mirror: bool = False,
    priorities: dict[str, str] | None = None,
    species: Species | None = None,
    hydrogens_added: list[str] | None = None,
) -> list[str]:
    """Write the input for one species, after the checks that must pass first.

    ``species`` is the species with its missing hydrogens put back, when that
    was done; it replaces the one read from the crystal.
    """
    completed = species is not None
    species = species if completed else reading.species[index]
    if mirror:
        species = mirrored(species)

    from .hydrogens import lacking

    if not completed and not allow_missing_hydrogens and lacking(reading, index):
        raise CrystalError(
            f"{reading.missing_hydrogens:g} hydrogen(s) per formula unit are missing "
            "from the CIF, so this would be a different molecule. Add them first, "
            "or pass --allow-missing-hydrogens if you really mean it."
        )
    if log.has_errors() and any(i.code == "crystal.overlap" for i in log):
        raise CrystalError("overlapping atoms (unresolved disorder); see the errors above")
    problem = metals.parity_problem(species.elements, charge, multiplicity)
    if problem:
        raise CrystalError(f"impossible electronic state: {problem}")

    stem = SAFE.sub("_", name or f"{species.formula}_{reading.cod_id or reading.path.stem}")
    source = f"COD {reading.cod_id}" if reading.cod_id else reading.path.name
    job = JobSpec(
        name=stem,
        title=(
            f"{species.formula} | experimental geometry from {source} "
            f"({reading.spacegroup}) | charge {charge} mult {multiplicity} | m2i {__version__}"
        ),
        elements=species.elements,
        coords=centred(species).tolist(),
        charge=charge,
        multiplicity=multiplicity,
        profile=profile,
        metadata={"formula": species.formula, "source": source},
    )
    output_dir = Path(output_dir)
    written = [str(get_writer(profile.program).write(job, output_dir / stem, log))]

    record = provenance(reading, index, charge, multiplicity, oxidation or {}, profile, log)
    # Describe what was written, which is the mirror image when asked for.
    centres = analyse(species)
    record["species"].update(
        formula=species.formula,
        atoms=species.n_atoms,
        labels=species.labels,
        coordination=[line for c in centres for line in c.describe()],
    )
    record["reconstruction"]["hydrogens_added"] = hydrogens_added or []
    record["chirality"] = {
        "mirror_image_written": mirror,
        "racemic_crystal": reading.racemic,
        "flack": reading.flack,
        "notes": chirality_notes(reading, centres, priorities)[0],
        "priorities_set_by_user": priorities or {},
    }
    record["written_files"] = list(written)
    path = validate.write_provenance(output_dir / f"{stem}.m2i.json", record)
    written.append(str(path))
    return written


def centred(species: Species):
    """Coordinates moved so the metal (or the centroid) sits at the origin.

    A pure translation: it changes nothing physical, but a molecule left where
    it sat in the unit cell reads badly and hides the metal among offsets.
    """
    metals_at = [i for i, e in enumerate(species.elements) if e in species.metals]
    origin = (
        species.coords[metals_at[0]] if len(metals_at) == 1 else species.coords.mean(axis=0)
    )
    return species.coords - origin


def provenance(reading, index, charge, multiplicity, oxidation, profile, log) -> dict:
    species = reading.species[index]
    return {
        "m2i_version": __version__,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python_version": platform.python_version(),
        "source": {
            "file": reading.path.name,
            "cod_id": reading.cod_id,
            "block": reading.block,
            "spacegroup": reading.spacegroup,
            "cell": list(reading.cell),
            "z": reading.z,
            "formula_sum": reading.formula_sum,
        },
        "reconstruction": {
            "disorder": reading.disorder,
            "hydrogens_normalised": reading.hydrogens_normalised,
            "missing_hydrogens_per_formula_unit": reading.missing_hydrogens,
            "extended_networks": reading.extended,
            "species_in_cell": [
                {"formula": s.formula, "copies": s.copies} for s in reading.species
            ],
        },
        "species": {
            "formula": species.formula,
            "atoms": species.n_atoms,
            "labels": species.labels,
            "coordination": [line for c in analyse(species) for line in c.describe()],
        },
        "electronic_state": {
            "charge": charge,
            "multiplicity": multiplicity,
            "oxidation_states": oxidation,
            "decided_by": "user",
        },
        "profile": profile.to_dict(),
        "issues": log.to_list(),
    }
