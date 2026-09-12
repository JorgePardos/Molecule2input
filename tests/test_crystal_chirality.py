"""Chirality of metal complexes: helicity, planar chirality, racemic crystals.

Every sign convention here is checked against COD structures whose authors
assign the label -- and against their mirror images, because a labelled set
that happens to be all one hand would also be passed by code that always
answers the same.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("gemmi", reason="the CIF route needs the [crystal] extra")

from m2i import config  # noqa: E402
from m2i.crystal import read_cif  # noqa: E402
from m2i.crystal import jobs as crystal_jobs  # noqa: E402
from m2i.crystal.chirality import _has_mirror, mirrored  # noqa: E402
from m2i.crystal.coordination import analyse, skew_angle  # noqa: E402
from m2i.types import IssueLog  # noqa: E402

COD = Path(__file__).parent / "data" / "cod"


def centres(cod: str, *, mirror: bool = False):
    reading = read_cif(COD / f"{cod}.cif", IssueLog())
    species = reading.species[0]
    return reading, analyse(mirrored(species) if mirror else species)


def helicities(found) -> set:
    return {c.helicity for c in found if c.helicity}


# -- Delta / Lambda --------------------------------------------------------------


@pytest.mark.parametrize("cod, label", [
    ("2010630", "Lambda"),  # [Co(en)3]3+
    ("2006181", "Lambda"),  # [Ru(phen)3]2+
    ("2016683", "Delta"),   # [Ru(phen)3]2+
])
def test_helicity_matches_the_authors(cod, label):
    assert helicities(centres(cod)[1]) == {label}


@pytest.mark.parametrize("cod", ["2010630", "2016683"])
def test_the_mirror_image_has_the_other_helicity(cod):
    original = helicities(centres(cod)[1])
    reflected = helicities(centres(cod, mirror=True)[1])
    assert original | reflected == {"Delta", "Lambda"}


def test_all_three_chelate_pairs_agree_in_a_tris_chelate():
    (cobalt,) = centres("2010630")[1]
    assert cobalt.helicity_detail.startswith("3 chelate pair(s)")


def test_the_skew_angle_belongs_to_the_pair_not_to_how_it_was_listed():
    rng = np.random.default_rng(3)
    a1, a2, b1, b2 = rng.normal(size=(4, 3))
    angle = skew_angle(a1, a2, b1, b2)
    assert skew_angle(b1, b2, a1, a2) == pytest.approx(angle)  # order of the lines
    assert skew_angle(a2, a1, b1, b2) == pytest.approx(angle)  # direction of one
    assert skew_angle(a2, a1, b2, b1) == pytest.approx(angle)  # of both


def test_no_helicity_without_chelates():
    (cobalt,) = analyse(read_cif(COD / "2106540.cif", IssueLog()).species[0])  # six NH3
    assert cobalt.helicity is None


# -- planar chirality --------------------------------------------------------------


def test_planar_chirality_matches_the_authors():
    (iron,) = centres("2219943")[1]
    (ring,) = iron.planar
    assert ring.descriptor == "Sp"
    assert ring.certain  # P outranks C at the first atom: no bond order involved
    assert ring.first.startswith("P")


def test_the_mirror_image_is_rp():
    (iron,) = centres("2219943", mirror=True)[1]
    assert [r.descriptor for r in iron.planar] == ["Rp"]


def test_an_unsubstituted_ferrocene_has_no_planar_chirality():
    (iron,) = centres("7033930")[1]
    assert iron.planar == []


@pytest.mark.parametrize("labels, achiral", [
    (["A", "H", "H", "H", "H"], True),   # monosubstituted
    (["A", "A", "H", "H", "H"], True),   # 1,2 with equal groups
    (["A", "H", "A", "H", "H"], True),   # 1,3 with equal groups
    (["A", "B", "H", "H", "H"], False),  # 1,2 with different groups
    (["A", "H", "B", "H", "H"], False),  # 1,3 with different groups
    (["A", "B", "H", "H", "H", "H"], False),  # an arene, 1,2
    (["A", "H", "H", "B", "H", "H"], True),   # an arene, 1,4: a mirror through both
])
def test_which_rings_have_a_mirror(labels, achiral):
    assert _has_mirror(labels) is achiral


def test_the_user_can_settle_an_uncertain_ranking():
    (iron,) = centres("2219943")[1]
    ring = iron.planar[0]
    flipped = ring.flipped()
    assert flipped.descriptor == "Rp" and flipped.first == ring.second
    resolved = crystal_jobs.resolve_planar(iron, {ring.ring_label: ring.second})
    assert resolved[0].descriptor == "Rp"


# -- which enantiomer is in the crystal ------------------------------------------


def test_a_racemic_crystal_says_so():
    reading, found = centres("2017108")
    assert reading.racemic
    notes, chiral = crystal_jobs.chirality_notes(reading, found)
    assert chiral
    assert any("racemic" in note for note in notes)


def test_an_enantiopure_crystal_reports_its_flack_parameter():
    reading, found = centres("2010630")
    assert not reading.racemic
    notes, _ = crystal_jobs.chirality_notes(reading, found)
    assert any("Flack parameter .07(5)" in note for note in notes)


def test_an_achiral_complex_gets_no_chirality_notes():
    reading, found = centres("7033930")
    assert crystal_jobs.chirality_notes(reading, found) == ([], False)


def test_the_mirror_image_can_be_written(tmp_path):
    reading = read_cif(COD / "2010630.cif", IssueLog())
    written = crystal_jobs.write(
        reading, 0, charge=3, multiplicity=1, profile=config.load_profile("xyz_only"),
        output_dir=tmp_path, log=IssueLog(), mirror=True,
    )
    record = json.loads(Path(written[1]).read_text(encoding="utf-8"))
    assert record["chirality"]["mirror_image_written"] is True
    assert any("Delta" in line for line in record["species"]["coordination"])  # Lambda reflected
