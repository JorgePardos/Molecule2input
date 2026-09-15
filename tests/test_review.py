"""When a structure needs a human look: only when it was read from a picture
and something points at a misreading."""

from __future__ import annotations

from m2i import pipeline
from m2i.report import review
from m2i.types import IssueLog, RecognitionResult


def assess(smiles, backend="decimer", confidence=0.97, molblock=None):
    log = IssueLog()
    reading = RecognitionResult(smiles=smiles, molblock=molblock, confidence=confidence, backend=backend)
    molecule = pipeline.prepare_molecule(reading, log)
    return review.assess(molecule, log)


def test_typed_and_drawn_structures_are_never_asked_about():
    assert not assess("C[C@H](N)C(=O)O", backend="manual").needed
    assert not assess("C[C@H](N)C(=O)O", backend="cdxml").needed


def test_a_clear_reading_without_stereochemistry_is_accepted():
    decision = assess("c1ccccc1O")
    assert not decision.needed
    assert "97%" in decision.summary()


def test_low_confidence_needs_a_look():
    decision = assess("c1ccccc1O", confidence=0.6)
    assert decision.needed and "60%" in decision.summary()


def test_no_confidence_at_all_needs_a_look():
    assert assess("c1ccccc1O", confidence=None).needed


def test_stereochemistry_from_a_picture_always_needs_a_look():
    """The errors a confident DECIMER hides best are lost or invented
    stereodescriptors: at 0.93-0.99 in the calibration."""
    assert assess("C[C@H](N)C(=O)O", confidence=0.99).needed
    assert assess("C/C=C/C(=O)O", confidence=0.99).needed


def test_an_undefined_stereocentre_needs_a_look():
    assert assess("CC(N)Cc1ccccc1", confidence=0.99).needed


def test_what_m2i_warned_about_needs_a_look():
    decision = assess("CCO.[Na+].[Cl-]", confidence=0.99)  # fragments dropped
    assert decision.needed
