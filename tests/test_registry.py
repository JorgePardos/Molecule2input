"""Backend selection and reporting, with stand-in backends."""

from __future__ import annotations

import pytest
from conftest import codes

from m2i.recognition import registry
from m2i.recognition.base import BackendError
from m2i.types import RecognitionResult

S_ENANTIOMER = "F[C@@H](Cl)Br"


class FakeBackend:
    def __init__(self, name, smiles, *, molblock=None, confidence=None, fails=False):
        self.name = name
        self.description = f"fake {name}"
        self.strength = f"whatever {name} is good at"
        self.returns_molblock = molblock is not None
        self._result = RecognitionResult(
            smiles=smiles, molblock=molblock, confidence=confidence, backend=name
        )
        self._fails = fails

    def available(self):
        return True, "fake"

    def recognize(self, image_path, **options):
        if self._fails:
            raise RuntimeError("model exploded")
        return self._result


def molblock_for(smiles: str) -> str:
    from rdkit import Chem
    from rdkit.Chem import rdDepictor

    mol = Chem.MolFromSmiles(smiles)
    rdDepictor.Compute2DCoords(mol)
    Chem.WedgeMolBonds(mol, mol.GetConformer())
    return Chem.MolToMolBlock(mol)


# -- reading ---------------------------------------------------------------------


def test_a_failing_backend_is_an_error_with_its_reason(log):
    with pytest.raises(BackendError, match="model exploded"):
        registry.recognize(None, [FakeBackend("broken", "", fails=True)], log)


def test_no_backend_at_all_explains_how_to_get_one(log):
    with pytest.raises(BackendError) as excinfo:
        registry.recognize(None, [], log)
    message = str(excinfo.value)
    assert "m2i setup decimer" in message
    assert "ChemDraw" in message  # the file needs no model


def test_low_confidence_is_surfaced(log):
    registry.recognize(None, [FakeBackend("a", "CCO", confidence=0.4)], log)
    assert "recognition.low_confidence" in codes(log)


def test_high_confidence_is_not_a_warning(log):
    registry.recognize(None, [FakeBackend("a", "CCO", confidence=0.95)], log)
    assert "recognition.low_confidence" not in codes(log)


def test_sequence_model_gets_a_stereochemistry_caveat(log):
    registry.recognize(None, [FakeBackend("a", S_ENANTIOMER, confidence=0.9)], log)
    message = next(i for i in log if i.code == "recognition.sequence_model").message
    assert "measured" in message


def test_a_molblock_gets_no_caveat(log):
    backend = FakeBackend("a", S_ENANTIOMER, molblock=molblock_for(S_ENANTIOMER), confidence=0.9)
    registry.recognize(None, [backend], log)
    assert "recognition.sequence_model" not in codes(log)


# -- selection -----------------------------------------------------------------------


def test_the_installed_model_is_used(monkeypatch, log):
    fake = FakeBackend("decimer", "CCO")
    monkeypatch.setattr(registry, "installed_vision_backends", lambda: [fake])
    assert registry.resolve_backends(None, log) == [fake]
    assert "backend.auto" in codes(log)


def test_nothing_installed_selects_nothing(monkeypatch, log):
    monkeypatch.setattr(registry, "installed_vision_backends", list)
    assert registry.resolve_backends(None, log) == []


def test_an_unknown_backend_name_is_reported(log):
    assert registry.resolve_backends(["nonesuch"], log) == []
    assert "backend.unknown" in codes(log)


def test_screenshots_of_chemdraw_have_no_model_of_their_own(log):
    """MolScribe was dropped: the ChemDraw file is read exactly instead."""
    assert registry.resolve_backends(["molscribe"], log) == []
    assert "backend.unknown" in codes(log)


def test_an_uninstalled_backend_name_says_how_to_install_it(log):
    assert registry.resolve_backends(["decimer"], log) == []
    message = next(i for i in log if i.code == "backend.unavailable").message
    assert "m2i setup decimer" in message


def test_manual_is_always_available():
    rows = {row["name"]: row for row in registry.describe_backends()}
    assert rows["manual"]["available"]
    assert set(rows) == {"manual", "decimer"}
