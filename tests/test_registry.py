"""Backend selection and multi-backend consensus, with stand-in backends."""

from __future__ import annotations

import pytest
from conftest import codes

from m2i.recognition import registry
from m2i.recognition.base import BackendError
from m2i.types import IssueLog, RecognitionResult

# Two readings of the same skeleton that differ only in the stereocentre.
S_ENANTIOMER = "F[C@@H](Cl)Br"
R_ENANTIOMER = "F[C@H](Cl)Br"


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


# -- consensus -----------------------------------------------------------


def test_full_agreement_is_the_strongest_signal(log):
    backends = [
        FakeBackend("a", S_ENANTIOMER, confidence=0.8),
        FakeBackend("b", S_ENANTIOMER, confidence=0.7),
    ]
    result = registry.recognize(None, backends, log)
    assert result.smiles == S_ENANTIOMER
    assert "consensus.agree" in codes(log)
    assert not log.warnings


def test_same_skeleton_different_stereochemistry_warns(log):
    backends = [
        FakeBackend("a", S_ENANTIOMER, confidence=0.9),
        FakeBackend("b", R_ENANTIOMER, confidence=0.5),
    ]
    registry.recognize(None, backends, log)
    assert "consensus.stereo_disagreement" in codes(log)
    message = next(i for i in log if i.code == "consensus.stereo_disagreement").message
    assert "wedge" in message


def test_different_molecules_warn_loudly(log):
    backends = [
        FakeBackend("a", "CCO", confidence=0.9),
        FakeBackend("b", "c1ccccc1", confidence=0.4),
    ]
    registry.recognize(None, backends, log)
    assert "consensus.disagreement" in codes(log)


def test_a_molblock_beats_a_more_confident_smiles(log):
    """Measured stereochemistry beats generated stereochemistry."""
    graph = FakeBackend("graph", S_ENANTIOMER, molblock=molblock_for(S_ENANTIOMER), confidence=0.5)
    sequence = FakeBackend("sequence", R_ENANTIOMER, confidence=0.99)
    result = registry.recognize(None, [graph, sequence], log)
    assert result.backend == "graph"


def test_highest_confidence_wins_among_equals(log):
    weak = FakeBackend("weak", "CCO", confidence=0.3)
    strong = FakeBackend("strong", "c1ccccc1", confidence=0.95)
    assert registry.recognize(None, [weak, strong], log).backend == "strong"


def test_a_failing_backend_does_not_sink_the_run(log):
    backends = [FakeBackend("broken", "", fails=True), FakeBackend("ok", "CCO")]
    assert registry.recognize(None, backends, log).smiles == "CCO"
    assert "backend.failed" in codes(log)


def test_all_backends_failing_is_an_error(log):
    with pytest.raises(BackendError):
        registry.recognize(None, [FakeBackend("broken", "", fails=True)], log)


def test_no_backend_at_all_explains_how_to_get_one(log):
    with pytest.raises(BackendError) as excinfo:
        registry.recognize(None, [], log)
    assert "m2i setup" in str(excinfo.value)


# -- single-backend reporting -------------------------------------------


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


def test_graph_model_gets_no_caveat(log):
    backend = FakeBackend("a", S_ENANTIOMER, molblock=molblock_for(S_ENANTIOMER), confidence=0.9)
    registry.recognize(None, [backend], log)
    assert "recognition.sequence_model" not in codes(log)


# -- selection -----------------------------------------------------------


def install_fakes(monkeypatch, *names):
    fakes = {n: FakeBackend(n, "CCO") for n in names}
    monkeypatch.setattr(registry, "installed_vision_backends", lambda: list(fakes.values()))
    return fakes


def test_hand_drawn_images_prefer_decimer(monkeypatch, log):
    install_fakes(monkeypatch, "molscribe", "decimer")
    chosen = registry.resolve_backends(None, log, style="hand_drawn")
    assert [b.name for b in chosen] == ["decimer"]
    assert "backend.auto" in codes(log)


def test_clean_images_prefer_molscribe(monkeypatch, log):
    install_fakes(monkeypatch, "molscribe", "decimer")
    assert [b.name for b in registry.resolve_backends(None, log, style="clean")] == [
        "molscribe"
    ]


def test_the_preferred_backend_falls_back_to_what_is_installed(monkeypatch, log):
    install_fakes(monkeypatch, "molscribe")
    assert [b.name for b in registry.resolve_backends(None, log, style="hand_drawn")] == [
        "molscribe"
    ]


def test_all_runs_every_installed_backend(monkeypatch, log):
    install_fakes(monkeypatch, "molscribe", "decimer")
    chosen = registry.resolve_backends(["all"], log)
    assert len(chosen) == 2
    assert "backend.consensus" in codes(log)


def test_nothing_installed_selects_nothing(monkeypatch, log):
    monkeypatch.setattr(registry, "installed_vision_backends", list)
    assert registry.resolve_backends(None, log, style="clean") == []


def test_an_unknown_backend_name_is_reported(log):
    assert registry.resolve_backends(["nonesuch"], log) == []
    assert "backend.unknown" in codes(log)


def test_an_uninstalled_backend_name_says_how_to_install_it(log):
    assert registry.resolve_backends(["decimer"], log) == []
    message = next(i for i in log if i.code == "backend.unavailable").message
    assert "m2i setup decimer" in message


def test_manual_is_always_available():
    rows = {row["name"]: row for row in registry.describe_backends()}
    assert rows["manual"]["available"]
    assert set(rows) == {"manual", "molscribe", "decimer"}
