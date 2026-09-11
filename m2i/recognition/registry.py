"""Backend discovery, automatic selection and multi-backend consensus.

No single OCSR model wins on both inputs: on hand-drawn structures the DECIMER
family is far ahead, while on clean ChemDraw-style depictions the image-to-graph
models win because they measure wedges instead of generating them. So m2i picks
by what the image looks like, and when both are installed it can run them
together and compare InChIKeys -- turning two unreliable readings into a
trustworthy confidence signal.
"""

from __future__ import annotations

from pathlib import Path

from rdkit import Chem

from ..types import IssueLog, RecognitionResult
from . import _subprocess
from .base import BackendError, OCSRBackend
from .manual import ManualBackend

#: Which backend to prefer for each drawing style, best first.
STYLE_PREFERENCE = {
    "hand_drawn": ("decimer", "molscribe"),
    "clean": ("molscribe", "decimer"),
}


class SubprocessBackend:
    """A vision model living in its own virtual environment."""

    def __init__(self, spec) -> None:
        self.spec = spec
        self.name = spec.name
        self.description = spec.description
        self.returns_molblock = spec.returns_molblock
        self.strength = spec.strength

    def available(self) -> tuple[bool, str]:
        if _subprocess.is_installed(self.name):
            marker = _subprocess.read_marker(self.name) or {}
            when = marker.get("installed_utc", "unknown date")
            return True, f"installed {when}"
        if _subprocess.venv_python(self.name).is_file():
            return False, f"environment exists but is incomplete; run: m2i setup {self.name} --force"
        return False, f"not installed; run: m2i setup {self.name}  ({self.spec.download_size})"

    def recognize(self, image_path: Path, **options) -> RecognitionResult:
        ok, reason = self.available()
        if not ok:
            raise BackendError(f"{self.name} is not usable: {reason}")
        timeout = options.pop("timeout", _subprocess.DEFAULT_TIMEOUT)
        return _subprocess.recognize_with_worker(
            self.name,
            self.spec.worker,
            image_path,
            {"mode": "recognize", **options},
            timeout=timeout,
        )


def available_backends() -> dict[str, OCSRBackend]:
    """Every known backend, installed or not."""
    from ..backends import SPECS

    backends: dict[str, OCSRBackend] = {"manual": ManualBackend()}
    for name, spec in SPECS.items():
        backends[name] = SubprocessBackend(spec)
    return backends


def describe_backends() -> list[dict]:
    rows = []
    for name, backend in available_backends().items():
        ok, reason = backend.available()
        rows.append(
            {
                "name": name,
                "available": ok,
                "reason": reason,
                "description": backend.description,
                "returns_molblock": backend.returns_molblock,
                "strength": getattr(backend, "strength", "any input you type yourself"),
            }
        )
    return rows


def installed_vision_backends() -> list[OCSRBackend]:
    return [
        backend
        for name, backend in available_backends().items()
        if name != "manual" and backend.available()[0]
    ]


def resolve_backends(
    names: list[str] | None, log: IssueLog, *, style: str | None = None
) -> list[OCSRBackend]:
    """Turn user-requested names into usable backend objects.

    ``auto`` (the default) picks the single backend best suited to the drawing
    style; ``all`` runs every installed one so their answers can be compared.
    """
    known = available_backends()

    if not names or names == ["auto"]:
        return _auto_select(style, log)
    if names == ["all"]:
        chosen = installed_vision_backends()
        if len(chosen) > 1:
            log.info(
                "backend.consensus",
                f"Running {len(chosen)} backends and comparing their answers.",
            )
        return chosen

    chosen: list[OCSRBackend] = []
    for name in names:
        backend = known.get(name)
        if backend is None:
            log.warn(
                "backend.unknown",
                f"Unknown backend {name!r}; known: {', '.join(known)}.",
            )
            continue
        ok, reason = backend.available()
        if not ok:
            log.warn("backend.unavailable", f"Backend {name!r} unavailable: {reason}")
            continue
        chosen.append(backend)
    return chosen


def _auto_select(style: str | None, log: IssueLog) -> list[OCSRBackend]:
    installed = {b.name: b for b in installed_vision_backends()}
    if not installed:
        return []

    order = STYLE_PREFERENCE.get(style or "", ())
    for name in order:
        if name in installed:
            if style:
                log.info(
                    "backend.auto",
                    f"The drawing looks {style.replace('_', '-')}, so {name} was "
                    f"chosen ({installed[name].strength}). Use --backend to override, "
                    "or --backend all to run every installed model and compare.",
                )
            return [installed[name]]

    # No preference matched (unknown style, or only the other model installed).
    backend = next(iter(installed.values()))
    log.info("backend.auto", f"Using {backend.name} ({backend.strength}).")
    return [backend]


def recognize(
    image_path: Path | None,
    backends: list[OCSRBackend],
    log: IssueLog,
    **options,
) -> RecognitionResult:
    """Run the backends and reconcile their answers."""
    if not backends:
        raise BackendError(
            "no recognition backend available. Supply the structure yourself with "
            "--smiles/--molfile, or install one with `m2i setup decimer` "
            "(hand-drawn) or `m2i setup molscribe` (ChemDraw-style)."
        )

    results: list[RecognitionResult] = []
    for backend in backends:
        try:
            results.append(backend.recognize(image_path, **options))
        except Exception as exc:
            log.warn("backend.failed", f"Backend {backend.name} failed: {exc}")

    if not results:
        raise BackendError("every recognition backend failed; see the warnings above")
    if len(results) == 1:
        _report_single(results[0], log)
        return results[0]
    return _reconcile(results, log)


def _report_single(result: RecognitionResult, log: IssueLog) -> None:
    if result.confidence is not None:
        log.info(
            "recognition.confidence",
            f"{result.backend} reported a confidence of {result.confidence:.2f}.",
        )
    if result.confidence is not None and result.confidence < 0.7:
        log.warn(
            "recognition.low_confidence",
            f"{result.backend} is not confident about this reading "
            f"({result.confidence:.2f}). Check every atom and wedge.",
        )
    if not result.molblock:
        log.info(
            "recognition.sequence_model",
            f"{result.backend} returned only a SMILES, so its stereochemistry was "
            "generated rather than measured off the drawing. The stereocentres "
            "deserve a closer look than the connectivity.",
        )


def _reconcile(results: list[RecognitionResult], log: IssueLog) -> RecognitionResult:
    keys = {r.backend: _inchikey(r) for r in results}
    full = {k for k in keys.values() if k}
    skeletons = {k.split("-")[0] for k in full}

    if len(full) == 1:
        log.info(
            "consensus.agree",
            f"{len(results)} backends agree on the full structure, stereochemistry "
            "included. That is the strongest signal m2i can give you.",
        )
    elif len(skeletons) == 1:
        detail = "; ".join(f"{b}: {k}" for b, k in keys.items())
        log.warn(
            "consensus.stereo_disagreement",
            "The backends agree on the connectivity but disagree on the "
            f"stereochemistry ({detail}). Check every wedge before trusting the "
            "generated input.",
        )
    else:
        detail = "; ".join(f"{r.backend}: {r.smiles}" for r in results)
        log.warn(
            "consensus.disagreement",
            f"The backends read different molecules ({detail}). The reading with "
            "the highest confidence was kept, but this drawing needs manual review.",
        )

    chosen = max(results, key=_ranking)
    log.info("consensus.chosen", f"Kept the reading from {chosen.backend}.")
    return chosen


def _ranking(result: RecognitionResult) -> tuple[int, float]:
    # A molblock beats a bare SMILES: measured stereochemistry beats generated.
    return (1 if result.molblock else 0, result.confidence or 0.0)


def _inchikey(result: RecognitionResult) -> str:
    try:
        mol = (
            Chem.MolFromMolBlock(result.molblock)
            if result.molblock
            else Chem.MolFromSmiles(result.smiles)
        )
        return Chem.MolToInchiKey(mol) if mol else ""
    except Exception:
        return ""
