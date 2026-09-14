"""Backend discovery and selection for reading pictures.

One vision model is supported: DECIMER, the one that reads hand-drawn
structures. Drawings made in ChemDraw do not need a model at all -- the file
itself holds the bonds and the wedges, and is read exactly -- so a second model
for screenshots of them was dropped: uploading the .cdx is always better.
"""

from __future__ import annotations

from pathlib import Path

from ..types import IssueLog, RecognitionResult
from . import _subprocess
from .base import BackendError, OCSRBackend
from .manual import ManualBackend

#: Below this, a reading is reported as uncertain in so many words.
LOW_CONFIDENCE = 0.7


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


def resolve_backends(names: list[str] | None, log: IssueLog) -> list[OCSRBackend]:
    """The model to read a picture with: the installed one, or the one named."""
    if not names or names == ["auto"]:
        installed = installed_vision_backends()
        if installed:
            log.info("backend.auto", f"Reading the picture with {installed[0].name}.")
        return installed[:1]

    known = available_backends()
    chosen: list[OCSRBackend] = []
    for name in names:
        backend = known.get(name)
        if backend is None or name == "manual":
            log.warn(
                "backend.unknown",
                f"Unknown backend {name!r}; known: {', '.join(n for n in known if n != 'manual')}.",
            )
            continue
        ok, reason = backend.available()
        if not ok:
            log.warn("backend.unavailable", f"Backend {name!r} unavailable: {reason}")
            continue
        chosen.append(backend)
    return chosen[:1]


def recognize(
    image_path: Path | None,
    backends: list[OCSRBackend],
    log: IssueLog,
    **options,
) -> RecognitionResult:
    """Read the picture, and say how much the reading can be trusted."""
    if not backends:
        raise BackendError(
            "no recognition backend available. Supply the structure yourself with "
            "--smiles or --molfile (a ChemDraw file is read exactly), or install the "
            "model for hand-drawn structures with `m2i setup decimer`."
        )
    failures = []
    for backend in backends:
        try:
            result = backend.recognize(image_path, **options)
        except Exception as exc:  # noqa: BLE001 - reported below, with its reason
            failures.append(f"{backend.name}: {exc}")
            continue
        _report(result, log)
        return result
    raise BackendError("the picture could not be read (" + "; ".join(failures) + ")")


def _report(result: RecognitionResult, log: IssueLog) -> None:
    if result.confidence is not None:
        log.info(
            "recognition.confidence",
            f"{result.backend} reported a confidence of {result.confidence:.2f}.",
        )
    if result.confidence is not None and result.confidence < LOW_CONFIDENCE:
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
