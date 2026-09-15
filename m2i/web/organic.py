"""Organic molecules: give a structure, check it when in doubt, write the input.

The routes behind screens 2a-2d and 2g. Ports of ``picture_input``,
``file_input``, ``smiles_input``, ``check_step`` and ``output_step``.
"""

from __future__ import annotations

import hashlib
import io
from pathlib import Path

from fastapi import APIRouter, Depends, HTTPException, Query, UploadFile
from fastapi.responses import FileResponse, Response
from pydantic import BaseModel
from rdkit import Chem

from .. import config, pipeline
from ..chem.conformers import ConformerOptions
from ..organometallic import has_metal, read_raw
from ..preprocess import ImageError, estimate_drawing_style, prepare_image
from ..recognition import BackendError, recognize
from ..recognition.manual import STRUCTURE_FILE_SUFFIXES, from_structure_file
from ..recognition.registry import resolve_backends
from ..report import depict, review
from ..types import IssueLog, RecognitionResult
from . import recipes, shapes
from .session import Session, Source, current

router = APIRouter(prefix="/api")

IMAGE_SUFFIXES = (".png", ".jpg", ".jpeg", ".bmp", ".tif", ".tiff", ".webp")
CIF_SUFFIXES = (".cif",)


def error(message: str, status: int = 400):
    raise HTTPException(status_code=status, detail=message)


def save_upload(session: Session, upload: UploadFile, allowed) -> tuple[Path, str]:
    """Store an upload under a name m2i chooses.

    The browser-supplied file name never becomes part of a path: only its
    extension is kept, checked against the types this route accepts, and the
    file is named by its own hash. The size is capped by the server's upload
    limit, checked again here.
    """
    suffix = Path(upload.filename or "").suffix.lower()
    if suffix not in allowed:
        error(f"{suffix or 'That file'} is not accepted here: {', '.join(allowed)}")
    data = upload.file.read(max_upload_bytes() + 1)
    if len(data) > max_upload_bytes():
        error(f"That file is larger than {max_upload_bytes() // 2**20} MB.")
    digest = hashlib.sha256(data).hexdigest()[:16]
    path = session.folder / f"{digest}{suffix}"
    if not path.exists():
        path.write_bytes(data)
    return path, digest


def max_upload_bytes() -> int:
    import os

    return int(os.environ.get("M2I_MAX_UPLOAD_MB") or 20) * 2**20


def manual(smiles: str) -> RecognitionResult:
    return RecognitionResult(smiles=smiles, backend="manual", confidence=1.0)


def display_smiles(recognition: RecognitionResult) -> str:
    """What to put in the correction box: the reading itself, as a SMILES."""
    if recognition.molblock:
        mol = Chem.MolFromMolBlock(recognition.molblock)
        if mol is not None:
            return Chem.MolToSmiles(mol)
    return recognition.smiles or ""


def remember(session: Session, source: Source) -> dict:
    session.sources[source.key] = source
    return {
        "key": source.key,
        "kind": source.kind,
        "shown_name": source.shown_name,
        "style": source.style,
        "read": source.recognition is not None,
        "preview_url": f"/api/source/{source.key}/preview" if source.path and source.kind == "picture" else None,
        "issues": shapes.issues(source.issues),
    }


# -- step 1: structure --------------------------------------------------------------


@router.post("/source/image")
def source_image(file: UploadFile, session: Session = Depends(current)):
    with session.lock:
        path, digest = save_upload(session, file, IMAGE_SUFFIXES)
        log = IssueLog()
        try:
            image = prepare_image(path, log)
        except ImageError as exc:
            error(f"That file could not be opened as a picture: {exc}")
        style, _ = estimate_drawing_style(image)
        source = Source(kind="picture", key=f"picture:{digest}", path=path, image=image, style=style,
                        shown_name=Path(file.filename or "").name, issues=list(log))
        existing = session.sources.get(source.key)
        if existing is not None and existing.recognition is not None:
            source = existing  # the same picture again: keep its reading
        return remember(session, source)


@router.post("/source/{key}/read")
def source_read(key: str, session: Session = Depends(current)):
    with session.lock:
        source = _source(session, key)
        if source.kind != "picture":
            error("Only a picture needs reading.")
        log = IssueLog(list(source.issues))
        backends = resolve_backends(None, log)
        try:
            recognition = recognize(source.path, backends, log, hand_drawn=(source.style == "hand_drawn"))
        except BackendError as exc:
            error(str(exc))
        source.recognition = recognition
        source.issues = list(log)
        source.display_smiles = display_smiles(recognition)
        return remember(session, source)


@router.get("/source/{key}/preview")
def source_preview(key: str, session: Session = Depends(current)):
    source = _source(session, key)
    if source.image is None:
        error("No picture for this structure.", 404)
    buffer = io.BytesIO()
    source.image.convert("RGB").save(buffer, "PNG")
    return Response(buffer.getvalue(), media_type="image/png", headers={"Cache-Control": "private, max-age=3600"})


@router.post("/source/file")
def source_file(file: UploadFile, session: Session = Depends(current)):
    with session.lock:
        path, digest = save_upload(session, file, tuple(STRUCTURE_FILE_SUFFIXES))
        shown = Path(file.filename or "").name
        try:
            metal = any(has_metal(m) for m in read_raw(path))
        except Exception:  # noqa: BLE001 - the organic reader reports a bad file properly
            metal = False
        if metal:
            key = f"om:{digest}"
            session.sources[key] = Source(kind="om", key=key, path=path, shown_name=shown)
            return {"route": "organometallic", "key": key, "shown_name": shown}
        try:
            recognition = from_structure_file(path).recognize(None)
        except Exception as exc:  # noqa: BLE001 - RDKit raises assorted types on bad files
            error(f"That file could not be read: {exc}")
        source = Source(kind="file", key=f"file:{digest}", path=path, shown_name=shown,
                        recognition=recognition, display_smiles=display_smiles(recognition))
        return {"route": "organic", **remember(session, source), **_file_counts(recognition)}


class SmilesIn(BaseModel):
    smiles: str


@router.post("/source/smiles")
def source_smiles(body: SmilesIn, session: Session = Depends(current)):
    text = body.smiles.strip()
    if not text:
        error("Type a SMILES.")
    with session.lock:
        key = "smiles:" + hashlib.sha256(text.encode()).hexdigest()[:16]
        recognition = manual(text)
        source = Source(kind="smiles", key=key, recognition=recognition, display_smiles=text)
        return remember(session, source)


def _file_counts(recognition: RecognitionResult) -> dict:
    mol = Chem.MolFromMolBlock(recognition.molblock, sanitize=False) if recognition.molblock else None
    if mol is None:
        return {}
    wedge_dirs = (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH)
    wedges = sum(1 for b in mol.GetBonds() if b.GetBondDir() in wedge_dirs
                 or (b.HasProp("_MolFileBondStereo") and b.GetIntProp("_MolFileBondStereo") in (1, 6)))
    return {"atoms": mol.GetNumAtoms(), "bonds": mol.GetNumBonds(), "wedges": wedges}


# -- step 2: check ------------------------------------------------------------------


class Overrides(BaseModel):
    charge: int | None = None
    multiplicity: int | None = None
    keep_all_fragments: bool = False


def overrides(
    charge: int | None = Query(None), multiplicity: int | None = Query(None),
    keep_all_fragments: bool = Query(False),
) -> Overrides:
    return Overrides(charge=charge, multiplicity=multiplicity, keep_all_fragments=keep_all_fragments)


def _source(session: Session, key: str) -> Source:
    source = session.sources.get(key)
    if source is None:
        error("That structure is no longer in this session; give it again.", 404)
    return source


def evaluate(session: Session, key: str, extra: Overrides):
    """The structure as it stands: reading or correction, molecule, decision."""
    source = _source(session, key)
    if source.recognition is None:
        error("This picture has not been read yet.", 409)
    corrected = session.corrections.get(key)
    recognition = corrected or source.recognition
    log = IssueLog() if corrected else IssueLog(source.issues)
    if corrected:
        log.info("recognition.corrected",
                 f"Structure corrected by hand. The original reading "
                 f"({source.recognition.backend}) was {source.display_smiles}.")
    molecule, problem = None, None
    try:
        molecule = pipeline.prepare_molecule(
            recognition, log, charge=extra.charge, multiplicity=extra.multiplicity,
            keep_all_fragments=extra.keep_all_fragments,
        )
    except Exception as exc:  # noqa: BLE001 - shown to the user, who can fix it
        problem = str(exc)
    decision = review.assess(molecule, log) if molecule is not None else None
    return source, recognition, corrected, molecule, problem, log, decision


@router.get("/check/{key}")
def check(key: str, extra: Overrides = Depends(overrides), session: Session = Depends(current)):
    with session.lock:
        source, recognition, corrected, molecule, problem, log, decision = evaluate(session, key, extra)
        confirmed = molecule is not None and (key, molecule.smiles) in session.confirmed
        return {
            "key": key,
            "kind": source.kind,
            "shown_name": source.shown_name,
            "molecule": shapes.molecule(molecule) if molecule else None,
            "error": problem,
            "depiction_url": f"/api/depiction/{key}.png?smiles={_quote(molecule.smiles)}"
            f"&keep_all_fragments={str(extra.keep_all_fragments).lower()}" if molecule else None,
            "preview_url": f"/api/source/{key}/preview" if source.image is not None else None,
            "source": {
                "backend": recognition.backend,
                "confidence": recognition.confidence,
                "display_smiles": source.display_smiles,
                "corrected": corrected is not None,
                "correction": corrected.smiles if corrected else None,
            },
            "gate": {
                "needed": bool(decision and decision.needed) or molecule is None,
                "reasons": decision.reasons if decision else [],
                "summary": decision.summary() if decision else "",
                "confirmed": confirmed,
            },
            "errors": log.has_errors(),
            "issues": shapes.issues(log),
        }


def _quote(text: str) -> str:
    from urllib.parse import quote

    return quote(text, safe="")


class CorrectionIn(BaseModel):
    smiles: str


@router.post("/check/{key}/correct")
def correct(key: str, body: CorrectionIn, session: Session = Depends(current)):
    text = body.smiles.strip()
    with session.lock:
        source = _source(session, key)
        applied = session.corrections.get(key)
        applied_text = applied.smiles if applied else source.display_smiles
        # Compared with what was last applied, not with the canonical SMILES:
        # a correction typed in non-canonical form would otherwise differ from
        # its own canonicalisation every time.
        if text and text != applied_text:
            session.corrections[key] = manual(text)
        return {"ok": True}


@router.delete("/check/{key}/correct")
def uncorrect(key: str, session: Session = Depends(current)):
    with session.lock:
        session.corrections.pop(key, None)
        return {"ok": True}


class ConfirmIn(BaseModel):
    smiles: str
    confirmed: bool = True


@router.post("/check/{key}/confirm")
def confirm(key: str, body: ConfirmIn, session: Session = Depends(current)):
    """Confirmation belongs to this exact structure: a correction resets it."""
    with session.lock:
        _source(session, key)
        entry = (key, body.smiles)
        if body.confirmed:
            session.confirmed.add(entry)
        else:
            session.confirmed.discard(entry)
        return {"ok": True}


@router.get("/depiction/{key}.png")
def depiction(key: str, keep_all_fragments: bool = False, session: Session = Depends(current)):
    with session.lock:
        _, _, _, molecule, problem, _, _ = evaluate(session, key, Overrides(keep_all_fragments=keep_all_fragments))
        if molecule is None:
            error(problem or "This structure cannot be drawn.")
        png = depict.draw_molecule(molecule.mol, stereo=molecule.stereo)
    return Response(png, media_type="image/png", headers={"Cache-Control": "private, max-age=86400"})


# -- step 3: output -----------------------------------------------------------------


class ConformersIn(BaseModel):
    n_confs: int | None = None
    keep: int = 1
    seed: int = 0xF00D
    force_field: str = "mmff94s"


class GenerateIn(BaseModel):
    key: str
    settings: recipes.Settings
    conformers: ConformersIn = ConformersIn()
    overrides: Overrides = Overrides()


@router.post("/generate")
def generate(body: GenerateIn, session: Session = Depends(current)):
    with session.lock:
        source, _, _, molecule, problem, log, decision = evaluate(session, body.key, body.overrides)
        if molecule is None:
            error(problem or "This structure cannot be used.")
        if log.has_errors():
            error("Resolve the errors above before generating an input.")
        if decision.needed and (body.key, molecule.smiles) not in session.confirmed:
            error("Confirm the structure, or correct its SMILES, to generate an input.", 409)
        try:
            profile = recipes.build(body.settings)
        except config.ProfileError as exc:
            error(str(exc))
        review.record(log, decision, confirmed=True if decision.needed else None)

        options = ConformerOptions(
            n_confs=body.conformers.n_confs or None, keep=max(1, body.conformers.keep),
            seed=body.conformers.seed, force_field=body.conformers.force_field,
        )
        cache_key = (molecule.smiles, repr(options))
        if cache_key not in session.embeddings:
            embed_log = IssueLog()
            try:
                embedding = pipeline.embed(molecule, options, embed_log)
            except Exception as exc:  # noqa: BLE001
                error(f"Could not build a 3D structure: {exc}")
            session.embeddings[cache_key] = (embedding, list(embed_log))
        embedding, embed_issues = session.embeddings[cache_key]

        run_log = IssueLog(list(log) + embed_issues)
        try:
            result = pipeline.write_inputs(
                molecule, embedding,
                pipeline.PipelineOptions(profile=profile, output_dir=session.output,
                                         conformers=options, source_image=source.image),
                run_log,
            )
        except Exception as exc:  # noqa: BLE001
            error(f"Could not write the input: {exc}")

        written = [Path(p) for p in result.written_files]
        inputs = [p for p in written if not p.name.endswith(("_check.png", ".m2i.json"))]
        files = []
        for conformer, path in zip(result.conformers, inputs):
            entry = session.offer([path])[0]
            entry.update({
                "kind": "input",
                "preview": path.read_text(encoding="utf-8", errors="replace"),
                "conformer": {
                    "index": conformer.index, "relative_energy": conformer.relative_energy,
                    "force_field": conformer.force_field,
                },
            })
            files.append(entry)
        extras = []
        for path in written:
            if path in inputs:
                continue
            entry = session.offer([path])[0]
            entry["kind"] = "provenance" if path.name.endswith(".m2i.json") else "check_image"
            extras.append(entry)
        return {
            "format": body.settings.format,
            "format_label": recipes.FORMATS[body.settings.format],
            "jobs": " ".join(profile.jobs) or "single point",
            "files": files,
            "extras": extras,
            "conformers_kept": len(result.conformers),
            "issues": shapes.new_issues(run_log, log),
        }


@router.get("/download/{name}")
def download(name: str, session: Session = Depends(current)):
    path = session.files.get(name)
    if path is None or not path.is_file():
        error("That file is not available in this session.", 404)
    return FileResponse(path, filename=path.name)
