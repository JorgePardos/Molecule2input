"""MolScribe worker: image -> molecular graph -> molfile + SMILES.

This is the backend that makes stereochemistry *measured* rather than
generated: it predicts atoms with 2D coordinates and bond types including
solid and hashed wedges, so m2i can read R/S and E/Z off the layout exactly as
a chemist would.

Runs inside its own venv (numpy<2, opencv 4.5.5, torch). Never import m2i here.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from _protocol import package_versions, serve  # noqa: E402

CHECKPOINT_REPO = "yujieq/MolScribe"
CHECKPOINT_FILE = "swin_base_char_aux_1m.pth"

_MODEL = None


def checkpoint_path(request: dict) -> str:
    explicit = request.get("checkpoint")
    if explicit and Path(explicit).is_file():
        return explicit
    from huggingface_hub import hf_hub_download

    return hf_hub_download(CHECKPOINT_REPO, CHECKPOINT_FILE)


def load_model(request: dict):
    global _MODEL
    if _MODEL is None:
        import torch
        from molscribe import MolScribe

        device = request.get("device") or (
            "cuda" if os.environ.get("M2I_USE_CUDA") else "cpu"
        )
        _MODEL = MolScribe(checkpoint_path(request), device=torch.device(device))
    return _MODEL


def handle(request: dict) -> dict:
    mode = request.get("mode", "recognize")
    versions = package_versions("MolScribe", "torch", "numpy", "opencv-python")

    if mode == "warmup":
        # Downloading the checkpoint is the slow part of a first run; do it at
        # install time so the first real recognition is not a long wait. Then
        # actually load the model and read something: a broken dependency
        # resolution leaves the checkpoint on disk just the same.
        path = checkpoint_path(request)
        probe = None
        if request.get("image_path"):
            probe = load_model(request).predict_image_file(request["image_path"])[
                "smiles"
            ]
        return {"checkpoint": path, "versions": versions, "probe": probe}

    model = load_model(request)
    output = model.predict_image_file(
        request["image_path"], return_atoms_bonds=True, return_confidence=True
    )

    atoms = output.get("atoms") or []
    bonds = output.get("bonds") or []
    return {
        "smiles": output.get("smiles", ""),
        # MolScribe calls it 'molfile'; m2i calls the same thing a molblock.
        "molblock": output.get("molfile"),
        "confidence": output.get("confidence"),
        "raw": {
            "n_atoms": len(atoms),
            "n_bonds": len(bonds),
            "min_atom_confidence": _minimum(atoms, "atom_confidence"),
            "min_bond_confidence": _minimum(bonds, "bond_confidence"),
            "versions": versions,
        },
    }


def _minimum(items, key):
    values = [i.get(key) for i in items if isinstance(i, dict) and i.get(key) is not None]
    return min(values) if values else None


if __name__ == "__main__":
    raise SystemExit(serve(handle))
