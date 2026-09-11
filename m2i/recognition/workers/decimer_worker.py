"""DECIMER worker: image -> SMILES.

The strongest backend for hand-drawn structures by a wide margin, at the cost
of returning a token sequence rather than a graph: its stereochemistry is
generated with the SMILES, not measured off the drawing. m2i records that
distinction so the confidence signal stays honest.

Runs inside its own venv (TensorFlow). Never import m2i here.
"""

from __future__ import annotations

import inspect
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from _protocol import package_versions, serve  # noqa: E402

_PREDICT = None


def load_predict():
    """Import DECIMER lazily; the import itself fetches the model weights."""
    global _PREDICT
    if _PREDICT is None:
        from DECIMER import predict_SMILES

        _PREDICT = predict_SMILES
    return _PREDICT


def _supported_arguments(function) -> set:
    try:
        return set(inspect.signature(function).parameters)
    except (TypeError, ValueError):
        return set()


def handle(request: dict) -> dict:
    mode = request.get("mode", "recognize")
    versions = package_versions("decimer", "tensorflow")

    predict = load_predict()
    supported = _supported_arguments(predict)

    if mode == "warmup":
        # Importing DECIMER downloads the weights. If a probe image was given,
        # run one real prediction: that is the only proof the install works.
        probe = None
        if request.get("image_path"):
            probe, _ = _split(_predict(predict, supported, request))
        return {"versions": versions, "signature": sorted(supported), "probe": probe}

    prediction = _predict(predict, supported, request)
    smiles, confidence = _split(prediction)
    return {
        "smiles": smiles,
        # DECIMER is image-to-sequence: there is no 2D layout to read wedges from.
        "molblock": None,
        "confidence": confidence,
        "raw": {"versions": versions, "hand_drawn": _wants_hand_drawn(supported, request)},
    }


def _predict(predict, supported: set, request: dict):
    kwargs = {}
    if "hand_drawn" in supported:
        kwargs["hand_drawn"] = _wants_hand_drawn(supported, request)
    if "confidence" in supported and request.get("want_confidence", True):
        kwargs["confidence"] = True
    try:
        return predict(request["image_path"], **kwargs)
    except TypeError:
        # Signature introspection can fail on decorated functions; fall back to
        # the one argument every version of DECIMER has accepted.
        return predict(request["image_path"])


def _wants_hand_drawn(supported: set, request: dict) -> bool:
    return bool(request.get("hand_drawn", False)) and "hand_drawn" in supported


def _split(prediction):
    """DECIMER returns either a SMILES or a (SMILES, confidence) pair."""
    if isinstance(prediction, str):
        return prediction, None
    if isinstance(prediction, (tuple, list)) and prediction:
        smiles = prediction[0]
        rest = prediction[1] if len(prediction) > 1 else None
        return str(smiles), _mean_confidence(rest)
    return str(prediction), None


def _mean_confidence(value):
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    try:
        # Some versions return a per-token list of (token, probability) pairs.
        numbers = [
            float(item[1]) if isinstance(item, (tuple, list)) else float(item)
            for item in value
        ]
        return sum(numbers) / len(numbers) if numbers else None
    except (TypeError, ValueError, IndexError):
        return None


if __name__ == "__main__":
    raise SystemExit(serve(handle))
