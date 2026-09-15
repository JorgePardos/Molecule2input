"""Loading the photo model once per server, and saying how far it has got.

A model takes most of a minute to load and a couple of seconds to read a
picture once loaded, so it is loaded when the server starts, in the
background. The front end polls ``/api/health`` and shows a splash with this
state while the routes that need no model are already usable.

TensorFlow gives no progress while it loads, so none is invented: the state is
waiting, loading (with how long it has taken and how long it usually takes),
ready, failed, or absent when no model is installed.
"""

from __future__ import annotations

import threading
import time

from ..recognition import _subprocess
from ..recognition.base import BackendError
from ..recognition.registry import installed_vision_backends

#: About how long loading DECIMER takes on a CPU, for the progress display.
TYPICAL_SECONDS = 40

STARTED = time.time()
_STATE = {"state": "waiting", "backend": None, "since": None, "finished": None, "error": None}
_LOCK = threading.Lock()


def start() -> None:
    with _LOCK:
        if _STATE["state"] != "waiting":
            return
        _STATE["state"] = "loading"
        _STATE["since"] = time.time()
    threading.Thread(target=_load, name="m2i-warmup", daemon=True).start()


def _load() -> None:
    backends = installed_vision_backends()
    if not backends:
        _set(state="absent", finished=time.time())
        return
    backend = backends[0]
    _set(backend=backend.name)
    try:
        _subprocess.persistent_worker(backend.name, backend.spec.worker).ask(
            {"mode": "warmup"}, _subprocess.DEFAULT_TIMEOUT
        )
    except BackendError as exc:
        _set(state="failed", error=str(exc), finished=time.time())
        return
    _set(state="ready", finished=time.time())


def _set(**values) -> None:
    with _LOCK:
        _STATE.update(values)


def status() -> dict:
    with _LOCK:
        state = dict(_STATE)
    now = time.time()
    model = state["state"]
    elapsed = round((state["finished"] or now) - state["since"]) if state["since"] else 0
    return {
        # The site is usable either way; "starting" only means photos wait.
        "status": "starting" if model in ("waiting", "loading") else "ready",
        "uptime": round(now - STARTED),
        "model": {
            "state": model,
            "backend": state["backend"],
            "elapsed": elapsed,
            "typical": TYPICAL_SECONDS,
            "error": state["error"],
        },
        "routes_available": ["smiles", "file", "cif"] + (["picture"] if model == "ready" else []),
    }
