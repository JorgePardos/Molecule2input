"""Bridge to a recognition model living in its own virtual environment.

The vision backends cannot share an environment with m2i or with each other:
MolScribe pins ``numpy<2``, DECIMER pulls in TensorFlow, and m2i itself runs on
a current RDKit and NumPy 2.x. Each backend therefore gets its own venv and is
driven as a subprocess.

The protocol deliberately does **not** use stdout for the payload. These models
print progress bars, TensorFlow banners and CUDA warnings; anything sharing the
channel with the result would corrupt it. The request and the response are
files, and stdout/stderr are captured only to explain failures.
"""

from __future__ import annotations

import atexit
import json
import os
import subprocess
import sys
import tempfile
import threading
import time
import uuid
from pathlib import Path

from ..types import RecognitionResult
from .base import BackendError

#: Long enough for a cold start that loads a large checkpoint from disk.
DEFAULT_TIMEOUT = 600
WORKER_DIR = Path(__file__).parent / "workers"


def backend_home() -> Path:
    """Where the per-backend virtual environments live."""
    override = os.environ.get("M2I_BACKEND_HOME")
    if override:
        return Path(override)
    return Path.home() / ".m2i" / "backends"


def venv_dir(name: str) -> Path:
    return backend_home() / name


def venv_python(name: str) -> Path:
    """Interpreter inside the backend's venv."""
    root = venv_dir(name)
    if sys.platform == "win32":
        return root / "Scripts" / "python.exe"
    return root / "bin" / "python"


def marker_path(name: str) -> Path:
    """Written only after a successful install; its absence means 'not ready'."""
    return venv_dir(name) / "m2i-backend.json"


def read_marker(name: str) -> dict | None:
    path = marker_path(name)
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def is_installed(name: str) -> bool:
    return venv_python(name).is_file() and read_marker(name) is not None


def worker_path(filename: str | Path) -> Path:
    """Resolve a worker name against the bundled directory, or accept a path."""
    candidate = Path(filename)
    if candidate.is_absolute():
        if not candidate.is_file():
            raise BackendError(f"worker script missing: {candidate}")
        return candidate
    path = WORKER_DIR / filename
    if not path.is_file():
        raise BackendError(f"worker script missing: {path}")
    return path


def run_worker(
    name: str,
    worker: str,
    request: dict,
    *,
    timeout: int = DEFAULT_TIMEOUT,
    python: Path | None = None,
) -> dict:
    """Run a worker in the backend's venv and return its JSON response."""
    interpreter = python or venv_python(name)
    if not interpreter.is_file():
        raise BackendError(
            f"{name} is not installed (no interpreter at {interpreter}). "
            f"Run: m2i setup {name}"
        )

    with tempfile.TemporaryDirectory(prefix=f"m2i-{name}-") as workspace:
        request_file = Path(workspace) / "request.json"
        response_file = Path(workspace) / "response.json"
        request = dict(request, response_path=str(response_file))
        request_file.write_text(json.dumps(request), encoding="utf-8")

        env = _worker_env()

        try:
            completed = subprocess.run(
                [str(interpreter), str(worker_path(worker)), str(request_file)],
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=timeout,
                env=env,
            )
        except subprocess.TimeoutExpired as exc:
            raise BackendError(
                f"{name} timed out after {timeout}s. A first run downloads the model "
                "weights and can be slow; try again, or raise the timeout."
            ) from exc

        if not response_file.is_file():
            raise BackendError(
                f"{name} produced no result (exit code {completed.returncode}).\n"
                f"{_tail(completed.stderr) or _tail(completed.stdout)}"
            )

        try:
            payload = json.loads(response_file.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            raise BackendError(f"{name} returned a malformed response: {exc}") from exc

    if not payload.get("ok"):
        raise BackendError(f"{name}: {payload.get('error', 'unknown failure')}")
    return payload


class PersistentWorker:
    """One model kept loaded in its own process, answering requests in turn.

    Starting a worker per image re-imports TensorFlow and reloads the weights
    every time: most of a minute for DECIMER, for a prediction that takes a
    few seconds. A server makes that unbearable, so the process stays up.
    Requests go in one JSON line at a time on stdin; answers come back as
    files, as in a one-shot run. A lock serialises callers -- the models are
    not thread-safe, and on a two-core machine running them side by side would
    not be faster anyway.

    The model's output goes to a log file rather than a pipe: a pipe nobody
    reads fills up, and the worker would then block forever mid-banner.
    """

    def __init__(self, name: str, interpreter: Path, script: Path) -> None:
        self.name = name
        self.closed = False
        self.lock = threading.Lock()
        self.workspace = Path(tempfile.mkdtemp(prefix=f"m2i-{name}-"))
        self.log_path = self.workspace / "worker.log"
        self._log = open(self.log_path, "w", encoding="utf-8", errors="replace")  # noqa: SIM115
        self.process = subprocess.Popen(
            [str(interpreter), str(script), "--serve"],
            stdin=subprocess.PIPE,
            stdout=self._log,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            env=_worker_env(),
        )

    def alive(self) -> bool:
        return not self.closed and self.process.poll() is None

    def ask(self, request: dict, timeout: int) -> dict:
        with self.lock:
            response = self.workspace / f"{uuid.uuid4().hex}.json"
            try:
                self.process.stdin.write(json.dumps(dict(request, response_path=str(response))) + "\n")
                self.process.stdin.flush()
            except (BrokenPipeError, OSError, ValueError) as exc:
                raise BackendError(f"{self.name} stopped unexpectedly.\n{self._tail()}") from exc
            deadline = time.monotonic() + timeout
            while not response.is_file():
                if not self.alive():
                    raise BackendError(
                        f"{self.name} produced no result (exit code {self.process.returncode}).\n"
                        f"{self._tail()}"
                    )
                if time.monotonic() > deadline:
                    self.close()  # a stuck model is not given the next request too
                    raise BackendError(
                        f"{self.name} timed out after {timeout}s. A first run downloads the "
                        "model weights and can be slow; try again, or raise the timeout."
                    )
                time.sleep(0.05)
            try:
                payload = json.loads(response.read_text(encoding="utf-8"))
            except json.JSONDecodeError as exc:
                raise BackendError(f"{self.name} returned a malformed response: {exc}") from exc
            finally:
                response.unlink(missing_ok=True)
        if not payload.get("ok"):
            raise BackendError(f"{self.name}: {payload.get('error', 'unknown failure')}")
        return payload

    def close(self) -> None:
        if self.closed:
            return
        self.closed = True
        try:
            self.process.stdin.close()  # a healthy worker exits on end of input
            self.process.wait(timeout=5)
        except Exception:  # noqa: BLE001 - a stuck one is stopped regardless
            self.process.kill()
            self.process.wait()
        self._log.close()

    def _tail(self) -> str:
        try:
            if not self._log.closed:
                self._log.flush()
            return _tail(self.log_path.read_text(encoding="utf-8", errors="replace"))
        except (OSError, ValueError):
            return ""


_WORKERS: dict[tuple[str, str], PersistentWorker] = {}
_WORKERS_LOCK = threading.Lock()


def persistent_worker(name: str, worker: str | Path, python: Path | None = None) -> PersistentWorker:
    """The running worker for this model, started (or restarted) if needed."""
    interpreter = python or venv_python(name)
    if not interpreter.is_file():
        raise BackendError(
            f"{name} is not installed (no interpreter at {interpreter}). "
            f"Run: m2i setup {name}"
        )
    script = worker_path(worker)
    key = (str(interpreter), str(script))
    with _WORKERS_LOCK:
        current = _WORKERS.get(key)
        if current is None or not current.alive():
            if current is not None:
                current.close()
            current = _WORKERS[key] = PersistentWorker(name, interpreter, script)
        return current


@atexit.register
def stop_workers() -> None:
    with _WORKERS_LOCK:
        for worker in _WORKERS.values():
            worker.close()
        _WORKERS.clear()


def _worker_env() -> dict:
    # UTF-8 so error text survives the trip, and nothing imported from the
    # caller's environment.
    env = dict(os.environ)
    env["PYTHONIOENCODING"] = "utf-8"
    env["PYTHONNOUSERSITE"] = "1"
    env.pop("PYTHONPATH", None)
    env.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")
    env.setdefault("TRANSFORMERS_VERBOSITY", "error")
    return env


def recognize_with_worker(
    name: str,
    worker: str | Path,
    image_path: Path,
    request: dict,
    *,
    timeout: int,
    python: Path | None = None,
) -> RecognitionResult:
    payload = persistent_worker(name, worker, python).ask(
        dict(request, image_path=str(Path(image_path).resolve())), timeout
    )
    smiles = payload.get("smiles") or ""
    if not smiles and not payload.get("molblock"):
        raise BackendError(f"{name} returned an empty structure")
    return RecognitionResult(
        smiles=smiles,
        molblock=payload.get("molblock"),
        confidence=payload.get("confidence"),
        backend=name,
        raw=payload.get("raw"),
    )


def _tail(text: str | None, lines: int = 12) -> str:
    if not text:
        return ""
    return "\n".join(text.strip().splitlines()[-lines:])
