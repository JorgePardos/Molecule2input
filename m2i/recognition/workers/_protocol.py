"""Shared request/response handling for the backend workers.

Copied into each venv is not necessary -- the workers are executed by path and
import this module as a sibling, which works because the whole directory is on
sys.path when a script inside it is run.
"""

from __future__ import annotations

import json
import sys
import traceback
from pathlib import Path


def serve(handler) -> int:
    """Read the request file named in argv[1], run `handler`, write the response.

    Never raises: a crash inside the model must come back as a readable error
    rather than as a truncated stream the caller has to guess about.
    """
    if len(sys.argv) < 2:
        print("usage: worker.py <request.json>", file=sys.stderr)
        return 2

    try:
        request = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    except Exception as exc:  # noqa: BLE001 - the caller only sees the message
        print(f"unreadable request: {exc}", file=sys.stderr)
        return 2

    response_path = request.get("response_path")
    try:
        payload = handler(request) or {}
        payload["ok"] = True
    except BaseException as exc:  # noqa: BLE001 - includes SystemExit from models
        payload = {
            "ok": False,
            "error": f"{type(exc).__name__}: {exc}",
            "traceback": traceback.format_exc(limit=8),
        }

    if response_path:
        try:
            Path(response_path).write_text(
                json.dumps(payload, default=str), encoding="utf-8"
            )
        except Exception as exc:  # noqa: BLE001
            print(f"could not write response: {exc}", file=sys.stderr)
            return 3
    return 0 if payload.get("ok") else 1


def package_versions(*names: str) -> dict:
    """Best-effort version report, for the install marker and for debugging."""
    versions = {}
    for name in names:
        try:
            from importlib.metadata import version

            versions[name] = version(name)
        except Exception:  # noqa: BLE001
            versions[name] = "unknown"
    return versions
