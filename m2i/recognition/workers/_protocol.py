"""Shared request/response handling for the backend workers.

Copied into each venv is not necessary -- the workers are executed by path and
import this module as a sibling, which works because the whole directory is on
sys.path when a script inside it is run.
"""

from __future__ import annotations

import json
import os
import sys
import traceback
from pathlib import Path

#: Passed instead of a request file: keep the model loaded and take requests
#: one per line on stdin until it closes.
SERVE_FLAG = "--serve"


def serve(handler) -> int:
    """Answer one request file named in argv[1], or many with ``--serve``.

    Never raises: a crash inside the model must come back as a readable error
    rather than as a truncated stream the caller has to guess about.
    """
    if len(sys.argv) < 2:
        print("usage: worker.py <request.json> | --serve", file=sys.stderr)
        return 2
    if sys.argv[1] == SERVE_FLAG:
        return _serve_forever(handler)

    try:
        request = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    except Exception as exc:  # noqa: BLE001 - the caller only sees the message
        print(f"unreadable request: {exc}", file=sys.stderr)
        return 2
    payload = _answer(handler, request)
    if not _write(request.get("response_path"), payload):
        return 3
    return 0 if payload.get("ok") else 1


def _serve_forever(handler) -> int:
    """Loading a model takes most of a minute; answering takes seconds.

    So the process stays up, reading one JSON request per line. The payload
    still travels in a file -- stdout belongs to the model's chatter -- and is
    renamed into place only when complete, so the caller never reads half.
    """
    while True:
        line = sys.stdin.readline()
        if not line:
            return 0  # the caller closed the pipe: done
        if not line.strip():
            continue
        try:
            request = json.loads(line)
        except json.JSONDecodeError as exc:
            print(f"unreadable request: {exc}", file=sys.stderr)
            continue
        _write(request.get("response_path"), _answer(handler, request))


def _answer(handler, request: dict) -> dict:
    try:
        payload = handler(request) or {}
        payload["ok"] = True
    except BaseException as exc:  # noqa: BLE001 - includes SystemExit from models
        payload = {
            "ok": False,
            "error": f"{type(exc).__name__}: {exc}",
            "traceback": traceback.format_exc(limit=8),
        }
    return payload


def _write(response_path, payload: dict) -> bool:
    if not response_path:
        return True
    try:
        partial = Path(f"{response_path}.part")
        partial.write_text(json.dumps(payload, default=str), encoding="utf-8")
        os.replace(partial, response_path)
    except Exception as exc:  # noqa: BLE001
        print(f"could not write response: {exc}", file=sys.stderr)
        return False
    return True


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
