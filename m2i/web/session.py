"""Per-visitor state, kept on the server and found by a signed cookie.

What ``st.session_state`` held in the Streamlit app: the visitor's private
folder, what they gave (by kind, so a picture never leaks into the SMILES box),
corrections, confirmations, and the caches that make switching format cheap.

Everything lives in memory: a restart forgets the sessions, and their folders
are swept a day later like any other. Each session has its own lock, so two
requests from the same page cannot interleave on its state.
"""

from __future__ import annotations

import os
import secrets
import threading
import time
from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path

from fastapi import Request, Response
from itsdangerous import BadSignature, URLSafeSerializer

from ..gui import hosting

COOKIE = "m2i_session"
#: Kept in memory this long after the last request; the folder follows a day
#: later through hosting.sweep_sessions.
IDLE_LIMIT = hosting.SESSION_MAX_AGE

#: A fixed secret keeps sessions valid across restarts of several workers;
#: without one, a random secret per process is just as safe for one process.
_SIGNER = URLSafeSerializer(os.environ.get("M2I_SECRET") or secrets.token_hex(32), salt="m2i")


class LimitedCache(OrderedDict):
    """A dict that forgets its oldest entries past ``size``."""

    def __init__(self, size: int) -> None:
        super().__init__()
        self.size = size

    def __setitem__(self, key, value) -> None:
        super().__setitem__(key, value)
        self.move_to_end(key)
        while len(self) > self.size:
            self.popitem(last=False)


@dataclass
class Source:
    """Something a visitor gave: a picture, a file or a SMILES."""

    kind: str  # picture | file | smiles
    key: str
    path: Path | None = None
    image: object = None  # the prepared PIL image of a picture
    style: str | None = None
    shown_name: str | None = None
    recognition: object = None  # RecognitionResult, once read
    issues: list = field(default_factory=list)
    display_smiles: str = ""


@dataclass
class Session:
    id: str
    folder: Path
    lock: threading.RLock = field(default_factory=threading.RLock)
    last_seen: float = field(default_factory=time.time)
    sources: dict[str, Source] = field(default_factory=dict)  # by key
    corrections: dict = field(default_factory=dict)  # key -> RecognitionResult
    confirmed: set = field(default_factory=set)  # (key, smiles)
    embeddings: LimitedCache = field(default_factory=lambda: LimitedCache(6))
    crystals: LimitedCache = field(default_factory=lambda: LimitedCache(4))
    om_drawings: LimitedCache = field(default_factory=lambda: LimitedCache(4))
    om_built: LimitedCache = field(default_factory=lambda: LimitedCache(4))
    files: dict[str, Path] = field(default_factory=dict)  # downloadable, by name

    @property
    def output(self) -> Path:
        folder = self.folder / "output"
        folder.mkdir(parents=True, exist_ok=True)
        return folder

    def touch(self) -> None:
        self.last_seen = time.time()
        # Swept away while its visitor was gone: start it again empty.
        self.folder.mkdir(parents=True, exist_ok=True)

    def offer(self, paths) -> list[dict]:
        """Make written files downloadable by this session, and describe them."""
        out = []
        for path in map(Path, paths):
            self.files[path.name] = path
            out.append({"name": path.name, "download_url": f"/api/download/{path.name}"})
        return out


_SESSIONS: dict[str, Session] = {}
_LOCK = threading.Lock()


def current(request: Request, response: Response) -> Session:
    """The visitor's session, created (with its cookie) on the first request."""
    session_id = None
    token = request.cookies.get(COOKIE)
    if token:
        try:
            session_id = _SIGNER.loads(token)
        except BadSignature:
            session_id = None
    with _LOCK:
        _forget_idle()
        session = _SESSIONS.get(session_id) if session_id else None
        if session is None:
            session = Session(id=secrets.token_urlsafe(18), folder=hosting.new_session_folder())
            _SESSIONS[session.id] = session
            response.set_cookie(
                COOKIE, _SIGNER.dumps(session.id), httponly=True, samesite="lax",
                secure=request.url.scheme == "https", max_age=int(IDLE_LIMIT),
            )
    session.touch()
    return session


def _forget_idle() -> None:
    limit = time.time() - IDLE_LIMIT
    for key in [k for k, s in _SESSIONS.items() if s.last_seen < limit]:
        del _SESSIONS[key]
