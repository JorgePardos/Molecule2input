"""Whether the interface is served to other people, and what that takes.

Set ``M2I_HOSTED=1`` (the Docker image does). A few things only make sense
locally and are left out when hosted: a folder of the server to save into,
where the files were written on its disk, and ``m2i setup`` instructions
nobody visiting the page can follow.

A server also runs for months, so what each visitor leaves behind -- their
uploads and generated files, in a folder of their own -- has to go away on
its own. Streamlit gives no signal when a session ends, so folders nobody has
touched for a day are removed whenever a new session starts.
"""

from __future__ import annotations

import os
import shutil
import tempfile
import time
from pathlib import Path

HOSTED = os.environ.get("M2I_HOSTED", "").strip().lower() in ("1", "true", "yes")

SESSION_PREFIX = "m2i-session-"
#: A session folder untouched for this long belongs to someone who has left.
SESSION_MAX_AGE = 24 * 3600


def new_session_folder() -> Path:
    sweep_sessions()
    return Path(tempfile.mkdtemp(prefix=SESSION_PREFIX))


def sweep_sessions(root: Path | None = None, max_age: float = SESSION_MAX_AGE, now: float | None = None) -> int:
    """Remove session folders nobody has written to for ``max_age`` seconds."""
    root = Path(root or tempfile.gettempdir())
    now = time.time() if now is None else now
    removed = 0
    for folder in root.glob(f"{SESSION_PREFIX}*"):
        try:
            if folder.is_dir() and now - _last_written(folder) > max_age:
                shutil.rmtree(folder, ignore_errors=True)
                removed += 1
        except OSError:
            continue  # another session removing the same folder: fine
    return removed


def _last_written(folder: Path) -> float:
    """The newest modification time in the folder, however deep."""
    newest = folder.stat().st_mtime
    for path in folder.rglob("*"):
        newest = max(newest, path.stat().st_mtime)
    return newest
