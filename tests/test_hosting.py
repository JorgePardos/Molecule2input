"""A server runs for months: what visitors leave behind has to go away on its own."""

from __future__ import annotations

import os

from m2i.gui import hosting


def make_session(root, name, age, now):
    folder = root / f"{hosting.SESSION_PREFIX}{name}"
    (folder / "output").mkdir(parents=True)
    written = folder / "output" / "molecule.gjf"
    written.write_text("x")
    for path in (folder, folder / "output", written):
        os.utime(path, (now - age, now - age))
    return folder


def test_sessions_untouched_for_a_day_are_removed(tmp_path):
    now = 1_000_000_000.0
    old = make_session(tmp_path, "old", age=2 * 24 * 3600, now=now)
    recent = make_session(tmp_path, "recent", age=3600, now=now)
    unrelated = tmp_path / "somebody-elses-folder"
    unrelated.mkdir()
    os.utime(unrelated, (0, 0))

    assert hosting.sweep_sessions(tmp_path, now=now) == 1
    assert not old.exists()
    assert recent.exists()
    assert unrelated.exists()


def test_a_recent_file_deep_inside_keeps_the_session(tmp_path):
    now = 1_000_000_000.0
    folder = make_session(tmp_path, "busy", age=3 * 24 * 3600, now=now)
    os.utime(folder / "output" / "molecule.gjf", (now - 60, now - 60))
    assert hosting.sweep_sessions(tmp_path, now=now) == 0
    assert folder.exists()
