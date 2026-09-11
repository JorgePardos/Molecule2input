"""The backend installer: specs, status and removal.

Nothing here installs anything -- the real installation is exercised by hand,
and by `m2i doctor` afterwards.
"""

from __future__ import annotations

import json

import pytest

from m2i import backends
from m2i.recognition import _subprocess


@pytest.fixture(autouse=True)
def isolated_backend_home(monkeypatch, tmp_path):
    """Never touch the user's real ~/.m2i while testing."""
    monkeypatch.setenv("M2I_BACKEND_HOME", str(tmp_path))
    return tmp_path


def test_every_spec_is_complete():
    assert set(backends.SPECS) == {"molscribe", "decimer"}
    for spec in backends.SPECS.values():
        assert spec.install_steps and all(spec.install_steps)
        assert spec.download_size and spec.strength and spec.description
        assert backends.worker_exists(spec), f"{spec.name} has no worker script"


def test_exactly_one_backend_returns_a_molblock():
    """The asymmetry is the point: one measures stereochemistry, one writes it."""
    graph = [s for s in backends.SPECS.values() if s.returns_molblock]
    sequence = [s for s in backends.SPECS.values() if not s.returns_molblock]
    assert [s.name for s in graph] == ["molscribe"]
    assert [s.name for s in sequence] == ["decimer"]


def test_molscribe_pins_are_repeated_in_every_step():
    """albumentations pulls a newer NumPy and the headless opencv; if the pins
    are not repeated, the second step silently undoes the first."""
    steps = backends.SPECS["molscribe"].install_steps
    numpy_steps = [s for s in steps if any(p.startswith("numpy") for p in s)]
    assert len(numpy_steps) >= 2

    opencv_step = next(s for s in steps if any("albumentations" in p for p in s))
    assert any(p.startswith("opencv-python==") for p in opencv_step)
    assert any(p.startswith("opencv-python-headless==") for p in opencv_step)


def test_status_of_a_backend_that_is_not_installed():
    state = backends.status("decimer")
    assert state["known"]
    assert not state["installed"]
    assert not state["python_exists"]
    assert state["marker"] is None


def test_status_of_an_unknown_backend():
    assert backends.status("nonesuch")["known"] is False


def test_status_needs_the_marker_not_just_the_venv(isolated_backend_home):
    python = _subprocess.venv_python("decimer")
    python.parent.mkdir(parents=True)
    python.write_text("")
    assert not backends.status("decimer")["installed"]

    _subprocess.marker_path("decimer").write_text(json.dumps({"name": "decimer"}))
    assert backends.status("decimer")["installed"]


def test_installing_an_unknown_backend_is_refused():
    with pytest.raises(backends.SetupError) as excinfo:
        backends.install("nonesuch")
    assert "molscribe" in str(excinfo.value)


def test_removing_something_that_is_not_there_is_not_an_error():
    assert backends.remove("decimer") is False


def test_remove_deletes_only_that_backend(isolated_backend_home):
    for name in ("decimer", "molscribe"):
        (isolated_backend_home / name).mkdir(parents=True)
        (isolated_backend_home / name / "file").write_text("x")

    assert backends.remove("decimer") is True
    assert not (isolated_backend_home / "decimer").exists()
    assert (isolated_backend_home / "molscribe").exists()


def test_describe_tells_the_user_what_they_are_choosing():
    text = backends.describe(backends.SPECS["decimer"])
    assert "hand-drawn" in text
    assert "GB" in text  # the download size is stated before anything is fetched


def test_all_status_covers_every_spec():
    assert {row["name"] for row in backends.all_status()} == set(backends.SPECS)
