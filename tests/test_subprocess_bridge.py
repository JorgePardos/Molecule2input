"""The venv bridge. Uses throwaway worker scripts run by the current
interpreter, so none of this needs a model installed."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from m2i.recognition import _subprocess
from m2i.recognition.base import BackendError

HERE = Path(_subprocess.__file__).parent


def write_worker(tmp_path: Path, body: str) -> Path:
    """A worker that uses the real protocol helper, like the shipped ones do."""
    script = tmp_path / "fake_worker.py"
    script.write_text(
        "import sys\n"
        f"sys.path.insert(0, {str(HERE / 'workers')!r})\n"
        "from _protocol import serve\n"
        f"{body}\n"
        "raise SystemExit(serve(handle))\n",
        encoding="utf-8",
    )
    return script


def run(script: Path, request: dict, **kwargs):
    return _subprocess.run_worker(
        "fake", script, request, python=Path(sys.executable), **kwargs
    )


def test_round_trip(tmp_path):
    script = write_worker(
        tmp_path,
        "def handle(request):\n"
        "    return {'smiles': 'CCO', 'echo': request.get('marker')}\n",
    )
    payload = run(script, {"mode": "recognize", "marker": 42})
    assert payload["ok"] is True
    assert payload["smiles"] == "CCO"
    assert payload["echo"] == 42


def test_model_chatter_on_stdout_does_not_corrupt_the_result(tmp_path):
    """The whole reason the payload travels in a file: these models print
    progress bars, TensorFlow banners and CUDA warnings."""
    script = write_worker(
        tmp_path,
        "def handle(request):\n"
        "    print('2026-07-31 oneDNN custom operations are on')\n"
        "    print('{\"not\": \"the response\"}')\n"
        "    import sys; print('W tensorflow/core: something', file=sys.stderr)\n"
        "    return {'smiles': 'c1ccccc1'}\n",
    )
    assert run(script, {})["smiles"] == "c1ccccc1"


def test_worker_exception_becomes_a_readable_error(tmp_path):
    script = write_worker(
        tmp_path,
        "def handle(request):\n    raise ValueError('checkpoint is corrupt')\n",
    )
    with pytest.raises(BackendError) as excinfo:
        run(script, {})
    assert "checkpoint is corrupt" in str(excinfo.value)
    assert "ValueError" in str(excinfo.value)


def test_worker_that_dies_without_responding_is_reported(tmp_path):
    script = tmp_path / "crash.py"
    script.write_text("import os\nos._exit(3)\n", encoding="utf-8")
    with pytest.raises(BackendError) as excinfo:
        run(script, {})
    assert "no result" in str(excinfo.value)


def test_timeout_explains_the_first_run_cost(tmp_path):
    script = tmp_path / "slow.py"
    script.write_text("import time\ntime.sleep(30)\n", encoding="utf-8")
    with pytest.raises(BackendError) as excinfo:
        run(script, {}, timeout=1)
    assert "timed out" in str(excinfo.value)
    assert "weights" in str(excinfo.value)


def test_missing_interpreter_points_at_setup(tmp_path):
    with pytest.raises(BackendError) as excinfo:
        _subprocess.run_worker(
            "decimer", "decimer_worker.py", {}, python=tmp_path / "nope.exe"
        )
    assert "m2i setup decimer" in str(excinfo.value)


def test_missing_worker_script_is_reported():
    with pytest.raises(BackendError):
        _subprocess.worker_path("does_not_exist.py")


def test_recognize_wrapper_builds_a_recognition_result(tmp_path):
    script = write_worker(
        tmp_path,
        "def handle(request):\n"
        "    return {'smiles': 'CCO', 'molblock': 'block', 'confidence': 0.9,\n"
        "            'raw': {'n_atoms': 3}}\n",
    )
    result = _subprocess.recognize_with_worker(
        "fake", script, tmp_path, {"mode": "recognize"}, timeout=60,
        python=Path(sys.executable),
    )
    assert result.smiles == "CCO"
    assert result.molblock == "block"
    assert result.confidence == 0.9
    assert result.backend == "fake"
    assert result.raw == {"n_atoms": 3}


def test_empty_structure_is_rejected(tmp_path):
    script = write_worker(tmp_path, "def handle(request):\n    return {'smiles': ''}\n")
    with pytest.raises(BackendError) as excinfo:
        _subprocess.recognize_with_worker(
            "fake", script, tmp_path, {}, timeout=60, python=Path(sys.executable)
        )
    assert "empty" in str(excinfo.value)


def test_image_path_is_made_absolute(tmp_path):
    script = write_worker(
        tmp_path,
        "def handle(request):\n"
        "    return {'smiles': 'C', 'raw': {'path': request['image_path']}}\n",
    )
    result = _subprocess.recognize_with_worker(
        "fake", script, Path("some/relative/drawing.png"), {}, timeout=60,
        python=Path(sys.executable),
    )
    assert Path(result.raw["path"]).is_absolute()


def test_backend_home_is_configurable(monkeypatch, tmp_path):
    monkeypatch.setenv("M2I_BACKEND_HOME", str(tmp_path))
    assert _subprocess.backend_home() == tmp_path
    assert _subprocess.venv_dir("decimer") == tmp_path / "decimer"
    assert _subprocess.venv_python("decimer").parent.parent == tmp_path / "decimer"


def test_installation_requires_both_interpreter_and_marker(monkeypatch, tmp_path):
    monkeypatch.setenv("M2I_BACKEND_HOME", str(tmp_path))
    assert not _subprocess.is_installed("decimer")

    python = _subprocess.venv_python("decimer")
    python.parent.mkdir(parents=True)
    python.write_text("")
    assert not _subprocess.is_installed("decimer")  # interpreter but no marker

    _subprocess.marker_path("decimer").write_text(json.dumps({"name": "decimer"}))
    assert _subprocess.is_installed("decimer")


def test_corrupt_marker_is_treated_as_not_installed(monkeypatch, tmp_path):
    monkeypatch.setenv("M2I_BACKEND_HOME", str(tmp_path))
    python = _subprocess.venv_python("decimer")
    python.parent.mkdir(parents=True)
    python.write_text("")
    _subprocess.marker_path("decimer").write_text("{not json")
    assert not _subprocess.is_installed("decimer")


# -- a worker kept running between images ---------------------------------------


def ask(script: Path, request: dict, timeout: int = 60) -> dict:
    return _subprocess.persistent_worker("fake", script, Path(sys.executable)).ask(request, timeout)


def test_the_model_stays_loaded_between_requests(tmp_path):
    """Reloading DECIMER costs most of a minute per image; answering costs seconds."""
    script = write_worker(
        tmp_path,
        "import os\n"
        "LOADS = []\n"
        "def handle(request):\n"
        "    LOADS.append(1)\n"
        "    return {'pid': os.getpid(), 'requests_seen': len(LOADS)}\n",
    )
    first, second = ask(script, {}), ask(script, {})
    assert first["pid"] == second["pid"]
    assert second["requests_seen"] == 2
    _subprocess.stop_workers()


def test_chatter_and_errors_do_not_end_the_session(tmp_path):
    script = write_worker(
        tmp_path,
        "def handle(request):\n"
        "    print('W tensorflow/core: banner ' * 5000)\n"  # more than a pipe buffer
        "    if request.get('fail'):\n"
        "        raise ValueError('unreadable picture')\n"
        "    return {'smiles': 'CCO'}\n",
    )
    with pytest.raises(BackendError, match="unreadable picture"):
        ask(script, {"fail": True})
    assert ask(script, {})["smiles"] == "CCO"
    _subprocess.stop_workers()


def test_a_worker_that_died_is_started_again(tmp_path):
    script = write_worker(
        tmp_path,
        "import os\n"
        "def handle(request):\n"
        "    if request.get('crash'):\n"
        "        os._exit(9)\n"
        "    return {'pid': os.getpid()}\n",
    )
    before = ask(script, {})["pid"]
    with pytest.raises(BackendError, match="no result"):
        ask(script, {"crash": True})
    assert ask(script, {})["pid"] != before
    _subprocess.stop_workers()


def test_a_stuck_worker_times_out_and_is_replaced(tmp_path):
    script = write_worker(
        tmp_path,
        "import os, time\n"
        "def handle(request):\n"
        "    if request.get('hang'):\n"
        "        time.sleep(60)\n"
        "    return {'pid': os.getpid()}\n",
    )
    before = ask(script, {})["pid"]
    with pytest.raises(BackendError, match="timed out"):
        ask(script, {"hang": True}, timeout=1)
    assert ask(script, {})["pid"] != before
    _subprocess.stop_workers()


def test_concurrent_callers_are_answered_in_turn(tmp_path):
    from concurrent.futures import ThreadPoolExecutor

    script = write_worker(
        tmp_path,
        "def handle(request):\n"
        "    return {'echo': request['n']}\n",
    )
    with ThreadPoolExecutor(4) as pool:
        answers = list(pool.map(lambda n: ask(script, {"n": n})["echo"], range(12)))
    assert answers == list(range(12))
    _subprocess.stop_workers()
