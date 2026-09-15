"""Installing the vision backends into their own virtual environments.

The recognition model (DECIMER, for hand-drawn structures) brings TensorFlow,
which cannot share an environment with m2i itself, so `m2i setup <name>` builds a dedicated venv,
installs a curated requirement list into it, downloads the model weights, and
writes a marker file recording exactly what landed there.

The requirement lists are curated rather than copied from each project's
``requirements.txt``: those pin git commits and versions that no longer build
on current Python. Where a pin was replaced, the reason is in the comment.
"""

from __future__ import annotations

import json
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from .recognition._subprocess import (
    marker_path,
    read_marker,
    venv_dir,
    venv_python,
    worker_path,
)

#: pip can take a long time on a cold cache; tensorflow is large.
INSTALL_TIMEOUT = 3600
WARMUP_TIMEOUT = 1800


class SetupError(RuntimeError):
    """The backend could not be installed."""


@dataclass(frozen=True)
class BackendSpec:
    name: str
    description: str
    returns_molblock: bool
    worker: str
    #: pip invocations run in order; split so a resolver conflict is localised.
    install_steps: tuple[tuple[str, ...], ...]
    #: Approximate download size, so the user can decide before starting.
    download_size: str
    #: What this backend is actually good at, shown by `m2i doctor`.
    strength: str
    warmup_request: dict = field(default_factory=lambda: {"mode": "warmup"})
    notes: str = ""


DECIMER = BackendSpec(
    name="decimer",
    description=(
        "image-to-sequence; a transformer decoder that writes the SMILES "
        "directly, with a model fine-tuned on hand-drawn structures"
    ),
    strength="hand-drawn structures, where it is far ahead of every alternative",
    returns_molblock=False,
    worker="decimer_worker.py",
    download_size="about 1.5 GB (tensorflow, then the model weights)",
    install_steps=(("decimer",),),
    notes=(
        "Returns only a SMILES: its stereochemistry comes from the generated "
        "token sequence, not from measuring the drawing."
    ),
)

SPECS: dict[str, BackendSpec] = {DECIMER.name: DECIMER}


# -- status --------------------------------------------------------------


def status(name: str) -> dict:
    spec = SPECS.get(name)
    marker = read_marker(name)
    python = venv_python(name)
    return {
        "name": name,
        "known": spec is not None,
        "venv": str(venv_dir(name)),
        "python_exists": python.is_file(),
        "installed": python.is_file() and marker is not None,
        "marker": marker,
    }


def all_status() -> list[dict]:
    return [status(name) for name in SPECS]


# -- installation --------------------------------------------------------


def install(
    name: str,
    *,
    base_python: Path | None = None,
    force: bool = False,
    on_output=print,
) -> dict:
    """Create the venv, install the requirements, warm the model up."""
    spec = SPECS.get(name)
    if spec is None:
        raise SetupError(f"unknown backend {name!r}; known: {', '.join(SPECS)}")

    root = venv_dir(name)
    python = venv_python(name)

    if python.is_file() and not force:
        if read_marker(name):
            on_output(f"{name} is already installed at {root}")
            return read_marker(name)
        on_output(f"Reusing the existing environment at {root}")
    else:
        if force and root.exists():
            on_output(f"Removing {root}")
            remove(name)
        _create_venv(name, root, base_python or Path(sys.executable), on_output)

    _pip(python, ["install", "--upgrade", "pip", "setuptools", "wheel"], on_output)

    for step, packages in enumerate(spec.install_steps, start=1):
        on_output(f"\n[{step}/{len(spec.install_steps)}] installing {' '.join(packages)}")
        _pip(python, ["install", *packages], on_output)

    on_output("\nDownloading model weights (this is the slow part)...")
    warmup = _warmup(spec, on_output)

    marker = {
        "name": name,
        "installed_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python": str(python),
        "base_python": str(base_python or sys.executable),
        "install_steps": [list(s) for s in spec.install_steps],
        "warmup": warmup,
    }
    marker_path(name).write_text(json.dumps(marker, indent=2), encoding="utf-8")
    on_output(f"\n{name} is ready. Try: m2i from-image <drawing.png> --backend {name}")
    return marker


def remove(name: str) -> bool:
    """Delete a backend's environment. Nothing else on the system is touched."""
    import shutil

    root = venv_dir(name)
    if not root.exists():
        return False
    shutil.rmtree(root, ignore_errors=True)
    return not root.exists()


# -- internals -----------------------------------------------------------


def _create_venv(name: str, root: Path, base_python: Path, on_output) -> None:
    on_output(f"Creating a virtual environment at {root}")
    root.parent.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        [str(base_python), "-m", "venv", str(root)],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    if result.returncode != 0 or not venv_python(name).is_file():
        raise SetupError(
            f"could not create the virtual environment at {root}:\n"
            f"{result.stderr.strip() or result.stdout.strip()}"
        )


def _pip(python: Path, arguments: list[str], on_output) -> None:
    command = [str(python), "-m", "pip", *arguments, "--disable-pip-version-check"]
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
        bufsize=1,
    )
    tail: list[str] = []
    assert process.stdout is not None
    for line in process.stdout:
        line = line.rstrip()
        tail.append(line)
        del tail[:-40]
        if _worth_showing(line):
            on_output(f"    {line}")
    if process.wait(timeout=INSTALL_TIMEOUT) != 0:
        raise SetupError(
            "pip failed:\n  " + "\n  ".join(tail[-20:]) + f"\n\ncommand: {' '.join(command)}"
        )


def _worth_showing(line: str) -> bool:
    """pip is verbose; surface progress without drowning the terminal."""
    if not line.strip():
        return False
    return line.startswith(("Collecting", "Installing collected", "Successfully", "ERROR", "WARNING: Ignoring"))


def _warmup(spec: BackendSpec, on_output) -> dict:
    """Download the weights *and* run one real prediction.

    Downloading alone is not proof of anything: a dependency resolution that
    quietly broke the model would still leave the checkpoint on disk, and the
    backend would be marked ready until the first real drawing failed. So the
    warmup feeds it a generated probe image and only succeeds if a structure
    comes back.
    """
    import tempfile

    from .recognition._subprocess import run_worker

    with tempfile.TemporaryDirectory(prefix="m2i-warmup-") as workspace:
        probe = _render_probe(Path(workspace) / "probe.png")
        request = dict(spec.warmup_request)
        if probe is not None:
            request["image_path"] = str(probe)

        try:
            payload = run_worker(
                spec.name, spec.worker, request, timeout=WARMUP_TIMEOUT
            )
        except Exception as exc:  # noqa: BLE001 - reported, not swallowed
            raise SetupError(
                f"the packages installed but {spec.name} could not run its model:\n{exc}"
            ) from exc

    result = {k: v for k, v in payload.items() if k != "ok"}
    read = result.get("probe") or result.get("smiles")
    if probe is not None and not read:
        raise SetupError(
            f"{spec.name} loaded but returned nothing for the probe image; the "
            "installation is not usable"
        )
    on_output(f"    model ready (read the probe drawing as {read})")
    return result


def _render_probe(path: Path) -> Path | None:
    """A small depiction to prove the model runs, drawn with m2i's own RDKit."""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDepictor
        from rdkit.Chem.Draw import rdMolDraw2D

        mol = Chem.MolFromSmiles("CCO")
        rdDepictor.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DCairo(400, 300)
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
        drawer.FinishDrawing()
        path.write_bytes(drawer.GetDrawingText())
        return path
    except Exception:  # noqa: BLE001 - the warmup degrades to a load-only check
        return None


def describe(spec: BackendSpec) -> str:
    return (
        f"{spec.name}\n"
        f"  best at   {spec.strength}\n"
        f"  how       {spec.description}\n"
        f"  download  {spec.download_size}\n"
        f"  note      {spec.notes}"
    )


def worker_exists(spec: BackendSpec) -> bool:
    try:
        worker_path(spec.worker)
        return True
    except Exception:  # noqa: BLE001
        return False
