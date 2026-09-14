"""Browser interface: give it a structure, check the reading, pick a format.

Three steps, in the order a chemist actually works:

1. **Structure** -- a picture, a ChemDraw or molfile, or a SMILES.
2. **Check** -- the reading next to the original, editable. Nothing is
   generated until this has been looked at.
3. **Output** -- choose the program. Only the settings that program needs are
   shown, and what you download is the file for that format.

The 3D geometry is computed once per structure and reused. Switching from
Gaussian to ORCA rewrites a text file rather than re-embedding the molecule,
so both files carry identical coordinates.
"""

from __future__ import annotations

import hashlib
import io
import sys
import tempfile
from dataclasses import replace
from pathlib import Path

import streamlit as st

# Allow `streamlit run path/to/app.py` without the package being installed.
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from rdkit import Chem  # noqa: E402

from m2i import __version__, config, pipeline  # noqa: E402
from m2i.chem.conformers import ConformerOptions  # noqa: E402
from m2i.gui.hosting import HOSTED  # noqa: E402
from m2i.preprocess import ImageError, estimate_drawing_style, prepare_image  # noqa: E402
from m2i.recognition import BackendError, recognize  # noqa: E402
from m2i.recognition.manual import STRUCTURE_FILE_SUFFIXES, from_structure_file  # noqa: E402
from m2i.recognition.registry import (  # noqa: E402
    describe_backends,
    installed_vision_backends,
    resolve_backends,
)
from m2i.report import depict  # noqa: E402
from m2i.types import ERROR, INFO, WARNING, IssueLog, RecognitionResult  # noqa: E402

IMAGE_TYPES = ["png", "jpg", "jpeg", "bmp", "tif", "tiff", "webp"]
STRUCTURE_TYPES = [suffix.lstrip(".") for suffix in STRUCTURE_FILE_SUFFIXES]

PICTURE, STRUCTURE_FILE, SMILES, CRYSTAL = (
    "Picture", "ChemDraw / molfile", "SMILES", "Crystal (.cif)"
)

#: What the user chooses between, and what each choice produces.
FORMATS = {
    "gaussian": "Gaussian (.gjf)",
    "orca": "ORCA (.inp)",
    "xyz": "XYZ (.xyz)",
    "sdf": "SDF (.sdf)",
}
QM_PROGRAMS = ("gaussian", "orca")
DISPERSION = {"gaussian": config.GAUSSIAN_DISPERSION, "orca": config.ORCA_DISPERSION}
NO_DISPERSION = "none"
#: Geometries kept in memory per session; each holds every embedded conformer.
EMBEDDING_CACHE_SIZE = 6

st.set_page_config(page_title="m2i - molecule to input", page_icon="⚗️", layout="wide")


@st.cache_resource(show_spinner=False)
def _preload_recognition() -> bool:
    """Load the recognition models once per server, in the background.

    A model takes most of a minute to load and a couple of seconds to read a
    picture once loaded. Served to other people, that minute is spent as soon
    as someone opens the page -- while they pick their photo -- instead of on
    their first reading. Runs once for the whole process, not per visitor.
    """
    import threading

    from m2i.recognition import _subprocess

    def load() -> None:
        for backend in installed_vision_backends():
            try:
                _subprocess.persistent_worker(backend.name, backend.spec.worker).ask(
                    {"mode": "warmup"}, _subprocess.DEFAULT_TIMEOUT
                )
            except BackendError:
                pass  # the first real reading will say what is wrong

    threading.Thread(target=load, name="m2i-preload", daemon=True).start()
    return True


def main() -> None:
    if HOSTED:
        _preload_recognition()
    st.title("Molecule to calculation input")
    st.caption(
        f"m2i {__version__} - read a structure, check it, and generate the input "
        "for the program you use"
    )

    settings = sidebar()
    st.subheader("1. Structure")
    kind = st.segmented_control(
        "Start from", [PICTURE, STRUCTURE_FILE, SMILES, CRYSTAL], default=PICTURE,
        required=True, key="source_kind",
    )
    if kind == CRYSTAL:
        # A crystal is not drawn but rebuilt: its own page, sharing only the
        # safe upload handling and the per-format settings.
        from m2i.gui import crystal_page

        crystal_page.render(
            settings,
            save_upload=_save_upload,
            format_settings=format_settings,
            show_issues=_show_issues,
        )
        return

    source = structure_step(kind)
    metal_drawing = st.session_state.get("organometallic_upload")
    if kind == STRUCTURE_FILE and metal_drawing:
        # A drawn metal complex: the geometry has to be settled (which donor
        # is trans to which) before there is anything to check.
        from m2i.gui import organometallic_page

        organometallic_page.render(settings, format_settings=format_settings, **metal_drawing)
        return
    if source is None:
        if not _input_pending():
            st.info(
                "Start with a picture of the molecule (a photo, a scan or a screenshot), "
                "a ChemDraw or molfile, a SMILES, or a crystal structure."
            )
        return

    molecule, log = check_step(source, settings)
    if molecule is None:
        return
    output_step(molecule, log, source, settings)


# -- sidebar: everything that is not about the output format --------------


def sidebar() -> dict:
    st.sidebar.header("Electronic state")
    st.sidebar.caption("Derived from the structure unless overridden.")
    charge = None
    if st.sidebar.checkbox("Override the charge"):
        charge = int(st.sidebar.number_input("Charge", value=0, step=1, format="%d"))
    mult = None
    if st.sidebar.checkbox("Override the multiplicity"):
        mult = int(st.sidebar.number_input("Multiplicity", value=1, min_value=1, step=1))
    keep_all_fragments = st.sidebar.checkbox(
        "Keep every fragment",
        help="Salts and counter-ions are dropped by default, keeping the largest "
        "fragment. Tick this to calculate the whole assembly.",
    )

    st.sidebar.header("3D geometry")
    keep = st.sidebar.slider("Conformers to write (lowest in energy)", 1, 10, 1)
    n_confs = st.sidebar.number_input(
        "Conformers to embed (0 = automatic)", value=0, min_value=0, step=10
    )
    force_field = st.sidebar.selectbox(
        "Force field", ("mmff94s", "mmff94", "uff", "none")
    )
    seed = int(st.sidebar.number_input("Random seed", value=0xF00D, step=1))

    # A folder of this session, not one shared by everyone: served over the
    # web, two people generating the same molecule must not overwrite each
    # other. What you download is the same file either way.
    output_dir = str(_workspace() / "output")
    if not HOSTED:
        st.sidebar.header("Files")
        output_dir = st.sidebar.text_input("Also save to folder", output_dir, key="output_dir")

    with st.sidebar.expander("Recognition models"):
        for row in describe_backends():
            if row["name"] == "manual":
                continue
            mark = "installed" if row["available"] else "not installed"
            st.write(f"**{row['name']}** - {mark}")
            st.caption(f"Best at {row['strength']}.")
            if not row["available"] and not HOSTED:
                st.caption(f"`m2i setup {row['name']}`")

    return {
        "charge": charge,
        "mult": mult,
        "keep_all_fragments": keep_all_fragments,
        "conformers": ConformerOptions(
            n_confs=int(n_confs) or None, keep=keep, seed=seed, force_field=force_field
        ),
        "output_dir": Path(output_dir),
    }


# -- step 1: structure ----------------------------------------------------


def structure_step(kind: str) -> dict | None:
    """Show the input chosen by the user; return the structure it produced.

    Each kind of input keeps its own state, so flipping between them never
    lets a stale picture reading leak into the SMILES box or the other way
    round.
    """
    {PICTURE: picture_input, STRUCTURE_FILE: file_input, SMILES: smiles_input}[kind]()
    return _sources().get(kind)


def picture_input() -> None:
    upload = st.file_uploader(
        "Picture of the molecule - hand-drawn or from ChemDraw, photo or screenshot",
        type=IMAGE_TYPES,
        key="picture_upload",
    )
    if upload is None:
        _sources().pop(PICTURE, None)
        return

    path, digest = _save_upload(upload)
    prep_log = IssueLog()
    try:
        image = prepare_image(path, prep_log)
    except ImageError as exc:
        # A corrupt or mislabelled file must not take the page down with it.
        _sources().pop(PICTURE, None)
        st.error(f"That file could not be opened as a picture: {exc}")
        return
    style, _ = estimate_drawing_style(image)
    looks = "hand-drawn" if style == "hand_drawn" else "clean, software-drawn"
    suited = "decimer" if style == "hand_drawn" else "molscribe"

    installed = installed_vision_backends()
    if not installed:
        if HOSTED:
            st.warning(
                "This server does not read pictures: the recognition models are too large "
                "to host here. Type the structure you read, or upload the ChemDraw file, "
                "and the rest is automatic."
            )
        else:
            st.warning(
                f"No recognition model is installed. This picture looks {looks}; the "
                f"model suited to it is installed with `m2i setup {suited}`. Until "
                "then, type the structure you read and the rest is automatic."
            )
        typed = st.text_input("SMILES read from the picture", key="picture_smiles").strip()
        if typed:
            _set_source(PICTURE, f"picture:{digest}:{typed}", _manual(typed), image, prep_log)
        else:
            _sources().pop(PICTURE, None)
        return

    compare = False
    if len(installed) > 1:
        compare = st.checkbox(
            "Run both models and compare their answers",
            help="Slower, but when both read the same molecule (stereochemistry "
            "included) that is the strongest confidence signal available.",
        )
    else:
        st.caption(f"The picture looks {looks}; {installed[0].name} will read it.")

    key = f"picture:{digest}:{'all' if compare else 'auto'}"
    current = _sources().get(PICTURE)
    if current and current["key"] == key:
        return
    _sources().pop(PICTURE, None)  # a new picture: the old reading is stale

    if not st.button("Read the structure", type="primary"):
        return
    log = IssueLog(list(prep_log))
    backends = [
        b
        for b in resolve_backends(["all"] if compare else None, log, style=style)
        if b.name != "manual"
    ]
    with st.spinner("Reading the drawing. The first run loads the model and takes longer."):
        try:
            recognition = recognize(path, backends, log, hand_drawn=(style == "hand_drawn"))
        except BackendError as exc:
            st.error(str(exc))
            return
    _set_source(PICTURE, key, recognition, image, log)
    st.rerun()  # redraw without the button: this picture has been read


def file_input() -> None:
    upload = st.file_uploader(
        "ChemDraw (.cdxml, .cdx) or molfile (.mol, .sdf)",
        type=STRUCTURE_TYPES,
        key="structure_upload",
    )
    if upload is None:
        _sources().pop(STRUCTURE_FILE, None)
        st.session_state.pop("organometallic_upload", None)
        return
    st.caption(
        "Read straight from the file: bonds and wedges are exactly as drawn, with "
        "no recognition model involved."
    )

    path, digest = _save_upload(upload)
    if _draws_a_metal(path):
        _sources().pop(STRUCTURE_FILE, None)
        st.session_state["organometallic_upload"] = {
            "path": path, "digest": digest, "shown_name": Path(upload.name).name,
        }
        return
    st.session_state.pop("organometallic_upload", None)
    key = f"file:{digest}"
    current = _sources().get(STRUCTURE_FILE)
    if current and current["key"] == key:
        return
    try:
        recognition = from_structure_file(path).recognize(None)
    except Exception as exc:  # noqa: BLE001 - RDKit raises assorted types on bad files
        _sources().pop(STRUCTURE_FILE, None)
        st.error(f"That file could not be read: {exc}")
        return
    _set_source(STRUCTURE_FILE, key, recognition, None, IssueLog())


def _draws_a_metal(path: Path) -> bool:
    """Whether the drawing holds a metal, which needs the other route entirely.

    Cached: the answer decides which page is drawn, so it is asked on every
    rerun, and reading the file again each time would be wasteful.
    """
    from m2i.organometallic import has_metal, read_raw

    cache = st.session_state.setdefault("metal_drawings", {})
    key = str(path)
    if key not in cache:
        try:
            cache[key] = any(has_metal(m) for m in read_raw(path))
        except Exception:  # noqa: BLE001 - the organic reader reports a bad file properly
            cache[key] = False
    return cache[key]


def smiles_input() -> None:
    typed = st.text_input(
        "SMILES", key="smiles_text", placeholder="C[C@H](N)C(=O)O"
    ).strip()
    if not typed:
        _sources().pop(SMILES, None)
        return
    _set_source(SMILES, f"smiles:{typed}", _manual(typed), None, IssueLog())


# -- step 2: check ---------------------------------------------------------


def check_step(source: dict, settings: dict):
    """The verification gate. Returns (molecule, log), or (None, log)."""
    st.subheader("2. Check the structure")

    corrections = st.session_state.setdefault("corrections", {})
    corrected = corrections.get(source["key"])
    recognition = corrected or source["recognition"]

    log = IssueLog() if corrected else IssueLog(source["issues"])
    if corrected:
        log.info(
            "recognition.corrected",
            f"Structure corrected by hand. The original reading "
            f"({source['recognition'].backend}) was {source['display_smiles']}.",
        )

    molecule, error = None, None
    try:
        molecule = pipeline.prepare_molecule(
            recognition,
            log,
            charge=settings["charge"],
            multiplicity=settings["mult"],
            keep_all_fragments=settings["keep_all_fragments"],
        )
    except Exception as exc:  # noqa: BLE001 - shown to the user, who can fix it
        error = str(exc)

    left, right = st.columns(2)
    with left:
        if source["image"] is not None:
            st.image(source["image"], caption="Original", width="stretch")
        else:
            st.caption("No picture for this structure.")
    with right:
        if molecule is not None:
            png = depict.draw_molecule(molecule.mol, stereo=molecule.stereo)
            st.image(
                io.BytesIO(png),
                caption="Understood by m2i - atoms numbered as in the generated input",
                width="stretch",
            )
        else:
            st.error(f"This structure cannot be used: {error}")

    # Even when the reading is unusable, the way out is to correct it here.
    _correction_field(source, corrections)

    if molecule is None:
        return None, log

    columns = st.columns(4)
    columns[0].metric("Formula", molecule.formula)
    columns[1].metric("Charge", molecule.charge)
    columns[2].metric("Multiplicity", molecule.multiplicity)
    columns[3].metric("Stereocentres", len(molecule.stereo.centers))
    st.caption(f"InChIKey {molecule.inchikey or '(unavailable)'}")
    st.write(f"**Stereochemistry:** {molecule.stereo.describe()}")

    _show_issues(log)
    return molecule, log


def _correction_field(source: dict, corrections: dict) -> None:
    applied = corrections.get(source["key"])
    applied_text = applied.smiles if applied else source["display_smiles"]

    text = st.text_input(
        "SMILES - edit it to correct the reading",
        value=source["display_smiles"],
        key=f"correction:{source['key']}",
    ).strip()
    # Compare with what was last applied, not with the canonical SMILES: a
    # correction typed in non-canonical form would otherwise differ from its
    # own canonicalisation on every run and loop forever.
    if text and text != applied_text:
        corrections[source["key"]] = _manual(text)
        st.rerun()

    if applied and st.button("Back to the original reading"):
        corrections.pop(source["key"], None)
        st.session_state.pop(f"correction:{source['key']}", None)
        st.rerun()


# -- step 3: output --------------------------------------------------------


def output_step(molecule, log: IssueLog, source: dict, settings: dict) -> None:
    st.subheader("3. Output")
    if log.has_errors():
        st.error("Resolve the errors above before generating an input.")
        return

    fmt = st.segmented_control(
        "Format", list(FORMATS), format_func=FORMATS.get, default="gaussian",
        required=True, key="format",
    )
    profile = format_settings(fmt)
    if profile is None:
        return

    conformers = settings["conformers"]
    signature = (
        source["key"],
        molecule.smiles,
        molecule.charge,
        molecule.multiplicity,
        repr(profile),
        repr(conformers),
        str(settings["output_dir"]),
    )

    if st.button(f"Generate {FORMATS[fmt]}", type="primary"):
        _generate(molecule, log, source, settings, profile, signature)

    stored = st.session_state.get("result")
    if stored and stored["signature"] == signature:
        _show_result(stored, fmt, log)


def format_settings(fmt: str) -> config.JobProfile | None:
    """The settings this format needs, and nothing else."""
    if fmt not in QM_PROGRAMS:
        if fmt == "xyz":
            st.caption(
                "Geometry only: element symbols and Cartesian coordinates, readable "
                "by practically every program."
            )
        else:
            st.caption(
                "Geometry only, but it keeps bonds, charges and stereochemistry: "
                "the format to hand to another cheminformatics tool."
            )
        return config.JobProfile(
            name=f"{fmt}_geometry", program=fmt, method="", basis="", jobs=()
        )

    recipes = _profiles_for(fmt)
    if recipes:
        names = list(recipes)
        name = st.selectbox(
            "Recipe",
            names,
            format_func=lambda n: f"{n} - {recipes[n].description}"
            if recipes[n].description
            else n,
            key=f"recipe:{fmt}",
        )
        base = recipes[name]
    else:
        name, base = "default", config.JobProfile(program=fmt)

    # Keys include the recipe so choosing another recipe loads its defaults.
    scope = f"{fmt}:{name}"
    first, second, third = st.columns(3)
    method = first.text_input("Method", base.method, key=f"method:{scope}")
    basis = second.text_input("Basis set", base.basis, key=f"basis:{scope}")
    choices = [NO_DISPERSION, *DISPERSION[fmt]]
    current = (base.dispersion or NO_DISPERSION).lower()
    dispersion = third.selectbox(
        "Dispersion",
        choices,
        index=choices.index(current) if current in choices else 0,
        key=f"dispersion:{scope}",
    )

    first, second, third = st.columns(3)
    solvent = first.text_input(
        "Solvent (blank = gas phase)",
        base.solvent.name if base.solvent else "",
        key=f"solvent:{scope}",
    ).strip()
    nproc = second.number_input(
        "Cores", value=base.resources.nproc, min_value=1, step=1, key=f"nproc:{scope}"
    )
    mem = third.text_input("Memory", base.resources.mem, key=f"mem:{scope}")

    try:
        profile = replace(
            base,
            method=method.strip(),
            basis=basis.strip(),
            dispersion=None if dispersion == NO_DISPERSION else dispersion,
            solvent=config.Solvent(
                model=base.solvent.model if base.solvent else "", name=solvent
            )
            if solvent
            else None,
            resources=config.Resources(
                mem=mem.strip() or base.resources.mem,
                nproc=int(nproc),
                maxcore_mb=base.resources.maxcore_mb,
            ),
        )
        profile.validate()
    except config.ProfileError as exc:
        st.error(str(exc))
        return None

    jobs = " ".join(profile.jobs) or "single point"
    st.caption(f"Jobs from the recipe: {jobs}.")
    return profile


def _generate(molecule, log, source, settings, profile, signature) -> None:
    conformers = settings["conformers"]
    cache = st.session_state.setdefault("embeddings", {})
    key = (molecule.smiles, repr(conformers))

    if key not in cache:
        embed_log = IssueLog()
        with st.spinner(
            "Embedding and minimising conformers. This is done once per structure; "
            "other formats reuse the same geometry."
        ):
            try:
                embedding = pipeline.embed(molecule, conformers, embed_log)
            except Exception as exc:  # noqa: BLE001
                st.error(f"Could not build a 3D structure: {exc}")
                return
        while len(cache) >= EMBEDDING_CACHE_SIZE:
            cache.pop(next(iter(cache)))
        cache[key] = (embedding, list(embed_log))
    embedding, embed_issues = cache[key]

    run_log = IssueLog(list(log) + embed_issues)
    options = pipeline.PipelineOptions(
        profile=profile,
        output_dir=settings["output_dir"],
        conformers=conformers,
        source_image=source["image"],
    )
    try:
        result = pipeline.write_inputs(molecule, embedding, options, run_log)
    except Exception as exc:  # noqa: BLE001
        st.error(f"Could not write the input: {exc}")
        return
    st.session_state["result"] = {
        "signature": signature,
        "result": result,
        "issues": run_log,
    }


def _show_result(stored: dict, fmt: str, check_log: IssueLog) -> None:
    result = stored["result"]
    written = [Path(p) for p in result.written_files]
    inputs = [p for p in written if not p.name.endswith(("_check.png", ".m2i.json"))]
    extras = [p for p in written if p not in inputs]

    where = "" if HOSTED or not inputs else f" A copy is in {inputs[0].parent}."
    st.success(
        f"{len(inputs)} {FORMATS[fmt]} file(s) ready.{where}" if inputs else "Nothing was written."
    )

    for conformer, path in zip(result.conformers, inputs):
        data = path.read_bytes()
        label = f"conformer {conformer.index + 1}"
        if conformer.relative_energy is not None:
            label += f", +{conformer.relative_energy:.2f} kcal/mol ({conformer.force_field})"
        left, right = st.columns([1, 3])
        left.download_button(
            f"Download {path.name}",
            data,
            file_name=path.name,
            key=f"download:{path.name}",
            type="primary",
        )
        with right.expander(f"{path.name} - {label}"):
            st.code(data.decode("utf-8", errors="replace"), language="text")

    # Only what the 3D step or the writer added; the check step showed the rest.
    already = {(i.code, i.message) for i in check_log}
    for issue in stored["issues"]:
        if issue.level == WARNING and (issue.code, issue.message) not in already:
            st.warning(issue.message)

    if extras:
        with st.expander("Also written: comparison image and provenance record"):
            for path in extras:
                st.download_button(
                    f"Download {path.name}",
                    path.read_bytes(),
                    file_name=path.name,
                    key=f"download:{path.name}",
                )


# -- helpers ---------------------------------------------------------------


def _input_pending() -> bool:
    """Something has been given but not turned into a structure yet -- a
    picture waiting for its button, or a file that failed. The prompt to
    start would be the wrong thing to show then."""
    state = st.session_state
    return any(state.get(key) is not None for key in ("picture_upload", "structure_upload"))


def _sources() -> dict:
    return st.session_state.setdefault("sources", {})


def _set_source(kind, key, recognition, image, log) -> None:
    if _sources().get(kind, {}).get("key") == key:
        return
    _sources()[kind] = {
        "key": key,
        "recognition": recognition,
        "image": image,
        "issues": list(log),
        "display_smiles": _display_smiles(recognition),
    }


def _display_smiles(recognition: RecognitionResult) -> str:
    """What to put in the correction box: the reading itself, as a SMILES."""
    if recognition.molblock:
        mol = Chem.MolFromMolBlock(recognition.molblock)
        if mol is not None:
            return Chem.MolToSmiles(mol)
    return recognition.smiles or ""


def _manual(smiles: str) -> RecognitionResult:
    return RecognitionResult(smiles=smiles, backend="manual", confidence=1.0)


def _workspace() -> Path:
    """A private temporary folder for this session's uploads."""
    if "workspace" not in st.session_state:
        st.session_state["workspace"] = tempfile.mkdtemp(prefix="m2i-session-")
    return Path(st.session_state["workspace"])


def _save_upload(upload) -> tuple[Path, str]:
    """Store an upload under a name m2i chooses.

    The browser-supplied file name never becomes part of a path: only its
    extension is kept, and the uploader has already restricted that to the
    accepted types.
    """
    data = upload.getvalue()
    digest = hashlib.sha256(data).hexdigest()[:16]
    path = _workspace() / f"{digest}{Path(upload.name).suffix.lower()}"
    if not path.exists():
        path.write_bytes(data)
    return path, digest


def _profiles_for(program: str) -> dict[str, config.JobProfile]:
    found = {}
    for name in sorted(config.available_profiles()):
        try:
            profile = config.load_profile(name)
        except config.ProfileError:
            continue
        if profile.program == program:
            found[name] = profile
    return found


def _show_issues(log: IssueLog) -> None:
    for issue in log:
        if issue.level == ERROR:
            st.error(issue.message)
        elif issue.level == WARNING:
            st.warning(issue.message)
    notes = [issue for issue in log if issue.level == INFO]
    if notes:
        with st.expander("Details"):
            for issue in notes:
                st.caption(issue.message)


main()
