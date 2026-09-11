"""Browser interface: drop a drawing, check what m2i read, generate the input.

The layout follows the one rule that matters: nothing is written until the
structure has been shown next to the original and the chemist has looked at it.
The SMILES field is editable precisely because the recognition will sometimes
be wrong, and correcting it by hand must be the fastest path forward.
"""

from __future__ import annotations

import io
import sys
import tempfile
from pathlib import Path

import streamlit as st

# Allow `streamlit run path/to/app.py` without the package being installed.
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from m2i import __version__, config, pipeline  # noqa: E402
from m2i.chem.conformers import ConformerOptions  # noqa: E402
from m2i.preprocess import prepare_image  # noqa: E402
from m2i.recognition import recognize  # noqa: E402
from m2i.recognition.manual import ManualBackend  # noqa: E402
from m2i.recognition.registry import describe_backends, resolve_backends  # noqa: E402
from m2i.report import depict  # noqa: E402
from m2i.types import ERROR, INFO, WARNING, IssueLog, RecognitionResult  # noqa: E402
from m2i.writers import known_programs  # noqa: E402

st.set_page_config(page_title="m2i - molecule to input", page_icon="*", layout="wide")


def main() -> None:
    st.title("Molecule to calculation input")
    st.caption(
        f"m2i {__version__} - read a drawn structure, check it, generate the "
        "quantum-chemistry input"
    )

    settings = sidebar()
    recognition, source_image = structure_input()

    if recognition is None:
        st.info(
            "Give m2i a structure: upload a picture and type what you read, drop a "
            "ChemDraw .mol export, or paste a SMILES."
        )
        return

    molecule, log = parse(recognition)
    if molecule is None:
        return

    verification_panel(molecule, log, source_image)
    generation_panel(molecule, log, settings, source_image)


# -- panels --------------------------------------------------------------


def sidebar() -> dict:
    st.sidebar.header("Calculation")
    profiles = sorted(config.available_profiles())
    default = profiles.index("gaussian_opt_freq") if "gaussian_opt_freq" in profiles else 0
    profile_name = st.sidebar.selectbox("Profile", profiles, index=default)

    profile = config.load_profile(profile_name)
    if profile.description:
        st.sidebar.caption(profile.description)

    program = st.sidebar.selectbox(
        "Program",
        known_programs(),
        index=known_programs().index(profile.program)
        if profile.program in known_programs()
        else 0,
    )
    method = st.sidebar.text_input("Method", profile.method)
    basis = st.sidebar.text_input("Basis set", profile.basis)

    st.sidebar.header("Electronic state")
    charge = st.sidebar.number_input(
        "Charge (blank = derived)", value=0, step=1, format="%d"
    )
    use_charge = st.sidebar.checkbox("Force this charge", value=False)
    mult = st.sidebar.number_input("Multiplicity", value=1, min_value=1, step=1)
    use_mult = st.sidebar.checkbox("Force this multiplicity", value=False)

    st.sidebar.header("Conformers")
    keep = st.sidebar.slider("Files to write (lowest in energy)", 1, 10, 1)
    n_confs = st.sidebar.number_input(
        "Conformers to embed (0 = automatic)", value=0, min_value=0, step=10
    )
    force_field = st.sidebar.selectbox(
        "Force field", ("mmff94s", "mmff94", "uff", "none")
    )
    seed = st.sidebar.number_input("Random seed", value=0xF00D, step=1)

    st.sidebar.header("Output")
    output_dir = st.sidebar.text_input(
        "Folder", str(Path(tempfile.gettempdir()) / "m2i_output")
    )

    with st.sidebar.expander("Recognition backends"):
        for row in describe_backends():
            mark = "ready" if row["available"] else "not installed"
            st.write(f"**{row['name']}** - {mark}")
            st.caption(row["description"])

    return {
        "profile": profile.merged_with(program=program, method=method, basis=basis),
        "charge": int(charge) if use_charge else None,
        "mult": int(mult) if use_mult else None,
        "keep": keep,
        "n_confs": int(n_confs) or None,
        "force_field": force_field,
        "seed": int(seed),
        "output_dir": Path(output_dir),
    }


def structure_input() -> tuple[RecognitionResult | None, object | None]:
    """Returns (recognition, preprocessed source image)."""
    image_tab, molfile_tab, smiles_tab = st.tabs(
        ["From a drawing", "From a .mol / .sdf", "From a SMILES"]
    )

    with image_tab:
        upload = st.file_uploader(
            "Picture of the molecule (hand-drawn or ChemDraw)",
            type=["png", "jpg", "jpeg", "bmp", "tif", "tiff", "webp"],
        )
        source_image = None
        if upload is not None:
            temp = Path(tempfile.gettempdir()) / f"m2i_upload_{upload.name}"
            temp.write_bytes(upload.getvalue())
            source_image = prepare_image(temp, IssueLog())
            st.session_state["source_image"] = source_image

            backends = resolve_backends(None, IssueLog())
            vision = [b for b in backends if b.name != "manual"]
            if vision:
                if st.button("Recognise structure", type="primary"):
                    log = IssueLog()
                    st.session_state["recognition"] = recognize(temp, vision, log)
                    st.session_state["smiles_field"] = st.session_state[
                        "recognition"
                    ].smiles
            else:
                st.warning(
                    "No vision backend is installed yet (phase 2). Read the structure "
                    "yourself and type it below - everything after this point is "
                    "automatic."
                )
                typed = st.text_input("SMILES you read from the drawing", key="image_smiles")
                if typed.strip():
                    st.session_state["recognition"] = RecognitionResult(
                        smiles=typed.strip(), backend="manual", confidence=1.0
                    )
                    st.session_state["smiles_field"] = typed.strip()

    with molfile_tab:
        molupload = st.file_uploader("Molfile", type=["mol", "sdf"], key="molupload")
        if molupload is not None:
            text = molupload.getvalue().decode("utf-8", errors="replace")
            st.session_state["recognition"] = ManualBackend(
                molblock=text.split("$$$$")[0]
            ).recognize(None)
            st.caption(
                "Molblock loaded: stereochemistry comes from the wedge bonds and the "
                "2D layout."
            )

    with smiles_tab:
        typed = st.text_input(
            "SMILES", value=st.session_state.get("smiles_field", ""), key="smiles_direct"
        )
        if typed.strip():
            st.session_state["recognition"] = RecognitionResult(
                smiles=typed.strip(), backend="manual", confidence=1.0
            )

    return st.session_state.get("recognition"), st.session_state.get("source_image")


def parse(recognition: RecognitionResult):
    log = IssueLog()
    try:
        molecule = pipeline.prepare_molecule(recognition, log)
    except Exception as exc:
        st.error(f"Could not read that structure: {exc}")
        return None, log
    return molecule, log


def verification_panel(molecule, log: IssueLog, source_image) -> None:
    st.subheader("1. Check the structure")
    left, right = st.columns(2)

    with left:
        if source_image is not None:
            st.image(source_image, caption="Original drawing", use_container_width=True)
        else:
            st.caption("No source image for this structure.")

    with right:
        png = depict.draw_molecule(molecule.mol, stereo=molecule.stereo)
        st.image(
            io.BytesIO(png),
            caption="Understood by m2i - atoms numbered as in the generated input",
            use_container_width=True,
        )

    corrected = st.text_input(
        "SMILES (edit to correct the reading)", value=molecule.smiles, key="correction"
    )
    if corrected.strip() and corrected.strip() != molecule.smiles:
        st.session_state["recognition"] = RecognitionResult(
            smiles=corrected.strip(), backend="manual", confidence=1.0
        )
        st.rerun()

    columns = st.columns(4)
    columns[0].metric("Formula", molecule.formula)
    columns[1].metric("Charge", molecule.charge)
    columns[2].metric("Multiplicity", molecule.multiplicity)
    columns[3].metric("Stereocentres", len(molecule.stereo.centers))

    st.caption(f"InChIKey {molecule.inchikey or '(unavailable)'}")
    st.write(f"**Stereochemistry:** {molecule.stereo.describe()}")

    for issue in log:
        if issue.level == ERROR:
            st.error(issue.message)
        elif issue.level == WARNING:
            st.warning(issue.message)

    with st.expander("Details"):
        for issue in log:
            if issue.level == INFO:
                st.caption(issue.message)


def generation_panel(molecule, log: IssueLog, settings: dict, source_image) -> None:
    st.subheader("2. Generate the input")

    if log.has_errors():
        st.error("Resolve the errors above before generating an input file.")
        return

    profile = settings["profile"]
    st.caption(
        f"{profile.program} - {profile.method}/{profile.basis} "
        f"{' '.join(profile.jobs)} - {profile.resources.nproc} cores, "
        f"{profile.resources.mem}"
    )

    if not st.button("Generate", type="primary"):
        return

    options = pipeline.PipelineOptions(
        profile=profile,
        output_dir=settings["output_dir"],
        charge=settings["charge"],
        multiplicity=settings["mult"],
        conformers=ConformerOptions(
            n_confs=settings["n_confs"],
            keep=settings["keep"],
            seed=settings["seed"],
            force_field=settings["force_field"],
        ),
        source_image=source_image,
    )

    run_log = IssueLog(list(log))
    with st.spinner("Embedding conformers and writing files..."):
        try:
            result = pipeline.generate_inputs(molecule, options, run_log)
        except Exception as exc:
            st.error(f"Generation failed: {exc}")
            return

    st.success(f"{len(result.written_files)} file(s) written to {settings['output_dir']}")

    for conformer in result.conformers:
        label = f"conformer {conformer.index + 1}"
        if conformer.relative_energy is not None:
            label += f" (+{conformer.relative_energy:.2f} kcal/mol {conformer.force_field})"
        st.caption(label)

    for path_str in result.written_files:
        path = Path(path_str)
        if path.suffix.lower() == ".png":
            continue
        data = path.read_bytes()
        st.download_button(
            f"Download {path.name}", data, file_name=path.name, key=f"dl_{path.name}"
        )
        if path.suffix.lower() in (".gjf", ".inp", ".xyz"):
            with st.expander(f"Preview {path.name}"):
                st.code(data.decode("utf-8", errors="replace"), language="text")

    for issue in run_log.warnings:
        st.warning(issue.message)


main()
