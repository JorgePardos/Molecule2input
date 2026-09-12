"""The crystal (.cif) route in the browser.

Different from the other inputs in one fundamental way: nothing is drawn or
generated. The molecule is rebuilt from the crystal, so what has to be checked
is not "did m2i read my drawing" but "is this the species I meant, complete,
with the right coordination" -- hence a 3D view and the coordination sphere
instead of a 2D depiction, and questions instead of derived values for
everything a crystal cannot tell: oxidation state, charge, spin.
"""

from __future__ import annotations

import json
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

from m2i.chem import metals
from m2i.crystal import CrystalError, read_cif
from m2i.crystal import jobs as crystal_jobs
from m2i.crystal.coordination import analyse
from m2i.types import WARNING, IssueLog

#: A CIF has no bond orders, so an SDF would be a list of atoms pretending to
#: be a molecule file. The formats offered are the ones that need none.
FORMATS = {"gaussian": "Gaussian (.gjf)", "orca": "ORCA (.inp)", "xyz": "XYZ (.xyz)"}
#: Crystals kept per session; each holds every species of the structure.
CRYSTAL_CACHE_SIZE = 4
#: A local Streamlit component with 3Dmol.js bundled beside it, so the view
#: works without a network connection.
_VIEWER = components.declare_component(
    "m2i_viewer", path=str(Path(__file__).parent / "viewer")
)


def render(settings: dict, *, save_upload, format_settings, show_issues) -> None:
    upload = st.file_uploader(
        "Crystal structure (.cif) - from the CSD, the COD or your own refinement",
        type=["cif"],
        key="cif_upload",
    )
    if upload is None:
        st.info(
            "The geometry is taken from the crystal as measured, so stereochemistry "
            "and ligand shapes are kept exactly. Oxidation state, charge and spin are "
            "not in a CIF: you will be asked for them."
        )
        return

    normalise = st.checkbox(
        "Extend C-H, N-H and O-H to standard lengths",
        value=True,
        help="X-ray places hydrogens about 0.1 Å too close to their atom. "
        "Hydrides and other M-H are never moved.",
    )
    path, digest = save_upload(upload)
    reading, issues = _read(path, digest, normalise, Path(upload.name).name)
    if reading is None:
        return

    _crystal_summary(reading)
    index = _choose_species(reading, digest)
    original = reading.species[index]

    st.subheader("2. Check the structure")
    show_issues(IssueLog(issues))
    species, hydrogens_added, allow_missing = original, None, False
    if reading.missing_hydrogens:
        from m2i.crystal.hydrogens import lacking

        if lacking(reading, index):
            species, hydrogens_added, allow_missing = _hydrogens(reading, index, digest)
        else:
            st.caption(
                f"The hydrogens missing from the crystal belong to other species; "
                f"{original.formula} is complete."
            )
            allow_missing = True

    left, right = st.columns(2)
    with left:
        viewer(species)
    with right:
        coordination(original)
    priorities, mirror = _chirality(reading, original, digest)

    st.subheader("3. Electronic state")
    state = _electronic_state(reading, index, digest, species)

    st.subheader("4. Output")
    fmt = st.segmented_control(
        "Format", list(FORMATS), format_func=FORMATS.get, default="gaussian",
        required=True, key="cif_format",
    )
    profile = format_settings(fmt)
    if profile is None or state is None:
        return
    if reading.missing_hydrogens and hydrogens_added is None and not allow_missing:
        st.error("Hydrogens are missing from the CIF (see above); nothing will be written.")
        return

    charge, multiplicity, oxidation = state
    signature = (digest, normalise, index, charge, multiplicity, json.dumps(oxidation),
                 repr(profile), str(settings["output_dir"]), species.formula,
                 tuple(hydrogens_added or ()), mirror, json.dumps(priorities))
    if st.button(f"Generate {FORMATS[fmt]}", type="primary"):
        log = IssueLog(issues)
        try:
            written = crystal_jobs.write(
                reading, index, charge=charge, multiplicity=multiplicity,
                profile=profile, output_dir=settings["output_dir"], log=log,
                oxidation=oxidation, allow_missing_hydrogens=allow_missing,
                species=species if hydrogens_added is not None else None,
                hydrogens_added=hydrogens_added, mirror=mirror, priorities=priorities,
            )
        except CrystalError as exc:
            st.error(str(exc))
            return
        st.session_state["cif_result"] = {
            "signature": signature, "files": written, "issues": list(log),
        }

    stored = st.session_state.get("cif_result")
    if stored and stored["signature"] == signature:
        show_result(stored, fmt, issues)


# -- pieces ------------------------------------------------------------------


def _read(path: Path, digest: str, normalise: bool, shown_name: str):
    cache = st.session_state.setdefault("crystals", {})
    key = (digest, normalise)
    if key not in cache:
        log = IssueLog()
        try:
            reading = read_cif(path, log, normalise_hydrogens=normalise)
        except CrystalError as exc:
            st.error(str(exc))
            return None, []
        # The upload was stored under a hash; show and name files by the original.
        reading.path = Path(shown_name)
        while len(cache) >= CRYSTAL_CACHE_SIZE:  # a session must not grow without end
            cache.pop(next(iter(cache)))
        cache[key] = (reading, list(log))
    return cache[key]


def _crystal_summary(reading) -> None:
    source = f"COD {reading.cod_id}" if reading.cod_id else reading.path.name
    z = f" · Z = {reading.z:g}" if reading.z else ""
    a, b, c, alpha, beta, gamma = reading.cell
    st.caption(
        f"{source} · {reading.spacegroup}{z} · a {a:g}, b {b:g}, c {c:g} Å, "
        f"α {alpha:g}°, β {beta:g}°, γ {gamma:g}°"
    )
    rows = [
        {
            "Species": s.formula,
            "Copies in the cell": s.copies,
            "Atoms": s.n_atoms,
            "Metal": ", ".join(sorted(set(s.metals))) or "-",
        }
        for s in reading.species
    ]
    st.dataframe(rows, hide_index=True, width="stretch")


def _choose_species(reading, digest: str) -> int:
    if len(reading.species) == 1:
        return 0
    return st.radio(
        "Species to calculate",
        list(range(len(reading.species))),
        format_func=lambda i: f"{reading.species[i].formula} (x{reading.species[i].copies})",
        horizontal=True,
        key=f"species:{digest}",
    )


def _hydrogens(reading, index: int, digest: str):
    """Propose the missing hydrogens, let the user adjust them, place them.

    Returns (species to write, what was added or None, allow writing without).
    """
    from m2i.crystal import hydrogens

    proposal = hydrogens.plan(reading, index)
    st.markdown("**Missing hydrogens**")
    rows = [
        {
            "Atom": s.label,
            "H to add": s.count,
            "Why": s.reason,
            "Sure": "yes" if s.certain else "assumed",
        }
        for s in proposal.sites
    ]
    edited = st.data_editor(
        rows,
        key=f"hydrogens:{digest}:{index}",
        hide_index=True,
        width="stretch",
        disabled=["Atom", "Why", "Sure"],
        column_config={
            "H to add": st.column_config.NumberColumn(min_value=0, max_value=4, step=1),
        },
    )
    counts = {row["Atom"]: int(row["H to add"] or 0) for row in edited}
    changed = {label: n for label, n in counts.items()
               if n != next(s.count for s in proposal.sites if s.label == label)}
    final = proposal.with_counts(changed)

    if final.target is None:
        st.warning(
            f"{final.total} hydrogens proposed. The formula does not say how many this "
            "species needs, so check the list."
        )
    elif final.matches:
        st.success(f"{final.total} hydrogens, as the formula requires.")
    else:
        st.error(
            f"{final.total} hydrogens proposed, but the formula needs {final.target}. "
            "Adjust the rows marked 'assumed'."
        )

    add = st.checkbox(
        "Add these hydrogens", value=final.matches, key=f"add_h:{digest}:{index}",
    )
    if add:
        completed, added = hydrogens.apply(reading.species[index], final)
        return completed, added, False
    allow = st.checkbox(
        "Write the input without them - I know the hydrogens are missing",
        key=f"allow_missing:{digest}",
    )
    return reading.species[index], None, allow


def _chirality(reading, species, digest: str):
    """Notes on chirality, the questions it raises, and the mirror-image option."""
    centres = analyse(species)
    priorities = {}
    for centre in centres:
        for ring in centre.planar:
            if ring.certain:
                continue
            st.caption(f"Ring {ring.ring_label}: {ring.reason}.")
            priorities[ring.ring_label] = st.radio(
                f"Higher CIP priority in ring {ring.ring_label}",
                [ring.first, ring.second],
                horizontal=True,
                key=f"priority:{digest}:{ring.ring_label}",
            )
    notes, chiral = crystal_jobs.chirality_notes(reading, centres, priorities)
    mirror = False
    if chiral:
        st.markdown("**Chirality**")
        for note in notes:
            st.caption(note)
        mirror = st.checkbox(
            "Write the mirror image (the other enantiomer)", key=f"mirror:{digest}"
        )
    return priorities, mirror


def coordination(species) -> None:
    centres = analyse(species)
    if not centres:
        st.caption("No metal in this species.")
        return
    for centre in centres:
        st.markdown(f"**{centre.label}** - {centre.geometry}"
                    + (f" ({centre.detail})" if centre.detail else ""))
        # A table, not preformatted text: a narrow column would cut the
        # distances off, and they are the part that matters.
        st.dataframe(
            [
                {
                    "Donor": (f"η{d.hapticity}-" if d.hapticity > 1 else "") + d.label,
                    "Distance (Å)": round(d.distance, 3),
                    "Type": d.kind,
                }
                for d in centre.donors
            ],
            hide_index=True,
            width="stretch",
        )
        for isomer in centre.isomers:
            st.markdown(f"**{isomer}**")
        for a, b, angle in centre.trans_pairs:
            st.caption(f"trans: {a} / {b} ({angle:.1f}°)")


def _electronic_state(reading, index: int, digest: str, species):
    """Ask what the crystal cannot say. Returns (charge, mult, oxidation) or None.

    ``species`` is what will be written -- with any hydrogens put back, which
    matters: they change the electron count the multiplicity must fit.
    """
    asked = crystal_jobs.questions(reading, index)
    # The formula is part of the key: adding hydrogens changes which spin
    # states the electron count allows, and a stale choice must not survive.
    scope = f"{digest}:{index}:{species.formula}"

    oxidation: dict[str, int] = {}
    if asked.needs_answers:
        st.write(
            f"**{species.formula}** contains {', '.join(asked.metals)}. A crystal does "
            "not record oxidation states, charge or spin, so these are yours to set."
        )
        columns = st.columns(max(1, len(asked.metals)))
        for column, element in zip(columns, asked.metals):
            value = column.number_input(
                f"Oxidation state of {element}", value=None, step=1,
                placeholder="e.g. 2", key=f"ox:{scope}:{element}",
            )
            if value is not None:
                oxidation[element] = int(value)
        for line in crystal_jobs.spin_hint(species, oxidation):
            st.caption(line)

    left, right = st.columns(2)
    charge = left.number_input(
        "Charge", value=asked.suggested_charge, step=1,
        placeholder="required", key=f"charge:{scope}",
    )
    left.caption(f"Suggestion: {asked.charge_reason}.")
    if charge is None:
        right.number_input("Multiplicity", value=None, disabled=True, key=f"mult_off:{scope}")
        st.warning("Set the charge to continue.")
        return None
    charge = int(charge)

    mult = multiplicity(right, species, charge, oxidation, asked.needs_answers, scope)
    if mult is None:
        st.warning("Choose the spin state to continue.")
        return None

    problem = metals.parity_problem(species.elements, charge, mult)
    if problem:
        st.error(f"Impossible electronic state: {problem}.")
        return None
    return charge, mult, oxidation


def multiplicity(column, species, charge, oxidation, is_metal, scope):
    options = []
    if len(species.metals) == 1 and species.metals[0] in oxidation:
        options = [
            o for o in metals.spin_options(species.metals[0], oxidation[species.metals[0]])
            if metals.parity_problem(species.elements, charge, o.multiplicity) is None
        ]
    if options:
        # Only states the electron count allows are offered at all.
        choice = column.selectbox(
            "Spin state", options, index=None if is_metal else 0,
            format_func=lambda o: f"{o.multiplicity} - {o.label} ({o.unpaired} unpaired)",
            placeholder="choose", key=f"spin:{scope}:{charge}",
        )
        return choice.multiplicity if choice else None
    lowest = crystal_jobs.default_multiplicity(species, charge, oxidation)
    value = column.number_input(
        "Multiplicity", value=None if is_metal else lowest, min_value=1, step=1,
        placeholder=f"lowest allowed: {lowest}", key=f"mult:{scope}",
    )
    return int(value) if value is not None else None


def show_result(stored: dict, fmt: str, earlier: list) -> None:
    files = [Path(p) for p in stored["files"]]
    inputs = [p for p in files if not p.name.endswith(".m2i.json")]
    st.success(f"{FORMATS[fmt]} ready. A copy is in {inputs[0].parent}.")
    seen = {(i.code, i.message) for i in earlier}
    for issue in stored["issues"]:
        if issue.level == WARNING and (issue.code, issue.message) not in seen:
            st.warning(issue.message)
    for path in inputs:
        data = path.read_bytes()
        st.download_button(
            f"Download {path.name}", data, file_name=path.name,
            key=f"download:{path.name}", type="primary",
        )
        with st.expander(f"Preview {path.name}"):
            st.code(data.decode("utf-8", errors="replace"), language="text")
    for path in files:
        if path.name.endswith(".m2i.json"):
            st.download_button(
                f"Provenance record ({path.name})", path.read_bytes(),
                file_name=path.name, key=f"download:{path.name}",
            )


def viewer(species, height: int = 380, key: str = "crystal_viewer") -> None:
    """The molecule in 3D, metals labelled."""
    coords = crystal_jobs.centred(species)
    xyz = f"{species.n_atoms}\n{species.formula}\n" + "\n".join(
        f"{e} {x:.4f} {y:.4f} {z:.4f}" for e, (x, y, z) in zip(species.elements, coords)
    )
    labels = [
        {"text": species.labels[i], "x": float(x), "y": float(y), "z": float(z)}
        for i, (x, y, z) in enumerate(coords)
        if species.elements[i] in species.metals
    ]
    _VIEWER(xyz=xyz, labels=labels, height=height, key=key, default=None)
