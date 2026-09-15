"""A metal complex drawn in ChemDraw (or as a molfile), in the browser.

The drawing gives the connectivity but not the geometry: which donors are
trans to which is read from the bond directions on the page and the wedges,
offered as a list of possible arrangements, best fit first, and built in 3D
once chosen. What a drawing cannot say -- oxidation state, charge, spin -- is
asked, as on the crystal page.
"""

from __future__ import annotations

import json
from pathlib import Path

import streamlit as st

from m2i.chem import metals
from m2i.crystal import jobs as crystal_jobs
from m2i.gui import crystal_page
from m2i.organometallic import arrangements, build, mirror_partners, read
from m2i.organometallic import jobs as om_jobs
from m2i.recognition.base import BackendError
from m2i.types import IssueLog

FORMATS = crystal_page.FORMATS
#: Readings and built geometries kept per session. Each is a whole molecule,
#: and a browser session can walk through many.
CACHE_SIZE = 4


def render(settings: dict, *, path: Path, digest: str, shown_name: str, format_settings) -> None:
    base = _read(path, digest, {})
    if base is None:
        return
    substituents = _substituents(base, digest)
    drawing = _read(path, digest, substituents) if substituents else base
    if drawing is None:
        return
    st.markdown(
        f"**Metal complex:** {drawing.metal_symbol} with {len(drawing.donors)} donor atoms - "
        + ", ".join(drawing.names[d] for d in drawing.donors)
    )
    for level, _, message in drawing.notes:
        (st.warning if level == "warning" else st.caption)(message)

    options = arrangements(drawing.donors, drawing.drawn, drawing.chelated, labels=drawing.labels)
    if not options:
        st.error(
            f"{len(drawing.donors)} atoms are bound to {drawing.metal_symbol}; complexes with 2 "
            "to 6 donors can be built from a drawing. Is one of those bonds drawn by mistake? "
            "A bare line ending at the metal is read as a methyl."
        )
        return

    st.subheader("2. Check the structure")
    mirrors = mirror_partners(options, drawing.labels, drawing.chelated)

    def label(i: int) -> str:
        twin = f"; mirror image of {mirrors[i] + 1}" if i in mirrors else ""
        return f"{i + 1}. {om_jobs.describe(options[i], drawing)}{twin}  (fit {options[i].misfit:.2f})"

    if om_jobs.ambiguous(options):
        if om_jobs.mirror_only(options, drawing):
            st.warning(
                "The first two are mirror images, and the drawing does not tell them apart: "
                "only wedges and hashes on the bonds to the metal do. Choose the enantiomer."
            )
        else:
            st.warning("The drawing fits the first two about equally well: choose one.")
    shown = list(range(min(len(options), 12)))
    choice = st.radio(
        "Arrangement of the ligands (best fit to the drawing first)", shown,
        format_func=label, key=f"om_isomer:{digest}",
    )

    built = _build(drawing, options[choice], f"{digest}:{sorted(substituents.items())}",
                   choice, settings["conformers"].seed)
    if built is None:
        return
    species = built.species
    left, right = st.columns(2)
    with left:
        crystal_page.viewer(species, key="organometallic_viewer")
        st.caption(
            f"{species.formula}, built by distance geometry and pre-optimised with UFF: a "
            "starting geometry for your optimisation."
        )
    with right:
        crystal_page.coordination(species)
    for _, _, message in built.notes:
        st.warning(message)

    st.subheader("3. Electronic state")
    state = _electronic_state(drawing, species, f"{digest}:{choice}")

    st.subheader("4. Output")
    fmt = st.segmented_control(
        "Format", list(FORMATS), format_func=FORMATS.get, default="gaussian",
        required=True, key="om_format",
    )
    profile = format_settings(fmt)
    if profile is None or state is None:
        return
    charge, multiplicity, oxidation = state
    signature = (digest, choice, charge, multiplicity, json.dumps(oxidation), repr(profile),
                 str(settings["output_dir"]), settings["conformers"].seed,
                 json.dumps(substituents, sort_keys=True))
    if st.button(f"Generate {FORMATS[fmt]}", type="primary"):
        log = IssueLog()
        try:
            written = om_jobs.write(
                built, drawing, source=Path(shown_name), charge=charge,
                multiplicity=multiplicity, profile=profile,
                output_dir=settings["output_dir"], log=log, oxidation=oxidation,
            )
        except BackendError as exc:
            st.error(str(exc))
            return
        st.session_state["om_result"] = {"signature": signature, "files": written, "issues": list(log)}

    stored = st.session_state.get("om_result")
    if stored and stored["signature"] == signature:
        crystal_page.show_result(stored, fmt, [])


def _read(path: Path, digest: str, substituents: dict[str, str]):
    cache = st.session_state.setdefault("om_drawings", {})
    key = (digest, tuple(sorted(substituents.items())))
    if key not in cache:
        try:
            cache[key] = read(path, substituents)
        except BackendError as exc:
            st.error(str(exc))
            return None
        _trim(cache)
    return cache[key]


def _trim(cache: dict, keep: int = CACHE_SIZE) -> None:
    """Keep the last few readings only: a session must not grow without end."""
    for key in list(cache)[:-keep]:
        del cache[key]


def _substituents(drawing, digest: str) -> dict[str, str]:
    """The groups a drawing leaves out, asked for in a table, hydrogen by default.

    A sketch normally shows P rather than PiPr2, and hydrogen is what the
    valence needs but almost never what was meant, so the gaps are listed
    with the hydrogens filled in and can be written over.
    """
    gaps = {**drawing.assumed, **drawing.completed}
    if not gaps:
        return {}
    st.markdown("**What the drawing leaves out**")
    st.caption(
        "Write the real groups - iPr2, Ph2, Cy2, Me, OMe, CH2OH, or SMILES - one per "
        "missing bond. Left empty, the valence is filled with hydrogen, which is "
        "hardly ever what a scheme means."
    )
    rows = [
        {"Atom": name, "Bonds missing": len(added), "Groups": ""}
        for name, added in gaps.items()
    ]
    edited = st.data_editor(
        rows, key=f"om_subs:{digest}", hide_index=True, width="stretch",
        disabled=["Atom", "Bonds missing"],
        column_config={
            "Groups": st.column_config.TextColumn(
                "Groups", help="Blank leaves hydrogens", default="",
            ),
        },
    )
    return {
        row["Atom"]: str(row["Groups"]).strip()
        for row in edited
        if str(row["Groups"]).strip()
    }


def _build(drawing, arrangement, digest: str, choice: int, seed: int):
    cache = st.session_state.setdefault("om_built", {})
    key = (digest, choice, seed)
    if key not in cache:
        with st.spinner("Building the complex in 3D"):
            try:
                cache[key] = build(drawing, arrangement, seed=seed)
            except BackendError as exc:
                st.error(str(exc))
                return None
        _trim(cache)
    return cache[key]


def _electronic_state(drawing, species, scope: str):
    """Oxidation state (optional, for the spin options), charge and spin."""
    oxidation: dict[str, int] = {}
    symbol = drawing.metal_symbol
    value = st.number_input(
        f"Oxidation state of {symbol}", value=None, step=1, placeholder="e.g. 1",
        key=f"om_ox:{scope}",
        help="Only used to offer the spin states that make sense for the metal.",
    )
    if value is not None:
        oxidation[symbol] = int(value)
    for line in crystal_jobs.spin_hint(species, oxidation):
        st.caption(line)

    left, right = st.columns(2)
    charge = left.number_input(
        "Charge", value=drawing.drawn_charge, step=1, key=f"om_charge:{scope}",
    )
    left.caption(
        f"Suggestion: the formal charges drawn add up to {drawing.drawn_charge:+d}; a drawing "
        "often leaves charges out, so check it."
    )
    charge = int(charge)
    mult = crystal_page.multiplicity(right, species, charge, oxidation, True, f"om:{scope}")
    if mult is None:
        st.warning("Choose the spin state to continue.")
        return None
    problem = metals.parity_problem(species.elements, charge, mult)
    if problem:
        st.error(f"Impossible electronic state: {problem}.")
        return None
    return charge, mult, oxidation
