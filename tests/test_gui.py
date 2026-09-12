"""The browser interface, driven headless with Streamlit's AppTest.

File uploads cannot be simulated this way, so the picture and ChemDraw inputs
are covered by their own tests; everything from the structure onwards --
checking, correcting, choosing a format, generating -- is exercised here.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest  # noqa: E402

APP = Path(__file__).resolve().parents[1] / "m2i" / "gui" / "app.py"
TARGET = "C[C@H](O)/C=C/c1ccc(Cl)cc1"


@pytest.fixture
def app(tmp_path):
    at = AppTest.from_file(str(APP), default_timeout=300).run()
    at.text_input(key="output_dir").set_value(str(tmp_path)).run()
    return at


def enter_smiles(at, smiles):
    at.segmented_control(key="source_kind").set_value("SMILES").run()
    at.text_input(key="smiles_text").input(smiles).run()
    return at


def choose_format(at, fmt):
    at.segmented_control(key="format").set_value(fmt).run()
    return at


def generate(at):
    next(b for b in at.button if b.label.startswith("Generate")).click().run()
    return at


def downloads(at) -> list[str]:
    return [d.proto.label.removeprefix("Download ") for d in at.get("download_button")]


def metrics(at) -> dict:
    return {m.label: m.value for m in at.metric}


def no_exceptions(at):
    assert not at.exception, [e.value for e in at.exception]


# -- opening -------------------------------------------------------------


def test_the_app_opens_with_a_prompt_and_nothing_else(app):
    no_exceptions(app)
    assert any("Start with a picture" in i.value for i in app.info)
    assert not app.metric


# -- the check step ------------------------------------------------------


def test_a_smiles_is_checked_before_anything_is_generated(app):
    enter_smiles(app, TARGET)
    no_exceptions(app)
    assert metrics(app) == {
        "Formula": "C10H11ClO",
        "Charge": "0",
        "Multiplicity": "1",
        "Stereocentres": "1",
    }
    assert any("C2: S; 4=5: E" in m.value for m in app.markdown)
    assert downloads(app) == []  # nothing is written until asked


def test_a_non_canonical_correction_settles_instead_of_looping(app):
    """A correction is compared with what was last applied, not with its own
    canonical form -- otherwise 'OCC' would differ from 'CCO' forever."""
    enter_smiles(app, "c1ccccc1")
    app.text_input(key="correction:smiles:c1ccccc1").input("OCC").run()
    no_exceptions(app)
    assert metrics(app)["Formula"] == "C2H6O"
    assert any(b.label == "Back to the original reading" for b in app.button)


def test_going_back_to_the_original_reading(app):
    enter_smiles(app, "c1ccccc1")
    app.text_input(key="correction:smiles:c1ccccc1").input("OCC").run()
    next(b for b in app.button if b.label == "Back to the original reading").click().run()
    no_exceptions(app)
    assert metrics(app)["Formula"] == "C6H6"


def test_an_unusable_reading_can_still_be_corrected(app):
    """A reading that fails to parse must leave the way out open."""
    enter_smiles(app, "C1CC")  # unclosed ring
    no_exceptions(app)
    assert app.error, "the problem should be shown"
    assert not app.metric

    app.text_input(key="correction:smiles:C1CC").input("C1CC1").run()
    no_exceptions(app)
    assert metrics(app)["Formula"] == "C3H6"


def test_the_charge_can_be_overridden(app):
    enter_smiles(app, "CC(=O)O")
    next(c for c in app.sidebar.checkbox if c.label == "Override the charge").check().run()
    app.sidebar.number_input[0].set_value(-1).run()
    no_exceptions(app)
    assert metrics(app)["Charge"] == "-1"


# -- choosing a format ---------------------------------------------------


@pytest.mark.parametrize("fmt, suffix", [
    ("gaussian", ".gjf"), ("orca", ".inp"), ("xyz", ".xyz"), ("sdf", ".sdf"),
])
def test_you_get_the_format_you_chose(app, fmt, suffix):
    enter_smiles(app, TARGET)
    generate(choose_format(app, fmt))
    no_exceptions(app)

    names = downloads(app)
    inputs = [n for n in names if not n.endswith(("_check.png", ".m2i.json"))]
    assert inputs == [f"FBQMXUBIFKPRNV{suffix}"]


def test_qm_formats_ask_for_a_method_and_geometry_formats_do_not(app):
    enter_smiles(app, "CCO")
    for fmt in ("gaussian", "orca"):
        choose_format(app, fmt)
        labels = {w.label for w in app.text_input}
        assert {"Method", "Basis set", "Memory"} <= labels, fmt
    for fmt in ("xyz", "sdf"):
        choose_format(app, fmt)
        labels = {w.label for w in app.text_input}
        assert "Method" not in labels and "Basis set" not in labels, fmt


def test_dispersion_choices_belong_to_the_chosen_program(app):
    enter_smiles(app, "CCO")
    choose_format(app, "gaussian")
    gaussian = next(s for s in app.selectbox if s.label == "Dispersion").options
    choose_format(app, "orca")
    orca = next(s for s in app.selectbox if s.label == "Dispersion").options
    assert "gd3bj" in gaussian and "gd3bj" not in orca
    assert "d3bj" in orca and "d3bj" not in gaussian


def test_the_method_typed_reaches_the_file(app, tmp_path):
    enter_smiles(app, "CCO")
    choose_format(app, "gaussian")
    next(w for w in app.text_input if w.label == "Method").set_value("wb97xd").run()
    next(w for w in app.text_input if w.label == "Basis set").set_value("def2TZVP").run()
    generate(app)
    no_exceptions(app)
    route = (tmp_path / "LFQSCWFLJHTTHZ.gjf").read_text(encoding="utf-8")
    assert "wb97xd/def2TZVP" in route


def test_the_geometry_is_embedded_once_for_every_format(app, monkeypatch):
    """Switching format rewrites a text file; it must not re-embed."""
    from m2i import pipeline

    calls = []
    real = pipeline.embed
    monkeypatch.setattr(pipeline, "embed", lambda *a, **k: calls.append(1) or real(*a, **k))

    enter_smiles(app, TARGET)
    for fmt in ("gaussian", "orca", "xyz", "sdf", "gaussian"):
        generate(choose_format(app, fmt))
    no_exceptions(app)
    assert len(calls) == 1


def test_changing_a_setting_hides_files_that_no_longer_match_it(app):
    """A download button must never offer a file generated with settings the
    screen no longer shows."""
    enter_smiles(app, "CCO")
    generate(choose_format(app, "gaussian"))
    assert downloads(app)

    next(w for w in app.text_input if w.label == "Method").set_value("pbe0").run()
    assert downloads(app) == []

    generate(app)
    assert downloads(app)


def test_a_new_structure_does_not_show_the_previous_files(app):
    enter_smiles(app, "CCO")
    generate(choose_format(app, "xyz"))
    assert downloads(app)

    app.text_input(key="smiles_text").input("CCN").run()
    assert downloads(app) == []


def test_the_crystal_page_opens_and_explains_itself(app):
    """Uploads cannot be simulated headless; the page itself must still render,
    and say up front that oxidation state, charge and spin will be asked."""
    app.segmented_control(key="source_kind").set_value("Crystal (.cif)").run()
    no_exceptions(app)
    assert any("you will be asked" in i.value for i in app.info)
