"""The web interface's API, driven through FastAPI's test client.

The front end is plain JavaScript that only calls these routes, so what a
visitor can do is what is tested here; the browser itself was checked by hand.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402

from m2i.types import RecognitionResult  # noqa: E402
from m2i.web import main, organic, warmup  # noqa: E402

DATA = Path(__file__).parent / "data"


@pytest.fixture
def client(monkeypatch, tmp_path):
    monkeypatch.setattr(warmup, "start", lambda: None)  # no model loading in tests
    monkeypatch.setenv("TMPDIR", str(tmp_path))
    with TestClient(main.app) as c:
        yield c


def give_smiles(client, smiles="C[C@H](N)C(=O)O"):
    answer = client.post("/api/source/smiles", json={"smiles": smiles})
    assert answer.status_code == 200, answer.text
    return answer.json()["key"]


# -- the page and the server ---------------------------------------------------------


def test_the_page_is_served_with_its_security_headers(client):
    page = client.get("/")
    assert page.status_code == 200
    assert "m2i" in page.text
    assert "default-src 'self'" in page.headers["content-security-policy"]
    assert page.headers["x-content-type-options"] == "nosniff"
    assert client.get("/app.js").status_code == 200
    assert client.get("/viewer/3Dmol-min.js").status_code == 200


def test_health_always_answers(client):
    health = client.get("/api/health").json()
    assert health["status"] in ("starting", "ready")
    assert {"smiles", "file", "cif"} <= set(health["routes_available"])


def test_about_names_the_accepted_types_and_the_limit(client):
    about = client.get("/api/about").json()
    assert "cdx" in about["types"]["file"] and "png" in about["types"]["picture"]
    assert about["max_upload_mb"] == 20


# -- organic molecules -----------------------------------------------------------------------


def test_a_smiles_goes_straight_to_the_output(client):
    key = give_smiles(client)
    check = client.get(f"/api/check/{key}").json()
    assert check["molecule"]["formula"] == "C3H7NO2"
    assert check["gate"]["needed"] is False
    assert client.get(check["depiction_url"]).headers["content-type"] == "image/png"

    result = client.post("/api/generate", json={"key": key, "settings": {"format": "orca"}}).json()
    names = [f["name"] for f in result["files"]]
    assert names and names[0].endswith(".inp")
    assert "B3LYP" in result["files"][0]["preview"]
    assert {e["kind"] for e in result["extras"]} == {"provenance", "check_image"}
    assert client.get(result["files"][0]["download_url"]).status_code == 200


def test_switching_format_reuses_the_geometry(client, monkeypatch):
    from m2i import pipeline

    calls = []
    real = pipeline.embed
    monkeypatch.setattr(pipeline, "embed", lambda *a, **k: calls.append(1) or real(*a, **k))
    key = give_smiles(client, "CCO")
    for fmt in ("gaussian", "orca", "xyz"):
        assert client.post("/api/generate", json={"key": key, "settings": {"format": fmt}}).status_code == 200
    assert len(calls) == 1


def test_a_doubtful_picture_is_held_until_confirmed(client, monkeypatch, tmp_path):
    from PIL import Image

    image = tmp_path / "photo.png"
    Image.new("RGB", (120, 90), "white").save(image)
    reading = RecognitionResult(smiles="C[C@H](N)C(=O)O", confidence=0.97, backend="decimer")
    monkeypatch.setattr(organic, "resolve_backends", lambda names, log: ["decimer"])
    monkeypatch.setattr(organic, "recognize", lambda *a, **k: reading)
    with image.open("rb") as fh:
        key = client.post("/api/source/image", files={"file": ("photo.png", fh, "image/png")}).json()["key"]
    assert client.post(f"/api/source/{key}/read").status_code == 200

    check = client.get(f"/api/check/{key}").json()
    assert check["gate"]["needed"] and not check["gate"]["confirmed"]
    assert "written by the model" in " ".join(check["gate"]["reasons"])
    blocked = client.post("/api/generate", json={"key": key, "settings": {"format": "xyz"}})
    assert blocked.status_code == 409

    smiles = check["molecule"]["smiles"]
    client.post(f"/api/check/{key}/confirm", json={"smiles": smiles})
    assert client.get(f"/api/check/{key}").json()["gate"]["confirmed"]
    assert client.post("/api/generate", json={"key": key, "settings": {"format": "xyz"}}).status_code == 200

    # A correction is another structure: it has to be confirmed again.
    client.post(f"/api/check/{key}/correct", json={"smiles": "C[C@@H](N)C(=O)O"})
    corrected = client.get(f"/api/check/{key}").json()
    assert corrected["source"]["corrected"] and not corrected["gate"]["confirmed"]


def test_an_unusable_smiles_is_explained_and_still_correctable(client):
    key = give_smiles(client, "C1CC")
    check = client.get(f"/api/check/{key}").json()
    assert check["molecule"] is None and check["error"]
    client.post(f"/api/check/{key}/correct", json={"smiles": "C1CC1"})
    assert client.get(f"/api/check/{key}").json()["molecule"]["formula"] == "C3H6"


# -- uploads and sessions ----------------------------------------------------------------------


def test_uploads_of_the_wrong_type_are_refused(client):
    answer = client.post("/api/source/file", files={"file": ("evil.exe", b"MZ", "application/octet-stream")})
    assert answer.status_code == 400
    assert "not accepted" in answer.json()["error"]


def test_an_oversized_upload_is_refused_before_it_is_read(client, monkeypatch):
    monkeypatch.setenv("M2I_MAX_UPLOAD_MB", "1")
    big = b"x" * (3 * 2**20)
    answer = client.post("/api/source/file", files={"file": ("big.mol", big, "chemical/x-mdl-molfile")})
    assert answer.status_code == 413


def test_the_uploaded_name_never_becomes_a_path(client):
    mol = (Path(__file__).parents[1] / "examples" / "drawn_molecule.mol").read_bytes()
    answer = client.post("/api/source/file", files={"file": ("../../etc/passwd.mol", mol, "text/plain")}).json()
    assert answer["route"] == "organic"
    assert answer["key"].startswith("file:") and "/" not in answer["key"]


def test_one_visitor_cannot_download_anothers_files(client):
    key = give_smiles(client, "CCO")
    result = client.post("/api/generate", json={"key": key, "settings": {"format": "xyz"}}).json()
    url = result["files"][0]["download_url"]
    with TestClient(main.app) as stranger:
        assert stranger.get(url).status_code == 404
        assert stranger.get(f"/api/check/{key}").status_code == 404


# -- metal complexes and crystals ---------------------------------------------------------------


def test_a_drawn_complex_is_built_with_the_groups_given(client):
    with (DATA / "drawings" / "bare_phosphine.cdxml").open("rb") as fh:
        answer = client.post("/api/source/file", files={"file": ("sketch.cdxml", fh)}).json()
    assert answer["route"] == "organometallic"
    key = answer["key"]
    page = client.post(f"/api/om/{key}/state", json={"substituents": {"P4": "Cy3", "P5": "Ph3"}}).json()
    assert page["built"]["formula"] == "C43H65ClP3Pt"
    assert {g["atom"] for g in page["drawing"]["gaps"]} == {"P4", "P5"}
    assert page["arrangements"][0]["label"].startswith("square planar")

    impossible = client.post(f"/api/om/{key}/generate", json={
        "substituents": {"P4": "Cy3", "P5": "Ph3"}, "charge": 0, "multiplicity": 1, "settings": {"format": "xyz"}})
    assert impossible.status_code == 400 and "Impossible" in impossible.json()["error"]
    written = client.post(f"/api/om/{key}/generate", json={
        "substituents": {"P4": "Cy3", "P5": "Ph3"}, "charge": 1, "multiplicity": 1, "settings": {"format": "gaussian"}})
    assert written.status_code == 200, written.text
    assert written.json()["files"][0]["name"].endswith(".gjf")


def test_a_crystal_is_rebuilt_and_its_hydrogens_put_back(client):
    with (DATA / "cod" / "1538403.cif").open("rb") as fh:
        summary = client.post("/api/source/cif", files={"file": ("1538403.cif", fh)}).json()
    assert summary["source"] == "COD 1538403"
    key = summary["key"]
    page = client.post(f"/api/crystal/{key}/state", json={"species": 0, "oxidation": {"Pt": 2}}).json()
    assert page["hydrogens"]["needed"] and page["hydrogens"]["matches"]
    assert page["species"]["formula"] == "Cl2H6N2Pt"
    assert page["electronic"]["options"][0]["multiplicity"] == 1
    assert client.post(f"/api/crystal/{key}/generate", json={
        "species": 0, "oxidation": {"Pt": 2}, "charge": 0, "multiplicity": 1, "settings": {"format": "sdf"}}).status_code == 400
    written = client.post(f"/api/crystal/{key}/generate", json={
        "species": 0, "oxidation": {"Pt": 2}, "charge": 0, "multiplicity": 1, "settings": {"format": "orca"}})
    assert written.status_code == 200, written.text
