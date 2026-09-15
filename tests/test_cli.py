"""Command line behaviour, including the verification gate."""

from __future__ import annotations

import json

import pytest

from m2i import cli
from m2i.types import RecognitionResult


def run(*args) -> int:
    return cli.main(list(args))


def test_from_smiles_writes_the_expected_files(tmp_path, capsys):
    code = run("from-smiles", "C[C@H](N)C(=O)O", "-o", str(tmp_path), "--yes")
    assert code == 0

    gjf = list(tmp_path.glob("*.gjf"))
    assert len(gjf) == 1
    text = gjf[0].read_text()
    assert text.startswith("%chk=")
    assert "\n0 1\n" in text

    record = json.loads(list(tmp_path.glob("*.m2i.json"))[0].read_text())
    assert record["molecule"]["smiles"] == "C[C@H](N)C(=O)O"
    assert record["issues"]

    assert "C[C@H](N)C(=O)O" in capsys.readouterr().out


def test_a_typed_structure_is_written_without_asking(tmp_path, monkeypatch):
    """Nothing was recognised, so there is nothing to confirm."""
    monkeypatch.setattr("m2i.cli.stdin_is_terminal", lambda: True)
    monkeypatch.setattr("builtins.input", lambda _: pytest.fail("asked about a typed SMILES"))
    assert run("from-smiles", "C[C@H](N)C(=O)O", "-o", str(tmp_path)) == 0
    record = json.loads(next(tmp_path.glob("*.m2i.json")).read_text())
    assert "review.not_needed" in {issue["code"] for issue in record["issues"]}


def read_picture_as(monkeypatch, tmp_path, smiles, confidence):
    """from-image with DECIMER replaced by a reading chosen by the test."""
    from PIL import Image

    image = tmp_path / "photo.png"
    Image.new("RGB", (300, 300), "white").save(image)
    reading = RecognitionResult(smiles=smiles, confidence=confidence, backend="decimer")
    monkeypatch.setattr("m2i.cli.resolve_backends", lambda names, log: ["decimer"])
    monkeypatch.setattr("m2i.cli.recognize", lambda *args, **kwargs: reading)
    return image


def test_a_doubtful_picture_is_not_written_without_a_terminal(tmp_path, monkeypatch):
    monkeypatch.setattr("sys.stdin.isatty", lambda: False, raising=False)
    image = read_picture_as(monkeypatch, tmp_path, "CCO", confidence=0.5)
    assert run("from-image", str(image), "-o", str(tmp_path)) == 130
    assert list(tmp_path.glob("*.gjf")) == []


def test_stereochemistry_read_from_a_picture_is_always_shown(tmp_path, monkeypatch, capsys):
    """DECIMER writes stereocentres instead of measuring them: even a confident
    reading with one is asked about."""
    monkeypatch.setattr("m2i.cli.stdin_is_terminal", lambda: True)
    asked = []
    monkeypatch.setattr("builtins.input", lambda prompt: asked.append(prompt) or "y")
    image = read_picture_as(monkeypatch, tmp_path, "C[C@H](N)C(=O)O", confidence=0.99)
    assert run("from-image", str(image), "-o", str(tmp_path)) == 0
    assert asked
    assert "written by the model" in capsys.readouterr().out
    record = json.loads(next(tmp_path.glob("*.m2i.json")).read_text())
    assert "review.confirmed" in {issue["code"] for issue in record["issues"]}


def test_a_refusal_writes_nothing(tmp_path, monkeypatch):
    monkeypatch.setattr("m2i.cli.stdin_is_terminal", lambda: True)
    monkeypatch.setattr("builtins.input", lambda prompt: "n")
    image = read_picture_as(monkeypatch, tmp_path, "CCO", confidence=0.4)
    assert run("from-image", str(image), "-o", str(tmp_path)) == 130
    assert list(tmp_path.glob("*.gjf")) == []


def test_a_closed_stdin_is_no_answer_rather_than_a_crash(tmp_path, monkeypatch):
    """On Windows, input redirected from NUL claims to be a terminal."""
    monkeypatch.setattr("m2i.cli.stdin_is_terminal", lambda: True)

    def closed(_prompt):
        raise EOFError

    monkeypatch.setattr("builtins.input", closed)
    image = read_picture_as(monkeypatch, tmp_path, "CCO", confidence=0.4)
    assert run("from-image", str(image), "-o", str(tmp_path)) == 130
    assert list(tmp_path.glob("*.gjf")) == []


def test_a_confident_plain_reading_goes_straight_through(tmp_path, monkeypatch):
    monkeypatch.setattr("m2i.cli.stdin_is_terminal", lambda: True)
    monkeypatch.setattr("builtins.input", lambda _: pytest.fail("asked about a clear reading"))
    image = read_picture_as(monkeypatch, tmp_path, "c1ccccc1O", confidence=0.97)
    assert run("from-image", str(image), "-o", str(tmp_path)) == 0
    assert list(tmp_path.glob("*.gjf"))


def test_from_cif_asks_about_the_metal_and_uses_the_answers(tmp_path, monkeypatch, capsys):
    from pathlib import Path

    import pytest

    pytest.importorskip("gemmi")

    monkeypatch.setattr("m2i.cli.stdin_is_terminal", lambda: True)
    answers = iter(["3", "1", "2"])  # Fe oxidation state, charge, multiplicity
    prompts = []
    monkeypatch.setattr("builtins.input", lambda p: prompts.append(p) or next(answers))

    ferrocene = Path(__file__).parent / "data" / "cod" / "7033930.cif"
    assert run("from-cif", str(ferrocene), "-o", str(tmp_path)) == 0

    assert any("Oxidation state of Fe" in p for p in prompts)
    assert "Fe(+3) is d5" in capsys.readouterr().out
    gjf = next(tmp_path.glob("*.gjf")).read_text(encoding="utf-8")
    assert "\n1 2\n" in gjf  # ferrocenium, doublet


def test_program_and_basis_overrides(tmp_path):
    run(
        "from-smiles",
        "O",
        "-o",
        str(tmp_path),
        "--program",
        "orca",
        "--basis",
        "def2-TZVP",
        "--yes",
    )
    text = list(tmp_path.glob("*.inp"))[0].read_text()
    assert "def2-TZVP" in text
    assert "* xyz 0 1" in text


def test_keep_writes_several_numbered_inputs(tmp_path):
    run("from-smiles", "CCCCCCO", "-o", str(tmp_path), "--keep", "3", "--yes")
    assert len(list(tmp_path.glob("*_c0*.gjf"))) == 3


def test_charge_and_multiplicity_flags(tmp_path):
    run(
        "from-smiles",
        "[CH3]",
        "-o",
        str(tmp_path),
        "--charge",
        "0",
        "--mult",
        "2",
        "--yes",
    )
    assert "\n0 2\n" in list(tmp_path.glob("*.gjf"))[0].read_text()


def test_impossible_spin_state_refuses_to_write(tmp_path, capsys):
    code = run("from-smiles", "CC", "-o", str(tmp_path), "--mult", "2", "--yes")
    assert code == 1
    assert list(tmp_path.glob("*.gjf")) == []
    captured = capsys.readouterr()
    assert "impossible" in (captured.out + captured.err).lower()


def test_invalid_smiles_exits_with_an_error(tmp_path, capsys):
    assert run("from-smiles", "NOTASMILES!!", "-o", str(tmp_path), "--yes") == 1
    assert "error" in capsys.readouterr().err.lower()


def test_from_molfile_round_trips_stereochemistry(tmp_path):
    from rdkit import Chem
    from rdkit.Chem import rdDepictor

    mol = Chem.MolFromSmiles("F[C@@H](Cl)Br")
    rdDepictor.Compute2DCoords(mol)
    Chem.WedgeMolBonds(mol, mol.GetConformer())
    molfile = tmp_path / "drawn.mol"
    molfile.write_text(Chem.MolToMolBlock(mol))

    assert run("from-molfile", str(molfile), "-o", str(tmp_path), "--yes") == 0
    record = json.loads(list(tmp_path.glob("*.m2i.json"))[0].read_text())
    assert record["molecule"]["stereo"]["centers"][0]["label"] == "S"


def test_batch_over_a_smiles_list(tmp_path):
    listing = tmp_path / "list.smi"
    listing.write_text(
        "CCO ethanol\nc1ccccc1 benzene\nCC(=O)[O-] acetate\nBROKEN!! bad\n"
    )
    out = tmp_path / "out"
    code = run("batch", str(listing), "-o", str(out), "--yes")
    assert code == 0  # partial failure is not a total failure

    import csv

    rows = list(csv.DictReader((out / "manifest.csv").open(encoding="utf-8")))
    by_name = {r["name"]: r for r in rows}
    assert by_name["ethanol"]["status"] == "ok"
    assert by_name["acetate"]["charge"] == "-1"
    assert by_name["bad"]["status"] == "failed"
    assert (out / "ethanol.gjf").is_file()


def test_batch_over_a_folder_of_molfiles(tmp_path):
    from rdkit import Chem
    from rdkit.Chem import rdDepictor

    source = tmp_path / "drawings"
    source.mkdir()
    for name, smiles in (("a", "CCO"), ("b", "CCN")):
        mol = Chem.MolFromSmiles(smiles)
        rdDepictor.Compute2DCoords(mol)
        (source / f"{name}.mol").write_text(Chem.MolToMolBlock(mol))

    out = tmp_path / "out"
    assert run("batch", str(source), "-o", str(out), "--yes") == 0
    assert (out / "a.gjf").is_file() and (out / "b.gjf").is_file()


def test_from_image_without_a_backend_explains_itself(tmp_path, capsys):
    from PIL import Image

    image = tmp_path / "drawing.png"
    Image.new("RGB", (400, 400), "white").save(image)

    code = run("from-image", str(image), "-o", str(tmp_path), "--yes")
    assert code == 1
    assert "--smiles" in capsys.readouterr().err


def test_from_image_with_a_typed_smiles(tmp_path):
    from PIL import Image

    image = tmp_path / "drawing.png"
    Image.new("RGB", (400, 400), "white").save(image)

    code = run(
        "from-image", str(image), "--smiles", "CCO", "-o", str(tmp_path), "--yes"
    )
    assert code == 0
    assert list(tmp_path.glob("*.gjf"))


def test_setup_list_states_the_cost_before_downloading(capsys, monkeypatch, tmp_path):
    monkeypatch.setenv("M2I_BACKEND_HOME", str(tmp_path))
    assert run("setup", "--list") == 0
    out = capsys.readouterr().out
    assert "decimer" in out and "molscribe" not in out
    assert "GB" in out  # the download size is stated up front
    assert "hand-drawn" in out


def test_setup_without_an_argument_lists_instead_of_guessing(capsys, monkeypatch, tmp_path):
    monkeypatch.setenv("M2I_BACKEND_HOME", str(tmp_path))
    assert run("setup") == 0
    assert "Install with" in capsys.readouterr().out


def test_setup_rejects_an_unknown_backend(capsys, monkeypatch, tmp_path):
    monkeypatch.setenv("M2I_BACKEND_HOME", str(tmp_path))
    assert run("setup", "nonesuch") == 2
    assert "unknown backend" in capsys.readouterr().err


def test_setup_remove_is_safe_when_nothing_is_installed(capsys, monkeypatch, tmp_path):
    monkeypatch.setenv("M2I_BACKEND_HOME", str(tmp_path))
    assert run("setup", "decimer", "--remove") == 0
    assert "Nothing to remove" in capsys.readouterr().out


def test_from_image_without_a_model_points_at_decimer_and_the_file(tmp_path, capsys, monkeypatch):
    """With nothing installed, the way out is the model for photos, or the
    ChemDraw file, which needs none."""
    monkeypatch.setenv("M2I_BACKEND_HOME", str(tmp_path))
    from PIL import Image

    image = tmp_path / "clean.png"
    Image.new("RGB", (400, 400), "white").save(image)

    assert run("from-image", str(image), "-o", str(tmp_path), "--yes") == 1
    err = capsys.readouterr().err
    assert "m2i setup decimer" in err
    assert "ChemDraw" in err


def test_doctor_and_profiles_run(capsys):
    assert run("doctor") == 0
    assert "Recognition backends" in capsys.readouterr().out
    assert run("profiles") == 0
    assert "gaussian_opt_freq" in capsys.readouterr().out
    assert run("profiles", "orca_opt_freq") == 0
    assert "def2-SVP" in capsys.readouterr().out


def test_unknown_profile_is_a_configuration_error(tmp_path, capsys):
    code = run("from-smiles", "CCO", "-p", "nope", "-o", str(tmp_path), "--yes")
    assert code == 2
    assert "not found" in capsys.readouterr().err


def test_batch_writes_whether_each_structure_needs_a_look(tmp_path):
    import csv

    listing = tmp_path / "list.smi"
    listing.write_text("CCO ethanol\n", encoding="utf-8")
    out = tmp_path / "out"
    assert run("batch", str(listing), "-o", str(out), "--yes") == 0
    row = next(csv.DictReader(open(out / "manifest.csv", encoding="utf-8")))
    assert row["review"] == "not needed"  # typed, not read from a picture


def test_batch_refuses_a_metal_complex_instead_of_misreading_it(tmp_path):
    """Through the organic route a complex comes out as [CH3][Pt]... without a word."""
    import csv
    import shutil
    from pathlib import Path

    source = tmp_path / "drawings"
    source.mkdir()
    shutil.copy(Path(__file__).parent / "data" / "drawings" / "bare_phosphine.cdxml", source)
    out = tmp_path / "out"
    assert run("batch", str(source), "-o", str(out), "--yes") != 0
    row = next(csv.DictReader(open(out / "manifest.csv", encoding="utf-8")))
    assert row["status"] == "failed"
    assert "from-molfile" in row["warning_detail"]
    assert not list(out.glob("*.gjf"))
