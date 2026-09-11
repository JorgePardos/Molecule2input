"""Image preprocessing and the hand-drawn / clean heuristic that picks a model."""

from __future__ import annotations

import numpy as np
import pytest
from PIL import Image

from m2i.preprocess import estimate_drawing_style, prepare_image
from m2i.preprocess.image import MIN_SHORT_SIDE, ImageError
from m2i.types import IssueLog


def render_clean(path, smiles="C[C@H](O)/C=C/c1ccc(Cl)cc1", size=(600, 450)):
    """A software depiction: pure white background, no mid-tones."""
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D

    mol = Chem.MolFromSmiles(smiles)
    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DCairo(*size)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    path.write_bytes(drawer.GetDrawingText())
    return path


def simulate_photo(source: Image.Image, seed=0) -> Image.Image:
    """What a phone picture of a pencil drawing looks like to a histogram:
    paper that is nearly-but-never-exactly white, plus a broad band of grey."""
    rng = np.random.default_rng(seed)
    grey = np.asarray(source.convert("L"), dtype=float)
    grey = np.where(grey > 200, 244.0, grey)  # paper, not paper-white
    grey += rng.normal(0, 14, grey.shape)  # texture, shadow, pressure
    return Image.fromarray(np.clip(grey, 0, 255).astype("uint8")).convert("RGB")


def test_a_software_depiction_reads_as_clean(tmp_path):
    image = Image.open(render_clean(tmp_path / "clean.png"))
    style, metrics = estimate_drawing_style(image)
    assert style == "clean"
    assert metrics["background_purity"] > 0.9
    assert metrics["midtone_fraction"] < 0.1


def test_a_photographed_drawing_reads_as_hand_drawn(tmp_path):
    clean = Image.open(render_clean(tmp_path / "clean.png"))
    style, metrics = estimate_drawing_style(simulate_photo(clean))
    assert style == "hand_drawn"
    assert metrics["background_purity"] < 0.9


def test_the_heuristic_reports_its_evidence(tmp_path):
    """It is a guess, so it has to be inspectable rather than silent."""
    image = Image.open(render_clean(tmp_path / "clean.png"))
    _, metrics = estimate_drawing_style(image)
    assert set(metrics) >= {
        "pure_white_fraction",
        "midtone_fraction",
        "ink_fraction",
        "background_purity",
        "style",
    }


@pytest.mark.parametrize("seed", range(4))
def test_the_heuristic_is_stable_across_noise_samples(tmp_path, seed):
    clean = Image.open(render_clean(tmp_path / "clean.png"))
    assert estimate_drawing_style(simulate_photo(clean, seed=seed))[0] == "hand_drawn"


# -- preprocessing -------------------------------------------------------


def test_small_images_are_upscaled_with_a_warning(tmp_path):
    path = tmp_path / "tiny.png"
    Image.new("RGB", (120, 90), "white").save(path)
    log = IssueLog()
    result = prepare_image(path, log)
    assert min(result.size) >= MIN_SHORT_SIDE
    assert any(i.code == "image.upscaled" for i in log)


def test_whitespace_is_cropped_away(tmp_path):
    canvas = Image.new("RGB", (900, 900), "white")
    canvas.paste(Image.new("RGB", (100, 80), "black"), (400, 400))
    path = tmp_path / "sparse.png"
    canvas.save(path)

    result = prepare_image(path, IssueLog())
    assert result.width < 900 and result.height < 900


def test_transparency_is_flattened_onto_white(tmp_path):
    """A transparent PNG otherwise renders as black and swamps the drawing."""
    path = tmp_path / "alpha.png"
    canvas = Image.new("RGBA", (400, 400), (0, 0, 0, 0))
    canvas.paste(Image.new("RGBA", (60, 60), (0, 0, 0, 255)), (170, 170))
    canvas.save(path)

    result = prepare_image(path, IssueLog())
    assert result.mode == "RGB"
    corner = result.getpixel((1, 1))
    assert corner == (255, 255, 255)


def test_a_blank_image_does_not_crash_the_crop(tmp_path):
    path = tmp_path / "blank.png"
    Image.new("RGB", (500, 500), "white").save(path)
    assert prepare_image(path, IssueLog()).size == (500, 500)


def test_a_missing_image_is_reported(tmp_path):
    with pytest.raises(ImageError):
        prepare_image(tmp_path / "nope.png", IssueLog())


def test_a_file_that_is_not_an_image_is_reported(tmp_path):
    path = tmp_path / "notanimage.png"
    path.write_text("this is not a png")
    with pytest.raises(ImageError):
        prepare_image(path, IssueLog())
