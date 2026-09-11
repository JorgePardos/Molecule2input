"""Get an image into the shape OCSR models expect.

Deliberately conservative: crop the whitespace, fix the orientation, and make
sure the drawing is not tiny. Aggressive binarisation or thinning tends to eat
the thin end of a hash wedge, which is exactly the information that decides
R from S.
"""

from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageChops, ImageOps

from ..types import IssueLog

#: Below this, recognition accuracy degrades badly.
MIN_SHORT_SIDE = 300
#: Above this there is nothing to gain and inference gets slow.
MAX_LONG_SIDE = 2000
#: Whitespace kept around the drawing after autocropping, in pixels.
MARGIN = 12

SUPPORTED_SUFFIXES = {".png", ".jpg", ".jpeg", ".bmp", ".gif", ".tif", ".tiff", ".webp"}


class ImageError(ValueError):
    """The image could not be read."""


def load_image(path: Path) -> Image.Image:
    path = Path(path)
    if not path.is_file():
        raise ImageError(f"no such image: {path}")
    try:
        image = Image.open(path)
        image.load()
    except Exception as exc:
        raise ImageError(f"{path.name}: {exc}") from exc
    return image


def prepare_image(
    path: Path, log: IssueLog | None = None, *, save_to: Path | None = None
) -> Image.Image:
    """Return a cleaned-up RGB copy, optionally saving it for inspection."""
    image = load_image(path)
    image = ImageOps.exif_transpose(image)

    if image.mode in ("RGBA", "LA", "P"):
        # Flatten transparency onto white; a transparent background otherwise
        # renders as black and swamps the drawing.
        image = image.convert("RGBA")
        background = Image.new("RGBA", image.size, (255, 255, 255, 255))
        image = Image.alpha_composite(background, image)
    image = image.convert("RGB")

    image = _autocrop(image)
    image = _rescale(image, log)

    if save_to is not None:
        save_to = Path(save_to)
        save_to.parent.mkdir(parents=True, exist_ok=True)
        image.save(save_to)
    return image


def estimate_drawing_style(image: Image.Image) -> tuple[str, dict]:
    """Guess whether this is a clean vector-style depiction or a hand drawing.

    Which backend to reach for depends entirely on this, and the two look very
    different in a histogram. Software exports have a mathematically pure white
    background and almost no mid-tones: every pixel is either paper or ink.
    A photographed or scanned drawing has a background that is *nearly* white
    but never exactly, plus a wide band of grey from paper texture, shadows and
    pencil pressure.

    Returns ``("clean" | "hand_drawn", metrics)``. It is a heuristic and it is
    reported to the user, never applied silently -- ``--backend`` overrides it.
    """
    grey = ImageOps.grayscale(image)
    histogram = grey.histogram()
    total = sum(histogram) or 1

    pure_white = histogram[255] / total
    light = sum(histogram[200:]) / total
    midtones = sum(histogram[60:200]) / total
    ink = sum(histogram[:60]) / total
    # Of everything that looks like background, how much is *exactly* white?
    background_purity = histogram[255] / (sum(histogram[200:]) or 1)

    metrics = {
        "pure_white_fraction": round(pure_white, 4),
        "light_fraction": round(light, 4),
        "midtone_fraction": round(midtones, 4),
        "ink_fraction": round(ink, 4),
        "background_purity": round(background_purity, 4),
    }

    clean = background_purity > 0.9 and midtones < 0.1
    metrics["style"] = "clean" if clean else "hand_drawn"
    return metrics["style"], metrics


def _autocrop(image: Image.Image) -> Image.Image:
    grey = ImageOps.grayscale(image)
    # Compare against a white canvas; the difference is the drawing.
    background = Image.new("L", grey.size, 255)
    bbox = ImageChops.difference(grey, background).getbbox()
    if not bbox:
        return image  # blank image, nothing to crop
    left, top, right, bottom = bbox
    return image.crop(
        (
            max(0, left - MARGIN),
            max(0, top - MARGIN),
            min(image.width, right + MARGIN),
            min(image.height, bottom + MARGIN),
        )
    )


def _rescale(image: Image.Image, log: IssueLog | None) -> Image.Image:
    short = min(image.size)
    long = max(image.size)

    if short < MIN_SHORT_SIDE:
        factor = MIN_SHORT_SIDE / short
        new_size = (round(image.width * factor), round(image.height * factor))
        if log is not None:
            log.info(
                "image.upscaled",
                f"Image upscaled from {image.width}x{image.height} to "
                f"{new_size[0]}x{new_size[1]}; a low-resolution drawing is the most "
                "common cause of a misread structure.",
            )
        return image.resize(new_size, Image.LANCZOS)

    if long > MAX_LONG_SIDE:
        factor = MAX_LONG_SIDE / long
        new_size = (round(image.width * factor), round(image.height * factor))
        return image.resize(new_size, Image.LANCZOS)

    return image
