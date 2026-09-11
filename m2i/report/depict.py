"""The verification picture: what was drawn, next to what m2i understood.

Everything hinges on this image being checkable at a glance, so the annotations
are built here rather than left to RDKit's defaults: atoms are numbered from 1
to match both the warning messages and the atom order in the generated input
file, and stereocentres the drawing failed to define are highlighted instead of
quietly rendered as flat.
"""

from __future__ import annotations

import io
from pathlib import Path

from PIL import Image, ImageDraw
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from ..chem import stereo as stereo_mod
from ..types import StereoSummary

PANEL = (520, 460)
LABEL_HEIGHT = 26
GAP = 10
BACKGROUND = (255, 255, 255)
LABEL_COLOR = (40, 40, 40)
#: Undefined stereochemistry, the thing most worth noticing.
ALERT = (1.0, 0.45, 0.45)


def draw_molecule(
    mol: Chem.Mol,
    size: tuple[int, int] = PANEL,
    *,
    stereo: StereoSummary | None = None,
    number_atoms: bool = True,
) -> bytes:
    """Render an annotated 2D depiction as PNG bytes."""
    mol = _flat_copy(mol)
    summary = stereo if stereo is not None else stereo_mod.summarize(mol)
    highlight_atoms, highlight_bonds = _annotate(mol, summary, number_atoms)

    drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    options = drawer.drawOptions()
    # Annotations are set as atom/bond notes above; letting RDKit add its own
    # would overwrite them and reintroduce 0-based numbering.
    options.addStereoAnnotation = False
    options.addAtomIndices = False
    options.annotationFontScale = 0.75
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer,
        mol,
        highlightAtoms=highlight_atoms,
        highlightAtomColors={i: ALERT for i in highlight_atoms},
        highlightBonds=highlight_bonds,
        highlightBondColors={i: ALERT for i in highlight_bonds},
    )
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def comparison_image(
    mol: Chem.Mol,
    out_path: Path,
    *,
    source_image: Path | Image.Image | None = None,
    caption: str = "",
    stereo: StereoSummary | None = None,
) -> Path:
    """Write a side-by-side PNG: the original drawing vs the parsed structure."""
    parsed = Image.open(io.BytesIO(draw_molecule(mol, stereo=stereo))).convert("RGB")

    panels: list[tuple[str, Image.Image]] = []
    if source_image is not None:
        original = (
            source_image
            if isinstance(source_image, Image.Image)
            else Image.open(source_image).convert("RGB")
        )
        panels.append(("Original drawing", _fit(original, PANEL)))
    panels.append(("Understood by m2i (atoms numbered as in the input file)", _fit(parsed, PANEL)))

    caption_height = LABEL_HEIGHT if caption else 0
    width = sum(p.width for _, p in panels) + GAP * (len(panels) + 1)
    height = PANEL[1] + LABEL_HEIGHT + caption_height + GAP * 2

    canvas = Image.new("RGB", (width, height), BACKGROUND)
    draw = ImageDraw.Draw(canvas)

    x = GAP
    for label, panel in panels:
        draw.text((x + 4, GAP), label, fill=LABEL_COLOR)
        canvas.paste(panel, (x, GAP + LABEL_HEIGHT))
        x += panel.width + GAP

    if caption:
        draw.text(
            (GAP + 4, GAP + LABEL_HEIGHT + PANEL[1] + 4),
            caption[:220],
            fill=LABEL_COLOR,
        )

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(out_path)
    return out_path


# -- internals -----------------------------------------------------------


def _flat_copy(mol: Chem.Mol) -> Chem.Mol:
    """A 2D copy safe to annotate without touching the caller's molecule."""
    mol = Chem.Mol(mol)
    if mol.GetNumConformers() and mol.GetConformer().Is3D():
        # A 3D conformer projects into an unreadable tangle; lay it out flat.
        try:
            mol = Chem.RemoveHs(mol)
        except Exception:
            pass
        mol.RemoveAllConformers()
    if mol.GetNumConformers() == 0:
        rdDepictor.Compute2DCoords(mol)
        rdDepictor.StraightenDepiction(mol)
    return mol


def _annotate(
    mol: Chem.Mol, summary: StereoSummary, number_atoms: bool
) -> tuple[list[int], list[int]]:
    """Attach 1-based numbers and CIP labels; return what to highlight."""
    notes: dict[int, list[str]] = {}
    if number_atoms:
        for atom in mol.GetAtoms():
            notes[atom.GetIdx()] = [str(atom.GetIdx() + 1)]

    highlight_atoms: list[int] = []
    for center in summary.centers:
        if center.atom_index >= mol.GetNumAtoms():
            continue
        tag = f"({center.label})" if center.specified and center.label else "(?)"
        notes.setdefault(center.atom_index, []).append(tag)
        if not center.specified:
            highlight_atoms.append(center.atom_index)

    for idx, parts in notes.items():
        mol.GetAtomWithIdx(idx).SetProp("atomNote", " ".join(parts))

    highlight_bonds: list[int] = []
    for stereo_bond in summary.bonds:
        bond = mol.GetBondBetweenAtoms(stereo_bond.begin_index, stereo_bond.end_index)
        if bond is None:
            continue
        if stereo_bond.specified and stereo_bond.label:
            bond.SetProp("bondNote", f"({stereo_bond.label})")
        else:
            bond.SetProp("bondNote", "(?)")
            highlight_bonds.append(bond.GetIdx())

    return highlight_atoms, highlight_bonds


def _fit(image: Image.Image, box: tuple[int, int]) -> Image.Image:
    """Letterbox an image into `box` without distorting it."""
    fitted = image.copy()
    fitted.thumbnail(box, Image.LANCZOS)
    panel = Image.new("RGB", box, BACKGROUND)
    panel.paste(fitted, ((box[0] - fitted.width) // 2, (box[1] - fitted.height) // 2))
    return panel
