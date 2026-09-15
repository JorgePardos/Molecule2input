"""The backend where the chemist is the recogniser.

Always available, always the fallback, and the reference implementation of the
backend contract. It is also what makes m2i useful before any vision model is
installed: draw the molecule, read it yourself, type the SMILES.
"""

from __future__ import annotations

from pathlib import Path

from ..types import RecognitionResult
from .base import BackendError


class ManualBackend:
    name = "manual"
    description = "SMILES (or molblock) supplied by the user"
    returns_molblock = True  # if the user hands one over

    def __init__(
        self,
        smiles: str | None = None,
        molblock: str | None = None,
        notes: list[tuple[str, str, str]] | None = None,
    ) -> None:
        self.smiles = smiles
        self.molblock = molblock
        #: (level, code, message) about how the input had to be interpreted.
        self.notes = list(notes or [])

    def available(self) -> tuple[bool, str]:
        return True, "always available"

    def recognize(self, image_path: Path | None = None, **options) -> RecognitionResult:
        if self.molblock:
            return RecognitionResult(
                smiles=self.smiles or "",
                molblock=self.molblock,
                confidence=1.0,
                backend=self.name,
                raw={"notes": self.notes} if self.notes else None,
            )
        if not self.smiles:
            raise BackendError(
                "the manual backend needs a --smiles value (or a molfile via "
                "--molfile); no vision model is installed yet"
            )
        return RecognitionResult(
            smiles=self.smiles, molblock=None, confidence=1.0, backend=self.name
        )


#: Files that already contain the structure, so no vision model is involved.
MOLFILE_SUFFIXES = (".mol", ".sdf", ".mdl")
CHEMDRAW_SUFFIXES = (".cdxml", ".cdx")
STRUCTURE_FILE_SUFFIXES = MOLFILE_SUFFIXES + CHEMDRAW_SUFFIXES


def from_structure_file(path: Path, smiles: str | None = None) -> ManualBackend:
    """Build a manual backend from a drawing file: molfile, SDF or ChemDraw.

    These are the best inputs m2i can get. The wedges are already in the file,
    so the stereochemistry is read exactly as drawn, with no model in between.
    """
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in CHEMDRAW_SUFFIXES:
        molblock, notes = _read_chemdraw(path)
        return ManualBackend(smiles=smiles, molblock=molblock, notes=notes)
    if suffix in MOLFILE_SUFFIXES:
        text = path.read_text(encoding="utf-8", errors="replace")
        # An SDF may hold several records; only the first is used.
        return ManualBackend(smiles=smiles, molblock=text.split("$$$$")[0])
    raise BackendError(
        f"{path.name}: unsupported structure file; expected one of "
        f"{', '.join(STRUCTURE_FILE_SUFFIXES)}"
    )


#: Kept for callers written before ChemDraw support existed.
from_molfile = from_structure_file


def _read_chemdraw(path: Path) -> tuple[str, list[tuple[str, str, str]]]:
    """Every structure in a ChemDraw document, as one molblock, plus notes.

    Three things about real ChemDraw files shape this function:

    - The bytes are handed to RDKit, never the path. RDKit's file readers use
      narrow-character paths on Windows, so a folder called "Artículo" makes
      the file "not exist".
    - The format is stated explicitly. Left to guess, RDKit parses a binary
      .cdx as XML and fails.
    - The page is read unsanitised and each molecule is then handled on its
      own. Reading sanitised makes RDKit drop, without a word, every molecule
      that fails -- on an active-site figure with carboxylates drawn as
      delocalised, that is every residue. Here a molecule that cannot be
      interpreted is reported by name instead.

    Several molecules are combined as disconnected fragments so the
    pipeline's fragment handling applies: it keeps the largest and names what
    it dropped, or keeps them all when asked.
    """
    from rdkit import Chem

    from ..chem.sanitize import resolve_delocalised_bonds

    suffix = path.suffix.lower()
    if suffix == ".cdx" and not Chem.HasChemDrawCDXSupport():
        raise BackendError(
            f"{path.name}: this RDKit build cannot read binary .cdx files. In "
            "ChemDraw, use File > Save As and choose .cdxml or .mol instead."
        )

    params = Chem.CDXMLParserParams()
    params.format = Chem.CDXMLFormat.CDX if suffix == ".cdx" else Chem.CDXMLFormat.CDXML
    params.sanitize = False
    params.removeHs = False
    try:
        found = Chem.MolsFromCDXML(path.read_bytes(), params)
        raw = [mol for mol in found if mol is not None and mol.GetNumAtoms()]
    except Exception as exc:  # noqa: BLE001 - RDKit raises assorted types here
        raise BackendError(f"{path.name}: could not read the ChemDraw file ({exc})") from exc

    # RDKit sometimes returns several drawn molecules as one; handle each
    # fragment on its own so a single bad one cannot take the others with it.
    fragments = [
        frag
        for mol in raw
        for frag in Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    ]

    notes: list[tuple[str, str, str]] = []
    kept, anions, resolved_as_drawn, failures = [], [], 0, []
    for mol in fragments:
        editable = Chem.RWMol(mol)
        groups = resolve_delocalised_bonds(editable)
        try:
            Chem.SanitizeMol(editable)
            clean = Chem.RemoveHs(editable.GetMol())
        except Exception as exc:  # noqa: BLE001
            failure = f"{Chem.MolToSmiles(mol)} ({exc})"
            failures.append(failure)
            notes.append((
                "warning",
                "chemdraw.molecule_unreadable",
                f"A molecule on the page could not be interpreted and was left out: "
                f"{failure}.",
            ))
            continue
        kept.append(clean)
        for group in groups:
            if group.added_charge:
                anions.append(Chem.MolToSmiles(clean))
            else:
                resolved_as_drawn += 1

    if not kept and failures:
        raise BackendError(
            f"{path.name}: none of the drawn molecules could be interpreted: "
            + "; ".join(failures)
        )
    if not kept:
        raise BackendError(
            f"{path.name}: no structure found. Text, arrows and shapes are ignored; "
            "the page needs at least one drawn molecule."
        )

    if anions:
        net = sum(Chem.GetFormalCharge(m) for m in kept)
        notes.append((
            "warning",
            "chemdraw.delocalised_anion",
            f"{len(anions)} group(s) drawn with delocalised bonds and no charge were "
            f"read as anions, -1 each, in: {', '.join(anions)}. The page now carries "
            f"a net charge of {net:+d}. If any of them is meant to be protonated, "
            "correct the SMILES or override the charge.",
        ))
    if resolved_as_drawn:
        notes.append((
            "info",
            "chemdraw.delocalised_resolved",
            f"{resolved_as_drawn} group(s) drawn with delocalised bonds were written "
            "as a Lewis structure, keeping the charges as drawn.",
        ))

    combined = kept[0]
    for mol in kept[1:]:
        combined = Chem.CombineMols(combined, mol)
    return Chem.MolToMolBlock(combined), notes
