"""Turn whatever the recognition layer produced into a clean RDKit molecule.

The molblock is preferred over the SMILES whenever a backend provides one: it
carries 2D coordinates and wedge/hash bonds, which is how stereochemistry is
determined geometrically rather than guessed by a sequence model.
"""

from __future__ import annotations

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

from ..types import IssueLog, RecognitionResult

# RDKit prints parse failures to stderr; we surface them as Issues instead.
RDLogger.DisableLog("rdApp.*")


class MoleculeParseError(ValueError):
    """The input could not be turned into a valid molecule."""


def mol_from_smiles(smiles: str) -> Chem.Mol:
    if not smiles or not smiles.strip():
        raise MoleculeParseError("empty SMILES")
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        raise MoleculeParseError(_diagnose_smiles(smiles.strip()))
    return mol


def mol_from_molblock(molblock: str) -> Chem.Mol:
    if not molblock or not molblock.strip():
        raise MoleculeParseError("empty molblock")
    mol = Chem.MolFromMolBlock(molblock, sanitize=True, removeHs=True)
    if mol is None:
        # Retry unsanitised so we can report *which* atom is wrong.
        raw = Chem.MolFromMolBlock(molblock, sanitize=False, removeHs=False)
        if raw is None:
            raise MoleculeParseError("molblock could not be parsed at all")
        raise MoleculeParseError(_describe_problems(raw) or "molblock failed sanitisation")
    return mol


def mol_from_recognition(result: RecognitionResult, log: IssueLog) -> Chem.Mol:
    """Build a molecule from a recognition result, molblock first."""
    if result.molblock:
        try:
            mol = mol_from_molblock(result.molblock)
            log.info(
                "stereo.source.molblock",
                "Stereochemistry taken from the 2D layout and wedge bonds (molblock).",
            )
            return _reperceive_2d_stereo(mol)
        except MoleculeParseError as exc:
            log.warn(
                "parse.molblock_failed",
                f"Molblock from {result.backend} unusable ({exc}); falling back to the SMILES.",
            )
    mol = mol_from_smiles(result.smiles)
    if result.molblock is None and result.backend not in ("manual", "unknown"):
        log.info(
            "stereo.source.smiles",
            f"{result.backend} returned only a SMILES: stereochemistry comes from the "
            "generated token sequence, not from measuring the drawing.",
        )
    return mol


def _reperceive_2d_stereo(mol: Chem.Mol) -> Chem.Mol:
    """Re-derive double-bond stereo from the 2D coordinates when available.

    ``MolFromMolBlock`` already reads wedges into chiral tags, but a molblock
    written by an OCSR model sometimes carries crossed/undefined bond flags
    that are better resolved from the coordinates themselves.
    """
    if mol.GetNumConformers() == 0:
        return mol
    conf = mol.GetConformer()
    if conf.Is3D():
        return mol
    try:
        Chem.DetectBondStereochemistry(mol)
        Chem.AssignChiralTypesFromBondDirs(mol)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    except Exception:  # a malformed layout must not kill the run
        pass
    return mol


def sanitize(mol: Chem.Mol, log: IssueLog) -> Chem.Mol:
    """Sanitise in place-ish, reporting the offending atom when it fails."""
    problems = Chem.DetectChemistryProblems(mol)
    if problems:
        raise MoleculeParseError(_describe_problems(mol, problems))
    Chem.SanitizeMol(mol)
    return mol


def split_fragments(
    mol: Chem.Mol, log: IssueLog, *, keep_largest: bool = True
) -> Chem.Mol:
    """Handle salts and counter-ions.

    With ``keep_largest`` the biggest fragment is kept (the usual intent when a
    drawing includes a counter-ion); otherwise every fragment is kept and the
    calculation treats them as one supermolecular system.
    """
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if len(frags) <= 1:
        return mol

    smiles = [Chem.MolToSmiles(f) for f in frags]
    if not keep_largest:
        log.warn(
            "fragments.kept",
            f"{len(frags)} disconnected fragments kept as a single system: "
            f"{', '.join(smiles)}. Charge and multiplicity apply to the whole set.",
        )
        return mol

    largest = max(frags, key=lambda f: (f.GetNumHeavyAtoms(), f.GetNumAtoms()))
    dropped = [s for f, s in zip(frags, smiles) if f is not largest]
    log.warn(
        "fragments.dropped",
        f"{len(frags)} disconnected fragments found; keeping the largest "
        f"({Chem.MolToSmiles(largest)}) and dropping {', '.join(dropped)}. "
        "Use --keep-all-fragments to calculate the whole assembly instead.",
    )
    return largest


def canonical_smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)


def inchikey(mol: Chem.Mol, log: IssueLog | None = None) -> str:
    try:
        key = Chem.MolToInchiKey(mol)
    except Exception:
        key = ""
    if not key and log is not None:
        log.info("inchi.unavailable", "InChIKey could not be generated for this molecule.")
    return key


def formula(mol: Chem.Mol) -> str:
    return rdMolDescriptors.CalcMolFormula(mol)


def suggest_name(mol: Chem.Mol, fallback: str = "molecule") -> str:
    """A filesystem-safe identifier: the InChIKey skeleton when available."""
    key = inchikey(mol)
    return key.split("-")[0] if key else fallback


# -- diagnostics ---------------------------------------------------------


def _describe_problems(mol: Chem.Mol, problems=None) -> str:
    problems = problems if problems is not None else Chem.DetectChemistryProblems(mol)
    messages = []
    for problem in problems:
        message = problem.Message()
        try:
            idx = problem.GetAtomIdx()
            atom = mol.GetAtomWithIdx(idx)
            message = f"atom {idx + 1} ({atom.GetSymbol()}): {message}"
        except Exception:
            pass
        messages.append(message)
    return "; ".join(messages)


def _diagnose_smiles(smiles: str) -> str:
    """Give a more useful message than 'MolFromSmiles returned None'."""
    raw = Chem.MolFromSmiles(smiles, sanitize=False)
    if raw is None:
        return f"invalid SMILES syntax: {smiles!r}"
    described = _describe_problems(raw)
    return f"SMILES parsed but is chemically invalid: {described or smiles!r}"
