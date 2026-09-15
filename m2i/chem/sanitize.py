"""Turn whatever the recognition layer produced into a clean RDKit molecule.

The molblock is preferred over the SMILES whenever a backend provides one: it
carries 2D coordinates and wedge/hash bonds, which is how stereochemistry is
determined geometrically rather than guessed by a sequence model.
"""

from __future__ import annotations

from dataclasses import dataclass

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
    # Whatever produced the result may have had to interpret the input -- a
    # ChemDraw reader resolving delocalised bonds, say. It reports that here.
    for level, code, message in (result.raw or {}).get("notes", ()):
        log.add(level, code, message)

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
            f"{', '.join(smiles)}. Charge and multiplicity apply to the whole set. "
            "Their placement relative to each other in 3D is arbitrary: a 2D "
            "drawing says nothing about it. For a cluster model of an active "
            "site, the coordinates have to come from the crystal structure or an "
            "MD frame.",
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


# -- bonds drawn as delocalised ------------------------------------------

#: Terminal atoms that carry the charge of a delocalised anion.
_CHALCOGENS = {8, 16, 34}  # O, S, Se


@dataclass
class DelocalisedGroup:
    """One group whose bonds were drawn as delocalised (bond order 1.5)."""

    central: int
    added_charge: int


def resolve_delocalised_bonds(mol: Chem.RWMol) -> list[DelocalisedGroup]:
    """Rewrite bonds drawn as delocalised outside a ring as a Lewis structure.

    Mechanism drawings often show a carboxylate or a guanidinium with two
    dashed, "one and a half" bonds. ChemDraw stores those as order 1.5, which
    RDKit reads as aromatic -- and an aromatic bond outside a ring has no valid
    Lewis structure, so the molecule fails sanitisation and, left alone, would
    simply vanish from the page.

    Each group gets one double bond: to a terminal atom drawn with a positive
    charge if there is one (the =NH2+ of a guanidinium), otherwise to a
    terminal chalcogen. A terminal O/S left singly bonded, with no hydrogen and
    no charge drawn, can only be the anion -- a delocalised drawing of a
    carboxylic acid would make no sense -- so it becomes O-. That changes the
    net charge, which is why every such group is returned for reporting.
    """
    mol.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(mol)
    bonds = [
        b
        for b in mol.GetBonds()
        if b.GetBondType() == Chem.BondType.AROMATIC and not b.IsInRing()
    ]
    groups: list[DelocalisedGroup] = []
    for group in _connected_bonds(bonds):
        degree: dict[int, int] = {}
        for bond in group:
            for idx in (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()):
                degree[idx] = degree.get(idx, 0) + 1
            bond.SetBondType(Chem.BondType.SINGLE)
            bond.SetIsAromatic(False)
        for idx in degree:
            atom = mol.GetAtomWithIdx(idx)
            if not any(b.GetIsAromatic() for b in atom.GetBonds()):
                atom.SetIsAromatic(False)

        central = max(degree, key=lambda idx: (degree[idx], -idx))
        centre = mol.GetAtomWithIdx(central)
        spokes = [b for b in group if central in (b.GetBeginAtomIdx(), b.GetEndAtomIdx())]

        def preference(bond, centre=centre):
            other = bond.GetOtherAtom(centre)
            return (other.GetFormalCharge() > 0, other.GetAtomicNum() in _CHALCOGENS)

        double = max(spokes, key=preference)
        double.SetBondType(Chem.BondType.DOUBLE)

        added = 0
        for bond in spokes:
            if bond is double:
                continue
            other = bond.GetOtherAtom(centre)
            if (
                other.GetAtomicNum() in _CHALCOGENS
                and other.GetDegree() == 1
                and other.GetFormalCharge() == 0
                and other.GetNumExplicitHs() == 0
            ):
                other.SetFormalCharge(-1)
                other.SetNoImplicit(True)
                added -= 1
        groups.append(DelocalisedGroup(central=central, added_charge=added))
    return groups


def _connected_bonds(bonds: list) -> list[list]:
    """Split bonds into groups that share atoms."""
    remaining = list(bonds)
    groups = []
    while remaining:
        group = [remaining.pop()]
        atoms = {group[0].GetBeginAtomIdx(), group[0].GetEndAtomIdx()}
        grew = True
        while grew:
            grew = False
            for bond in list(remaining):
                ends = {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}
                if ends & atoms:
                    group.append(bond)
                    atoms |= ends
                    remaining.remove(bond)
                    grew = True
        groups.append(group)
    return groups
