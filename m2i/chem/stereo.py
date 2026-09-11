"""Stereochemistry: CIP labelling and detection of undefined stereocentres.

The important distinction this module makes is between *"the drawing did not
define this centre"* and *"the recognition lost it"*. Both look identical in a
SMILES string, and silently picking one enantiomer is the worst thing this
program could do.
"""

from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ..types import IssueLog, StereoBond, StereoCenter, StereoSummary

_ATOM_TYPES = {
    Chem.StereoType.Atom_Tetrahedral,
    Chem.StereoType.Atom_SquarePlanar,
    Chem.StereoType.Atom_TrigonalBipyramidal,
    Chem.StereoType.Atom_Octahedral,
}
_BOND_TYPES = {
    Chem.StereoType.Bond_Double,
    Chem.StereoType.Bond_Cumulene_Even,
    Chem.StereoType.Bond_Atropisomer,
}


def assign(mol: Chem.Mol) -> Chem.Mol:
    """Perceive stereochemistry and attach true CIP descriptors."""
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    try:
        rdCIPLabeler.AssignCIPLabels(mol)
    except Exception:
        # The rigorous CIP implementation can bail out on exotic systems; the
        # legacy labels from AssignStereochemistry are still attached.
        pass
    return mol


def assign_from_3d(mol: Chem.Mol) -> Chem.Mol:
    """Re-derive stereochemistry from a 3D conformer, then CIP-label it."""
    Chem.AssignStereochemistryFrom3D(mol)
    try:
        rdCIPLabeler.AssignCIPLabels(mol)
    except Exception:
        pass
    return mol


def summarize(mol: Chem.Mol) -> StereoSummary:
    """List every potentially stereogenic element and whether it is defined."""
    summary = StereoSummary()
    for element in Chem.FindPotentialStereo(mol):
        specified = element.specified == Chem.StereoSpecified.Specified
        if element.type in _ATOM_TYPES:
            atom = mol.GetAtomWithIdx(element.centeredOn)
            center = StereoCenter(
                atom_index=element.centeredOn,
                symbol=atom.GetSymbol(),
                label=_atom_label(atom) if specified else None,
                specified=specified,
            )
            if is_equivalent_oxoanion(atom):
                summary.ignored_centers.append(center)
            else:
                summary.centers.append(center)
        elif element.type in _BOND_TYPES:
            bond = mol.GetBondWithIdx(element.centeredOn)
            summary.bonds.append(
                StereoBond(
                    begin_index=bond.GetBeginAtomIdx(),
                    end_index=bond.GetEndAtomIdx(),
                    label=_bond_label(bond) if specified else None,
                    specified=specified,
                )
            )
    return summary


#: Centres where a chalcogen pair can make the atom *look* stereogenic.
_OXOANION_CENTERS = {15, 16, 33, 34}  # P, S, As, Se
_CHALCOGENS = {8, 16, 34}  # O, S, Se


def is_equivalent_oxoanion(atom: Chem.Atom) -> bool:
    """True if this centre only looks stereogenic because of a P=O / P-O(-) pair.

    A phosphodiester is written with one double-bonded oxygen and one anionic
    oxygen, so RDKit sees four different substituents and flags the phosphorus.
    Chemically the charge is delocalised over both: the two oxygens are the same
    thing, the centre is not stereogenic, and the DNA backbone is famously not
    chiral at phosphorus. Left unfiltered this warning fires on every nucleotide,
    ATP, phospholipid and phosphate ester -- teaching the user to ignore exactly
    the warnings this program exists to give.

    The discriminator is the element. In a phosphorothioate the pair is =S and
    O(-), which resonance does *not* interchange, the centre is genuinely
    stereogenic, and its Sp/Rp isomers matter -- so it is not filtered here.
    A hydroxyl counts alongside the anion: P(=O)(OH) differs only by where a
    rapidly exchanging proton sits, which is not a resolvable configuration.
    """
    if atom.GetAtomicNum() not in _OXOANION_CENTERS:
        return False

    doubly_bonded: set[int] = set()
    exchangeable: set[int] = set()
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if other.GetAtomicNum() not in _CHALCOGENS or other.GetDegree() != 1:
            continue
        if bond.GetBondType() == Chem.BondType.DOUBLE and other.GetFormalCharge() == 0:
            doubly_bonded.add(other.GetAtomicNum())
        elif bond.GetBondType() == Chem.BondType.SINGLE:
            if other.GetFormalCharge() == -1 or other.GetTotalNumHs() == 1:
                exchangeable.add(other.GetAtomicNum())

    # Only the *same* element on both sides makes the two positions equivalent.
    return bool(doubly_bonded & exchangeable)


def report(summary: StereoSummary, log: IssueLog) -> None:
    """Turn the summary into user-facing issues."""
    if summary.centers or summary.bonds:
        log.info("stereo.summary", f"Stereochemistry: {summary.describe()}")

    if summary.ignored_centers:
        where = ", ".join(
            f"{c.symbol}{c.atom_index + 1}" for c in summary.ignored_centers
        )
        log.info(
            "stereo.equivalent_oxoanion",
            f"{where}: not treated as a stereocentre. The doubly bonded and "
            "anionic oxygens are equivalent by resonance, so the configuration "
            "is not real (a phosphorothioate, where they are not equivalent, "
            "would still be flagged).",
        )

    unspecified_centers = summary.unspecified_centers
    if unspecified_centers:
        where = ", ".join(f"{c.symbol}{c.atom_index + 1}" for c in unspecified_centers)
        log.warn(
            "stereo.unspecified_center",
            f"{len(unspecified_centers)} stereocentre(s) are not defined ({where}). "
            "The drawing did not specify them, or the recognition lost the wedge. "
            "The 3D structure will pick one configuration arbitrarily -- fix the "
            "SMILES with @/@@ if a particular enantiomer is intended.",
        )

    unspecified_bonds = summary.unspecified_bonds
    if unspecified_bonds:
        where = ", ".join(
            f"{b.begin_index + 1}={b.end_index + 1}" for b in unspecified_bonds
        )
        log.warn(
            "stereo.unspecified_bond",
            f"{len(unspecified_bonds)} double bond(s) have undefined E/Z geometry "
            f"({where}). The 3D structure will pick one arbitrarily -- use / and \\ "
            "in the SMILES to fix it.",
        )


def compare(reference: StereoSummary, candidate: StereoSummary) -> list[str]:
    """Descriptors that changed between the 2D perception and the 3D geometry.

    An empty list means the embedding preserved the stereochemistry. Elements
    that were undefined in the reference are ignored: they had no configuration
    to preserve.
    """
    ref = reference.fingerprint()
    cand = candidate.fingerprint()
    mismatches = []
    for key, expected in ref.items():
        found = cand.get(key)
        if found is None:
            mismatches.append(f"{_pretty(key)}: {expected} -> lost")
        elif found != expected:
            mismatches.append(f"{_pretty(key)}: {expected} -> {found}")
    return mismatches


def _pretty(key: str) -> str:
    kind, _, rest = key.partition(":")
    if kind == "atom":
        return f"atom {int(rest) + 1}"
    begin, _, end = rest.partition("-")
    return f"bond {int(begin) + 1}={int(end) + 1}"


def _atom_label(atom: Chem.Atom) -> str | None:
    if atom.HasProp("_CIPCode"):
        return atom.GetProp("_CIPCode")
    return None


def _bond_label(bond: Chem.Bond) -> str | None:
    if bond.HasProp("_CIPCode"):
        return bond.GetProp("_CIPCode")
    # Fallback only reached if AssignCIPLabels failed. STEREOE/STEREOZ are
    # already CIP-based; STEREOCIS/STEREOTRANS are relative to RDKit's chosen
    # reference atoms and must not be reported as E/Z.
    stereo = bond.GetStereo()
    if stereo == Chem.BondStereo.STEREOE:
        return "E"
    if stereo == Chem.BondStereo.STEREOZ:
        return "Z"
    return None
