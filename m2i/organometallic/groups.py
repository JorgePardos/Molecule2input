"""Abbreviations: turning "PiPr2", "CH2OH" or "OTf" into atoms.

A drawing says what the chemist meant, not what the calculation needs. Three
things can be missing, and only the last of them needs asking:

- a nickname ChemDraw interpreted (CO, Ph, Ts): the expansion is stored in the
  file itself and RDKit builds it -- nothing to do here;
- a label ChemDraw did *not* interpret, which it keeps as plain text: the text
  is still exact, so it is expanded here rather than guessed;
- substituents never drawn at all (a phosphorus with a single line to the
  metal): no file holds them, so they are filled with hydrogen and said aloud,
  or given by the user.

The groups below are written with their attachment atom first, which is how a
label reads: PMe2 hangs two methyls on the P, CH2OH is a chain starting at C.
"""

from __future__ import annotations

import re
from functools import lru_cache

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Geometry import Point3D

from ..recognition.base import BackendError

#: name -> SMILES of the group, its first atom being the one that attaches.
GROUPS: dict[str, str] = {
    "h": "[H]",
    "me": "C", "et": "CC", "pr": "CCC", "npr": "CCC", "ipr": "C(C)C",
    "bu": "CCCC", "nbu": "CCCC", "ibu": "CC(C)C", "sbu": "C(C)CC", "tbu": "C(C)(C)C",
    "cy": "C1CCCCC1", "cyclohexyl": "C1CCCCC1", "cyclopentyl": "C1CCCC1",
    "cpr": "C1CC1", "cyclopropyl": "C1CC1", "ad": "C12CC3CC(CC(C3)C1)C2",
    "ph": "c1ccccc1", "bn": "Cc1ccccc1", "bz": "C(=O)c1ccccc1",
    "tol": "c1ccc(C)cc1", "ptol": "c1ccc(C)cc1", "mes": "c1c(C)cc(C)cc1C",
    "xyl": "c1c(C)cccc1C", "naph": "c1ccc2ccccc2c1", "c6f5": "c1c(F)c(F)c(F)c(F)c1F",
    "ch3": "C", "ch2": "C", "ch": "C", "c": "C",
    "vinyl": "C=C", "allyl": "CC=C", "cho": "C=O", "cn": "C#N",
    "ac": "C(C)=O", "cf3": "C(F)(F)F", "co2me": "C(=O)OC", "co2et": "C(=O)OCC",
    "cooh": "C(=O)O", "co2h": "C(=O)O",
    "oh": "O", "o": "O", "ome": "OC", "oet": "OCC", "oipr": "OC(C)C",
    "otbu": "OC(C)(C)C", "oph": "Oc1ccccc1", "obn": "OCc1ccccc1", "oac": "OC(C)=O",
    "otf": "OS(=O)(=O)C(F)(F)F", "ots": "OS(=O)(=O)c1ccc(C)cc1", "oms": "OS(C)(=O)=O",
    "nh2": "N", "nh": "N", "n": "N", "nme2": "N(C)C", "net2": "N(CC)CC",
    "nhme": "NC", "nhph": "Nc1ccccc1", "no2": "[N+](=O)[O-]", "n3": "N=[N+]=[N-]",
    "sh": "S", "sme": "SC", "set": "SCC", "sph": "Sc1ccccc1",
    "ts": "S(=O)(=O)c1ccc(C)cc1", "ms": "S(C)(=O)=O", "tf": "S(=O)(=O)C(F)(F)F",
    "sime3": "[Si](C)(C)C", "tms": "[Si](C)(C)C", "tbs": "[Si](C)(C)C(C)(C)C",
    "pph2": "P(c1ccccc1)c1ccccc1", "pme2": "P(C)C", "pipr2": "P(C(C)C)C(C)C",
    "f": "F", "cl": "Cl", "br": "Br", "i": "I",
}

#: How people actually write them, normalised to the keys above.
ALIASES = {
    "methyl": "me", "ethyl": "et", "propyl": "pr", "isopropyl": "ipr", "pri": "ipr",
    "isopr": "ipr", "butyl": "bu", "tertbutyl": "tbu", "tbutyl": "tbu",
    "phenyl": "ph", "benzyl": "bn", "benzoyl": "bz", "mesityl": "mes", "tolyl": "tol",
    "adamantyl": "ad", "hydrogen": "h", "hydride": "h",
    "methoxy": "ome", "ethoxy": "oet", "hydroxy": "oh", "triflate": "otf",
    "tosyl": "ts", "tosylate": "ots", "mesyl": "ms", "acetate": "oac", "acetyl": "ac",
    "trimethylsilyl": "tms", "cyano": "cn", "nitro": "no2", "azide": "n3",
    "pentafluorophenyl": "c6f5",
}

#: How a group is written back to the user.
DISPLAY = {
    "ipr": "iPr", "nbu": "nBu", "npr": "nPr", "ibu": "iBu", "sbu": "sBu", "tbu": "tBu",
    "ph": "Ph", "me": "Me", "et": "Et", "pr": "Pr", "bu": "Bu", "cy": "Cy", "bn": "Bn",
    "bz": "Bz", "mes": "Mes", "tol": "Tol", "ptol": "p-Tol", "xyl": "Xyl", "ad": "Ad",
    "naph": "Naph", "c6f5": "C6F5", "cf3": "CF3", "cn": "CN", "cho": "CHO", "ac": "Ac",
    "ch3": "CH3", "ch2": "CH2", "ch": "CH", "c": "C", "h": "H",
    "oh": "OH", "ome": "OMe", "oet": "OEt", "oipr": "OiPr", "otbu": "OtBu", "oph": "OPh",
    "obn": "OBn", "oac": "OAc", "otf": "OTf", "ots": "OTs", "oms": "OMs", "o": "O",
    "nh2": "NH2", "nh": "NH", "n": "N", "nme2": "NMe2", "net2": "NEt2", "nhme": "NHMe",
    "nhph": "NHPh", "no2": "NO2", "n3": "N3", "sh": "SH", "sme": "SMe", "set": "SEt",
    "sph": "SPh", "ts": "Ts", "ms": "Ms", "tf": "Tf", "tms": "TMS", "tbs": "TBS",
    "sime3": "SiMe3", "pph2": "PPh2", "pme2": "PMe2", "pipr2": "PiPr2",
    "f": "F", "cl": "Cl", "br": "Br", "i": "I", "vinyl": "vinyl", "allyl": "allyl",
}

#: No atom carries more than a handful of substituents. The caps are here
#: because this text comes from whoever is at the keyboard, and "Me99999999"
#: must be a clear refusal rather than a million atoms.
MOST_GROUPS = 12
LONGEST_SPEC = 200

_SEPARATOR = re.compile(r"[,;/+]|\s+")
#: Written with bond orders, brackets or branches that are not a count: SMILES,
#: not an abbreviation, so it is read as SMILES first (C(C)C is isopropyl, and
#: reading it letter by letter would make it three methyls).
_LOOKS_LIKE_SMILES = re.compile(r"[=#\[\]@]|\([^)]*\)(?!\d)")
_COUNT = re.compile(r"^(\d+)")


def names(spec: str, symbol: str | None = None) -> list[str]:
    """The groups a label stands for: "iPr2" -> two isopropyls.

    ``symbol`` is the element the label sits on ("P"), so the label as
    ChemDraw shows it ("PiPr2") can be passed straight in.
    """
    text = (spec or "").strip()
    if not text:
        return []
    if len(text) > LONGEST_SPEC:
        raise BackendError(f"that is a very long substituent: {LONGEST_SPEC} characters at most.")
    out: list[str] = []
    for token in _SEPARATOR.split(text):
        if token:
            out += _token_names(token, symbol)
        symbol = None  # only the first token can carry the element itself
    if len(out) > MOST_GROUPS:
        raise BackendError(
            f"{spec} means {len(out)} groups on one atom; {MOST_GROUPS} is already more "
            "than any atom carries."
        )
    return out


def complete(mol: Chem.Mol, plan: dict[int, list[str]]) -> Chem.Mol:
    """Hang the named groups on the given atoms, every other atom left alone.

    The coordinates of the drawing must survive untouched -- they are what the
    isomer is read from -- so new atoms are laid out around their own parent
    instead of by recomputing the whole depiction.
    """
    rw = Chem.RWMol(mol)
    conf = rw.GetConformer()
    for index, group_names in plan.items():
        for name in group_names:
            _attach(rw, conf, index, _group(name), _free_direction(rw, conf, index))
    out = rw.GetMol()
    out.UpdatePropertyCache(strict=False)
    return out


def chain(mol: Chem.Mol, index: int, group_names: list[str]) -> Chem.Mol:
    """Hang the groups one after another, as a written chain: CH2-OH."""
    rw = Chem.RWMol(mol)
    conf = rw.GetConformer()
    at = index
    for name in group_names:
        at = _attach(rw, conf, at, _group(name), _free_direction(rw, conf, at))
    out = rw.GetMol()
    out.UpdatePropertyCache(strict=False)
    return out


def replace(mol: Chem.Mol, index: int, group_names: list[str]) -> Chem.Mol:
    """Turn one atom into what its label says: a bare C labelled PhCH2OH.

    The atom keeps its place and its bonds -- it becomes the first atom of the
    first group -- and the rest of the label is chained onto it, which is how a
    label reads: the leftmost atom is the one the bond arrives at.
    """
    rw = Chem.RWMol(mol)
    conf = rw.GetConformer()
    at = _graft(rw, conf, index, _group(group_names[0]))
    for name in group_names[1:]:
        at = _attach(rw, conf, at, _group(name), _free_direction(rw, conf, at))
    out = rw.GetMol()
    out.UpdatePropertyCache(strict=False)
    return out


def pretty(group_names: list[str]) -> list[str]:
    """The groups as a chemist writes them: ["ipr", "ph"] -> ["iPr", "Ph"]."""
    shown = []
    for name in group_names:
        display = DISPLAY.get(name.lower())
        if display is None:  # a SMILES written by hand keeps its own spelling
            display = name.capitalize() if name.isalpha() else name
        shown.append(display)
    return shown


# -- reading a label -------------------------------------------------------------


def _token_names(token: str, symbol: str | None) -> list[str]:
    if _LOOKS_LIKE_SMILES.search(token) and Chem.MolFromSmiles(token) is not None:
        return [token]
    text = token.replace("(", "").replace(")", "").replace("-", "").replace(".", "")
    if symbol and len(text) > len(symbol) and text[: len(symbol)].lower() == symbol.lower():
        rest = text[len(symbol):]
        if rest and not rest[0].isdigit():
            text = rest
    out: list[str] = []
    while text:
        name, rest = _longest(text)
        if name is None:
            if Chem.MolFromSmiles(text) is not None:
                return out + [text]  # SMILES written by hand: kept as it is
            raise BackendError(
                f"m2i does not know the group {token}. Write it as Me, iPr, Ph, Cy, tBu, "
                "OMe..., as a combination such as CH2OH, or as SMILES."
            )
        count = _COUNT.match(rest)
        if count:
            if int(count.group(1)) > MOST_GROUPS:
                raise BackendError(
                    f"{token}: {count.group(1)} of the same group is more than any atom carries."
                )
            out += [name] * int(count.group(1))
            rest = rest[count.end():]
        else:
            out.append(name)
        text = rest
    return out


def _longest(text: str) -> tuple[str | None, str]:
    lowered = text.lower()
    for length in range(min(len(lowered), 16), 0, -1):
        piece = lowered[:length]
        if piece in ALIASES:
            return ALIASES[piece], text[length:]
        if piece in GROUPS:
            return piece, text[length:]
    return None, text


def _group(name: str) -> Chem.Mol:
    return Chem.Mol(_parsed(name))  # a copy: the caller lays coordinates on it


@lru_cache(maxsize=256)
def _parsed(name: str) -> Chem.Mol:
    key = name.lower()
    smiles = GROUPS.get(key) or GROUPS.get(ALIASES.get(key, ""), name)
    group = Chem.MolFromSmiles(smiles, sanitize=False)
    if group is None:
        raise BackendError(f"{name} is not a group m2i knows, nor a SMILES it can read.")
    group.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(group)
    return group


# -- putting the atoms on the page -------------------------------------------------


def _free_direction(rw: Chem.RWMol, conf, index: int) -> np.ndarray:
    """A 2D direction pointing away from everything already bonded here."""
    here = conf.GetAtomPosition(index)
    taken = []
    for neighbour in rw.GetAtomWithIdx(index).GetNeighbors():
        there = conf.GetAtomPosition(neighbour.GetIdx())
        v = np.array([there.x - here.x, there.y - here.y])
        if np.linalg.norm(v) > 1e-6:
            taken.append(v / np.linalg.norm(v))
    best, widest = np.array([1.0, 0.0]), -1.0
    for degrees in range(0, 360, 5):
        angle = np.radians(degrees)
        candidate = np.array([np.cos(angle), np.sin(angle)])
        gap = min((1.0 - float(candidate @ t) for t in taken), default=2.0)
        if gap > widest:
            best, widest = candidate, gap
    return best


def _graft(rw: Chem.RWMol, conf, index: int, group: Chem.Mol) -> int:
    """Make an existing atom the first atom of a group, adding the rest around it."""
    first = group.GetAtomWithIdx(0)
    atom = rw.GetAtomWithIdx(index)
    atom.SetAtomicNum(first.GetAtomicNum())
    atom.SetFormalCharge(first.GetFormalCharge())
    atom.SetIsAromatic(first.GetIsAromatic())
    atom.SetNoImplicit(False)
    atom.SetNumExplicitHs(0)
    rest = Chem.RWMol(group)
    rest.RemoveAtom(0)
    if rest.GetNumAtoms():
        # The group without its first atom, hung back on the atom it belongs to.
        mapping = _add_atoms(rw, conf, index, group)
        for bond in group.GetBonds():
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            rw.AddBond(mapping.get(a, index), mapping.get(b, index), bond.GetBondType())
    return index


def _add_atoms(rw: Chem.RWMol, conf, index: int, group: Chem.Mol) -> dict[int, int]:
    """Add every atom of the group but its first, laid out around ``index``."""
    rdDepictor.Compute2DCoords(group)
    local = group.GetConformer()
    origin = np.array([local.GetAtomPosition(0).x, local.GetAtomPosition(0).y])
    direction = _free_direction(rw, conf, index)
    angle = np.arctan2(direction[1], direction[0])
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    here = conf.GetAtomPosition(index)
    base = np.array([here.x, here.y])
    mapping = {}
    for atom in group.GetAtoms():
        if atom.GetIdx() == 0:
            continue
        mapping[atom.GetIdx()] = _copy_atom(rw, atom)
        p = local.GetAtomPosition(atom.GetIdx())
        xy = base + rotation @ (np.array([p.x, p.y]) - origin)
        conf.SetAtomPosition(mapping[atom.GetIdx()], Point3D(float(xy[0]), float(xy[1]), 0.0))
    return mapping


def _copy_atom(rw: Chem.RWMol, atom) -> int:
    new = Chem.Atom(atom.GetAtomicNum())
    new.SetFormalCharge(atom.GetFormalCharge())
    new.SetNoImplicit(atom.GetNoImplicit())
    new.SetNumExplicitHs(atom.GetNumExplicitHs())
    new.SetIsAromatic(atom.GetIsAromatic())
    return rw.AddAtom(new)


def _attach(rw: Chem.RWMol, conf, index: int, group: Chem.Mol, direction: np.ndarray) -> int:
    """Add the group, bond its first atom here, lay it out; return that atom."""
    rdDepictor.Compute2DCoords(group)
    local = group.GetConformer()
    origin = np.array([local.GetAtomPosition(0).x, local.GetAtomPosition(0).y])
    angle = np.arctan2(direction[1], direction[0])
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    here = conf.GetAtomPosition(index)
    base = np.array([here.x, here.y]) + 1.5 * direction

    mapping = {}
    for atom in group.GetAtoms():
        mapping[atom.GetIdx()] = _copy_atom(rw, atom)
        p = local.GetAtomPosition(atom.GetIdx())
        xy = base + rotation @ (np.array([p.x, p.y]) - origin)
        conf.SetAtomPosition(mapping[atom.GetIdx()], Point3D(float(xy[0]), float(xy[1]), 0.0))
    for bond in group.GetBonds():
        rw.AddBond(mapping[bond.GetBeginAtomIdx()], mapping[bond.GetEndAtomIdx()], bond.GetBondType())
        if bond.GetIsAromatic():
            rw.GetBondBetweenAtoms(
                mapping[bond.GetBeginAtomIdx()], mapping[bond.GetEndAtomIdx()]
            ).SetIsAromatic(True)
    rw.AddBond(index, mapping[0], Chem.BondType.SINGLE)
    return mapping[0]
