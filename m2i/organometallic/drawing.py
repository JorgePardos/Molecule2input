"""Reading a drawn metal complex: connectivity, hydrogens, and what the page hints.

ChemDraw, and molfiles, draw the bonds to a metal as ordinary lines (sometimes
as dative arrows). Read naively, a phosphine bound to Mn has four bonds and a
carbonyl carbon looks like a formyl group. So, before anything else:

- a carbon bound to the metal whose only other neighbour is a terminal O or N
  is a carbonyl, cyanide or isocyanide: it gets no hydrogen;
- every other donor with a label (P, N, O, S...) binds through a lone pair:
  its bond to the metal becomes dative and its hydrogens are those drawn;
- an unlabelled carbon bound to the metal is an alkyl or aryl ligand: that
  bond counts towards its valence, so a bare line ending at the metal is a
  methyl, as it would be on paper;
- a halogen or hydrogen on its own, and a labelled N, O or S one bond short
  (amido, alkoxide, thiolate), are X ligands: their bond to the metal is a
  covalent one.

What cannot be settled -- a phosphine with one substituent, a label ChemDraw
could not interpret -- is reported rather than guessed: those are signs of an
incomplete drawing, and the resulting molecule would not be the one meant.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from rdkit import Chem

from ..crystal.structure import METALS
from ..recognition.base import BackendError
from . import chemdraw, groups
from .geometry import drawn_directions

#: A labelled atom with fewer bonds (hydrogens included) than this, and no
#: negative charge, is short of substituents -- an incomplete drawing.
#: Carbon is left out: a carbonyl carbon is legitimately short once read as C#O.
USUAL_VALENCE = {"P": 3, "N": 3, "O": 2, "S": 2, "As": 3, "Se": 2, "B": 3}
HALOGENS = {"F", "Cl", "Br", "I"}
#: Donors that, drawn one bond short, are read as anions bound covalently.
#: A phosphine one bond short is far more often an unfinished drawing.
ANIONIC = {"N", "O", "S"}
INCOMPLETE_HINT = {
    "P": "A phosphine needs its three substituents: draw them, or write PMe2, PPh2, PCy2...",
    "N": "Count its hydrogens (NH, NH2) and any ring it belongs to.",
    "As": "An arsine needs its three substituents.",
}


@dataclass
class OrganometallicDrawing:
    mol: Chem.Mol
    metal: int
    donors: list[int]
    #: donor -> 3D direction guessed from the page
    drawn: dict[int, np.ndarray]
    #: donor pairs linked by a chelate ring small enough to force them cis
    chelated: list[tuple[int, int]]
    #: donor -> a name that is equal for chemically equivalent donors
    labels: dict[int, str]
    #: donor -> readable description, e.g. "P of C9H14N3P2"
    names: dict[int, str]
    notes: list[tuple[str, str, str]] = field(default_factory=list)
    #: sum of the formal charges as drawn: a starting point for the charge
    drawn_charge: int = 0
    #: atom name -> groups added to it, because the drawing left them out
    completed: dict[str, list[str]] = field(default_factory=dict)
    #: atom name -> hydrogens nobody drew, added by the usual convention
    assumed: dict[str, list[str]] = field(default_factory=dict)

    @property
    def metal_symbol(self) -> str:
        return self.mol.GetAtomWithIdx(self.metal).GetSymbol()


def has_metal(mol: Chem.Mol) -> bool:
    return any(a.GetSymbol() in METALS for a in mol.GetAtoms())


def read_raw(path: Path) -> list[Chem.Mol]:
    """Every molecule on the page, unsanitised, with 2D coordinates and wedges."""
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in (".cdx", ".cdxml"):
        params = Chem.CDXMLParserParams()
        params.format = Chem.CDXMLFormat.CDX if suffix == ".cdx" else Chem.CDXMLFormat.CDXML
        params.sanitize = False
        params.removeHs = False
        mols = Chem.MolsFromCDXML(path.read_bytes(), params)
    elif suffix in (".mol", ".sdf", ".mdl"):
        text = path.read_text(encoding="utf-8", errors="replace").split("$$$$")[0]
        mols = [Chem.MolFromMolBlock(text, sanitize=False, removeHs=False)]
    else:
        raise BackendError(f"{path.name}: not a drawing file m2i can read")
    return [m for m in mols if m is not None and m.GetNumAtoms()]


def read(path: Path, substituents: dict[str, str] | None = None) -> OrganometallicDrawing:
    """The metal complex drawn in ``path``, prepared for building in 3D.

    ``substituents`` fills in what the drawing leaves out, keyed by the atom
    names m2i prints: {"P6": "iPr2"}.
    """
    mols = read_raw(path)
    with_metal = [m for m in mols if has_metal(m)]
    if not with_metal:
        raise BackendError(f"{Path(path).name}: no metal in the drawing")
    notes: list[tuple[str, str, str]] = []
    mol = with_metal[0]
    mol = _read_labels(path, mol, notes)
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if len(mols) > 1 or len(fragments) > 1:
        # Counter-ions and solvent drawn beside the complex are not part of it.
        pieces = [f for f in fragments if has_metal(f)]
        others = len(mols) - 1 + len(fragments) - len(pieces)
        mol = pieces[0]
        if others:
            notes.append(("warning", "organometallic.other_fragments",
                          f"{others} other fragment(s) on the page were left out: only the "
                          "species containing the metal is built."))
    return prepare(mol, notes, substituents)


def _read_labels(path: Path, mol: Chem.Mol, notes: list) -> Chem.Mol:
    """Expand the labels ChemDraw kept as plain text and RDKit turned into carbon.

    A label ChemDraw did not interpret is still exact text, so it is expanded
    here. It is read from its first atom, which is where the bond arrives; when
    the metal is meant to bind through another atom of the label, only drawing
    that atom can say so, and the warning says as much.
    """
    for index, text in sorted(chemdraw.labels(path, mol).items()):
        try:
            group_names = groups.names(text)
        except BackendError as exc:
            notes.append(("warning", "organometallic.unreadable_label",
                          f"ChemDraw could not interpret the label {text}, and neither can "
                          f"m2i: {exc} It is left as the bare atom RDKit read."))
            continue
        mol = groups.replace(mol, index, group_names)
        notes.append(("warning", "organometallic.uninterpreted_label",
                      f"ChemDraw did not interpret the label {text} (it is plain text in the "
                      f"file, not a structure). m2i read it as {chr(45).join(groups.pretty(group_names))}, "
                      "bound to the metal through its first atom. If the metal binds through "
                      "another atom, draw that atom and its bond."))
    return mol


def prepare(mol: Chem.Mol, notes: list | None = None,
            substituents: dict[str, str] | None = None) -> OrganometallicDrawing:
    notes = list(notes or [])
    metals = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in METALS]
    if len(metals) > 1:
        raise BackendError(
            "complexes with more than one metal are not supported from drawings yet; "
            "start from a crystal structure (.cif) instead"
        )
    metal = metals[0]
    rw = Chem.RWMol(mol)
    rw.UpdatePropertyCache(strict=False)
    donors = sorted(n.GetIdx() for n in rw.GetAtomWithIdx(metal).GetNeighbors())
    if not donors:
        raise BackendError("the metal has no bonds drawn to it")

    haptic = [d for d in donors if any(n.GetIdx() in donors for n in rw.GetAtomWithIdx(d).GetNeighbors())]
    if haptic:
        raise BackendError(
            "ligands bound through several adjacent atoms (eta-2 alkenes, Cp, arenes) are "
            "not supported from drawings yet; start from a crystal structure (.cif)"
        )

    depth = _wedge_depths(rw, metal)
    for d in donors:
        _set_donor_bond(rw, metal, d, notes)

    fixed = rw.GetMol()
    fixed.UpdatePropertyCache(strict=False)
    fixed, completed = _fill(fixed, metal, substituents or {}, notes)
    try:
        Chem.SanitizeMol(fixed)
    except Exception as exc:  # whatever RDKit raises, reported with its reason
        raise BackendError(f"the drawn complex cannot be interpreted: {exc}") from exc

    conf = fixed.GetConformer() if fixed.GetNumConformers() else None
    if conf is None:
        raise BackendError("the drawing has no 2D coordinates")
    at = conf.GetAtomPosition
    drawn = drawn_directions(
        (at(metal).x, at(metal).y), {d: (at(d).x, at(d).y) for d in donors}, depth
    )
    implicit, assumed = _implicit_hydrogens(fixed, donors)
    notes += implicit
    labels, names, chelated = _ligand_facts(fixed, metal, donors)
    charge = sum(a.GetFormalCharge() for a in fixed.GetAtoms())
    return OrganometallicDrawing(
        mol=fixed, metal=metal, donors=donors, drawn=drawn, chelated=chelated,
        labels=labels, names=names, notes=notes, drawn_charge=charge, completed=completed,
        assumed=assumed,
    )


# -- the bonds to the metal ----------------------------------------------------


def _set_donor_bond(rw: Chem.RWMol, metal: int, donor: int, notes) -> None:
    atom = rw.GetAtomWithIdx(donor)
    bond = rw.GetBondBetweenAtoms(metal, donor)
    others = [n for n in atom.GetNeighbors() if n.GetIdx() != metal]

    if atom.GetSymbol() == "C" and len(others) == 1 and others[0].GetDegree() <= 2 \
            and others[0].GetSymbol() in ("O", "N"):
        # Carbonyl, cyanide or isocyanide: C triply bonded, no hydrogen.
        partner = others[0]
        rw.GetBondBetweenAtoms(donor, partner.GetIdx()).SetBondType(Chem.BondType.TRIPLE)
        atom.SetFormalCharge(-1)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(0)
        if partner.GetSymbol() == "O" or partner.GetDegree() == 2:
            partner.SetFormalCharge(1)  # C#O and C#N-R are neutral; cyanide keeps its -1
        _make_dative(rw, metal, donor, bond)
        return

    if atom.GetSymbol() == "C" and not atom.GetNoImplicit():
        # An unlabelled carbon: an alkyl or aryl ligand, its bond to the metal
        # counted like any other.
        bond.SetBondType(Chem.BondType.SINGLE)
        return

    if atom.GetFormalCharge() == 0 and (
        atom.GetSymbol() == "H" or (atom.GetSymbol() in HALOGENS and not others)
    ):
        # Hydride or halide: an X ligand whose one bond is the one to the metal.
        bond.SetBondType(Chem.BondType.SINGLE)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(0)
        return

    if atom.GetNoImplicit() and atom.GetFormalCharge() == 0 and atom.GetSymbol() in ANIONIC:
        # A labelled N, O or S drawn one bond short -- "N" in a pincer backbone,
        # the "O" of an acac or a carboxylate -- is an amido, alkoxide or
        # thiolate: its bond to the metal is a covalent one, as it was drawn.
        used = _bonds_besides(atom, metal) + atom.GetNumExplicitHs()
        if used == USUAL_VALENCE[atom.GetSymbol()] - 1:
            bond.SetBondType(Chem.BondType.SINGLE)
            kind = {"N": "an amido (N-)", "O": "an alkoxide/carboxylate (O-)", "S": "a thiolate (S-)"}
            notes.append(("info", "organometallic.anionic_donor",
                          f"{atom.GetSymbol()}{donor + 1} is drawn one bond short: read as "
                          f"{kind[atom.GetSymbol()]} donor, bound covalently. If it carries a "
                          "hydrogen, draw it (NH, OH)."))
            return

    _make_dative(rw, metal, donor, bond)


def _bonds_besides(atom, metal: int) -> float:
    return sum(
        b.GetBondTypeAsDouble() for b in atom.GetBonds()
        if b.GetOtherAtomIdx(atom.GetIdx()) != metal
    )


def _make_dative(rw: Chem.RWMol, metal: int, donor: int, bond) -> None:
    """A bond from the donor's lone pair: it does not count towards its valence."""
    rw.RemoveBond(metal, donor)
    rw.AddBond(donor, metal, Chem.BondType.DATIVE)


def _wedge_depths(mol: Chem.Mol, metal: int) -> dict[int, float]:
    """+1 for a donor drawn in front of the page, -1 behind, 0 in its plane."""
    depth = {}
    for bond in mol.GetAtomWithIdx(metal).GetBonds():
        other = bond.GetOtherAtomIdx(metal)
        direction = bond.GetBondDir()
        # An unsanitised molfile keeps its wedges only as the raw stereo field.
        stereo = bond.GetIntProp("_MolFileBondStereo") if bond.HasProp("_MolFileBondStereo") else 0
        sign = 0.0
        if direction == Chem.BondDir.BEGINWEDGE or stereo == 1:
            sign = 1.0
        elif direction == Chem.BondDir.BEGINDASH or stereo == 6:
            sign = -1.0
        if bond.GetBeginAtomIdx() != metal:
            sign = -sign  # the narrow end at the ligand: the metal is the one in front
        depth[other] = sign
    return depth


def _shortfalls(mol: Chem.Mol, metal: int, labelled: bool = True) -> list[tuple[int, str, int]]:
    """Atoms with fewer bonds drawn than usual: (atom, symbol, how many short).

    ``labelled`` picks the ones whose hydrogens the drawing fixed -- a real
    gap -- rather than the ones RDKit fills in by convention.
    """
    short = []
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if atom.GetIdx() == metal or symbol not in USUAL_VALENCE:
            continue
        if atom.GetNoImplicit() is not labelled:
            continue
        # A covalent bond to the metal (amido, alkoxide) counts; a dative one does not.
        used = sum(
            b.GetBondTypeAsDouble() for b in atom.GetBonds() if b.GetBondType() != Chem.BondType.DATIVE
        ) + atom.GetNumExplicitHs()
        expected = USUAL_VALENCE[symbol] + atom.GetFormalCharge() * (1 if symbol in ("N", "P", "O", "S") else 0)
        if used < expected:
            short.append((atom.GetIdx(), symbol, int(expected - used)))
    return short


def _fill(mol: Chem.Mol, metal: int, substituents: dict[str, str], notes: list):
    """Give the short atoms the substituents asked for, or hydrogens by default.

    Nothing here is a chemical guess: hydrogen is what the valence needs, and
    it is said out loud every time, because PH2 is not PiPr2 and the
    calculation would be of another molecule.
    """
    wanted = {key.strip().upper(): value for key, value in (substituents or {}).items()}
    short = _shortfalls(mol, metal)
    # Atoms whose hydrogens RDKit assumed (a molfile "P" is phosphine) are not
    # filled in, but they can still be given their real groups.
    assumed = [gap for gap in _shortfalls(mol, metal, labelled=False)
               if f"{gap[1]}{gap[0] + 1}".upper() in wanted]
    if not short and not assumed:
        return mol, {}
    mol = Chem.RWMol(mol)
    for index, _, _ in assumed:
        atom = mol.GetAtomWithIdx(index)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(0)
    mol = mol.GetMol()
    mol.UpdatePropertyCache(strict=False)

    plan, chains, filled = {}, {}, {}
    for index, symbol, missing in sorted(short + assumed):
        name = f"{symbol}{index + 1}"
        spec = wanted.pop(name.upper(), None)
        if spec is None:
            plan[index] = ["h"] * missing
            filled[name] = ["H"] * missing
            notes.append(("warning", "organometallic.incomplete",
                          f"{name} is {missing} bond(s) short of the usual "
                          f"{USUAL_VALENCE[symbol]}: filled with {missing} hydrogen(s). "
                          + INCOMPLETE_HINT.get(symbol, "")
                          + f" Give the real ones with --sub {name}=... "
                          + "(Me, iPr2, Ph2, Cy2, OMe...) or in the browser."))
            continue
        group_names, as_chain = _interpret(spec, symbol, missing, name)
        if as_chain:
            chains[index] = group_names
        else:
            plan[index] = group_names
        filled[name] = groups.pretty(group_names)
        notes.append(("info", "organometallic.completed",
                      f"{name}: {spec} as you asked, so "
                      + ", ".join(filled[name]) + " added."))
    if wanted:
        raise BackendError(
            "nothing is missing at " + ", ".join(sorted(wanted))
            + "; --sub only adds what a drawing leaves out."
        )
    out = groups.complete(mol, plan) if plan else mol
    for index, group_names in chains.items():
        out = groups.chain(out, index, group_names)
    return out, filled


def _interpret(spec: str, symbol: str, missing: int, name: str) -> tuple[list[str], bool]:
    """Read a substituent spec, with or without the element it sits on.

    "PiPr2" on a phosphorus means two isopropyls, but "Ph3" means three
    phenyls, not three hydrogens: the reading that fills the missing bonds is
    the one meant, and the plain one wins a tie. A spec that fills a single
    bond with several groups is a chain, the way CH2OH is written.
    """
    readings = []
    for element in (None, symbol):
        try:
            parsed = groups.names(spec, element)
        except BackendError:
            continue
        if parsed and parsed not in readings:
            readings.append(parsed)
    if not readings:
        groups.names(spec, symbol)  # raises, with its own explanation
    for parsed in readings:
        if len(parsed) == missing:
            return parsed, False
    if missing == 1:
        return readings[0], True
    raise BackendError(
        f"{name} is {missing} bond(s) short but {spec} gives {len(readings[0])} group(s); "
        "write one per missing bond, for instance iPr2."
    )


def _implicit_hydrogens(mol: Chem.Mol, donors: list[int]):
    """Say which hydrogens a donor got that nobody drew.

    A molfile "N" bound to a metal is ammonia by the usual convention, and a
    "P" is phosphine. That convention is kept -- but it is written down, since
    the same drawing is often meant as an amide or a substituted phosphine.
    """
    notes, assumed = [], {}
    for donor in donors:
        atom = mol.GetAtomWithIdx(donor)
        count = atom.GetTotalNumHs()
        if not count or atom.GetNoImplicit() or atom.GetSymbol() not in USUAL_VALENCE:
            continue
        assumed[f"{atom.GetSymbol()}{donor + 1}"] = ["H"] * count
        notes.append(("info", "organometallic.implicit_hydrogens",
                      f"{atom.GetSymbol()}{donor + 1} had no substituents drawn: read as "
                      f"{atom.GetSymbol()}H{count} by the usual convention. If it carries "
                      f"groups, give them with --sub {atom.GetSymbol()}{donor + 1}=iPr2."))
    return notes, assumed


def _ligand_facts(mol: Chem.Mol, metal: int, donors: list[int]):
    """Labels for equivalent donors, readable names, and which donors are chelated."""
    rw = Chem.RWMol(mol)
    for d in donors:
        rw.RemoveBond(d, metal)
    detached = rw.GetMol()
    detached.UpdatePropertyCache(strict=False)
    fragments = Chem.GetMolFrags(detached, asMols=False, sanitizeFrags=False)
    ranks = list(Chem.CanonicalRankAtoms(detached, breakTies=False))
    labels, names = {}, {}
    fragment_of = {}
    for number, fragment in enumerate(fragments):
        if fragment == (metal,):
            continue
        smiles = Chem.MolFragmentToSmiles(detached, atomsToUse=list(fragment), canonical=True)
        formula = _formula(detached, fragment)
        for atom in fragment:
            fragment_of[atom] = number
        for d in donors:
            if d in fragment:
                labels[d] = f"{smiles}|{ranks[d]}"
                names[d] = f"{mol.GetAtomWithIdx(d).GetSymbol()}{d + 1} of {formula}"
    chelated = []
    for i, a in enumerate(donors):
        for b in donors[i + 1:]:
            if fragment_of.get(a) != fragment_of.get(b):
                continue
            path = Chem.GetShortestPath(detached, a, b)
            if path and len(path) - 1 <= 5:  # a chelate ring of at most seven members
                chelated.append((a, b))
    return labels, names, chelated


def _formula(mol: Chem.Mol, atoms) -> str:
    from ..crystal.structure import hill_formula

    elements = []
    for index in atoms:
        atom = mol.GetAtomWithIdx(index)
        elements.append(atom.GetSymbol())
        elements += ["H"] * atom.GetTotalNumHs()
    return hill_formula(elements)
