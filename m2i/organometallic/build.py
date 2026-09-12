"""A metal complex in 3D, with the coordination arrangement that was chosen.

RDKit's distance geometry cannot build a metal complex on its own: the bounds
it derives around a metal are inconsistent (triangle smoothing fails), and it
has no idea of an octahedron. So the bounds are assembled here instead:

- each ligand's own bounds, from RDKit, which handles ligands well;
- the coordination sphere: metal-donor distances, and donor-donor distances
  for the angles of the chosen polyhedron -- as ranges, because a chelate's
  bite angle is nowhere near an ideal 90 degrees;
- the metal to each donor's neighbours, from the donor's hybridisation, so a
  carbonyl points straight out and a phosphine's substituents splay.

Distance geometry cannot tell an arrangement from its mirror image, so the
handedness is checked afterwards and put right. A UFF pre-optimisation then
tidies the result, with the metal-ligand bonds made ordinary single bonds
(UFF ignores dative ones, and would let the metal and its donors only repel).

The result is a starting geometry for a quantum-chemical optimisation, not a
final structure. Tested against DFT geometries of pincer complexes: the core
comes within about 0.2 Å.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem
from rdkit import DistanceGeometry as DG
from rdkit.Chem import AllChem, rdDistGeom

from ..crystal.structure import Species
from ..recognition.base import BackendError
from .drawing import OrganometallicDrawing
from .geometry import TRANS, Arrangement, is_chiral

#: Angular ranges (degrees) allowed for cis and trans donor pairs: wide for
#: two donors of one chelate, whose bite angle is nowhere near 90 degrees...
CIS_RANGE = (72.0, 108.0)
TRANS_RANGE = (158.0, 180.0)
#: ...narrower when one of them is in a chelate, and narrow between two
#: monodentate ligands. The narrow ones matter in the UFF step: RDKit's UFF
#: bends every angle at a metal towards the tetrahedral (its energy is highest
#: at 180 degrees), so the angles end up at the edge of whatever is allowed.
CIS_CHELATE = (80.0, 100.0)
TRANS_CHELATE = (165.0, 180.0)
CIS_MONO = (85.0, 95.0)
TRANS_MONO = (172.0, 180.0)
#: Metal-donor distances start from covalent radii, shortened a little.
RADIUS_SHIFT = -0.15
ATTEMPTS = 25
#: The hybridisation UFF expects of a metal in each polyhedron.
METAL_HYBRIDISATION = {
    "tetrahedral": Chem.HybridizationType.SP3,
    "square planar": Chem.HybridizationType.SP2D,
    "trigonal bipyramidal": Chem.HybridizationType.SP3D,
    "square pyramidal": Chem.HybridizationType.SP3D,
    "octahedral": Chem.HybridizationType.SP3D2,
}


@dataclass
class BuiltComplex:
    mol: Chem.Mol  # with hydrogens and a 3D conformer; metal bonds dative
    species: Species
    arrangement: Arrangement
    notes: list[tuple[str, str, str]] = field(default_factory=list)


def build(drawing: OrganometallicDrawing, arrangement: Arrangement, *, seed: int = 0xF00D) -> BuiltComplex:
    mol = Chem.AddHs(drawing.mol)
    metal = drawing.metal
    notes: list[tuple[str, str, str]] = []

    bounds, radius = _bounds(mol, metal, arrangement)
    chiral = is_chiral(arrangement, drawing.labels, drawing.chelated)
    reference = _reference_triple(arrangement)
    has_stereo = any(
        a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms() if a.GetIdx() != metal
    )

    conformer = None
    for attempt in range(ATTEMPTS):
        params = AllChem.ETKDGv3()
        params.useRandomCoords = True
        params.randomSeed = seed + attempt
        params.enforceChirality = True
        params.SetBoundsMat(bounds)
        if AllChem.EmbedMolecule(mol, params) < 0:
            continue
        positions = mol.GetConformer().GetPositions()
        if _trans_built(positions, metal, arrangement) != arrangement.trans_pairs():
            continue  # the bounds were only met approximately: another try
        if chiral and not _same_hand(positions, metal, arrangement, reference):
            if has_stereo:
                continue  # reflecting would invert the ligands' own stereocentres
            positions[:, 0] *= -1
            _set_positions(mol, positions)
        conformer = positions
        break
    if conformer is None:
        raise BackendError(
            f"could not build the {arrangement.polyhedron} arrangement in 3D; the "
            "chelates may not fit it (a rigid pincer cannot be fac, for instance)"
        )

    _preoptimise(mol, metal, arrangement, radius, drawing.chelated, notes)
    species = _species(mol)
    _check(species, metal, arrangement, reference, chiral, notes)
    return BuiltComplex(mol=mol, species=species, arrangement=arrangement, notes=notes)


# -- bounds -----------------------------------------------------------------------


def _bounds(mol: Chem.Mol, metal: int, arrangement: Arrangement):
    table = Chem.GetPeriodicTable()
    n = mol.GetNumAtoms()
    bm = np.zeros((n, n))
    for a in range(n):
        za = mol.GetAtomWithIdx(a).GetAtomicNum()
        for b in range(a + 1, n):
            zb = mol.GetAtomWithIdx(b).GetAtomicNum()
            bm[a][b] = 100.0
            bm[b][a] = 0.7 * (table.GetRvdw(za) + table.GetRvdw(zb))

    # each ligand's own bounds
    rw = Chem.RWMol(mol)
    for d in arrangement.sites:
        rw.RemoveBond(d, metal)
    detached = rw.GetMol()
    detached.UpdatePropertyCache(strict=False)
    mapping: list = []
    fragments = Chem.GetMolFrags(detached, asMols=True, sanitizeFrags=True, fragsMolAtomMapping=mapping)
    for fragment, atoms in zip(fragments, mapping):
        if atoms == (metal,):
            continue
        local = rdDistGeom.GetMoleculeBoundsMatrix(fragment)
        for x in range(len(atoms)):
            for y in range(x + 1, len(atoms)):
                a, b = sorted((atoms[x], atoms[y]))
                upper, lower = (local[x][y], local[y][x]) if atoms[x] < atoms[y] else (local[y][x], local[x][y])
                bm[a][b], bm[b][a] = max(upper, lower), min(upper, lower)

    def put(a, b, low, high):
        a, b = sorted((a, b))
        bm[a][b], bm[b][a] = high, low

    def chord(r1, r2, degrees):
        return float(np.sqrt(r1 * r1 + r2 * r2 - 2 * r1 * r2 * np.cos(np.radians(degrees))))

    zm = mol.GetAtomWithIdx(metal).GetAtomicNum()
    radius = {}
    for d in arrangement.sites:
        radius[d] = table.GetRcovalent(zm) + table.GetRcovalent(mol.GetAtomWithIdx(d).GetAtomicNum()) + RADIUS_SHIFT
        put(metal, d, radius[d] - 0.08, radius[d] + 0.08)

    donors = list(arrangement.sites)
    for i, a in enumerate(donors):
        for b in donors[i + 1:]:
            low, high = TRANS_RANGE if _is_trans(arrangement, a, b) else CIS_RANGE
            put(a, b, chord(radius[a], radius[b], low), chord(radius[a], radius[b], high))

    for d in donors:
        atom = mol.GetAtomWithIdx(d)
        hybrid = atom.GetHybridization()
        if hybrid == Chem.HybridizationType.SP:
            low, high = 165.0, 180.0
        elif hybrid in (Chem.HybridizationType.SP2, Chem.HybridizationType.SP2D):
            low, high = 108.0, 135.0
        else:
            low, high = 100.0, 128.0
        for neighbour in atom.GetNeighbors():
            x = neighbour.GetIdx()
            if x == metal:
                continue
            a, b = sorted((d, x))
            bond = 0.5 * (bm[a][b] + bm[b][a])
            put(metal, x, chord(radius[d], bond, low), chord(radius[d], bond, high))

    if not DG.DoTriangleSmoothing(bm):
        raise BackendError(
            "the chosen arrangement is geometrically impossible for these ligands "
            "(their own shape contradicts it)"
        )
    return bm, radius


def _is_trans(arrangement: Arrangement, a: int, b: int) -> bool:
    return frozenset((a, b)) in arrangement.trans_pairs()


# -- handedness -------------------------------------------------------------------


def _reference_triple(arrangement: Arrangement):
    donors = sorted(arrangement.sites)
    for i in range(len(donors)):
        for j in range(i + 1, len(donors)):
            for k in range(j + 1, len(donors)):
                trio = (donors[i], donors[j], donors[k])
                volume = np.linalg.det(np.array([arrangement.sites[t] for t in trio]))
                if abs(volume) > 0.1:
                    return trio, float(np.sign(volume))
    return None, 0.0


def _same_hand(positions, metal, arrangement, reference) -> bool:
    trio, sign = reference
    if trio is None:
        return True
    vectors = np.array([positions[t] - positions[metal] for t in trio])
    return np.sign(np.linalg.det(vectors)) == sign


def _set_positions(mol: Chem.Mol, positions) -> None:
    conf = mol.GetConformer()
    for i, p in enumerate(positions):
        conf.SetAtomPosition(i, p.tolist())


# -- pre-optimisation -------------------------------------------------------------


def _preoptimise(mol, metal, arrangement, radius, chelated, notes) -> None:
    kekulized = Chem.Mol(mol)
    Chem.Kekulize(kekulized, clearAromaticFlags=True)
    rw = Chem.RWMol(kekulized)
    for bond in rw.GetBonds():
        if bond.GetBondType() == Chem.BondType.DATIVE:
            bond.SetBondType(Chem.BondType.SINGLE)
    ffmol = rw.GetMol()
    ffmol.UpdatePropertyCache(strict=False)
    Chem.SetHybridization(ffmol)
    # UFF types a metal by its hybridisation; RDKit's guess from the bond
    # count calls a square-planar Pt sp3, for which UFF has no parameters.
    hybrid = METAL_HYBRIDISATION.get(arrangement.polyhedron)
    if hybrid is not None:
        ffmol.GetAtomWithIdx(metal).SetHybridization(hybrid)
    Chem.GetSymmSSSR(ffmol)
    if not AllChem.UFFHasAllMoleculeParams(ffmol):
        notes.append(("warning", "organometallic.no_uff",
                      "UFF lacks parameters for part of this complex: the geometry is the "
                      "distance-geometry one, not pre-optimised."))
        return
    field_ = AllChem.UFFGetMoleculeForceField(ffmol)
    donors = list(arrangement.sites)
    for d in donors:
        field_.UFFAddDistanceConstraint(metal, d, False, radius[d] - 0.35, radius[d] + 0.15, 100.0)
    pairs = {frozenset(pair) for pair in chelated}
    in_chelate = {d for pair in chelated for d in pair}
    for i, a in enumerate(donors):
        for b in donors[i + 1:]:
            trans = _is_trans(arrangement, a, b)
            if frozenset((a, b)) in pairs:
                low, high = TRANS_RANGE if trans else CIS_RANGE
            elif a in in_chelate or b in in_chelate:
                low, high = TRANS_CHELATE if trans else CIS_CHELATE
            else:
                low, high = TRANS_MONO if trans else CIS_MONO
            field_.UFFAddAngleConstraint(a, metal, b, False, low, high, 100.0)
    field_.Minimize(maxIts=5000)
    _set_positions(mol, ffmol.GetConformer().GetPositions())


# -- checks and conversion ----------------------------------------------------------


def _species(mol: Chem.Mol) -> Species:
    positions = mol.GetConformer().GetPositions()
    elements = [a.GetSymbol() for a in mol.GetAtoms()]
    # Numbered as in the drawing (hydrogens follow), so the names printed for
    # the built complex are the ones the isomer list used.
    labels = [f"{e}{i + 1}" for i, e in enumerate(elements)]
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    return Species(elements=elements, coords=np.array(positions), labels=labels, bonds=bonds)


def _trans_built(positions, metal, arrangement: Arrangement) -> set:
    donors = list(arrangement.sites)
    found = set()
    for i, a in enumerate(donors):
        for b in donors[i + 1:]:
            u, v = positions[a] - positions[metal], positions[b] - positions[metal]
            angle = np.degrees(np.arccos(np.clip(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v), -1, 1)))
            if angle >= TRANS:
                found.add(frozenset((a, b)))
    return found


def _check(species, metal, arrangement, reference, chiral, notes) -> None:
    """The built geometry must show the arrangement that was asked for."""
    positions = species.coords
    if _trans_built(positions, metal, arrangement) != arrangement.trans_pairs():
        raise BackendError("the 3D structure lost the chosen arrangement during pre-optimisation")
    if chiral and not _same_hand(positions, metal, arrangement, reference):
        raise BackendError("the 3D structure came out as the mirror image of the chosen arrangement")
