"""The coordination sphere of every metal: what binds, how, and in what shape.

CIP descriptors and InChI -- the checks m2i uses for organic molecules -- are
blind at a metal centre: they say nothing about cis/trans, fac/mer or
hapticity. What is reported here is what a chemist looks at instead, so that
the geometry going into the input can be checked against the structure the
chemist expects.

Donors that are bonded to each other form one haptic ligand (the five carbons
of a Cp are one eta-5 donor) and are represented by their centroid, which is
what the angles and the polyhedron refer to.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

import numpy as np

from .chirality import planar_chirality
from .structure import METALS, Species, hill_formula

#: Two donors with an angle at the metal above this are trans.
TRANS_ANGLE = 150.0


@dataclass
class Donor:
    atoms: list[int]
    label: str
    #: What kind of donor this is: the donor element(s) plus the ligand it is
    #: part of, so a PPh3 phosphorus and a PMe3 phosphorus are told apart.
    kind: str
    position: np.ndarray
    distance: float
    #: Which ligand molecule this donor is part of (-1 for a metal).
    ligand: int = -1

    @property
    def hapticity(self) -> int:
        return len(self.atoms)


@dataclass
class MetalCentre:
    index: int
    label: str
    element: str
    donors: list[Donor] = field(default_factory=list)
    geometry: str = ""
    detail: str = ""
    isomers: list[str] = field(default_factory=list)
    trans_pairs: list[tuple[str, str, float]] = field(default_factory=list)
    #: "Delta" or "Lambda" for an octahedral centre with two or more chelates.
    helicity: str | None = None
    helicity_detail: str = ""
    #: Planar chirality of each eta-5/eta-6 ring that has it.
    planar: list = field(default_factory=list)

    @property
    def coordination_number(self) -> int:
        return len(self.donors)

    def describe(self) -> list[str]:
        lines = [
            f"{self.label} ({self.element}): {self.geometry}"
            + (f", {self.detail}" if self.detail else "")
        ]
        for donor in self.donors:
            eta = f"eta{donor.hapticity}-" if donor.hapticity > 1 else ""
            lines.append(f"    {eta}{donor.label:<22} {donor.distance:6.3f} A   [{donor.kind}]")
        for a, b, angle in self.trans_pairs:
            lines.append(f"    trans: {a} / {b}  ({angle:.1f} deg)")
        for isomer in self.isomers:
            lines.append(f"    isomer: {isomer}")
        if self.helicity:
            lines.append(f"    helicity: {self.helicity} ({self.helicity_detail})")
        for ring in self.planar:
            lines.append(f"    {ring.describe()}")
        return lines


def analyse(species: Species) -> list[MetalCentre]:
    """Every metal atom in the species, with its donors and polyhedron."""
    ligand_of = _ligands(species)
    centres = []
    for index, element in enumerate(species.elements):
        if element not in METALS:
            continue
        centre = MetalCentre(index=index, label=species.labels[index], element=element)
        centre.donors = _donors(species, index, ligand_of)
        _classify(species, centre)
        _helicity(species, centre)
        centre.planar = [
            found
            for donor in centre.donors
            if donor.hapticity >= 5
            and (found := planar_chirality(species, index, donor.atoms)) is not None
        ]
        centres.append(centre)
    return centres


# -- donors ----------------------------------------------------------------


def _ligands(species: Species) -> dict[int, tuple[int, str]]:
    """(ligand id, ligand formula) for each non-metal atom, metals removed."""
    adjacency: dict[int, set[int]] = {i: set() for i in range(species.n_atoms)}
    for a, b in species.bonds:
        if species.elements[a] in METALS or species.elements[b] in METALS:
            continue
        adjacency[a].add(b)
        adjacency[b].add(a)

    ligand_of: dict[int, tuple[int, str]] = {}
    for start in range(species.n_atoms):
        if start in ligand_of or species.elements[start] in METALS:
            continue
        members, stack = {start}, [start]
        while stack:
            for other in adjacency[stack.pop()] - members:
                members.add(other)
                stack.append(other)
        formula = hill_formula(species.elements[m] for m in members)
        for member in members:
            ligand_of[member] = (start, formula)
    return ligand_of


def _donors(species: Species, metal: int, ligand_of: dict) -> list[Donor]:
    bound = species.neighbours(metal)
    bound_set = set(bound)
    # Donor atoms bonded to each other are one haptic donor.
    groups, placed = [], set()
    for atom in bound:
        if atom in placed:
            continue
        group, stack = {atom}, [atom]
        while stack:
            here = stack.pop()
            for other in species.neighbours(here):
                if other in bound_set and other not in group and species.elements[other] not in METALS:
                    group.add(other)
                    stack.append(other)
        if species.elements[atom] in METALS:
            group = {atom}  # a metal-metal bond is its own donor
        placed |= group
        groups.append(sorted(group))

    donors = []
    for group in groups:
        position = species.coords[group].mean(axis=0)
        elements = [species.elements[g] for g in group]
        if len(group) == 1:
            label = species.labels[group[0]]
        else:
            label = f"{hill_formula(elements)} ({species.labels[group[0]]}..)"
        ligand_id, ligand = ligand_of.get(group[0], (-1, species.elements[group[0]]))
        kind = f"{hill_formula(elements) if len(group) > 1 else elements[0]} of {ligand}"
        donors.append(
            Donor(
                atoms=group,
                label=label,
                kind=kind,
                position=position,
                distance=float(np.linalg.norm(position - species.coords[metal])),
                ligand=ligand_id,
            )
        )
    donors.sort(key=lambda d: (d.kind, d.distance))
    return donors


# -- polyhedron and isomers -----------------------------------------------


def _classify(species: Species, centre: MetalCentre) -> None:
    metal = species.coords[centre.index]
    donors = centre.donors
    n = len(donors)
    vectors = [d.position - metal for d in donors]

    angles = {}
    for i in range(n):
        for j in range(i + 1, n):
            angles[(i, j)] = _angle(vectors[i], vectors[j])
    ordered = sorted(angles.values(), reverse=True)

    centre.trans_pairs = [
        (donors[i].label, donors[j].label, angle)
        for (i, j), angle in sorted(angles.items())
        if angle >= TRANS_ANGLE
    ]

    if n == 0:
        centre.geometry = "no donor atoms within bonding distance"
    elif n == 1:
        centre.geometry = "one donor"
    elif n == 2:
        haptic = all(d.hapticity >= 5 for d in donors)
        angle = ordered[0]
        if haptic:
            centre.geometry = "sandwich" if angle >= 170 else "bent sandwich"
            centre.detail = f"centroid-M-centroid {angle:.1f} deg"
        else:
            centre.geometry = "linear" if angle >= 160 else "bent"
            centre.detail = f"{angle:.1f} deg"
    elif n == 3:
        total = sum(ordered)
        if total >= 350:
            centre.geometry = "T-shaped" if ordered[0] >= 150 else "trigonal planar"
        else:
            centre.geometry = "trigonal pyramidal"
        centre.detail = f"sum of angles {total:.1f} deg"
    elif n == 4:
        tau4 = (360.0 - (ordered[0] + ordered[1])) / 141.0
        centre.geometry = (
            "square planar" if tau4 < 0.25 else "tetrahedral" if tau4 > 0.75 else "seesaw"
        )
        centre.detail = f"tau4 = {tau4:.2f}"
    elif n == 5:
        tau5 = (ordered[0] - ordered[1]) / 60.0
        centre.geometry = (
            "square pyramidal" if tau5 < 0.3 else "trigonal bipyramidal" if tau5 > 0.7
            else "between square pyramid and trigonal bipyramid"
        )
        centre.detail = f"tau5 = {tau5:.2f}"
    elif n == 6:
        trans = len(centre.trans_pairs)
        centre.geometry = (
            "octahedral" if trans == 3 else "trigonal prismatic" if trans == 0
            else "distorted octahedral"
        )
        centre.detail = f"{trans} trans pair(s)"
    else:
        centre.geometry = f"coordination number {n}"

    centre.isomers = _isomers(centre, angles)


def _isomers(centre: MetalCentre, angles: dict) -> list[str]:
    """cis/trans for pairs of identical donors, fac/mer for triples.

    Only meaningful on square-planar and octahedral centres, where 'trans'
    has a definite meaning.
    """
    if centre.geometry not in ("square planar", "octahedral", "distorted octahedral"):
        return []
    donors = centre.donors
    kinds = Counter(d.kind for d in donors)
    found = []
    for kind, count in sorted(kinds.items()):
        if count == len(donors):
            continue
        members = [i for i, d in enumerate(donors) if d.kind == kind]
        pair_angles = [angles[tuple(sorted((a, b)))] for a in members for b in members if a < b]
        trans = sum(angle >= TRANS_ANGLE for angle in pair_angles)
        if count == 2:
            found.append(f"{'trans' if trans else 'cis'} ({kind} x2)")
        elif count == 3 and centre.coordination_number == 6:
            found.append(f"{'mer' if trans else 'fac'} ({kind} x3)")
    return found


#: Two chelate edges closer to parallel than this define no helicity.
PARALLEL_EDGES = 15.0


def _helicity(species: Species, centre: MetalCentre) -> None:
    """Delta or Lambda, from the skew-line convention (IUPAC IR-9.3.4.11).

    Each chelate spans an edge of the octahedron. Two edges belonging to
    different chelates and sharing no donor are skew lines, and the sense in
    which one turns onto the other about their common perpendicular is the
    helicity. A tris-chelate gives three such pairs, a cis-bis-chelate one,
    and a proper Delta or Lambda centre has them all agree.

    The sign was fixed empirically, not from memory: against five COD
    structures of resolved complexes whose authors assign the label (Lambda:
    2006181 [Ru(phen)3]2+, 2010630 [Co(en)3]3+; Delta: 2016683 [Ru(phen)3]2+,
    2006867 and 1572653), with ``skew_angle`` as defined below Delta comes out
    negative in every one. Getting this backwards would mislabel every complex.
    """
    if centre.geometry not in ("octahedral", "distorted octahedral"):
        return
    metal = species.coords[centre.index]
    edges = []
    for ligand in {d.ligand for d in centre.donors if d.ligand >= 0}:
        members = [d for d in centre.donors if d.ligand == ligand and d.hapticity == 1]
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                a, b = members[i], members[j]
                if _angle(a.position - metal, b.position - metal) < TRANS_ANGLE:
                    edges.append((ligand, a, b))
    signs = []
    for i in range(len(edges)):
        for j in range(i + 1, len(edges)):
            (la, a1, a2), (lb, b1, b2) = edges[i], edges[j]
            if la == lb or {id(a1), id(a2)} & {id(b1), id(b2)}:
                continue
            angle = skew_angle(a1.position, a2.position, b1.position, b2.position)
            if angle is not None and abs(angle) >= PARALLEL_EDGES:
                signs.append(angle)
    if not signs:
        return
    delta = sum(a < 0 for a in signs)
    if delta == len(signs):
        centre.helicity = "Delta"
    elif delta == 0:
        centre.helicity = "Lambda"
    else:
        centre.helicity = "mixed"
    centre.helicity_detail = (
        f"{len(signs)} chelate pair(s), skew angles "
        + ", ".join(f"{a:+.0f}" for a in signs)
        + " deg"
    )


def skew_angle(a1, a2, b1, b2) -> float | None:
    """Signed angle (deg, -90..90) turning line a1-a2 onto line b1-b2.

    Measured about the common perpendicular pointing from the first line to
    the second, positive by the right-hand rule. The result does not depend on
    the order of the two lines or on the direction in which either is
    traversed, which is what makes it a property of the pair of chelates
    rather than of how they were listed.
    """
    u, v = np.asarray(a2) - np.asarray(a1), np.asarray(b2) - np.asarray(b1)
    normal = np.cross(u, v)
    if np.linalg.norm(normal) < 1e-8:
        return None  # parallel
    axis = normal / np.linalg.norm(normal)
    if np.dot(axis, np.asarray(b1) - np.asarray(a1)) < 0:
        axis = -axis
    angle = float(np.degrees(np.arctan2(np.dot(normal, axis), np.dot(u, v))))
    if angle > 90:
        angle -= 180
    elif angle < -90:
        angle += 180
    return angle


def _angle(u: np.ndarray, v: np.ndarray) -> float:
    cosine = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))
