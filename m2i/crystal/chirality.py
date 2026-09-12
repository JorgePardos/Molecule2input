"""Chirality that CIP centre labels do not cover: planar chirality, and whether
the crystal holds one enantiomer or both.

Planar chirality (Schlögl's convention for metallocenes and arene complexes):
a eta-5/eta-6 ring bearing two different substituents in a 1,2 or 1,3
pattern has no mirror plane once one face is bound to the metal. The ring is
viewed from the side away from the metal; if the path from the substituent of
highest CIP priority to the next one runs clockwise, it is Rp, else Sp. Other
conventions exist -- naming the ring carbon as a stereocentre gives "(1R)"
style labels that need not agree -- so the convention is always stated.

What is decided here without doubt is *whether* a ring is planar chiral and
which way round its substituents sit. The label needs the two substituents
ranked by CIP, and a CIF has no bond orders: when a ranking depends on them it
is reported as uncertain, so the user can be asked.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .structure import METALS, Species

#: Typical number of neighbours (hydrogens included) of a saturated atom.
#: Fewer means a multiple bond is possible, whose duplicate atoms a CIF
#: cannot provide -- the source of doubt in a CIP ranking.
SATURATED_VALENCE = {"C": 4, "N": 3, "O": 2, "S": 2, "P": 3, "B": 3, "Si": 4}
#: How far along a substituent the CIP-like comparison looks.
RANKING_DEPTH = 6


@dataclass
class PlanarChirality:
    ring_label: str
    #: Labels of the substituent atoms, highest CIP priority first.
    first: str
    second: str
    descriptor: str  # "Rp" or "Sp"
    certain: bool
    reason: str

    def describe(self) -> str:
        doubt = "" if self.certain else f" -- uncertain: {self.reason}"
        return (
            f"planar chiral ring {self.ring_label}: {self.descriptor} (Schlögl convention; "
            f"{self.first} ranked above {self.second}){doubt}"
        )

    def flipped(self) -> "PlanarChirality":
        """The same ring if the user ranks the two substituents the other way."""
        return PlanarChirality(
            ring_label=self.ring_label,
            first=self.second,
            second=self.first,
            descriptor="Sp" if self.descriptor == "Rp" else "Rp",
            certain=True,
            reason="priority set by the user",
        )


def planar_chirality(species: Species, metal: int, ring: list[int]) -> PlanarChirality | None:
    """Rp/Sp for one eta-5/eta-6 ring bound to ``metal``, or None if achiral."""
    order = _ring_order(species, ring)
    if order is None:
        return None
    ring_set = set(order)
    exo = {atom: _exo_atom(species, atom, ring_set) for atom in order}
    labels = [_substituent_signature(species, exo[a], ring_set) for a in order]
    if _has_mirror(labels):
        return None

    substituted = [a for a in order if exo[a] is not None and species.elements[exo[a]] != "H"]
    if len(substituted) < 2:
        return None
    ranked, certain, reason = _rank(species, [exo[a] for a in substituted], ring_set)
    top, second = ranked[0], ranked[1]
    carbon_of = {exo[a]: a for a in substituted}

    centre = species.coords[order].mean(axis=0)
    away = centre - species.coords[metal]
    away /= np.linalg.norm(away)
    u = species.coords[carbon_of[top]] - centre
    v = species.coords[carbon_of[second]] - centre
    # Counterclockwise seen from the far side (right-hand rule about the axis
    # pointing away from the metal) is positive; clockwise is Rp.
    turn = float(np.dot(away, np.cross(u, v)))
    return PlanarChirality(
        ring_label=species.labels[order[0]],
        first=species.labels[top],
        second=species.labels[second],
        descriptor=RP_WHEN_CLOCKWISE if turn < 0 else SP_WHEN_COUNTERCLOCKWISE,
        certain=certain,
        reason=reason,
    )


# The two constants make the convention explicit in one place; the mapping is
# checked in the tests against resolved COD structures labelled by their
# authors, and against their mirror images.
RP_WHEN_CLOCKWISE = "Rp"
SP_WHEN_COUNTERCLOCKWISE = "Sp"


def is_racemic_crystal(operations) -> bool:
    """True if the space group has an improper operation (inversion, mirror,
    glide): the crystal then holds both enantiomers of any chiral molecule."""
    for op in operations:
        if np.linalg.det(np.array(op.rot, dtype=float)) < 0:
            return True
    return False


def mirrored(species: Species) -> Species:
    """The enantiomer: every x coordinate reflected."""
    coords = species.coords.copy()
    coords[:, 0] *= -1
    return Species(
        elements=list(species.elements),
        coords=coords,
        labels=list(species.labels),
        bonds=list(species.bonds),
        copies=species.copies,
    )


# -- rings and substituents ----------------------------------------------


def _ring_order(species: Species, ring: list[int]) -> list[int] | None:
    ring_set = set(ring)
    neighbours = {a: [b for b in species.neighbours(a) if b in ring_set] for a in ring}
    if any(len(n) != 2 for n in neighbours.values()):
        return None  # not a simple cycle (fused, or partly bound)
    order, previous = [ring[0]], None
    while len(order) < len(ring):
        here = order[-1]
        following = [b for b in neighbours[here] if b != previous]
        previous = here
        order.append(following[0])
    return order


def _exo_atom(species: Species, atom: int, ring: set[int]) -> int | None:
    outside = [
        b for b in species.neighbours(atom)
        if b not in ring and species.elements[b] not in METALS
    ]
    heavy = [b for b in outside if species.elements[b] != "H"]
    if heavy:
        return heavy[0]
    return outside[0] if outside else None


def _substituent_signature(species: Species, root: int | None, ring: set[int]) -> tuple:
    """Sphere-by-sphere composition of a substituent: equal only for equal groups."""
    if root is None:
        return ()
    spheres = []
    seen, frontier = {root} | ring, [root]
    for _ in range(RANKING_DEPTH):
        if not frontier:
            break
        spheres.append(tuple(sorted(species.elements[a] for a in frontier)))
        following = []
        for atom in frontier:
            for other in species.neighbours(atom):
                if other not in seen and species.elements[other] not in METALS:
                    seen.add(other)
                    following.append(other)
        frontier = following
    return tuple(spheres)


def _has_mirror(labels: list) -> bool:
    """With the metal face fixed, the only symmetry that can remove chirality
    is a mirror containing the ring normal: it reverses the ring sequence."""
    n = len(labels)
    for k in range(n):
        if all(labels[i] == labels[(k - i) % n] for i in range(n)):
            return True
    return False


def _rank(species: Species, roots: list[int], ring: set[int]):
    """Order substituent root atoms by an approximate CIP comparison.

    Branches are compared sphere by sphere on atomic numbers, which is exact
    whenever the first difference is in the heaviest atom of a sphere: the
    duplicate atoms of a multiple bond copy an atom already present, so they
    can never exceed that maximum. A difference found further down a sphere
    is reliable only if no atom above it could carry a multiple bond.
    """
    from rdkit import Chem

    table = Chem.GetPeriodicTable()

    def spheres(root):
        rows, parents_saturated = [], []
        seen, frontier = {root} | ring, [root]
        saturated = True
        for _ in range(RANKING_DEPTH):
            if not frontier:
                break
            rows.append(sorted((table.GetAtomicNumber(species.elements[a]) for a in frontier), reverse=True))
            parents_saturated.append(saturated)
            following = []
            for atom in frontier:
                partners = [b for b in species.neighbours(atom) if species.elements[b] not in METALS]
                expected = SATURATED_VALENCE.get(species.elements[atom])
                if expected is not None and len(partners) < expected:
                    saturated = False
                for other in partners:
                    if other not in seen:
                        seen.add(other)
                        following.append(other)
            frontier = following
        return rows, parents_saturated

    explored = {root: spheres(root) for root in roots}

    def compare(a, b):
        rows_a, sat_a = explored[a]
        rows_b, sat_b = explored[b]
        for level in range(max(len(rows_a), len(rows_b))):
            ra = rows_a[level] if level < len(rows_a) else []
            rb = rows_b[level] if level < len(rows_b) else []
            if ra == rb:
                continue
            width = max(len(ra), len(rb))
            ra, rb = ra + [0] * (width - len(ra)), rb + [0] * (width - len(rb))
            first = next(i for i in range(width) if ra[i] != rb[i])
            safe = first == 0 or (
                (level >= len(sat_a) or sat_a[level]) and (level >= len(sat_b) or sat_b[level])
            )
            return (1 if ra > rb else -1), safe
        return 0, False

    import functools

    doubts = []

    def key(a, b):
        result, safe = compare(a, b)
        if not safe:
            doubts.append((a, b))
        return -result

    ranked = sorted(roots, key=functools.cmp_to_key(key))
    top_two = set(ranked[:2])
    relevant = [pair for pair in doubts if set(pair) <= top_two]
    if relevant:
        a, b = ranked[0], ranked[1]
        return ranked, False, (
            f"which of {species.labels[a]} and {species.labels[b]} has the higher CIP "
            "priority depends on bond orders, which a CIF does not record"
        )
    return ranked, True, ""
