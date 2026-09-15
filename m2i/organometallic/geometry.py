"""Coordination polyhedra, their isomers, and which one a drawing shows.

A 2D drawing of a metal complex fixes the connectivity but only hints at the
arrangement: which ligand is trans to which, and -- for chiral arrangements --
which enantiomer. The hints are the directions of the metal-ligand bonds on
the page, plus wedges and hashes for ligands in front of or behind it.

So every distinct arrangement compatible with the chelates is enumerated, each
is scored by how well it can be rotated onto the drawn directions, and the
best is proposed. When two arrangements fit about equally well the drawing
does not decide, and the user is asked.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import permutations

import numpy as np

_S = 1 / np.sqrt(2)

#: Unit vectors of the ideal polyhedra.
POLYHEDRA = {
    "linear": np.array([[1, 0, 0], [-1, 0, 0]], float),
    "trigonal planar": np.array(
        [[1, 0, 0], [-0.5, np.sqrt(3) / 2, 0], [-0.5, -np.sqrt(3) / 2, 0]]
    ),
    "tetrahedral": np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float) / np.sqrt(3),
    "square planar": np.array([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]], float),
    "trigonal bipyramidal": np.array(
        [[0, 0, 1], [0, 0, -1], [1, 0, 0], [-0.5, np.sqrt(3) / 2, 0], [-0.5, -np.sqrt(3) / 2, 0]]
    ),
    "square pyramidal": np.array(
        [[0, 0, 1], [1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]], float
    ),
    "octahedral": np.array(
        [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float
    ),
}

#: The polyhedra that can hold each number of donors, most common first.
CANDIDATES = {
    2: ("linear",),
    3: ("trigonal planar",),
    4: ("tetrahedral", "square planar"),
    5: ("trigonal bipyramidal", "square pyramidal"),
    6: ("octahedral",),
}

#: Donors closer than this (in degrees) on the polyhedron are cis.
TRANS = 150.0
#: Two arrangements fitting the drawing within this RMS of each other are a tie.
AMBIGUOUS = 0.15


@dataclass
class Arrangement:
    polyhedron: str
    #: donor atom index -> unit vector of its site
    sites: dict[int, np.ndarray]
    #: RMS misfit to the drawn directions (0 = perfect)
    misfit: float

    def trans_pairs(self) -> set[frozenset]:
        donors = list(self.sites)
        return {
            frozenset((a, b))
            for i, a in enumerate(donors)
            for b in donors[i + 1:]
            if _angle(self.sites[a], self.sites[b]) >= TRANS
        }

    def handedness(self) -> float:
        """Sign of the arrangement: flips for its mirror image (0 if achiral-looking)."""
        return _handedness(self.sites)


def arrangements(
    donors: list[int],
    drawn: dict[int, np.ndarray],
    chelated: list[tuple[int, int]],
    polyhedron: str | None = None,
    labels: dict[int, str] | None = None,
) -> list[Arrangement]:
    """Every distinct arrangement, best fit to the drawing first.

    ``chelated`` lists donor pairs joined by a short chelate ring: they must
    be cis. Arrangements that differ only by a rotation are the same isomer
    and appear once; mirror images are kept apart, since they are the two
    enantiomers of a chiral arrangement. ``labels`` names chemically
    equivalent donors alike (two CO ligands), so swapping them does not make
    a new isomer; the arrangement kept is the one that fits the drawing best.
    """
    labels = labels or {d: str(d) for d in donors}
    shapes = [polyhedron] if polyhedron else list(CANDIDATES.get(len(donors), ()))
    found: list[Arrangement] = []
    for shape in shapes:
        vertices = POLYHEDRA[shape]
        if len(vertices) != len(donors):
            continue
        # Every permutation of six donors is tried, so the inner loop is
        # arithmetic on vertex numbers: the geometry (which rotations map the
        # polyhedron onto itself, which vertex pairs are trans) is worked out
        # once per polyhedron and cached.
        maps = _vertex_maps(shape)
        trans = _trans_vertices(shape)
        seen: dict = {}
        for perm in permutations(range(len(donors))):
            vertex_of = dict(zip(donors, perm))
            if any((vertex_of[a], vertex_of[b]) in trans for a, b in chelated):
                continue
            key = _key(vertex_of, maps, labels, chelated)
            sites = {d: vertices[v] for d, v in vertex_of.items()}
            candidate = Arrangement(shape, sites, _fit(sites, drawn))
            if key in seen:
                # The same isomer with equivalent donors swapped: keep the
                # labelling that matches the drawing best.
                if candidate.misfit < found[seen[key]].misfit:
                    found[seen[key]] = candidate
                continue
            seen[key] = len(found)
            found.append(candidate)
    found.sort(key=lambda a: a.misfit)
    return found


def is_chiral(arrangement: Arrangement, labels: dict[int, str], chelated=()) -> bool:
    """Whether the mirror image is a different arrangement of these donors.

    The chelates count: [Co(en)3]3+ has six equal donors, and is chiral only
    because of which pairs the three rings join.
    """
    vertices = POLYHEDRA[arrangement.polyhedron]
    mirror = {d: v * np.array([-1.0, 1.0, 1.0]) for d, v in arrangement.sites.items()}
    return _signature(arrangement.sites, vertices, labels, chelated) != _signature(
        mirror, vertices, labels, chelated
    )


def mirror_partners(options: list[Arrangement], labels: dict[int, str], chelated=()) -> dict[int, int]:
    """Position in ``options`` -> position of its mirror image, for the chiral ones."""
    keys = {}
    for number, option in enumerate(options):
        vertices = POLYHEDRA[option.polyhedron]
        keys[(option.polyhedron, _signature(option.sites, vertices, labels, chelated))] = number
    partners = {}
    for number, option in enumerate(options):
        vertices = POLYHEDRA[option.polyhedron]
        mirror = {d: v * np.array([-1.0, 1.0, 1.0]) for d, v in option.sites.items()}
        other = keys.get((option.polyhedron, _signature(mirror, vertices, labels, chelated)))
        if other is not None and other != number:
            partners[number] = other
    return partners


def drawn_directions(metal_xy, donor_xy: dict[int, tuple], depth: dict[int, float]) -> dict[int, np.ndarray]:
    """3D guesses of the bond directions from a drawing.

    ``depth`` is +1 for a wedge (towards the viewer), -1 for a hash, 0 for a
    plain line; a wedged bond is drawn foreshortened, so its out-of-plane
    component is as large as its in-plane one.
    """
    out = {}
    for donor, (x, y) in donor_xy.items():
        v = np.array([x - metal_xy[0], y - metal_xy[1], 0.0])
        length = np.linalg.norm(v[:2]) or 1.0
        v[2] = depth.get(donor, 0.0) * length
        out[donor] = v / np.linalg.norm(v)
    return out


# -- internals -----------------------------------------------------------------


def _angle(u, v) -> float:
    return float(np.degrees(np.arccos(np.clip(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)), -1, 1))))


def _fit(sites: dict[int, np.ndarray], drawn: dict[int, np.ndarray]) -> float:
    """RMS distance between the drawn directions and the best proper rotation of the sites."""
    donors = [d for d in sites if d in drawn]
    if len(donors) < 2:
        return 0.0
    A = np.array([sites[d] for d in donors])
    B = np.array([drawn[d] for d in donors])
    U, _, Vt = np.linalg.svd(A.T @ B)
    D = np.eye(3)
    D[2, 2] = np.sign(np.linalg.det(U @ Vt)) or 1.0
    R = U @ D @ Vt
    return float(np.sqrt(((A @ R - B) ** 2).sum(axis=1).mean()))


def _signature(sites: dict[int, np.ndarray], vertices: np.ndarray, labels: dict[int, str], chelated=()):
    """The key of an arrangement given as vectors (the stored ``sites``)."""
    shape = next(name for name, corners in POLYHEDRA.items() if corners is vertices)
    vertex_of = {d: int(np.argmin(np.linalg.norm(vertices - v, axis=1))) for d, v in sites.items()}
    return _key(vertex_of, _vertex_maps(shape), labels, chelated)


def _key(vertex_of: dict[int, int], maps, labels: dict[int, str], chelated):
    """The same key for arrangements related by a proper rotation of the
    polyhedron, or by swapping equivalent donors: the label found at each
    vertex, and the vertex pairs the chelates join, minimised over the
    rotations."""
    best = None
    for mapping in maps:
        at_vertex = [""] * len(mapping)
        for d, v in vertex_of.items():
            at_vertex[mapping[v]] = labels[d]
        rings = tuple(sorted(
            tuple(sorted((mapping[vertex_of[a]], mapping[vertex_of[b]]))) for a, b in chelated
        ))
        key = (tuple(at_vertex), rings)
        if best is None or key < best:
            best = key
    return best


@lru_cache(maxsize=None)
def _vertex_maps(shape: str) -> tuple[tuple[int, ...], ...]:
    """Each proper rotation of the polyhedron, as a vertex -> vertex map."""
    vertices = POLYHEDRA[shape]
    maps = []
    for rotation in _rotations(vertices):
        turned = vertices @ rotation.T
        maps.append(tuple(
            int(np.argmin(np.linalg.norm(vertices - v, axis=1))) for v in turned
        ))
    return tuple(maps)


@lru_cache(maxsize=None)
def _trans_vertices(shape: str) -> frozenset:
    """Vertex pairs too far apart to be cis, in both orders."""
    vertices = POLYHEDRA[shape]
    pairs = set()
    for a in range(len(vertices)):
        for b in range(len(vertices)):
            if a != b and _angle(vertices[a], vertices[b]) >= TRANS:
                pairs.add((a, b))
    return frozenset(pairs)


_ROTATION_CACHE: dict[bytes, list[np.ndarray]] = {}


def _rotations(vertices: np.ndarray) -> list[np.ndarray]:
    """Proper rotations mapping the polyhedron onto itself."""
    key = vertices.tobytes()
    if key in _ROTATION_CACHE:
        return _ROTATION_CACHE[key]
    rotations = []
    n = len(vertices)
    for a in range(n):
        for b in range(n):
            if a == b:
                continue
            for i in range(n):
                for j in range(n):
                    if i == j:
                        continue
                    if abs(np.dot(vertices[a], vertices[b]) - np.dot(vertices[i], vertices[j])) > 1e-6:
                        continue
                    R = _frame(vertices[i], vertices[j]) @ _frame(vertices[a], vertices[b]).T
                    if np.linalg.det(R) < 0.5:
                        continue
                    mapped = vertices @ R.T
                    if all(np.min(np.linalg.norm(vertices - m, axis=1)) < 1e-6 for m in mapped)                             and not any(np.allclose(R, r) for r in rotations):
                        rotations.append(R)
    _ROTATION_CACHE[key] = rotations
    return rotations


def _frame(u, v) -> np.ndarray:
    e1 = u / np.linalg.norm(u)
    w = v - np.dot(v, e1) * e1
    if np.linalg.norm(w) < 1e-9:
        w = np.cross(e1, [1.0, 0.0, 0.0])
        if np.linalg.norm(w) < 1e-9:
            w = np.cross(e1, [0.0, 1.0, 0.0])
    e2 = w / np.linalg.norm(w)
    return np.column_stack([e1, e2, np.cross(e1, e2)])


def _handedness(sites: dict[int, np.ndarray]) -> float:
    """A chirality sign from the first three non-coplanar donors (in index order)."""
    donors = sorted(sites)
    for i in range(len(donors)):
        for j in range(i + 1, len(donors)):
            for k in range(j + 1, len(donors)):
                volume = np.linalg.det(np.array([sites[donors[i]], sites[donors[j]], sites[donors[k]]]))
                if abs(volume) > 0.1:
                    return float(np.sign(volume))
    return 0.0
