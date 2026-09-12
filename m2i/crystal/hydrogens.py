"""Put back the hydrogens an X-ray experiment did not locate.

Placing a hydrogen once its count and the hybridisation of its atom are known
is standard geometry (it is what SHELXL's HFIX does). The hard part is the
count, because a CIF records neither bond orders nor protonation:

- for carbon, and for nitrogen in organic parts, the heavy-atom geometry
  settles it well -- angles and planarity tell sp3 from sp2 from sp, and bond
  lengths do so for terminal atoms;
- for atoms bound to a metal it often does not. An O bound only to the metal
  may be an oxo, a hydroxo or an aqua ligand; the M-O distance points to one,
  but the ranges overlap -- and the choice changes the charge.

So every count is proposed with the reason behind it, the doubtful ones are
marked, and the total is checked against the hydrogens the formula declares.
When the proposal matches the formula it can be applied as it is; when it
does not, the doubtful atoms are exactly the questions to ask.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .structure import METALS, CrystalReading, Species

#: X-H bond lengths used for placed hydrogens (neutron-normalised).
XH = {"C": 1.089, "N": 1.015, "O": 0.993, "S": 1.338, "B": 1.19, "Si": 1.48, "P": 1.42}
#: Lengths of single bonds from carbon, to tell them from shorter, multiple ones.
SINGLE_FROM_C = {"C": 1.50, "N": 1.46, "O": 1.41, "S": 1.80, "P": 1.83}
TETRAHEDRAL = 109.47


@dataclass
class HydrogenSite:
    atom: int
    label: str
    element: str
    count: int
    #: 'tetrahedral', 'trigonal', 'linear', 'bent' or 'water': how to place them.
    mode: str
    certain: bool
    reason: str
    options: tuple[int, ...] = ()


@dataclass
class HydrogenPlan:
    sites: list[HydrogenSite] = field(default_factory=list)
    #: Hydrogens one copy of the species needs according to the formula,
    #: when that can be worked out; None otherwise.
    target: int | None = None

    @property
    def total(self) -> int:
        return sum(s.count for s in self.sites)

    @property
    def matches(self) -> bool:
        return self.target is not None and self.total == self.target

    @property
    def doubtful(self) -> list[HydrogenSite]:
        return [s for s in self.sites if not s.certain]

    def with_counts(self, counts: dict[str, int]) -> "HydrogenPlan":
        """The same plan with some counts set by the user (by atom label)."""
        sites = []
        for site in self.sites:
            if site.label in counts:
                site = HydrogenSite(**{**site.__dict__, "count": int(counts[site.label]),
                                       "certain": True, "reason": "set by the user"})
            sites.append(site)
        return HydrogenPlan(sites=sites, target=self.target)


# -- deciding how many --------------------------------------------------------


def plan(reading: CrystalReading, index: int) -> HydrogenPlan:
    """Where the missing hydrogens of one species probably go."""
    species = reading.species[index]
    proposal = HydrogenPlan(sites=_sites(species))
    if not reading.missing_hydrogens or not reading.z:
        return proposal

    missing_in_cell = round(reading.missing_hydrogens * reading.z)
    needing = {
        i: sum(s.count for s in _sites(other)) for i, other in enumerate(reading.species)
    }
    needing = {i: n for i, n in needing.items() if n}
    if sum(reading.species[i].copies * n for i, n in needing.items()) == missing_in_cell:
        # The proposals, species by species, account for exactly the missing
        # hydrogens: each is then consistent with the formula.
        proposal.target = needing.get(index, 0)
    elif set(needing) == {index} and not reading.extended:
        # Only this species can hold them, so it must hold them all -- and a
        # proposal that falls short is then visibly wrong.
        per_copy, remainder = divmod(missing_in_cell, species.copies)
        proposal.target = per_copy if not remainder else None
    # Otherwise the missing hydrogens cannot be attributed: they may sit on a
    # network, or on atoms whose geometry proposed none. No target is given
    # rather than a wrong one.
    return proposal


def lacking(reading: CrystalReading, index: int) -> bool:
    """Whether *this* species is short of hydrogens.

    The formula check is crystal-wide. Often the hydrogens that were not
    located belong to lattice water or another solvent, and the complex
    itself is complete: it must not be held back for them.
    """
    if not reading.missing_hydrogens:
        return False
    proposal = plan(reading, index)
    return proposal.total > 0 or bool(proposal.target)


def _sites(species: Species) -> list[HydrogenSite]:
    haptic = _haptic_atoms(species)
    sites = []
    for atom, element in enumerate(species.elements):
        if element in METALS or element == "H":
            continue
        site = _infer(species, atom, haptic)
        if site is not None and (site.count or not site.certain):
            sites.append(site)
    return sites


def _haptic_atoms(species: Species) -> set[int]:
    """Atoms bonded to a metal together with a neighbour bound to the same metal."""
    found = set()
    for atom, element in enumerate(species.elements):
        if element in METALS:
            continue
        for metal in (n for n in species.neighbours(atom) if species.elements[n] in METALS):
            if any(
                other in species.neighbours(metal)
                for other in species.neighbours(atom)
                if species.elements[other] not in METALS
            ):
                found.add(atom)
    return found


def _infer(species: Species, atom: int, haptic: set[int]) -> HydrogenSite | None:
    """How many hydrogens this atom lacks, from the geometry of its neighbours.

    Every rule states the *total* number of neighbours the atom should have
    (a sigma-bound metal included, a face-on ring's metal not); the hydrogens
    to add are that total minus the neighbours already there.
    """
    element = species.elements[atom]
    neighbours = species.neighbours(atom)
    have_h = sum(species.elements[n] == "H" for n in neighbours)
    metals = [n for n in neighbours if species.elements[n] in METALS]
    heavy = [n for n in neighbours if species.elements[n] != "H" and species.elements[n] not in METALS]
    sigma = heavy + ([] if atom in haptic else metals)

    def site(total, mode, certain, reason, options=()):
        count = max(0, total - len(sigma) - have_h)
        return HydrogenSite(atom, species.labels[atom], element, count, mode, certain,
                            reason, tuple(options))

    if element == "C":
        return _carbon(species, atom, sigma, metals, haptic, site)
    if element == "N":
        return _nitrogen(species, atom, sigma, heavy, metals, site)
    if element == "O":
        return _oxygen(species, atom, sigma, metals, site)
    if element == "B":
        return _boron(species, atom, sigma, site)
    if element == "S" and len(sigma) == 1 and not metals and species.elements[sigma[0]] == "C":
        length = _distance(species, atom, sigma[0])
        if length >= 1.76:
            return site(2, "bent", False, f"C-S {length:.2f} A: a thiol SH", (1, 0))
    return None  # halogens, P, B, S otherwise: nothing added by default


def _carbon(species, atom, sigma, metals, haptic, site):
    k = len(sigma)
    if atom in haptic:
        return site(3, "trigonal", True, "carbon of a ring bound face-on to the metal: sp2")
    if k >= 4:
        return None
    if k == 3:
        total = _angle_sum(species, atom, sigma)
        if total >= 350:
            return None
        return site(4, "tetrahedral", True, f"three neighbours, pyramidal ({total:.0f} deg): CH")
    if k == 2:
        angle = _angle(species, sigma[0], atom, sigma[1])
        if angle >= 170:
            return None  # sp, no hydrogen
        # An sp2 CH needs two pieces of evidence, not one: a wide angle *and*
        # a bond shorter than single. Either alone misleads -- thermal motion
        # shortens the bonds of an ethyl group or an ethylenediamine chelate,
        # and a strained chelate opens the angle of a plain CH2. In a planar
        # ring (thiophene, pyrrole) a short bond is enough, since there the
        # angle is small by construction.
        short = _has_short_bond(species, atom, sigma)
        ring = _smallest_ring(species, atom)
        planar = ring is not None and _is_planar(species, ring)
        if short and (angle >= 116 or planar):
            return site(3, "trigonal", angle >= 118 or planar,
                        f"two neighbours at {angle:.0f} deg, a short bond"
                        + (f", in a planar {len(ring)}-ring" if planar else "") + ": sp2 CH",
                        (1, 2))
        if angle >= 116:
            single = all(
                _distance(species, atom, n) >= SINGLE_FROM_C.get(species.elements[n], 1.50) - 0.03
                for n in sigma
                if species.elements[n] not in METALS
            )
            return site(4, "tetrahedral", single,
                        f"two neighbours at {angle:.0f} deg but single bonds: a CH2 in a "
                        "strained chain", (2, 1))
        return site(4, "tetrahedral", angle <= 114,
                    f"two neighbours at {angle:.0f} deg: sp3 CH2", (2, 1))
    if k == 1:
        other = sigma[0]
        partner = species.elements[other]
        length = _distance(species, atom, other)
        if partner in METALS:
            return site(4, "tetrahedral", False,
                        f"bound only to {partner} ({length:.2f} A): a methyl, or a carbene",
                        (3, 2, 1))
        single = SINGLE_FROM_C.get(partner, 1.50)
        if length >= single - 0.04:
            return site(4, "tetrahedral", True, f"terminal, C-{partner} {length:.2f} A: CH3")
        if length >= single - 0.14 and _is_sp3(species, other, atom):
            # A terminal carbon librates, and X-ray then shortens its bond. A
            # vinyl's inner carbon would be trigonal; this one is tetrahedral.
            return site(4, "tetrahedral", True,
                        f"terminal, C-{partner} {length:.2f} A (short, typical of thermal "
                        f"motion) on an sp3 {partner}: CH3")
        if length >= single - 0.22:
            return site(3, "trigonal", False, f"terminal, C-{partner} {length:.2f} A: =CH2", (2, 3))
        return site(2, "linear", False, f"terminal, C-{partner} {length:.2f} A: a triple bond", (1, 0))
    return site(4, "tetrahedral", False, "an isolated carbon", (4,))


def _boron(species, atom, sigma, site):
    k = len(sigma)
    if k >= 4:
        return None
    if k == 3:
        total = _angle_sum(species, atom, sigma)
        if total >= 350:
            return None  # a trigonal borane BR3
        return site(4, "tetrahedral", False, f"pyramidal boron ({total:.0f} deg): BH", (1, 0))
    if k == 1:
        partner = species.elements[sigma[0]]
        return site(4, "tetrahedral", False, f"bound only to {partner}: a borane adduct, BH3", (3, 2))
    if k == 2:
        return site(4, "tetrahedral", False, "two neighbours: BH2", (2, 1))
    return None


def _is_linear(species, atom) -> bool:
    around = [n for n in species.neighbours(atom) if species.elements[n] != "H"]
    return len(around) == 2 and _angle(species, around[0], atom, around[1]) >= 165


def _is_sp3(species, atom, excluding) -> bool:
    """Whether a neighbour's own geometry is tetrahedral."""
    around = [n for n in species.neighbours(atom) if species.elements[n] != "H"]
    if len(around) >= 4:
        return True
    if len(around) == 3:
        return _angle_sum(species, atom, around) < 345
    if len(around) == 2:
        return _angle(species, around[0], atom, around[1]) < 114
    return False


def _nitrogen(species, atom, sigma, heavy, metals, site):
    k = len(sigma)
    bound = bool(metals)
    if k >= 4:
        return None
    if k == 3:
        total = _angle_sum(species, atom, sigma)
        if total >= 350:
            return None  # planar: amide, aniline or coordinated pyridine-type N
        if bound:
            return site(4, "tetrahedral", True, f"coordinated and pyramidal ({total:.0f} deg): NH")
        return site(3, "tetrahedral", False, f"pyramidal amine ({total:.0f} deg), or ammonium", (0, 1))
    if k == 2:
        angle = _angle(species, sigma[0], atom, sigma[1])
        if angle >= 165:
            return None  # nitrile or isocyanide
        ring = _ring_size(species, atom)
        if bound:
            carbon = [n for n in heavy if species.elements[n] == "C"]
            length = _distance(species, atom, carbon[0]) if carbon else 1.47
            if ring:
                return None  # a ring nitrogen bound to the metal: pyridine-, imidazole-type
            if length >= 1.44:
                return site(4, "tetrahedral", True, f"coordinated amine, N-C {length:.2f} A: NH2")
            return site(3, "trigonal", False, f"coordinated, N-C {length:.2f} A: an imine NH", (1, 0))
        if ring == 6 and angle >= 114:
            return None  # pyridine-type
        if ring == 5:
            return site(3, "trigonal", False,
                        f"in a five-membered ring ({angle:.0f} deg): a pyrrole-type NH, "
                        "or an imidazole-type N", (1, 0))
        if _has_short_bond(species, atom, sigma, margin=0.14):
            return None  # an imine C=N-R
        return site(3, "trigonal" if angle >= 116 else "tetrahedral", True,
                    f"two neighbours at {angle:.0f} deg: NH")
    if k == 1:
        other = sigma[0]
        partner = species.elements[other]
        length = _distance(species, atom, other)
        if partner in METALS:
            if length < 1.75:
                return site(1, "linear", False, f"{partner}-N {length:.2f} A: nitrido or imido", (0, 1))
            return site(4, "tetrahedral", False,
                        f"bound only to {partner} ({length:.2f} A): an ammine NH3, or an amido",
                        (3, 2))
        if length < 1.20 or (partner == "C" and _is_linear(species, other)):
            # A nitrile. Its C-N length can come out distorted in a crystal
            # (a disordered acetonitrile, say), but a linear carbon cannot.
            return None
        # An amide NH2 has C-N about 1.32 A and an imine =NH about 1.28 A:
        # close enough that the band between them stays a question.
        if length < 1.28:
            return site(2, "trigonal", False, f"N-{partner} {length:.2f} A: =NH", (1, 2))
        return site(3, "trigonal" if length < 1.42 else "tetrahedral", length >= 1.31,
                    f"terminal, N-{partner} {length:.2f} A: NH2", (2, 1))
    return site(3, "tetrahedral", False, "an isolated nitrogen: ammonia, or ammonium", (3, 4))


def _oxygen(species, atom, sigma, metals, site):
    k = len(sigma)
    if k == 0:
        return site(2, "water", True, "an isolated oxygen: water", (2, 1))
    if k == 1:
        other = sigma[0]
        partner = species.elements[other]
        length = _distance(species, atom, other)
        if partner in METALS:
            if length < 1.75:
                return site(1, "bent", False, f"{partner}-O {length:.2f} A: an oxo", (0, 1, 2))
            if length < 1.98:
                return site(2, "bent", False, f"{partner}-O {length:.2f} A: a hydroxo", (1, 0, 2))
            return site(3, "trigonal", False, f"{partner}-O {length:.2f} A: an aqua", (2, 1))
        if partner == "C":
            if length >= 1.36:
                return site(2, "bent", True, f"C-O {length:.2f} A: a hydroxyl")
            if length <= 1.28:
                return None  # carbonyl or carboxylate
            return site(2 if length >= 1.32 else 1, "bent", False,
                        f"C-O {length:.2f} A: an acid OH, or a carboxylate/phenolate O", (1, 0))
        if partner in ("P", "S") and length >= 1.55:
            return site(2, "bent", False, f"{partner}-O {length:.2f} A: {partner}-OH", (1, 0))
        return None
    if k == 2 and metals:
        if len(metals) == 2:
            shortest = min(_distance(species, atom, m) for m in metals)
            return site(2 if shortest < 1.85 else 3, "bent", False,
                        f"bridging two metals (M-O {shortest:.2f} A): mu-oxo or mu-hydroxo", (0, 1))
        metal_o = _distance(species, atom, metals[0])
        # A neutral alcohol binds more weakly, and so further away, than an
        # alkoxide: a long M-O points to the alcohol.
        if metal_o >= 2.15:
            return site(3, "bent", False,
                        f"bound to a metal ({metal_o:.2f} A, long) and a carbon: a coordinated "
                        "alcohol", (1, 0))
        return site(2, "bent", False,
                    f"bound to a metal ({metal_o:.2f} A) and a carbon: an alkoxide", (0, 1))
    return None


# -- placing them --------------------------------------------------------------


def apply(species: Species, proposal: HydrogenPlan) -> tuple[Species, list[str]]:
    """A copy of the species with the proposed hydrogens placed."""
    elements = list(species.elements)
    coords = [np.array(c) for c in species.coords]
    labels = list(species.labels)
    bonds = list(species.bonds)
    added = []

    for site in proposal.sites:
        if site.count <= 0:
            continue
        directions = _directions(species, site, np.array(coords))
        length = XH.get(site.element, 1.0)
        for n, direction in enumerate(directions):
            position = coords[site.atom] + direction * length
            elements.append("H")
            coords.append(position)
            labels.append(f"H{site.label}{chr(ord('A') + n)}")
            bonds.append((site.atom, len(elements) - 1))
        added.append(f"{site.label}: +{site.count} H ({site.reason})")

    return (
        Species(elements=elements, coords=np.array(coords), labels=labels,
                bonds=bonds, copies=species.copies),
        added,
    )


def _directions(species: Species, site: HydrogenSite, coords: np.ndarray) -> list[np.ndarray]:
    here = coords[site.atom]
    units = [
        (coords[n] - here) / np.linalg.norm(coords[n] - here)
        for n in species.neighbours(site.atom)
    ]
    n = site.count

    if site.mode == "water" or not units:
        return _free(n)
    if site.mode == "linear":
        return [-units[0]]
    if len(units) >= 2 and n == 1:
        opposite = -np.sum(units, axis=0)
        if np.linalg.norm(opposite) < 1e-6:
            opposite = np.cross(units[0], units[1])
        return [opposite / np.linalg.norm(opposite)]
    if len(units) == 2 and n == 2:
        bisector = -(units[0] + units[1])
        bisector /= np.linalg.norm(bisector)
        normal = np.cross(units[0], units[1])
        normal /= np.linalg.norm(normal)
        half = np.radians(TETRAHEDRAL / 2)
        return [bisector * np.cos(half) + s * normal * np.sin(half) for s in (1, -1)]
    if len(units) == 1:
        neighbour = species.neighbours(site.atom)[0]
        if site.mode == "trigonal":
            # =CH2, an amide NH2: in the plane of the neighbour's substituents.
            plane = _plane_reference(species, site.atom, neighbour, coords)
            return _cone(units[0], 120.0, n, 180.0 if n == 2 else 0.0, plane) if plane is not None \
                else _best_rotor(species, site.atom, coords, units[0], 120.0, n, 180.0)
        angle = TETRAHEDRAL if site.mode == "tetrahedral" else 109.0
        return _best_rotor(species, site.atom, coords, units[0], angle, n, 360.0 / n)
    return _free(n)


def _plane_reference(species, atom, neighbour, coords):
    """A direction in the plane of the neighbour's other substituents."""
    axis = coords[neighbour] - coords[atom]
    axis /= np.linalg.norm(axis)
    for other in species.neighbours(neighbour):
        if other == atom:
            continue
        v = coords[other] - coords[neighbour]
        v = v - np.dot(v, axis) * axis
        if np.linalg.norm(v) > 1e-3:
            return v / np.linalg.norm(v)
    return None


def _cone(axis, angle, n, step, reference):
    other = np.cross(axis, reference)
    theta = np.radians(angle)
    result = []
    for k in range(n):
        phi = np.radians(k * step)
        radial = reference * np.cos(phi) + other * np.sin(phi)
        result.append(axis * np.cos(theta) + radial * np.sin(theta))
    return result


def _best_rotor(species, atom, coords, axis, angle, n, step):
    """Hydrogens on a cone around one bond, turned to stay clear of other atoms."""
    reference = _perpendicular(axis)
    other = np.cross(axis, reference)
    here = coords[atom]
    bonded = set(species.neighbours(atom)) | {atom}
    others = np.array([c for i, c in enumerate(coords) if i not in bonded]) if len(coords) > len(bonded) else None
    best, best_score = None, -np.inf
    for offset in range(0, 360, 10):
        directions = []
        for k in range(n):
            phi = np.radians(offset + k * step)
            radial = reference * np.cos(phi) + other * np.sin(phi)
            theta = np.radians(angle)
            directions.append(axis * np.cos(theta) + radial * np.sin(theta))
        if others is None or not len(others):
            return directions
        tips = [here + d for d in directions]
        score = min(np.min(np.linalg.norm(others - tip, axis=1)) for tip in tips)
        if score > best_score:
            best, best_score = directions, score
    return best


def _free(n: int) -> list[np.ndarray]:
    """Directions for an atom with no neighbour to orient by (e.g. lattice water)."""
    half = np.radians(104.5 / 2)
    if n == 1:
        return [np.array([0.0, 0.0, 1.0])]
    if n == 2:
        return [np.array([np.sin(half), 0.0, np.cos(half)]), np.array([-np.sin(half), 0.0, np.cos(half)])]
    corners = [np.array(v, dtype=float) for v in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))]
    return [v / np.linalg.norm(v) for v in corners[:n]]


# -- geometry helpers ---------------------------------------------------------


def _perpendicular(v: np.ndarray) -> np.ndarray:
    trial = np.array([1.0, 0.0, 0.0]) if abs(v[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    p = np.cross(v, trial)
    return p / np.linalg.norm(p)


def _distance(species, a, b) -> float:
    return float(np.linalg.norm(species.coords[a] - species.coords[b]))


def _angle(species, a, centre, b) -> float:
    u = species.coords[a] - species.coords[centre]
    v = species.coords[b] - species.coords[centre]
    cosine = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
    return float(np.degrees(np.arccos(np.clip(cosine, -1, 1))))


def _angle_sum(species, atom, neighbours) -> float:
    return sum(
        _angle(species, neighbours[i], atom, neighbours[j])
        for i in range(len(neighbours))
        for j in range(i + 1, len(neighbours))
    )


def _has_short_bond(species, atom, neighbours, margin: float = 0.06) -> bool:
    """A bond clearly shorter than single: conjugation or a multiple bond."""
    if species.elements[atom] != "C":
        return any(
            species.elements[n] == "C"
            and _distance(species, atom, n) < SINGLE_FROM_C.get(species.elements[atom], 1.46) - margin
            for n in neighbours
        )
    return any(
        _distance(species, atom, n) < SINGLE_FROM_C.get(species.elements[n], 1.50) - margin
        for n in neighbours
        if species.elements[n] not in METALS
    )


def _smallest_ring(species, atom, limit: int = 7) -> list[int] | None:
    """Atoms of the smallest ring through ``atom`` (non-metal, non-H bonds)."""
    def usable(n):
        return species.elements[n] not in METALS and species.elements[n] != "H"

    best = None
    for first in [n for n in species.neighbours(atom) if usable(n)]:
        parents = {first: atom}
        frontier = [first]
        for _ in range(limit - 1):
            following = []
            for here in frontier:
                for n in species.neighbours(here):
                    if not usable(n) or n == parents[here]:
                        continue
                    if n == atom and here != first:
                        path, node = [], here
                        while node != atom:
                            path.append(node)
                            node = parents[node]
                        ring = [atom] + path[::-1]
                        if best is None or len(ring) < len(best):
                            best = ring
                    elif n not in parents and n != atom:
                        parents[n] = here
                        following.append(n)
            frontier = following
    return best


def _is_planar(species, ring: list[int], tolerance: float = 0.1) -> bool:
    """All ring atoms within ``tolerance`` Å of their best plane."""
    points = species.coords[ring] - species.coords[ring].mean(axis=0)
    normal = np.linalg.svd(points)[2][-1]
    return float(np.max(np.abs(points @ normal))) < tolerance


def _ring_size(species, atom, limit: int = 7) -> int | None:
    """Size of the smallest ring through ``atom`` (non-metal bonds only)."""
    start_neighbours = [
        n for n in species.neighbours(atom)
        if species.elements[n] not in METALS and species.elements[n] != "H"
    ]
    best = None
    for first in start_neighbours:
        frontier, seen, depth = [first], {atom, first}, 1
        while frontier and depth < limit:
            depth += 1
            following = []
            for here in frontier:
                for n in species.neighbours(here):
                    if species.elements[n] in METALS or species.elements[n] == "H":
                        continue
                    if n == atom and depth > 2:
                        best = depth if best is None else min(best, depth)
                    elif n not in seen:
                        seen.add(n)
                        following.append(n)
            frontier = following
    return best
