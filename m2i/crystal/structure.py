"""Discrete molecules out of a crystal structure (CIF).

A CIF does not describe a molecule. It describes the asymmetric unit of a
periodic crystal: a fraction of the unit cell that the symmetry operations
complete. A complex sitting on an inversion centre is stored as half a
complex; a molecule straddling the cell edge is stored in pieces. So before a
single coordinate can go into an input file, the crystal has to be rebuilt:

1. keep one component of any disorder (the major one);
2. apply the symmetry to fill the unit cell, merging atoms that symmetry
   places on the same special position;
3. find the bonds, across cell boundaries, from covalent radii;
4. walk each bonded network to reassemble whole molecules in one piece --
   or discover that the network never ends (a polymer or a framework), which
   has no molecule to calculate and is reported as such.

What comes out is experimental geometry. Nothing is generated, so the
stereochemistry and the shape of every ligand are exactly those measured.
"""

from __future__ import annotations

import itertools
import re
from collections import Counter, deque
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..types import IssueLog

#: Added to the sum of covalent radii when deciding whether two atoms bond.
BOND_TOLERANCE = 0.45  # Å
#: Two symmetry-generated atoms closer than this are the same site.
SAME_SITE = 0.3  # Å
#: Two atoms of one molecule closer than this betray unresolved disorder.
OVERLAP = 0.6  # Å

#: Standard X-H lengths (neutron-normalised, as used by the CSD and Mercury).
#: X-ray places H at the centre of its electron density, ~0.1 Å too short.
XH_LENGTHS = {"C": 1.089, "N": 1.015, "O": 0.993, "B": 1.190}

METALS = frozenset(
    """
    Li Be Na Mg Al K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Rb Sr Y Zr Nb Mo Tc Ru
    Rh Pd Ag Cd In Sn Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta
    W Re Os Ir Pt Au Hg Tl Pb Bi Po Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md
    No Lr
    """.split()
)


class CrystalError(ValueError):
    """The file cannot be turned into molecules."""


@dataclass
class Species:
    """One kind of discrete molecule or ion found in the crystal."""

    elements: list[str]
    #: Cartesian coordinates in Å, the whole molecule in one piece.
    coords: np.ndarray
    labels: list[str]
    bonds: list[tuple[int, int]]
    #: How many copies the unit cell holds.
    copies: int = 1

    @property
    def n_atoms(self) -> int:
        return len(self.elements)

    @property
    def formula(self) -> str:
        return hill_formula(self.elements)

    @property
    def metals(self) -> list[str]:
        return [e for e in self.elements if e in METALS]

    @property
    def has_hydrogen(self) -> bool:
        return "H" in self.elements

    def neighbours(self, index: int) -> list[int]:
        return [b if a == index else a for a, b in self.bonds if index in (a, b)]


@dataclass
class CrystalReading:
    path: Path
    block: str
    spacegroup: str
    cell: tuple[float, ...]
    z: float | None
    formula_sum: str | None
    cod_id: str | None
    species: list[Species] = field(default_factory=list)
    #: Formulas of bonded networks that never close: polymers, frameworks.
    extended: list[str] = field(default_factory=list)
    #: Per formula unit, hydrogens the CIF declares but does not contain.
    missing_hydrogens: float = 0.0
    disorder: str | None = None
    hydrogens_normalised: int = 0
    #: The space group has an inversion, mirror or glide: both enantiomers of
    #: any chiral molecule are in the crystal.
    racemic: bool = False
    #: Flack parameter as written in the CIF, e.g. "0.02(3)".
    flack: str | None = None

    @property
    def main(self) -> int:
        """Index of the species to calculate by default: the metal complex."""
        return 0


# -- reading ---------------------------------------------------------------


def read_cif(path: Path, log: IssueLog, *, normalise_hydrogens: bool = True) -> CrystalReading:
    """A CIF file in, the discrete molecules of the crystal out."""
    try:
        import gemmi
    except ImportError as exc:
        raise CrystalError(
            "reading CIF files needs gemmi, which is not installed: "
            'pip install gemmi  (or pip install -e ".[crystal]")'
        ) from exc

    path = Path(path)
    try:
        # Bytes, not the path: gemmi, like RDKit, is happier without
        # non-ASCII folder names on Windows.
        doc = gemmi.cif.read_string(path.read_bytes().decode("utf-8", errors="replace"))
    except Exception as exc:  # noqa: BLE001 - gemmi raises assorted types
        raise CrystalError(f"{path.name}: not a readable CIF ({exc})") from exc

    block = _structure_block(doc)
    if block is None:
        raise CrystalError(
            f"{path.name}: no atom coordinates found (no _atom_site_fract_x loop)."
        )
    structure = gemmi.make_small_structure_from_block(block)
    if not structure.sites:
        raise CrystalError(f"{path.name}: the atom list is empty.")

    reading = CrystalReading(
        path=path,
        block=block.name,
        spacegroup=structure.spacegroup_hm or "unknown",
        cell=tuple(
            round(v, 4)
            for v in (
                structure.cell.a, structure.cell.b, structure.cell.c,
                structure.cell.alpha, structure.cell.beta, structure.cell.gamma,
            )
        ),
        z=_number(block.find_value("_cell_formula_units_Z")),
        formula_sum=_text(block.find_value("_chemical_formula_sum")),
        cod_id=_text(block.find_value("_cod_database_code")),
        flack=_text(block.find_value("_refine_ls_abs_structure_Flack")),
    )

    sites = _select_disorder(list(structure.sites), reading, log)
    operations = _operations(structure, log)
    from .chirality import is_racemic_crystal

    reading.racemic = is_racemic_crystal(operations)
    labels, elements, fractional = _fill_cell(sites, operations, structure.cell)
    matrix = np.array(structure.cell.orth.mat.tolist())

    edges = _bonds(elements, fractional, matrix)
    components = _assemble(elements, fractional, edges)
    _collect_species(reading, components, labels, elements, fractional, matrix, log)
    _check_formula(reading, elements, log)
    if normalise_hydrogens:
        for species in reading.species:
            reading.hydrogens_normalised += normalise_xh(species)
    _report(reading, log)
    return reading


def _structure_block(doc):
    for block in doc:
        if block.find_value("_atom_site_fract_x") is not None or block.find_loop(
            "_atom_site_fract_x"
        ):
            return block
    return None


def _operations(structure, log: IssueLog):
    """Symmetry operations: the ones listed in the file, else the space group's."""
    import gemmi

    if structure.symops:
        return [gemmi.Op(text) for text in structure.symops]
    if structure.spacegroup is not None:
        return list(structure.spacegroup.operations())
    log.warn(
        "crystal.no_symmetry",
        "The CIF lists no symmetry operations and no recognisable space group; "
        "only the asymmetric unit is used, so molecules may be incomplete.",
    )
    return [gemmi.Op("x,y,z")]


def _select_disorder(sites, reading: CrystalReading, log: IssueLog):
    """Keep ordered atoms plus the disorder component with the most occupancy.

    Choosing 'group 1' blindly is wrong surprisingly often: in the ferrocene
    structure used by the tests the major component is group 2.
    """
    groups: dict[int, list] = {}
    for site in sites:
        if site.disorder_group > 0:
            groups.setdefault(site.disorder_group, []).append(site)

    if not groups:
        partial = [s.label for s in sites if s.occ < 0.99]
        if partial:
            log.warn(
                "crystal.partial_occupancy",
                f"{len(partial)} atom(s) have partial occupancy but no disorder "
                f"group ({', '.join(partial[:6])}{'...' if len(partial) > 6 else ''}). "
                "All are kept; check the structure for overlapping atoms.",
            )
        return sites

    occupancy = {g: float(np.mean([s.occ for s in members])) for g, members in groups.items()}
    chosen = max(occupancy, key=lambda g: (occupancy[g], -g))
    dropped = sorted(g for g in groups if g != chosen)
    reading.disorder = (
        f"disorder group {chosen} kept (occupancy {occupancy[chosen]:.2f}); "
        + ", ".join(f"group {g} ({occupancy[g]:.2f})" for g in dropped)
        + " dropped"
    )
    log.warn(
        "crystal.disorder",
        f"The structure is disordered. Kept the major component: {reading.disorder}.",
    )
    negative = [s.label for s in sites if s.disorder_group < 0]
    if negative:
        log.warn(
            "crystal.special_disorder",
            f"{len(negative)} atom(s) are disordered across a symmetry element "
            f"({', '.join(negative[:6])}); check the result for overlapping atoms.",
        )
    return [s for s in sites if s.disorder_group <= 0 or s.disorder_group == chosen]


def _fill_cell(sites, operations, cell):
    """Apply every operation, wrap into the cell, merge coincident sites."""
    matrix = np.array(cell.orth.mat.tolist())
    labels: list[str] = []
    elements: list[str] = []
    fractional: list[np.ndarray] = []
    seen: Counter = Counter()

    for op_index, op in enumerate(operations):
        for site in sites:
            frac = np.array(op.apply_to_xyz(list(site.fract.tolist()))) % 1.0
            element = _element(site)
            if _already_there(frac, element, elements, fractional, matrix):
                continue
            seen[site.label] += 1
            copy = seen[site.label]
            labels.append(site.label if copy == 1 else f"{site.label}_{copy}")
            elements.append(element)
            fractional.append(frac)
    return labels, elements, np.array(fractional)


def _element(site) -> str:
    symbol = site.element.name
    return "H" if symbol == "D" else symbol


def _already_there(frac, element, elements, fractional, matrix) -> bool:
    if not fractional:
        return False
    existing = np.array(fractional)
    delta = existing - frac
    delta -= np.round(delta)
    distance = np.linalg.norm(delta @ matrix.T, axis=1)
    close = np.where(distance < SAME_SITE)[0]
    return any(elements[i] == element for i in close)


def _bonds(elements, fractional, matrix):
    """Bonded pairs across the periodic boundary: (i, j, lattice shift of j)."""
    from rdkit import Chem

    table = Chem.GetPeriodicTable()
    radii = np.array([table.GetRcovalent(table.GetAtomicNumber(e)) for e in elements])
    cutoff = radii[:, None] + radii[None, :] + BOND_TOLERANCE
    hydrogen = np.array([e == "H" for e in elements])
    not_hh = ~(hydrogen[:, None] & hydrogen[None, :])
    cart = fractional @ matrix.T

    edges = []
    for shift in itertools.product((-1, 0, 1), repeat=3):
        shifted = (fractional + np.array(shift)) @ matrix.T
        distance = np.linalg.norm(cart[:, None, :] - shifted[None, :, :], axis=2)
        bonded = (distance < cutoff) & (distance > 0.1) & not_hh
        for i, j in zip(*np.nonzero(bonded)):
            edges.append((int(i), int(j), shift))
    return edges


def _assemble(elements, fractional, edges):
    """Connected components, each with the lattice image of every atom.

    An atom reached twice through different images means the network closes
    on itself only after crossing the cell: it is infinite.
    """
    adjacency: dict[int, list] = {i: [] for i in range(len(elements))}
    for i, j, shift in edges:
        adjacency[i].append((j, np.array(shift)))

    image: dict[int, np.ndarray] = {}
    components = []
    for start in range(len(elements)):
        if start in image:
            continue
        image[start] = np.zeros(3, dtype=int)
        members, periodic = [start], False
        queue = deque([start])
        while queue:
            i = queue.popleft()
            for j, shift in adjacency[i]:
                wanted = image[i] + shift
                if j not in image:
                    image[j] = wanted
                    members.append(j)
                    queue.append(j)
                elif not np.array_equal(image[j], wanted):
                    periodic = True
        components.append((sorted(members), periodic, {m: image[m] for m in members}))
    return components


def _collect_species(reading, components, labels, elements, fractional, matrix, log):
    by_formula: dict[str, list] = {}
    for members, periodic, images in components:
        formula = hill_formula([elements[i] for i in members])
        if periodic:
            if formula not in reading.extended:
                reading.extended.append(formula)
            continue
        by_formula.setdefault(formula, []).append((members, images))

    for formula, copies in by_formula.items():
        # The copy holding the lowest-numbered atom contains the asymmetric
        # unit as published, completed by symmetry.
        members, images = min(copies, key=lambda c: c[0][0])
        index = {atom: k for k, atom in enumerate(members)}
        coords = np.array([(fractional[a] + images[a]) @ matrix.T for a in members])
        bonds = sorted(
            {
                tuple(sorted((index[a], index[b])))
                for a in members
                for b in _bonded_partners(a, members, coords, index, elements)
            }
        )
        species = Species(
            elements=[elements[a] for a in members],
            coords=coords,
            labels=[labels[a] for a in members],
            bonds=bonds,
            copies=len(copies),
        )
        _check_overlaps(species, log)
        reading.species.append(species)

    reading.species.sort(key=lambda s: (-len(s.metals), -sum(e != "H" for e in s.elements)))
    if reading.extended:
        log.warn(
            "crystal.extended",
            f"Part of this crystal is an extended network with no discrete molecule: "
            f"{', '.join(reading.extended)} repeats without end across the cell "
            "(a coordination polymer or framework). It cannot be calculated as a "
            "molecule; a cluster cut from it would need its own capping choices.",
        )
    if not reading.species:
        raise CrystalError(
            f"{reading.path.name}: no discrete molecule in this crystal -- "
            "everything belongs to an extended network."
        )


def _bonded_partners(atom, members, coords, index, elements):
    from rdkit import Chem

    table = Chem.GetPeriodicTable()
    here = coords[index[atom]]
    r_here = table.GetRcovalent(table.GetAtomicNumber(elements[atom]))
    for other in members:
        if other == atom or (elements[atom] == "H" and elements[other] == "H"):
            continue
        r_other = table.GetRcovalent(table.GetAtomicNumber(elements[other]))
        if np.linalg.norm(coords[index[other]] - here) < r_here + r_other + BOND_TOLERANCE:
            yield other


def _check_overlaps(species: Species, log: IssueLog) -> None:
    coords = species.coords
    distance = np.linalg.norm(coords[:, None] - coords[None, :], axis=2)
    np.fill_diagonal(distance, np.inf)
    i, j = np.unravel_index(np.argmin(distance), distance.shape)
    if distance[i, j] < OVERLAP:
        log.error(
            "crystal.overlap",
            f"In {species.formula}, {species.labels[i]} and {species.labels[j]} are "
            f"{distance[i, j]:.2f} Å apart: two positions of a disordered atom were "
            "both kept. This geometry cannot be calculated as it is.",
        )


# -- checks and corrections ------------------------------------------------


def _check_formula(reading: CrystalReading, elements: list[str], log: IssueLog) -> None:
    """Compare what was rebuilt with the formula the authors declared.

    It is the single best check on everything above -- symmetry, disorder,
    bond detection -- and it is how missing hydrogens are found: old X-ray
    structures often never located them.
    """
    if not reading.formula_sum or not reading.z:
        log.info(
            "crystal.formula_unchecked",
            "The CIF gives no _chemical_formula_sum or Z; the reconstruction cannot "
            "be checked against the declared formula.",
        )
        return

    declared = parse_formula(reading.formula_sum)
    found = {e: n / reading.z for e, n in Counter(elements).items()}
    missing_h = declared.get("H", 0) - found.get("H", 0)
    others = {
        e
        for e in set(declared) | set(found)
        if e != "H" and abs(declared.get(e, 0) - found.get(e, 0)) > 0.01
    }

    if missing_h > 0.01:
        reading.missing_hydrogens = missing_h
        # A warning, not an error: m2i proposes where they go, and the input
        # is only held back if the species to be calculated is the one short.
        log.warn(
            "crystal.missing_hydrogens",
            f"The formula declares {declared.get('H', 0):g} H per formula unit, but "
            f"the file contains {found.get('H', 0):g}: {missing_h:g} hydrogen(s) were "
            "never located in the experiment. Without them the molecule would have "
            "another electron count and charge; m2i proposes where they go.",
        )
    if others:
        detail = ", ".join(
            f"{e}: declared {declared.get(e, 0):g}, found {found.get(e, 0):g}"
            for e in sorted(others)
        )
        log.warn(
            "crystal.formula_mismatch",
            f"The rebuilt crystal does not match _chemical_formula_sum per formula "
            f"unit ({detail}). Solvent left out of the model, or disorder, can do "
            "this; check the species list.",
        )
    if missing_h <= 0.01 and not others:
        log.info(
            "crystal.formula_ok",
            f"Rebuilt contents match the declared formula ({reading.formula_sum}, Z = {reading.z:g}).",
        )


def normalise_xh(species: Species) -> int:
    """Move each terminal H along its bond to the standard X-H length."""
    moved = 0
    for h in range(species.n_atoms):
        if species.elements[h] != "H":
            continue
        partners = species.neighbours(h)
        if len(partners) != 1:
            continue
        x = partners[0]
        length = XH_LENGTHS.get(species.elements[x])
        if length is None:
            continue  # hydrides and other X-H: leave as measured
        vector = species.coords[h] - species.coords[x]
        norm = np.linalg.norm(vector)
        if norm > 0 and abs(norm - length) > 1e-3:
            species.coords[h] = species.coords[x] + vector / norm * length
            moved += 1
    return moved


def _report(reading: CrystalReading, log: IssueLog) -> None:
    if reading.hydrogens_normalised:
        log.info(
            "crystal.xh_normalised",
            f"{reading.hydrogens_normalised} C-H/N-H/O-H distance(s) extended to "
            "standard neutron values (X-ray places H about 0.1 Å too close). Hydrides "
            "and other M-H were left as measured.",
        )


# -- formulas --------------------------------------------------------------


def hill_formula(elements) -> str:
    counts = Counter(elements)
    order = []
    if "C" in counts:
        order = ["C"] + (["H"] if "H" in counts else [])
    order += sorted(e for e in counts if e not in order)
    return "".join(f"{e}{counts[e] if counts[e] > 1 else ''}" for e in order)


def parse_formula(text: str) -> dict[str, float]:
    counts: dict[str, float] = {}
    for symbol, number in re.findall(r"([A-Z][a-z]?)\s*([0-9]*\.?[0-9]*)", text):
        counts[symbol] = counts.get(symbol, 0) + (float(number) if number else 1.0)
    return counts


def _number(value) -> float | None:
    text = _text(value)
    if not text:
        return None
    try:
        return float(re.sub(r"\(.*\)", "", text))
    except ValueError:
        return None


def _text(value) -> str | None:
    if value is None:
        return None
    text = str(value).strip().strip("'\"")
    return None if text in ("", "?", ".") else text
