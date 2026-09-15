"""Electronic state of metal-containing species: what to ask, and what to check.

For an organic molecule m2i derives the charge from the formal charges in
the structure and the multiplicity from the electron count. Neither works for
a metal complex. The oxidation state is not in the coordinates, and for most
d-electron counts the spin state -- high, intermediate or low -- is a genuine
chemical question that only the chemist can answer. So these are asked, and
this module turns the answers into sensible options and catches the one
mistake that can be checked with certainty: a multiplicity whose parity is
impossible for the number of electrons, which Gaussian and ORCA would reject
only after the job has been queued.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from rdkit import Chem

#: Group number of the d-block elements.
_D_BLOCK = {
    **dict.fromkeys(("Sc", "Y", "La", "Ac"), 3),
    **dict.fromkeys(("Ti", "Zr", "Hf"), 4),
    **dict.fromkeys(("V", "Nb", "Ta"), 5),
    **dict.fromkeys(("Cr", "Mo", "W"), 6),
    **dict.fromkeys(("Mn", "Tc", "Re"), 7),
    **dict.fromkeys(("Fe", "Ru", "Os"), 8),
    **dict.fromkeys(("Co", "Rh", "Ir"), 9),
    **dict.fromkeys(("Ni", "Pd", "Pt"), 10),
    **dict.fromkeys(("Cu", "Ag", "Au"), 11),
    **dict.fromkeys(("Zn", "Cd", "Hg"), 12),
}
_LANTHANIDES = "Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu".split()
_ACTINIDES = "Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr".split()

#: Ions whose charge is not in doubt, keyed by Hill formula. Used only to
#: *suggest* the charge of the complex they balance.
COMMON_IONS = {
    "F": -1, "Cl": -1, "Br": -1, "I": -1,
    "BF4": -1, "F6P": -1, "F6Sb": -1, "AsF6": -1, "ClO4": -1, "NO3": -1,
    "CF3O3S": -1,  # triflate
    "C24H20B": -1,  # BPh4-
    "C32H12BF24": -1,  # BArF24-
    "O4S": -2, "CO3": -2,
    "Li": 1, "Na": 1, "K": 1, "Rb": 1, "Cs": 1,
    "H4N": 1,  # ammonium
    "C4H12N": 1,  # NMe4+
    "C8H20N": 1,  # NEt4+
    "C16H36N": 1,  # NBu4+
    "C24H20P": 1,  # PPh4+
}
#: Common crystallisation solvents, which carry no charge.
NEUTRAL_SOLVENTS = {
    "H2O", "CH4O", "C2H6O", "C3H6O", "C2H3N", "CH2Cl2", "CHCl3", "C4H8O",
    "C4H10O", "C6H6", "C7H8", "C6H14", "C5H12", "C2H6OS", "C3H7NO", "C4H8O2",
}


@dataclass
class SpinOption:
    multiplicity: int
    unpaired: int
    label: str


def valence_electron_count(symbol: str, oxidation_state: int) -> tuple[str, int] | None:
    """('d', n) or ('f', n) for a metal in a given oxidation state."""
    table = Chem.GetPeriodicTable()
    if symbol in _D_BLOCK:
        count = _D_BLOCK[symbol] - oxidation_state
        return ("d", count) if 0 <= count <= 10 else None
    if symbol in _LANTHANIDES:
        count = table.GetAtomicNumber(symbol) - 54 - oxidation_state
        return ("f", count) if 0 <= count <= 14 else None
    if symbol in _ACTINIDES:
        count = table.GetAtomicNumber(symbol) - 86 - oxidation_state
        return ("f", count) if 0 <= count <= 14 else None
    return None


def spin_options(symbol: str, oxidation_state: int) -> list[SpinOption]:
    """Every spin state a single metal centre can take, low to high."""
    shell = valence_electron_count(symbol, oxidation_state)
    if shell is None:
        return []
    kind, n = shell
    capacity = 10 if kind == "d" else 14
    highest = min(n, capacity - n)
    options = []
    for unpaired in range(highest % 2, highest + 1, 2):
        if highest == highest % 2:
            label = "only possibility"
        elif unpaired == highest % 2:
            label = "low spin"
        elif unpaired == highest:
            label = "high spin"
        else:
            label = "intermediate spin"
        options.append(SpinOption(multiplicity=unpaired + 1, unpaired=unpaired, label=label))
    return options


def electron_count(elements, charge: int) -> int:
    table = Chem.GetPeriodicTable()
    return sum(table.GetAtomicNumber(e) for e in elements) - charge


def parity_problem(elements, charge: int, multiplicity: int) -> str | None:
    """Why this charge/multiplicity pair is impossible, or None if it is fine."""
    if multiplicity < 1:
        return "the multiplicity must be at least 1"
    electrons = electron_count(elements, charge)
    if electrons < 0:
        return f"a charge of {charge:+d} leaves no electrons"
    if (electrons % 2 == 0) != (multiplicity % 2 == 1):
        parity = "even" if electrons % 2 == 0 else "odd"
        allowed = "odd (1, 3, 5...)" if electrons % 2 == 0 else "even (2, 4, 6...)"
        return (
            f"with charge {charge:+d} the species has {electrons} electrons ({parity}), "
            f"so the multiplicity must be {allowed}; {multiplicity} is impossible"
        )
    if multiplicity - 1 > electrons:
        return f"{multiplicity - 1} unpaired electrons is more than the {electrons} there are"
    return None


def suggest_charge(main_formula: str, main_copies: int, others) -> tuple[int | None, str]:
    """A charge for the main species from the ions that balance it, if possible.

    ``others`` is a sequence of (formula, copies in the cell). The suggestion
    is offered only when every other species is a known ion or a neutral
    solvent; anything unrecognised means the charge could be anything.
    """
    balance, unknown, used = 0, [], []
    for formula, copies in others:
        if formula in COMMON_IONS:
            balance += COMMON_IONS[formula] * copies
            used.append(f"{copies} x {formula} ({COMMON_IONS[formula]:+d})")
        elif formula not in NEUTRAL_SOLVENTS:
            unknown.append(formula)
    if unknown:
        return None, (
            f"cannot suggest a charge: the crystal also holds {', '.join(unknown)}, "
            "whose charge m2i does not know"
        )
    if balance == 0:
        return 0, "no counter-ions in the crystal, so the species is presumably neutral"
    if balance % main_copies:
        return None, f"the counter-ions ({'; '.join(used)}) do not divide evenly among the {main_copies} copies"
    charge = -balance // main_copies
    return charge, f"balances the counter-ions in the cell: {'; '.join(used)} for {main_copies} x {main_formula}"


def metal_summary(elements) -> dict[str, int]:
    from ..crystal.structure import METALS

    return dict(Counter(e for e in elements if e in METALS))
