"""Charge and spin multiplicity.

Getting these two integers wrong is the most common way a quantum-chemistry
input fails, so they are derived explicitly and any manual override that
disagrees with the derivation is reported rather than silently accepted.
"""

from __future__ import annotations

from rdkit import Chem

from ..types import IssueLog


def net_charge(mol: Chem.Mol) -> int:
    """Sum of the formal charges."""
    return Chem.GetFormalCharge(mol)


def total_electrons(mol: Chem.Mol, charge: int | None = None) -> int:
    """Electron count including implicit hydrogens."""
    if charge is None:
        charge = net_charge(mol)
    electrons = 0
    for atom in mol.GetAtoms():
        electrons += atom.GetAtomicNum()
        electrons += atom.GetTotalNumHs()  # implicit/explicit H are 1 electron each
    return electrons - charge


def radical_electrons(mol: Chem.Mol) -> int:
    return sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())


def multiplicity(mol: Chem.Mol, charge: int | None = None) -> int:
    """Spin multiplicity 2S+1.

    Radical electrons flagged on the molecule win; otherwise the parity of the
    electron count decides between a closed-shell singlet and a doublet.
    """
    radicals = radical_electrons(mol)
    if radicals:
        return radicals + 1
    return 1 if total_electrons(mol, charge) % 2 == 0 else 2


def resolve(
    mol: Chem.Mol,
    log: IssueLog,
    *,
    charge_override: int | None = None,
    multiplicity_override: int | None = None,
) -> tuple[int, int]:
    """Return the (charge, multiplicity) to write, reporting disagreements."""
    derived_charge = net_charge(mol)
    charge = derived_charge if charge_override is None else charge_override
    if charge_override is not None and charge_override != derived_charge:
        log.warn(
            "electronic.charge_override",
            f"Charge forced to {charge_override} but the structure adds up to "
            f"{derived_charge}. Check the formal charges in the drawing.",
        )

    derived_mult = multiplicity(mol, charge)
    mult = derived_mult if multiplicity_override is None else multiplicity_override
    if multiplicity_override is not None and multiplicity_override != derived_mult:
        log.info(
            "electronic.multiplicity_override",
            f"Multiplicity forced to {multiplicity_override} (structure suggests "
            f"{derived_mult}).",
        )

    electrons = total_electrons(mol, charge)
    unpaired = mult - 1
    if (electrons - unpaired) % 2 != 0:
        log.error(
            "electronic.impossible",
            f"Charge {charge} with multiplicity {mult} is impossible: "
            f"{electrons} electrons cannot leave {unpaired} unpaired. "
            "Adjust --charge or --mult.",
        )
    elif mult > 1 and multiplicity_override is None:
        log.warn(
            "electronic.open_shell",
            f"Open-shell species (multiplicity {mult}); the calculation will need "
            "an unrestricted reference.",
        )

    radicals = radical_electrons(mol)
    if radicals:
        log.info(
            "electronic.radicals",
            f"{radicals} radical electron(s) detected in the structure.",
        )

    return charge, mult
