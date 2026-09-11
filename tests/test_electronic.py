"""Charge and multiplicity derivation."""

from __future__ import annotations

import pytest
from conftest import codes

from m2i.chem import electronic, sanitize


@pytest.mark.parametrize(
    "smiles, charge, multiplicity",
    [
        ("CC", 0, 1),  # ethane
        ("CC(=O)[O-]", -1, 1),  # acetate
        ("C[NH3+]", 1, 1),  # methylammonium
        ("[CH3]", 0, 2),  # methyl radical
        ("[O-]C(=O)C(=O)[O-]", -2, 1),  # oxalate
        ("[NH4+]", 1, 1),
        ("c1ccccc1", 0, 1),
    ],
)
def test_derived_state(smiles, charge, multiplicity, log):
    mol = sanitize.mol_from_smiles(smiles)
    assert electronic.resolve(mol, log) == (charge, multiplicity)


def test_electron_count_includes_implicit_hydrogens():
    mol = sanitize.mol_from_smiles("CC")  # C2H6: 2*6 + 6 = 18
    assert electronic.total_electrons(mol) == 18


def test_forced_charge_that_contradicts_the_drawing_warns(log):
    mol = sanitize.mol_from_smiles("CC(=O)[O-]")
    charge, _ = electronic.resolve(mol, log, charge_override=0)
    assert charge == 0
    assert "electronic.charge_override" in codes(log)


def test_impossible_spin_state_is_an_error(log):
    mol = sanitize.mol_from_smiles("CC")  # 18 electrons cannot leave 1 unpaired
    electronic.resolve(mol, log, multiplicity_override=2)
    assert "electronic.impossible" in codes(log)
    assert log.has_errors()


def test_triplet_oxygen_can_be_forced_without_complaint(log):
    """O=O parses as a singlet; the chemist has to say otherwise, and may."""
    mol = sanitize.mol_from_smiles("O=O")
    charge, multiplicity = electronic.resolve(mol, log, multiplicity_override=3)
    assert (charge, multiplicity) == (0, 3)
    assert not log.has_errors()


def test_radical_cation_combines_charge_and_spin(log):
    mol = sanitize.mol_from_smiles("[CH3+]")
    charge, multiplicity = electronic.resolve(mol, log)
    assert charge == 1
    assert multiplicity == 1  # methyl cation: 8 electrons, closed shell
