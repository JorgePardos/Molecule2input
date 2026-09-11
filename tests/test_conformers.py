"""3D embedding, pruning, and the stereochemistry safety net."""

from __future__ import annotations

import math

import pytest
from conftest import codes
from rdkit import Chem

from m2i.chem import conformers as conf_mod
from m2i.chem import sanitize
from m2i.chem import stereo as stereo_mod

CHIRAL_SMILES = [
    "F[C@@H](Cl)Br",
    "C[C@H](N)C(=O)O",
    "C[C@H](O)/C=C/C(F)Cl",
    r"C/C=C\C",
    "O[C@H]1CC[C@@H](O)CC1",
    "C[C@@H]1CC[C@H](C)CC1",
    "N[C@@H](Cc1ccccc1)C(=O)O",
    "C[C@H](O)[C@@H](N)C(=O)O",
]


def embed(smiles: str, log, **kwargs):
    mol = sanitize.mol_from_smiles(smiles)
    stereo_mod.assign(mol)
    reference = stereo_mod.summarize(mol)
    options = conf_mod.ConformerOptions(**kwargs)
    mol_h, conformers = conf_mod.generate(mol, log, reference, options)
    return mol, reference, mol_h, conformers


@pytest.mark.parametrize("smiles", CHIRAL_SMILES)
def test_every_written_conformer_keeps_the_intended_configuration(smiles, log):
    _, reference, mol_h, conformers = embed(smiles, log, keep=5)
    assert conformers

    for conformer in conformers:
        single = Chem.Mol(mol_h)
        single.RemoveAllConformers()
        single.AddConformer(
            Chem.Conformer(mol_h.GetConformer(conformer.index)), assignId=True
        )
        stereo_mod.assign_from_3d(single)
        assert stereo_mod.compare(reference, stereo_mod.summarize(single)) == []


def test_conformers_come_back_ordered_by_energy(log):
    _, _, _, conformers = embed("CC(C)CC(C)CC(C)CO", log, keep=5)
    energies = [c.energy for c in conformers]
    assert energies == sorted(energies)
    assert conformers[0].relative_energy == pytest.approx(0.0)
    assert all(c.relative_energy >= -1e-9 for c in conformers)


def test_returned_molecule_holds_exactly_the_selected_conformers(log):
    """Index i of the list must be conformer i of the molecule, or the SDF and
    the .gjf would describe different geometries."""
    _, _, mol_h, conformers = embed("CC(C)CC(C)CC(C)CO", log, keep=3)
    assert mol_h.GetNumConformers() == len(conformers) == 3

    for conformer in conformers:
        positions = mol_h.GetConformer(conformer.index).GetPositions()
        for written, actual in zip(conformer.coords, positions):
            assert written == pytest.approx(list(actual))


def test_distinct_conformers_are_separated_by_the_rms_threshold(log):
    threshold = 0.5
    _, _, mol_h, conformers = embed("CCCCCCO", log, keep=6, prune_rms=threshold)
    assert len(conformers) > 1

    probe = conf_mod._rmsd_probe(mol_h)
    for i in range(len(conformers)):
        for j in range(i + 1, len(conformers)):
            assert conf_mod._rms(probe, i, j, True) > threshold


def test_polar_hydrogen_rotamers_survive_pruning(log):
    """Pruning on heavy atoms alone would collapse ethanol to one conformer."""
    _, _, _, conformers = embed("CCO", log, keep=5)
    assert len(conformers) >= 2


def test_stereochemistry_violation_is_caught_and_refused(log):
    """If the requested configuration cannot be built, refuse rather than write
    the wrong enantiomer."""
    mol = sanitize.mol_from_smiles("F[C@@H](Cl)Br")
    stereo_mod.assign(mol)
    reference = stereo_mod.summarize(mol)
    assert reference.centers[0].label == "S"
    reference.centers[0].label = "R"  # ask for the mirror image of the geometry

    with pytest.raises(conf_mod.EmbeddingError):
        conf_mod.generate(mol, log, reference, conf_mod.ConformerOptions(keep=1))
    assert "conformers.stereo_violation" in codes(log)


def test_hydrogens_are_added_and_heavy_atom_order_is_preserved(log):
    mol, _, mol_h, conformers = embed("C[C@H](N)C(=O)O", log, keep=1)
    heavy = mol.GetNumAtoms()
    assert mol_h.GetNumAtoms() > heavy
    for index in range(heavy):
        assert (
            mol_h.GetAtomWithIdx(index).GetSymbol()
            == mol.GetAtomWithIdx(index).GetSymbol()
        )
    assert conformers[0].n_atoms == mol_h.GetNumAtoms()


def test_uff_fallback_when_mmff_has_no_parameters(log):
    """Boron is outside MMFF94; the run must degrade instead of failing."""
    _, _, _, conformers = embed("B(O)(O)c1ccccc1", log, keep=1)
    assert conformers
    if conformers[0].force_field == "UFF":
        assert "conformers.no_mmff_params" in codes(log)


def test_geometry_is_reproducible_for_a_fixed_seed(log):
    from m2i.types import IssueLog

    _, _, _, first = embed("CC(C)CC(C)CO", log, keep=1, seed=1234)
    _, _, _, second = embed("CC(C)CC(C)CO", IssueLog(), keep=1, seed=1234)
    for a, b in zip(first[0].coords, second[0].coords):
        assert a == pytest.approx(list(b), abs=1e-6)


def test_no_force_field_still_produces_a_geometry(log):
    _, _, _, conformers = embed("CCO", log, keep=1, force_field="none")
    assert conformers[0].force_field == "none"
    assert conformers[0].energy is None
    assert all(math.isfinite(v) for xyz in conformers[0].coords for v in xyz)
