"""Stereochemistry perception: CIP labels and undefined-centre detection."""

from __future__ import annotations

import pytest
from conftest import codes

from m2i.chem import sanitize
from m2i.chem import stereo as stereo_mod


def summarize(smiles: str):
    mol = sanitize.mol_from_smiles(smiles)
    stereo_mod.assign(mol)
    return mol, stereo_mod.summarize(mol)


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("F[C@@H](Cl)Br", "S"),
        ("F[C@H](Cl)Br", "R"),
        ("C[C@@H](N)C(=O)O", "R"),  # D-alanine as written (N before the carboxyl)
        ("C[C@H](N)C(=O)O", "S"),  # L-alanine
    ],
)
def test_tetrahedral_cip_labels(smiles, expected):
    _, summary = summarize(smiles)
    assert len(summary.centers) == 1
    assert summary.centers[0].specified
    assert summary.centers[0].label == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        (r"C/C=C/C", "E"),
        (r"C/C=C\C", "Z"),
        (r"C/C=C/Cl", "E"),
    ],
)
def test_double_bond_cip_labels(smiles, expected):
    _, summary = summarize(smiles)
    assert len(summary.bonds) == 1
    assert summary.bonds[0].specified
    assert summary.bonds[0].label == expected


def test_undefined_centre_is_reported_not_guessed(log):
    _, summary = summarize("CC(F)Cl")
    assert len(summary.centers) == 1
    assert not summary.centers[0].specified
    assert summary.centers[0].label is None
    assert not summary.is_complete

    stereo_mod.report(summary, log)
    assert "stereo.unspecified_center" in codes(log)


def test_undefined_double_bond_is_reported(log):
    _, summary = summarize("CC=CC")
    assert len(summary.bonds) == 1
    assert not summary.bonds[0].specified

    stereo_mod.report(summary, log)
    assert "stereo.unspecified_bond" in codes(log)


def test_mixed_specified_and_unspecified():
    _, summary = summarize("C[C@H](O)/C=C/C(F)Cl")
    assert len(summary.centers) == 2
    assert len(summary.bonds) == 1
    assert len(summary.unspecified_centers) == 1
    assert summary.unspecified_bonds == []
    assert not summary.is_complete


def test_achiral_molecule_has_nothing_to_lose(log):
    _, summary = summarize("c1ccccc1")
    assert summary.centers == []
    assert summary.bonds == []
    assert summary.is_complete
    stereo_mod.report(summary, log)
    assert not log.warnings


def test_fingerprint_only_covers_specified_elements():
    _, summary = summarize("C[C@H](O)/C=C/C(F)Cl")
    fingerprint = summary.fingerprint()
    # the undefined centre must not appear: it has no configuration to preserve
    assert len(fingerprint) == 2
    assert set(fingerprint.values()) == {"S", "E"}


def test_compare_detects_inversion_and_loss():
    _, summary = summarize("F[C@@H](Cl)Br")
    _, inverted = summarize("F[C@H](Cl)Br")
    _, flat = summarize("FC(Cl)Br")

    assert stereo_mod.compare(summary, summary) == []
    assert stereo_mod.compare(summary, inverted) == ["atom 2: S -> R"]
    assert stereo_mod.compare(summary, flat) == ["atom 2: S -> lost"]


def test_molblock_wedges_drive_the_stereochemistry(log):
    """A molblock with a wedge must produce a defined centre."""
    from m2i.types import RecognitionResult

    reference = sanitize.mol_from_smiles("F[C@@H](Cl)Br")
    from rdkit import Chem
    from rdkit.Chem import rdDepictor

    rdDepictor.Compute2DCoords(reference)
    Chem.WedgeMolBonds(reference, reference.GetConformer())
    molblock = Chem.MolToMolBlock(reference)

    result = RecognitionResult(smiles="", molblock=molblock, backend="test")
    mol = sanitize.mol_from_recognition(result, log)
    stereo_mod.assign(mol)
    summary = stereo_mod.summarize(mol)

    assert summary.centers[0].specified
    assert summary.centers[0].label == "S"
    assert "stereo.source.molblock" in codes(log)


# -- phosphorus that only looks stereogenic ------------------------------

# dTDP-beta-L-rhamnose dianion; stereochemistry checked against PubChem CID 121966.
DTDP_RHAMNOSE = (
    "Cc1cn([C@H]2C[C@H](O)[C@@H](COP(=O)([O-])OP(=O)([O-])O[C@H]3O[C@@H](C)"
    "[C@H](O)[C@@H](O)[C@H]3O)O2)c(=O)[nH]c1=O"
)


def test_phosphodiester_phosphorus_is_not_a_stereocentre(log):
    """The DNA backbone is not chiral at phosphorus: the P=O and the P-O(-) are
    the same position with the charge delocalised over both."""
    mol, summary = summarize("CO[P](=O)([O-])OCC")
    stereo_mod.report(summary, log)

    assert summary.unspecified_centers == []
    assert "stereo.unspecified_center" not in codes(log)
    assert [c.symbol for c in summary.ignored_centers] == ["P"]


def test_a_real_nucleotide_produces_no_phantom_warning(log):
    """dTDP-L-rhamnose: eight genuine carbon centres, two phantom phosphorus
    ones. Before this filter the warning fired on every nucleotide drawn."""
    mol, summary = summarize(DTDP_RHAMNOSE)
    stereo_mod.report(summary, log)

    assert len(summary.centers) == 8
    assert all(c.symbol == "C" and c.specified for c in summary.centers)
    assert [c.symbol for c in summary.ignored_centers] == ["P", "P"]
    assert "stereo.unspecified_center" not in codes(log)
    assert "stereo.equivalent_oxoanion" in codes(log)
    assert "P" not in summary.describe()


def test_phosphorothioate_is_still_flagged(log):
    """Sulfur breaks the equivalence, so Sp/Rp is real and must be reported --
    this is the case the filter must not swallow."""
    mol, summary = summarize("CO[P](=S)([O-])OCC")
    stereo_mod.report(summary, log)

    assert [c.symbol for c in summary.unspecified_centers] == ["P"]
    assert summary.ignored_centers == []
    assert "stereo.unspecified_center" in codes(log)


def test_the_sulfur_can_be_on_either_side():
    """P(=O)(S-) is the same compound written the other way round."""
    _, summary = summarize("CO[P](=O)([S-])OCC")
    assert [c.symbol for c in summary.unspecified_centers] == ["P"]


def test_a_neutral_phosphate_is_filtered_too():
    """P(=O)(OH) differs from P(-O(-)) only by where a mobile proton sits."""
    _, summary = summarize("CO[P](=O)(O)OCC")
    assert summary.unspecified_centers == []
    assert [c.symbol for c in summary.ignored_centers] == ["P"]


def test_a_phosphine_oxide_keeps_its_stereocentre():
    """Three different carbons and one oxygen: nothing to make equivalent."""
    _, summary = summarize("CC[P](=O)(C)c1ccccc1")
    assert [c.symbol for c in summary.unspecified_centers] == ["P"]
    assert summary.ignored_centers == []


def test_a_sulfoxide_keeps_its_stereocentre():
    """Sulfoxides are resolvable and their configuration is real."""
    _, summary = summarize("CS(=O)c1ccccc1")
    assert [c.symbol for c in summary.unspecified_centers] == ["S"]
    assert summary.ignored_centers == []


def test_carbon_is_never_filtered():
    """The rule must not reach beyond P/S-like centres."""
    _, summary = summarize("OC(N)(F)Cl")
    assert summary.ignored_centers == []
    assert [c.symbol for c in summary.unspecified_centers] == ["C"]


def test_ignored_centres_reach_the_provenance_record():
    _, summary = summarize(DTDP_RHAMNOSE)
    assert summary.to_dict()["ignored_centers"] == [
        {"atom_index": 11, "symbol": "P"},
        {"atom_index": 15, "symbol": "P"},
    ]


def test_the_filter_does_not_disturb_the_3d_comparison():
    """A filtered centre has no descriptor, so it cannot report a false
    inversion when the 2D and 3D summaries are compared."""
    _, summary = summarize(DTDP_RHAMNOSE)
    assert not any(key.startswith("atom:11") for key in summary.fingerprint())
    assert stereo_mod.compare(summary, summary) == []
