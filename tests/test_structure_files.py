"""Drawing files -- ChemDraw and molfiles -- where no vision model is needed."""

from __future__ import annotations

import pytest
from conftest import codes
from rdkit import Chem
from rdkit.Chem import rdDepictor

from m2i import pipeline
from m2i.recognition.base import BackendError
from m2i.recognition.manual import from_molfile, from_structure_file

TARGET = "C[C@H](O)/C=C/c1ccc(Cl)cc1"


def depicted(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    rdDepictor.Compute2DCoords(mol)
    return mol


def write_cdxml(path, mol):
    path.write_text(Chem.MolToCDXMLBlock(mol), encoding="utf-8")
    return path


def read(path, log):
    return pipeline.prepare_molecule(from_structure_file(path).recognize(None), log)


def test_cdxml_keeps_the_stereochemistry_as_drawn(tmp_path, log):
    path = write_cdxml(tmp_path / "drawing.cdxml", depicted(TARGET))
    molecule = read(path, log)
    assert molecule.inchikey == Chem.MolToInchiKey(Chem.MolFromSmiles(TARGET))
    assert molecule.stereo.describe() == "C2: S; 4=5: E"


def test_cdxml_is_read_through_its_wedges(tmp_path, log):
    """The file becomes a molblock, so the stereochemistry is measured from the
    drawing -- the same rule that applies to a MolScribe reading."""
    path = write_cdxml(tmp_path / "drawing.cdxml", depicted(TARGET))
    result = from_structure_file(path).recognize(None)
    assert result.molblock and not result.smiles
    read(path, log)
    assert "stereo.source.molblock" in codes(log)


def test_a_page_with_several_molecules_keeps_the_largest_and_says_so(tmp_path, log):
    page = Chem.CombineMols(depicted(TARGET), depicted("[Na+].[Cl-]"))
    molecule = read(write_cdxml(tmp_path / "page.cdxml", page), log)
    assert molecule.formula == "C10H11ClO"
    dropped = next(i for i in log if i.code == "fragments.dropped")
    assert "[Na+]" in dropped.message and "[Cl-]" in dropped.message


def test_a_page_with_no_molecule_is_an_error_not_an_empty_input(tmp_path):
    path = tmp_path / "scheme.cdxml"
    path.write_text(
        '<?xml version="1.0" encoding="UTF-8" ?>\n'
        '<CDXML><page id="1"><t id="2" p="10 10"><s>Scheme 1</s></t></page></CDXML>\n',
        encoding="utf-8",
    )
    with pytest.raises(BackendError) as excinfo:
        from_structure_file(path)
    assert "no structure found" in str(excinfo.value)


def test_a_molfile_still_works(tmp_path, log):
    mol = depicted(TARGET)
    path = tmp_path / "drawing.mol"
    path.write_text(Chem.MolToMolBlock(mol), encoding="utf-8")
    assert read(path, log).stereo.describe() == "C2: S; 4=5: E"


def test_only_the_first_record_of_an_sdf_is_used(tmp_path, log):
    path = tmp_path / "two.sdf"
    writer = Chem.SDWriter(str(path))
    writer.write(depicted(TARGET))
    writer.write(depicted("CCO"))
    writer.close()
    assert read(path, log).formula == "C10H11ClO"


def test_an_unsupported_file_is_refused(tmp_path):
    path = tmp_path / "drawing.skc"
    path.write_text("whatever")
    with pytest.raises(BackendError) as excinfo:
        from_structure_file(path)
    assert ".cdxml" in str(excinfo.value)


def test_the_old_name_still_works():
    assert from_molfile is from_structure_file


# -- what real ChemDraw pages throw at the reader ------------------------


def draw_page(path, smiles):
    """A ChemDraw page drawn from an unsanitised SMILES, so that bonds written
    ':' outside a ring come out as ChemDraw's delocalised, order-1.5 bonds."""
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol.UpdatePropertyCache(strict=False)
    rdDepictor.Compute2DCoords(mol)
    path.write_text(Chem.MolToCDXMLBlock(mol), encoding="utf-8")
    return path


def test_a_folder_with_an_accent_is_not_a_missing_file(tmp_path, log):
    """RDKit's file readers use narrow paths on Windows: 'Artículo' made the
    file 'not exist'. The reader hands RDKit bytes instead."""
    folder = tmp_path / "Artículo"
    folder.mkdir()
    path = write_cdxml(folder / "dibujo.cdxml", depicted(TARGET))
    assert read(path, log).stereo.describe() == "C2: S; 4=5: E"


def test_a_delocalised_carboxylate_is_read_as_the_anion_and_said_so(tmp_path, log):
    path = draw_page(tmp_path / "glu.cdxml", "CCC(:O):O")
    molecule = read(path, log)
    assert molecule.smiles == "CCC(=O)[O-]"
    assert molecule.charge == -1

    warning = next(i for i in log if i.code == "chemdraw.delocalised_anion")
    assert warning.level == "warning"
    assert "-1" in warning.message and "protonated" in warning.message


def test_a_delocalised_guanidinium_keeps_the_charge_drawn(tmp_path, log):
    """The double bond goes to the N drawn with the '+', and nothing is added."""
    path = draw_page(tmp_path / "arg.cdxml", "CNC(:[NH2]):[NH2+]")
    molecule = read(path, log)
    assert molecule.smiles == "CNC(N)=[NH2+]"
    assert molecule.charge == 1
    assert "chemdraw.delocalised_anion" not in codes(log)


def test_the_net_charge_of_a_page_is_stated(tmp_path, log):
    path = draw_page(tmp_path / "site.cdxml", "CCC(:O):O.CCC(:O):O.CNC(:[NH2]):[NH2+]")
    molecule = pipeline.prepare_molecule(
        from_structure_file(path).recognize(None), log, keep_all_fragments=True
    )
    assert molecule.charge == -1
    warning = next(i for i in log if i.code == "chemdraw.delocalised_anion")
    assert "2 group(s)" in warning.message
    assert "net charge of -1" in warning.message


def test_one_bad_molecule_does_not_take_the_page_with_it(tmp_path, log):
    """Reading sanitised made RDKit drop failing molecules without a word; on a
    real active-site figure that was every residue. Now they are named."""
    path = draw_page(tmp_path / "page.cdxml", "CCO.CC(C)(C)(C)(C)C")
    molecule = read(path, log)
    assert molecule.formula == "C2H6O"
    warning = next(i for i in log if i.code == "chemdraw.molecule_unreadable")
    assert "CC(C)(C)(C)(C)C" in warning.message


def test_a_page_where_nothing_can_be_read_says_why(tmp_path):
    path = draw_page(tmp_path / "bad.cdxml", "CC(C)(C)(C)(C)C")
    with pytest.raises(BackendError) as excinfo:
        from_structure_file(path)
    assert "valence" in str(excinfo.value)


def test_several_molecules_kept_together_come_with_a_3d_caveat(tmp_path, log):
    """A 2D page says nothing about where the molecules sit relative to each
    other, which matters for a cluster model of an active site."""
    page = Chem.CombineMols(depicted(TARGET), depicted("O"))
    path = write_cdxml(tmp_path / "cluster.cdxml", page)
    pipeline.prepare_molecule(
        from_structure_file(path).recognize(None), log, keep_all_fragments=True
    )
    kept = next(i for i in log if i.code == "fragments.kept")
    assert "arbitrary" in kept.message and "crystal structure" in kept.message


@pytest.mark.skipif(
    not __import__("os").environ.get("M2I_TEST_CDX"),
    reason="set M2I_TEST_CDX to a real binary .cdx file to run this check",
)
def test_a_real_binary_cdx_file(log):
    """Binary .cdx cannot be generated from Python, and real drawings are
    usually unpublished work that must not live in the repository, so this
    check runs against a file you point it at."""
    from pathlib import Path
    import os

    path = Path(os.environ["M2I_TEST_CDX"])
    molecule = pipeline.prepare_molecule(
        from_structure_file(path).recognize(None), log, keep_all_fragments=True
    )
    assert molecule.mol.GetNumAtoms() > 0
    assert "chemdraw.molecule_unreadable" not in codes(log)
