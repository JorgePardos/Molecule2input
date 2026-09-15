"""End-to-end pipeline behaviour."""

from __future__ import annotations

import json

import numpy as np
import pytest
from conftest import codes
from rdkit import Chem

from m2i import config, pipeline
from m2i.chem.conformers import ConformerOptions
from m2i.types import IssueLog, RecognitionResult


def options(tmp_path, profile="gaussian_opt_freq", **kwargs):
    kwargs.setdefault("conformers", ConformerOptions(keep=1))
    return pipeline.PipelineOptions(
        profile=config.load_profile(profile), output_dir=tmp_path, **kwargs
    )


def manual(smiles: str) -> RecognitionResult:
    return RecognitionResult(smiles=smiles, backend="manual", confidence=1.0)


def test_full_run_writes_input_image_and_provenance(tmp_path, log):
    result = pipeline.run(manual("C[C@H](N)C(=O)O"), options(tmp_path), log)

    suffixes = {p.rsplit(".", 1)[-1] for p in result.written_files}
    assert suffixes == {"gjf", "png", "json"}
    assert result.ok

    record = json.loads((tmp_path / f"{result.molecule.name}.m2i.json").read_text())
    assert record["molecule"]["smiles"] == "C[C@H](N)C(=O)O"
    assert record["molecule"]["stereo"]["centers"][0]["label"] == "S"
    assert record["profile"]["program"] == "gaussian"
    assert record["rdkit_version"]
    assert record["conformers"]


def test_multiple_conformers_get_numbered_files(tmp_path, log):
    result = pipeline.run(
        manual("CCCCCCO"), options(tmp_path, conformers=ConformerOptions(keep=3)), log
    )
    inputs = sorted(p for p in result.written_files if p.endswith(".gjf"))
    assert len(inputs) == 3
    assert inputs[0].endswith("_c01.gjf")
    assert inputs[2].endswith("_c03.gjf")


def test_single_conformer_gets_an_unnumbered_file(tmp_path, log):
    result = pipeline.run(manual("O"), options(tmp_path), log)
    inputs = [p for p in result.written_files if p.endswith(".gjf")]
    assert len(inputs) == 1
    assert "_c01" not in inputs[0]


def test_geometry_is_identical_across_output_formats(tmp_path, log):
    """The .sdf and the .gjf for conformer N must describe the same geometry."""
    recognition = manual("CCCCO")
    base = config.load_profile("gaussian_opt_freq")
    for program in ("gaussian", "sdf"):
        pipeline.run(
            recognition,
            pipeline.PipelineOptions(
                profile=base.merged_with(program=program),
                output_dir=tmp_path,
                name="probe",
                conformers=ConformerOptions(keep=2),
            ),
            IssueLog(),
        )

    for index in (1, 2):
        gjf = (tmp_path / f"probe_c{index:02d}.gjf").read_text().splitlines()
        start = gjf.index("0 1") + 1
        from_gjf = np.array(
            [[float(v) for v in line.split()[1:4]] for line in gjf[start:] if line.strip()]
        )
        mol = Chem.MolFromMolFile(
            str(tmp_path / f"probe_c{index:02d}.sdf"), removeHs=False
        )
        from_sdf = mol.GetConformer().GetPositions()
        assert np.abs(from_gjf - from_sdf).max() < 1e-3


def test_salt_keeps_the_largest_fragment_and_says_so(tmp_path, log):
    molecule = pipeline.prepare_molecule(manual("CC(=O)[O-].[Na+]"), log)
    assert molecule.smiles == "CC(=O)[O-]"
    assert molecule.charge == -1
    assert "fragments.dropped" in codes(log)


def test_keeping_all_fragments_is_possible(tmp_path, log):
    molecule = pipeline.prepare_molecule(
        manual("CC(=O)[O-].[Na+]"), log, keep_all_fragments=True
    )
    assert molecule.charge == 0
    assert "fragments.kept" in codes(log)


def test_name_is_made_filesystem_safe(tmp_path, log):
    molecule = pipeline.prepare_molecule(manual("CCO"), log, name="my molecule/2 <x>")
    assert "/" not in molecule.name and " " not in molecule.name
    assert molecule.name.startswith("my_molecule")


def test_default_name_is_the_inchikey_skeleton(tmp_path, log):
    molecule = pipeline.prepare_molecule(manual("CCO"), log)
    assert molecule.name == molecule.inchikey.split("-")[0]


def test_title_is_a_single_line(tmp_path, log):
    """A blank line in the title section would truncate a Gaussian input."""
    result = pipeline.run(manual(r"C/C=C\C"), options(tmp_path), log)
    text = next(p for p in result.written_files if p.endswith(".gjf"))
    lines = open(text).read().split("\n")
    assert lines[4] == "" and lines[6] == ""
    assert lines[5].strip()


def test_invalid_smiles_raises_with_a_useful_message(log):
    from m2i.chem.sanitize import MoleculeParseError

    with pytest.raises(MoleculeParseError) as excinfo:
        pipeline.prepare_molecule(manual("C1CC"), log)
    assert "SMILES" in str(excinfo.value)


def test_charge_override_reaches_the_input_file(tmp_path, log):
    result = pipeline.run(
        manual("CC(=O)O"), options(tmp_path, charge=-1, multiplicity=1), log
    )
    text = open(next(p for p in result.written_files if p.endswith(".gjf"))).read()
    assert "\n-1 1\n" in text


# -- geometry computed once, written in any format -----------------------


def test_every_format_carries_the_same_geometry(tmp_path, log):
    """The GUI embeds once and writes whichever format is asked for; the
    files must describe identical coordinates."""
    from m2i import config
    from m2i.chem.conformers import ConformerOptions

    molecule = pipeline.prepare_molecule(
        RecognitionResult(smiles="C[C@H](O)/C=C/c1ccc(Cl)cc1", backend="manual"), log
    )
    embedding = pipeline.embed(molecule, ConformerOptions(keep=1), log)

    written = {}
    for program in ("gaussian", "orca", "xyz"):
        profile = config.JobProfile(program=program, method="b3lyp", basis="def2-SVP")
        options = pipeline.PipelineOptions(
            profile=profile,
            output_dir=tmp_path / program,
            write_comparison=False,
            write_provenance=False,
        )
        result = pipeline.write_inputs(molecule, embedding, options, log)
        written[program] = _cartesian_block(result.written_files[0])

    assert written["gaussian"] == written["orca"] == written["xyz"]
    assert len(written["xyz"]) == 23


def test_generate_inputs_is_embed_then_write(tmp_path, log):
    from m2i import config

    molecule = pipeline.prepare_molecule(
        RecognitionResult(smiles="CCO", backend="manual"), log
    )
    options = pipeline.PipelineOptions(
        profile=config.JobProfile(program="xyz", method="", basis="", jobs=()),
        output_dir=tmp_path,
    )
    result = pipeline.generate_inputs(molecule, options, log)
    assert result.conformers
    assert any(path.endswith(".xyz") for path in result.written_files)


def _cartesian_block(path) -> list[tuple]:
    import re
    from pathlib import Path

    rows = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) == 4 and re.fullmatch(r"[A-Z][a-z]?", parts[0]):
            try:
                rows.append((parts[0], *(round(float(v), 5) for v in parts[1:])))
            except ValueError:
                continue
    return rows
