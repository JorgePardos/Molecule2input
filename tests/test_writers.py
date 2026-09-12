"""Input-file layout. Coordinates are fixed here so the assertions are about
the file format, not about the conformer generator."""

from __future__ import annotations

import pytest
from conftest import codes

from m2i import config
from m2i.types import IssueLog, JobSpec
from m2i.writers import WriterError, get_writer

WATER = (
    ["O", "H", "H"],
    [(0.0, 0.0, 0.117), (0.0, 0.757, -0.469), (0.0, -0.757, -0.469)],
)


def make_job(profile: config.JobProfile, charge: int = 0, multiplicity: int = 1):
    elements, coords = WATER
    return JobSpec(
        name="water",
        title="O | XLYOFNOQVPJJNP-UHFFFAOYSA-N | m2i test",
        elements=elements,
        coords=coords,
        charge=charge,
        multiplicity=multiplicity,
        profile=profile,
        metadata={"smiles": "O"},
    )


def render(profile_name: str, **overrides) -> tuple[str, IssueLog]:
    profile = config.load_profile(profile_name).merged_with(**overrides)
    log = IssueLog()
    job = make_job(profile)
    return get_writer(profile.program).render(job, log), log


# -- Gaussian ------------------------------------------------------------


def test_gaussian_block_order():
    text, _ = render("gaussian_opt_freq")
    lines = text.split("\n")

    assert lines[0] == "%chk=water.chk"
    assert lines[1] == "%mem=8GB"
    assert lines[2] == "%nprocshared=8"
    assert lines[3].startswith("#p ")
    assert lines[4] == ""  # blank before the title
    assert lines[5] == "O | XLYOFNOQVPJJNP-UHFFFAOYSA-N | m2i test"
    assert lines[6] == ""  # blank after the title
    assert lines[7] == "0 1"
    assert lines[8].split() == ["O", "0.00000000", "0.00000000", "0.11700000"]
    assert len([l for l in lines[8:11] if l.strip()]) == 3


def test_gaussian_ends_with_a_blank_line():
    """Gaussian truncates the geometry without the terminating blank line."""
    text, _ = render("gaussian_opt_freq")
    assert text.endswith("\n\n")


def test_gaussian_route_contents():
    text, _ = render("gaussian_opt_freq_solvent")
    route = text.split("\n")[3]
    assert route.startswith("#p ")
    assert "opt" in route and "freq" in route
    assert "b3lyp/6-31G(d,p)" in route
    assert "empiricaldispersion=gd3bj" in route
    assert "scrf=(smd,solvent=chloroform)" in route


def test_gaussian_charge_and_multiplicity_line():
    profile = config.load_profile("gaussian_opt_freq")
    log = IssueLog()
    job = make_job(profile, charge=-2, multiplicity=3)
    text = get_writer("gaussian").render(job, log)
    assert "\n-2 3\n" in text


def test_gaussian_link1_chains_and_reads_the_checkpoint():
    profile = config.load_profile("gaussian_opt_freq").merged_with(
        link_jobs=({"jobs": ["nmr"], "method": "b3lyp", "basis": "6-311+G(2d,p)"},)
    )
    text = get_writer("gaussian").render(make_job(profile), IssueLog())

    assert text.count("--Link1--") == 1
    second = text.split("--Link1--")[1]
    assert "geom=check" in second and "guess=read" in second
    assert "nmr" in second
    assert "6-311+G(2d,p)" in second
    # the follow-up step must not repeat the geometry
    assert " O " not in second


def test_gen_basis_without_a_basis_block_is_an_error():
    _, log = render("gaussian_opt_freq", basis="gen")
    assert "gaussian.gen_without_block" in codes(log)
    assert log.has_errors()


def test_gen_basis_with_a_block_is_accepted():
    _, log = render(
        "gaussian_opt_freq",
        basis="gen",
        extra_sections=("O H 0\n6-31G(d)\n****",),
    )
    assert not log.has_errors()


def test_a_heavy_element_outside_the_basis_gets_an_ecp():
    """6-31G(d) stops at Kr, and Gaussian would stop at the W. It gets LANL2DZ."""
    profile = config.load_profile("gaussian_opt_freq")
    job = JobSpec(
        name="wcl",
        title="t",
        elements=["W", "Cl"],
        coords=[(0.0, 0.0, 0.0), (0.0, 0.0, 2.3)],
        charge=0,
        multiplicity=1,
        profile=profile,
    )
    log = IssueLog()
    text = get_writer("gaussian").render(job, log)
    assert "b3lyp/genecp" in text
    tail = text[text.index("W   ") :]
    assert "Cl 0\n6-31G(d)\n****\nW 0\nLANL2DZ\n****\n\nW 0\nLANL2DZ\n" in tail
    assert text.endswith("LANL2DZ\n\n")  # Gaussian needs the blank line after the ECP
    assert "basis.auto_ecp" in codes(log)


def test_def2_basis_covers_heavy_elements_without_warning():
    profile = config.load_profile("gaussian_opt_freq").merged_with(basis="def2TZVP")
    job = JobSpec(
        name="wcl",
        title="t",
        elements=["W", "Cl"],
        coords=[(0.0, 0.0, 0.0), (0.0, 0.0, 2.3)],
        charge=0,
        multiplicity=1,
        profile=profile,
    )
    log = IssueLog()
    text = get_writer("gaussian").render(job, log)
    assert "basis.auto_ecp" not in codes(log)
    assert "genecp" not in text


# -- ORCA ----------------------------------------------------------------


def test_orca_layout():
    text, _ = render("orca_opt_freq")
    lines = text.split("\n")

    assert lines[0].startswith("# ")
    assert lines[1].startswith("! ")
    assert "B3LYP" in lines[1] and "def2-SVP" in lines[1]
    assert "D3BJ" in lines[1] and "Opt" in lines[1] and "Freq" in lines[1]
    assert "%pal nprocs 8 end" in text
    assert "%maxcore" in text

    assert "* xyz 0 1" in text
    body = text.split("* xyz 0 1\n")[1]
    assert body.rstrip().endswith("*")
    assert len([l for l in body.split("\n") if l.strip() and l.strip() != "*"]) == 3


def test_orca_cpcm_is_a_keyword_and_smd_is_a_block():
    cpcm, _ = render("orca_opt_freq", solvent={"model": "cpcm", "name": "water"})
    assert "CPCM(water)" in cpcm.split("\n")[1]
    assert "%cpcm" not in cpcm

    smd, _ = render("orca_opt_freq", solvent={"model": "smd", "name": "water"})
    assert "CPCM(water)" not in smd
    assert "smd true" in smd
    assert 'SMDsolvent "water"' in smd


def test_orca_maxcore_is_per_core():
    text, _ = render("orca_opt_freq")
    maxcore = int([l for l in text.split("\n") if l.startswith("%maxcore")][0].split()[1])
    assert 256 <= maxcore <= 8 * 1024


# -- geometry ------------------------------------------------------------


def test_xyz_header_and_atom_count():
    text, _ = render("xyz_only")
    lines = text.split("\n")
    assert lines[0] == "3"
    assert "O" in lines[1]  # comment carries the SMILES/title
    assert len([l for l in lines[2:] if l.strip()]) == 3


def test_unknown_program_is_rejected():
    with pytest.raises(WriterError):
        get_writer("gamess")


def test_writer_forces_the_right_extension(tmp_path):
    profile = config.load_profile("gaussian_opt_freq")
    written = get_writer("gaussian").write(
        make_job(profile), tmp_path / "water.txt", IssueLog()
    )
    assert written.suffix == ".gjf"
    assert written.read_text().startswith("%chk=")


def test_written_file_uses_unix_newlines(tmp_path):
    profile = config.load_profile("gaussian_opt_freq")
    written = get_writer("gaussian").write(
        make_job(profile), tmp_path / "water", IssueLog()
    )
    assert b"\r\n" not in written.read_bytes()
