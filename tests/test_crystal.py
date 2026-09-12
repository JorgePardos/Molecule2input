"""The CIF route: crystal in, experimental geometry and a checked input out.

The three COD structures in tests/data/cod each stress something different;
the README there says what. Small synthetic P1 crystals cover the cases no
real file here happens to show.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from conftest import codes

pytest.importorskip("gemmi", reason="the CIF route needs the [crystal] extra")

from m2i import config
from m2i.chem import metals
from m2i.cli import main as cli
from m2i.crystal import CrystalError, read_cif
from m2i.crystal import jobs as crystal_jobs
from m2i.crystal.coordination import analyse
from m2i.types import IssueLog

COD = Path(__file__).parent / "data" / "cod"
FERROCENE, TRANSPLATIN, COBALT = COD / "7033930.cif", COD / "1538403.cif", COD / "2106540.cif"


def read(path, log=None):
    return read_cif(path, log if log is not None else IssueLog())


def p1_cif(tmp_path, atoms, *, a=10.0, formula=None, z=None, name="test.cif"):
    """A P1 crystal from (label, element, x, y, z[, occupancy]) in fractions."""
    lines = [
        "data_test",
        f"_cell_length_a {a}", "_cell_length_b 10", "_cell_length_c 10",
        "_cell_angle_alpha 90", "_cell_angle_beta 90", "_cell_angle_gamma 90",
        "_symmetry_space_group_name_H-M 'P 1'",
    ]
    if formula:
        lines.append(f"_chemical_formula_sum '{formula}'")
    if z:
        lines.append(f"_cell_formula_units_Z {z}")
    lines += [
        "loop_", "_symmetry_equiv_pos_as_xyz", "'x, y, z'",
        "loop_", "_atom_site_label", "_atom_site_type_symbol",
        "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z",
        "_atom_site_occupancy",
    ]
    for atom in atoms:
        label, element, x, y, zz, *occ = atom
        lines.append(f"{label} {element} {x} {y} {zz} {occ[0] if occ else 1}")
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


FORMALDEHYDE = [
    ("C1", "C", 0.5, 0.5, 0.5),
    ("O1", "O", 0.621, 0.5, 0.5),
    ("H1", "H", 0.445, 0.5953, 0.5),
    ("H2", "H", 0.445, 0.4047, 0.5),
]


# -- rebuilding molecules from the crystal --------------------------------


def test_half_a_ferrocene_in_the_file_becomes_a_whole_one():
    reading = read(FERROCENE)
    assert [s.formula for s in reading.species] == ["C10H10Fe"]
    ferrocene = reading.species[0]
    assert ferrocene.n_atoms == 21
    assert ferrocene.copies == 2  # Z = 2


def test_the_major_disorder_component_is_kept_even_when_it_is_group_2():
    log = IssueLog()
    reading = read(FERROCENE, log)
    assert "group 2 kept (occupancy 0.69)" in reading.disorder
    assert "crystal.disorder" in codes(log)


def test_three_atoms_on_an_inversion_centre_give_trans_platin():
    reading = read(TRANSPLATIN)
    platin = reading.species[0]
    assert platin.formula == "Cl2N2Pt"
    assert platin.n_atoms == 5  # Pt, Cl, N in the file; Cl and N doubled by symmetry


def test_independent_copies_and_counter_ions_are_counted():
    """Four symmetry-independent Co in the file, twelve complexes per cell."""
    reading = read(COBALT)
    assert [(s.formula, s.copies) for s in reading.species] == [("CoN6", 12), ("Cl", 36)]


def test_an_extended_network_is_named_not_cut(tmp_path):
    chain = p1_cif(
        tmp_path,
        [("Cu1", "Cu", 0.0, 0.0, 0.0), ("Cl1", "Cl", 0.5, 0.0, 0.0)] + FORMALDEHYDE,
        a=4.6,
    )
    log = IssueLog()
    reading = read(chain, log)
    assert reading.extended == ["ClCu"]
    assert [s.formula for s in reading.species] == ["CH2O"]
    assert "crystal.extended" in codes(log)


def test_a_crystal_that_is_only_a_network_has_nothing_to_calculate(tmp_path):
    chain = p1_cif(tmp_path, [("Cu1", "Cu", 0.0, 0.0, 0.0), ("Cl1", "Cl", 0.5, 0.0, 0.0)], a=4.6)
    with pytest.raises(CrystalError) as excinfo:
        read(chain)
    assert "extended network" in str(excinfo.value)


def test_overlapping_atoms_block_the_input(tmp_path):
    """Two disorder positions both kept would put atoms 0.45 Å apart."""
    path = p1_cif(tmp_path, FORMALDEHYDE + [("C2", "C", 0.545, 0.5, 0.5, 0.5)])
    log = IssueLog()
    reading = read(path, log)
    assert "crystal.overlap" in codes(log)
    assert "crystal.partial_occupancy" in codes(log)
    with pytest.raises(CrystalError):
        crystal_jobs.write(
            reading, 0, charge=0, multiplicity=1, output_dir=tmp_path,
            profile=config.load_profile("xyz_only"), log=log,
        )


# -- checks on what was rebuilt --------------------------------------------


def test_the_rebuilt_crystal_matches_the_declared_formula():
    log = IssueLog()
    read(FERROCENE, log)
    assert "crystal.formula_ok" in codes(log)


def test_hydrogens_never_located_are_caught():
    log = IssueLog()
    reading = read(TRANSPLATIN, log)
    assert reading.missing_hydrogens == 6
    # A warning: m2i proposes the hydrogens itself, and holds the input back
    # only if the species being written is the one short of them.
    assert next(i for i in log if i.code == "crystal.missing_hydrogens").level == "warning"


def test_a_formula_mismatch_other_than_hydrogen_is_reported(tmp_path):
    path = p1_cif(tmp_path, FORMALDEHYDE, formula="C2 H2 O", z=1)
    log = IssueLog()
    read(path, log)
    assert "crystal.formula_mismatch" in codes(log)


def test_xray_hydrogens_are_moved_to_standard_lengths():
    reading = read(FERROCENE)
    ferrocene = reading.species[0]
    lengths = [
        np.linalg.norm(ferrocene.coords[h] - ferrocene.coords[ferrocene.neighbours(h)[0]])
        for h, e in enumerate(ferrocene.elements) if e == "H"
    ]
    assert reading.hydrogens_normalised == 10
    assert np.allclose(lengths, 1.089, atol=1e-3)


def test_xray_hydrogens_can_be_kept_as_measured():
    reading = read_cif(FERROCENE, IssueLog(), normalise_hydrogens=False)
    assert reading.hydrogens_normalised == 0


def test_a_folder_with_an_accent_is_fine(tmp_path):
    folder = tmp_path / "Artículo"
    folder.mkdir()
    path = folder / "ferroceno.cif"
    path.write_bytes(FERROCENE.read_bytes())
    assert read(path).species[0].formula == "C10H10Fe"


def test_a_file_that_is_not_a_cif_is_refused(tmp_path):
    path = tmp_path / "nothing.cif"
    path.write_text("data_empty\n_cell_length_a 10\n")
    with pytest.raises(CrystalError):
        read(path)


# -- the coordination sphere ----------------------------------------------


def test_ferrocene_is_a_sandwich_of_two_eta5_rings():
    (iron,) = analyse(read(FERROCENE).species[0])
    assert iron.geometry == "sandwich"
    assert [d.hapticity for d in iron.donors] == [5, 5]
    assert all(abs(d.distance - 1.654) < 0.01 for d in iron.donors)  # Fe-centroid


def test_transplatin_is_square_planar_and_trans():
    (platinum,) = analyse(read(TRANSPLATIN).species[0])
    assert platinum.geometry == "square planar"
    assert platinum.isomers == ["trans (Cl of Cl x2)", "trans (N of N x2)"]


def test_the_cobalt_complex_is_octahedral():
    (cobalt,) = analyse(read(COBALT).species[0])
    assert cobalt.geometry == "octahedral"
    assert len(cobalt.trans_pairs) == 3
    assert all(1.95 < d.distance < 1.98 for d in cobalt.donors)


def test_cis_is_told_from_trans(tmp_path):
    """A synthetic cis-[PtCl2(NH3)2] (heavy atoms only), to check the other answer."""
    atoms = [
        ("Pt1", "Pt", 0.5, 0.5, 0.5),
        ("Cl1", "Cl", 0.732, 0.5, 0.5), ("Cl2", "Cl", 0.5, 0.732, 0.5),
        ("N1", "N", 0.295, 0.5, 0.5), ("N2", "N", 0.5, 0.295, 0.5),
    ]
    (platinum,) = analyse(read(p1_cif(tmp_path, atoms)).species[0])
    assert platinum.geometry == "square planar"
    assert platinum.isomers == ["cis (Cl of Cl x2)", "cis (N of N x2)"]


# -- what to ask, and what to check ---------------------------------------


@pytest.mark.parametrize("element, state, expected", [
    ("Fe", 2, ("d", 6)), ("Co", 3, ("d", 6)), ("Pt", 2, ("d", 8)),
    ("Cu", 2, ("d", 9)), ("Zn", 2, ("d", 10)), ("Gd", 3, ("f", 7)), ("U", 4, ("f", 2)),
])
def test_electron_counts(element, state, expected):
    assert metals.valence_electron_count(element, state) == expected


def test_spin_options_for_d5():
    options = metals.spin_options("Fe", 3)
    assert [(o.multiplicity, o.label) for o in options] == [
        (2, "low spin"), (4, "intermediate spin"), (6, "high spin"),
    ]


def test_d10_has_one_option():
    assert [o.multiplicity for o in metals.spin_options("Zn", 2)] == [1]


def test_an_impossible_multiplicity_is_explained():
    problem = metals.parity_problem(["Fe", "C", "C"], 0, 2)  # 38 electrons
    assert "even" in problem and "odd" in problem
    assert metals.parity_problem(["Fe", "C", "C"], 0, 1) is None


def test_the_charge_is_suggested_from_the_counter_ions():
    charge, reason = metals.suggest_charge("CoN6", 12, [("Cl", 36)])
    assert charge == 3
    assert "36 x Cl" in reason


def test_no_suggestion_when_a_species_is_unknown():
    charge, reason = metals.suggest_charge("CoN6", 1, [("C6H5NO2", 1)])
    assert charge is None and "C6H5NO2" in reason


def test_solvent_does_not_count_towards_the_charge():
    assert metals.suggest_charge("C10H10Fe", 2, [("H2O", 4)])[0] == 0


def test_default_multiplicity_is_the_lowest_the_metal_allows():
    ferrocene = read(FERROCENE).species[0]
    assert crystal_jobs.default_multiplicity(ferrocene, 0, {"Fe": 2}) == 1
    assert crystal_jobs.default_multiplicity(ferrocene, 1, {"Fe": 3}) == 2


# -- writing -----------------------------------------------------------------


def test_the_gaussian_input_keeps_the_crystal_geometry(tmp_path):
    reading = read(FERROCENE)
    written = crystal_jobs.write(
        reading, 0, charge=0, multiplicity=1, oxidation={"Fe": 2}, output_dir=tmp_path,
        profile=config.load_profile("gaussian_opt_freq"), log=IssueLog(),
    )
    gjf = Path(written[0]).read_text(encoding="utf-8")
    assert "\n0 1\n" in gjf
    rows = [line.split() for line in gjf.splitlines() if len(line.split()) == 4]
    elements = [r[0] for r in rows]
    coords = np.array([[float(v) for v in r[1:]] for r in rows])
    assert sorted(elements) == sorted(reading.species[0].elements)
    iron = coords[elements.index("Fe")]
    assert np.allclose(iron, 0.0)  # centred on the metal
    fe_c = [np.linalg.norm(c - iron) for c, e in zip(coords, elements) if e == "C"]
    assert 2.0 < min(fe_c) and max(fe_c) < 2.1


def test_the_provenance_records_the_choices(tmp_path):
    import json

    written = crystal_jobs.write(
        read(FERROCENE), 0, charge=0, multiplicity=1, oxidation={"Fe": 2},
        output_dir=tmp_path, profile=config.load_profile("orca_opt_freq"), log=IssueLog(),
    )
    record = json.loads(Path(written[1]).read_text(encoding="utf-8"))
    assert record["source"]["cod_id"] == "7033930"
    assert record["electronic_state"] == {
        "charge": 0, "multiplicity": 1, "oxidation_states": {"Fe": 2}, "decided_by": "user",
    }
    assert "group 2 kept" in record["reconstruction"]["disorder"]
    assert any("sandwich" in line for line in record["species"]["coordination"])


def test_missing_hydrogens_stop_the_input_unless_allowed(tmp_path):
    reading = read(TRANSPLATIN)
    kwargs = dict(charge=0, multiplicity=1, output_dir=tmp_path, log=IssueLog(),
                  profile=config.load_profile("gaussian_opt_freq"))
    with pytest.raises(CrystalError) as excinfo:
        crystal_jobs.write(reading, 0, **kwargs)
    assert "missing" in str(excinfo.value)
    assert crystal_jobs.write(reading, 0, allow_missing_hydrogens=True, **kwargs)


def test_an_impossible_state_stops_the_input(tmp_path):
    with pytest.raises(CrystalError) as excinfo:
        crystal_jobs.write(
            read(FERROCENE), 0, charge=0, multiplicity=2, output_dir=tmp_path,
            profile=config.load_profile("gaussian_opt_freq"), log=IssueLog(),
        )
    assert "impossible" in str(excinfo.value)


def test_platinum_gets_an_ecp_in_orca_too(tmp_path):
    profile = config.load_profile("orca_opt_freq").merged_with(basis="6-31G(d)")
    written = crystal_jobs.write(
        read(TRANSPLATIN), 0, charge=0, multiplicity=1, profile=profile,
        output_dir=tmp_path, log=IssueLog(), allow_missing_hydrogens=True,
    )
    text = Path(written[0]).read_text(encoding="utf-8")
    assert 'NewGTO Pt "def2-TZVP" end' in text and 'NewECP Pt "def2-ECP" end' in text


# -- the command line ---------------------------------------------------------


def test_cli_asks_for_the_charge_of_a_metal_complex(tmp_path, capsys):
    assert cli(["from-cif", str(FERROCENE), "-o", str(tmp_path)]) == 2
    assert "pass --charge (suggested: +0)" in capsys.readouterr().err


def test_cli_asks_for_the_spin_of_a_metal_complex(tmp_path, capsys):
    assert cli(["from-cif", str(FERROCENE), "-o", str(tmp_path), "--charge", "0"]) == 2
    assert "spin state" in capsys.readouterr().err


def test_cli_writes_when_told(tmp_path, capsys):
    code = cli([
        "from-cif", str(FERROCENE), "-o", str(tmp_path),
        "--charge", "0", "--mult", "1", "--oxidation", "Fe=2",
    ])
    out = capsys.readouterr().out
    assert code == 0
    assert "Fe(+2) is d6" in out and "sandwich" in out
    assert (tmp_path / "C10H10Fe_7033930.gjf").is_file()


def test_cli_needs_no_questions_for_an_organic_crystal(tmp_path):
    path = p1_cif(tmp_path, FORMALDEHYDE, formula="C H2 O", z=1)
    assert cli(["from-cif", str(path), "-o", str(tmp_path), "--program", "xyz"]) == 0
    assert (tmp_path / "CH2O_test.xyz").is_file()
