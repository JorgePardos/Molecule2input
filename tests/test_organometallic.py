"""Metal complexes from drawings: isomer counting, reading the page, building in 3D.

Every drawing here is synthetic, written as a molfile by the helpers below: a
metal at the origin, donors placed on the page, wedges and hashes on the bonds
that leave the plane -- what a chemist would draw in ChemDraw.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from m2i.cli import main
from m2i.crystal.coordination import analyse
from m2i.crystal.structure import Species
from m2i.organometallic import (
    arrangements,
    build,
    groups,
    is_chiral,
    mirror_partners,
    read,
    read_raw,
)
from m2i.organometallic import jobs as om_jobs
from m2i.organometallic.drawing import prepare
from m2i.organometallic.geometry import POLYHEDRA
from m2i.recognition.base import BackendError

WEDGE, HASH = 1, 6
CHARGE_CODE = {1: 3, 2: 2, 3: 1, -1: 5, -2: 6, -3: 7}


def molblock(atoms, bonds) -> str:
    """V2000 text. atoms: (symbol, x, y[, charge]); bonds: (i, j, order[, stereo])."""
    lines = ["", "  m2i-test          2D", ""]
    lines.append(f"{len(atoms):3d}{len(bonds):3d}  0  0  0  0  0  0  0  0999 V2000")
    for symbol, x, y, *_ in atoms:
        lines.append(f"{x:10.4f}{y:10.4f}{0.0:10.4f} {symbol:<3} 0  0  0  0  0  0  0  0  0  0  0  0")
    for i, j, order, *rest in bonds:
        lines.append(f"{i + 1:3d}{j + 1:3d}{order:3d}{(rest[0] if rest else 0):3d}")
    charged = [(k + 1, a[3]) for k, a in enumerate(atoms) if len(a) > 3 and a[3]]
    if charged:
        lines.append("M  CHG" + f"{len(charged):3d}" + "".join(f"{k:4d}{q:4d}" for k, q in charged))
    lines.append("M  END")
    return "\n".join(lines) + "\n"


class Page:
    """A drawing being assembled: the metal at the origin, ligands around it."""

    def __init__(self, metal: str, charge: int = 0):
        self.atoms = [(metal, 0.0, 0.0, charge)]
        self.bonds = []

    def add(self, symbol, x, y, charge=0):
        self.atoms.append((symbol, float(x), float(y), charge))
        return len(self.atoms) - 1

    def bond(self, i, j, order=1, stereo=0):
        self.bonds.append((i, j, order, stereo))

    def to_metal(self, donor, depth=0):
        """A metal-donor bond drawn from the metal: +1 wedge, -1 hash, 0 plain."""
        self.bond(0, donor, 1, {1: WEDGE, -1: HASH, 0: 0}[depth])

    def carbonyl(self, direction, depth=0):
        u = np.asarray(direction, float) / np.linalg.norm(direction)
        c = self.add("C", *(1.5 * u))
        o = self.add("O", *(2.7 * u))
        self.bond(c, o, 2)
        self.to_metal(c, depth)
        return c

    def phosphine(self, direction, depth=0):
        u = np.asarray(direction, float) / np.linalg.norm(direction)
        n = np.array([-u[1], u[0]])
        p = self.add("P", *(1.8 * u))
        for offset in (-0.8, 0.0, 0.8):
            c = self.add("C", *(2.8 * u + offset * n))
            self.bond(p, c)
        self.to_metal(p, depth)
        return p

    def simple(self, symbol, direction, depth=0, charge=0):
        u = np.asarray(direction, float) / np.linalg.norm(direction)
        d = self.add(symbol, *(1.8 * u), charge=charge)
        self.to_metal(d, depth)
        return d

    def write(self, path):
        path.write_text(molblock(self.atoms, self.bonds), encoding="utf-8")
        return path


def square_planar(tmp_path, trans: bool):
    page = Page("Pt")
    if trans:
        page.simple("Cl", (1, 0)), page.simple("Cl", (-1, 0))
        page.simple("N", (0, 1)), page.simple("N", (0, -1))
    else:
        page.simple("Cl", (1, 0)), page.simple("Cl", (0, 1))
        page.simple("N", (-1, 0)), page.simple("N", (0, -1))
    return page.write(tmp_path / ("trans.mol" if trans else "cis.mol"))


#: An octahedron as usually drawn: four bonds in the page, one wedge towards
#: the viewer (lower left) and one hash away from it (upper right).
OCTAHEDRON = {
    "up": ((0, 1), 0), "down": ((0, -1), 0), "left": ((-1, 0), 0), "right": ((1, 0), 0),
    "front": ((-0.7, -0.7), 1), "back": ((0.7, 0.7), -1),
}


def mo_co3_pme3_3(tmp_path, isomer: str):
    phosphines = {"mer": ("up", "down", "left"), "fac": ("up", "right", "front")}[isomer]
    page = Page("Mo")
    for site, (direction, depth) in OCTAHEDRON.items():
        if site in phosphines:
            page.phosphine(direction, depth)
        else:
            page.carbonyl(direction, depth)
    return page.write(tmp_path / f"{isomer}.mol")


# -- tris(ethylenediamine)cobalt(III), drawn from a 3D model -------------------

#: Three cis edges of the octahedron, one per chelate, covering all six vertices.
EN_EDGES = [((1, 0, 0), (0, 1, 0)), ((-1, 0, 0), (0, 0, 1)), ((0, -1, 0), (0, 0, -1))]


def _rotation():
    a, b = np.radians(25), np.radians(35)
    rx = np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
    ry = np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0], [-np.sin(b), 0, np.cos(b)]])
    return ry @ rx


def co_en3_model(mirror: bool):
    """Heavy atoms of [Co(en)3]3+ in 3D: Co, then N, C, C, N for each chelate."""
    coords, bonds = [np.zeros(3)], []
    for a, b in EN_EDGES:
        a, b = np.array(a, float), np.array(b, float)
        mid = (a + b) / np.linalg.norm(a + b)
        n1, n2 = 2.0 * a, 2.0 * b
        c1, c2 = n1 + 0.9 * mid + 0.4 * (b - a), n2 + 0.9 * mid - 0.4 * (b - a)
        start = len(coords)
        coords += [n1, c1, c2, n2]
        bonds += [(0, start), (0, start + 3), (start, start + 1), (start + 1, start + 2), (start + 2, start + 3)]
    coords = np.array(coords) @ _rotation().T
    if mirror:
        coords[:, 0] *= -1
    return coords, bonds


def co_en3_drawing(tmp_path, mirror: bool):
    coords, bonds = co_en3_model(mirror)
    page = Page("Co", charge=3)
    for position in coords[1:]:
        page.add("N" if len(page.atoms) % 4 in (0, 1) else "C", position[0], position[1])
    for i, j in bonds:
        if i == 0:
            z = coords[j][2]
            page.to_metal(j, 1 if z > 0.6 else -1 if z < -0.6 else 0)
        else:
            page.bond(i, j)
    return page.write(tmp_path / ("lambda.mol" if mirror else "delta.mol"))


def helicity_of_model(mirror: bool) -> str:
    coords, bonds = co_en3_model(mirror)
    elements = ["Co"] + ["N" if k % 4 in (1, 0) else "C" for k in range(1, len(coords))]
    species = Species(elements=elements, coords=coords, labels=[f"{e}{k + 1}" for k, e in enumerate(elements)],
                      bonds=bonds)
    return analyse(species)[0].helicity


# -- counting isomers ------------------------------------------------------------


def _count(n_donors, labels, chelated=(), polyhedron=None):
    donors = list(range(n_donors))
    drawn = {d: np.array([1.0, 0.0, 0.0]) for d in donors}
    return arrangements(donors, drawn, list(chelated), polyhedron=polyhedron,
                        labels=dict(enumerate(labels)))


@pytest.mark.parametrize(
    "labels, polyhedron, expected",
    [
        ("AAAABB", "octahedral", 2),       # cis / trans
        ("AAABBB", "octahedral", 2),       # fac / mer
        ("ABCDEF", "octahedral", 30),      # 15 pairs of enantiomers
        ("AABB", "square planar", 2),
        ("ABCD", "square planar", 3),
        ("ABCD", "tetrahedral", 2),        # a pair of enantiomers
        ("AABC", "tetrahedral", 1),
    ],
)
def test_isomer_counts(labels, polyhedron, expected):
    assert len(_count(len(labels), labels, polyhedron=polyhedron)) == expected


def test_chelate_forces_cis():
    # M(AA)B4 with the chelate spanning a trans pair is impossible: one isomer.
    options = _count(6, "AABBBB", chelated=[(0, 1)], polyhedron="octahedral")
    assert len(options) == 1
    assert (0, 1) not in {tuple(sorted(p)) for p in options[0].trans_pairs()}


def test_mixed_octahedron_is_chiral_but_trans_one_is_not():
    options = _count(6, "AABBCC", polyhedron="octahedral")
    chiral = [o for o in options if is_chiral(o, dict(enumerate("AABBCC")))]
    assert len(options) == 6 and len(chiral) == 2  # all-cis is the chiral pair


def test_polyhedra_are_unit_vectors():
    for vertices in POLYHEDRA.values():
        assert np.allclose(np.linalg.norm(np.asarray(vertices), axis=1), 1.0)


# -- reading the drawing -----------------------------------------------------------


def test_trans_platin_is_read_and_placed(tmp_path):
    drawing = read(square_planar(tmp_path, trans=True))
    assert drawing.metal_symbol == "Pt" and len(drawing.donors) == 4
    options = arrangements(drawing.donors, drawing.drawn, drawing.chelated, labels=drawing.labels)
    best = options[0]
    assert best.polyhedron == "square planar"
    cl = {d for d in drawing.donors if drawing.mol.GetAtomWithIdx(d).GetSymbol() == "Cl"}
    assert {tuple(sorted(p)) for p in best.trans_pairs()} >= {tuple(sorted(cl))}
    assert not om_jobs.ambiguous(options)


def test_hydrogens_on_donors(tmp_path):
    drawing = read(square_planar(tmp_path, trans=False))
    for d in drawing.donors:
        atom = drawing.mol.GetAtomWithIdx(d)
        assert atom.GetTotalNumHs() == (3 if atom.GetSymbol() == "N" else 0)
    assert drawing.drawn_charge == 0


def test_carbonyl_carbon_gets_no_hydrogen(tmp_path):
    drawing = read(mo_co3_pme3_3(tmp_path, "mer"))
    carbons = [d for d in drawing.donors if drawing.mol.GetAtomWithIdx(d).GetSymbol() == "C"]
    assert len(carbons) == 3
    assert all(drawing.mol.GetAtomWithIdx(c).GetTotalNumHs() == 0 for c in carbons)
    assert drawing.drawn_charge == 0
    assert not [n for n in drawing.notes if n[0] == "warning"]


def test_bare_phosphorus_is_reported_as_incomplete(tmp_path):
    page = Page("Mn")
    for site, (direction, depth) in OCTAHEDRON.items():
        page.carbonyl(direction, depth) if site != "up" else page.simple("P", direction, depth)
    (mol,) = read_raw(page.write(tmp_path / "bare.mol"))
    # A ChemDraw label "P" fixes its hydrogens: none.
    phosphorus = next(a for a in mol.GetAtoms() if a.GetSymbol() == "P")
    phosphorus.SetNoImplicit(True)
    drawing = prepare(mol)
    assert any(code == "organometallic.incomplete" for _, code, _ in drawing.notes)


def test_amido_and_halide_bind_covalently(tmp_path):
    page = Page("Pd")
    n = page.simple("N", (1, 0))
    for direction in ((1, 1), (1, -1)):
        c = page.add("C", *(2.8 * np.asarray(direction) / np.sqrt(2) + [1.0, 0]))
        page.bond(n, c)
    page.simple("Cl", (-1, 0)), page.simple("Br", (0, 1)), page.simple("I", (0, -1))
    (mol,) = read_raw(page.write(tmp_path / "amido.mol"))
    mol.GetAtomWithIdx(n).SetNoImplicit(True)  # "N" as a ChemDraw label: no H
    drawing = prepare(mol)
    assert any(code == "organometallic.anionic_donor" for _, code, _ in drawing.notes)
    assert all(drawing.mol.GetAtomWithIdx(d).GetTotalNumHs() == 0 for d in drawing.donors)
    assert not [note for note in drawing.notes if note[0] == "warning"]


def test_haptic_ligand_is_refused(tmp_path):
    page = Page("Pt")
    a, b = page.add("C", 1.5, 0.7), page.add("C", 1.5, -0.7)
    page.bond(a, b, 2)
    page.to_metal(a), page.to_metal(b)
    page.simple("Cl", (-1, 0))
    with pytest.raises(BackendError, match="adjacent atoms"):
        read(page.write(tmp_path / "ethylene.mol"))


def test_two_metals_are_refused(tmp_path):
    page = Page("Mn")
    other = page.add("Mn", 3.0, 0.0)
    page.bond(0, other)
    page.carbonyl((-1, 0))
    with pytest.raises(BackendError, match="more than one metal"):
        read(page.write(tmp_path / "dimer.mol"))


# -- building in 3D -------------------------------------------------------------------


def _trans_symbols(built, drawing):
    species = built.species
    centre = analyse(species)[0]
    return sorted(tuple(sorted(x.rstrip("0123456789") for x in (a, b))) for a, b, _ in centre.trans_pairs)


@pytest.mark.parametrize("trans", [True, False])
def test_platin_builds_with_its_isomer(tmp_path, trans):
    drawing = read(square_planar(tmp_path, trans))
    options = arrangements(drawing.donors, drawing.drawn, drawing.chelated, labels=drawing.labels)
    built = build(drawing, options[0])
    assert built.species.formula == "Cl2H6N2Pt"
    centre = analyse(built.species)[0]
    assert "square planar" in centre.geometry
    pairs = _trans_symbols(built, drawing)
    assert pairs == ([("Cl", "Cl"), ("N", "N")] if trans else [("Cl", "N"), ("Cl", "N")])


@pytest.mark.parametrize("isomer", ["mer", "fac"])
def test_mer_and_fac_are_read_from_the_wedges(tmp_path, isomer):
    drawing = read(mo_co3_pme3_3(tmp_path, isomer))
    options = arrangements(drawing.donors, drawing.drawn, drawing.chelated, labels=drawing.labels)
    assert len(options) == 2
    built = build(drawing, options[0])
    assert built.species.formula == "C12H27MoO3P3"
    pairs = _trans_symbols(built, drawing)
    expected = [("C", "P"), ("C", "P"), ("C", "P")] if isomer == "fac" else [("C", "C"), ("C", "P"), ("P", "P")]
    assert pairs == expected


@pytest.mark.parametrize("mirror", [False, True])
def test_helicity_survives_drawing_and_building(tmp_path, mirror):
    expected = helicity_of_model(mirror)
    assert expected in ("Delta", "Lambda")
    drawing = read(co_en3_drawing(tmp_path, mirror))
    assert drawing.drawn_charge == 3
    options = arrangements(drawing.donors, drawing.drawn, drawing.chelated, labels=drawing.labels)
    assert len(options) == 2  # Delta and Lambda: the chelates leave nothing else
    built = build(drawing, options[0])
    assert built.species.formula == "C6H24CoN6"
    assert analyse(built.species)[0].helicity == expected


def test_mirror_images_differ():
    assert {helicity_of_model(False), helicity_of_model(True)} == {"Delta", "Lambda"}


# -- the command line ------------------------------------------------------------------


def test_cli_writes_a_gaussian_input_for_a_drawn_complex(tmp_path, capsys):
    path = square_planar(tmp_path, trans=True)
    out = tmp_path / "out"
    code = main(["from-molfile", str(path), "--charge", "0", "--mult", "1", "--yes",
                 "-p", "gaussian_opt_freq_def2svp", "-o", str(out)])
    printed = capsys.readouterr().out
    assert code == 0, printed
    assert "square planar" in printed
    gjf = next(out.glob("*.gjf"))
    text = gjf.read_text()
    assert "def2svp" in text.lower() and "Pt" in text
    assert (out / (gjf.stem + ".m2i.json")).is_file()


def test_cli_needs_the_charge_when_nobody_can_be_asked(tmp_path, capsys):
    code = main(["from-molfile", str(square_planar(tmp_path, trans=True)), "--yes",
                 "-o", str(tmp_path)])
    assert code == 2
    assert "--charge" in capsys.readouterr().err


def test_cli_rejects_an_impossible_spin(tmp_path, capsys):
    code = main(["from-molfile", str(square_planar(tmp_path, trans=True)), "--charge", "0",
                 "--mult", "2", "--yes", "-o", str(tmp_path)])
    assert code == 1
    assert "impossible" in capsys.readouterr().err


def test_mirror_images_are_paired_in_the_list(tmp_path):
    # Without wedges the two enantiomers of [Co(en)3]3+ fit a flat drawing equally.
    drawing = read(co_en3_drawing(tmp_path, mirror=False))
    options = arrangements(drawing.donors, drawing.drawn, drawing.chelated, labels=drawing.labels)
    assert mirror_partners(options, drawing.labels, drawing.chelated) == {0: 1, 1: 0}
    trans = read(square_planar(tmp_path, trans=True))
    flat = arrangements(trans.donors, trans.drawn, trans.chelated, labels=trans.labels)
    assert mirror_partners(flat, trans.labels, trans.chelated) == {}


def test_cli_reports_a_haptic_ligand_as_an_error(tmp_path, capsys):
    page = Page("Pt")
    a, b = page.add("C", 1.5, 0.7), page.add("C", 1.5, -0.7)
    page.bond(a, b, 2)
    page.to_metal(a), page.to_metal(b)
    code = main(["from-molfile", str(page.write(tmp_path / "zeise.mol")), "--yes", "--charge", "0",
                 "--mult", "1", "-o", str(tmp_path)])
    assert code == 1
    assert "crystal structure" in capsys.readouterr().err


# -- the browser page ---------------------------------------------------------------


def _page(path, out):
    from pathlib import Path

    from m2i import config
    from m2i.chem.conformers import ConformerOptions
    from m2i.gui import organometallic_page

    settings = {"conformers": ConformerOptions(seed=0xF00D), "output_dir": Path(out)}
    organometallic_page.render(
        settings, path=Path(path), digest="test", shown_name=Path(path).name,
        format_settings=lambda fmt: config.load_profile("gaussian_opt_freq_def2svp"),
    )


def test_the_page_offers_the_isomers_and_writes_the_input(tmp_path):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest

    path = mo_co3_pme3_3(tmp_path, "mer")
    at = AppTest.from_function(_page, kwargs={"path": str(path), "out": str(tmp_path / "out")},
                               default_timeout=300).run()
    assert not at.exception, [e.value for e in at.exception]
    radio = at.radio(key="om_isomer:test")
    assert len(radio.options) == 2 and "trans" in radio.options[0]
    # Mo(0), d6: the low-spin singlet is offered once the oxidation state is given.
    at.number_input(key=next(n.key for n in at.number_input if n.key.startswith("om_ox:"))).set_value(0).run()
    spin = next(s for s in at.selectbox if s.key.startswith("spin:om:"))
    spin.select_index(next(i for i, o in enumerate(spin.options) if o.startswith("1 "))).run()
    next(b for b in at.button if b.label.startswith("Generate")).click().run()
    assert not at.exception, [e.value for e in at.exception]
    written = list((tmp_path / "out").glob("*.gjf"))
    assert len(written) == 1
    assert "C12H27MoO3P3" in written[0].name
    assert "0 1" in written[0].read_text().splitlines()


# -- what a drawing leaves out ---------------------------------------------------------

SKETCH = Path(__file__).parent / "data" / "drawings" / "bare_phosphine.cdxml"


@pytest.mark.parametrize(
    "spec, symbol, expected",
    [
        ("iPr2", "P", ["ipr", "ipr"]),
        ("PiPr2", "P", ["ipr", "ipr"]),
        ("Ph, Me", None, ["ph", "me"]),
        ("CH2OH", None, ["ch2", "oh"]),
        ("PhCH2OH", None, ["ph", "ch2", "oh"]),
        ("OMe", None, ["ome"]),
        ("C(C)(C)C", None, None),  # SMILES: one group, whatever it is called
    ],
)
def test_labels_are_read(spec, symbol, expected):
    got = groups.names(spec, symbol)
    assert got == expected if expected is not None else len(got) >= 1


def test_an_unknown_group_is_refused_with_advice():
    with pytest.raises(BackendError, match="does not know"):
        groups.names("Qx3")


def test_a_sketch_is_completed_with_hydrogen_and_says_so():
    drawing = read(SKETCH)
    # The label ChemDraw left as text is expanded; the bare phosphorus is not invented.
    assert drawing.completed == {"P4": ["H", "H", "H"]}
    assert drawing.assumed == {"P5": ["H", "H", "H"]}
    codes = {code for _, code, _ in drawing.notes}
    assert "organometallic.uninterpreted_label" in codes
    assert "organometallic.incomplete" in codes
    assert Chem.MolToSmiles(drawing.mol).count("C(C)C") == 2  # the two iPr of the label


def test_the_real_groups_replace_the_hydrogens():
    drawing = read(SKETCH, {"P4": "Cy3", "P5": "Ph3"})
    assert drawing.completed == {"P4": ["Cy"] * 3, "P5": ["Ph"] * 3}
    formula = rdMolDescriptors.CalcMolFormula(Chem.AddHs(drawing.mol))
    assert formula == "C43H65ClP3Pt"
    assert not [note for note in drawing.notes if note[1] == "organometallic.incomplete"]


def test_a_wrong_number_of_groups_is_refused():
    with pytest.raises(BackendError, match="one per missing bond"):
        read(SKETCH, {"P4": "Cy2"})


def test_substituents_for_an_atom_that_needs_none_are_refused():
    with pytest.raises(BackendError, match="nothing is missing"):
        read(SKETCH, {"Pt1": "Me"})


def test_the_cli_takes_the_groups_and_refuses_to_invent_them(tmp_path, capsys):
    # A cation, [PtCl(PR3)3]+: with the groups it asks for, the electron count works out.
    code = main(["from-molfile", str(SKETCH), "--charge", "1", "--mult", "1", "-o", str(tmp_path)])
    assert code == 1
    assert "--sub" in capsys.readouterr().err

    out = tmp_path / "out"
    code = main(["from-molfile", str(SKETCH), "--sub", "P4=Cy3", "--sub", "P5=Ph3",
                 "--charge", "1", "--mult", "1", "--yes", "-o", str(out)])
    printed = capsys.readouterr().out
    assert code == 0, printed
    assert "Cy, Cy, Cy" in printed
    written = next(out.glob("*.gjf"))
    assert "C43H65ClP3Pt" in written.name
    record = json.loads(next(out.glob("*.m2i.json")).read_text())
    assert record["drawing"]["substituents_added"]["P4"] == ["Cy", "Cy", "Cy"]


def test_the_page_asks_for_the_missing_groups(tmp_path):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest

    at = AppTest.from_function(_page, kwargs={"path": str(SKETCH), "out": str(tmp_path / "out")},
                               default_timeout=300).run()
    assert not at.exception, [e.value for e in at.exception]
    assert any("What the drawing leaves out" in m.value for m in at.markdown)
    tables = [e for e in at.main if type(e).__name__ == "Dataframe"]
    assert tables, "the gaps in the drawing should be offered as a table"
    assert set(tables[0].value["Atom"]) == {"P4", "P5"}
    assert list(tables[0].value["Groups"]) == ["", ""]  # blank keeps the hydrogens


def test_a_substituent_given_as_smiles_is_accepted():
    drawing = read(SKETCH, {"P4": "C(C)(C)C, Me, Ph", "P5": "Ph3"})
    assert drawing.completed["P4"] == ["C(C)(C)C", "Me", "Ph"]


@pytest.mark.parametrize("spec", ["Me99999", "Me" * 200])
def test_an_absurd_substituent_is_refused_rather_than_built(spec):
    # The text comes from whoever is at the keyboard: it must not turn into
    # a million atoms on the server.
    with pytest.raises(BackendError):
        read(SKETCH, {"P4": spec})
