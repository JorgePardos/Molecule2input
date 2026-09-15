"""Putting back the hydrogens an X-ray experiment did not locate.

The strongest check is the round trip: take structures whose hydrogens *were*
located, remove them, and ask m2i where they go. On the fixtures below every
atom must come back right, or be one m2i marked as doubtful.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("gemmi", reason="the CIF route needs the [crystal] extra")

from m2i.cli import main as cli  # noqa: E402
from m2i.crystal import Species, read_cif  # noqa: E402
from m2i.crystal import hydrogens  # noqa: E402
from m2i.types import IssueLog  # noqa: E402

COD = Path(__file__).parent / "data" / "cod"
TRANSPLATIN, COBALT = COD / "1538403.cif", COD / "2106540.cif"


def read(name):
    return read_cif(COD / name, IssueLog())


def without_hydrogens(species: Species) -> tuple[Species, dict[int, int]]:
    """The species stripped of H, and how many each heavy atom had."""
    keep = [i for i, e in enumerate(species.elements) if e != "H"]
    index = {old: new for new, old in enumerate(keep)}
    stripped = Species(
        elements=[species.elements[i] for i in keep],
        coords=species.coords[keep],
        labels=[species.labels[i] for i in keep],
        bonds=[(index[a], index[b]) for a, b in species.bonds if a in index and b in index],
    )
    truth = {index[i]: sum(species.elements[j] == "H" for j in species.neighbours(i)) for i in keep}
    return stripped, truth


# -- the round trip ----------------------------------------------------------------


@pytest.mark.parametrize("cod", ["7033930", "2010630", "2016683", "2017108", "2006181", "2219943"])
def test_hydrogens_removed_from_a_real_structure_come_back(cod):
    reading = read(f"{cod}.cif")
    for species in reading.species:
        if "H" not in species.elements or species.n_atoms < 4:
            continue
        stripped, truth = without_hydrogens(species)
        sites = {s.atom: s for s in hydrogens._sites(stripped)}
        for atom, expected in truth.items():
            site = sites.get(atom)
            found = site.count if site else 0
            if found != expected:
                assert site is not None and not site.certain, (
                    f"{species.formula} {stripped.labels[atom]}: {found} H instead of "
                    f"{expected}, and not marked as doubtful"
                )


def test_a_face_on_ring_keeps_one_hydrogen_per_carbon():
    stripped, _ = without_hydrogens(read("7033930.cif").species[0])
    sites = hydrogens._sites(stripped)
    assert len(sites) == 10
    assert all(s.count == 1 and s.certain and s.element == "C" for s in sites)


# -- when the formula settles it ------------------------------------------------------


def test_the_ammines_of_transplatin_are_completed_and_match_the_formula():
    reading = read("1538403.cif")
    proposal = hydrogens.plan(reading, 0)
    assert proposal.target == 6 and proposal.total == 6 and proposal.matches
    completed, added = hydrogens.apply(reading.species[0], proposal)
    assert completed.formula == "Cl2H6N2Pt"
    assert len(added) == 2


def test_placed_hydrogens_have_standard_geometry():
    reading = read("1538403.cif")
    completed, _ = hydrogens.apply(reading.species[0], hydrogens.plan(reading, 0))
    platinum = completed.elements.index("Pt")
    for n, element in enumerate(completed.elements):
        if element != "N":
            continue
        for h in (j for j in completed.neighbours(n) if completed.elements[j] == "H"):
            nh = completed.coords[h] - completed.coords[n]
            nm = completed.coords[platinum] - completed.coords[n]
            assert np.linalg.norm(nh) == pytest.approx(1.015, abs=1e-3)
            angle = np.degrees(np.arccos(np.dot(nh, nm) / np.linalg.norm(nh) / np.linalg.norm(nm)))
            assert angle == pytest.approx(109.47, abs=0.5)


def test_eighteen_hydrogens_for_the_hexaammine():
    proposal = hydrogens.plan(read("2106540.cif"), 0)
    assert proposal.target == 18 and proposal.matches


def test_counts_set_by_the_user_replace_the_proposal():
    proposal = hydrogens.plan(read("1538403.cif"), 0).with_counts({"N1": 2})
    assert proposal.total == 5 and not proposal.matches
    changed = next(s for s in proposal.sites if s.label == "N1")
    assert changed.certain and changed.reason == "set by the user"


# -- the questions a crystal cannot answer ------------------------------------------


def test_an_oxygen_bound_only_to_a_metal_is_a_question():
    """Oxo, hydroxo or aqua: the distance proposes one, the chemist decides."""
    from test_crystal import p1_cif  # the same small-crystal helper

    import tempfile

    with tempfile.TemporaryDirectory() as folder:
        path = p1_cif(Path(folder), [("Zn1", "Zn", 0.5, 0.5, 0.5), ("O1", "O", 0.71, 0.5, 0.5)])
        site = next(s for s in hydrogens._sites(read_cif(path, IssueLog()).species[0]))
    assert site.label == "O1" and not site.certain
    assert site.count == 2 and "aqua" in site.reason  # Zn-O 2.1 Å


# -- the command line --------------------------------------------------------------


def test_cli_completes_hydrogens_that_match_the_formula(tmp_path, capsys):
    code = cli(["from-cif", str(TRANSPLATIN), "-o", str(tmp_path), "--charge", "0", "--mult", "1"])
    assert code == 0
    assert "Placed 6 hydrogens" in capsys.readouterr().out
    gjf = next(tmp_path.glob("*.gjf")).read_text(encoding="utf-8")
    assert sum(1 for line in gjf.splitlines() if line.split()[:1] == ["H"]) == 6


def test_cli_stops_and_names_the_doubtful_atoms_when_the_formula_disagrees(tmp_path, capsys):
    code = cli([
        "from-cif", str(TRANSPLATIN), "-o", str(tmp_path),
        "--charge", "0", "--mult", "1", "--h", "N1=2",
    ])
    assert code == 2
    err = capsys.readouterr().err
    assert "do not match the formula" in err and "--accept-hydrogens" in err


def test_cli_accepts_a_proposal_when_told(tmp_path):
    code = cli([
        "from-cif", str(TRANSPLATIN), "-o", str(tmp_path),
        "--charge", "0", "--mult", "2", "--h", "N1=2", "--accept-hydrogens",
    ])
    assert code == 0  # Cl2H5N2Pt: an odd electron count, hence the doublet


def test_cli_can_still_write_without_them(tmp_path):
    code = cli([
        "from-cif", str(TRANSPLATIN), "-o", str(tmp_path),
        "--charge", "0", "--mult", "1", "--allow-missing-hydrogens",
    ])
    assert code == 0
    assert next(tmp_path.glob("*.gjf")).name.startswith("Cl2N2Pt")


def test_a_complete_complex_is_not_held_back_by_water_elsewhere(tmp_path, capsys):
    """In this crystal the 12 unlocated H per formula unit belong to a
    Na/water/Cl network; the [Co(en)3]3+ cations have all theirs."""
    reading = read("2010630.cif")
    assert reading.missing_hydrogens == 12 and reading.extended
    assert not hydrogens.lacking(reading, 0)
    assert hydrogens.plan(reading, 0).target is None  # no wrong attribution

    code = cli(["from-cif", str(COD / "2010630.cif"), "-o", str(tmp_path),
                "--charge", "3", "--mult", "1"])
    assert code == 0
    assert "belong to other species" in capsys.readouterr().out
