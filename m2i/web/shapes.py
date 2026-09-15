"""Engine objects turned into the JSON the front end reads. Presentation only."""

from __future__ import annotations

from ..crystal import jobs as crystal_jobs
from ..crystal.coordination import analyse
from ..types import IssueLog


def issues(log) -> list[dict]:
    return [{"level": i.level, "code": i.code, "message": i.message} for i in log]


def new_issues(log: IssueLog, earlier) -> list[dict]:
    """Warnings a later step added, leaving out what an earlier one showed."""
    seen = {(i.code, i.message) for i in earlier}
    return [x for x in issues(log) if x["level"] != "info" and (x["code"], x["message"]) not in seen]


def molecule(mol) -> dict:
    stereo = mol.stereo
    return {
        "formula": mol.formula,
        "smiles": mol.smiles,
        "inchikey": mol.inchikey,
        "charge": mol.charge,
        "multiplicity": mol.multiplicity,
        "stereo": {
            "centers": len(stereo.centers),
            "bonds": len(stereo.bonds),
            "describe": stereo.describe(),
        },
    }


def coordination(species) -> list[dict]:
    return [
        {
            "label": centre.label,
            "element": centre.element,
            "geometry": centre.geometry,
            "detail": centre.detail,
            "donors": [
                {
                    "label": (f"η{d.hapticity}-" if d.hapticity > 1 else "") + d.label,
                    "distance": round(float(d.distance), 3),
                    "kind": d.kind,
                }
                for d in centre.donors
            ],
            "trans_pairs": [
                {"a": a, "b": b, "angle": round(float(angle), 1)} for a, b, angle in centre.trans_pairs
            ],
            "isomers": list(centre.isomers),
            "helicity": centre.helicity,
            "helicity_detail": centre.helicity_detail,
        }
        for centre in analyse(species)
    ]


def viewer(species) -> dict:
    """What the 3D view needs: an XYZ block and the metals to label."""
    coords = crystal_jobs.centred(species)
    xyz = f"{species.n_atoms}\n{species.formula}\n" + "\n".join(
        f"{e} {x:.4f} {y:.4f} {z:.4f}" for e, (x, y, z) in zip(species.elements, coords)
    )
    labels = [
        {"text": species.labels[i], "x": float(x), "y": float(y), "z": float(z)}
        for i, (x, y, z) in enumerate(coords)
        if species.elements[i] in species.metals
    ]
    return {"xyz": xyz, "labels": labels}
