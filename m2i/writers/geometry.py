"""Program-neutral geometry output: .xyz and .sdf."""

from __future__ import annotations

from ..types import IssueLog, JobSpec
from .base import BaseWriter, format_geometry


class XyzWriter(BaseWriter):
    program = "xyz"
    extension = ".xyz"

    def render(self, job: JobSpec, log: IssueLog) -> str:
        comment = job.title or job.name
        energy = job.metadata.get("energy")
        if energy is not None:
            comment = f"{comment} | E({job.metadata.get('force_field', 'FF')})={energy:.4f} kcal/mol"
        body = format_geometry(job.elements, job.coords, width=16, decimals=8)
        return f"{job.n_atoms}\n{comment}\n{body}\n"


class SdfWriter(BaseWriter):
    """SDF keeps the connectivity, so it round-trips back into RDKit cleanly."""

    program = "sdf"
    extension = ".sdf"

    def render(self, job: JobSpec, log: IssueLog) -> str:
        from rdkit import Chem

        mol = job.mol
        if mol is None:
            log.warn(
                "sdf.no_connectivity",
                "No molecule object available; writing an .xyz-style SDF without bonds.",
            )
            return XyzWriter().render(job, log)

        mol = Chem.Mol(mol)
        mol.SetProp("_Name", job.name)
        for key in ("smiles", "inchikey", "charge", "multiplicity", "energy"):
            value = job.metadata.get(key)
            if value is not None:
                mol.SetProp(f"m2i_{key}", str(value))
        return Chem.MolToMolBlock(mol, confId=-1) + "$$$$\n"
