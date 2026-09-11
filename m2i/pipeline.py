"""Stage orchestration: recognition result in, input files out.

Split into two halves on purpose. ``prepare_molecule`` is cheap and is what the
CLI and the GUI call to build the verification view; ``generate_inputs`` is the
expensive half and only runs once a human has confirmed the structure.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from rdkit import Chem

from . import __version__
from .chem import conformers as conf_mod
from .chem import electronic, sanitize
from .chem import stereo as stereo_mod
from .config import JobProfile
from .report import depict, validate
from .types import (
    Conformer,
    IssueLog,
    JobSpec,
    MoleculeSpec,
    PipelineResult,
    RecognitionResult,
)
from .writers import get_writer

SAFE_NAME = re.compile(r"[^A-Za-z0-9._-]+")


@dataclass
class PipelineOptions:
    profile: JobProfile
    output_dir: Path = Path(".")
    name: str | None = None
    charge: int | None = None
    multiplicity: int | None = None
    conformers: conf_mod.ConformerOptions = field(
        default_factory=conf_mod.ConformerOptions
    )
    keep_all_fragments: bool = False
    write_provenance: bool = True
    write_comparison: bool = True
    source_image: Path | None = None


# -- first half: structure ------------------------------------------------


def prepare_molecule(
    recognition: RecognitionResult,
    log: IssueLog,
    *,
    name: str | None = None,
    charge: int | None = None,
    multiplicity: int | None = None,
    keep_all_fragments: bool = False,
) -> MoleculeSpec:
    """Recognition result -> validated, stereo-annotated molecule."""
    mol = sanitize.mol_from_recognition(recognition, log)
    mol = sanitize.split_fragments(mol, log, keep_largest=not keep_all_fragments)
    sanitize.sanitize(mol, log)

    stereo_mod.assign(mol)
    summary = stereo_mod.summarize(mol)
    stereo_mod.report(summary, log)

    net_charge, mult = electronic.resolve(
        mol, log, charge_override=charge, multiplicity_override=multiplicity
    )

    return MoleculeSpec(
        mol=mol,
        smiles=sanitize.canonical_smiles(mol),
        inchikey=sanitize.inchikey(mol, log),
        formula=sanitize.formula(mol),
        charge=net_charge,
        multiplicity=mult,
        stereo=summary,
        name=safe_name(name or sanitize.suggest_name(mol)),
        source=recognition,
    )


# -- second half: geometry and files -------------------------------------


def generate_inputs(
    molecule: MoleculeSpec, options: PipelineOptions, log: IssueLog
) -> PipelineResult:
    """Molecule -> 3D conformers -> input files on disk."""
    result = PipelineResult(molecule=molecule, issues=log, profile=options.profile)

    mol_h, conformers = conf_mod.generate(
        molecule.mol, log, molecule.stereo, options.conformers
    )
    result.conformers = conformers

    output_dir = Path(options.output_dir)
    writer = get_writer(options.profile.program)
    multiple = len(conformers) > 1

    for conformer in conformers:
        stem = molecule.name + (f"_c{conformer.index + 1:02d}" if multiple else "")
        job = JobSpec(
            name=stem,
            title=_title(molecule, conformer, multiple),
            elements=conformer.elements,
            coords=conformer.coords,
            charge=molecule.charge,
            multiplicity=molecule.multiplicity,
            profile=options.profile,
            metadata={
                "smiles": molecule.smiles,
                "inchikey": molecule.inchikey,
                "charge": molecule.charge,
                "multiplicity": molecule.multiplicity,
                "energy": conformer.energy,
                "force_field": conformer.force_field,
            },
            mol=_single_conformer_mol(mol_h, conformer.index),
        )
        path = writer.write(job, output_dir / stem, log)
        result.written_files.append(str(path))

    if options.write_comparison and molecule.mol is not None:
        try:
            image_path = depict.comparison_image(
                molecule.mol,
                output_dir / f"{molecule.name}_check.png",
                source_image=options.source_image,
                stereo=molecule.stereo,
                caption=f"{molecule.smiles}   |   {molecule.stereo.describe()}",
            )
            result.written_files.append(str(image_path))
        except Exception as exc:
            log.info("report.no_image", f"Comparison image not written: {exc}")

    if options.write_provenance:
        record = validate.provenance(result)
        path = validate.write_provenance(
            output_dir / f"{molecule.name}.m2i.json", record
        )
        result.written_files.append(str(path))

    return result


def run(
    recognition: RecognitionResult,
    options: PipelineOptions,
    log: IssueLog | None = None,
) -> PipelineResult:
    """The whole pipeline, without the human verification gate."""
    log = log or IssueLog()
    molecule = prepare_molecule(
        recognition,
        log,
        name=options.name,
        charge=options.charge,
        multiplicity=options.multiplicity,
        keep_all_fragments=options.keep_all_fragments,
    )
    return generate_inputs(molecule, options, log)


# -- helpers -------------------------------------------------------------


def safe_name(name: str) -> str:
    cleaned = SAFE_NAME.sub("_", (name or "molecule").strip()).strip("_.")
    return cleaned[:60] or "molecule"


def _title(molecule: MoleculeSpec, conformer: Conformer, multiple: bool) -> str:
    """Single-line title; a blank line here would truncate a Gaussian input."""
    parts = [molecule.smiles, molecule.inchikey or molecule.formula]
    if multiple:
        relative = conformer.relative_energy
        marker = f"conf {conformer.index + 1}"
        if relative is not None:
            marker += f" (+{relative:.2f} kcal/mol {conformer.force_field})"
        parts.append(marker)
    parts.append(f"m2i {__version__}")
    return " | ".join(parts).replace("\n", " ")


def _single_conformer_mol(mol_h: Chem.Mol, index: int) -> Chem.Mol:
    """A copy carrying only the requested conformer, for the SDF writer."""
    conf_ids = [c.GetId() for c in mol_h.GetConformers()]
    if not conf_ids:
        return mol_h
    single = Chem.Mol(mol_h)
    single.RemoveAllConformers()
    target = conf_ids[index] if index < len(conf_ids) else conf_ids[0]
    single.AddConformer(Chem.Conformer(mol_h.GetConformer(target)), assignId=True)
    return single
