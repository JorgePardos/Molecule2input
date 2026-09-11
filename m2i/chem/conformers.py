"""3D embedding: ETKDGv3 + MMFF94s, with a stereochemistry safety net.

The last step is the important one. ETKDG can, on strained or crowded systems,
produce a geometry whose configuration does not match the input -- and a
force-field optimisation can push a marginal centre through planarity. Every
surviving conformer is therefore re-labelled from its own 3D coordinates and
compared against the 2D perception; anything that inverted is discarded.
"""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdMolDescriptors

from ..types import Conformer, IssueLog, StereoSummary
from . import stereo as stereo_mod

DEFAULT_SEED = 0xF00D
DEFAULT_PRUNE_RMS = 0.5
# Symmetry-aware RMSD is combinatorial; fall back to the plain version above this.
BEST_RMS_ATOM_LIMIT = 80


class EmbeddingError(RuntimeError):
    """No usable 3D geometry could be produced."""


@dataclass
class ConformerOptions:
    n_confs: int | None = None  # None -> heuristic from rotatable bonds
    keep: int = 1  # how many lowest-energy conformers to write
    seed: int = DEFAULT_SEED
    prune_rms: float = DEFAULT_PRUNE_RMS
    max_iters: int = 2000
    force_field: str = "mmff94s"  # mmff94s | mmff94 | uff | none


def suggest_n_confs(mol: Chem.Mol) -> int:
    """More rotatable bonds, more starting points."""
    rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rot <= 3:
        return 10
    if rot <= 7:
        return 50
    if rot <= 12:
        return 150
    return 300


def generate(
    mol: Chem.Mol,
    log: IssueLog,
    reference_stereo: StereoSummary,
    options: ConformerOptions | None = None,
) -> tuple[Chem.Mol, list[Conformer]]:
    """Embed, optimise, prune and validate. Returns (mol_with_Hs, conformers)."""
    options = options or ConformerOptions()
    mol_h = Chem.AddHs(mol, addCoords=True)

    n_confs = options.n_confs or suggest_n_confs(mol)
    conf_ids = _embed(mol_h, n_confs, options, log)
    if not conf_ids:
        raise EmbeddingError(
            "ETKDG could not generate any 3D geometry. The structure may be "
            "over-constrained (impossible ring strain or stereochemistry)."
        )

    energies, field_name = _optimise(mol_h, options, log)

    order = sorted(conf_ids, key=lambda cid: energies.get(cid, float("inf")))
    kept = _validate_stereo(mol_h, order, reference_stereo, log)
    if not kept:
        raise EmbeddingError(
            "Every conformer changed the stereochemistry during embedding or "
            "optimisation. The requested configuration may be impossible."
        )

    kept = _prune_by_rms(mol_h, kept, options.prune_rms, log)
    selected = kept[: max(1, options.keep)]
    if len(kept) > len(selected):
        log.info(
            "conformers.selected",
            f"{len(kept)} distinct conformers found; writing the {len(selected)} "
            "lowest in energy.",
        )

    lowest = energies.get(selected[0], None)
    conformers: list[Conformer] = []
    for position, cid in enumerate(selected):
        energy = energies.get(cid)
        conformers.append(
            Conformer(
                index=position,
                elements=[a.GetSymbol() for a in mol_h.GetAtoms()],
                coords=[tuple(xyz) for xyz in mol_h.GetConformer(cid).GetPositions()],
                energy=energy,
                force_field=field_name,
                relative_energy=(
                    None if energy is None or lowest is None else energy - lowest
                ),
            )
        )

    # Return a molecule holding only the selected conformers, in the same order
    # as `conformers`, so callers can index the two together without bookkeeping.
    return _keep_conformers(mol_h, selected), conformers


def _keep_conformers(mol_h: Chem.Mol, conf_ids: list[int]) -> Chem.Mol:
    trimmed = Chem.Mol(mol_h)
    trimmed.RemoveAllConformers()
    for cid in conf_ids:
        trimmed.AddConformer(Chem.Conformer(mol_h.GetConformer(cid)), assignId=True)
    return trimmed


# -- steps ---------------------------------------------------------------


def _embed(
    mol_h: Chem.Mol, n_confs: int, options: ConformerOptions, log: IssueLog
) -> list[int]:
    params = AllChem.ETKDGv3()
    params.randomSeed = options.seed
    params.useSmallRingTorsions = True
    params.enforceChirality = True
    # No pruning here: ETKDG prunes on heavy atoms *before* minimisation, which
    # both discards O-H/N-H rotamers and compares structures that have not yet
    # fallen into their basins. All pruning happens after optimisation instead.
    params.pruneRmsThresh = -1.0
    params.numThreads = 0  # all cores

    conf_ids = list(AllChem.EmbedMultipleConfs(mol_h, numConfs=n_confs, params=params))
    if conf_ids:
        return conf_ids

    log.warn(
        "conformers.retry_random",
        "Distance-geometry embedding failed; retrying with random starting "
        "coordinates (common for macrocycles and cage systems).",
    )
    params.useRandomCoords = True
    params.enforceChirality = True
    return list(AllChem.EmbedMultipleConfs(mol_h, numConfs=n_confs, params=params))


def _optimise(
    mol_h: Chem.Mol, options: ConformerOptions, log: IssueLog
) -> tuple[dict[int, float], str]:
    """Force-field optimise every conformer; returns {conf_id: energy}."""
    requested = options.force_field.lower()
    if requested == "none":
        return {}, "none"

    variant = None
    if requested.startswith("mmff"):
        if AllChem.MMFFHasAllMoleculeParams(mol_h):
            variant = "MMFF94s" if requested == "mmff94s" else "MMFF94"
        else:
            log.warn(
                "conformers.no_mmff_params",
                "MMFF94 has no parameters for this molecule (unusual element or "
                "valence); falling back to UFF, which is cruder.",
            )

    if variant:
        results = AllChem.MMFFOptimizeMoleculeConfs(
            mol_h, maxIters=options.max_iters, mmffVariant=variant, numThreads=0
        )
        field_name = variant
    else:
        results = AllChem.UFFOptimizeMoleculeConfs(
            mol_h, maxIters=options.max_iters, numThreads=0
        )
        field_name = "UFF"

    conf_ids = [c.GetId() for c in mol_h.GetConformers()]
    energies: dict[int, float] = {}
    not_converged = 0
    for cid, (flag, energy) in zip(conf_ids, results):
        energies[cid] = energy
        if flag != 0:
            not_converged += 1
    if not_converged:
        log.info(
            "conformers.not_converged",
            f"{not_converged}/{len(conf_ids)} conformers hit the {options.max_iters}-step "
            "force-field limit. Harmless: the QM optimisation starts from here anyway.",
        )
    return energies, field_name


def _validate_stereo(
    mol_h: Chem.Mol,
    conf_ids: list[int],
    reference: StereoSummary,
    log: IssueLog,
) -> list[int]:
    """Drop conformers whose 3D geometry contradicts the intended configuration."""
    if not reference.fingerprint():
        return conf_ids  # nothing was defined, nothing can be violated

    kept: list[int] = []
    first_report: list[str] | None = None
    for cid in conf_ids:
        single = Chem.Mol(mol_h)
        single.RemoveAllConformers()
        single.AddConformer(Chem.Conformer(mol_h.GetConformer(cid)), assignId=True)
        stereo_mod.assign_from_3d(single)
        mismatches = stereo_mod.compare(reference, stereo_mod.summarize(single))
        if mismatches:
            if first_report is None:
                first_report = mismatches
        else:
            kept.append(cid)

    dropped = len(conf_ids) - len(kept)
    if dropped:
        detail = "; ".join(first_report or [])
        log.warn(
            "conformers.stereo_violation",
            f"{dropped}/{len(conf_ids)} conformers were discarded because the 3D "
            f"geometry did not reproduce the intended configuration ({detail}).",
        )
    return kept


def _prune_by_rms(
    mol_h: Chem.Mol, conf_ids: list[int], threshold: float, log: IssueLog
) -> list[int]:
    """Keep energy-ordered conformers separated by more than `threshold` RMSD."""
    if threshold <= 0 or len(conf_ids) < 2:
        return conf_ids

    probe = _rmsd_probe(mol_h)
    use_best = probe.GetNumAtoms() <= BEST_RMS_ATOM_LIMIT

    kept: list[int] = []
    for cid in conf_ids:
        if all(_rms(probe, cid, other, use_best) > threshold for other in kept):
            kept.append(cid)
    return kept


def _rmsd_probe(mol_h: Chem.Mol) -> Chem.Mol:
    """Throwaway copy for RMSD: heavy atoms plus hydrogens on heteroatoms.

    Pruning on heavy atoms alone -- the usual shortcut -- collapses O-H and N-H
    rotamers, which are genuinely different structures for a QM calculation.
    Carbon-bound hydrogens are dropped because their permutational symmetry
    (methyl groups) makes the symmetry-aware RMSD explode combinatorially.
    """
    editable = Chem.RWMol(mol_h)
    doomed = [
        atom.GetIdx()
        for atom in mol_h.GetAtoms()
        if atom.GetAtomicNum() == 1
        and all(nbr.GetAtomicNum() == 6 for nbr in atom.GetNeighbors())
        and atom.GetDegree() > 0
    ]
    for idx in sorted(doomed, reverse=True):
        editable.RemoveAtom(idx)
    probe = editable.GetMol()
    try:
        Chem.SanitizeMol(probe)
    except Exception:
        return Chem.RemoveHs(Chem.Mol(mol_h))
    return probe


def _rms(probe: Chem.Mol, cid_a: int, cid_b: int, use_best: bool) -> float:
    try:
        if use_best:
            return rdMolAlign.GetBestRMS(probe, probe, prbId=cid_a, refId=cid_b)
        return AllChem.GetConformerRMS(probe, cid_a, cid_b, prealigned=False)
    except Exception:
        return float("inf")  # incomparable -> treat as distinct
