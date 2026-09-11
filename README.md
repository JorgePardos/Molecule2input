# m2i — Molecule2Input

From a drawn organic molecule to a quantum-chemistry input file: read the
structure (stereochemistry included), check it, embed it in 3D, and write the
Gaussian `.gjf` / ORCA `.inp` / `.xyz` / `.sdf`.

```
drawing ──▶ recognition ──▶ [ YOU CHECK IT ] ──▶ CIP + charge/spin ──▶ 3D ──▶ input file
           (uncertain)                              (deterministic RDKit)
```

## Why there is a verification step

The recognition layer is the only statistically uncertain part of the program,
and no single model covers both use cases: on **hand-drawn** structures the
DECIMER family reaches ~73% exact accuracy while image-to-graph models sit
below 11%; on **clean ChemDraw-style** depictions MolScribe and MolNexTR win,
because they predict the molecular graph *with* 2D coordinates and wedge bonds
and therefore derive R/S and E/Z geometrically instead of generating them as
text.

Everything after the recognition — CIP labelling, charge and multiplicity,
ETKDG embedding, MMFF optimisation, file writing — is deterministic. So the
human check sits exactly on the boundary where the errors are.

## Install

RDKit, NumPy, Pillow and PyYAML are the only requirements for the core.

```bash
pip install -e ".[gui,dev]"
```

If `m2i` does not end up on your PATH (common with the Microsoft Store build of
Python), use `python -m m2i.cli` everywhere below.

## Quick start

Install a recognition model (once), then point it at a drawing:

```bash
python -m m2i.cli setup molscribe
```
```bash
python -m m2i.cli from-image drawing.png -o out/
```

Or skip the models entirely and give it the structure yourself:

```bash
python -m m2i.cli from-smiles "C[C@H](N)C(=O)O" -o out/
```

It prints what it understood, writes a check image, and asks for confirmation
before writing anything:

```
Structure understood by m2i:
  SMILES        C[C@H](N)C(=O)O
  Formula       C3H7NO2
  InChIKey      QNAYBMKLOCPYGJ-REOHCLBHSA-N
  Charge/mult   0 / 1
  Stereo        C2: S

Check the structure: out/QNAYBMKLOCPYGJ_check.png
Is this the molecule you drew? [y/N]
```

The check image shows the original next to what m2i parsed, with **atoms
numbered exactly as in the generated input file**, CIP descriptors on every
defined stereocentre, and undefined ones highlighted in red.

Add `--yes` to skip the prompt in scripts.

## Commands

```bash
python -m m2i.cli from-smiles "C[C@H](O)/C=C/C(F)Cl" -o out/ --keep 3
```
```bash
python -m m2i.cli from-molfile drawing.mol -o out/ -p orca_opt_freq
```
```bash
python -m m2i.cli from-image photo.png -o out/
```
```bash
python -m m2i.cli from-image photo.png --backend all -o out/
```
```bash
python -m m2i.cli batch structures.smi -o out/ --manifest out/manifest.csv
```
```bash
python -m m2i.cli gui
```
```bash
python -m m2i.cli doctor
```

`batch` accepts a `.smi`/`.csv` list (`SMILES name` per line) or a folder of
`.mol`/`.sdf` files, and writes a manifest recording every warning.

A `.mol` exported from ChemDraw remains the most reliable input of all: the
wedge bonds are already there, with no model in between.

## Recognition backends

Neither model is installed by default, because their dependencies are mutually
incompatible and neither can share your environment. Each gets its own
virtual environment under `~/.m2i/backends` (override with `M2I_BACKEND_HOME`):

```bash
python -m m2i.cli setup --list
```
```bash
python -m m2i.cli setup decimer
```
```bash
python -m m2i.cli setup molscribe
```

| | MolScribe | DECIMER |
|---|---|---|
| approach | image → graph | image → SMILES sequence |
| best at | clean, ChemDraw-style depictions | hand-drawn structures |
| returns | molblock: 2D layout + wedges | a SMILES string |
| stereochemistry | **measured** off the drawing | generated with the tokens |
| download | ~2.5 GB | ~1.5 GB |

`m2i` picks between them by looking at the image: a software export has a
mathematically pure white background and almost no mid-tones, while a
photographed drawing does not. The choice is always reported, never silent, and
`--backend <name>` overrides it.

With both installed, `--backend all` runs them together and compares InChIKeys:

- **identical** → the strongest confidence signal m2i can give you;
- **same skeleton, different stereochemistry** → a loud warning naming both readings;
- **different molecules** → the higher-confidence reading is kept and flagged for review.

When they disagree only on stereochemistry, the reading that came with a
molblock wins, because a measured wedge beats a generated one.

Remove a backend with `python -m m2i.cli setup <name> --remove`; nothing else
on the system is touched.

## Profiles

A profile is the calculation recipe, kept out of the command line:

```yaml
name: gaussian_opt_freq
program: gaussian
method: b3lyp
basis: 6-31G(d)
dispersion: gd3bj
solvent: {model: smd, name: chloroform}
jobs: [opt, freq]
resources: {mem: 8GB, nproc: 8}
```

Profiles are searched in `$M2I_PROFILE_PATH`, then `./profiles`, then
`~/.m2i/profiles`, then the built-in ones — so a local file shadows a built-in
of the same name. `python -m m2i.cli profiles` lists them.

Anything in a profile can be overridden per run: `--program`, `--method`,
`--basis`, `--charge`, `--mult`, `--conformers`, `--keep`, `--force-field`,
`--seed`.

## What it protects you from

- **Undefined stereocentres.** `FindPotentialStereo` separates "the drawing
  never defined this centre" from "the recognition lost it". Both look
  identical in a SMILES; neither is ever resolved silently.
- **Phantom stereocentres at phosphorus.** A phosphodiester is written with one
  `P=O` and one `P–O⁻`, so RDKit sees four different substituents and calls the
  phosphorus potentially stereogenic. The charge is delocalised over both
  oxygens, the centre is not stereogenic, and the DNA backbone is famously not
  chiral there. Unfiltered, that warning fires on every nucleotide, ATP and
  phospholipid — which teaches you to ignore exactly the warnings this program
  exists to give. The discriminator is the element: in a phosphorothioate the
  pair is `=S` and `O⁻`, resonance does not interchange them, Sp/Rp is real,
  and the warning still fires.
- **Silent inversion during embedding.** Every conformer is re-labelled from
  its own 3D coordinates and compared against the 2D perception. Anything that
  inverted is discarded; if all of them invert, the run fails rather than
  writing the wrong enantiomer.
- **Collapsed rotamers.** Pruning happens *after* the force-field
  minimisation and keeps hydrogens on heteroatoms, so O–H and N–H rotamers
  survive — they are genuinely different structures for a QM calculation.
- **Broken input files.** Gaussian's blank-line layout, `gen`/`genecp` without
  a basis block, basis sets that do not cover a heavy element, impossible
  charge/multiplicity combinations.
- **Losing track.** Every input gets a `<name>.m2i.json` sidecar with the
  SMILES, InChIKey, CIP labels, conformer energies, random seed, versions and
  the full list of warnings.

## Output

For `python -m m2i.cli from-smiles "CCCCCCO" --keep 3 -o out/`:

```
out/
  ZSIAUFGUXNUGDI_c01.gjf     lowest-energy conformer
  ZSIAUFGUXNUGDI_c02.gjf
  ZSIAUFGUXNUGDI_c03.gjf
  ZSIAUFGUXNUGDI_check.png   the verification image
  ZSIAUFGUXNUGDI.m2i.json    provenance
```

## Status

Complete: the pipeline, the CLI, the GUI, and both vision backends.

Verified on Windows 11 with the Microsoft Store build of Python 3.11.9. Both
models install and run on CPU, and on a test depiction of
`C[C@H](O)/C=C/c1ccc(Cl)cc1` both returned exactly that structure —
stereocentre and double-bond geometry included, same InChIKey.

Two things worth knowing if you port this elsewhere:

- Upstream MolScribe pins `opencv-python==4.5.5.64`, which looks like it rules
  out Python 3.11. It does not: that release ships an **abi3** wheel, so the
  `cp36` tag covers every later CPython.
- Its `albumentations` dependency pulls in both a newer NumPy and the
  *headless* opencv distribution. Two opencv packages in one environment write
  the same `cv2` module and whichever installs last wins, so the pins are
  repeated in the second install step rather than stated once.

## Tests

```bash
python -m pytest -q
```
