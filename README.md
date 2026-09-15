---
title: m2i - Molecule to input
emoji: ⚗️
colorFrom: indigo
colorTo: blue
sdk: docker
app_port: 7860
pinned: false
license: mit
short_description: From a drawn molecule to a Gaussian or ORCA input
---

# m2i — Molecule2Input

From a molecule — typed, drawn in ChemDraw, photographed from paper, or taken
from a crystal structure — to the input file of a quantum-chemistry
calculation. m2i works out the connectivity, stereochemistry, charge and spin,
builds a 3D geometry that keeps the stereochemistry, and writes a Gaussian
`.gjf`, ORCA `.inp`, `.xyz` or `.sdf`, with a record of everything it did.

```
SMILES, .cdx, .mol ─────────────────────────────┐
photo ──▶ DECIMER ──▶ [ a look, if in doubt ] ──┴─▶ CIP + charge/spin ──▶ 3D ──▶ input file
.cif ──▶ rebuilt from the crystal ──▶ [ species, hydrogens, spin ] ──────────────▶ input file
```

**Documentation:** [user manual](https://github.com/JorgePardos/Molecule2input/blob/main/docs/MANUAL.md) ·
[what changed in 0.2.0](https://github.com/JorgePardos/Molecule2input/blob/main/CHANGELOG.md)

## What it reads

| Input | How it is read | You are asked |
|---|---|---|
| **SMILES** | Exactly | Nothing |
| **ChemDraw** (`.cdx`, `.cdxml`) or **molfile** (`.mol`, `.sdf`) | Exactly, wedges included | Nothing |
| **Metal complex** drawn in ChemDraw | Bonds to the metal read as chemists draw them; every isomer listed, best fit to the drawing first; built in 3D and verified | The substituents the drawing leaves out, the isomer if the drawing does not decide, oxidation state, charge and spin |
| **Photo** of a hand-drawn molecule | DECIMER, a model for hand-drawn structures | To confirm the reading — only when it is in doubt |
| **Crystal structure** (`.cif`) | Rebuilt with its measured geometry; missing hydrogens put back; coordination and chirality at the metal reported | The species, doubtful hydrogens, oxidation state, charge and spin |

If you have the ChemDraw file, use it rather than a picture of it: the file
holds the bonds, wedges and labels exactly.

## When it asks you to check the structure

Only a photo can be misread. A reading is shown for confirmation when the
model is less than 90% confident, when the molecule has any stereochemistry
(DECIMER writes it as text instead of measuring it off the drawing), or when
anything looked wrong while reading it. Everything else goes straight through.

The rules were set by measurement, not guessed: on 120 DECIMER readings of
known molecules, 4 of the 7 wrong ones came back at 93–99% confidence, so the
confidence alone would have let them through; with the three rules together,
none of the 7 would have been accepted unseen. Accepting and confirming are
both recorded, with the reasons. Details in the
[manual](https://github.com/JorgePardos/Molecule2input/blob/main/docs/MANUAL.md#4-photos-of-hand-drawn-molecules).

## Install

Python 3.10 or newer.

```bash
pip install -e ".[gui,crystal]"
```

For photos, add the model once (about 1.5 GB, in an environment of its own):

```bash
m2i setup decimer
```

On Windows, if `m2i` is not found, use `python -m m2i.cli` instead.

## Quick start

```bash
m2i from-smiles "C[C@H](N)C(=O)O" -o out/
```
```bash
m2i from-molfile drawing.cdx -o out/ -p orca_opt_freq
```
```bash
m2i from-image photo.jpg -o out/
```
```bash
m2i from-molfile complex.cdx --sub P6=iPr2 -p gaussian_opt_freq_def2svp -o out/
```
```bash
m2i from-cif structure.cif -o out/
```
```bash
m2i batch structures.smi -o out/
```

Each run prints what it understood and writes the input and a `.m2i.json`
provenance record; for organic molecules, also a `_check.png` with the atoms
numbered as in the input.
`m2i profiles` lists the calculation recipes; `m2i doctor` shows what is
installed.

## Browser interface

```bash
pip install -e ".[web,crystal]"
```
```bash
m2i gui
```

Opens at http://localhost:8501. The same routes, one page: give the structure
(Picture, ChemDraw / molfile, SMILES or Crystal), look at it only when it
deserves a look, choose Gaussian, ORCA, XYZ or SDF and download the file.
Metal complexes and crystals get their own steps, with a 3D view that works
offline. It works on a phone too.

It is a small FastAPI application (`m2i/web/`) serving a plain HTML and
JavaScript front end (`m2i/web/static/`), with no build step and nothing
loaded from other sites. `m2i gui --streamlit` still opens the previous
Streamlit interface.

## Run it as a website

m2i runs as a website from a Linux machine you already have — a lab
workstation, a spare PC — reached from anywhere through a free Cloudflare
tunnel: no open ports, no public IP, HTTPS, and photos included. On that
machine, with Docker installed:

```bash
sh deploy/lab/m2i.sh check
```
```bash
sh deploy/lab/m2i.sh start
```
```bash
sh deploy/lab/m2i.sh url
```

The first start builds the image (15–25 minutes); `url` prints the address to
share. The step-by-step guide — installing Docker, a fixed address on your own
domain, restricting access to your group, troubleshooting — is in
[docs/DEPLOY.md](https://github.com/JorgePardos/Molecule2input/blob/main/docs/DEPLOY.md).

The machine needs Linux (x86_64 or arm64), 4 GB of memory for m2i (the photo
model alone uses 2.4 GB), 10 GB of disk, and outbound internet access. Each
visitor works in a private session; uploads (up to 20 MB) and generated files
stay on that machine, and are removed a day after they were last used.

The same image runs on any Docker host (`docker build -t m2i .`). Hugging
Face Spaces also run it, but Docker Spaces need a paid plan since July 2026;
the block at the top of this file is their configuration and
`deploy/push_space.sh` publishes to one.

## What it protects you from

- **Misread photos** — see above.
- **Undefined stereocentres**, reported by name and never assigned at random.
- **Phantom stereocentres at phosphorus**: a phosphodiester's P is not
  stereogenic, however it is written; a phosphorothioate's is, and is flagged.
- **Silent inversion in 3D**: every conformer's stereochemistry is compared
  with the structure's, and inverted ones are discarded.
- **Invented substituents**: a phosphine drawn as a bare P is filled with
  hydrogen, said out loud, and the real groups asked for; unattended, the
  command line will not write it without `--yes`.
- **Wrong isomers of complexes**: the built complex must have the trans pairs
  and the handedness chosen, or nothing is written.
- **Impossible electronic states**: a multiplicity the electron count forbids
  is refused.
- **Broken input files**: Gaussian's blank-line layout, basis sets that do not
  cover an element (an ECP is added), `gen` without a basis block.
- **Losing track**: every input comes with a provenance record — versions,
  structure, source and confidence, profile, conformers, and every warning and
  decision.

## Status

Verified on Windows 11 with Python 3.11 and DECIMER on CPU, from the command
line and in the browser, locally and in hosted mode. The Docker image has not
been built on a Linux machine yet; every package it pins has Linux wheels for
both x86_64 and arm64.

```bash
python -m pytest -q
```

## Licence

MIT. The bundled 3D viewer, [3Dmol.js](https://3dmol.csb.pitt.edu), is
BSD-3-Clause (`m2i/gui/viewer/`). The crystal structures used in the tests
come from the [Crystallography Open Database](https://www.crystallography.net)
and are public domain (CC0).
