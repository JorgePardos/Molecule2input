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
pip install -e ".[gui,crystal,dev]"
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
python -m m2i.cli from-molfile drawing.cdxml -o out/ -p orca_opt_freq
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
pictures and drawing files, and writes a manifest recording every warning.

A ChemDraw file (`.cdxml`, `.cdx`) or a `.mol`/`.sdf` export is the most
reliable input of all: the wedge bonds are already there, with no model in
between. A ChemDraw page holding several molecules — a counter-ion drawn
apart, a reagent — is read as fragments: the largest is kept, and the warning
names what was dropped (`--keep-all-fragments` keeps the whole assembly).

## Browser interface

```bash
python -m m2i.cli gui
```

Three steps, in the order you would work:

1. **Structure** — a picture (photo, scan or screenshot), a ChemDraw or
   molfile, or a SMILES. For pictures, the model is picked from how the image
   looks, and with both installed you can have them compare their answers.
2. **Check** — the reading next to the original, with atom numbers and CIP
   labels. The SMILES is editable: correcting a bad reading by hand is the
   fastest way forward, and it works even when the reading could not be parsed.
3. **Output** — choose Gaussian, ORCA, XYZ or SDF. Only the settings that
   program needs are shown (a recipe, method, basis, dispersion, solvent,
   cores, memory), and what you download is the file for that format.

The 3D geometry is computed once per structure. Switching format afterwards
rewrites a text file without re-embedding, so every format you download
carries the same coordinates.

### Serving it to other people

The repository is ready to run as a [Hugging Face Space](https://huggingface.co/docs/hub/spaces-sdks-docker):
the block at the top of this README is the Space configuration, and the
`Dockerfile` builds the image (pinned versions in `requirements-web.txt`).
Create a Space with the Docker SDK (huggingface.co/new-space), then publish
the last commit to it:

```bash
sh deploy/push_space.sh https://huggingface.co/spaces/<user>/<space>
```

The Space gets a snapshot of what the image needs, not this repository's
history: Spaces refuse histories holding binary files outside Git LFS, and the
tests and examples have no business on a public server. git asks for your
Hugging Face username and, as the password, an access token with write
permission.

The same image runs anywhere else:

```bash
docker build -t m2i . && docker run -p 7860:7860 m2i
```

It is a plain Streamlit app with no database and no state of its own, so one
process serves everybody. Everything a visitor does belongs to their own
session: the upload is stored under a name m2i chooses in a private temporary
folder, readings and geometries are cached there (a few at a time), and the
files are written to that session's folder, so two people generating the same
molecule cannot overwrite each other. The image sets `M2I_HOSTED=1`, which
hides what only makes sense on your own machine (a folder to save into, paths
on the server's disk, `m2i setup` instructions). Uploads are capped at 20 MB.
The 3D viewer is bundled, so the page needs no network of its own.

The recognition models are the only heavy part, and only pictures need them:
SMILES, ChemDraw, molfiles and CIFs work without them, and the hosted image
leaves them out.

## Crystal structures (.cif)

For metal complexes and anything else where the 3D shape cannot be trusted to
a force field, start from the crystal. Nothing is generated: the molecule is
rebuilt from the crystal, so stereochemistry, hapticity and the shape of every
ligand are exactly those measured.

```bash
pip install -e ".[crystal]"
```
```bash
python -m m2i.cli from-cif structure.cif -o out/
```

A CIF describes the asymmetric unit of a periodic crystal, not a molecule, so
m2i first rebuilds it: it keeps the major component of any disorder (by
occupancy, not by group number), applies the symmetry, merges atoms on special
positions, finds bonds across cell boundaries and walks each network back into
one piece. A network that never closes — a coordination polymer, a framework —
is reported as such rather than cut arbitrarily.

Then it checks the result against `_chemical_formula_sum`, which is also how
it catches hydrogens the experiment never located. X-ray C–H, N–H and O–H
lengths (about 0.1 Å short) are extended to standard values; hydrides are left
as measured.

**Missing hydrogens are put back.** How many each atom lacks is read from the
heavy-atom geometry — angles and planarity tell sp3 from sp2 from sp, bond
lengths do so for terminal atoms, and an sp2 CH needs both a wide angle *and*
a short bond, because thermal motion shortens bonds and strained chelates open
angles. Atoms bound to a metal are often genuinely ambiguous (oxo, hydroxo or
aqua?), so they are proposed from the M–X distance and marked as assumed. The
total is checked against the formula: when it matches, the hydrogens are placed
at standard lengths, rotors turned clear of other atoms; when it does not, the
assumed atoms are asked about (`--h ATOM=N` without a terminal). Missing
hydrogens that belong to lattice water or a network do not hold back a
complete complex.

Checked by removing the hydrogens from COD structures that had them and
putting them back: on 25 structures used to write the rules every atom came
back right or was flagged as doubtful; on 19 new ones, 628 of 653 atoms came
back right, 23 of the rest were flagged, and the 2 wrong but unflagged atoms
belonged to a badly modelled solvent.

For every metal it reports the coordination sphere — donors, hapticity,
distances, polyhedron (τ4, τ5, trans pairs) and cis/trans or fac/mer isomers —
because CIP labels and InChI, the checks used for organic molecules, say
nothing at a metal centre. It also reports chirality at the metal:

- **Δ/Λ** for octahedral centres with two or more chelates, from the IUPAC
  skew-line convention. The sign was fixed against five COD structures of
  resolved complexes labelled by their authors — implemented from memory it
  came out inverted for all five.
- **Planar chirality** (Rp/Sp, Schlögl's convention) for η5/η6 rings with
  different substituents in a 1,2 or 1,3 pattern, checked against thirteen
  labelled rings and their mirror images. The label needs the substituents
  ranked by CIP, which without bond orders is sometimes undecidable; m2i says
  so and asks which ranks higher.
- Whether the crystal is **racemic** (inversion, mirror or glide in the space
  group: both enantiomers present) or enantiopure, with its Flack parameter.
  `--mirror` writes the other enantiomer.

**What a crystal does not say is asked**: the oxidation state of each metal,
from which the d- or f-electron count and the possible spin states follow; the
charge, suggested from the counter-ions in the cell when they are recognised;
and the multiplicity. A multiplicity whose parity the electron count forbids
is rejected before any file is written. Without a terminal, `--charge` and
`--mult` are required for a metal complex; m2i does not pick a spin state for
you.

Basis sets that stop short of an element (a Pople basis ends at Kr) get an
effective core potential automatically: `genecp` with LANL2DZ in Gaussian,
`def2-TZVP` with `def2-ECP` in ORCA; `ecp_basis` in a profile chooses another.

The browser interface has the same route under **Crystal (.cif)**, with a 3D
view of the species in place of a 2D depiction, an editable table of the
proposed hydrogens, and the chirality notes. The viewer is 3Dmol.js
(BSD-3-Clause), bundled in `m2i/gui/viewer/` so it works offline.

## Metal complexes drawn in ChemDraw

When there is no crystal structure, draw the complex. `from-molfile` (and the
**ChemDraw or molfile** input of the browser interface) notices the metal and
takes a different route, because at a metal a drawing gives the connectivity
but not the geometry:

```bash
python -m m2i.cli from-molfile complex.cdx -p gaussian_opt_freq_def2svp -o out/
```

1. **Reading the bonds to the metal.** A C bound to the metal with only a
   terminal O or N is a carbonyl, cyanide or isocyanide (no hydrogen). A
   labelled P, N, O, S binds through a lone pair (a dative bond, its hydrogens
   those drawn). An unlabelled C is an alkyl or aryl. A lone halogen or H is a
   halide or hydride, and a labelled N, O or S one bond short (amido,
   alkoxide) is bound covalently. What does not add up is reported instead
   of guessed: a phosphine with fewer than three substituents, a label
   ChemDraw did not interpret, or a bare line ending at the metal (read as a
   methyl) — see *What a drawing leaves out*, below.
2. **Which donor is trans to which.** Every distinct arrangement of the donors
   on the polyhedron is listed: octahedral, square planar, tetrahedral,
   trigonal bipyramidal and so on, with chelate rings kept cis and equivalent
   ligands counted once. The list is sorted by how well each fits the bond
   directions on the page, wedges and hashes included. Mirror images are paired
   in the list. Without wedges on the bonds to the metal the drawing cannot
   tell two enantiomers apart; m2i says so, and asks (`--isomer N` without a
   terminal).
3. **Building it in 3D.** RDKit's distance geometry cannot build a metal
   complex on its own, so m2i assembles the distance bounds itself: each
   ligand as RDKit sees it, plus the coordination sphere of the chosen
   arrangement. The handedness is checked, then UFF tidies up the structure,
   with restraints that hold the arrangement in place. The result is
   verified: its trans pairs and its handedness must be the ones chosen.
   Tested against DFT geometries of pincer complexes: right isomer and
   enantiomer, core within about 0.2 Å.
4. **The electronic state is asked**, as for a crystal: oxidation state (for
   the spin options), charge (the formal charges drawn are only a suggestion)
   and multiplicity, checked against the electron count.

Complexes with one metal and 2 to 6 donor atoms are supported. η-bound ligands
(alkenes, Cp, arenes) and polynuclear complexes are not yet supported from
drawings; start from a CIF for those.

### What a drawing leaves out

Schemes are drawn to be read by people: a pincer is drawn with a bare P, and a
ligand is written as a label rather than as atoms. Three different things hide
behind that, and only the last of them is a question for you:

- **A nickname ChemDraw interpreted** (CO, Ph, Ts, Cy). The file stores the
  structure it stands for, so there is nothing to guess — it is simply built.
- **A label ChemDraw did not interpret.** Those stay as plain text on a node
  with no structure, and RDKit turns them into a bare carbon. m2i reads the
  text out of the file and expands it — `PiPr2`, `CH2OH`, `OTf`, `PhCH2OH`,
  or any SMILES — from its first atom, which is where the bond arrives, and
  says so, since a label like `PhCH2OH` does not say which atom binds the
  metal.
- **Substituents nobody drew.** No file holds them: a P with one line to the
  metal is a phosphine whose groups exist only in the chemist's head, and no
  model can recover what was never written — it could only invent it. So the
  valence is filled with hydrogen, which is said out loud every time, and the
  real groups are asked for: `--sub P6=iPr2` (repeatable), a question at the
  terminal, or a table in the browser. Without a terminal, m2i refuses to
  write an input with invented hydrogens unless you pass `--yes`.

Groups are written as chemists write them (`Me`, `iPr2`, `Ph3`, `Cy2`, `tBu`,
`OMe`, `NMe2`, `CF3`, `OTf`, `CH2OH`, `SiMe3`...), one per missing bond, or as
SMILES for anything else. Whatever was added is listed in the run and recorded
in the `.m2i.json`, so the file always says which atoms came from the drawing
and which from you.

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
