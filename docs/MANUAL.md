# m2i user manual

Version 0.2.0. For what changed since the previous version, see
[CHANGELOG.md](../CHANGELOG.md).

m2i turns a molecule into the input file of a quantum-chemistry calculation.
You give it the structure in whatever form you have it; it works out the
connectivity, stereochemistry, charge and spin, builds a 3D geometry that
keeps the stereochemistry, and writes a Gaussian, ORCA, XYZ or SDF file,
together with a record of everything it did.

Contents

1. [Installing](#1-installing)
2. [Which input to use](#2-which-input-to-use)
3. [Organic molecules](#3-organic-molecules)
4. [Photos of hand-drawn molecules](#4-photos-of-hand-drawn-molecules)
5. [Metal complexes drawn in ChemDraw](#5-metal-complexes-drawn-in-chemdraw)
6. [Crystal structures](#6-crystal-structures)
7. [Calculation settings](#7-calculation-settings)
8. [What you get](#8-what-you-get)
9. [The browser interface](#9-the-browser-interface)
10. [Command reference](#10-command-reference)
11. [Troubleshooting](#11-troubleshooting)
12. [Limits](#12-limits)

---

## 1. Installing

Python 3.10 or newer. From the repository folder:

```bash
pip install -e ".[web,crystal]"
```

- The core needs RDKit, NumPy, Pillow and PyYAML.
- `web` adds FastAPI and Uvicorn, for the browser interface (`gui`, the older Streamlit one, still works).
- `crystal` adds gemmi, for CIF files.

On Windows, if `m2i` is not found after installing (common with the Microsoft
Store build of Python), write `python -m m2i.cli` wherever this manual says
`m2i`.

**Reading photos** needs one more step: DECIMER, the model for hand-drawn
structures. It brings TensorFlow, which cannot share an environment with m2i,
so it is installed into its own, under `~/.m2i/backends` (or wherever
`M2I_BACKEND_HOME` points):

```bash
m2i setup decimer
```

About 1.5 GB is downloaded. The installation is only marked ready once the
model has read a probe drawing correctly. Nothing else needs it: SMILES,
ChemDraw files, molfiles and crystal structures work without.

Check what is installed at any time:

```bash
m2i doctor
```

## 2. Which input to use

| What you have | Command | Browser tab | Asked to check it? |
|---|---|---|---|
| A SMILES | `m2i from-smiles` | SMILES | Never |
| A ChemDraw file (`.cdx`, `.cdxml`) or molfile (`.mol`, `.sdf`) | `m2i from-molfile` | ChemDraw / molfile | Never |
| A metal complex drawn in ChemDraw | `m2i from-molfile` | ChemDraw / molfile | You choose the isomer and the electronic state |
| A photo or scan of a hand-drawn molecule | `m2i from-image` | Picture | Only when the reading is in doubt |
| A crystal structure (`.cif`) | `m2i from-cif` | Crystal (.cif) | You choose the species and the electronic state |
| Many structures at once | `m2i batch` | — | Never asked; the manifest says which need a look |

**If you have the ChemDraw file, use it** — never a screenshot of it. The file
holds the bonds, the wedges and the labels exactly; a picture of it can only
be read approximately.

## 3. Organic molecules

### From a SMILES

```bash
m2i from-smiles "C[C@H](N)C(=O)O" -o out/
```

```
Structure understood by m2i:
  SMILES        C[C@H](N)C(=O)O
  Formula       C3H7NO2
  InChIKey      QNAYBMKLOCPYGJ-REOHCLBHSA-N
  Charge/mult   0 / 1
  Stereo        C2: S
  Source        manual
  Review        not needed: the structure was typed in, not recognised from a picture

Wrote 3 file(s):
  out/QNAYBMKLOCPYGJ.gjf
  out/QNAYBMKLOCPYGJ_check.png
  out/QNAYBMKLOCPYGJ.m2i.json
```

The charge comes from the formal charges in the structure, and the
multiplicity from its electron count (the lowest the count allows). Override
either with `--charge` and `--mult`.

### From a ChemDraw file or molfile

```bash
m2i from-molfile drawing.cdx -o out/ -p orca_opt_freq
```

The stereochemistry is read from the wedges as drawn. A page with several
molecules — a counter-ion drawn apart, a reagent — keeps the largest one and
names what was dropped; `--keep-all-fragments` keeps them all. If the drawing
contains a metal, the metal-complex route is taken instead
([section 5](#5-metal-complexes-drawn-in-chemdraw)).

### What is reported rather than fixed

- **A stereocentre the structure never defined** is reported by name (`C4:
  unspecified`) and never assigned at random. Define it in the structure, or
  accept that the input will carry an arbitrary configuration.
- **A phosphate or phosphonate** is not flagged as a stereocentre just
  because it is written with one P=O and one P–O⁻: that phosphorus is not
  stereogenic. A phosphorothioate still is, and is still flagged.
- **A valence m2i had to repair** is named, with what was changed.

Errors stop the run; warnings are printed and recorded.

### Many at once

```bash
m2i batch structures.smi -o out/
```

`batch` takes a `.smi`, `.smiles`, `.txt` or `.csv` list (`SMILES name` on
each line), or a folder of drawing files and pictures. It never stops to ask;
instead it writes `manifest.csv`, one row per structure, with the files
written, every warning, and a `review` column saying whether the structure
needs a look and why. Metal-complex drawings are refused in a batch — they
need the questions of [section 5](#5-metal-complexes-drawn-in-chemdraw) — and
marked as failed with the reason.

## 4. Photos of hand-drawn molecules

```bash
m2i from-image photo.jpg -o out/
```

DECIMER reads the picture and writes the SMILES it sees. The first picture
after starting takes most of a minute, because the model has to be loaded;
after that, a couple of seconds each.

### When you are asked to check the reading

A reading is shown for confirmation — with the reasons — whenever anything
points at a mistake:

1. **The model is less than 90% confident.**
2. **The molecule has any stereochemistry**, defined or not. DECIMER writes
   stereocentres and double-bond geometry as part of the text it generates;
   it does not measure them off the drawing.
3. **m2i warned about something** while reading it: a valence it repaired, a
   fragment it dropped, an undefined stereocentre.

```
Worth a look before anything is written:
  - its stereochemistry (C2: S) was written by the model, not measured off the drawing

Check the structure: out/QNAYBMKLOCPYGJ_check.png
Is this the molecule you drew? [y/N]
```

Open the check image: the photo next to the structure as understood, with the
atoms numbered as in the input file, CIP labels on every stereocentre and
undefined ones in red. Answer `y` to write the input; anything else writes
nothing.

A reading with none of the three is accepted without a question, and m2i says
so:

```
  Review        not needed: decimer is 97% confident, the molecule has no stereochemistry, and nothing looked wrong
```

**Why these rules and not the confidence alone.** 60 known molecules were
rendered by RDKit, clean and degraded like a photo (rotation, blur, noise,
JPEG compression), and read by DECIMER through m2i. 113 of the 120 readings
were exact. Of the 7 wrong ones, 4 came back at 93–99% confidence — a
penicillin with its stereocentres wrong, a resveratrol that lost its E, an
amphetamine with a stereocentre invented — so no confidence threshold catches
them. With the three rules together, none of the 7 would have been accepted
unseen, and 70 of the 113 right readings went through without a question.
Rendered depictions are not hand drawings, and real photos are misread more
often; expect more questions with them.

### When the reading is wrong

- At the terminal, answer `n`, then give the structure yourself:
  `m2i from-image photo.jpg --smiles "the right SMILES"`. The photo still
  appears in the check image.
- In the browser, edit the SMILES under the reading; the correction replaces
  it, and "Back to the original reading" undoes it.

### Without a terminal

In scripts, a reading that needs a look is not written: the run ends with
exit code 130. `--yes` accepts it anyway; the reasons are still recorded in
the provenance file.

### Taking better photos

Dark lines on light paper, the whole molecule in the frame, taken from above
rather than at an angle. Write stereocentres with clear wedges and hashes, and
atom labels in capitals.

## 5. Metal complexes drawn in ChemDraw

```bash
m2i from-molfile complex.cdx -p gaussian_opt_freq_def2svp -o out/
```

At a metal, a drawing fixes which atoms are bound but not the geometry, so
this route goes through four steps.

### 5.1 How the bonds to the metal are read

| Drawn as | Read as |
|---|---|
| C bound to the metal whose only other neighbour is a terminal O or N | Carbonyl, cyanide or isocyanide (no hydrogen) |
| A labelled P, N, O, S… | A donor through its lone pair (dative bond); its hydrogens are those written in the label |
| An unlabelled C | An alkyl or aryl ligand, its bond to the metal counted |
| A halogen or H on its own | Halide or hydride |
| A labelled N, O or S one bond short (`N` in a pincer backbone, the `O` of an acac) | Amido, alkoxide or thiolate, bound covalently |

### 5.2 What the drawing leaves out

Schemes are drawn for people: a pincer with a bare P, a ligand written as a
label. m2i separates three cases:

- **A nickname ChemDraw interpreted** (`CO`, `Ph`, `Cy`, `Ts`…). The file
  stores the structure behind it, and it is built as such.
- **A label ChemDraw did not interpret** — it stays as plain text, and
  ChemDraw shows it in red. m2i reads the text out of the file and expands it:
  `PiPr2`, `PPh2`, `NMe2`, `OTf`, `CH2OH`, `PhCH2OH`… A label is read from its
  first atom, which is where the bond arrives; if the metal binds through
  another atom of the label, draw that atom and its bond.
- **Substituents nobody drew** — a P with a single line to the metal. No file
  holds them, so no program can recover them. The valence is filled with
  hydrogen, you are told, and asked for the real groups.

Give the groups on the command line, one `--sub` per atom, with the atom names
m2i prints:

```bash
m2i from-molfile complex.cdx --sub P6=iPr2 --sub P7=iPr2
```

Or answer when asked at the terminal, or fill in the table in the browser.
Write one group per missing bond, as chemists write them: `Me`, `Et`, `iPr`,
`tBu`, `Cy`, `Ph`, `Bn`, `Mes`, `Tol`, `CF3`, `OMe`, `OEt`, `OPh`, `OAc`,
`OTf`, `NMe2`, `SMe`, `SiMe3`, `CH2OH`, and counts (`iPr2`, `Ph3`), lists
(`Me, Ph`), or SMILES for anything else (`C(C)(C)C`). The element may be
included or not: `P6=PiPr2` and `P6=iPr2` are the same.

Without a terminal, m2i will not write an input with invented hydrogens: give
the groups with `--sub`, or accept the hydrogens with `--yes`.

### 5.3 Which donor is trans to which

Every distinct arrangement of the donors on the polyhedron is listed —
linear, trigonal planar, tetrahedral, square planar, trigonal bipyramidal,
square pyramidal, octahedral — with chelates kept cis and equivalent ligands
counted once. They are sorted by how well each fits the directions of the
bonds on the page, wedges and hashes included:

```
Arrangements, best fit to the drawing first:
  [1] fit 0.00  square planar; trans: Cl3 / P2, P4 / P5
  [2] fit 0.61  tetrahedral (chiral)
  [3] fit 0.61  tetrahedral (chiral)
  [4] fit 1.00  square planar; trans: Cl3 / P5, P2 / P4
```

The best fit is used. When the first two fit about equally well, you are
asked; without a terminal, choose with `--isomer N`. Mirror images are paired
in the list (`mirror image of [2]`): only wedges and hashes on the bonds to
the metal can tell them apart, so **draw at least one wedge on a bond to the
metal** if the enantiomer matters.

### 5.4 Building it, and the electronic state

The complex is built by distance geometry with bounds assembled for its
coordination sphere, pre-optimised with UFF, and checked: its trans pairs and
its handedness must be the ones chosen, or nothing is written. Against DFT
geometries of pincer complexes, the right isomer and enantiomer came out, with
the core within about 0.2 Å. It is a starting geometry for your
optimisation, not a final structure.

Then you are asked for the **oxidation state** of the metal (only used to
list the spin states that make sense: `Co(+3) is d6: multiplicity 1 (low
spin), 3 (intermediate spin), 5 (high spin)`), the **charge** (the formal
charges in the drawing are offered as a suggestion — drawings often leave
charges out) and the **multiplicity**. A multiplicity the electron count
forbids is refused. Without a terminal, `--charge` and `--mult` are required;
`--oxidation Mn=1` gives the oxidation state.

Supported: one metal with 2 to 6 donor atoms. Not yet: η-bound ligands
(alkenes, Cp, arenes) and more than one metal — start from a crystal
structure for those.

## 6. Crystal structures

```bash
m2i from-cif structure.cif -o out/
```

Nothing is generated: the molecule is taken from the crystal with its measured
geometry, so its stereochemistry, hapticity and ligand shapes are exactly
those of the experiment. CIFs from the CSD, the COD or your own refinement
work.

### What m2i does with the file

- Rebuilds the molecules from the asymmetric unit: keeps the major component
  of any disorder, applies the symmetry, merges atoms on special positions,
  joins molecules split across cell boundaries. A network that never closes (a
  coordination polymer, a framework) is reported, not cut.
- Lists the **species** in the cell (`[1] C6H24CoN6 x2 … metal: Co`). When
  there are several you are asked which; without a terminal, choose with
  `--species N`, or the first is used.
- Checks everything against `_chemical_formula_sum`, and extends X-ray C–H,
  N–H and O–H bonds (about 0.1 Å short) to standard lengths
  (`--keep-xray-hydrogens` leaves them).

### Missing hydrogens

X-ray experiments often do not locate every hydrogen. m2i proposes how many
each atom lacks from the heavy-atom geometry, marks the proposals it is not
sure of (atoms bound to a metal: oxo, hydroxo or aqua?), and checks the total
against the formula. When it matches, the hydrogens are placed; when it does
not, you are asked about the doubtful atoms. Without a terminal:

- `--h O1=2` sets the hydrogens on an atom (repeatable);
- `--accept-hydrogens` accepts the proposal even if it disagrees with the formula;
- `--allow-missing-hydrogens` writes the input without them (not recommended:
  the electron count changes).

### Coordination and chirality

For every metal: donors, distances, polyhedron (with τ4 or τ5), trans pairs,
cis/trans and fac/mer isomers, **Δ/Λ** for octahedral centres with two or more
chelates, **planar chirality** (Rp/Sp) of substituted η5/η6 rings, and
whether the crystal is racemic or enantiopure (with its Flack parameter).
`--mirror` writes the other enantiomer. When the CIP ranking needed for a
planar-chirality label cannot be settled from the crystal, you are asked which
substituent ranks higher.

### The electronic state

As for drawn complexes: oxidation state, charge (suggested from the
counter-ions in the cell when they are recognised) and multiplicity, all
yours to decide; `--oxidation`, `--charge` and `--mult` without a terminal.

## 7. Calculation settings

### Profiles

A profile is a calculation recipe. The built-in ones:

| Profile | Program | Recipe |
|---|---|---|
| `gaussian_opt_freq` (default) | Gaussian | Opt + freq, B3LYP-D3(BJ)/6-31G(d) |
| `gaussian_opt_freq_def2svp` | Gaussian | Opt (maxcycle=200) + freq (noraman), B3LYP-D3(BJ)/def2-SVP, ultrafine grid |
| `gaussian_opt_freq_solvent` | Gaussian | Opt + freq in SMD chloroform, B3LYP-D3(BJ)/6-31G(d,p) |
| `gaussian_sp_dft` | Gaussian | Single point, M06-2X/def2-TZVP |
| `orca_opt_freq` | ORCA | Opt + freq, B3LYP-D3(BJ)/def2-SVP, RIJCOSX |
| `xyz_only` | — | Geometry only |

```bash
m2i profiles
m2i profiles orca_opt_freq
```

Choose one with `-p NAME`, or `-p path/to/recipe.yaml`. Write your own:

```yaml
name: my_recipe
description: What it is for
program: gaussian          # gaussian, orca, xyz, sdf
method: b3lyp
basis: def2svp
dispersion: gd3bj
solvent: {model: smd, name: water}
jobs: [opt, freq]
extra_keywords: [int=grid=ultrafine]
resources: {mem: 16GB, nproc: 8}
```

Profiles are looked for in `$M2I_PROFILE_PATH`, `./profiles`,
`~/.m2i/profiles`, then the built-in ones; a file of the same name earlier in
that list wins.

Any of it can be overridden per run: `--program`, `--method`, `--basis`,
`--charge`, `--mult`.

### Heavy elements

When a basis set stops short of an element (6-31G(d) ends at Kr), an effective
core potential is added automatically — `genecp` with LANL2DZ in Gaussian,
def2-ECP in ORCA — and you are told. `ecp_basis` in a profile chooses another.

### Conformers (organic molecules)

m2i embeds several 3D conformers, optimises them with a force field and keeps
the lowest in energy:

- `--conformers N` — how many to embed (by default 10 to 300, more with more rotatable bonds);
- `--keep N` — how many to write (default 1);
- `--force-field` — `mmff94s` (default), `mmff94`, `uff` or `none`;
- `--rms-threshold` — below this RMSD (Å, default 0.5) two conformers are the same;
- `--seed` — for reproducible geometries.

Every conformer's stereochemistry is compared with the structure's; one that
inverted is discarded, and if all did, nothing is written.

## 8. What you get

For an organic molecule, in the output folder (`-o`, default the current one),
named by InChIKey unless `--name` is given:

| File | What it is |
|---|---|
| `NAME.gjf` / `.inp` / `.xyz` / `.sdf` | The input; `NAME_c01`, `NAME_c02`… when several conformers are kept |
| `NAME_check.png` | The structure as understood, atoms numbered as in the input, next to the photo when there was one |
| `NAME.m2i.json` | Provenance |

For complexes and crystals, the input and the provenance file, named by
formula and source file (`C6H24CoN6_2010630.gjf`).

The **provenance file** records the m2i, RDKit and Python versions, the
molecule (SMILES, InChIKey, formula, charge, multiplicity, stereochemistry,
where it came from and with what confidence), the profile, the conformers and
their energies, the files written, and every issue — including whether the
structure was accepted without review, confirmed or waived, and why. For
complexes, also the arrangement chosen, the substituents you added and the
hydrogens assumed; for crystals, the species, hydrogens added and chirality.

## 9. The browser interface

```bash
m2i gui
```

Opens at http://localhost:8501 (`--port` for another). The same routes as the
command line, on one page, with the steps listed on the left:

1. **Structure** — choose Picture, ChemDraw / molfile, SMILES or Crystal
   (.cif), and give it: drop the file or type the SMILES. A photo is read as
   soon as it is given.
2. **Check** — for a structure that needs no check, one line ("no check
   needed") with its depiction a click away. For a photo reading in doubt, the
   photo next to the reading, the reasons, and two buttons: "Not this
   molecule" (then correct the SMILES) and "This is the molecule I drew",
   which opens the output. The SMILES can be edited in both cases.
3. **Output** — choose Gaussian, ORCA, XYZ or SDF. Only that program's
   settings are shown (recipe, method, basis, dispersion, solvent, cores,
   memory), and **File name** says what the files will be called — leave it
   blank and m2i names them after the structure, or type your own
   (`pincer_Mn_opt`); the line underneath shows the name as it will be
   written, spaces and dots turned into underscores. Press Generate, preview
   the file and download it with its provenance record.

**Advanced**, under the steps, overrides the charge and multiplicity and sets
the conformer search. On a phone the steps become a progress bar at the top.

A ChemDraw file with a metal opens the complex page: the table of what the
drawing leaves out, the list of arrangements, a 3D view of the built complex
with its coordination sphere, the electronic state, and the output. A CIF
opens the crystal page: the species, the missing-hydrogens table, the 3D view,
coordination and chirality, the electronic state, and the output.

The 3D viewer (3Dmol.js) is bundled with m2i and works without a network.

### Hosted

m2i also runs as a website for other people, from a Linux machine of your own
behind a Cloudflare tunnel; see [DEPLOY.md](DEPLOY.md). What differs from
running it on your own computer:

- each visitor works in a private session; uploads and outputs stay in it,
  and are removed a day after they were last used;
- there is no "save to folder" option and no server paths are shown: what you
  download is the file;
- uploads are limited to 20 MB;
- the photo model is loaded once for the whole server, in the background, as
  soon as someone opens the page; right after the server starts, the first
  photo can take a minute.

## 10. Command reference

Options common to every `from-*` command and `batch`:

| Option | Meaning |
|---|---|
| `-o, --output-dir DIR` | Where to write (default: current folder) |
| `-p, --profile NAME` | Profile name or `.yaml` path (default: `gaussian_opt_freq`) |
| `--program`, `--method`, `--basis` | Override the profile |
| `--charge N`, `--mult N` | Override the charge and multiplicity |
| `--name NAME` | Base name of the output files |
| `--conformers N`, `--keep N`, `--force-field FF`, `--rms-threshold X`, `--seed N` | Conformer search ([section 7](#conformers-organic-molecules)) |
| `--keep-all-fragments` | Keep counter-ions and other fragments |
| `-y, --yes` | Do not ask about a doubtful photo reading or invented hydrogens; the reasons are still recorded |
| `-v, --verbose` | Show informational notes too |

Commands:

| Command | Arguments and options |
|---|---|
| `m2i from-smiles SMILES` | — |
| `m2i from-molfile FILE` | `--sub ATOM=GROUPS` (repeatable), `--isomer N`, `--oxidation EL=N` |
| `m2i from-image PICTURE` | `--smiles SMILES` or `--molfile FILE` to give the structure yourself; `--backend NAME` |
| `m2i from-cif FILE` | `--species N`, `--oxidation EL=N` (repeatable), `--h ATOM=N` (repeatable), `--accept-hydrogens`, `--allow-missing-hydrogens`, `--keep-xray-hydrogens`, `--mirror` |
| `m2i batch SOURCE` | `--manifest FILE`, `--backend NAME` |
| `m2i setup [decimer]` | `--list`, `--remove`, `--force`, `--python PATH`, `--yes` |
| `m2i profiles [NAME]` | — |
| `m2i doctor` | — |
| `m2i gui` | `--port N` |

Exit codes: 0 written; 1 an error (the message says which); 2 an invalid
option value; 130 not confirmed, nothing written.

## 11. Troubleshooting

**"no vision backend is installed yet"** — run `m2i setup decimer`, or give
the structure with `--smiles`, or use the ChemDraw file with `from-molfile`.

**"Not an interactive terminal; re-run with --yes once you have checked the
structure."** — a photo reading needed a look and nobody could be asked.
Check the `_check.png` written, then run again with `--yes`.

**The photo was read as the wrong molecule** — correct it with `--smiles`,
or edit the SMILES in the browser. Photos with faint lines, shadows or a
strong angle are read worst.

**A phosphine comes out as PH2 or PH3** — the drawing has a bare P. Give its
groups: `--sub P6=iPr2`.

**A label is red in ChemDraw and its atoms are wrong in m2i** — ChemDraw did
not interpret it. m2i expands it from its first atom; if the metal binds
through another atom (the O of `PhCH2OH`), draw that atom and its bond.

**"7 atoms are bound to Mn"** — a line in the drawing ends at the metal
without a label and was read as a methyl. Delete it or label it.

**Two arrangements are "mirror images" and it asks which** — the drawing has
no wedge or hash on the bonds to the metal. Draw one, or choose with
`--isomer`.

**"impossible electronic state"** — the multiplicity does not fit the number
of electrons for that charge (an even count needs an odd multiplicity, and the
other way round). Check the charge first: drawings often leave it out.

**"UFF lacks parameters for part of this complex"** — the geometry is the
distance-geometry one, not pre-optimised. It is still a valid starting point;
expect your optimisation to take a few more steps.

**Crystal: "an extended network"** — that part of the crystal is a polymer or
framework with no discrete molecule; choose another species.

## 12. Limits

- Photo reading is statistical. DECIMER reaches about 73% exact readings on
  hand-drawn structures in its authors' benchmarks; m2i's review rules catch
  the errors they were measured on, not every possible one.
- Drawn complexes: one metal, 2–6 donors, no η-bound ligands.
- Geometries of drawn complexes and organic conformers are starting points
  from force fields, not optimised structures.
- m2i does not choose spin states or oxidation states for metals; it lists
  what is possible and refuses what is impossible.
- Batch mode handles organic molecules; complexes and crystals are built one
  at a time.
