# Changelog

## 0.2.0 — 2026-09-14

This version reads metal complexes, runs as a web application, and stops
asking you to confirm structures that could not have been misread.

### Added

- **Metal complexes from crystal structures** (`from-cif`, and the *Crystal
  (.cif)* page). The species is rebuilt from the asymmetric unit with its
  measured geometry; missing hydrogens are put back from the heavy-atom
  geometry and checked against the formula; the coordination sphere is
  reported (polyhedron, trans pairs, cis/trans and fac/mer, Δ/Λ, planar
  chirality, racemic or enantiopure crystal). Oxidation state, charge and spin
  are asked, never guessed.
- **Metal complexes from ChemDraw drawings** (`from-molfile`, and the
  *ChemDraw / molfile* page, which notices the metal). Bonds to the metal are
  read as chemists draw them; every arrangement of the donors is listed, best
  fit to the drawing first; the complex is built in 3D and its trans pairs and
  handedness verified.
- **What a drawing leaves out**: labels ChemDraw kept as plain text (`PiPr2`,
  `PhCH2OH`) are read out of the file and expanded; substituents nobody drew
  are filled with hydrogen, said out loud, and asked for (`--sub P6=iPr2`, a
  question at the terminal, a table in the browser).
- **Web deployment** from a Linux machine of your own through a free
  Cloudflare tunnel: `deploy/lab/m2i.sh` (check, start, url, update…),
  `docs/DEPLOY.md`, a `Dockerfile` that runs on any Docker host, and a hosted
  mode (`M2I_HOSTED=1`) that keeps every visitor in their own session and
  clears sessions unused for a day. `deploy/push_space.sh` publishes to a
  Hugging Face Space instead, which needs a paid plan since July 2026.
- **Photos on the web**: the image includes DECIMER, the model for hand-drawn
  structures, kept loaded between pictures.
- `docs/MANUAL.md`, a full user manual; this changelog.
- Profile `gaussian_opt_freq_def2svp` (B3LYP-D3(BJ)/def2-SVP, ultrafine grid).
- Automatic effective core potentials when a basis set stops short of an
  element.

### Changed

- **You are asked to check a structure only when it is in doubt.** A SMILES,
  a ChemDraw file or a molfile is never asked about. A photo reading is shown
  for confirmation when the model is less than 90% confident, when the
  molecule has any stereochemistry, or when anything was warned about while
  reading it. The rules were set on 120 DECIMER readings of known molecules:
  none of the 7 wrong ones would have gone through unseen, although 4 of them
  came with 93–99% confidence. The decision and its reasons are recorded in
  the provenance file; `batch` writes them into the manifest.
- In the browser, a structure that needs no check is summarised in one line,
  its depiction a click away; one that does holds the Generate button until
  "This is the molecule I drew" is ticked.
- Recognition models are kept loaded between pictures: a reading takes a
  couple of seconds after the first, instead of most of a minute each time.
- `batch` refuses metal-complex drawings with an explanation instead of
  reading them as organic molecules.

### Removed

- **MolScribe**, the model for screenshots of ChemDraw drawings, with the
  consensus between two models and `--backend all`. Uploading the `.cdx` or
  `.cdxml` is always better: the file is read exactly, wedges and labels
  included. A leftover environment can be deleted from
  `~/.m2i/backends/molscribe`.

## 0.1.0

First version: organic molecules from a SMILES, a ChemDraw file or molfile, or
a picture (MolScribe or DECIMER), checked by the user every time, embedded in
3D with stereochemistry verified, and written as Gaussian, ORCA, XYZ or SDF
input, from the command line or a local browser interface.
