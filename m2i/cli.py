"""Command line interface.

A picture can be misread, so a reading with anything against it -- low
confidence, stereochemistry the model wrote rather than measured, a warning --
is shown and confirmed before anything is written (see m2i.report.review).
A SMILES, a ChemDraw file or a molfile says exactly what it contains and is
not asked about. --yes waives the question; batch mode waives it and writes
the reasons into the manifest instead.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from . import __version__, config, pipeline
from .chem.conformers import ConformerOptions, EmbeddingError
from .chem.sanitize import MoleculeParseError
from .recognition import BackendError, describe_backends, recognize
from .recognition.manual import STRUCTURE_FILE_SUFFIXES, ManualBackend, from_molfile
from .recognition.registry import resolve_backends
from .report import review, validate
from .types import IssueLog, RecognitionResult
from .writers import WriterError, known_programs

IMAGE_SUFFIXES = {".png", ".jpg", ".jpeg", ".bmp", ".gif", ".tif", ".tiff", ".webp"}
MOLFILE_SUFFIXES = set(STRUCTURE_FILE_SUFFIXES)
SMILES_LIST_SUFFIXES = {".smi", ".smiles", ".txt", ".csv"}


# -- argument parsing ----------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="m2i",
        description="Turn a drawn molecule into a quantum-chemistry input file.",
    )
    parser.add_argument("--version", action="version", version=f"m2i {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "-o", "--output-dir", type=Path, default=Path("."), help="where to write (default: .)"
    )
    common.add_argument(
        "-p", "--profile", default="gaussian_opt_freq", help="profile name or .yaml path"
    )
    common.add_argument(
        "--program", choices=known_programs(), help="override the profile's program"
    )
    common.add_argument("--method", help="override the profile's method")
    common.add_argument("--basis", help="override the profile's basis set")
    common.add_argument("--charge", type=int, help="override the derived net charge")
    common.add_argument("--mult", type=int, help="override the derived multiplicity")
    common.add_argument("--name", help="base name for the output files")
    common.add_argument(
        "--conformers", type=int, help="how many conformers to embed (default: automatic)"
    )
    common.add_argument(
        "--keep", type=int, default=1, help="how many lowest-energy conformers to write"
    )
    common.add_argument(
        "--force-field",
        choices=("mmff94s", "mmff94", "uff", "none"),
        default="mmff94s",
        help="force field for the pre-optimisation",
    )
    common.add_argument("--seed", type=int, default=0xF00D, help="random seed (reproducibility)")
    common.add_argument(
        "--rms-threshold",
        type=float,
        default=0.5,
        help="RMSD below which two conformers count as the same",
    )
    common.add_argument(
        "--keep-all-fragments",
        action="store_true",
        help="calculate salts/counter-ions together instead of keeping the largest fragment",
    )
    common.add_argument(
        "-y", "--yes", action="store_true",
        help="do not ask about a doubtful picture reading (the reasons are still recorded)"
    )
    common.add_argument("-v", "--verbose", action="store_true", help="show info messages too")

    smiles_parser = subparsers.add_parser(
        "from-smiles", parents=[common], help="build inputs from a SMILES string"
    )
    smiles_parser.add_argument("smiles")

    image_parser = subparsers.add_parser(
        "from-image", parents=[common], help="build inputs from a picture of a molecule"
    )
    image_parser.add_argument("image", type=Path)
    image_parser.add_argument(
        "--smiles", help="structure read by you (required until a vision backend is installed)"
    )
    image_parser.add_argument(
        "--molfile", type=Path, help="structure as a .mol/.sdf/.cdxml/.cdx file"
    )
    image_parser.add_argument(
        "--backend",
        action="append",
        help="recognition backend (default: the installed one)",
    )

    mol_parser = subparsers.add_parser(
        "from-molfile",
        parents=[common],
        help="build inputs from a ChemDraw (.cdxml/.cdx) or molfile (.mol/.sdf)",
    )
    mol_parser.add_argument("molfile", type=Path)
    mol_parser.add_argument(
        "--isomer",
        type=int,
        help="for a metal complex: which arrangement, as numbered in the list (1 = best fit)",
    )
    mol_parser.add_argument(
        "--sub",
        action="append",
        default=[],
        metavar="ATOM=GROUPS",
        help="substituents a drawing leaves out, e.g. P6=iPr2; repeatable",
    )
    mol_parser.add_argument(
        "--oxidation",
        action="append",
        default=[],
        metavar="EL=N",
        help="for a metal complex: oxidation state of the metal, e.g. Mn=1",
    )

    cif_parser = subparsers.add_parser(
        "from-cif",
        parents=[common],
        help="build inputs from a crystal structure (.cif), keeping its geometry",
    )
    cif_parser.add_argument("cif", type=Path)
    cif_parser.add_argument(
        "--species", type=int, help="which species of the crystal (1 = the listed first)"
    )
    cif_parser.add_argument(
        "--oxidation",
        action="append",
        default=[],
        metavar="EL=N",
        help="oxidation state of a metal, e.g. Fe=2; repeatable",
    )
    cif_parser.add_argument(
        "--keep-xray-hydrogens",
        action="store_true",
        help="do not extend C-H/N-H/O-H to standard lengths",
    )
    cif_parser.add_argument(
        "--allow-missing-hydrogens",
        action="store_true",
        help="write the input without the hydrogens the CIF lacks (not recommended)",
    )
    cif_parser.add_argument(
        "--h",
        action="append",
        default=[],
        metavar="ATOM=N",
        help="hydrogens to put on an atom, overriding the proposal, e.g. O1=2; repeatable",
    )
    cif_parser.add_argument(
        "--accept-hydrogens",
        action="store_true",
        help="accept the proposed hydrogens even when they do not match the formula",
    )
    cif_parser.add_argument(
        "--mirror",
        action="store_true",
        help="write the mirror image (the other enantiomer)",
    )

    batch_parser = subparsers.add_parser(
        "batch", parents=[common], help="process a folder or a list of SMILES"
    )
    batch_parser.add_argument("source", type=Path, help="directory, or .smi/.csv list")
    batch_parser.add_argument(
        "--manifest", type=Path, help="CSV summary (default: <output-dir>/manifest.csv)"
    )
    batch_parser.add_argument(
        "--backend",
        action="append",
        help="recognition backend for images (default: the installed one)",
    )

    setup_parser = subparsers.add_parser(
        "setup", help="install a recognition backend into its own environment"
    )
    setup_parser.add_argument("backend", nargs="?", help="decimer")
    setup_parser.add_argument(
        "--list", action="store_true", help="show what each backend is for"
    )
    setup_parser.add_argument(
        "--remove", action="store_true", help="delete the backend's environment"
    )
    setup_parser.add_argument(
        "--force", action="store_true", help="reinstall from scratch"
    )
    setup_parser.add_argument(
        "--python", type=Path, help="base interpreter to build the venv from"
    )
    setup_parser.add_argument(
        "--yes", action="store_true", help="do not ask before downloading"
    )

    subparsers.add_parser("doctor", help="report the environment and available backends")

    profiles_parser = subparsers.add_parser("profiles", help="list or show job profiles")
    profiles_parser.add_argument("name", nargs="?", help="profile to display")

    gui_parser = subparsers.add_parser("gui", help="launch the browser interface")
    gui_parser.add_argument("--port", type=int, default=8501)
    gui_parser.add_argument(
        "--streamlit", action="store_true",
        help="the previous Streamlit interface instead of the web one",
    )

    return parser


# -- commands ------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "from-smiles":
            return cmd_from_smiles(args)
        if args.command == "from-image":
            return cmd_from_image(args)
        if args.command == "from-molfile":
            return cmd_from_molfile(args)
        if args.command == "from-cif":
            return cmd_from_cif(args)
        if args.command == "batch":
            return cmd_batch(args)
        if args.command == "setup":
            return cmd_setup(args)
        if args.command == "doctor":
            return cmd_doctor(args)
        if args.command == "profiles":
            return cmd_profiles(args)
        if args.command == "gui":
            return cmd_gui(args)
    except config.ProfileError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except (MoleculeParseError, BackendError, EmbeddingError, WriterError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 2


def cmd_from_smiles(args) -> int:
    recognition = RecognitionResult(
        smiles=args.smiles, backend="manual", confidence=1.0
    )
    return _process_one(recognition, args)


def cmd_from_image(args) -> int:
    if not args.image.is_file():
        print(f"error: no such image: {args.image}", file=sys.stderr)
        return 1

    log = IssueLog()
    from .preprocess import estimate_drawing_style, prepare_image

    image = prepare_image(args.image, log)

    if args.molfile:
        backend = from_molfile(args.molfile, smiles=args.smiles)
        recognition = backend.recognize(args.image)
    elif args.smiles:
        recognition = ManualBackend(smiles=args.smiles).recognize(args.image)
    else:
        style, metrics = estimate_drawing_style(image)
        if args.verbose:
            print(f"Drawing style: {style} {metrics}")
        # The manual backend cannot read an image, so it does not count here.
        backends = resolve_backends(args.backend, log)
        if not backends:
            print(
                "error: no vision backend is installed yet, so the structure has to "
                "be supplied with --smiles or --molfile.\n"
                "       A ChemDraw file is read exactly, with no model: use from-molfile.\n"
                "       For photos of hand-drawn structures: m2i setup decimer",
                file=sys.stderr,
            )
            _print_issues(log, args.verbose)
            return 1
        recognition = recognize(
            args.image, backends, log, hand_drawn=(style == "hand_drawn")
        )

    return _process_one(recognition, args, log=log, source_image=image)


def cmd_from_molfile(args) -> int:
    if not args.molfile.is_file():
        print(f"error: no such file: {args.molfile}", file=sys.stderr)
        return 1
    from .organometallic import has_metal, read_raw

    try:
        drawn = read_raw(args.molfile)
    except Exception:  # noqa: BLE001 - the organic reader gives the proper error
        drawn = []
    if any(has_metal(m) for m in drawn):
        return cmd_organometallic(args)
    recognition = from_molfile(args.molfile).recognize(None)
    return _process_one(recognition, args)


def _read_drawing(path: Path, substituents: dict[str, str]):
    from .organometallic import read

    try:
        return read(path, substituents)
    except BackendError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return None


def _show_drawing(drawing, path: Path) -> None:
    print(f"\nMetal complex drawn in {path.name}: {drawing.metal_symbol} with "
          f"{len(drawing.donors)} donor atoms")
    for donor in drawing.donors:
        print(f"  {drawing.names[donor]}")
    for level, _, message in drawing.notes:
        print(f"  {'!' if level == 'warning' else '-'} {message}")


def cmd_organometallic(args) -> int:
    """A drawn metal complex: read it, settle the isomer, build it, ask the rest.

    The isomer is read from the drawing -- bond directions on the page, wedges
    and hashes -- and proposed; when the drawing does not decide, it is asked.
    Oxidation state, charge and spin are asked as for a crystal.
    """
    from .crystal import jobs as crystal_jobs
    from .crystal.coordination import analyse
    from .organometallic import arrangements, build, mirror_partners
    from .organometallic import jobs as om_jobs

    interactive = stdin_is_terminal() and not args.yes
    substituents = {}
    for item in args.sub:
        atom, _, spec = item.partition("=")
        if not atom.strip() or not spec.strip():
            print(f"error: --sub expects ATOM=GROUPS, got {item!r}", file=sys.stderr)
            return 2
        substituents[atom.strip()] = spec.strip()

    drawing = _read_drawing(args.molfile, substituents)
    if drawing is None:
        return 1
    _show_drawing(drawing, args.molfile)

    # What a drawing leaves out is filled with hydrogen, and hydrogen is rarely
    # what was meant: it is the first thing worth asking about.
    guessed = [name for name, added in drawing.completed.items()
               if name not in substituents and set(added) == {"H"}]
    if guessed and interactive:
        print("\nThose atoms are unfinished in the drawing. Write the real groups "
              "(iPr2, Ph2, Cy2, Me, OMe...), or leave blank to keep the hydrogens.")
        for name in guessed:
            answer = _ask(f"  {name}: ")
            if answer:
                substituents[name] = answer
        if any(name in substituents for name in guessed):
            drawing = _read_drawing(args.molfile, substituents)
            if drawing is None:
                return 1
            _show_drawing(drawing, args.molfile)
    elif guessed and not args.yes:
        print("error: hydrogens would be used where the drawing says nothing ("
              + ", ".join(guessed) + "). Give the real groups with --sub P6=iPr2, "
              "or accept the hydrogens with --yes.", file=sys.stderr)
        return 1

    options = arrangements(drawing.donors, drawing.drawn, drawing.chelated, labels=drawing.labels)
    if not options:
        print(f"error: {len(drawing.donors)} atoms are bound to {drawing.metal_symbol}; m2i builds "
              "2- to 6-coordinate complexes from drawings. Is one of those bonds drawn by mistake? "
              "(a bare line ending at the metal is read as a methyl)", file=sys.stderr)
        return 1
    mirrors = mirror_partners(options, drawing.labels, drawing.chelated)
    print("\nArrangements, best fit to the drawing first:")
    for number, option in enumerate(options[:6], start=1):
        twin = f"; mirror image of [{mirrors[number - 1] + 1}]" if number - 1 in mirrors else ""
        print(f"  [{number}] fit {option.misfit:.2f}  {om_jobs.describe(option, drawing)}{twin}")
    if len(options) > 6:
        print(f"  ... {len(options) - 6} more (see --isomer)")

    choice = 1
    if args.isomer:
        if not 1 <= args.isomer <= len(options):
            print(f"error: --isomer must be between 1 and {len(options)}", file=sys.stderr)
            return 2
        choice = args.isomer
    elif om_jobs.ambiguous(options):
        why = ("they are mirror images, and only wedges and hashes on the bonds to the metal "
               "tell those apart" if om_jobs.mirror_only(options, drawing)
               else "they fit it about equally well")
        if interactive:
            print(f"  The drawing does not decide between [1] and [2]: {why}.")
            choice = _ask_int("Which arrangement? ", default=1, lowest=1, highest=len(options))
        else:
            print(f"  ! The drawing does not decide between [1] and [2] ({why}); using [1]. "
                  "Choose another with --isomer.")
    chosen = options[choice - 1]

    try:
        built = build(drawing, chosen, seed=args.seed)
    except BackendError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    species = built.species
    print(f"\nBuilt {species.formula} in 3D (distance geometry + UFF: a starting geometry "
          "for your optimisation):")
    for centre in analyse(species):
        print("  " + "\n  ".join(centre.describe()))
    for level, _, message in built.notes:
        print(f"  ! {message}")

    oxidation = {}
    for item in args.oxidation:
        element, _, value = item.partition("=")
        try:
            oxidation[element.strip()] = int(value)
        except ValueError:
            print(f"error: --oxidation expects EL=N, got {item!r}", file=sys.stderr)
            return 2
    symbol = drawing.metal_symbol
    if interactive and symbol not in oxidation:
        value = _ask_int(f"Oxidation state of {symbol} (blank to skip): ", default=None)
        if value is not None:
            oxidation[symbol] = value
    for line in crystal_jobs.spin_hint(species, oxidation):
        print("  " + line)

    suggestion = drawing.drawn_charge
    print(f"  Charge: the formal charges drawn add up to {suggestion:+d} (a drawing often "
          "leaves charges out: check it).")
    charge = args.charge
    if charge is None and interactive:
        charge = _ask_int("Charge of the complex: ", default=suggestion)
    if charge is None:
        print(f"error: pass --charge for a metal complex (the drawing adds up to {suggestion:+d})",
              file=sys.stderr)
        return 2
    lowest = crystal_jobs.default_multiplicity(species, charge, oxidation)
    mult = args.mult
    if mult is None and interactive:
        mult = _ask_int("Multiplicity: ", default=lowest, lowest=1)
    if mult is None:
        print(f"error: pass --mult for a metal complex (the lowest the electron count "
              f"allows is {lowest})", file=sys.stderr)
        return 2

    profile = config.load_profile(args.profile).merged_with(
        program=args.program, method=args.method, basis=args.basis
    )
    log = IssueLog()
    try:
        written = om_jobs.write(
            built, drawing, source=args.molfile, charge=charge, multiplicity=mult,
            profile=profile, output_dir=Path(args.output_dir), log=log,
            oxidation=oxidation, name=args.name,
        )
    except BackendError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    for issue in log:
        if issue.code.startswith("basis."):
            print(f"  ! {issue.message}")
    print(f"\nWrote {len(written)} file(s):")
    for path in written:
        print(f"  {path}")
    return 0


def cmd_batch(args) -> int:
    entries = _collect_batch_entries(args.source, args)
    if not entries:
        print(f"error: nothing to process in {args.source}", file=sys.stderr)
        return 1

    output_dir = Path(args.output_dir)
    manifest_path = args.manifest or output_dir / "manifest.csv"
    rows = []
    failures = 0

    print(f"Processing {len(entries)} structure(s) into {output_dir}\n")
    for name, load in entries:
        log = IssueLog()
        row = {"name": name, "status": "ok"}
        try:
            # Recognition happens inside the loop so that its warnings land in
            # this entry's log, and therefore in the manifest.
            recognition = load(log)
            row["smiles"] = recognition.smiles
            options = _pipeline_options(args, name=name)
            result = pipeline.run(recognition, options, log)
            molecule = result.molecule
            row.update(
                {
                    "smiles": molecule.smiles,
                    "inchikey": molecule.inchikey,
                    "formula": molecule.formula,
                    "charge": molecule.charge,
                    "multiplicity": molecule.multiplicity,
                    "stereo": molecule.stereo.describe(),
                    "conformers": len(result.conformers),
                    "files": "; ".join(Path(f).name for f in result.written_files),
                    "warnings": len(log.warnings),
                    "warning_detail": " | ".join(i.message for i in log.warnings),
                }
            )
            decision = review.assess(molecule, log)
            row["review"] = f"needed: {decision.summary()}" if decision.needed else "not needed"
            flag = f"({len(log.warnings)} warning(s))" if log.warnings else ""
            if decision.needed:
                flag += " - check it"
            print(f"  ok   {name:<28} {molecule.smiles} {flag}")
        except Exception as exc:
            failures += 1
            row.update({"status": "failed", "warning_detail": str(exc)})
            print(f"  FAIL {name:<28} {exc}")
        rows.append(row)

    _write_manifest(manifest_path, rows)
    print(f"\nManifest: {manifest_path}")
    if failures:
        print(f"{failures} of {len(entries)} structures failed.", file=sys.stderr)
    return 1 if failures == len(entries) else 0


def cmd_from_cif(args) -> int:
    """A crystal structure in, an input with its experimental geometry out.

    Nothing is generated: the molecule is rebuilt from the crystal, so the
    stereochemistry and the shape of every ligand are the measured ones. What
    the crystal cannot say -- oxidation state, charge, spin -- is asked.
    """
    from .crystal import CrystalError, read_cif
    from .crystal import jobs as crystal_jobs

    if not args.cif.is_file():
        print(f"error: no such file: {args.cif}", file=sys.stderr)
        return 1

    log = IssueLog()
    try:
        reading = read_cif(args.cif, log, normalise_hydrogens=not args.keep_xray_hydrogens)
    except CrystalError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    interactive = stdin_is_terminal() and not args.yes

    source = f" (COD {reading.cod_id})" if reading.cod_id else ""
    z = f", Z = {reading.z:g}" if reading.z else ""
    print(f"\nCrystal: {reading.path.name}{source}, {reading.spacegroup}{z}")
    print("Species in the unit cell:")
    for number, species in enumerate(reading.species, start=1):
        metal = f"   metal: {', '.join(sorted(set(species.metals)))}" if species.metals else ""
        print(f"  [{number}] {species.formula:<18} x{species.copies:<3} {species.n_atoms:>4} atoms{metal}")
    for formula in reading.extended:
        print(f"  [-] {formula:<18} extended network (not a molecule)")

    index = 0
    if args.species:
        if not 1 <= args.species <= len(reading.species):
            print(f"error: --species must be between 1 and {len(reading.species)}", file=sys.stderr)
            return 2
        index = args.species - 1
    elif interactive and len(reading.species) > 1:
        index = _ask_int("Which species? ", default=1, lowest=1, highest=len(reading.species)) - 1
    species = reading.species[index]
    asked = crystal_jobs.questions(reading, index)

    if asked.centres:
        print(f"\nCoordination in {species.formula}:")
        for centre in asked.centres:
            print("  " + "\n  ".join(centre.describe()))
    _print_issues(log, args.verbose, header="\nNotes:")

    priorities = _ask_planar_priorities(asked.centres, interactive)
    notes, chiral = crystal_jobs.chirality_notes(reading, asked.centres, priorities)
    if chiral:
        print("\nChirality:")
        for note in notes:
            print(f"  {note}")
        if args.mirror:
            print("  Writing the mirror image, as asked (--mirror).")

    from .crystal.hydrogens import lacking

    completed, hydrogens_added = None, None
    if reading.missing_hydrogens and not lacking(reading, index):
        print(
            f"\nThe hydrogens missing from the crystal belong to other species; "
            f"{species.formula} is complete."
        )
    elif reading.missing_hydrogens and not args.allow_missing_hydrogens:
        outcome = _complete_hydrogens(reading, index, args, interactive)
        if isinstance(outcome, int):
            return outcome
        completed, hydrogens_added = outcome
        species = completed

    oxidation = {}
    for item in args.oxidation:
        element, _, value = item.partition("=")
        try:
            oxidation[element.strip()] = int(value)
        except ValueError:
            print(f"error: --oxidation expects EL=N, got {item!r}", file=sys.stderr)
            return 2

    if asked.needs_answers:
        print(f"\n{species.formula} contains {', '.join(asked.metals)}: the crystal does not "
              "tell the oxidation state, charge or spin, so they are up to you.")
        if interactive:
            for element in asked.metals:
                if element not in oxidation:
                    value = _ask_int(f"Oxidation state of {element} (blank to skip): ", default=None)
                    if value is not None:
                        oxidation[element] = value
        for line in crystal_jobs.spin_hint(species, oxidation):
            print("  " + line)

    suggestion = asked.suggested_charge
    print(f"  Charge: {'suggested ' + format(suggestion, '+d') if suggestion is not None else 'no suggestion'} ({asked.charge_reason})")

    charge = args.charge
    if charge is None and interactive:
        charge = _ask_int("Charge of the species: ", default=suggestion)
    if charge is None and not asked.needs_answers and suggestion is not None:
        charge = suggestion
    if charge is None:
        print(
            "error: the charge is not known; pass --charge"
            + (f" (suggested: {suggestion:+d})" if suggestion is not None else ""),
            file=sys.stderr,
        )
        return 2

    mult = args.mult
    default_mult = crystal_jobs.default_multiplicity(species, charge, oxidation)
    if mult is None and interactive:
        mult = _ask_int("Multiplicity: ", default=default_mult, lowest=1)
    if mult is None and not asked.needs_answers:
        mult = default_mult
    if mult is None:
        print(
            f"error: a metal complex needs its spin state stated; pass --mult "
            f"(the lowest the electron count allows is {default_mult})",
            file=sys.stderr,
        )
        return 2

    profile = config.load_profile(args.profile).merged_with(
        program=args.program, method=args.method, basis=args.basis
    )
    try:
        written = crystal_jobs.write(
            reading,
            index,
            charge=charge,
            multiplicity=mult,
            profile=profile,
            output_dir=Path(args.output_dir),
            log=log,
            oxidation=oxidation,
            name=args.name,
            allow_missing_hydrogens=args.allow_missing_hydrogens,
            mirror=args.mirror,
            priorities=priorities,
            species=completed,
            hydrogens_added=hydrogens_added,
        )
    except CrystalError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    later = [i for i in log if i.code.startswith("basis.")]
    for issue in later:
        print(f"  ! {issue.message}")
    print(f"\nWrote {len(written)} file(s):")
    for path in written:
        print(f"  {path}")
    return 0


def _complete_hydrogens(reading, index, args, interactive):
    """Put back the missing hydrogens: (species, what was added), or an exit code.

    The proposal is applied on its own when it matches the formula. When it
    does not, the doubtful atoms are asked about -- or, without a terminal,
    named, so they can be set with --h.
    """
    from .crystal import hydrogens

    overrides = {}
    for item in args.h:
        atom, _, value = item.partition("=")
        try:
            overrides[atom.strip()] = int(value)
        except ValueError:
            print(f"error: --h expects ATOM=N, got {item!r}", file=sys.stderr)
            return 2

    proposal = hydrogens.plan(reading, index).with_counts(overrides)
    target = f", the formula needs {proposal.target}" if proposal.target is not None else ""
    print(f"\nMissing hydrogens: {proposal.total} proposed{target}.")
    for site in proposal.sites:
        mark = "" if site.certain else "  (assumed)"
        print(f"  {site.label:<8} +{site.count}  {site.reason}{mark}")

    if not proposal.matches and not args.accept_hydrogens:
        if not interactive:
            doubtful = ", ".join(f"{s.label} (+{s.count})" for s in proposal.doubtful) or "none"
            print(
                "error: the proposed hydrogens do not match the formula, so where they "
                f"go is not settled. Doubtful atoms: {doubtful}. Set them with --h ATOM=N, "
                "accept the proposal with --accept-hydrogens, or skip them with "
                "--allow-missing-hydrogens.",
                file=sys.stderr,
            )
            return 2
        answers = {}
        for site in proposal.doubtful:
            options = "/".join(str(o) for o in site.options) if site.options else ""
            hint = f" ({options})" if options else ""
            answers[site.label] = _ask_int(
                f"  {site.label}: {site.reason}. Hydrogens{hint}?", default=site.count, lowest=0
            )
        proposal = proposal.with_counts(answers)
        if not proposal.matches:
            needed = proposal.target if proposal.target is not None else "an unknown number"
            answer = _ask(
                f"  That places {proposal.total} H; the formula needs {needed}. Continue? [y/N] "
            ).lower()
            if answer not in ("y", "yes", "s", "si"):
                print("Aborted; nothing was written.")
                return 130

    completed, added = hydrogens.apply(reading.species[index], proposal)
    print(f"  Placed {sum(s.count for s in proposal.sites)} hydrogens: {completed.formula}.")
    return completed, added


def _ask_planar_priorities(centres, interactive) -> dict[str, str]:
    """Ask which substituent ranks higher where m2i could not decide."""
    priorities = {}
    for centre in centres:
        for ring in centre.planar:
            if ring.certain or not interactive:
                continue
            print(f"\n  Ring {ring.ring_label}: {ring.reason}.")
            answer = _ask(
                f"  Which has the higher CIP priority, {ring.first} or {ring.second}? "
                f"[{ring.first}] "
            )
            priorities[ring.ring_label] = answer if answer in (ring.first, ring.second) else ring.first
    return priorities


def stdin_is_terminal() -> bool:
    """Whether someone can actually answer a question typed at them.

    ``isatty()`` alone is not enough on Windows: input redirected from NUL
    reports itself as a terminal, and the first question then dies with
    EOFError. A real console is one whose mode can be queried.
    """
    try:
        if not sys.stdin or not sys.stdin.isatty():
            return False
    except (AttributeError, ValueError):
        return False
    if sys.platform != "win32":
        return True
    try:
        import ctypes
        import msvcrt

        handle = msvcrt.get_osfhandle(sys.stdin.fileno())
        mode = ctypes.c_uint()
        return bool(ctypes.windll.kernel32.GetConsoleMode(handle, ctypes.byref(mode)))
    except Exception:  # noqa: BLE001 - if in doubt, do not ask
        return False


def _ask(prompt: str) -> str:
    """input() that treats a closed stdin as no answer rather than a crash."""
    try:
        return input(prompt).strip()
    except EOFError:
        print()
        return ""


def _ask_int(prompt: str, *, default=None, lowest=None, highest=None):
    """Ask for an integer; blank takes the default (None when there is none)."""
    shown = f"{prompt.rstrip()} [{default}] " if default is not None else prompt
    while True:
        answer = _ask(shown)
        if not answer:
            return default
        try:
            value = int(answer)
        except ValueError:
            print("  a whole number, please")
            continue
        if lowest is not None and value < lowest:
            print(f"  at least {lowest}, please")
            continue
        if highest is not None and value > highest:
            print(f"  at most {highest}, please")
            continue
        return value


def cmd_setup(args) -> int:
    from . import backends

    if args.list or not args.backend:
        print("Recognition backends m2i can install:\n")
        for spec in backends.SPECS.values():
            state = backends.status(spec.name)
            mark = "installed" if state["installed"] else "not installed"
            print(backends.describe(spec))
            print(f"  state     {mark}  ({state['venv']})\n")
        print(
            "It gets its own virtual environment: its dependencies cannot share yours.\n"
            "Drawings made in ChemDraw need no model: pass the .cdx or .cdxml itself.\n"
            "Install with: m2i setup <name>"
        )
        return 0

    name = args.backend
    if name not in backends.SPECS:
        print(
            f"error: unknown backend {name!r}; known: {', '.join(backends.SPECS)}",
            file=sys.stderr,
        )
        return 2

    if args.remove:
        if backends.remove(name):
            print(f"Removed {backends.venv_dir(name)}")
        else:
            print(f"Nothing to remove for {name}")
        return 0

    spec = backends.SPECS[name]
    if not args.yes and stdin_is_terminal() and not backends.status(name)["installed"]:
        print(backends.describe(spec))
        answer = _ask(f"\nDownload and install {name}? [y/N] ").lower()
        if answer not in ("y", "yes", "s", "si", "sí"):
            print("Aborted.")
            return 130

    try:
        backends.install(
            name, base_python=args.python, force=args.force, on_output=print
        )
    except backends.SetupError as exc:
        print(f"\nerror: {exc}", file=sys.stderr)
        return 1
    return 0


def cmd_doctor(args) -> int:
    import platform

    from rdkit import rdBase

    from .recognition._subprocess import backend_home

    print(f"m2i {__version__}")
    print(f"  python        {platform.python_version()} ({sys.executable})")
    print(f"  rdkit         {rdBase.rdkitVersion}")
    try:
        import streamlit  # noqa: F401

        gui_state = "installed"
    except ImportError:
        gui_state = "not installed  (pip install streamlit)"
    print(f"  gui           {gui_state}")
    print(f"  backend home  {backend_home()}")

    print("\nRecognition backends:")
    for row in describe_backends():
        mark = "available" if row["available"] else "unavailable"
        stereo_note = (
            "returns a molblock: stereochemistry is measured"
            if row["returns_molblock"]
            else "returns a SMILES: stereochemistry is generated"
        )
        print(f"  {row['name']:<12} {mark:<12} best at: {row['strength']}")
        print(f"               {stereo_note}")
        if not row["available"]:
            print(f"               -> {row['reason']}")

    profiles = config.available_profiles()
    print(f"\nProfiles ({len(profiles)}):")
    for name, path in sorted(profiles.items()):
        print(f"  {name:<28} {path}")
    return 0


def cmd_profiles(args) -> int:
    profiles = config.available_profiles()
    if not args.name:
        print(f"{len(profiles)} profile(s):\n")
        for name, path in sorted(profiles.items()):
            try:
                profile = config.load_profile(name)
                print(f"  {name:<28} [{profile.program}] {profile.description}")
            except config.ProfileError as exc:
                print(f"  {name:<28} INVALID: {exc}")
        print(f"\nSearched: {', '.join(str(p) for p in config.profile_search_path())}")
        return 0

    profile = config.load_profile(args.name)
    import yaml

    print(yaml.safe_dump(profile.to_dict(), sort_keys=False, allow_unicode=True))
    return 0


def cmd_gui(args) -> int:
    if args.streamlit:
        try:
            from streamlit.web import cli as stcli
        except ImportError:
            print("error: the Streamlit interface needs streamlit.\n       pip install streamlit",
                  file=sys.stderr)
            return 1
        app = Path(__file__).parent / "gui" / "app.py"
        sys.argv = ["streamlit", "run", str(app), "--server.port", str(args.port)]
        return stcli.main()

    try:
        import fastapi  # noqa: F401
        import uvicorn
    except ImportError:
        print('error: the browser interface needs its web extra.\n       pip install -e ".[web]"',
              file=sys.stderr)
        return 1
    print(f"m2i is running at http://localhost:{args.port}  (Ctrl+C to stop)")
    # One process: the sessions and the loaded photo model live in its memory.
    uvicorn.run("m2i.web.main:app", host="127.0.0.1", port=args.port, workers=1, log_level="warning")
    return 0


# -- shared machinery ----------------------------------------------------


def _process_one(
    recognition: RecognitionResult,
    args,
    *,
    log: IssueLog | None = None,
    source_image=None,
) -> int:
    log = log or IssueLog()
    options = _pipeline_options(args, source_image=source_image)

    molecule = pipeline.prepare_molecule(
        recognition,
        log,
        name=args.name,
        charge=args.charge,
        multiplicity=args.mult,
        keep_all_fragments=args.keep_all_fragments,
    )

    print("\nStructure understood by m2i:")
    print(f"  SMILES        {molecule.smiles}")
    print(f"  Formula       {molecule.formula}")
    print(f"  InChIKey      {molecule.inchikey or '(unavailable)'}")
    print(f"  Charge/mult   {molecule.charge} / {molecule.multiplicity}")
    print(f"  Stereo        {molecule.stereo.describe()}")
    print(f"  Source        {recognition.backend}")
    _print_issues(log, args.verbose)

    if log.has_errors():
        print("\nRefusing to write inputs while errors are unresolved.", file=sys.stderr)
        return 1

    decision = review.assess(molecule, log)
    if not decision.needed:
        print(f"  Review        not needed: {decision.summary()}")
        review.record(log, decision, confirmed=None)
    elif args.yes:
        review.record(log, decision, confirmed=False)
    else:
        print("\nWorth a look before anything is written:")
        for reason in decision.reasons:
            print(f"  - {reason}")
        if not _confirm(molecule, options, log):
            print("Aborted; nothing was written.")
            return 130
        review.record(log, decision, confirmed=True)

    result = pipeline.generate_inputs(molecule, options, log)

    print(f"\nWrote {len(result.written_files)} file(s):")
    for path in result.written_files:
        print(f"  {path}")
    _print_issues(log, args.verbose, header="\nNotes:")
    return 0


def _confirm(molecule, options: pipeline.PipelineOptions, log: IssueLog) -> bool:
    """Show the structure and ask. This is the point of the program."""
    from .report import depict

    preview = Path(options.output_dir) / f"{molecule.name}_check.png"
    try:
        depict.comparison_image(
            molecule.mol,
            preview,
            source_image=options.source_image,
            stereo=molecule.stereo,
            caption=f"{molecule.smiles}   |   {molecule.stereo.describe()}",
        )
        print(f"\nCheck the structure: {preview}")
    except Exception as exc:
        print(f"\n(could not render the check image: {exc})")

    if not stdin_is_terminal():
        print(
            "Not an interactive terminal; re-run with --yes once you have checked "
            "the structure.",
            file=sys.stderr,
        )
        return False

    answer = _ask("Is this the molecule you drew? [y/N] ").lower()
    return answer in ("y", "yes", "s", "si", "sí")


def _pipeline_options(args, *, name: str | None = None, source_image=None):
    profile = config.load_profile(args.profile).merged_with(
        program=args.program, method=args.method, basis=args.basis
    )
    return pipeline.PipelineOptions(
        profile=profile,
        output_dir=Path(args.output_dir),
        name=name if name is not None else args.name,
        charge=args.charge,
        multiplicity=args.mult,
        conformers=ConformerOptions(
            n_confs=args.conformers,
            keep=args.keep,
            seed=args.seed,
            prune_rms=args.rms_threshold,
            force_field=args.force_field,
        ),
        keep_all_fragments=args.keep_all_fragments,
        source_image=source_image,
    )


def _print_issues(log: IssueLog, verbose: bool, header: str = "") -> None:
    text = validate.format_issues(log, verbose=verbose)
    if text:
        if header:
            print(header)
        print(text)


def _collect_batch_entries(source: Path, args) -> list[tuple[str, object]]:
    """Read a directory of structures, or a list of SMILES.

    Each entry is ``(name, loader)`` where the loader takes the entry's own
    IssueLog. Recognising an image is deferred until then so the warnings it
    produces -- low confidence, backends disagreeing -- reach the manifest
    instead of being written to a log nobody reads.
    """
    source = Path(source)
    entries: list[tuple[str, object]] = []

    if source.is_file() and source.suffix.lower() in SMILES_LIST_SUFFIXES:
        for number, line in enumerate(
            source.read_text(encoding="utf-8").splitlines(), start=1
        ):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.replace(",", " ").split()
            smiles = parts[0]
            name = parts[1] if len(parts) > 1 else f"{source.stem}_{number:03d}"
            entries.append((name, _smiles_loader(smiles)))
        return entries

    if source.is_dir():
        for path in sorted(source.iterdir()):
            suffix = path.suffix.lower()
            if suffix in MOLFILE_SUFFIXES:
                entries.append((path.stem, _molfile_loader(path)))
            elif suffix in IMAGE_SUFFIXES:
                entries.append((path.stem, _image_loader(path, args)))
        return entries

    return entries


def _smiles_loader(smiles: str):
    def load(log: IssueLog) -> RecognitionResult:
        return RecognitionResult(smiles=smiles, backend="manual", confidence=1.0)

    return load


def _molfile_loader(path: Path):
    def load(log: IssueLog) -> RecognitionResult:
        from .organometallic import has_metal, read_raw

        try:
            metal = any(has_metal(m) for m in read_raw(path))
        except Exception:  # noqa: BLE001 - the organic reader reports a bad file
            metal = False
        if metal:
            # Read as an organic molecule, a complex comes out as nonsense
            # ([CH3][Pt]...), and a batch cannot ask what a complex needs.
            raise BackendError(
                "a metal complex: build it with from-molfile, which settles the isomer "
                "and asks for the oxidation state, charge and spin"
            )
        return from_molfile(path).recognize(None)

    return load


def _image_loader(path: Path, args):
    def load(log: IssueLog) -> RecognitionResult:
        from .preprocess import estimate_drawing_style, prepare_image

        image = prepare_image(path, log)
        style, _ = estimate_drawing_style(image)
        backends = resolve_backends(getattr(args, "backend", None), log)
        if not backends:
            raise BackendError(
                f"{path.name} needs a vision backend and none is installed "
                "(m2i setup decimer), or give the ChemDraw file instead"
            )
        return recognize(path, backends, log, hand_drawn=(style == "hand_drawn"))

    return load


def _write_manifest(path: Path, rows: list[dict]) -> None:
    columns = [
        "name",
        "status",
        "smiles",
        "inchikey",
        "formula",
        "charge",
        "multiplicity",
        "stereo",
        "conformers",
        "files",
        "warnings",
        "warning_detail",
        "review",
    ]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


if __name__ == "__main__":
    raise SystemExit(main())
