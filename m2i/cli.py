"""Command line interface.

The verification step is the point of the whole program, so it is on by default
and has to be waived explicitly with --yes. Batch mode waives it implicitly and
records every warning in the manifest instead.
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
from .recognition.manual import ManualBackend, from_molfile
from .recognition.registry import resolve_backends
from .report import validate
from .types import IssueLog, RecognitionResult
from .writers import WriterError, known_programs

IMAGE_SUFFIXES = {".png", ".jpg", ".jpeg", ".bmp", ".gif", ".tif", ".tiff", ".webp"}
MOLFILE_SUFFIXES = {".mol", ".sdf", ".mdl"}
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
        "-y", "--yes", action="store_true", help="skip the structure verification prompt"
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
    image_parser.add_argument("--molfile", type=Path, help="structure as a .mol/.sdf file")
    image_parser.add_argument(
        "--backend",
        action="append",
        help="recognition backend; repeatable, or 'all'/'auto' (default: auto)",
    )

    mol_parser = subparsers.add_parser(
        "from-molfile", parents=[common], help="build inputs from a .mol/.sdf (e.g. ChemDraw export)"
    )
    mol_parser.add_argument("molfile", type=Path)

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
        help="recognition backend for images; repeatable, or 'all'/'auto'",
    )

    setup_parser = subparsers.add_parser(
        "setup", help="install a recognition backend into its own environment"
    )
    setup_parser.add_argument("backend", nargs="?", help="molscribe or decimer")
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
        backends = [
            b
            for b in resolve_backends(args.backend, log, style=style)
            if b.name != "manual"
        ]
        if not backends:
            print(
                "error: no vision backend is installed yet, so the structure has to "
                "be supplied with --smiles or --molfile.\n"
                f"       This drawing looks {style.replace('_', '-')}; install the "
                f"model suited to it:\n"
                f"         m2i setup {'decimer' if style == 'hand_drawn' else 'molscribe'}\n"
                "       Run 'm2i setup --list' to see both.",
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
    recognition = from_molfile(args.molfile).recognize(None)
    return _process_one(recognition, args)


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
            flag = f"({len(log.warnings)} warning(s))" if log.warnings else ""
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
            "Each one gets its own virtual environment: their dependencies are "
            "mutually\nincompatible and none of them can share yours.\n"
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
    if not args.yes and sys.stdin.isatty() and not backends.status(name)["installed"]:
        print(backends.describe(spec))
        answer = input(f"\nDownload and install {name}? [y/N] ").strip().lower()
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
    try:
        from streamlit.web import cli as stcli
    except ImportError:
        print(
            "error: the GUI needs streamlit.\n       pip install streamlit",
            file=sys.stderr,
        )
        return 1

    app = Path(__file__).parent / "gui" / "app.py"
    sys.argv = ["streamlit", "run", str(app), "--server.port", str(args.port)]
    return stcli.main()


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

    if not args.yes and not _confirm(molecule, options, log):
        print("Aborted; nothing was written.")
        return 130

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

    if not sys.stdin.isatty():
        print(
            "Not an interactive terminal; re-run with --yes once you have checked "
            "the structure.",
            file=sys.stderr,
        )
        return False

    answer = input("Is this the molecule you drew? [y/N] ").strip().lower()
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
        return from_molfile(path).recognize(None)

    return load


def _image_loader(path: Path, args):
    def load(log: IssueLog) -> RecognitionResult:
        from .preprocess import estimate_drawing_style, prepare_image

        image = prepare_image(path, log)
        style, _ = estimate_drawing_style(image)
        backends = [
            b
            for b in resolve_backends(getattr(args, "backend", None), log, style=style)
            if b.name != "manual"
        ]
        if not backends:
            raise BackendError(
                f"{path.name} needs a vision backend; none is installed "
                f"(this drawing looks {style.replace('_', '-')}: "
                f"m2i setup {'decimer' if style == 'hand_drawn' else 'molscribe'})"
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
