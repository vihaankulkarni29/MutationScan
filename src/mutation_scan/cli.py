"""Command-line interface.

    mutationscan run --genomes DIR --references DIR --out DIR

One command in, one folder out. The gene list is not an argument: it comes from
the reference files present in ``--references``, so pointing the tool at a
different set of proteins is the whole of "configuring" it.

The heavy lifting lives in the stage modules; this file only wires them together
in order and hands them paths from a validated :class:`~mutation_scan.config.Config`.
For large cohorts prefer the Snakemake workflow, which runs the same functions
with checkpointing and cluster support.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

from . import __version__
from .config import Config, ConfigError, load_config
from .cooccurrence import cooccurrence
from .diagnostics import diagnose_empty_extraction, diagnose_empty_variants
from .epistasis import epistasis_networks
from .extract import STATUS_EXTRACTED, extract
from .preflight import preflight
from .references import prepare_references
from .summary import write_run_summary
from .variants import VARIANT_STATUS_CALLED, call_variants

logger = logging.getLogger(__name__)

__all__ = [
    "EXIT_EMPTY_RESULT",
    "EXIT_INPUTS_MISSING",
    "EXIT_INTERRUPTED",
    "EXIT_OK",
    "EXIT_USAGE",
    "build_parser",
    "main",
    "run",
]

#: Exit codes. A caller scripting this tool can branch on these: 2 means the run
#: never started, 3 means it started and found nothing.
EXIT_OK = 0
EXIT_INPUTS_MISSING = 1  # config-check only: config valid, data absent
EXIT_USAGE = 2  # bad usage, bad config, bad environment, failed preflight
EXIT_EMPTY_RESULT = 3  # the run completed but produced no usable result
EXIT_INTERRUPTED = 130

_LOG_FORMAT = "%(asctime)s %(levelname)-7s %(message)s"
_LOG_DATEFMT = "%H:%M:%S"


def _configure_logging(verbose: bool = False, quiet: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.WARNING if quiet else logging.INFO
    logging.basicConfig(level=level, format=_LOG_FORMAT, datefmt=_LOG_DATEFMT, stream=sys.stderr)


def _split_list(value: str | None) -> tuple[str, ...] | None:
    """``"a,b"`` or ``"a b"`` -> ``("a", "b")``. ``None`` stays ``None``."""
    if value is None:
        return None
    items = [item.strip() for item in value.replace(",", " ").split()]
    return tuple(item for item in items if item)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mutationscan",
        description=(
            "Find amino-acid substitutions in target proteins across a cohort of "
            "genomes, and count how often pairs of proteins are mutated together. "
            "No organism, gene or drug is assumed: targets are the reference "
            "files you supply."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "example:\n"
            "  mutationscan run --genomes data/genomes --references data/references \\\n"
            "                   --out data/output --job-name my_run --threads 4\n"
            "\n"
            "Name references <gene>_WT.faa (or <gene>.faa); each filename becomes a target."
        ),
    )
    parser.add_argument("--version", action="version", version=f"mutationscan {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser(
        "run", help="run the full pipeline: extract -> call variants -> count co-occurrence"
    )
    run_parser.add_argument("--genomes", metavar="DIR", help="directory of genome assemblies")
    run_parser.add_argument(
        "--references", metavar="DIR", help="directory of reference proteins; defines the targets"
    )
    run_parser.add_argument(
        "--out", metavar="DIR", help="output root; results go in <out>/<job-name>"
    )
    run_parser.add_argument("--job-name", metavar="NAME", help="subfolder name for this run")
    run_parser.add_argument(
        "--targets",
        metavar="LIST",
        help="comma-separated subset of genes to analyse (default: every reference found)",
    )
    run_parser.add_argument(
        "--taxid",
        metavar="ID",
        help="NCBI taxonomy id; if set, a missing reference is fetched from UniProt",
    )
    run_parser.add_argument(
        "--min-identity", type=float, metavar="PCT", help="reject a protein below this %% identity"
    )
    run_parser.add_argument(
        "--min-coverage", type=float, metavar="PCT", help="reject a protein below this %% coverage"
    )
    run_parser.add_argument(
        "--cooccurrence-level",
        metavar="LIST",
        help="'gene', 'variant', or 'gene,variant' (default: gene)",
    )
    run_parser.add_argument(
        "--min-count", type=int, metavar="N", help="only report pairs co-mutated in >= N genomes"
    )
    run_parser.add_argument("--threads", type=int, metavar="N", help="concurrent tblastn processes")
    run_parser.add_argument(
        "--config", metavar="FILE", help="YAML config to start from; flags override it"
    )
    run_parser.add_argument(
        "--allow-empty",
        action="store_true",
        help=(
            "exit 0 even when a stage produces nothing. Without this a run that "
            "extracts no protein, or passes no pair through the QC gates, exits 3"
        ),
    )
    run_parser.add_argument("-v", "--verbose", action="store_true", help="debug logging")
    run_parser.add_argument("-q", "--quiet", action="store_true", help="warnings and errors only")

    check_parser = subparsers.add_parser(
        "config-check", help="validate a config and print the run layout it resolves to"
    )
    check_parser.add_argument("--config", metavar="FILE", default="config/config.yaml")
    check_parser.add_argument("--json", action="store_true", help="machine-readable output")

    epi_parser = subparsers.add_parser(
        "epistasis",
        help=(
            "epistasis analysis: genome mutation map + ranked epistatic networks "
            "from an existing mutations table"
        ),
    )
    epi_parser.add_argument(
        "--mutations", metavar="CSV", required=True, help="path to mutations.csv"
    )
    epi_parser.add_argument(
        "--qc", metavar="CSV", help="path to variant_qc.csv (optional but recommended)"
    )
    epi_parser.add_argument(
        "--out", metavar="DIR", required=True, help="output directory for epistasis results"
    )
    epi_parser.add_argument(
        "--fdr-threshold",
        type=float,
        metavar="PCT",
        default=0.05,
        help="FDR significance threshold (default: 0.05)",
    )
    epi_parser.add_argument(
        "--min-count",
        type=int,
        metavar="N",
        default=1,
        help="minimum co-occurrence count to retain a pair (default: 1)",
    )
    epi_parser.add_argument("--threads", type=int, metavar="N", default=1)
    epi_parser.add_argument("-v", "--verbose", action="store_true", help="debug logging")
    epi_parser.add_argument("-q", "--quiet", action="store_true", help="warnings and errors only")

    return parser


def _config_from_args(args: argparse.Namespace) -> Config:
    """Config file first (if given), then command-line overrides on top."""
    if args.config:
        config = load_config(args.config)
    else:
        config = Config.from_mapping({})

    levels = _split_list(getattr(args, "cooccurrence_level", None))
    cooccurrence_config = config.cooccurrence
    if levels is not None or args.min_count is not None:
        from .config import VALID_COOCCURRENCE_LEVELS, CooccurrenceConfig

        chosen = levels if levels is not None else cooccurrence_config.levels
        unknown = [level for level in chosen if level not in VALID_COOCCURRENCE_LEVELS]
        if unknown:
            raise ConfigError(
                f"Unknown co-occurrence level(s): {', '.join(unknown)}. "
                f"Valid: {', '.join(VALID_COOCCURRENCE_LEVELS)}."
            )
        cooccurrence_config = CooccurrenceConfig(
            levels=tuple(chosen),
            min_count=(
                args.min_count if args.min_count is not None else cooccurrence_config.min_count
            ),
        )

    return config.with_overrides(
        job_name=args.job_name,
        genomes_dir=Path(args.genomes).expanduser().resolve() if args.genomes else None,
        references_dir=(Path(args.references).expanduser().resolve() if args.references else None),
        output_root=Path(args.out).expanduser().resolve() if args.out else None,
        targets=_split_list(args.targets),
        uniprot_taxid=args.taxid,
        min_identity_percent=args.min_identity,
        min_coverage_percent=args.min_coverage,
        threads=args.threads,
        cooccurrence=cooccurrence_config,
    )


def _report_empty(lines: Sequence[str], allow_empty: bool) -> int:
    """Log an empty-result diagnostic and return the exit code it warrants."""
    log = logger.warning if allow_empty else logger.error
    for line in lines:
        log("%s", line)
    if allow_empty:
        logger.warning("Continuing anyway: --allow-empty was given.")
        return EXIT_OK
    logger.error("Pass --allow-empty to accept an empty result as success.")
    return EXIT_EMPTY_RESULT


def run(config: Config, allow_empty: bool = False) -> int:
    """Execute all four stages in order. Returns a process exit code."""
    config.require_inputs()

    logger.info("MutationScan %s -- run '%s'", __version__, config.job_name)
    logger.info("  genomes:    %s", config.genomes_dir)
    logger.info("  references: %s", config.references_dir)
    logger.info("  output:     %s", config.run_dir)

    # Validate everything before spending a single alignment on it. A swapped
    # --genomes/--references pair costs seconds here and an hour without this.
    report = preflight(config)
    report.log()
    if not report.ok:
        logger.error("Preflight found %d problem(s); nothing was run.", len(report.errors))
        return EXIT_USAGE

    config.run_dir.mkdir(parents=True, exist_ok=True)

    # 1. Targets come from the reference filenames -- this is the only place the
    #    gene list is decided, and it is decided by the user's own directory.
    targets = report.targets
    logger.info("Targets (%d): %s", len(targets), ", ".join(targets))

    references = prepare_references(
        config.references_dir, config.refs_dir, targets, config.uniprot_taxid
    )
    if references.missing:
        logger.error(
            "No reference available for: %s. Supply <gene>_WT.faa in %s, or set "
            "--taxid to fetch from UniProt.",
            ", ".join(references.missing),
            config.references_dir,
        )
        return EXIT_USAGE

    # 2. Extract every target protein from every genome.
    manifest = extract(
        config.genomes_dir,
        config.refs_dir,
        config.proteins_dir,
        references,
        manifest_csv=config.manifest_csv,
        tblastn_binary=config.tblastn_binary,
        threads=config.threads,
    )

    extracted = (
        int((manifest["Status"] == STATUS_EXTRACTED).sum())
        if not manifest.empty and "Status" in manifest.columns
        else 0
    )
    if extracted == 0:
        code = _report_empty(diagnose_empty_extraction(manifest, config), allow_empty)
        if code != EXIT_OK:
            return code

    # 3. Call substitutions against the reference.
    variants = call_variants(
        config.proteins_dir,
        config.refs_dir,
        manifest=manifest,
        min_identity_percent=config.min_identity_percent,
        min_coverage_percent=config.min_coverage_percent,
        mutations_csv=config.mutations_csv,
        qc_csv=config.variant_qc_csv,
    )

    called = (
        int((variants.qc["Status"] == VARIANT_STATUS_CALLED).sum())
        if not variants.qc.empty and "Status" in variants.qc.columns
        else 0
    )
    if called == 0:
        code = _report_empty(diagnose_empty_variants(variants.qc, config), allow_empty)
        if code != EXIT_OK:
            return code

    # 4. Count co-occurrence. Descriptive only.
    #
    #    No emptiness gate here on purpose: a single-target cohort yields zero
    #    pairs, and so does a cohort where nothing happens to be mutated. Both
    #    are legitimate results, not failures.
    counts = cooccurrence(
        variants.mutations,
        variants.qc,
        levels=config.cooccurrence.levels,
        min_count=config.cooccurrence.min_count,
        genes_csv=config.cooccurrence_genes_csv,
        matrix_csv=config.cooccurrence_matrix_csv,
        variants_csv=(
            config.cooccurrence_variants_csv if "variant" in config.cooccurrence.levels else None
        ),
    )

    # 5. Epistasis -- genome mutation maps + statistical pair networks.
    #    Only runs when there are at least two mutations to compare.
    epistasis_networks(
        mutations_csv=config.mutations_csv,
        qc_csv=config.variant_qc_csv,
        output_dir=config.run_dir,
        fdr_threshold=config.epistasis.fdr_threshold,
        min_count=config.epistasis.min_count,
        threads=config.epistasis.threads,
    )

    write_run_summary(
        config,
        references=references,
        manifest=manifest,
        qc=variants.qc,
        mutations=variants.mutations,
        cooccurrence=counts,
        tblastn_version=report.tblastn_version,
    )

    logger.info("Done. Results in %s", config.run_dir)
    return EXIT_OK


def _config_check(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    inputs_present = {
        "genomes_dir": config.genomes_dir.is_dir(),
        "references_dir": config.references_dir.is_dir(),
    }

    # Only worth validating contents once both directories exist -- otherwise
    # preflight would just restate what inputs_present already says.
    report = preflight(config) if all(inputs_present.values()) else None

    layout = {
        "config_file": str(Path(args.config).resolve()),
        "settings": config.to_dict(),
        "run_dir": str(config.run_dir),
        "expected_outputs": [str(path) for path in config.outputs()],
        "inputs_present": inputs_present,
        "preflight": {
            "ran": report is not None,
            "ok": report.ok if report else False,
            "errors": list(report.errors) if report else [],
            "warnings": list(report.warnings) if report else [],
            "genomes": report.n_genomes if report else 0,
            "references": report.n_references if report else 0,
            "targets": list(report.targets) if report else [],
            "tblastn": report.tblastn_version if report else "",
        },
    }

    ready = all(inputs_present.values()) and report is not None and report.ok

    if args.json:
        print(json.dumps(layout, indent=2))
        return EXIT_OK if ready else EXIT_INPUTS_MISSING

    print(f"config: {layout['config_file']}  (valid)")
    print(f"run dir: {config.run_dir}")
    for label, present in inputs_present.items():
        path = getattr(config, label)
        print(f"  {label}: {path} {'[found]' if present else '[MISSING]'}")

    if report is not None:
        print(f"  genomes: {report.n_genomes}")
        print(
            f"  targets ({len(report.targets)}): "
            f"{', '.join(report.targets) if report.targets else '(none found)'}"
        )
        if report.tblastn_version:
            print(f"  tblastn: {report.tblastn_version}")
        for warning in report.warnings:
            print(f"  WARNING: {warning}")
        for error in report.errors:
            print(f"  ERROR: {error}")

    print("expected outputs:")
    for path in config.outputs():
        print(f"  {path}")

    print("ready to run" if ready else "NOT ready to run")
    return EXIT_OK if ready else EXIT_INPUTS_MISSING


def _run_epistasis(args: argparse.Namespace) -> int:
    _configure_logging(verbose=getattr(args, "verbose", False), quiet=getattr(args, "quiet", False))
    try:
        epistasis_networks(
            mutations_csv=args.mutations,
            qc_csv=args.qc,
            output_dir=args.out,
            fdr_threshold=args.fdr_threshold,
            min_count=args.min_count,
            threads=args.threads,
        )
    except (ValueError, OSError) as exc:
        logger.error("Epistasis failed: %s", exc)
        return EXIT_USAGE

    logger.info("Results in %s", args.out)
    return EXIT_OK


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    _configure_logging(verbose=getattr(args, "verbose", False), quiet=getattr(args, "quiet", False))

    try:
        if args.command == "run":
            return run(_config_from_args(args), allow_empty=args.allow_empty)
        if args.command == "config-check":
            return _config_check(args)
        if args.command == "epistasis":
            return _run_epistasis(args)
    except ConfigError as exc:
        logger.error("Configuration error: %s", exc)
        return EXIT_USAGE
    except KeyboardInterrupt:  # pragma: no cover - interactive
        logger.warning("Interrupted")
        return EXIT_INTERRUPTED

    return EXIT_USAGE  # pragma: no cover - argparse enforces a command


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
