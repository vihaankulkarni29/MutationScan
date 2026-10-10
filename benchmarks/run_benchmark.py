"""Command-line entry point for the AMR benchmark.

    python benchmarks/run_benchmark.py \
        --tools mutationscan,amrfinderplus \
        --dataset synthetic --references data/references \
        --output benchmarks/results

The entry point owns the workflow -- pick a dataset,
build the ground truth, pick the tools, run the
harness, write the report -- and delegates the
mechanics to :mod:`amr_benchmark`, :mod:`download_data`
and :mod:`report`. It is the only benchmark module
that a user is expected to invoke directly.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from .amr_benchmark import BenchmarkHarness, GroundTruth
from .download_data import generate_synthetic
from .report import write_report
from .tools import build_tool

__all__ = ["main"]

logger = logging.getLogger(__name__)

_ALL_TOOLS = ("mutationscan", "amrfinderplus", "mykrobe", "kleborate", "staramr")


def _load_ground_truth(path: Path) -> GroundTruth:
    """Read a ``ground_truth.json`` written by the
    synthetic generator or curated by hand."""
    with path.open(encoding="utf-8") as handle:
        document = json.load(handle)

    mutations = {tuple(record) for record in document.get("mutations", [])}
    epistasis = {frozenset(record) for record in document.get("epistasis", [])}
    dn_ds: dict[str, float] = {}
    for gene, value in document.get("dn_ds", {}).items():
        dn_ds[gene] = float("inf") if value == "inf" else float(value)
    lineages = dict(document.get("lineages", {}))
    known = {tuple(record) for record in document.get("known_mutations", [])}
    negatives = {tuple(record) for record in document.get("negatives", [])}
    return GroundTruth(
        mutations=mutations,
        epistasis=epistasis,
        dn_ds=dn_ds,
        lineages=lineages,
        known_mutations=known,
        negatives=negatives,
    )


def _build_inputs(
    args: argparse.Namespace,
    output_dir: Path,
) -> dict[str, Any]:
    inputs: dict[str, Any] = {}
    if args.genomes:
        inputs["genomes"] = Path(args.genomes)
    if args.references:
        inputs["references"] = Path(args.references)
    if args.lineage_manifest:
        inputs["lineage_manifest"] = Path(args.lineage_manifest)
    if args.threads:
        inputs["threads"] = args.threads
    if args.organism:
        inputs["organism"] = args.organism
    if args.species:
        inputs["species"] = args.species
    return inputs


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Benchmark MutationScan against AMR tools.")
    parser.add_argument(
        "--tools",
        default=",".join(_ALL_TOOLS),
        help=f"comma-separated tools to run (default: all: {','.join(_ALL_TOOLS)})",
    )
    parser.add_argument(
        "--dataset",
        default="synthetic",
        help="'synthetic', or a path to a ground_truth.json",
    )
    parser.add_argument(
        "--references",
        default="data/references",
        help="directory of <gene>_WT.faa reference proteins",
    )
    parser.add_argument(
        "--genomes",
        help="directory of genome assemblies (real data)",
    )
    parser.add_argument(
        "--lineage-manifest",
        help="TSV/CSV mapping Accession -> Lineage",
    )
    parser.add_argument(
        "--output",
        default="benchmarks/results",
        help="where results and the report are written",
    )
    parser.add_argument("--genomes-count", type=int, default=50)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--threads", type=int)
    parser.add_argument("--organism", help="AMRFinderPlus organism")
    parser.add_argument("--species", help="Mykrobe species")
    parser.add_argument("-v", "--verbose", action="store_true", help="debug logging")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)-7s %(message)s",
    )

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Ground truth: generate a synthetic cohort, or load
    #    a curated one.
    genomes_dir: Path | None = None
    lineage_manifest: Path | None = None
    if args.dataset == "synthetic":
        synthetic_dir = output_dir / "synthetic"
        generate_synthetic(
            args.references,
            synthetic_dir,
            n_genomes=args.genomes_count,
            seed=args.seed,
        )
        truth_path = synthetic_dir / "ground_truth.json"
        genomes_dir = synthetic_dir / "genomes"
        lineage_manifest = synthetic_dir / "lineage_manifest.tsv"
    else:
        truth_path = Path(args.dataset)
        genomes_dir = Path(args.genomes) if args.genomes else None
        lineage_manifest = Path(args.lineage_manifest) if args.lineage_manifest else None

    ground_truth = _load_ground_truth(truth_path)
    logger.info(
        "ground truth: %d mutations, %d epistatic pairs, %d genes with dN/dS, %d lineages",
        len(ground_truth.mutations),
        len(ground_truth.epistasis),
        len(ground_truth.dn_ds),
        len(ground_truth.lineages),
    )

    # 2. Inputs each tool needs.
    inputs = _build_inputs(args, output_dir)
    if genomes_dir is not None:
        inputs["genomes"] = genomes_dir
    if lineage_manifest is not None and "lineage_manifest" not in inputs:
        inputs["lineage_manifest"] = lineage_manifest

    # 3. Tools. A tool whose binary is absent is skipped,
    #    not fatal: the benchmark still runs for the rest,
    #    and a missing tool is never scored as a zero.
    tools = []
    for name in args.tools.split(","):
        name = name.strip()
        if not name:
            continue
        try:
            tool = build_tool(name)
        except ValueError as exc:
            logger.error("%s", exc)
            return 2
        if not tool.available:
            logger.warning("skipping %s: %s is not on PATH", tool.name, tool.binary)
            continue
        tools.append(tool)

    if not tools:
        logger.error("no tools available to benchmark")
        return 2

    # 4. Run and score.
    harness = BenchmarkHarness(ground_truth, tools, inputs, work_dir=output_dir / "runs")
    results = harness.run()

    # 5. Report.
    paths = write_report(results, output_dir)
    for result in results:
        logger.info(
            "%s: sensitivity=%s precision=%s epistasis F1=%s dN/dS r=%s lineage=%s",
            result.tool,
            f"{result.mutations.sensitivity:.3f}",
            f"{result.mutations.precision:.3f}",
            f"{result.epistasis.f1:.3f}",
            f"{result.dn_ds.pearson:.3f}",
            f"{result.lineage.concordance:.3f}",
        )
    print(f"report: {paths['markdown']}")
    print(f"report: {paths['html']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
