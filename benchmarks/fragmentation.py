"""Fragmentation stress test: mutation recall vs assembly N50.

A gene split across contig boundaries is found only on the
contig carrying the best tblastn HSP; mutations on other
fragments are lost. This is an annotation-independent limit
of tblastn-based extraction.

Run as a module::

    python -m benchmarks.fragmentation --references data/references

Produces a CSV, JSON and a two-panel figure.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import random
import sys
import time
from collections.abc import Sequence
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from mutation_scan.fasta import read_protein_sequence  # noqa: E402

from .amr_benchmark import GroundTruth  # noqa: E402
from .download_data import generate_synthetic  # noqa: E402
from .tools import MutationScanTool  # noqa: E402

logger = logging.getLogger(__name__)

__all__ = ["FragmentPoint", "fragment_genomes", "run_fragmentation", "main"]

#: Inter-contig gap (N bases) written between pieces so tblastn
#: cannot join them across the boundary.
_CONTIG_GAP = "N" * 100

#: Default fragmentation levels: number of contigs per gene.
_DEFAULT_FRAGMENTS = [1, 2, 4, 8]


class FragmentPoint:
    """One (fragments per gene, mutation recall) measurement."""

    def __init__(
        self,
        fragments_per_gene: int,
        runtime_seconds: float,
        peak_memory_mb: float,
        sensitivity: float,
        precision: float,
        f1: float,
        novel_recall: float | None,
    ) -> None:
        self.fragments_per_gene = fragments_per_gene
        self.runtime_seconds = runtime_seconds
        self.peak_memory_mb = peak_memory_mb
        self.sensitivity = sensitivity
        self.precision = precision
        self.f1 = f1
        self.novel_recall = novel_recall

    def as_row(self) -> dict[str, Any]:
        return {
            "fragments_per_gene": self.fragments_per_gene,
            "runtime_seconds": round(self.runtime_seconds, 3),
            "peak_memory_mb": round(self.peak_memory_mb, 1),
            "sensitivity": round(self.sensitivity, 3),
            "precision": round(self.precision, 3),
            "f1": round(self.f1, 3),
            "novel_recall": (
                round(self.novel_recall, 3) if self.novel_recall is not None else "N/A"
            ),
        }


def _read_reference(path: Path) -> str:
    record = read_protein_sequence(path, fallback_id=path.stem)
    return str(record.seq)


def _write_fasta_record(handle: Any, header: str, seq: str) -> None:
    handle.write(f">{header}\n")
    for i in range(0, len(seq), 60):
        handle.write(seq[i : i + 60] + "\n")


def fragment_genomes(
    source_dir: Path,
    output_dir: Path,
    fragments_per_gene: int,
    seed: int,
) -> None:
    """Fragment every genome in ``source_dir`` into ``fragments_per_gene`` pieces.

    Each gene's CDS is cut at random internal positions. The resulting
    contigs are named ``<accession>_c{idx}`` and joined with N gaps
    so tblastn cannot merge them. Ground-truth mutations are NOT
    modified; they remain on the original coordinate system.
    """
    rng = random.Random(seed)
    output_dir.mkdir(parents=True, exist_ok=True)

    for src_path in sorted(source_dir.glob("*.fna")):
        accession = src_path.stem
        # Parse the original multi-record FASTA: one record per gene.
        gene_cds: dict[str, str] = {}
        with src_path.open() as handle:
            current_gene: str | None = None
            current_seq: list[str] = []
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_gene is not None:
                        gene_cds[current_gene] = "".join(current_seq)
                    current_gene = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_gene is not None:
                gene_cds[current_gene] = "".join(current_seq)

        # Write fragmented version.
        dst_path = output_dir / f"{accession}.fna"
        with dst_path.open("w") as handle:
            for gene, cds in gene_cds.items():
                if fragments_per_gene == 1:
                    _write_fasta_record(handle, f"{accession}_{gene}", cds)
                else:
                    # Cut at random internal points (not at codons, so we
                    # split a gene mid-sequence). This models a contig break
                    # anywhere in the CDS.
                    cut_points = sorted(rng.sample(range(1, len(cds)), fragments_per_gene - 1))
                    pieces = []
                    start = 0
                    for cut in cut_points:
                        pieces.append(cds[start:cut])
                        start = cut
                    pieces.append(cds[start:])

                    for idx, piece in enumerate(pieces, 1):
                        header = f"{accession}_{gene}_c{idx}"
                        _write_fasta_record(handle, header, piece)
                        if idx < len(pieces):
                            handle.write(f"{_CONTIG_GAP}\n")


def _load_ground_truth(path: Path) -> GroundTruth:
    """Read a ``ground_truth.json`` written by the synthetic generator."""
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


def run_fragmentation(
    references_dir: Path,
    work_dir: Path,
    fragment_levels: Sequence[int],
    seed: int,
    genomes: int,
    threads: int,
) -> list[FragmentPoint]:
    """Measure recall at each fragmentation level.

    Steps:
    1. Generate a base cohort with ``generate_synthetic``.
    2. For each fragmentation level, re-write genomes with ``fragment_genomes``.
    3. Run MutationScan on the fragmented cohort and score against the
       original ground truth.
    """
    # 1. Base cohort.
    logger.info("generating base synthetic cohort (%d genomes)", genomes)
    base_dir = work_dir / "base"
    generate_synthetic(
        references_dir,
        base_dir,
        n_genomes=genomes,
        seed=seed,
    )
    truth_path = base_dir / "ground_truth.json"
    ground_truth = _load_ground_truth(truth_path)
    logger.info(
        "ground truth: %d mutations, %d epistatic pairs",
        len(ground_truth.mutations),
        len(ground_truth.epistasis),
    )

    # 2. Fragment and benchmark each level.
    points: list[FragmentPoint] = []
    tool = MutationScanTool()

    for frag in fragment_levels:
        logger.info("fragmenting to %d contigs/gene", frag)
        frag_dir = work_dir / f"frag_{frag}"
        fragment_genomes(base_dir / "genomes", frag_dir / "genomes", frag, seed)

        inputs = {
            "genomes": frag_dir / "genomes",
            "references": references_dir,
            "lineage_manifest": base_dir / "lineage_manifest.tsv",
            "threads": threads,
        }
        output_dir = frag_dir / "run"

        started = time.perf_counter()
        output = tool.run(inputs, output_dir)
        runtime = time.perf_counter() - started

        metrics = compute_mutation_metrics(ground_truth, output)
        point = FragmentPoint(
            fragments_per_gene=frag,
            runtime_seconds=runtime,
            peak_memory_mb=output.peak_memory_mb,
            sensitivity=metrics.sensitivity,
            precision=metrics.precision,
            f1=metrics.f1,
            novel_recall=metrics.novel_recall,
        )
        points.append(point)
        logger.info(
            "  -> sens=%.3f prec=%.3f f1=%.3f time=%.1fs mem=%.0fMiB",
            point.sensitivity,
            point.precision,
            point.f1,
            point.runtime_seconds,
            point.peak_memory_mb,
        )

    return points


# Re-export for inline use without importing amr_benchmark.
from .amr_benchmark import compute_mutation_metrics  # noqa: E402


def _write_csv(points: Sequence[FragmentPoint], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "fragments_per_gene",
                "runtime_seconds",
                "peak_memory_mb",
                "sensitivity",
                "precision",
                "f1",
                "novel_recall",
            ],
        )
        writer.writeheader()
        for point in points:
            writer.writerow(point.as_row())


def _write_json(points: Sequence[FragmentPoint], path: Path) -> None:
    document = [point.as_row() for point in points]
    with path.open("w", encoding="utf-8") as handle:
        json.dump(document, handle, indent=2)


def _render_figure(
    points: Sequence[FragmentPoint],
    fragment_levels: Sequence[int],
    path: Path,
) -> None:
    """Render a recall-vs-fragmentation figure."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib unavailable; skipping the figure")
        return

    fig, (ax_recall, ax_time) = plt.subplots(1, 2, figsize=(11, 4.2))

    # Panel A: fragments per gene vs sensitivity/F1.
    x = [point.fragments_per_gene for point in points]
    ax_recall.plot(
        x,
        [point.sensitivity for point in points],
        marker="o",
        linewidth=2,
        label="sensitivity",
    )
    ax_recall.plot(
        x,
        [point.f1 for point in points],
        marker="s",
        linewidth=2,
        label="F1",
        color="tab:orange",
    )
    ax_recall.set_xscale("log", base=2)
    ax_recall.set_xlabel("contigs per gene")
    ax_recall.set_ylabel("score")
    ax_recall.set_title("Panel A: recall vs fragmentation")
    ax_recall.grid(True, alpha=0.3)
    ax_recall.legend()

    # Panel B: fragments per gene vs wall-clock.
    ax_time.plot(
        x,
        [point.runtime_seconds for point in points],
        marker="^",
        linewidth=2,
        color="tab:green",
    )
    ax_time.set_xscale("log", base=2)
    ax_time.set_xlabel("contigs per gene")
    ax_time.set_ylabel("wall-clock time (s)")
    ax_time.set_title("Panel B: runtime vs fragmentation")
    ax_time.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--references",
        default="data/references",
        help="directory of <gene>_WT.faa reference proteins",
    )
    parser.add_argument(
        "--output",
        default="benchmarks/results/fragmentation",
        help="where the CSV, JSON and figure are written",
    )
    parser.add_argument(
        "--fragments",
        type=int,
        nargs="+",
        default=_DEFAULT_FRAGMENTS,
        help="contigs per gene to test (default: 1 2 4 8)",
    )
    parser.add_argument("--genomes", type=int, default=20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)-7s %(message)s",
    )

    references_dir = Path(args.references)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    points = run_fragmentation(
        references_dir,
        output_dir,
        args.fragments,
        args.seed,
        args.genomes,
        args.threads,
    )

    _write_csv(points, output_dir / "fragmentation.csv")
    _write_json(points, output_dir / "fragmentation.json")
    _render_figure(points, args.fragments, output_dir / "fragmentation.png")

    print(f"wrote {output_dir / 'fragmentation.csv'}")
    print(f"wrote {output_dir / 'fragmentation.json'}")
    print(f"wrote {output_dir / 'fragmentation.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
