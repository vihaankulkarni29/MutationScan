"""Resource-scaling benchmark: wall-clock and peak-RSS curves.

Two figures for the manuscript, both produced from one
run:

* **Panel A** -- cohort size vs. wall-clock time, at a
  fixed thread count. The claim is linear ``O(N)``
  scaling: doubling the cohort doubles the time.
* **Panel B** -- CPU cores vs. peak RSS, at a fixed
  cohort size. The claim is that memory stays flat
  under 2 GiB regardless of thread count, because the
  parallelism is process-level (tblastn) rather than
  in-memory.

Both curves are measured, not asserted. The script
prints the numbers, writes a CSV, and renders the
figure, so a reviewer can reproduce either panel
independently.

Run it as a module::

    python -m benchmarks.scaling --references data/references

The default cohort sizes are small enough to finish in
a couple of minutes; pass ``--sizes 10 50 100 250 500
1000`` for the full manuscript curve (an hour or more,
since tblastn dominates).
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
import time
from collections.abc import Sequence
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from mutation_scan.fasta import read_protein_sequence  # noqa: E402

from .download_data import generate_synthetic  # noqa: E402
from .tools import MutationScanTool  # noqa: E402

logger = logging.getLogger(__name__)

__all__ = ["ScalingPoint", "measure_point", "run_scaling", "main"]


def _read_reference(path: Path) -> str:
    record = read_protein_sequence(path, fallback_id=path.stem)
    return str(record.seq)


class ScalingPoint:
    """One (cohort size, threads) measurement."""

    def __init__(
        self,
        genomes: int,
        threads: int,
        runtime_seconds: float,
        peak_memory_mb: float,
    ) -> None:
        self.genomes = genomes
        self.threads = threads
        self.runtime_seconds = runtime_seconds
        self.peak_memory_mb = peak_memory_mb

    @property
    def genomes_per_hour(self) -> float:
        if self.runtime_seconds <= 0:
            return 0.0
        return self.genomes / self.runtime_seconds * 3600.0

    def as_row(self) -> dict[str, Any]:
        return {
            "genomes": self.genomes,
            "threads": self.threads,
            "runtime_seconds": round(self.runtime_seconds, 3),
            "peak_memory_mb": round(self.peak_memory_mb, 1),
            "genomes_per_hour": round(self.genomes_per_hour, 1),
        }


def measure_point(
    references_dir: Path,
    work_dir: Path,
    genomes: int,
    threads: int,
    seed: int,
) -> ScalingPoint:
    """Generate a synthetic cohort of ``genomes`` and time MutationScan.

    The cohort is regenerated at each size so the measurement
    reflects a real run end to end -- generation included in
    neither the time nor the memory, which is the honest cost
    of the pipeline itself.
    """
    cohort_dir = work_dir / f"cohort_{genomes}"
    generate_synthetic(
        references_dir,
        cohort_dir,
        n_genomes=genomes,
        seed=seed,
    )

    tool = MutationScanTool()
    inputs = {
        "genomes": cohort_dir / "genomes",
        "references": references_dir,
        "lineage_manifest": cohort_dir / "lineage_manifest.tsv",
        "threads": threads,
    }
    output_dir = work_dir / f"run_{genomes}_{threads}"

    started = time.perf_counter()
    output = tool.run(inputs, output_dir)
    runtime = time.perf_counter() - started

    return ScalingPoint(
        genomes=genomes,
        threads=threads,
        runtime_seconds=runtime,
        peak_memory_mb=output.peak_memory_mb,
    )


def run_scaling(
    references_dir: Path,
    work_dir: Path,
    sizes: Sequence[int],
    thread_counts: Sequence[int],
    seed: int,
) -> list[ScalingPoint]:
    """Measure every (size, threads) combination.

    Panel A uses the first thread count across all sizes;
    Panel B uses the first size across all thread counts.
    Every combination is measured so both panels come from
    one run and the CSV is complete.
    """
    points: list[ScalingPoint] = []
    for genomes in sizes:
        for threads in thread_counts:
            logger.info("measuring: %d genomes, %d threads", genomes, threads)
            point = measure_point(references_dir, work_dir, genomes, threads, seed)
            points.append(point)
            logger.info(
                "  -> %.1fs, %.0f MiB peak RSS",
                point.runtime_seconds,
                point.peak_memory_mb,
            )
    return points


def _write_csv(points: Sequence[ScalingPoint], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "genomes",
                "threads",
                "runtime_seconds",
                "peak_memory_mb",
                "genomes_per_hour",
            ],
        )
        writer.writeheader()
        for point in points:
            writer.writerow(point.as_row())


def _write_json(points: Sequence[ScalingPoint], path: Path) -> None:
    document = [point.as_row() for point in points]
    with path.open("w", encoding="utf-8") as handle:
        json.dump(document, handle, indent=2)


def _render_figure(
    points: Sequence[ScalingPoint],
    sizes: Sequence[int],
    thread_counts: Sequence[int],
    path: Path,
) -> None:
    """Render the two-panel manuscript figure."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib unavailable; skipping the figure")
        return

    by_size_threads = {(point.genomes, point.threads): point for point in points}

    fig, (ax_time, ax_mem) = plt.subplots(1, 2, figsize=(11, 4.2))

    # Panel A: cohort size vs wall-clock, at the first thread count.
    panel_threads = thread_counts[0]
    x = list(sizes)
    y = [
        by_size_threads[(size, panel_threads)].runtime_seconds
        if (size, panel_threads) in by_size_threads
        else float("nan")
        for size in x
    ]
    ax_time.plot(x, y, marker="o", linewidth=2)
    ax_time.set_xlabel("cohort size (genomes)")
    ax_time.set_ylabel("wall-clock time (s)")
    ax_time.set_title(f"Panel A: scaling at {panel_threads} threads")
    ax_time.grid(True, alpha=0.3)

    # Panel B: cores vs peak RSS, at the first cohort size.
    panel_size = sizes[0]
    x = list(thread_counts)
    y = [
        by_size_threads[(panel_size, threads)].peak_memory_mb
        if (panel_size, threads) in by_size_threads
        else float("nan")
        for threads in x
    ]
    ax_mem.plot(x, y, marker="s", linewidth=2, color="tab:orange")
    ax_mem.set_xlabel("CPU cores")
    ax_mem.set_ylabel("peak RSS (MiB)")
    ax_mem.set_title(f"Panel B: memory at {panel_size} genomes")
    ax_mem.grid(True, alpha=0.3)

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
        default="benchmarks/results/scaling",
        help="where the CSV, JSON and figure are written",
    )
    parser.add_argument(
        "--sizes",
        type=int,
        nargs="+",
        default=[10, 20, 40, 80],
        help="cohort sizes for Panel A (default: 10 20 40 80)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        nargs="+",
        default=[2, 4, 8],
        help="thread counts for Panel B (default: 2 4 8)",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)-7s %(message)s",
    )

    references_dir = Path(args.references)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    points = run_scaling(
        references_dir,
        output_dir,
        args.sizes,
        args.threads,
        args.seed,
    )

    _write_csv(points, output_dir / "scaling.csv")
    _write_json(points, output_dir / "scaling.json")
    _render_figure(points, args.sizes, args.threads, output_dir / "scaling.png")

    print(f"wrote {output_dir / 'scaling.csv'}")
    print(f"wrote {output_dir / 'scaling.json'}")
    print(f"wrote {output_dir / 'scaling.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
