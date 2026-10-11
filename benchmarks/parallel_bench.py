"""Measure the parallel layer, so its thresholds are evidence and not folklore.

    python benchmarks/parallel_bench.py

Reports, for a synthetic CPU-bound workload:

* pool start-up cost, isolated;
* a scaling curve, so the point where parallelism starts paying is visible;
* whether the parallel result equals the serial one, at every size.

Two things this exists to protect:

1. **Reproducibility.** Every size is checked for output equality, not just
   timed. A speedup on a path that quietly reorders or drops results is not a
   speedup.
2. **Honest thresholds.** ``resolve_chunk_size`` and
   ``epistasis._SERIAL_PAIR_LIMIT`` are empirical. The numbers here are what
   justifies re-tuning them on new hardware, and ``--verify-thresholds`` fails
   loudly when the machine can no longer support them.

Timings are inherently noisy -- a shared or thermally throttled machine will
report nonsense -- so this prints the numbers and states a verdict rather than
asserting a speedup. Only the correctness checks and the threshold check are
hard failures.
"""

from __future__ import annotations

import argparse
import statistics
import sys
import time
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from mutation_scan.parallel import (  # noqa: E402
    WORKLOAD_PROCESS,
    available_cores,
    default_mp_context,
    parallel_map,
    resolve_chunk_size,
)

#: Inner-loop count for the synthetic task. Sized so one item costs a few
#: milliseconds -- in the same ballpark as one Biopython global alignment.
DEFAULT_WORK = 400_000

SIZES: tuple[int, ...] = (8, 32, 64, 128, 256, 512)


# Module level: a process pool must be able to import these by name.


def burn(value: int, work: int = DEFAULT_WORK) -> int:
    """CPU-bound Python, the kind a thread pool cannot help with."""
    total = 0
    for index in range(work):
        total += index % 7
    return total + value


def trivial(value: int) -> int:
    return value


@dataclass(frozen=True)
class Measurement:
    label: str
    serial: float
    parallel: float

    @property
    def speedup(self) -> float:
        return self.serial / self.parallel if self.parallel > 0 else float("inf")

    @property
    def parallel_slower(self) -> bool:
        return self.parallel > self.serial


def _time(
    fn: Callable[[int], int],
    size: int,
    threads: int,
    workload: str,
    repeats: int,
) -> tuple[float, list[int]]:
    """Best-of-N wall time, plus the results, so correctness can be checked too.

    Best-of rather than mean: the contention this measures is exactly the kind
    that produces occasional slow runs, and the best run is the one that
    represents the code rather than the scheduler.
    """
    best = float("inf")
    results: list[int] = []
    for _ in range(max(1, repeats)):
        start = time.perf_counter()
        results = list(
            parallel_map(fn, range(size), threads=threads, workload=workload, total=size)
        )
        best = min(best, time.perf_counter() - start)
    return best, results


def measure_startup(repeats: int) -> float:
    """Cost of bringing a pool up, with no real work in it.

    This is the fixed price that every parallel stage pays, and the reason a
    small cohort can be slower in parallel than serial.
    """
    samples = []
    for _ in range(max(1, repeats)):
        start = time.perf_counter()
        list(parallel_map(trivial, range(8), threads=4, workload=WORKLOAD_PROCESS, total=8))
        samples.append(time.perf_counter() - start)
    return statistics.median(samples)


def run(sizes: Sequence[int], workers: int, repeats: int) -> list[Measurement]:
    measurements = []
    for size in sizes:
        serial_time, serial = _time(burn, size, 1, WORKLOAD_PROCESS, repeats)
        parallel_time, par = _time(burn, size, workers, WORKLOAD_PROCESS, repeats)
        if serial != par:
            index = next(
                (i for i, (a, b) in enumerate(zip(serial, par, strict=True)) if a != b),
                None,
            )
            raise SystemExit(
                f"CORRECTNESS FAILURE at n={size}: parallel output differs from serial. "
                f"First difference at index {index}."
            )
        measurements.append(Measurement(f"n={size}", serial_time, parallel_time))
    return measurements


def print_report(measurements: Sequence[Measurement], workers: int, startup: float) -> None:
    print(f"cores            : {available_cores()}")
    print(f"start method     : {default_mp_context().get_start_method()}")
    print(f"workers compared : {workers}")
    print(f"pool start-up    : {startup * 1000:.0f} ms (4 workers, negligible work)")
    print(f"start-up / task  : {startup / max(1, workers) * 1000:.0f} ms per worker")
    print()
    print(f"{'workload':>10} {'serial':>10} {'parallel':>10} {'speedup':>9}  verdict")
    print("-" * 56)
    for measurement in measurements:
        if measurement.parallel_slower:
            verdict = "slower - below start-up crossover"
        elif measurement.speedup >= workers * 0.8:
            verdict = "scaling well"
        else:
            verdict = "partial"
        print(
            f"{measurement.label:>10} "
            f"{measurement.serial:>9.2f}s "
            f"{measurement.parallel:>9.2f}s "
            f"{measurement.speedup:>8.2f}x  {verdict}"
        )
    print()
    print("All sizes produced output identical to the serial run.")


def measure_fisher_cost(repeats: int = 200) -> float:
    """Cheapest plausible cost of one epistasis pair.

    ``fisher_exact`` is the dominant term in a pair test, and its cost varies by
    an order of magnitude with table shape (~160us balanced and small, ~1.5ms for
    a rare co-mutation in a large cohort). The *cheapest* measurement is used,
    which is the least favourable case for the serial threshold being set high
    enough -- if the threshold works here, it works everywhere.
    """
    from scipy.stats import fisher_exact

    table = [[125, 125], [125, 125]]
    fisher_exact(table)  # warm up
    start = time.perf_counter()
    for _ in range(max(1, repeats)):
        fisher_exact(table)
    return (time.perf_counter() - start) / repeats


def verify_thresholds(measurements: Sequence[Measurement], startup: float) -> int:
    """Check that the empirical thresholds still hold on this machine.

    Guards a regression that is invisible in output and only shows up as
    "parallelism is mysteriously slow": chunk sizing that collapses a small
    cohort into one batch, stranding every worker but one.
    """
    failures = []

    print("threshold check")
    print(f"  pool start-up        : {startup * 1000:.0f} ms")
    print(f"  serial/pool crossover: {startup:.3f}s of work")

    smallest_useful = next((m for m in measurements if not m.parallel_slower), None)
    if smallest_useful is not None:
        print(f"  break-even observed  : {smallest_useful.label} ({smallest_useful.serial:.2f}s)")
    else:
        print("  break-even observed  : not reached in these sizes")

    for measurement in measurements:
        size = int(measurement.label.split("=")[1])
        chunk = resolve_chunk_size(total=size, workers=4)
        batches = -(-size // chunk)
        if batches < 2 and measurement.parallel_slower:
            failures.append(
                f"chunk sizing strands workers at n={size}: chunk size {chunk} "
                f"yields {batches} batch(es)"
            )

    # Epistasis used to be guarded by a serial/pool pair threshold. That check
    # is gone with the pool: the stage now resolves counts with one Gram matrix
    # per universe and runs one test per distinct 2x2 table, which took the
    # 570-genome run from 146,070 tests to ~1,700. What is worth guarding now is
    # that no pool creeps back in -- it was measurably a net loss, and nothing in
    # the output would reveal it if someone reintroduced one.
    import mutation_scan.epistasis as epistasis

    fisher = measure_fisher_cost()
    print(f"  fisher_exact         : {fisher * 1e6:.0f}us per distinct table")

    if hasattr(epistasis, "parallel_map"):
        failures.append(
            "epistasis imports parallel_map again: the stage resolved every pair "
            "with one Fisher test per distinct table, so a pool only adds spawn cost"
        )
    else:
        print("  epistasis            : no process pool (exact counts + memoised tests)")

    for failure in failures:
        print(f"  FAIL: {failure}")
    if failures:
        return 1
    print("  OK: empirical thresholds hold on this machine")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=4, help="worker count to compare against")
    parser.add_argument("--repeats", type=int, default=2, help="best-of-N per measurement")
    parser.add_argument(
        "--verify-thresholds",
        action="store_true",
        help="exit non-zero if the empirical thresholds no longer hold",
    )
    args = parser.parse_args(argv)

    print("warming up (pool start-up is measured cold on the first run)...\n")
    startup = measure_startup(2)
    measurements = run(SIZES, args.workers, args.repeats)
    print_report(measurements, args.workers, startup)
    print()

    if args.verify_thresholds:
        return verify_thresholds(measurements, startup)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
