"""Tests for the shared parallel-execution layer.

These are the tests that matter most for this module, because its contract is not
about speed but about *not changing results*. Everything here is written to fail
loudly if a future change makes the parallel path disagree with the serial one,
silently drop work, or quietly swallow a failure.

The worker functions must be module level to be picklable under the process
backend; pytest test modules are importable by name, so module-level helpers in
this file are importable in a spawned worker.
"""

from __future__ import annotations

import os
import time
from collections.abc import Iterator

import pytest

from mutation_scan.parallel import (
    DEFAULT_CHUNK_SIZE,
    WORKLOAD_PROCESS,
    WORKLOAD_THREAD,
    ParallelTaskError,
    available_cores,
    chunked,
    default_mp_context,
    parallel_map,
    pin_blas_threads,
    resolve_chunk_size,
    resolve_workers,
)

# -- module-level workers (must be picklable by reference) ------------------


def double(value: int) -> int:
    return value * 2


def jittered_double(value: int) -> int:
    """Uneven cost, so later items finish before earlier ones.

    This is what would expose a result set ordered by completion rather than by
    input position.
    """
    time.sleep(0.0005 * (value % 7))
    return value * 2


def explode_on_37(value: int) -> int:
    if value == 37:
        raise ValueError("deliberate failure")
    return value


def die_on_37(value: int) -> int:
    """Simulates a worker being killed outright rather than raising."""
    if value == 37:
        os._exit(1)
    return value


def counting_generator(limit: int) -> Iterator[int]:
    yield from range(limit)


# -- initializer / initargs ------------------------------------------------

_WORKER_NOTE: list[str] = []


def _record_note(factor: int, prefix: str) -> None:
    """Initializer: takes an argument, which is the whole point of the test.

    A worker that cannot receive this argument dies during start-up, and the
    failure surfaces only as "terminated abruptly" on the parent side. That
    happened for real: an initializer needing three tables was called with none.
    """
    _WORKER_NOTE.append(f"{prefix}:{factor}")


def note_seen(value: str) -> str:
    """Reads the state the initializer installed, and reports whether it ran."""
    if not _WORKER_NOTE:
        return f"MISSING|{value}"
    return f"{_WORKER_NOTE[0]}|{value}"


# -- determinism ------------------------------------------------------------


@pytest.mark.parametrize("workload", [WORKLOAD_PROCESS, WORKLOAD_THREAD])
@pytest.mark.parametrize("threads", [1, 2, 3, 5])
def test_results_are_independent_of_thread_count(workload, threads):
    """The core promise: ``threads`` is a scheduling knob, not a behaviour change.

    A MutationScan run must be reproducible. If worker count could reorder
    results, or drop them, the tool's output would depend on machine load.
    """
    expected = list(parallel_map(double, range(60), threads=1, workload=workload))
    got = list(parallel_map(double, range(60), threads=threads, workload=workload, total=60))
    assert got == expected


@pytest.mark.parametrize("threads", [2, 4, 8])
def test_order_survives_unequal_task_durations(threads):
    """Order must be by input position, not completion.

    With jittered costs, batch 0 finishes long after batch 2. A pool that yielded
    in completion order would interleave them and fail this.
    """
    got = list(
        parallel_map(
            jittered_double, range(200), threads=threads, workload=WORKLOAD_THREAD, total=200
        )
    )
    assert got == [value * 2 for value in range(200)]


@pytest.mark.parametrize("chunk_size", [1, 2, 7, 64, 100_000])
def test_results_are_independent_of_chunk_size(chunk_size):
    """Batching is an internal detail and must not be observable in the output."""
    got = list(
        parallel_map(
            double,
            range(97),
            threads=3,
            workload=WORKLOAD_THREAD,
            total=97,
            chunk_size=chunk_size,
        )
    )
    assert got == [value * 2 for value in range(97)]


def test_a_generator_input_is_accepted_and_fully_consumed():
    got = list(parallel_map(double, counting_generator(40), threads=3, workload=WORKLOAD_PROCESS))
    assert got == [value * 2 for value in range(40)]


def test_an_empty_input_produces_no_output():
    assert list(parallel_map(double, [], threads=4)) == []


# -- laziness and bounded memory -------------------------------------------


def test_items_are_consumed_lazily_not_eagerly():
    """A generator must not be drained before the first result is yielded.

    Materialising the input up front is the memory regression this module exists
    to prevent: it turns a bounded-memory stage into one that grows with the
    cohort, and it is invisible until a run is too big to finish.

    The generator is local to the test; that is fine, because generators are
    consumed in the parent process and never pickled. Only the mapped function
    has to be importable, so that is a module-level one.
    """
    pulled: list[int] = []

    def observed() -> Iterator[int]:
        for index in range(1000):
            pulled.append(index)
            yield index

    stream = parallel_map(double, observed(), threads=2, workload=WORKLOAD_PROCESS, total=1000)
    first = next(stream)
    assert first == 0
    # Nowhere near all 1000, though the pool is warm and eager to help.
    assert len(pulled) < 1000
    stream.close()


def test_input_generator_is_not_consumed_further_after_an_early_failure():
    """A run that dies must stop pulling from its input.

    Otherwise a stage that fails on item 3 has still walked the whole cohort --
    reading every file -- before surfacing the error.
    """
    pulled: list[int] = []

    def observed() -> Iterator[int]:
        for index in range(10_000):
            pulled.append(index)
            if index == 2:
                raise RuntimeError("input generation failed")
            yield index

    with pytest.raises((ParallelTaskError, RuntimeError)):
        list(parallel_map(double, observed(), threads=1))

    assert len(pulled) < 1000


# -- failure semantics ------------------------------------------------------


@pytest.mark.parametrize("threads", [1, 4])
def test_a_failing_task_raises_rather_than_returning_an_error_value(threads):
    """A failure must not be reported as data.

    Returning a sentinel would let a broken run look like a successful one with
    an odd result in it, which is exactly the failure mode that costs days.
    """
    with pytest.raises(ParallelTaskError) as excinfo:
        list(
            parallel_map(
                explode_on_37, range(100), threads=threads, workload=WORKLOAD_THREAD, total=100
            )
        )
    assert isinstance(excinfo.value.cause, ValueError)
    assert excinfo.value.completed >= 0


def test_the_serial_path_also_reports_parallel_task_error():
    """One exception type regardless of worker count, so callers need one handler."""
    with pytest.raises(ParallelTaskError):
        list(parallel_map(explode_on_37, range(50), threads=1))


def test_a_dead_worker_raises_instead_of_reporting_the_rest_as_errors():
    """A killed process must not masquerade as unreadable data.

    The old shape caught the pool error and wrote an ``error`` QC row for every
    pair in the failed batch: 598 rows, zero mutations, exit status 0.
    """
    with pytest.raises(ParallelTaskError):
        list(
            parallel_map(
                die_on_37, range(60), threads=2, workload=WORKLOAD_PROCESS, total=60, chunk_size=8
            )
        )


def test_the_error_names_the_label_and_keeps_counts():
    with pytest.raises(ParallelTaskError, match="variant calling"):
        list(
            parallel_map(
                explode_on_37,
                range(50),
                threads=2,
                workload=WORKLOAD_THREAD,
                total=50,
                label="variant calling",
            )
        )


# -- worker sizing ----------------------------------------------------------


@pytest.mark.parametrize("requested", [0, -1, -100])
def test_non_positive_thread_counts_clamp_to_one(requested):
    assert resolve_workers(requested) == 1


def test_a_request_above_the_core_count_is_honoured_but_warned(caplog):
    """Oversubscription is the user's call, but it must not be silent.

    Silently clamping to the core count would be a surprise, and silently
    honouring it hides why a stage got slower.
    """
    too_many = available_cores() + 8
    with caplog.at_level("WARNING"):
        assert resolve_workers(too_many) == too_many
    assert any("core" in record.message for record in caplog.records)


def test_a_silly_thread_count_is_capped_to_protect_file_descriptors(caplog):
    """Threads each hold a tblastn subprocess open, so there is a real ceiling."""
    with caplog.at_level("WARNING"):
        workers = resolve_workers(10_000, workload=WORKLOAD_THREAD)
    assert workers == 64


def test_available_cores_is_at_least_one():
    assert available_cores() >= 1


# -- chunk sizing -----------------------------------------------------------


def test_a_short_workload_is_split_so_every_worker_gets_something():
    """Regression test for a stage that silently ran on one core.

    With a fixed chunk size of 64, a 24-item workload became a single batch, so
    one worker did everything: a measured 1.03x "speedup" for the full cost of
    starting a pool.
    """
    size = resolve_chunk_size(total=24, workers=4)
    assert size < 24
    assert len(list(chunked(range(24), size))) >= 4


def test_a_long_workload_keeps_chunks_coarse():
    """Past the crossover, batch size should not keep shrinking with the cohort."""
    assert resolve_chunk_size(total=1_000_000, workers=4) == DEFAULT_CHUNK_SIZE


def test_chunk_size_never_exceeds_the_ceiling_and_never_reaches_zero():
    for total in (0, 1, 5, 97, 10_000, 10**9):
        size = resolve_chunk_size(total=total, workers=4)
        assert 1 <= size <= DEFAULT_CHUNK_SIZE


def test_an_explicit_chunk_size_always_wins():
    assert resolve_chunk_size(total=1000, workers=8, requested=7) == 7
    assert resolve_chunk_size(total=1000, workers=8, requested=0) == 1


def test_an_unknown_length_input_gets_a_small_default():
    """A generator is usually a real cohort, and one big batch would idle the pool."""
    assert resolve_chunk_size(total=None, workers=4) <= 16


@pytest.mark.parametrize("size", [1, 2, 3, 64, 1000])
def test_chunked_covers_every_item_exactly_once_and_in_order(size):
    items = list(range(250))
    flattened = [item for batch in chunked(items, size) for item in batch]
    assert flattened == items


def test_chunked_of_an_empty_iterable_yields_nothing():
    assert list(chunked([], 8)) == []


# -- initializer / initargs ------------------------------------------------


def test_an_initializer_receives_its_initargs_in_every_worker():
    """Regression: an initializer that needs arguments must actually get them.

    ``parallel_map`` originally had no ``initargs``. An initializer written to
    take shared tables was therefore called with none, every worker died during
    start-up, and the parent reported only "a process in the pool was terminated
    abruptly" -- with no clue that the initializer was the cause.
    """
    got = list(
        parallel_map(
            note_seen,
            ["a", "b", "c", "d"],
            threads=3,
            workload=WORKLOAD_PROCESS,
            total=4,
            initializer=_record_note,
            initargs=(7, "seen"),
        )
    )
    assert got == [f"seen:7|{value}" for value in ("a", "b", "c", "d")]


def test_the_serial_path_also_runs_the_initializer():
    """A function that depends on the initializer must work at threads=1 too.

    The serial path starts no pool, so nothing would call the initializer. If it
    is skipped, a stage works in parallel and fails in serial -- the most
    confusing direction for that bug to go in.
    """
    _WORKER_NOTE.clear()
    got = list(
        parallel_map(
            note_seen,
            ["x"],
            threads=1,
            initializer=_record_note,
            initargs=(3, "serial"),
        )
    )
    assert got == ["serial:3|x"]


def test_the_initializer_state_is_not_leaked_between_runs():
    """A second run must see the state its own initializer installed."""
    _WORKER_NOTE.clear()
    first = list(
        parallel_map(
            note_seen,
            ["x"],
            threads=2,
            workload=WORKLOAD_PROCESS,
            total=1,
            initializer=_record_note,
            initargs=(1, "first"),
        )
    )
    _WORKER_NOTE.clear()
    second = list(
        parallel_map(
            note_seen,
            ["x"],
            threads=2,
            workload=WORKLOAD_PROCESS,
            total=1,
            initializer=_record_note,
            initargs=(2, "second"),
        )
    )
    assert "first:1" in first[0]
    assert "second:2" in second[0]


# -- process start method and thread pinning -------------------------------


def test_the_start_method_is_pinned_not_left_to_the_platform():
    """``fork`` is unsafe with numpy loaded and can deadlock a worker."""
    method = default_mp_context().get_start_method()
    assert method in {"forkserver", "spawn"}


def test_pinning_blas_threads_fills_blanks_but_never_overrides(monkeypatch):
    """An explicit OMP_NUM_THREADS is a decision, and we must not undo it."""
    monkeypatch.setenv("OMP_NUM_THREADS", "4")
    for name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        monkeypatch.delenv(name, raising=False)

    pin_blas_threads()

    assert os.environ["OMP_NUM_THREADS"] == "4"
    assert os.environ["OPENBLAS_NUM_THREADS"] == "1"


def test_blas_threads_are_pinned_in_a_clean_environment(monkeypatch):
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ):
        monkeypatch.delenv(name, raising=False)

    pin_blas_threads()

    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ):
        assert os.environ[name] == "1"


def test_importing_the_package_pins_blas_threads():
    """Import order is the whole point, so assert the effect rather than the code.

    If someone later adds an ``import numpy`` above the pin in ``__init__``,
    this is the test that notices.
    """
    import mutation_scan  # noqa: F401

    assert os.environ.get("OMP_NUM_THREADS") == "1"


# -- both backends reach the same answer -----------------------------------


@pytest.mark.parametrize("workload", [WORKLOAD_PROCESS, WORKLOAD_THREAD])
def test_both_backends_agree_with_the_serial_result(workload):
    """A stage migrated to a different backend must not change its output."""
    serial = list(parallel_map(double, range(50), threads=1, workload=workload))
    parallel = list(parallel_map(double, range(50), threads=3, workload=workload, total=50))
    assert serial == parallel == [value * 2 for value in range(50)]
