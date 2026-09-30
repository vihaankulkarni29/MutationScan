"""One parallel-execution policy for every MutationScan stage.

Every stage that can run more than one thing at once goes through
:func:`parallel_map`. Centralising it means the awkward decisions -- which
backend, how many workers, how work is batched, what happens when a worker
dies -- are made once, in one place, and are the same for extraction, variant
calling and epistasis.

Why the backend is chosen per stage
-----------------------------------
Two kinds of work appear in this pipeline, and they want opposite things:

``WORKLOAD_THREAD``
    The work spends its time blocked on a subprocess (``tblastn``) or inside a C
    extension that releases the GIL. Threads are the right tool: no interpreter
    to respawn, no data to pickle, and the GIL is not the constraint.

``WORKLOAD_PROCESS``
    The work is CPU-bound Python -- Biopython alignment, ``scipy.stats``. The
    GIL serialises it, so threads buy nothing and a process pool is the only way
    to use more than one core.

Choosing wrongly is not a small mistake: running CPU-bound Python on threads
gives a speedup near 1.0 and still pays the coordination cost.

The correctness rules this module enforces
------------------------------------------
**Determinism.** Parallel results are yielded in input order, not completion
order, so ``threads`` can never change the output of a run. That is what lets a
tool like this promise a result is reproducible.

**Bounded memory.** *items* is consumed lazily and at most
``_MAX_PENDING_PER_WORKER`` batches per worker are in flight, so peak memory
does not grow with cohort size. Materialising every task up front is fine for
600 genomes and a slow memory leak for 60,000.

**Loud failure.** A worker that dies raises :class:`ParallelTaskError` naming
the failing batch. Silently converting a crashed worker into a per-item "error"
result makes a broken run look like a data problem, which is how a crash gets
misdiagnosed for a long time.

Notes for callers
-----------------
On the ``spawn`` and ``forkserver`` start methods (Windows, and the POSIX
default here) an imported module is re-executed in each worker. Any script that
drives MutationScan from its own ``__main__`` must therefore guard its entry
point::

    if __name__ == "__main__":
        main()

Without that guard each worker re-runs the script and recursively spawns more
workers until the machine falls over.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Callable, Iterable, Iterator
from concurrent.futures import (
    FIRST_COMPLETED,
    Executor,
    Future,
    ProcessPoolExecutor,
    ThreadPoolExecutor,
    wait,
)
from multiprocessing import get_all_start_methods, get_context
from multiprocessing.context import BaseContext
from pathlib import Path
from typing import Final, TypeVar

logger = logging.getLogger(__name__)

__all__ = [
    "DEFAULT_CHUNK_SIZE",
    "WORKLOAD_PROCESS",
    "WORKLOAD_THREAD",
    "ParallelTaskError",
    "available_cores",
    "chunked",
    "default_mp_context",
    "parallel_map",
    "pin_blas_threads",
    "resolve_chunk_size",
    "resolve_workers",
]

T = TypeVar("T")
R = TypeVar("R")

#: CPU-bound Python: needs separate processes to use more than one core.
WORKLOAD_PROCESS: Final = "process"
#: Subprocess-bound or GIL-releasing work: threads are cheaper and sufficient.
WORKLOAD_THREAD: Final = "thread"

#: Items per submitted batch, as an upper bound. Adaptive sizing will happily go
#: below this; it is the ceiling, not the value.
DEFAULT_CHUNK_SIZE: Final = 64

#: Chunk size used when the caller cannot say how many items there are. Chosen
#: small deliberately: an unknown-length input is usually a generator over a real
#: cohort, and a large chunk means the pool sits idle with one batch in flight
#: before anything is emitted.
_DEFAULT_CHUNK_UNKNOWN_TOTAL: Final = 8

#: Batches in flight per worker. Bounds the memory held by queued work and by
#: results not yet emitted.
_MAX_PENDING_PER_WORKER: Final = 4

#: When the caller gives no total, report progress this often instead.
_PROGRESS_INTERVAL: Final = 500

#: Threads are cheap to create, but each can hold a subprocess open, so past
#: this the file-descriptor budget rather than the CPU is the limiting resource.
_THREAD_CEILING: Final = 64

#: Every threading layer that could otherwise fan out inside our own workers.
_THREAD_ENV_VARS: Final = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)


class ParallelTaskError(RuntimeError):
    """A worker raised, or died, while processing part of the work.

    Attributes:
        completed: how many items finished successfully before the failure.
        pending: how many items were in the failing batch.
        cause: the original exception, or ``None`` if the worker died outright
            without raising (an out-of-memory kill, for instance).
    """

    def __init__(self, message: str, *, completed: int, pending: int, cause: BaseException | None):
        super().__init__(message)
        self.completed = completed
        self.pending = pending
        self.cause = cause


# -- core detection ---------------------------------------------------------


def _cgroup_cpu_limit() -> float | None:
    """CPU quota imposed by the cgroup, if any.

    Containers are the case that matters: ``os.cpu_count()`` reports the *host*
    core count inside a container, so a tool that sizes a process pool from it
    will happily oversubscribe a 2-core container tenfold. This package ships a
    docker-compose file, so that is a supported deployment, not a hypothetical.
    """
    # cgroup v2: "<quota> <period>", with quota "max" when unlimited.
    try:
        v2 = Path("/sys/fs/cgroup/cpu.max")
        if v2.is_file():
            quota_text, _, period_text = v2.read_text().strip().partition(" ")
            if quota_text and quota_text != "max" and period_text:
                return int(quota_text) / int(period_text)
    except (OSError, ValueError, ZeroDivisionError):
        pass

    # cgroup v1: separate quota and period files.
    try:
        quota_file = Path("/sys/fs/cgroup/cpu/cpu.cfs_quota_us")
        period_file = Path("/sys/fs/cgroup/cpu/cpu.cfs_period_us")
        if quota_file.is_file() and period_file.is_file():
            quota_us = int(quota_file.read_text().strip())
            period_us = int(period_file.read_text().strip())
            if quota_us > 0 and period_us > 0:
                return quota_us / period_us
    except (OSError, ValueError, ZeroDivisionError):
        pass

    return None


def available_cores() -> int:
    """Usable CPU count for this process, honouring affinity and cgroup limits.

    Order of preference: the cgroup CPU quota (if containerised and limited),
    then the scheduler affinity mask (which reflects ``taskset``/``cpuset``),
    then the raw core count. Never returns less than 1.
    """
    quota = _cgroup_cpu_limit()
    if quota is not None:
        return max(1, int(quota))

    getaffinity = getattr(os, "sched_getaffinity", None)
    if getaffinity is not None:
        try:
            return max(1, len(getaffinity(0)))
        except OSError:  # pragma: no cover - platform dependent
            pass

    return max(1, os.cpu_count() or 1)


# -- process start method ---------------------------------------------------


def _supported_methods() -> frozenset[str]:
    """Start methods this platform can actually use.

    Probed rather than hardcoded: ``fork``-family methods are absent on Windows
    and restricted in some embedded interpreters, and asking for an unsupported
    method only fails at pool-construction time -- i.e. mid-run.
    """
    try:
        return frozenset(get_all_start_methods())
    except (AttributeError, ValueError):  # pragma: no cover - very old runtimes
        return frozenset({"spawn"})


def default_mp_context() -> BaseContext:
    """The multiprocessing context used for process-backed work.

    Deliberately not the platform default. ``fork`` is Python's default on POSIX
    and is unsafe here: the parent has already imported numpy/pandas and may hold
    threads, and forking a multi-threaded process can deadlock the child inside
    an allocator or a lock held by a thread that does not exist post-fork.
    ``forkserver`` forks from a clean single-threaded server, keeping most of
    ``fork``'s speed without the hazard. Windows has no such choice and gets
    ``spawn``.

    Pinned rather than defaulted, so a run cannot change behaviour because a
    different Python version changed its mind about the default.
    """
    available = _supported_methods()
    return get_context("forkserver" if "forkserver" in available else "spawn")


# -- worker sizing ----------------------------------------------------------


def resolve_workers(threads: int, workload: str = WORKLOAD_PROCESS) -> int:
    """Turn a user-supplied thread count into a safe worker count.

    ``threads <= 0`` is treated as "not specified" and clamped to 1 rather than
    rejected: a bad number should not cost the user their run. An explicit
    request above the available core count is honoured -- a user on a loaded
    machine may knowingly oversubscribe for throughput -- but is reported, since
    oversubscription makes a stage slower rather than faster.
    """
    requested = max(1, int(threads))
    cores = available_cores()

    if requested > cores:
        logger.warning(
            "Requested %d worker(s) but only %d core(s) are available; "
            "oversubscribing will slow this stage down, not speed it up.",
            requested,
            cores,
        )
    if workload == WORKLOAD_THREAD and requested > _THREAD_CEILING:
        logger.warning(
            "Capping %d threads at %d to stay within file-descriptor limits.",
            requested,
            _THREAD_CEILING,
        )
        return _THREAD_CEILING
    return requested


# -- batching ---------------------------------------------------------------


def chunked(items: Iterable[T], size: int = DEFAULT_CHUNK_SIZE) -> Iterator[list[T]]:
    """Yield lists of at most *size* items, consuming *items* lazily.

    A generator stays a generator: nothing is buffered beyond one batch, which
    is what keeps a cohort of any size at constant memory.
    """
    size = max(1, int(size))
    batch: list[T] = []
    for item in items:
        batch.append(item)
        if len(batch) >= size:
            yield batch
            batch = []
    if batch:
        yield batch


def _map_batch(fn: Callable[[T], R], batch: list[T]) -> list[R]:
    """Apply *fn* to every item in *batch*.

    Runs inside a worker, so *fn* must be a module-level function that pickles by
    reference -- not a lambda, closure, or locally-defined function.
    """
    return [fn(item) for item in batch]


def resolve_chunk_size(
    total: int | None,
    workers: int,
    requested: int | None = None,
) -> int:
    """Choose a batch size that actually spreads work across the workers.

    A fixed chunk size is the wrong tool. Too large and a small cohort collapses
    into a single batch that runs on one worker while the rest sit idle -- the
    stage is then no faster than running serially, with none of the
    observability to suggest why. Too small and per-batch pickling and scheduling
    dominate the actual work.

    So: aim for :data:`_MAX_PENDING_PER_WORKER` batches per worker, capped at
    :data:`DEFAULT_CHUNK_SIZE` so a big cohort keeps batches coarse enough for
    pickling to be negligible. With no known length, fall back to a small batch:
    an unknown-length input is almost always a generator over a real cohort, and
    we would rather start parallel early than discover one batch too late.
    """
    if requested is not None:
        return max(1, int(requested))
    if total is None:
        return _DEFAULT_CHUNK_UNKNOWN_TOTAL

    workers = max(1, int(workers))
    target_batches = max(1, workers * _MAX_PENDING_PER_WORKER)
    # Ceiling division, so we never end up with fewer batches than we aimed for.
    needed = -(-max(1, int(total)) // target_batches)
    return max(1, min(DEFAULT_CHUNK_SIZE, needed))


def _make_executor(
    workers: int,
    workload: str,
    initializer: Callable[..., None] | None,
    initargs: tuple[object, ...],
) -> Executor:
    if workload == WORKLOAD_THREAD:
        return ThreadPoolExecutor(max_workers=workers, thread_name_prefix="mscan")
    # mypy infers the executor's initargs type from the initializer, and a
    # variadic `Callable[..., None]` cannot pin it down, so it defaults to the
    # empty tuple and rejects any real arguments. The runtime behaviour that
    # this cannot express is covered by the initargs tests instead.
    return ProcessPoolExecutor(
        max_workers=workers,
        mp_context=default_mp_context(),
        initializer=initializer,
        initargs=initargs,  # type: ignore[arg-type]
    )


def _log_progress(*, completed: int, total: int | None, label: str) -> None:
    if not label or completed == 0:
        return
    if total:
        if completed % max(1, total // 20) and completed != total:
            return
        logger.info("  %s %d/%d (%.0f%%)", label, completed, total, completed / total * 100)
    elif completed % _PROGRESS_INTERVAL == 0:
        logger.info("  %s %d item(s) done", label, completed)


# -- the public map ---------------------------------------------------------


def parallel_map(
    fn: Callable[[T], R],
    items: Iterable[T],
    *,
    threads: int = 1,
    workload: str = WORKLOAD_PROCESS,
    total: int | None = None,
    chunk_size: int | None = None,
    initializer: Callable[..., None] | None = None,
    initargs: tuple[object, ...] = (),
    label: str = "",
) -> Iterator[R]:
    """Map *fn* over *items*, preserving input order.

    Args:
        fn: a module-level callable. For :data:`WORKLOAD_PROCESS` it must be
            importable by name in the worker, which rules out lambdas, closures
            and locally-defined functions.
        items: any iterable, including a generator. Consumed lazily.
        threads: worker count. ``1`` runs in the calling process, with no pool and
            no pickling.
        workload: :data:`WORKLOAD_PROCESS` or :data:`WORKLOAD_THREAD`.
        total: how many items there are. Pass it whenever you know: it is what
            lets :func:`resolve_chunk_size` spread a short input across the
            workers, and it is what makes the progress log show a percentage.
        chunk_size: items per batch. ``None`` picks a size from *total* and
            *threads*; pass a number only if you have a reason to override that.
        initializer: run once per worker process before any task -- the place to
            build something expensive, or to install shared state, that would
            otherwise be rebuilt per task.
        initargs: arguments passed to *initializer*. Sent to each worker once, so
            this is far cheaper than putting the same data in every task.
        label: short name for the progress log, e.g. ``"variant calling"``.

    Yields:
        One result per item, **in input order**, whatever ``threads`` is set to.

    Raises:
        ParallelTaskError: a worker raised or died. The serial path re-raises the
            original error this way too, so callers get one exception type
            regardless of worker count.
    """
    workers = resolve_workers(threads, workload)

    if workers == 1:
        # The serial path runs no pool, so nothing would call the initializer.
        # Run it here so shared state is set up identically either way, and so a
        # function that needs the initializer works at threads=1 too.
        if initializer is not None:
            initializer(*initargs)
        yield from _run_serial(fn, items, total, label)
        return

    yield from _run_pooled(
        fn,
        items,
        workers=workers,
        workload=workload,
        total=total,
        chunk_size=resolve_chunk_size(total, workers, chunk_size),
        initializer=initializer,
        initargs=initargs,
        label=label,
    )


def _run_serial(
    fn: Callable[[T], R], items: Iterable[T], total: int | None, label: str
) -> Iterator[R]:
    completed = 0
    for index, item in enumerate(items, start=1):
        try:
            result = fn(item)
        except Exception as exc:  # noqa: BLE001 - re-raised with context
            raise ParallelTaskError(
                f"{label or 'task'} failed on item {index}: {exc}",
                completed=completed,
                pending=1,
                cause=exc,
            ) from exc
        completed += 1
        yield result
        _log_progress(completed=completed, total=total, label=label)


def _run_pooled(
    fn: Callable[[T], R],
    items: Iterable[T],
    *,
    workers: int,
    workload: str,
    total: int | None,
    chunk_size: int,
    initializer: Callable[..., None] | None,
    initargs: tuple[object, ...],
    label: str,
) -> Iterator[R]:
    """Drive a pool with a bounded window of in-flight batches.

    The window is what bounds memory: at most ``workers * _MAX_PENDING_PER_WORKER``
    batches exist at once, so a generator input is never fully materialised.

    Results are buffered and emitted in input order even though batches complete
    out of order. That is the property that makes the output independent of
    ``threads``, and it is worth the buffering: without it, a run's output would
    depend on scheduler timing and results would not be reproducible.
    """
    max_pending = max(1, workers * _MAX_PENDING_PER_WORKER)
    batches = chunked(items, chunk_size)

    pool = _make_executor(workers, workload, initializer, initargs)

    # future -> (batch index, batch size). The size is recorded up front because
    # a worker that dies cannot tell us how much work it was holding.
    pending: dict[Future[list[R]], tuple[int, int]] = {}
    ready: dict[int, list[R]] = {}  # batch index -> results
    next_batch = 0  # index of the next batch to submit
    next_to_emit = 0  # lowest batch index not yet yielded
    completed = 0
    exhausted = False

    def _fill() -> None:
        nonlocal next_batch, exhausted
        while not exhausted and len(pending) < max_pending:
            try:
                batch = next(batches)
            except StopIteration:
                exhausted = True
                return
            future = pool.submit(_map_batch, fn, batch)
            pending[future] = (next_batch, len(batch))
            next_batch += 1

    try:
        _fill()
        while pending:
            done, _ = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                index, batch_size = pending.pop(future)
                try:
                    outcome = future.result()
                except Exception as exc:  # noqa: BLE001 - converted to a typed error
                    raise ParallelTaskError(
                        f"{label or 'task'} worker failed on a batch of "
                        f"{batch_size} item(s): {exc}",
                        completed=completed,
                        pending=batch_size,
                        cause=exc,
                    ) from exc
                ready[index] = outcome
                completed += len(outcome)

            _fill()

            while next_to_emit in ready:
                yield from ready.pop(next_to_emit)
                next_to_emit += 1
                _log_progress(completed=completed, total=total, label=label)
    finally:
        # Also covers the caller abandoning the iterator early: a generator that
        # is closed raises GeneratorExit here, and we must not leave workers
        # running or queued batches holding references to data.
        for future in pending:
            future.cancel()
        pool.shutdown(wait=True, cancel_futures=True)


# -- BLAS thread pinning ----------------------------------------------------


def pin_blas_threads(threads: int = 1) -> None:
    """Stop numpy/scipy/BLAS from spawning their own thread pool per worker.

    Without this, ``--variant-threads 8`` on a machine whose BLAS defaults to 8
    threads is really 64 threads, and the contention makes the stage *slower*
    than not parallelising it at all.

    Must run before numpy is imported: BLAS reads these variables once, at load
    time. That is why :mod:`mutation_scan` calls this from its ``__init__``
    before touching any submodule. A caller who imported numpy themselves
    first has already missed the window, and the pin will not take effect.

    Values the user set explicitly are respected; this only fills in blanks, so
    it never overrides a deliberate ``OMP_NUM_THREADS`` in the environment.
    """
    for name in _THREAD_ENV_VARS:
        if os.environ.get(name) is None:
            os.environ[name] = str(threads)
