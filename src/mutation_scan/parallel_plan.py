"""The one place that answers "how many workers will each stage actually use?".

There used to be two wirings of the same numbers: the CLI read ``--threads`` and
``--variant-threads`` from flags, while the Snakefile read ``CONFIG.threads`` and
``CONFIG.variant_threads``. Nothing checked that they agreed, and one of them
already disagreed about epistasis -- ``--threads 8`` widened extraction and
variant calling while epistasis silently stayed on its own default of 1, with no
way to change it from the command line.

So: one plan, derived from the :class:`~mutation_scan.config.Config`, used by the
run summary, by ``config-check``, and available to the workflow. Every place that
reports on a run's parallelism reads this rather than re-deriving it.

The distinction that matters is *requested* versus *effective*. A user who asks
for 8 workers and silently gets 1 has been lied to, and nothing else in the run's
output would reveal it. That is why the plan is computed here rather than left to
each stage.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .config import Config
from .parallel import WORKLOAD_PROCESS, WORKLOAD_THREAD, available_cores, resolve_workers

__all__ = ["Stage", "describe", "format_plan"]


@dataclass(frozen=True)
class Stage:
    """One stage's parallelism: what was asked for, and what will run."""

    name: str
    workload: str
    requested: int
    effective: int
    #: Set when the stage will fall back to in-process work below a size
    #: threshold, which is the case a reader is most likely to misjudge.
    serial_below: int | None = None
    #: Set when the stage is deliberately not parallelised, so that "1 worker"
    #: is not misread as a setting the user chose.
    not_parallelised: str = ""

    @property
    def kind(self) -> str:
        return "process" if self.workload == WORKLOAD_PROCESS else "thread"

    @property
    def is_serial(self) -> bool:
        return self.effective <= 1

    @property
    def clamped(self) -> bool:
        """The request was reduced: the stage will run narrower than asked."""
        return self.effective < self.requested

    @property
    def oversubscribed(self) -> bool:
        """Honoured as asked, but wider than the machine has cores.

        Distinct from :attr:`clamped` and worth separating: an oversubscribed
        stage runs at full width and finishes *slower*, which reads very
        differently in a timing report from one that quietly ran at a tenth of
        the requested width.
        """
        return self.effective > available_cores()

    def to_dict(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "stage": self.name,
            "workload": self.kind,
            "requested_workers": self.requested,
            "effective_workers": self.effective,
            "clamped": self.clamped,
            "oversubscribed": self.oversubscribed,
        }
        if self.serial_below is not None:
            payload["serial_below_items"] = self.serial_below
        if self.not_parallelised:
            payload["not_parallelised"] = self.not_parallelised
        return payload


#: Mirrors epistasis._SERIAL_PAIR_LIMIT. Duplicated deliberately: importing it
#: from epistasis would pull scipy into every `config-check`, which must stay
#: fast and dependency-light. test_parallel_plan_agrees_with_the_epistasis_limit
#: fails the moment the two disagree.
_EPIGASTIS_SERIAL_PAIR_LIMIT = 2048


def describe(config: Config) -> list[Stage]:
    """The resolved parallelism for every stage of a run, in run order."""
    return [
        Stage(
            "extraction",
            WORKLOAD_THREAD,
            config.threads,
            resolve_workers(config.threads, WORKLOAD_THREAD),
        ),
        Stage(
            "variant calling",
            WORKLOAD_PROCESS,
            config.variant_threads,
            resolve_workers(config.variant_threads, WORKLOAD_PROCESS),
        ),
        # Co-occurrence is O(genes^2) over 15-odd genes: microseconds, and
        # parallelising it would cost more in coordination than it saves. It is
        # listed so the plan accounts for every stage rather than only the
        # parallel ones, and marked so its single worker is not misread as a
        # setting someone chose.
        Stage(
            "co-occurrence",
            WORKLOAD_THREAD,
            1,
            1,
            not_parallelised="O(genes^2) over a small panel; not worth a pool",
        ),
        Stage(
            "epistasis",
            WORKLOAD_PROCESS,
            config.epistasis.threads,
            resolve_workers(config.epistasis.threads, WORKLOAD_PROCESS),
            serial_below=_EPIGASTIS_SERIAL_PAIR_LIMIT,
        ),
    ]


def format_plan(stages: list[Stage]) -> str:
    """Render the plan for a terminal, flagging anything unexpectedly serial."""
    width = max(len(stage.name) for stage in stages)
    lines = []
    for stage in stages:
        if stage.not_parallelised:
            detail = f"not parallelised ({stage.not_parallelised})"
        else:
            detail = f"{stage.kind}, {stage.requested} requested, {stage.effective} effective"
            if stage.serial_below is not None:
                detail += f"; in-process below {stage.serial_below} pairs"
            if stage.clamped:
                detail += "  <-- clamped"
            elif stage.oversubscribed:
                detail += f"  <-- oversubscribed vs {available_cores()} cores"
        lines.append(f"  {stage.name.ljust(width)}  {detail}")
    return "\n".join(lines)
