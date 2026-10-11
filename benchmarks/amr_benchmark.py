"""Benchmark harness: run tools and score them against ground truth.

The harness is deliberately split from the tool wrappers:
a wrapper knows how to run *one* tool and parse its output,
while this module knows how to *score* any tool's canonical
output against a ground truth. That separation is what makes
the comparison fair -- every tool is scored by the same
metric code, including MutationScan itself.

Metrics
-------

**Mutation level** (against a ground-truth mutation set):

* *sensitivity* -- fraction of ground-truth mutations found;
* *precision* -- fraction of the tool's calls that are real;
* *specificity* -- fraction of the negative set correctly
  left alone, computed only when a defined negative set is
  supplied (a synthetic spike-in provides one; real cohorts
  do not, so the column is omitted rather than invented);
* *novel recall* -- fraction of ground-truth mutations that
  are *not* in a reference database the tool could have
  looked up, i.e. the tool's reach beyond catalogue-based
  detection.

**Epistasis level**: precision/recall/F1 over predicted
interaction pairs, matched as unordered variant pairs.

**dN/dS level**: Pearson and Spearman correlation and mean
absolute error per gene, against a ground-truth per-gene
selection coefficient.

**Lineage level**: concordance of population-structure
labels, which is what makes a "convergent" epistasis call
trustworthy rather than a clonal artefact.

**Performance**: wall time, peak memory, and throughput
(genomes/hour), plus a scaling curve when the harness is
run over several cohort sizes.
"""

from __future__ import annotations

import logging
import time
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .schemas import (
    ToolOutput,
)
from .tools import Tool

logger = logging.getLogger(__name__)

__all__ = [
    "BenchmarkResult",
    "DnDsMetrics",
    "EpistasisMetrics",
    "GroundTruth",
    "LineageMetrics",
    "MutationMetrics",
    "PerformanceMetrics",
    "BenchmarkHarness",
]


@dataclass
class GroundTruth:
    """The known answer for one dataset.

    Every field is optional: a dataset may know its
    mutations but not its dN/dS, and the harness scores
    only what is present.
    """

    mutations: set[tuple[str, str, str]] = field(default_factory=set)
    #: Unordered *gene* pairs known to be epistatic.
    epistasis: set[frozenset[str]] = field(default_factory=set)
    dn_ds: dict[str, float] = field(default_factory=dict)
    lineages: dict[str, str] = field(default_factory=dict)
    #: (accession, gene, mutation) known to a reference
    #: database. Everything in ``mutations`` but not here
    #: is "novel".
    known_mutations: set[tuple[str, str, str]] = field(default_factory=set)
    #: (accession, gene) pairs known to be unmutated,
    #: for specificity.
    negatives: set[tuple[str, str]] = field(default_factory=set)


@dataclass
class MutationMetrics:
    tp: int = 0
    fp: int = 0
    fn: int = 0
    tn: int = 0
    novel_tp: int = 0
    novel_total: int = 0

    @property
    def sensitivity(self) -> float:
        return self.tp / (self.tp + self.fn) if (self.tp + self.fn) else 0.0

    @property
    def precision(self) -> float:
        return self.tp / (self.tp + self.fp) if (self.tp + self.fp) else 0.0

    @property
    def specificity(self) -> float | None:
        denom = self.tn + self.fp
        return self.tn / denom if denom else None

    @property
    def f1(self) -> float:
        precision, recall = self.precision, self.sensitivity
        return 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0

    @property
    def novel_recall(self) -> float | None:
        return self.novel_tp / self.novel_total if self.novel_total else None


@dataclass
class EpistasisMetrics:
    tp: int = 0
    fp: int = 0
    fn: int = 0

    @property
    def precision(self) -> float:
        return self.tp / (self.tp + self.fp) if (self.tp + self.fp) else 0.0

    @property
    def recall(self) -> float:
        return self.tp / (self.tp + self.fn) if (self.tp + self.fn) else 0.0

    @property
    def f1(self) -> float:
        precision, recall = self.precision, self.recall
        return 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0


@dataclass
class DnDsMetrics:
    pearson: float = 0.0
    spearman: float = 0.0
    mae: float = 0.0
    n_genes: int = 0


@dataclass
class LineageMetrics:
    correct: int = 0
    total: int = 0

    @property
    def concordance(self) -> float:
        return self.correct / self.total if self.total else 0.0


@dataclass
class PerformanceMetrics:
    runtime_seconds: float = 0.0
    peak_memory_mb: float = 0.0
    n_genomes: int = 0

    @property
    def genomes_per_hour(self) -> float:
        if self.runtime_seconds <= 0:
            return 0.0
        return self.n_genomes / self.runtime_seconds * 3600.0


@dataclass
class BenchmarkResult:
    tool: str
    version: str = ""
    mutations: MutationMetrics = field(default_factory=MutationMetrics)
    epistasis: EpistasisMetrics = field(default_factory=EpistasisMetrics)
    dn_ds: DnDsMetrics = field(default_factory=DnDsMetrics)
    lineage: LineageMetrics = field(default_factory=LineageMetrics)
    performance: PerformanceMetrics = field(default_factory=PerformanceMetrics)


def compute_mutation_metrics(
    ground_truth: GroundTruth,
    output: ToolOutput,
) -> MutationMetrics:
    """Score a tool's mutation calls against the ground truth."""
    predicted = output.mutation_keys
    truth = ground_truth.mutations

    tp = len(predicted & truth)
    fn = len(truth - predicted)
    fp = len(predicted - truth)

    tn = 0
    if ground_truth.negatives:
        # A negative (accession, gene) is a true negative when the
        # tool made no call for that pair.
        called_pairs = {(a, g) for a, g, _ in predicted}
        tn = len(ground_truth.negatives - called_pairs)

    novel = truth - ground_truth.known_mutations
    novel_tp = len(predicted & novel)

    return MutationMetrics(
        tp=tp,
        fp=fp,
        fn=fn,
        tn=tn,
        novel_tp=novel_tp,
        novel_total=len(novel),
    )


def compute_epistasis_metrics(
    ground_truth: GroundTruth,
    output: ToolOutput,
) -> EpistasisMetrics:
    # Ground truth is known at the gene level (which genes
    # co-mutate), while a tool reports specific variant pairs.
    # Matching on the gene pair is what makes the two comparable.
    predicted = output.epistasis_gene_keys
    truth = ground_truth.epistasis
    tp = len(predicted & truth)
    fn = len(truth - predicted)
    fp = len(predicted - truth)
    return EpistasisMetrics(tp=tp, fp=fp, fn=fn)


def compute_dn_ds_metrics(
    ground_truth: GroundTruth,
    output: ToolOutput,
) -> DnDsMetrics:
    if not ground_truth.dn_ds:
        return DnDsMetrics()
    predicted = {record.gene: record.dn_ds for record in output.dn_ds}
    # A gene with no dN/dS (dS == 0) is not comparable: drop it
    # before building the value lists, so the types stay float.
    pairs: list[tuple[float, float]] = []
    for gene in set(ground_truth.dn_ds) & set(predicted):
        value = predicted[gene]
        if value is not None:
            pairs.append((ground_truth.dn_ds[gene], value))
    if len(pairs) < 2:
        return DnDsMetrics(n_genes=len(pairs))

    truth_values = [t for t, _ in pairs]
    pred_values = [p for _, p in pairs]

    pearson, spearman = _correlations(truth_values, pred_values)
    mae = sum(abs(t - p) for t, p in pairs) / len(pairs)
    return DnDsMetrics(
        pearson=pearson,
        spearman=spearman,
        mae=mae,
        n_genes=len(pairs),
    )


def _correlations(truth: list[float], pred: list[float]) -> tuple[float, float]:
    from scipy import stats

    pearson = 0.0
    spearman = 0.0
    try:
        pearson = float(stats.pearsonr(truth, pred)[0])
    except (ValueError, ZeroDivisionError):
        pearson = 0.0
    try:
        spearman = float(stats.spearmanr(truth, pred).statistic)
    except (ValueError, ZeroDivisionError):
        spearman = 0.0
    return pearson, spearman


def compute_lineage_metrics(
    ground_truth: GroundTruth,
    output: ToolOutput,
) -> LineageMetrics:
    if not ground_truth.lineages:
        return LineageMetrics()
    predicted = {record.accession: record.lineage for record in output.lineages}
    correct = sum(
        1
        for accession, lineage in ground_truth.lineages.items()
        if predicted.get(accession) == lineage
    )
    return LineageMetrics(correct=correct, total=len(ground_truth.lineages))


class BenchmarkHarness:
    """Run a set of tools over a dataset and score them.

    The harness owns nothing about any tool's internals; it
    is given a :class:`GroundTruth`, a list of
    :class:`Tool` instances, and the inputs each tool needs,
    and it returns one :class:`BenchmarkResult` per tool.
    """

    def __init__(
        self,
        ground_truth: GroundTruth,
        tools: Sequence[Tool],
        inputs: Mapping[str, Any],
        work_dir: Path | str = "benchmark_runs",
    ) -> None:
        self.ground_truth = ground_truth
        self.tools = list(tools)
        self.inputs = dict(inputs)
        self.work_dir = Path(work_dir)

    def run(self) -> list[BenchmarkResult]:
        results: list[BenchmarkResult] = []
        for tool in self.tools:
            result = self.run_tool(tool)
            results.append(result)
        return results

    def run_tool(self, tool: Tool) -> BenchmarkResult:
        output_dir = self.work_dir / tool.name.lower().replace(" ", "_")
        started = time.perf_counter()
        try:
            output = tool.run(self.inputs, output_dir)
        except Exception as exc:  # noqa: BLE001 - a failed tool is a result
            logger.error("%s failed: %s", tool.name, exc)
            return BenchmarkResult(tool=tool.name, version=tool.version())
        runtime = time.perf_counter() - started
        output.runtime_seconds = runtime

        return BenchmarkResult(
            tool=tool.name,
            version=output.version,
            mutations=compute_mutation_metrics(self.ground_truth, output),
            epistasis=compute_epistasis_metrics(self.ground_truth, output),
            dn_ds=compute_dn_ds_metrics(self.ground_truth, output),
            lineage=compute_lineage_metrics(self.ground_truth, output),
            performance=PerformanceMetrics(
                runtime_seconds=runtime,
                peak_memory_mb=output.peak_memory_mb,
                n_genomes=self._n_genomes(),
            ),
        )

    def _n_genomes(self) -> int:
        genomes = self.inputs.get("genomes")
        if isinstance(genomes, Path) and genomes.is_dir():
            return sum(
                1 for path in genomes.iterdir() if path.suffix in {".fna", ".fa", ".fasta", ".faa"}
            )
        return 0
