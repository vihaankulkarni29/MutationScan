"""Unified output schemas for the AMR benchmarking harness.

Each external tool speaks its own dialect: AMRFinderPlus emits a
tabular protein-variant report, Mykrobe a JSON genotype call,
Kleborate a tabular summary, StarAMR a ResFinder-style TSV. The
schemas here are the *canonical* form every tool is normalised into
before any metric is computed, so a sensitivity number means the
same thing regardless of which tool produced it.

The canonical mutation record mirrors what MutationScan itself
writes to ``mutations.csv`` (Accession, Gene, Mutation, Ref_AA,
Position, Alt_AA), so MutationScan is its own first test case: the
harness can score a tool against MutationScan's output and against
the ground truth with the same code path.

Nothing here imports the tool wrappers or the harness; the
dependency runs the other way (tools and harness import these),
which keeps the schemas importable from a bare test without
pulling in subprocess plumbing.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any


class MutationType(str, Enum):
    """The kind of change a record describes.

    ``SNV`` is a single-residue substitution. ``LoF`` collapses the
    loss-of-function events MutationScan reports (frameshift,
    nonsense/truncation knockout) into one class, matching the
    pipeline's own ``LOF_ALLELE`` convention. ``indel`` covers
    in-frame insertions/deletions that are not knockouts.
    """

    SNV = "snv"
    LoF = "lof"
    INDEL = "indel"
    OTHER = "other"


@dataclass(frozen=True)
class MutationRecord:
    """One called mutation in one genome, canonicalised.

    ``gene`` and ``mutation`` are the join keys: MutationScan's own
    tables use ``Gene`` + ``Mutation`` (e.g. ``gyrA`` + ``S83L``),
    and the ground-truth and tool records are normalised to the same
    pair so a match is a string comparison, not a fuzzy one.
    """

    accession: str
    gene: str
    mutation: str
    mutation_type: MutationType = MutationType.OTHER
    ref_aa: str = ""
    position: int | None = None
    alt_aa: str = ""

    @property
    def key(self) -> tuple[str, str, str]:
        """The identity of a mutation, independent of its annotation."""
        return (self.accession, self.gene, self.mutation)


@dataclass(frozen=True)
class EpistasisRecord:
    """One statistically tested mutation pair.

    ``variant_a``/``variant_b`` are ``Gene:Mutation`` strings, the
    form MutationScan writes to ``epistasis_networks.csv``. ``q_value``
    and ``sign`` carry the significance and direction; a benchmark
    compares these against a ground-truth interaction set.
    """

    variant_a: str
    variant_b: str
    gene_a: str = ""
    gene_b: str = ""
    n_both: int = 0
    q_value: float = 1.0
    sign: str = "neutral"

    @property
    def key(self) -> frozenset[str]:
        """Unordered pair identity, so (A,B) and (B,A) are the same pair."""
        return frozenset((self.variant_a, self.variant_b))

    @property
    def gene_key(self) -> frozenset[str]:
        """Unordered *gene* pair, aggregating variant-level calls.

        Ground truth is usually known at the gene level (which
        genes co-mutate), while a tool reports specific variant
        pairs. Matching on the gene pair is what makes the two
        comparable.
        """
        return frozenset((self.gene_a, self.gene_b))


@dataclass(frozen=True)
class DnDsRecord:
    """Per-gene selection pressure, as MutationScan writes ``dn_ds.csv``."""

    gene: str
    n_sites: int = 0
    s_sites: int = 0
    dn: float = 0.0
    ds: float = 0.0
    dn_ds: float | None = None
    n_genomes: int = 0


@dataclass(frozen=True)
class LineageRecord:
    """A genome's population-structure label (ST, phylogroup, ...)."""

    accession: str
    lineage: str


@dataclass
class ToolOutput:
    """Everything one tool run produced, in canonical form.

    ``metadata`` carries tool name, version and wall time so the
    report can attribute a number to a specific build. Empty lists
    are legitimate: a tool that found nothing is a result, not an
    error, and the metrics must be able to score it.
    """

    tool: str
    version: str = ""
    runtime_seconds: float = 0.0
    peak_memory_mb: float = 0.0
    mutations: list[MutationRecord] = field(default_factory=list)
    epistasis: list[EpistasisRecord] = field(default_factory=list)
    dn_ds: list[DnDsRecord] = field(default_factory=list)
    lineages: list[LineageRecord] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def mutation_keys(self) -> set[tuple[str, str, str]]:
        """Set of ``(accession, gene, mutation)`` for fast set metrics."""
        return {record.key for record in self.mutations}

    @property
    def epistasis_keys(self) -> set[frozenset[str]]:
        """Set of unordered variant-pair keys."""
        return {record.key for record in self.epistasis}

    @property
    def epistasis_gene_keys(self) -> set[frozenset[str]]:
        """Set of unordered gene-pair keys, aggregating variants."""
        return {record.gene_key for record in self.epistasis}


__all__ = [
    "DnDsRecord",
    "EpistasisRecord",
    "LineageRecord",
    "MutationRecord",
    "MutationType",
    "ToolOutput",
]
