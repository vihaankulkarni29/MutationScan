"""Descriptive co-occurrence of mutations across genomes.

The question this stage answers, and the only one: **in the genomes where two
genes were both evaluated, how often are both of them mutated?**

That is a count, not a claim. There is no p-value, no enrichment, no odds
ratio, no severity weight and no network score anywhere in this module. Two
genes appearing together frequently may reflect shared ancestry, sampling bias
in the cohort, or nothing at all -- deciding which is the reader's work, not
the tool's.

What the tool *does* guarantee is an honest denominator. A gene pair is counted
only over genomes where **both** genes passed variant calling, so "gene B is
not mutated here" is never confused with "gene B was never observed here". A
gene missing from half the cohort will not manufacture co-occurrence with the
half that has it.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path
from typing import Iterable, Mapping, Optional, Sequence

import pandas as pd

from .variants import VARIANT_STATUS_CALLED

logger = logging.getLogger(__name__)

__all__ = [
    "GENE_COLUMNS",
    "VARIANT_COLUMNS",
    "CooccurrenceResult",
    "cooccurrence",
    "eligibility_from_qc",
    "load_result",
    "read_cooccurrence",
]

#: Long-form gene-pair table. Counts first, then the two plain fractions.
GENE_COLUMNS = (
    "Gene_A",
    "Gene_B",
    "N_Both",
    "N_A",
    "N_B",
    "N_Eligible",
    "Frac_B_given_A",
    "Frac_A_given_B",
)

#: Same shape at mutation resolution; pair eligibility still comes from genes.
VARIANT_COLUMNS = (
    "Variant_A",
    "Variant_B",
    "Gene_A",
    "Gene_B",
    "N_Both",
    "N_A",
    "N_B",
    "N_Eligible",
    "Frac_B_given_A",
    "Frac_A_given_B",
)


@dataclass
class CooccurrenceResult:
    """The three tables this stage can produce.

    Attributes:
        genes: one row per unordered gene pair.
        genes_matrix: gene x gene counts; cell is ``N_Both``, diagonal is the
            number of genomes mutating that gene.
        variants: one row per unordered ``Gene:Mutation`` pair, or ``None`` when
            the variant level was not requested.
        gene_totals: genomes mutating each gene (the matrix diagonal).
        eligible_per_gene: genomes where each gene was evaluated.
    """

    genes: pd.DataFrame
    genes_matrix: pd.DataFrame
    variants: Optional[pd.DataFrame] = None
    gene_totals: dict[str, int] = field(default_factory=dict)
    eligible_per_gene: dict[str, int] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.genes)


def eligibility_from_qc(qc: pd.DataFrame) -> dict[str, set[str]]:
    """Genomes in which each gene was successfully evaluated.

    Read from the *variant* QC table rather than the extraction manifest: a gene
    can be extracted and still be rejected by the identity or coverage gate, and
    a rejected gene was not evaluated for mutations, so it must not sit in
    anyone's denominator.

    Returns:
        ``{gene: {accession, ...}}``.
    """
    if qc is None or qc.empty:
        return {}

    missing = [
        column for column in ("Accession", "Gene", "Status") if column not in qc.columns
    ]
    if missing:
        raise ValueError(f"Variant QC is missing column(s): {', '.join(missing)}")

    called = qc[qc["Status"] == VARIANT_STATUS_CALLED]
    eligible: dict[str, set[str]] = {}
    for row in called.itertuples():
        eligible.setdefault(str(row.Gene), set()).add(str(row.Accession))
    return eligible


def _mutated_sets(
    mutations: pd.DataFrame, key: str
) -> dict[str, set[str]]:
    """``{key value: {accession where it is mutated}}``."""
    if mutations.empty:
        return {}
    grouped: dict[str, set[str]] = {}
    for value, group in mutations.groupby(key, sort=True):
        grouped[str(value)] = {str(accession) for accession in group["Accession"]}
    return grouped


def _pair_row(
    label_a: str,
    label_b: str,
    mutated_a: set[str],
    mutated_b: set[str],
    eligible: set[str],
) -> Optional[dict]:
    """Counts for one unordered pair, restricted to the shared denominator."""
    if not eligible:
        return None

    in_a = mutated_a & eligible
    in_b = mutated_b & eligible
    both = in_a & in_b

    return {
        "N_Both": len(both),
        "N_A": len(in_a),
        "N_B": len(in_b),
        "N_Eligible": len(eligible),
        "Frac_B_given_A": round(len(both) / len(in_a), 4) if in_a else 0.0,
        "Frac_A_given_B": round(len(both) / len(in_b), 4) if in_b else 0.0,
    }


def _gene_pairs(
    mutated: Mapping[str, set[str]],
    eligible: Mapping[str, set[str]],
    genes: Sequence[str],
    min_count: int,
) -> pd.DataFrame:
    rows = []
    for gene_a, gene_b in combinations(genes, 2):
        shared = eligible.get(gene_a, set()) & eligible.get(gene_b, set())
        counts = _pair_row(
            gene_a, gene_b, mutated.get(gene_a, set()), mutated.get(gene_b, set()), shared
        )
        if counts is None:
            logger.debug(
                "%s/%s: no genome evaluated both genes; pair omitted", gene_a, gene_b
            )
            continue
        if counts["N_Both"] < min_count:
            continue
        rows.append({"Gene_A": gene_a, "Gene_B": gene_b, **counts})

    frame = pd.DataFrame(rows, columns=list(GENE_COLUMNS))
    if not frame.empty:
        frame = frame.sort_values(
            ["N_Both", "Gene_A", "Gene_B"], ascending=[False, True, True], kind="stable"
        ).reset_index(drop=True)
    return frame


def _variant_pairs(
    mutations: pd.DataFrame,
    eligible: Mapping[str, set[str]],
    min_count: int,
) -> pd.DataFrame:
    """Pairs at ``Gene:Mutation`` resolution.

    Eligibility is still per gene: two variants of genes A and B are compared
    over the genomes where A and B were both evaluated. Two variants of the
    *same* gene are compared over that gene's own eligible set.
    """
    if mutations.empty:
        return pd.DataFrame(columns=list(VARIANT_COLUMNS))

    labelled = mutations.assign(
        _variant=mutations["Gene"].astype(str) + ":" + mutations["Mutation"].astype(str)
    )
    mutated = _mutated_sets(labelled, "_variant")
    gene_of = dict(zip(labelled["_variant"], labelled["Gene"].astype(str)))

    rows = []
    for variant_a, variant_b in combinations(sorted(mutated), 2):
        gene_a, gene_b = gene_of[variant_a], gene_of[variant_b]
        shared = (
            eligible.get(gene_a, set())
            if gene_a == gene_b
            else eligible.get(gene_a, set()) & eligible.get(gene_b, set())
        )
        counts = _pair_row(variant_a, variant_b, mutated[variant_a], mutated[variant_b], shared)
        if counts is None or counts["N_Both"] < min_count:
            continue
        rows.append(
            {
                "Variant_A": variant_a,
                "Variant_B": variant_b,
                "Gene_A": gene_a,
                "Gene_B": gene_b,
                **counts,
            }
        )

    frame = pd.DataFrame(rows, columns=list(VARIANT_COLUMNS))
    if not frame.empty:
        frame = frame.sort_values(
            ["N_Both", "Variant_A", "Variant_B"],
            ascending=[False, True, True],
            kind="stable",
        ).reset_index(drop=True)
    return frame


def _matrix(
    mutated: Mapping[str, set[str]],
    eligible: Mapping[str, set[str]],
    genes: Sequence[str],
) -> pd.DataFrame:
    """Gene x gene counts. Diagonal = genomes mutating that gene.

    Off-diagonal cells use the pair's own shared denominator, so the matrix is
    symmetric but its cells are not all comparable to one another -- read
    ``N_Eligible`` in the long table before comparing two pairs.
    """
    matrix = pd.DataFrame(0, index=list(genes), columns=list(genes), dtype="int64")
    for gene in genes:
        matrix.loc[gene, gene] = len(mutated.get(gene, set()) & eligible.get(gene, set()))
    for gene_a, gene_b in combinations(genes, 2):
        shared = eligible.get(gene_a, set()) & eligible.get(gene_b, set())
        both = len(mutated.get(gene_a, set()) & mutated.get(gene_b, set()) & shared)
        matrix.loc[gene_a, gene_b] = both
        matrix.loc[gene_b, gene_a] = both
    matrix.index.name = "Gene"
    return matrix


def cooccurrence(
    mutations: pd.DataFrame,
    qc: pd.DataFrame | None = None,
    levels: Iterable[str] = ("gene",),
    min_count: int = 1,
    genes_csv: Path | str | None = None,
    matrix_csv: Path | str | None = None,
    variants_csv: Path | str | None = None,
) -> CooccurrenceResult:
    """Count how often pairs of genes (or variants) are mutated in the same genome.

    Args:
        mutations: the variant stage's mutation table.
        qc: the variant stage's QC table, which defines eligibility. Without it
            the denominator falls back to the genomes present in *mutations*,
            which **overstates** the fractions -- a warning is logged.
        levels: ``"gene"`` and/or ``"variant"``.
        min_count: drop pairs co-mutated in fewer than this many genomes.
        genes_csv / matrix_csv / variants_csv: optional output paths.

    Returns:
        :class:`CooccurrenceResult`. Counts only; nothing here is a claim about
        interaction, selection or causation.
    """
    levels = tuple(str(level).strip().lower() for level in levels)
    unknown = [level for level in levels if level not in ("gene", "variant")]
    if unknown:
        raise ValueError(f"Unknown co-occurrence level(s): {', '.join(unknown)}")
    min_count = max(0, int(min_count))

    required = ("Accession", "Gene", "Mutation")
    missing = [column for column in required if column not in mutations.columns]
    if missing and not mutations.empty:
        raise ValueError(f"Mutations table is missing column(s): {', '.join(missing)}")
    if missing:
        mutations = pd.DataFrame(columns=list(required))

    mutations = mutations.astype({"Accession": str, "Gene": str, "Mutation": str})

    eligible = eligibility_from_qc(qc) if qc is not None else {}
    if not eligible:
        logger.warning(
            "No variant QC table supplied: eligibility falls back to the genomes "
            "that carry at least one mutation, so fractions will be overstated"
        )
        observed = {str(accession) for accession in mutations["Accession"]}
        eligible = {gene: set(observed) for gene in mutations["Gene"].unique()}

    genes = sorted(set(eligible) | set(mutations["Gene"].unique()))
    mutated_by_gene = _mutated_sets(mutations, "Gene")

    genes_long = (
        _gene_pairs(mutated_by_gene, eligible, genes, min_count)
        if "gene" in levels
        else pd.DataFrame(columns=list(GENE_COLUMNS))
    )
    matrix = _matrix(mutated_by_gene, eligible, genes)
    variants_long = (
        _variant_pairs(mutations, eligible, min_count) if "variant" in levels else None
    )

    result = CooccurrenceResult(
        genes=genes_long,
        genes_matrix=matrix,
        variants=variants_long,
        gene_totals={gene: int(matrix.loc[gene, gene]) for gene in genes},
        eligible_per_gene={gene: len(eligible.get(gene, set())) for gene in genes},
    )

    _write(result, genes_csv, matrix_csv, variants_csv)
    _log_summary(result)
    return result


def _write(
    result: CooccurrenceResult,
    genes_csv: Path | str | None,
    matrix_csv: Path | str | None,
    variants_csv: Path | str | None,
) -> None:
    for frame, path, index in (
        (result.genes, genes_csv, False),
        (result.genes_matrix, matrix_csv, True),
        (result.variants, variants_csv, False),
    ):
        if path is None or frame is None:
            continue
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(path, index=index)
        logger.info("Wrote %s (%d rows)", path, len(frame))


def _log_summary(result: CooccurrenceResult) -> None:
    if not result.gene_totals:
        logger.warning("No genes were evaluated; co-occurrence is empty")
        return

    logger.info(
        "Co-occurrence over %d gene(s): %s",
        len(result.gene_totals),
        ", ".join(
            f"{gene} mutated in {count}/{result.eligible_per_gene.get(gene, 0)}"
            for gene, count in result.gene_totals.items()
        ),
    )
    logger.info("  %d gene pair(s) reported", len(result.genes))
    if result.variants is not None:
        logger.info("  %d variant pair(s) reported", len(result.variants))


def read_cooccurrence(path: Path | str) -> pd.DataFrame:
    """Read a long-form co-occurrence table, keeping gene names as strings."""
    frame = pd.read_csv(path)
    for column in ("Gene_A", "Gene_B", "Variant_A", "Variant_B"):
        if column in frame.columns:
            frame[column] = frame[column].astype(str)
    return frame


def load_result(
    genes_csv: Path | str,
    matrix_csv: Path | str,
    variants_csv: Path | str | None = None,
    qc: pd.DataFrame | None = None,
) -> CooccurrenceResult:
    """Rebuild a :class:`CooccurrenceResult` from files a previous stage wrote.

    Lets a later step (the run summary) report on co-occurrence without recomputing
    it, which is what keeps the Snakemake wrappers free of logic.
    """
    matrix = pd.read_csv(matrix_csv, index_col=0)
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)

    variants = (
        read_cooccurrence(variants_csv)
        if variants_csv is not None and Path(variants_csv).is_file()
        else None
    )
    eligible = eligibility_from_qc(qc) if qc is not None else {}

    return CooccurrenceResult(
        genes=read_cooccurrence(genes_csv),
        genes_matrix=matrix,
        variants=variants,
        gene_totals={gene: int(matrix.loc[gene, gene]) for gene in matrix.index},
        eligible_per_gene={
            gene: len(eligible.get(gene, set())) for gene in matrix.index
        },
    )
