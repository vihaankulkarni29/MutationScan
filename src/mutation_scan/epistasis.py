"""Epistasis detection and genome mutation mapping.

Two outputs:

1. **Genome mutation map** — a minimalistic per-genome view of every mutation.
2. **Epistasis networks** — statistically tested mutation pairs ranked by
   biochemical severity.

No structural modelling, no docking, no PDB parsing, no Grantham scoring.
Statistical testing uses Fisher's exact test with Benjamini-Hochberg FDR.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from .scoring import build_severity_index, pair_severity
from .variants import LOF_ALLELE, with_loss_of_function

logger = logging.getLogger(__name__)

__all__ = [
    "EPISTASIS_COLUMNS",
    "GENOME_MAP_COLUMNS",
    "EpistasisResult",
    "detect_epistasis",
    "epistasis_networks",
    "genome_mutation_map",
    "read_epistasis",
]

GENOME_MAP_COLUMNS = (
    "Accession",
    "N_Mutated_Genes",
    "Mutated_Genes",
    "N_Mutations",
    "N_LoF",
    "Mutation_Map",
)

EPISTASIS_COLUMNS = (
    "Variant_A",
    "Variant_B",
    "Gene_A",
    "Gene_B",
    "N_Both",
    "N_A_only",
    "N_B_only",
    "N_Neither",
    "N_Eligible",
    "Observed",
    "Expected",
    "Odds_Ratio",
    "P_value",
    "Q_value",
    "Sign",
    "Pair_Severity",
)


@dataclass
class EpistasisResult:
    """Container for epistasis outputs."""

    genome_map: pd.DataFrame
    networks: pd.DataFrame


def genome_mutation_map(
    mutations: pd.DataFrame,
    qc: pd.DataFrame,
) -> pd.DataFrame:
    """Build a minimalistic per-genome mutation map.

    One row per genome. Lists every mutation as ``Gene:Mutation`` in a
    single pipe-separated string so the table stays narrow and grep-friendly.

    Args:
        mutations: mutations table (Accession, Gene, Mutation, Position)
        qc: variant QC table (used to restrict to evaluated (gene, genome) pairs)

    Returns:
        One-row-per-genome DataFrame.
    """
    mutations = with_loss_of_function(mutations, qc)
    if mutations.empty:
        return pd.DataFrame(columns=list(GENOME_MAP_COLUMNS))

    called_pairs: set[tuple[str, str]] = set()
    if not qc.empty and "Status" in qc.columns:
        called = qc[qc["Status"] == "called"]
        called_pairs = set(
            zip(
                called["Accession"].astype(str),
                called["Gene"].astype(str),
                strict=True,
            )
        )

    rows = []
    for accession, group in mutations.groupby("Accession", sort=True):
        filtered = group
        if called_pairs:
            # The membership test is pure Python either way; a row-wise apply
            # just adds pandas' per-row overhead on top of it. Once the pairwise
            # stage stopped looping over pairs this filter became the most
            # expensive thing left in the module.
            keep = [
                (acc, gene) in called_pairs
                for acc, gene in zip(
                    group["Accession"].astype(str),
                    group["Gene"].astype(str),
                    strict=True,
                )
            ]
            filtered = group[keep]
        if filtered.empty:
            continue
        genes = sorted(filtered["Gene"].astype(str).unique())
        ordered = filtered.sort_values(["Gene", "Position"])
        mutation_list = [f"{row.Gene}:{row.Mutation}" for row in ordered.itertuples()]
        rows.append(
            {
                "Accession": accession,
                "N_Mutated_Genes": len(genes),
                "Mutated_Genes": ";".join(genes),
                "N_Mutations": len(filtered),
                "N_LoF": int((filtered["Mutation"].astype(str) == LOF_ALLELE).sum()),
                "Mutation_Map": " | ".join(mutation_list),
            }
        )

    frame = pd.DataFrame(rows, columns=list(GENOME_MAP_COLUMNS))
    frame = frame.sort_values("Accession", kind="stable").reset_index(drop=True)
    return frame


def _benjamini_hochberg(p_values: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg FDR correction. Returns q-values in original order."""
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    if n == 0:
        return np.array([], dtype=float)

    order = np.argsort(p)
    sorted_p = p[order]
    ranks = np.arange(1, n + 1)
    q = sorted_p * n / ranks
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)

    q_values = np.empty(n, dtype=float)
    q_values[order] = q
    return q_values


#: Per-pair epistasis work is resolved in two exact stages -- a Gram matrix for
#: the co-occurrence counts, and one Fisher test per *distinct* 2x2 table -- so
#: the stage no longer spawns a process pool. `threads` is still accepted so
#: existing scripts and config files keep working; it no longer changes what runs.
#:
#: The escalation path for much larger cohorts is to pool over distinct tables
#: rather than pairs, since that is now the parallelisable unit and it has a
#: flat per-task payload. At 541 variants over 570 genomes there are 146,070
#: pairs but only ~1,700 distinct tables after min_count, so the vectorised path
#: is not close to needing it.


@dataclass(frozen=True)
class _PairTable:
    """One pair's 2x2 contingency counts, resolved but not yet tested."""

    index: int
    variant_a: str
    variant_b: str
    gene_a: str
    gene_b: str
    n11: int
    n10: int
    n01: int
    n00: int


def _universe_for(
    gene_a: str,
    gene_b: str,
    eligible: dict[str, set[str]],
) -> frozenset[str] | None:
    """The genome set a pair is scored against, or ``None`` if it is empty.

    Depends only on the two genes, so the number of distinct universes is
    bounded by the panel size rather than by the number of pairs.
    """
    if gene_a == gene_b:
        shared = eligible.get(gene_a, set())
    else:
        shared = eligible.get(gene_a, set()) & eligible.get(gene_b, set())
    return frozenset(shared) if shared else None


def _resolve_pair_tables(
    variants: list[str],
    pairs: list[tuple[str, str]],
    mutated: dict[str, set[str]],
    gene_of: dict[str, str],
    eligible: dict[str, set[str]],
) -> list[_PairTable]:
    """Build every pair's 2x2 counts, batching each universe into one matmul.

    The counts are pure arithmetic once the co-occurrence numbers are known:
    ``n11`` is a set intersection, and ``n10``/``n01``/``n00`` follow from the
    marginals. ``n11`` for *all* pairs in a universe is the Gram matrix of that
    universe's 0/1 presence matrix, which BLAS computes in one pass instead of
    once per pair.

    Pairs are returned in the order given, which the caller relies on: the final
    sort is stable, so equal ``(Q_value, Pair_Severity, N_Both)`` triples fall
    back to this enumeration order to break the tie.
    """
    variant_index = {variant: i for i, variant in enumerate(variants)}
    gene_of_pair = [(gene_of[a], gene_of[b]) for a, b in pairs]

    # A pair's universe depends only on its two genes, so memoising on the gene
    # pair turns 146k set intersections into at most panel_size**2 of them.
    universe_cache: dict[tuple[str, str], frozenset[str] | None] = {}
    by_universe: dict[frozenset[str], list[int]] = {}
    universe_at: list[frozenset[str] | None] = []
    for position, (gene_a, gene_b) in enumerate(gene_of_pair):
        gene_pair = (gene_a, gene_b)
        if gene_pair not in universe_cache:
            universe_cache[gene_pair] = _universe_for(gene_a, gene_b, eligible)
        shared = universe_cache[gene_pair]
        universe_at.append(shared)
        if shared is not None:
            by_universe.setdefault(shared, []).append(position)

    # One global presence matrix, sliced per universe. Building it per universe
    # instead would re-walk every variant's accession set once per universe.
    # Eligible genomes that carry no mutation still need a row: they belong to
    # universes and contribute an all-zero row, which is exactly their weight.
    accessions = sorted(
        {a for carriers in mutated.values() for a in carriers}
        | {a for gene_set in eligible.values() for a in gene_set}
    )
    acc_index = {accession: i for i, accession in enumerate(accessions)}
    presence_all = np.zeros((len(accessions), len(variants)), dtype=np.float64)
    for variant, carriers in mutated.items():
        presence_all[[acc_index[a] for a in carriers], variant_index[variant]] = 1.0

    n11_at: dict[int, int] = {}
    n10_at: dict[int, int] = {}
    n01_at: dict[int, int] = {}
    for shared, positions in by_universe.items():
        rows = [acc_index[a] for a in sorted(shared)]
        presence = presence_all[rows]

        # float64 is exact here (counts are bounded by the cohort size, far
        # below 2**53), so the Gram matrix is bit-identical whatever thread
        # count BLAS happens to use. G.T @ G is the variant-by-variant Gram
        # matrix, i.e. the co-occurrence count for every pair at once.
        cooccurrence = presence.T @ presence

        # Column sums give each variant's carrier count within this universe, so
        # n10 and n01 fall out of the intersection without further set work.
        marginals = presence.sum(axis=0)
        for position in positions:
            variant_a, variant_b = pairs[position]
            i, j = variant_index[variant_a], variant_index[variant_b]
            both = int(cooccurrence[i, j])
            n11_at[position] = both
            n10_at[position] = int(marginals[i]) - both
            n01_at[position] = int(marginals[j]) - both

    tables: list[_PairTable] = []
    for position, (variant_a, variant_b) in enumerate(pairs):
        shared = universe_at[position]
        if shared is None:
            continue
        n11, n10, n01 = n11_at[position], n10_at[position], n01_at[position]
        tables.append(
            _PairTable(
                index=position,
                variant_a=variant_a,
                variant_b=variant_b,
                gene_a=gene_of_pair[position][0],
                gene_b=gene_of_pair[position][1],
                n11=n11,
                n10=n10,
                n01=n01,
                n00=len(shared) - n11 - n10 - n01,
            )
        )

    return tables


def _test_distinct_tables(
    distinct: list[tuple[int, int, int, int]],
) -> dict[tuple[int, int, int, int], tuple[float, float]]:
    """Run Fisher's exact test once per *distinct* 2x2 table.

    The test is a pure function of the integer counts, so pairs sharing a table
    necessarily share a p-value and odds ratio. Memoising on the tuple is
    therefore exact, and at cohort scale it collapses the number of tests by
    more than an order of magnitude.
    """
    results: dict[tuple[int, int, int, int], tuple[float, float]] = {}
    for table in distinct:
        if table in results:
            continue
        n11, n10, n01, n00 = table
        odds_ratio, p_value = fisher_exact([[n11, n10], [n01, n00]], alternative="two-sided")
        if np.isnan(odds_ratio):
            # Haldane-Anscombe correction for structural zeros.
            c11, c10, c01, c00 = n11 + 0.5, n10 + 0.5, n01 + 0.5, n00 + 0.5
            _, p_value = fisher_exact([[c11, c10], [c01, c00]], alternative="two-sided")
            denom = c10 * c01
            odds_ratio = (c11 * c00) / denom if denom else float("inf")
        results[table] = (float(odds_ratio), float(p_value))
    return results


def detect_epistasis(
    mutations: pd.DataFrame,
    qc: pd.DataFrame,
    severity: pd.DataFrame | None = None,
    fdr_threshold: float = 0.05,
    min_count: int = 1,
    threads: int = 1,
) -> pd.DataFrame:
    """Detect epistatic mutation pairs and rank by biochemical severity.

    For every unordered pair of mutations:

    1. Restrict to genomes where **both** mutations were evaluated (the same
       honest denominator used by :func:`mutation_scan.cooccurrence`).
    2. Build a 2×2 contingency table and run Fisher's exact test.
    3. Apply Benjamini-Hochberg FDR correction.
    4. Classify as *positive*, *negative*, or *neutral*.
    5. Rank by ``Q_value`` → ``Pair_Severity`` ↓ → ``N_Both`` ↓.

    The pairwise work is resolved in three stages rather than as a loop over
    pairs. A pair's denominator depends only on its two genes, so the pairs are
    grouped into a handful of universes and each universe's co-occurrence counts
    come from a single matrix product. The counts are then keyed by their 2×2
    tuple so that Fisher's exact test runs once per *distinct* table rather than
    once per pair. Both steps are exact rearrangements, not approximations: see
    :func:`_resolve_pair_tables` and :func:`_test_distinct_tables`.

    Args:
        mutations: mutations table with Accession, Gene, Mutation, Position.
        qc: variant QC table with Accession, Gene, Status.
        severity: per-mutation severity scores from :func:`scoring.mutation_severity`.
            If ``None``, pairs are ranked only by count and p-value.
        fdr_threshold: significance cut after FDR correction.
        min_count: drop pairs with ``N_Both`` below this threshold.
        threads: retained for interface compatibility. The stage no longer
            spawns a process pool; see the note in the module plan.

    Returns:
        Ranked epistasis network DataFrame.
    """
    required = {"Accession", "Gene", "Mutation"}
    if not mutations.empty and not required.issubset(mutations.columns):
        missing = sorted(required - set(mutations.columns))
        raise ValueError(f"mutations table missing column(s): {', '.join(missing)}")

    # LoF alleles collapse into the substitution stream here so a direct
    # caller of detect_epistasis (one passing raw mutations + a QC table
    # that carries event columns) still sees LoF pairs. epistasis_networks
    # merges once for mutation_severity; with_loss_of_function is idempotent,
    # so the merge here is a no-op when the frame already carries LoF rows.
    mutations = with_loss_of_function(mutations, qc)

    if mutations.empty or len(mutations) < 2:
        return pd.DataFrame(columns=list(EPISTASIS_COLUMNS))

    # Eligibility sets from QC: gene -> set of accessions where it was called
    eligible: dict[str, set[str]] = {}
    if not qc.empty and "Status" in qc.columns:
        called = qc[qc["Status"] == "called"]
        for accession, gene in zip(
            called["Accession"].astype(str), called["Gene"].astype(str), strict=True
        ):
            eligible.setdefault(gene, set()).add(accession)

    # Mutation presence: variant_key -> set of accessions
    labelled = mutations.copy()
    labelled["_variant"] = labelled["Gene"].astype(str) + ":" + labelled["Mutation"].astype(str)
    mutated: dict[str, set[str]] = {}
    for variant, group in labelled.groupby("_variant", sort=True):
        mutated[variant] = set(group["Accession"].astype(str))

    gene_of = dict(zip(labelled["_variant"], labelled["Gene"].astype(str), strict=False))

    # Whether severity scoring was requested at all. Deliberately not a lookup
    # dict: scoring itself goes through build_severity_index, which averages a
    # variant's severity across every genome carrying it. Building a per-variant
    # dict here would have kept whichever row came last instead, and cost an
    # iterrows pass over the whole severity table to be thrown away.
    score_pairs = severity is not None and not severity.empty

    variants = sorted(mutated.keys())
    pairs = list(combinations(variants, 2))

    logger.info("Computing epistasis for %d variant pairs", len(pairs))

    tables = _resolve_pair_tables(variants, pairs, mutated, gene_of, eligible)
    testable = len(tables)

    # Apply min_count here rather than after building the frame. The filter runs
    # before BH either way, so this is the same set of rows -- but a pair with
    # N_Both == 0 can never reach the output, so there is no reason to test it.
    # At cohort scale that is the majority of pairs, and the saving is exact.
    tables = [table for table in tables if table.n11 >= min_count]

    distinct = sorted({(t.n11, t.n10, t.n01, t.n00) for t in tables})
    logger.info(
        "Resolved %d testable pair(s), %d above min_count, across %d distinct 2x2 table(s)",
        testable,
        len(tables),
        len(distinct),
    )

    tested = _test_distinct_tables(distinct)

    results: list[dict] = []
    for table in tables:
        odds_ratio, p_value = tested[(table.n11, table.n10, table.n01, table.n00)]
        n = table.n11 + table.n10 + table.n01 + table.n00
        expected = ((table.n11 + table.n10) * (table.n11 + table.n01)) / n if n else 0.0
        results.append(
            {
                "Variant_A": table.variant_a,
                "Variant_B": table.variant_b,
                "Gene_A": table.gene_a,
                "Gene_B": table.gene_b,
                "N_Both": table.n11,
                "N_A_only": table.n10,
                "N_B_only": table.n01,
                "N_Neither": table.n00,
                "N_Eligible": n,
                "Observed": table.n11,
                "Expected": round(expected, 4),
                "Odds_Ratio": round(odds_ratio, 4),
                "P_value": p_value,
            }
        )

    if not results:
        return pd.DataFrame(columns=list(EPISTASIS_COLUMNS))

    frame = pd.DataFrame(results)
    frame = frame[frame["N_Both"] >= min_count].copy()

    if frame.empty:
        return frame

    p_values = frame["P_value"].values
    frame["Q_value"] = _benjamini_hochberg(p_values)

    frame["Sign"] = "neutral"
    frame.loc[(frame["Odds_Ratio"] > 1) & (frame["Q_value"] < fdr_threshold), "Sign"] = "positive"
    frame.loc[(frame["Odds_Ratio"] < 1) & (frame["Q_value"] < fdr_threshold), "Sign"] = "negative"

    frame["Pair_Severity"] = 0.0
    if score_pairs:
        # Build the (gene, mutation) -> mean severity lookup once. Scoring each
        # pair with its own scan of the severity table is quadratic, and at
        # cohort scale that scan outweighed the Fisher tests themselves.
        index = build_severity_index(severity)
        frame["Pair_Severity"] = [
            pair_severity(variant_a, variant_b, severity, index)
            for variant_a, variant_b in zip(frame["Variant_A"], frame["Variant_B"], strict=True)
        ]

    frame = frame.sort_values(
        ["Q_value", "Pair_Severity", "N_Both"],
        ascending=[True, False, False],
        kind="stable",
    ).reset_index(drop=True)

    return frame[[c for c in EPISTASIS_COLUMNS if c in frame.columns]]


def epistasis_networks(
    mutations_csv: Path | str,
    qc_csv: Path | str | None = None,
    output_dir: Path | str | None = None,
    fdr_threshold: float = 0.05,
    min_count: int = 1,
    threads: int = 1,
    lineage_manifest: Path | str | None = None,
) -> dict:
    """Run the full epistasis pipeline.

    1. Build a per-genome mutation map.
    2. Compute per-mutation biochemical severity.
    3. Detect epistatic pairs and rank them.
    4. *(optional)* annotate with lineage information.

    Args:
        mutations_csv: path to ``mutations.csv``.
        qc_csv: path to ``variant_qc.csv`` (optional but recommended).
        output_dir: destination folder. Created if absent.
        fdr_threshold: Benjamini-Hochberg significance threshold.
        min_count: minimum ``N_Both`` to retain a pair.
        threads: parallel threads for pairwise tests.
        lineage_manifest: optional path to a TSV/CSV mapping Accession ->
            Lineage. When provided, the genome map gains a ``Lineage``
            column and the networks table gains ``Lineage_Distribution``
            and ``Clonality`` columns.

    Returns:
        Dictionary with paths to the output files.
    """
    mutations = pd.read_csv(mutations_csv)
    qc = pd.read_csv(qc_csv) if qc_csv and Path(qc_csv).is_file() else pd.DataFrame()

    mutations = with_loss_of_function(mutations, qc)
    logger.info(
        "Loaded %d mutations for epistasis analysis (%d loss-of-function)",
        len(mutations),
        int((mutations["Mutation"].astype(str) == LOF_ALLELE).sum()) if not mutations.empty else 0,
    )

    genome_map = genome_mutation_map(mutations, qc)

    from .scoring import mutation_severity

    severity = mutation_severity(mutations)
    networks = detect_epistasis(mutations, qc, severity, fdr_threshold, min_count, threads)

    if lineage_manifest:
        from .lineage import (
            annotate_epistasis,
            annotate_genome_map,
            read_lineage_manifest,
        )

        lineage_map = read_lineage_manifest(lineage_manifest)
        if lineage_map:
            genome_map = annotate_genome_map(genome_map, lineage_map)
            networks = annotate_epistasis(networks, mutations, lineage_map)
            logger.info(
                "Annotated epistasis outputs with lineage data (%d genomes mapped)",
                len(lineage_map),
            )

    output_dir = Path(output_dir) if output_dir else Path(".")
    output_dir.mkdir(parents=True, exist_ok=True)

    genome_path = output_dir / "genome_mutation_map.csv"
    severity_path = output_dir / "mutation_severity.csv"
    networks_path = output_dir / "epistasis_networks.csv"

    genome_map.to_csv(genome_path, index=False)
    severity.to_csv(severity_path, index=False)
    networks.to_csv(networks_path, index=False)

    logger.info("Wrote genome map: %s (%d genomes)", genome_path, len(genome_map))
    logger.info("Wrote severity: %s (%d mutations)", severity_path, len(severity))
    logger.info("Wrote networks: %s (%d pairs)", networks_path, len(networks))

    return {
        "genome_map": str(genome_path),
        "severity": str(severity_path),
        "networks": str(networks_path),
    }


def read_epistasis(path: Path | str) -> pd.DataFrame:
    """Read an epistasis networks CSV, preserving types."""
    frame = pd.read_csv(path)
    for col in ("Variant_A", "Variant_B", "Gene_A", "Gene_B", "Sign"):
        if col in frame.columns:
            frame[col] = frame[col].astype(str)
    return frame
