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

from .parallel import WORKLOAD_PROCESS, parallel_map
from .scoring import build_severity_index, pair_severity

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
    if mutations.empty:
        return pd.DataFrame(columns=list(GENOME_MAP_COLUMNS))

    called_pairs: set[tuple[str, str]] = set()
    if not qc.empty and "Status" in qc.columns:
        called_pairs = {
            (str(row["Accession"]), str(row["Gene"]))
            for _, row in qc[qc["Status"] == "called"].iterrows()
        }

    rows = []
    for accession, group in mutations.groupby("Accession", sort=True):
        filtered = group
        if called_pairs:
            filtered = group[
                group.apply(
                    lambda r: (str(r["Accession"]), str(r["Gene"])) in called_pairs,
                    axis=1,
                )
            ]
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


def _pair_epistasis(n11: int, n10: int, n01: int, n00: int) -> dict:
    """Fisher's exact test for one mutation pair.

    Handles degenerate tables (zero cells) by computing Fisher's exact on
    the integer table when possible, and falling back to a Haldane-Anscombe
    correction (+0.5 to every cell) for the p-value when the table has a
    structural zero. The odds ratio is computed directly from the (possibly
    corrected) counts so it never returns NaN.
    """
    n = n11 + n10 + n01 + n00

    # Try the raw table first
    table = [[n11, n10], [n01, n00]]
    odds_ratio, p_value = fisher_exact(table, alternative="two-sided")

    if np.isnan(odds_ratio):
        # Haldane-Anscombe correction for structural zeros
        c11, c10, c01, c00 = n11 + 0.5, n10 + 0.5, n01 + 0.5, n00 + 0.5
        _, p_value = fisher_exact([[c11, c10], [c01, c00]], alternative="two-sided")
        denom = c10 * c01
        odds_ratio = (c11 * c00) / denom if denom else float("inf")

    expected = ((n11 + n10) * (n11 + n01)) / n if n else 0.0
    return {
        "N_Both": n11,
        "N_A_only": n10,
        "N_B_only": n01,
        "N_Neither": n00,
        "N_Eligible": n,
        "Observed": n11,
        "Expected": round(expected, 4),
        "Odds_Ratio": round(float(odds_ratio), 4),
        "P_value": float(p_value),
    }


#: Below this many pairs, test them in-process.
#:
#: A Fisher exact test is not the cheap call it looks like: measured on this
#: workload it costs ~160us for a balanced small table and ~1.5ms for a rare
#: co-mutation in a large cohort (``benchmarks/parallel_bench.py --verify-thresholds``
#: re-checks this). Even at the cheap end, 2048 pairs is ~0.3s of work against a
#: ~0.2s pool start-up under the ``spawn`` start method that Windows forces, so
#: this threshold sits just past break-even even in the most favourable case, and
#: well past it for anything but a small balanced cohort.
#:
#: Erring low is the safe direction: running serially when a pool might have
#: given 1.5x costs almost nothing, whereas starting a pool that then gains
#: nothing costs the same 0.2s every time.
_SERIAL_PAIR_LIMIT = 2048

#: Per-process epistasis state, installed once by the pool initializer.
#:
#: A closure would be simpler but cannot be pickled, and shipping the three
#: tables to every task would pickle the whole cohort once per pair. Handing them
#: to the initializer instead costs one copy per *worker*.
_EPISTASIS_STATE: dict[str, object] = {}


def _init_epistasis_worker(
    eligible: dict[str, set[str]],
    mutated: dict[str, set[str]],
    gene_of: dict[str, str],
) -> None:
    """Pool initializer: install the shared tables in this worker."""
    _EPISTASIS_STATE["eligible"] = eligible
    _EPISTASIS_STATE["mutated"] = mutated
    _EPISTASIS_STATE["gene_of"] = gene_of


def _test_variant_pair(pair: tuple[str, str]) -> dict | None:
    """Build the 2x2 table for one variant pair and run Fisher's exact test.

    Module level, and therefore picklable, so a process pool can call it.
    Returns ``None`` for a pair with no shared denominator, which is the same
    signal the serial path used.
    """
    eligible = _EPISTASIS_STATE["eligible"]
    mutated = _EPISTASIS_STATE["mutated"]
    gene_of = _EPISTASIS_STATE["gene_of"]
    assert isinstance(eligible, dict) and isinstance(mutated, dict)
    assert isinstance(gene_of, dict)

    variant_a, variant_b = pair
    gene_a = gene_of[variant_a]
    gene_b = gene_of[variant_b]
    if gene_a == gene_b:
        shared = eligible.get(gene_a, set())
    else:
        shared = eligible.get(gene_a, set()) & eligible.get(gene_b, set())
    if not shared:
        return None

    in_a = mutated[variant_a] & shared
    in_b = mutated[variant_b] & shared
    n11 = len(in_a & in_b)
    n10 = len(in_a - in_b)
    n01 = len(in_b - in_a)
    n00 = len(shared) - n11 - n10 - n01

    result = _pair_epistasis(n11, n10, n01, n00)
    result["Variant_A"] = variant_a
    result["Variant_B"] = variant_b
    result["Gene_A"] = gene_a
    result["Gene_B"] = gene_b
    return result


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

    Args:
        mutations: mutations table with Accession, Gene, Mutation, Position.
        qc: variant QC table with Accession, Gene, Status.
        severity: per-mutation severity scores from :func:`scoring.mutation_severity`.
            If ``None``, pairs are ranked only by count and p-value.
        fdr_threshold: significance cut after FDR correction.
        min_count: drop pairs with ``N_Both`` below this threshold.
        threads: parallel worker processes for the pairwise tests.

    Returns:
        Ranked epistasis network DataFrame.
    """
    required = {"Accession", "Gene", "Mutation"}
    if not mutations.empty and not required.issubset(mutations.columns):
        missing = sorted(required - set(mutations.columns))
        raise ValueError(f"mutations table missing column(s): {', '.join(missing)}")

    if mutations.empty or len(mutations) < 2:
        return pd.DataFrame(columns=list(EPISTASIS_COLUMNS))

    # Eligibility sets from QC: gene -> set of accessions where it was called
    eligible: dict[str, set[str]] = {}
    if not qc.empty and "Status" in qc.columns:
        called = qc[qc["Status"] == "called"]
        for _, row in called.iterrows():
            eligible.setdefault(str(row["Gene"]), set()).add(str(row["Accession"]))

    # Mutation presence: variant_key -> set of accessions
    labelled = mutations.copy()
    labelled["_variant"] = labelled["Gene"].astype(str) + ":" + labelled["Mutation"].astype(str)
    mutated: dict[str, set[str]] = {}
    for variant, group in labelled.groupby("_variant", sort=True):
        mutated[variant] = set(group["Accession"].astype(str))

    gene_of = dict(zip(labelled["_variant"], labelled["Gene"].astype(str), strict=False))

    # Pre-build severity lookup if provided
    sev_lookup: dict[str, float] = {}
    if severity is not None and not severity.empty:
        for _, row in severity.iterrows():
            key = f"{row['Gene']}:{row['Mutation']}"
            sev_lookup[key] = float(row["severity_score"])

    variants = sorted(mutated.keys())
    pairs = list(combinations(variants, 2))

    use_pool = threads > 1 and len(pairs) > _SERIAL_PAIR_LIMIT
    if threads > 1 and not use_pool:
        logger.info(
            "Testing %d variant pair(s) in-process: below the %d-pair threshold "
            "where pool start-up pays for itself",
            len(pairs),
            _SERIAL_PAIR_LIMIT,
        )
    else:
        logger.info("Computing epistasis for %d variant pairs on %d thread(s)", len(pairs), threads)

    results: list[dict] = []
    # Processes, not threads: scipy's exact tests and the set arithmetic are
    # CPU-bound Python, and the GIL would serialise a thread pool. The shared
    # tables ride along in initargs, so they are sent once per worker rather than
    # once per task.
    for result in parallel_map(
        _test_variant_pair,
        pairs,
        threads=threads if use_pool else 1,
        workload=WORKLOAD_PROCESS,
        total=len(pairs),
        initializer=_init_epistasis_worker,
        initargs=(eligible, mutated, gene_of),
        label="epistasis",
    ):
        if result is not None:
            results.append(result)

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
    if sev_lookup:
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
) -> dict:
    """Run the full epistasis pipeline.

    1. Build a per-genome mutation map.
    2. Compute per-mutation biochemical severity.
    3. Detect epistatic pairs and rank them.

    Args:
        mutations_csv: path to ``mutations.csv``.
        qc_csv: path to ``variant_qc.csv`` (optional but recommended).
        output_dir: destination folder. Created if absent.
        fdr_threshold: Benjamini-Hochberg significance threshold.
        min_count: minimum ``N_Both`` to retain a pair.
        threads: parallel threads for pairwise tests.

    Returns:
        Dictionary with paths to the three output files.
    """
    mutations = pd.read_csv(mutations_csv)
    qc = pd.read_csv(qc_csv) if qc_csv and Path(qc_csv).is_file() else pd.DataFrame()

    logger.info("Loaded %d mutations for epistasis analysis", len(mutations))

    genome_map = genome_mutation_map(mutations, qc)

    from .scoring import mutation_severity

    severity = mutation_severity(mutations)
    networks = detect_epistasis(mutations, qc, severity, fdr_threshold, min_count, threads)

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
