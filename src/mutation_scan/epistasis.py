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
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from .scoring import pair_severity

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
        threads: parallel threads for pairwise tests.

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

    logger.info(
        "Computing epistasis for %d variant pairs on %d thread(s)",
        len(pairs),
        threads,
    )

    def _work(va: str, vb: str) -> dict | None:
        ga = gene_of[va]
        gb = gene_of[vb]
        if ga == gb:
            shared = eligible.get(ga, set())
        else:
            shared = eligible.get(ga, set()) & eligible.get(gb, set())
        if not shared:
            return None

        in_a = mutated[va] & shared
        in_b = mutated[vb] & shared
        n11 = len(in_a & in_b)
        n10 = len(in_a - in_b)
        n01 = len(in_b - in_a)
        n00 = len(shared) - n11 - n10 - n01

        res = _pair_epistasis(n11, n10, n01, n00)
        res["Variant_A"] = va
        res["Variant_B"] = vb
        res["Gene_A"] = ga
        res["Gene_B"] = gb
        return res

    results: list[dict] = []
    if threads > 1 and len(pairs) > 200:
        with ThreadPoolExecutor(max_workers=threads) as pool:
            futures = {pool.submit(_work, va, vb): (va, vb) for va, vb in pairs}
            for future in as_completed(futures):
                try:
                    res = future.result()
                    if res is not None:
                        results.append(res)
                except Exception as exc:
                    va, vb = futures[future]
                    logger.error("epistasis pair %s / %s failed: %s", va, vb, exc)
    else:
        for va, vb in pairs:
            res = _work(va, vb)
            if res is not None:
                results.append(res)

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
        pair_sevs = []
        for _, row in frame.iterrows():
            pair_sevs.append(pair_severity(row["Variant_A"], row["Variant_B"], severity))
        frame["Pair_Severity"] = pair_sevs

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
