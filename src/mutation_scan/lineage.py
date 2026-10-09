"""Lineage tagging for population-structure control.

When a cohort is dominated by a single clone (e.g. an ST131 outbreak),
genomes share mutations by descent, which inflates co-occurrence counts
and Fisher-exact p-values. Mapping each genome to a lineage (sequence
type, phylogroup, etc.) lets an epistatic pair be classified as:

* **Clonal** -- every co-mutated genome shares one lineage, so the signal
  may be descent rather than selection;
* **Convergent** -- the variants appear across two or more distinct
  lineages, the gold-standard evidence for true biological epistasis.

The lineage manifest is an external TSV/CSV supplied by the user and is
never assumed: the column names are ``Accession`` and ``Lineage`` (case-
insensitive). A genome absent from the manifest is tagged ``"unassigned"``.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from pathlib import Path

import pandas as pd

__all__ = [
    "CLONALITY_CLONAL",
    "CLONALITY_CONVERGENT",
    "CLONALITY_UNKNOWN",
    "classify_clonality",
    "annotate_epistasis",
    "annotate_genome_map",
    "lineage_distribution",
    "lineage_distribution_str",
    "read_lineage_manifest",
]

CLONALITY_UNKNOWN = "Unknown"
CLONALITY_CLONAL = "Clonal"
CLONALITY_CONVERGENT = "Convergent"

_UNASSIGNED = "unassigned"


def read_lineage_manifest(path: Path | str | None) -> dict[str, str]:
    """Read a lineage manifest mapping accession to a lineage label.

    TSV or CSV with columns ``Accession`` and ``Lineage`` (matched
    case-insensitively so ``accession``/``ACCESSION`` both work). Whitespace
    is stripped from both fields. A missing or empty file returns ``{}``.
    """
    if path is None:
        return {}
    path = Path(path)
    if not path.is_file():
        return {}
    try:
        # Try to infer separator from file extension, default to comma.
        if path.suffix.lower() in {".tsv", ".tab"}:
            frame = pd.read_csv(path, dtype=str, sep="\t")
        else:
            frame = pd.read_csv(path, dtype=str)
    except (OSError, pd.errors.ParserError):
        return {}

    cols = {c.lower(): c for c in frame.columns}
    acc_col = cols.get("accession")
    lin_col = cols.get("lineage")
    if acc_col is None or lin_col is None:
        return {}

    mapping: dict[str, str] = {}
    for accession, lineage in zip(
        frame[acc_col].astype(str), frame[lin_col].astype(str), strict=True
    ):
        acc = accession.strip()
        if acc:
            mapping[acc] = lineage.strip() if lineage.strip() else _UNASSIGNED
    return mapping


def lineage_distribution(
    accessions: Sequence[str], lineage_map: Mapping[str, str]
) -> dict[str, int]:
    """Count lineages among a set of accessions.

    Accessions missing from *lineage_map* are bucketed as ``"unassigned"``.
    Returns ``{}`` when the map is empty (no lineage data supplied).
    """
    if not lineage_map:
        return {}
    return dict(Counter(lineage_map.get(acc, _UNASSIGNED) for acc in accessions))


def classify_clonality(accessions: Sequence[str], lineage_map: Mapping[str, str]) -> str:
    """Classify a group by lineage diversity.

    * Two or more distinct known lineages -> ``"Convergent"`` (independent
      evolution, strong evidence for true epistasis).
    * Exactly one known lineage -> ``"Clonal"`` (shared descent, requires
      further scrutiny).
    * No lineage information at all -> ``"Unknown"``.
    """
    if not lineage_map:
        return CLONALITY_UNKNOWN
    dist = lineage_distribution(accessions, lineage_map)
    known = {
        lineage: count for lineage, count in dist.items() if lineage != _UNASSIGNED and count > 0
    }
    if not known:
        return CLONALITY_UNKNOWN
    if len(known) >= 2:
        return CLONALITY_CONVERGENT
    return CLONALITY_CLONAL


def annotate_genome_map(genome_map: pd.DataFrame, lineage_map: Mapping[str, str]) -> pd.DataFrame:
    """Add a ``Lineage`` column to the genome mutation map."""
    if genome_map.empty or "Accession" not in genome_map.columns:
        return genome_map
    if not lineage_map:
        return genome_map
    out = genome_map.copy()
    out["Lineage"] = [lineage_map.get(str(acc), _UNASSIGNED) for acc in out["Accession"]]
    return out


def lineage_distribution_str(accessions: Sequence[str], lineage_map: Mapping[str, str]) -> str:
    """Render lineage distribution as a compact ``ST3:5,ST11:3`` string."""
    if not lineage_map:
        return ""
    dist = lineage_distribution(accessions, lineage_map)
    parts = [f"{lineage}:{count}" for lineage, count in sorted(dist.items())]
    return ",".join(parts)


def _variant_carriers(mutations: pd.DataFrame, variant_label: str) -> set[str]:
    """Accessions carrying a specific ``Gene:Mutation`` variant."""
    if mutations.empty:
        return set()
    gene, _, mutation = variant_label.partition(":")
    mask = (mutations["Gene"].astype(str) == gene) & (mutations["Mutation"].astype(str) == mutation)
    return set(mutations.loc[mask, "Accession"].astype(str))


def annotate_epistasis(
    networks: pd.DataFrame,
    mutations: pd.DataFrame,
    lineage_map: Mapping[str, str],
) -> pd.DataFrame:
    """Add lineage distribution and clonality columns to epistasis networks.

    Returns the frame unchanged when *lineage_map* is empty.
    """
    if networks.empty or not lineage_map:
        return networks

    distributions: list[str] = []
    clonalities: list[str] = []
    for row in networks.itertuples(index=False):
        both = _variant_carriers(mutations, row.Variant_A) & _variant_carriers(
            mutations, row.Variant_B
        )
        dist_str = lineage_distribution_str(sorted(both), lineage_map)
        distributions.append(dist_str)
        clonalities.append(classify_clonality(sorted(both), lineage_map))

    out = networks.copy()
    out["Lineage_Distribution"] = distributions
    out["Clonality"] = clonalities
    return out
