"""Biochemical severity scoring for amino-acid substitutions.

No thermodynamics, no docking, no PDB parsing, no Grantham scoring.
Scores are derived from:

* BLOSUM62 evolutionary substitution matrix
* Charge change at pH 7
* Kyte-Doolittle hydropathy shift
* Amino-acid volume class shift

The composite severity is a transparent weighted combination on a 0-100 scale.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping

import pandas as pd
from Bio.Align import substitution_matrices

from .variants import LOF_ALLELE

logger = logging.getLogger(__name__)

__all__ = [
    "CHARGE",
    "HYDROPATHY",
    "VOLUME",
    "LOF_SEVERITY",
    "compute_severity",
    "mutation_severity",
    "pair_severity",
]

_BLOSUM62 = substitution_matrices.load("BLOSUM62")

#: Loss-of-function is a categorical class, not a substitution, so it carries
#: no per-residue biochemical delta. It scores at the top of the 0-100 scale:
#: a knockout or frameshift destroys (or near-destroys) the protein's
#: structure, which is the maximal functional impact the severity model can
#: express. A pair's severity is the mean of its two variants' severity
#: scores, so a LoF paired with anything inherits most of that impact.
LOF_SEVERITY = 100.0

CHARGE: Mapping[str, int] = {
    "A": 0,
    "C": 0,
    "D": -1,
    "E": -1,
    "F": 0,
    "G": 0,
    "H": 1,
    "I": 0,
    "K": 1,
    "L": 0,
    "M": 0,
    "N": 0,
    "P": 0,
    "Q": 0,
    "R": 1,
    "S": 0,
    "T": 0,
    "V": 0,
    "W": 0,
    "Y": 0,
}

HYDROPATHY: Mapping[str, float] = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}

VOLUME: Mapping[str, int] = {
    "A": 1,
    "C": 2,
    "D": 2,
    "E": 3,
    "F": 4,
    "G": 0,
    "H": 3,
    "I": 4,
    "K": 4,
    "L": 4,
    "M": 3,
    "N": 2,
    "P": 3,
    "Q": 3,
    "R": 5,
    "S": 1,
    "T": 2,
    "V": 3,
    "W": 4,
    "Y": 4,
}


def _safe_lookup(
    table: Mapping[str, float | int], key: str, default: float | int = 0
) -> float | int:
    return table.get(key.upper(), default)


def _blosum62_score(ref: str, alt: str) -> int:
    return int(_BLOSUM62.get(ref.upper(), {}).get(alt.upper(), 0))


def compute_severity(ref_aa: str, alt_aa: str) -> dict:
    """Return biochemical metrics for one substitution."""
    ref = ref_aa.upper()
    alt = alt_aa.upper()

    blosum = _blosum62_score(ref, alt)
    charge_delta = _safe_lookup(CHARGE, alt) - _safe_lookup(CHARGE, ref)
    kd_delta = float(_safe_lookup(HYDROPATHY, alt)) - float(_safe_lookup(HYDROPATHY, ref))
    vol_delta = int(_safe_lookup(VOLUME, alt)) - int(_safe_lookup(VOLUME, ref))

    blosum_norm = max(0.0, min(1.0, (4.0 - blosum) / 8.0))
    charge_norm = min(1.0, abs(charge_delta) / 2.0)
    kd_norm = min(1.0, abs(kd_delta) / 5.0)
    vol_norm = min(1.0, abs(vol_delta) / 3.0)

    severity = 100.0 * (0.40 * blosum_norm + 0.25 * charge_norm + 0.20 * kd_norm + 0.15 * vol_norm)

    return {
        "blosum62_score": blosum,
        "charge_delta": int(charge_delta),
        "kd_delta": round(kd_delta, 2),
        "volume_delta": vol_delta,
        "severity_score": round(float(severity), 2),
    }


def mutation_severity(mutations: pd.DataFrame) -> pd.DataFrame:
    """Attach biochemical severity to every row of a mutations table."""
    if mutations.empty:
        return pd.DataFrame(
            columns=[
                "Accession",
                "Gene",
                "Mutation",
                "Ref_AA",
                "Alt_AA",
                "blosum62_score",
                "charge_delta",
                "kd_delta",
                "volume_delta",
                "severity_score",
            ]
        )

    # compute_severity is a pure function of the two residues, so there are at
    # most a few hundred distinct answers however many mutation rows there are.
    # Memoising it turns ~66k calls at cohort scale into a few hundred, which is
    # the same exactness argument used for the epistasis contingency tables.
    # iterrows is avoided for the same reason: it boxes every row into a Series.
    cache: dict[tuple[str, str], dict] = {}
    accessions = mutations["Accession"].tolist() if "Accession" in mutations.columns else []
    genes = mutations["Gene"].tolist() if "Gene" in mutations.columns else []
    mutations_col = mutations["Mutation"].tolist() if "Mutation" in mutations.columns else []

    records: list[dict] = []
    for accession, gene, raw in zip(accessions, genes, mutations_col, strict=True):
        mut = str(raw)
        if mut == LOF_ALLELE:
            # A loss-of-function allele is not a substitution, so it has no
            # Ref_Alt residues to score with BLOSUM62. Without this guard the
            # "L" and "F" of the literal "LoF" parse as a Leu->Phe change and
            # emit a bogus biochemical score. Severity is categorical.
            records.append(
                {
                    "blosum62_score": 0,
                    "charge_delta": 0,
                    "kd_delta": 0.0,
                    "volume_delta": 0,
                    "severity_score": LOF_SEVERITY,
                    "Accession": accession,
                    "Gene": gene,
                    "Mutation": LOF_ALLELE,
                    "Ref_AA": "",
                    "Alt_AA": "",
                }
            )
            continue
        if len(mut) < 2:
            continue
        ref_aa = mut[0]
        alt_aa = mut[-1]
        key = (ref_aa, alt_aa)
        metrics = cache.get(key)
        if metrics is None:
            metrics = compute_severity(ref_aa, alt_aa)
            cache[key] = metrics
        records.append(
            {
                **metrics,
                "Accession": accession,
                "Gene": gene,
                "Mutation": mut,
                "Ref_AA": ref_aa,
                "Alt_AA": alt_aa,
            }
        )

    frame = pd.DataFrame(records)
    cols = [
        "Accession",
        "Gene",
        "Mutation",
        "Ref_AA",
        "Alt_AA",
        "blosum62_score",
        "charge_delta",
        "kd_delta",
        "volume_delta",
        "severity_score",
    ]
    return frame[cols]


def build_severity_index(severity: pd.DataFrame) -> dict[tuple[str, str], float]:
    """Mean severity per ``(gene, mutation)``, built once for many lookups.

    Pair scoring needs this same lookup once per candidate pair. Computing it
    per pair means re-scanning the whole severity table each time, which is
    quadratic: at cohort scale it dominated the epistasis stage so completely
    that the run never finished.
    """
    if severity.empty:
        return {}
    grouped = (
        severity.assign(
            _gene=severity["Gene"].astype(str),
            _mutation=severity["Mutation"].astype(str),
        )
        .groupby(["_gene", "_mutation"])["severity_score"]
        .mean()
    )
    return {(str(gene), str(mutation)): float(value) for (gene, mutation), value in grouped.items()}


def pair_severity(
    variant_a: str,
    variant_b: str,
    severity: pd.DataFrame,
    index: dict[tuple[str, str], float] | None = None,
) -> float:
    """Average severity of two variants across all genomes that carry them.

    Pass *index* (from :func:`build_severity_index`) when scoring many pairs;
    without it each call rescans *severity*.
    """
    if severity.empty:
        return 0.0
    if index is None:
        index = build_severity_index(severity)

    def _mean(variant: str) -> float:
        gene, mut = variant.split(":", 1)
        return index.get((gene, mut), 0.0)

    return round((_mean(variant_a) + _mean(variant_b)) / 2.0, 2)
