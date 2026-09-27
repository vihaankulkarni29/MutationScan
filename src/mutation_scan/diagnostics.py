"""Explaining an empty result.

A run can complete without error and still produce nothing: every alignment
missed, or every pair failed a QC gate. The tables are then all valid and all
empty, which looks identical to "this cohort genuinely has no mutations". These
functions tell the two apart and say which one happened.

Every message is built from the stage's own output -- the dominant status, the
observed identity range, the actual ``tblastn`` stderr -- rather than from a
fixed script. A canned "something went wrong" would be no better than the silent
empty run it replaces.

Pure functions over DataFrames: no I/O, no logging, no exit codes. The caller
decides what to do with the lines; the tests can assert on them directly.
"""

from __future__ import annotations

import pandas as pd

from .config import Config
from .extract import (
    STATUS_ERROR,
    STATUS_EXTRACTED,
    STATUS_NO_HIT,
    STATUS_NO_REFERENCE,
)
from .variants import (
    VARIANT_STATUS_CALLED,
    VARIANT_STATUS_ERROR,
    VARIANT_STATUS_LOW_COVERAGE,
    VARIANT_STATUS_LOW_IDENTITY,
    VARIANT_STATUS_MISSING_PROTEIN,
    VARIANT_STATUS_NO_REFERENCE,
)

__all__ = [
    "diagnose_empty_extraction",
    "diagnose_empty_variants",
    "dominant_status",
    "format_causes",
]


def dominant_status(frame: pd.DataFrame, exclude: str = "") -> tuple[str, int]:
    """The most common value in ``Status``, and how many rows carry it.

    Returns ``("", 0)`` when there is nothing to count.
    """
    if frame is None or frame.empty or "Status" not in frame.columns:
        return "", 0
    counts = frame["Status"].value_counts()
    if exclude:
        counts = counts.drop(labels=[exclude], errors="ignore")
    if counts.empty:
        return "", 0
    return str(counts.index[0]), int(counts.iloc[0])


def format_causes(causes: list[str]) -> list[str]:
    """Number a cause list, most likely first."""
    if not causes:
        return []
    return ["Likely causes, most common first:"] + [
        f"    {index}. {cause}" for index, cause in enumerate(causes, start=1)
    ]


def _range_of(frame: pd.DataFrame, column: str) -> tuple[float, float] | None:
    """``(min, max)`` of a numeric column, ignoring blanks. ``None`` if empty."""
    if column not in frame.columns:
        return None
    values = pd.to_numeric(frame[column], errors="coerce").dropna()
    if values.empty:
        return None
    return float(values.min()), float(values.max())


def _sample_details(frame: pd.DataFrame, status: str, limit: int = 2) -> list[str]:
    """Distinct ``Detail`` strings for rows with *status* -- the real error text."""
    if "Detail" not in frame.columns:
        return []
    subset = frame[frame["Status"] == status]
    details = [
        str(detail).strip() for detail in subset["Detail"].dropna().unique() if str(detail).strip()
    ]
    return details[:limit]


def diagnose_empty_extraction(manifest: pd.DataFrame, config: Config) -> list[str]:
    """Explain why extraction produced no proteins. One list entry per line."""
    attempted = 0 if manifest is None else len(manifest)

    if attempted == 0:
        return [
            "Extraction attempted 0 (genome, gene) pairs.",
            f"No work was scheduled. Checked {config.genomes_dir} for assemblies "
            f"and {config.references_dir} for reference proteins.",
        ]

    genomes = int(manifest["Accession"].nunique()) if "Accession" in manifest else 0
    genes = int(manifest["Gene"].nunique()) if "Gene" in manifest else 0
    status, count = dominant_status(manifest, exclude=STATUS_EXTRACTED)

    lines = [f"Extraction produced 0 proteins from {attempted} attempts."]

    if status == STATUS_NO_HIT:
        lines.append(f"All {count} attempts returned no_hit.")
        lines.extend(
            format_causes(
                [
                    "--genomes contains protein FASTA, not nucleotide assemblies",
                    "the references are from too distant an organism",
                    "--references and --genomes are swapped",
                    "the assemblies are heavily fragmented or largely masked",
                ]
            )
        )
    elif status == STATUS_ERROR:
        lines.append(f"All {count} attempts failed with an error.")
        for detail in _sample_details(manifest, STATUS_ERROR):
            lines.append(f"    tblastn said: {detail}")
        lines.extend(
            format_causes(
                [
                    "the BLAST+ installation is broken or a wrong binary is configured",
                    "a genome file is truncated or not valid FASTA",
                    "the per-alignment timeout is too low for these assemblies",
                ]
            )
        )
    elif status == STATUS_NO_REFERENCE:
        lines.append(f"{count} attempt(s) had no usable reference to align.")
        lines.append(
            f"References were prepared into {config.refs_dir}; check that every "
            "target resolved to a readable protein FASTA."
        )
    elif status:
        lines.append(f"Dominant status was '{status}' on {count} attempt(s).")

    lines.append(
        f"Checked: {config.genomes_dir} ({genomes} genome(s)), "
        f"{config.references_dir} ({genes} reference(s))."
    )
    return lines


def diagnose_empty_variants(qc: pd.DataFrame, config: Config) -> list[str]:
    """Explain why variant calling evaluated nothing. One list entry per line."""
    considered = 0 if qc is None else len(qc)

    if considered == 0:
        return [
            "Variant calling evaluated 0 (genome, gene) pairs.",
            "Extraction reported proteins but none reached the caller. This is an "
            f"internal inconsistency; inspect {config.manifest_csv}.",
        ]

    status, count = dominant_status(qc, exclude=VARIANT_STATUS_CALLED)
    lines = [f"Variant calling passed 0 of {considered} pair(s) through the QC gates."]

    if status == VARIANT_STATUS_LOW_IDENTITY:
        observed = _range_of(qc[qc["Status"] == status], "Identity_pct")
        measured = f" scoring {observed[0]:.1f}-{observed[1]:.1f}% identity" if observed else ""
        lines.append(
            f"{count} pair(s) were rejected by the identity gate{measured}, against "
            f"a threshold of {config.min_identity_percent:.1f}%."
        )
        lines.extend(
            format_causes(
                [
                    "the references are from a different species than the cohort",
                    f"--min-identity {config.min_identity_percent:.0f} is too strict "
                    "for this level of divergence",
                    "tblastn recovered a paralogue rather than the intended gene",
                ]
            )
        )
    elif status == VARIANT_STATUS_LOW_COVERAGE:
        observed = _range_of(qc[qc["Status"] == status], "Coverage_pct")
        measured = (
            f" covering {observed[0]:.1f}-{observed[1]:.1f}% of the reference" if observed else ""
        )
        lines.append(
            f"{count} pair(s) were rejected by the coverage gate{measured}, against "
            f"a threshold of {config.min_coverage_percent:.1f}%."
        )
        lines.extend(
            format_causes(
                [
                    "the gene is fragmented across contigs in these assemblies",
                    "the assemblies are drafts with the gene only partially present",
                    f"--min-coverage {config.min_coverage_percent:.0f} is too strict "
                    "for draft-quality data",
                ]
            )
        )
    elif status == VARIANT_STATUS_MISSING_PROTEIN:
        lines.append(
            f"{count} extracted protein(s) could not be found on disk. Extraction "
            f"and variant calling disagree; inspect {config.proteins_dir}."
        )
    elif status == VARIANT_STATUS_NO_REFERENCE:
        lines.append(f"{count} pair(s) had no usable reference in {config.refs_dir}.")
    elif status == VARIANT_STATUS_ERROR:
        lines.append(f"{count} pair(s) failed with an error.")
        for detail in _sample_details(qc, VARIANT_STATUS_ERROR):
            lines.append(f"    {detail}")
    elif status:
        lines.append(f"Dominant status was '{status}' on {count} pair(s).")

    lines.append(
        f"Gates in force: identity >= {config.min_identity_percent:.1f}%, "
        f"coverage >= {config.min_coverage_percent:.1f}%. Full detail in "
        f"{config.variant_qc_csv}."
    )
    return lines
