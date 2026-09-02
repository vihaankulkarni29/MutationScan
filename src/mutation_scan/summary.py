"""Run provenance: what was run, on what, with which settings, and what came out.

``run_summary.json`` is written by both entry points (the CLI and the Snakemake
workflow) through this one function, so a run is described identically however it
was launched. It is the file to attach to a manuscript or hand to a collaborator
who asks "what exactly did you do?".

It records counts and settings only. Like every other output in this package it
makes no claim about what the numbers mean.
"""

from __future__ import annotations

import json
import logging
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from . import __version__
from .config import Config
from .cooccurrence import CooccurrenceResult
from .extract import STATUS_EXTRACTED
from .references import ReferenceSet
from .variants import VARIANT_STATUS_CALLED

logger = logging.getLogger(__name__)

__all__ = ["build_run_summary", "write_run_summary"]


def _status_counts(frame: Optional[pd.DataFrame]) -> dict[str, int]:
    """``{status: n}`` for any table carrying a Status column."""
    if frame is None or frame.empty or "Status" not in frame.columns:
        return {}
    counts = frame["Status"].value_counts()
    return {str(status): int(count) for status, count in counts.items()}


def _n_unique(frame: Optional[pd.DataFrame], column: str) -> int:
    if frame is None or frame.empty or column not in frame.columns:
        return 0
    return int(frame[column].nunique())


def build_run_summary(
    config: Config,
    references: Optional[ReferenceSet] = None,
    manifest: Optional[pd.DataFrame] = None,
    qc: Optional[pd.DataFrame] = None,
    mutations: Optional[pd.DataFrame] = None,
    cooccurrence: Optional[CooccurrenceResult] = None,
    tblastn_version: str = "",
) -> dict[str, Any]:
    """Assemble the provenance record for one run.

    Every argument beyond *config* is optional so a partial run (or a dry run)
    still produces a readable summary rather than failing.
    """
    extracted = (
        manifest[manifest["Status"] == STATUS_EXTRACTED]
        if manifest is not None and not manifest.empty and "Status" in manifest.columns
        else None
    )
    called = (
        qc[qc["Status"] == VARIANT_STATUS_CALLED]
        if qc is not None and not qc.empty and "Status" in qc.columns
        else None
    )

    summary: dict[str, Any] = {
        "mutation_scan_version": __version__,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "tblastn": tblastn_version,
        },
        "config": config.to_dict(),
        "references": {
            "genes": list(references.genes) if references else list(config.targets),
            "missing": list(references.missing) if references else [],
            "sources": dict(references.sources) if references else {},
            "lengths": dict(references.lengths) if references else {},
        },
        "extraction": {
            "genomes": _n_unique(manifest, "Accession"),
            "pairs_attempted": int(len(manifest)) if manifest is not None else 0,
            "status_counts": _status_counts(manifest),
            "genomes_with_any_protein": _n_unique(extracted, "Accession"),
        },
        "variants": {
            "pairs_evaluated": int(len(qc)) if qc is not None else 0,
            "status_counts": _status_counts(qc),
            "pairs_called": int(len(called)) if called is not None else 0,
            "mutations": int(len(mutations)) if mutations is not None else 0,
            "genomes_with_any_mutation": _n_unique(mutations, "Accession"),
            "gates": {
                "min_identity_percent": config.min_identity_percent,
                "min_coverage_percent": config.min_coverage_percent,
            },
        },
        "cooccurrence": {
            "levels": list(config.cooccurrence.levels),
            "min_count": config.cooccurrence.min_count,
            "gene_pairs_reported": int(len(cooccurrence.genes)) if cooccurrence else 0,
            "variant_pairs_reported": (
                int(len(cooccurrence.variants))
                if cooccurrence is not None and cooccurrence.variants is not None
                else 0
            ),
            # Per gene: mutated in N of M genomes where the gene was evaluated.
            # These are the numerators and denominators behind every pair count.
            "mutated_per_gene": dict(cooccurrence.gene_totals) if cooccurrence else {},
            "eligible_per_gene": dict(cooccurrence.eligible_per_gene) if cooccurrence else {},
        },
        "outputs": [str(path) for path in config.outputs()],
        "note": (
            "Co-occurrence counts are descriptive. They are not tests of "
            "interaction, selection or causation."
        ),
    }
    return summary


def write_run_summary(
    config: Config,
    path: Path | str | None = None,
    **kwargs: Any,
) -> dict[str, Any]:
    """Build the summary and write it as JSON. Returns the summary dict."""
    summary = build_run_summary(config, **kwargs)
    target = Path(path) if path is not None else config.run_summary_json
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    logger.info("Wrote run summary: %s", target)
    return summary
