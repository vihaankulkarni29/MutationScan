"""Input validation that runs before any alignment work.

Every check here exists because its absence produces a *quiet* failure: a run
that completes, exits zero, and writes six empty tables. The classic case is
putting protein FASTA in the genomes directory (or swapping the two directories
outright), which costs a full ``tblastn`` sweep and yields nothing but
``no_hit`` rows, with no output saying why.

Two principles:

* **Cheap.** Only the first record of each file is read, capped at a few
  thousand residues, so validating a thousand 5 Mb assemblies takes about as
  long as one alignment.
* **Complete.** Every problem is collected and reported together. A user who
  has both swapped their directories *and* left a corrupt reference on disk
  should learn both facts from one run, not from two.

Nothing here writes to disk, so ``config-check`` can call it as a true dry run.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path

from .config import Config
from .extract import GENOME_EXTENSIONS, TblastnError, discover_genomes, verify_tblastn
from .fasta import (
    looks_like_nucleotide,
    read_first_record_sequence,
    read_protein_sequence,
)
from .references import REFERENCE_EXTENSIONS, discover_targets, resolve_reference

logger = logging.getLogger(__name__)

__all__ = ["PreflightReport", "preflight"]

#: How many offending filenames to name before summarising the rest. Enough to
#: recognise the pattern, few enough that a 500-genome mistake stays readable.
_MAX_NAMED = 5


def _summarise(paths: list[Path], limit: int = _MAX_NAMED) -> str:
    """``"a.faa, b.faa and 3 more"`` -- bounded regardless of list length."""
    names = [path.name for path in paths[:limit]]
    remainder = len(paths) - len(names)
    listed = ", ".join(names)
    return f"{listed} and {remainder} more" if remainder else listed


@dataclass
class PreflightReport:
    """What validation found. Empty ``errors`` means the run may proceed."""

    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    n_genomes: int = 0
    n_references: int = 0
    targets: list[str] = field(default_factory=list)
    tblastn_version: str = ""

    @property
    def ok(self) -> bool:
        return not self.errors

    def log(self) -> None:
        """Emit every finding at its proper level."""
        for warning in self.warnings:
            logger.warning("%s", warning)
        for error in self.errors:
            logger.error("%s", error)


def _check_tblastn(config: Config, report: PreflightReport) -> None:
    try:
        report.tblastn_version = verify_tblastn(config.tblastn_binary)
    except TblastnError as exc:
        report.errors.append(str(exc))


def _check_genomes(config: Config, report: PreflightReport) -> None:
    """Genomes must exist, be readable, and actually be nucleotide sequences."""
    genomes = discover_genomes(config.genomes_dir)
    report.n_genomes = len(genomes)

    if not genomes:
        present = sorted(
            {path.suffix.lower() for path in config.genomes_dir.iterdir() if path.is_file()}
            if config.genomes_dir.is_dir()
            else set()
        )
        detail = (
            f" The directory holds {', '.join(present)} files."
            if present
            else " The directory is empty."
        )
        report.errors.append(
            f"No genome assemblies in {config.genomes_dir}.{detail} Expected one "
            f"nucleotide FASTA per genome with an extension of "
            f"{', '.join(sorted(GENOME_EXTENSIONS))}."
        )
        return

    protein_like: list[Path] = []
    unreadable: list[Path] = []
    for path in genomes:
        sequence = read_first_record_sequence(path)
        if not sequence:
            unreadable.append(path)
        elif not looks_like_nucleotide(sequence):
            protein_like.append(path)

    if protein_like:
        report.errors.append(
            f"{len(protein_like)} of {len(genomes)} file(s) in {config.genomes_dir} "
            f"look like protein, not nucleotide assemblies: {_summarise(protein_like)}. "
            "Are --genomes and --references swapped?"
        )
    if unreadable:
        report.errors.append(
            f"{len(unreadable)} genome file(s) in {config.genomes_dir} contain no "
            f"readable sequence: {_summarise(unreadable)}."
        )


def _check_references(config: Config, report: PreflightReport) -> None:
    """References define the targets, so an unusable one silently drops a gene."""
    targets = discover_targets(config.references_dir, config.targets)
    report.targets = list(targets)

    if not targets:
        report.errors.append(
            f"No reference proteins in {config.references_dir}. Add one FASTA per "
            f"target protein, named <gene>_WT.faa (or <gene>.faa); each filename "
            f"becomes a target. Accepted extensions: "
            f"{', '.join(sorted(REFERENCE_EXTENSIONS))}."
        )
        return

    missing: list[str] = []
    unparsable: list[str] = []
    nucleotide_like: list[Path] = []

    for gene in targets:
        path = resolve_reference(config.references_dir, gene)
        if path is None:
            # Only reachable via explicit --targets naming a gene with no file.
            missing.append(gene)
            continue
        try:
            record = read_protein_sequence(path, fallback_id=gene)
        except (ValueError, OSError) as exc:
            unparsable.append(f"{path.name} ({exc})")
            continue
        if looks_like_nucleotide(str(record.seq)):
            nucleotide_like.append(path)

    report.n_references = len(targets) - len(missing) - len(unparsable)

    if missing:
        report.errors.append(
            f"No reference file for target(s): {', '.join(missing)}. Supply "
            f"<gene>_WT.faa in {config.references_dir}, or set uniprot_taxid to "
            "fetch from UniProt."
        )
    if unparsable:
        # An error, not a warning: a dropped reference silently removes a gene
        # from every downstream table, including the co-occurrence denominators.
        report.errors.append(
            f"{len(unparsable)} reference file(s) in {config.references_dir} could "
            f"not be read: {'; '.join(unparsable[:_MAX_NAMED])}."
        )
    if nucleotide_like:
        report.errors.append(
            f"{len(nucleotide_like)} reference file(s) in {config.references_dir} "
            f"look like nucleotide sequences, not proteins: "
            f"{_summarise(nucleotide_like)}. References must be amino-acid FASTA. "
            "Are --genomes and --references swapped?"
        )


def _check_output(config: Config, report: PreflightReport) -> None:
    """Writability of the nearest existing ancestor -- no directory is created."""
    ancestor = config.run_dir
    while not ancestor.exists() and ancestor != ancestor.parent:
        ancestor = ancestor.parent

    if not ancestor.is_dir():
        report.errors.append(f"Output path is not a directory: {ancestor}")
    elif not os.access(ancestor, os.W_OK):
        report.errors.append(
            f"Output directory is not writable: {ancestor} (results would go to {config.run_dir})"
        )


def preflight(config: Config, check_tblastn: bool = True) -> PreflightReport:
    """Validate a run's inputs before any alignment is attempted.

    Args:
        config: the resolved run configuration.
        check_tblastn: verify the BLAST+ binary. Disabled by tests that have no
            BLAST+ installed but still want the data checks.

    Returns:
        A :class:`PreflightReport`. Every problem found is collected, so the
        caller can show them all at once; ``report.ok`` is the go/no-go.
    """
    report = PreflightReport()

    if check_tblastn:
        _check_tblastn(config, report)
    _check_genomes(config, report)
    _check_references(config, report)
    _check_output(config, report)

    if report.ok:
        logger.info(
            "Preflight passed: %d genome(s), %d reference(s), targets: %s",
            report.n_genomes,
            report.n_references,
            ", ".join(report.targets),
        )
    return report
