"""Protein extraction from genomic DNA with tblastn.

``tblastn`` is a translating aligner: it aligns a protein query against all six
reading frames of nucleotide subject sequences and reports the *translated*
subject segment. Taking that segment (``sseq``) rather than translating a
predicted ORF ourselves is what keeps frameshift artefacts out of the extracted
proteins.

Every (genome, gene) pair attempted produces exactly one manifest row, whether
it succeeded or not. That manifest is the run's record of what was actually
observed, and downstream stages read it instead of guessing from filenames --
so "this gene is not mutated here" is never confused with "this gene was never
seen here".
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence

import pandas as pd

from .fasta import sanitize_protein_sequence, write_protein_record
from .references import ReferenceSet, resolve_reference

logger = logging.getLogger(__name__)

__all__ = [
    "DEFAULT_EVALUE",
    "GENOME_EXTENSIONS",
    "MANIFEST_COLUMNS",
    "STATUS_ERROR",
    "STATUS_EXTRACTED",
    "STATUS_NO_HIT",
    "STATUS_NO_REFERENCE",
    "TblastnError",
    "discover_genomes",
    "extract",
    "protein_path",
    "read_manifest",
    "verify_tblastn",
]

#: Assembly file extensions treated as genomes. One file per genome.
GENOME_EXTENSIONS = frozenset({".fna", ".fa", ".fasta", ".fas", ".fsa"})

#: Standard homology cutoff. Loose enough for diverged orthologues, tight
#: enough to reject chance alignments; identity/coverage gating happens later.
DEFAULT_EVALUE = 1e-5

#: Per-alignment wall-clock ceiling, so one pathological genome cannot stall a run.
TBLASTN_TIMEOUT = 180

STATUS_EXTRACTED = "extracted"
STATUS_NO_HIT = "no_hit"
STATUS_NO_REFERENCE = "no_reference"
STATUS_ERROR = "error"

MANIFEST_COLUMNS = (
    "Accession",
    "Gene",
    "Status",
    "Pident",
    "Aln_Length",
    "Ref_Length",
    "Ref_Coverage_pct",
    "Extracted_Length",
    "Internal_Stops",
    "Contig",
    "Detail",
)

# tblastn tabular fields we request, in order. sseq last so that a stray tab in
# an earlier field cannot silently shift the sequence column.
_OUTFMT_FIELDS = (
    "qseqid sseqid pident length mismatch gapopen qstart qend sstart send "
    "evalue bitscore qlen slen sseq"
)
_FIELD_COUNT = len(_OUTFMT_FIELDS.split())


class TblastnError(RuntimeError):
    """tblastn could not be run, or failed on a specific alignment."""


def verify_tblastn(binary: str = "tblastn") -> str:
    """Confirm the tblastn binary works. Returns its version banner.

    Raises:
        TblastnError: the binary is absent, not executable, or exits non-zero.
    """
    if shutil.which(binary) is None and not Path(binary).is_file():
        raise TblastnError(
            f"tblastn not found: {binary!r}. Install BLAST+ and ensure tblastn is on "
            "PATH, or set 'tblastn_binary' in the config to its full path."
        )

    try:
        result = subprocess.run(
            [binary, "-version"], capture_output=True, text=True, timeout=30
        )
    except OSError as exc:
        raise TblastnError(f"Could not execute {binary!r}: {exc}") from exc
    except subprocess.TimeoutExpired as exc:
        raise TblastnError(f"{binary} -version timed out") from exc

    if result.returncode != 0:
        # The old implementation fell through this case and reported success.
        raise TblastnError(
            f"{binary} -version exited {result.returncode}: "
            f"{(result.stderr or result.stdout or '').strip()}"
        )

    banner = (result.stdout or "").strip().splitlines()
    version = banner[0] if banner else "unknown version"
    logger.info("tblastn verified: %s", version)
    return version


def discover_genomes(genomes_dir: Path | str) -> list[Path]:
    """Genome assemblies in *genomes_dir*, sorted for reproducible runs."""
    genomes_dir = Path(genomes_dir)
    if not genomes_dir.is_dir():
        return []
    return sorted(
        path
        for path in genomes_dir.iterdir()
        if path.is_file()
        and path.suffix.lower() in GENOME_EXTENSIONS
        and path.stat().st_size > 0
    )


def protein_path(proteins_dir: Path | str, accession: str, gene: str) -> Path:
    """Where the extracted protein for one (genome, gene) pair is written."""
    return Path(proteins_dir) / f"{accession}_{gene}.faa"


@dataclass(frozen=True)
class TblastnHit:
    """The best HSP of one reference protein against one genome."""

    sequence: str
    pident: float
    aln_length: int
    ref_length: int
    ref_coverage_pct: float
    contig: str
    internal_stops: int


def _parse_hit(line: str) -> Optional[TblastnHit]:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < _FIELD_COUNT:
        raise TblastnError(
            f"Unexpected tblastn output: {len(fields)} fields, expected {_FIELD_COUNT}"
        )

    contig = fields[1]
    pident = float(fields[2])
    aln_length = int(fields[3])
    qstart, qend = int(fields[6]), int(fields[7])
    ref_length = int(fields[12])
    aligned_seq = fields[14]

    sequence = sanitize_protein_sequence(aligned_seq.replace("-", ""))
    if not sequence:
        return None

    covered = max(0, qend - qstart + 1)
    coverage = (covered / ref_length * 100.0) if ref_length else 0.0

    return TblastnHit(
        sequence=sequence,
        pident=pident,
        aln_length=aln_length,
        ref_length=ref_length,
        ref_coverage_pct=round(coverage, 2),
        contig=contig,
        internal_stops=sequence.count("*"),
    )


def run_tblastn(
    ref_faa: Path,
    genome_fna: Path,
    binary: str = "tblastn",
    evalue: float = DEFAULT_EVALUE,
    timeout: int = TBLASTN_TIMEOUT,
) -> Optional[TblastnHit]:
    """Align one reference protein against one genome; return the best HSP.

    Returns ``None`` when there is no hit above *evalue*.

    Raises:
        TblastnError: tblastn failed, timed out, or produced unreadable output.
    """
    command = [
        binary,
        "-query", str(ref_faa),
        "-subject", str(genome_fna),
        "-outfmt", f"6 {_OUTFMT_FIELDS}",
        "-evalue", repr(evalue),
        "-max_target_seqs", "1",
        "-max_hsps", "1",
    ]

    try:
        result = subprocess.run(
            command, capture_output=True, text=True, timeout=timeout, check=True
        )
    except subprocess.TimeoutExpired as exc:
        raise TblastnError(f"tblastn timed out after {timeout}s") from exc
    except subprocess.CalledProcessError as exc:
        raise TblastnError(
            f"tblastn exited {exc.returncode}: {(exc.stderr or '').strip()}"
        ) from exc
    except OSError as exc:
        raise TblastnError(f"Could not execute tblastn: {exc}") from exc

    output = (result.stdout or "").strip()
    if not output:
        return None
    return _parse_hit(output.splitlines()[0])


@dataclass(frozen=True)
class _Task:
    accession: str
    genome: Path
    gene: str
    reference: Optional[Path]


def _row(
    accession: str,
    gene: str,
    status: str,
    detail: str = "",
    hit: Optional[TblastnHit] = None,
) -> dict:
    return {
        "Accession": accession,
        "Gene": gene,
        "Status": status,
        "Pident": hit.pident if hit else pd.NA,
        "Aln_Length": hit.aln_length if hit else pd.NA,
        "Ref_Length": hit.ref_length if hit else pd.NA,
        "Ref_Coverage_pct": hit.ref_coverage_pct if hit else pd.NA,
        "Extracted_Length": len(hit.sequence) if hit else pd.NA,
        "Internal_Stops": hit.internal_stops if hit else pd.NA,
        "Contig": hit.contig if hit else "",
        "Detail": detail,
    }


def _run_task(task: _Task, proteins_dir: Path, binary: str, evalue: float) -> dict:
    if task.reference is None:
        return _row(task.accession, task.gene, STATUS_NO_REFERENCE, "no reference available")

    try:
        hit = run_tblastn(task.reference, task.genome, binary=binary, evalue=evalue)
    except TblastnError as exc:
        logger.warning("%s/%s: %s", task.accession, task.gene, exc)
        return _row(task.accession, task.gene, STATUS_ERROR, str(exc))

    if hit is None:
        return _row(task.accession, task.gene, STATUS_NO_HIT, "no alignment above e-value")

    try:
        write_protein_record(
            protein_path(proteins_dir, task.accession, task.gene),
            hit.sequence,
            seq_id=f"{task.accession}_{task.gene}",
            description=(
                f"gene={task.gene} accession={task.accession} "
                f"contig={hit.contig} pident={hit.pident}"
            ),
        )
    except (OSError, ValueError) as exc:
        logger.error("%s/%s: could not write protein: %s", task.accession, task.gene, exc)
        return _row(task.accession, task.gene, STATUS_ERROR, f"write failed: {exc}", hit)

    if hit.internal_stops:
        logger.info(
            "%s/%s: %d internal stop codon(s) in the extracted protein",
            task.accession,
            task.gene,
            hit.internal_stops,
        )

    return _row(task.accession, task.gene, STATUS_EXTRACTED, "", hit)


def extract(
    genomes_dir: Path | str,
    refs_dir: Path | str,
    proteins_out: Path | str,
    targets: Sequence[str] | ReferenceSet,
    manifest_csv: Path | str | None = None,
    tblastn_binary: str = "tblastn",
    threads: int = 1,
    evalue: float = DEFAULT_EVALUE,
    genomes: Iterable[Path] | None = None,
) -> pd.DataFrame:
    """Extract every target protein from every genome.

    Args:
        genomes_dir: directory of nucleotide assemblies (one file per genome).
        refs_dir: directory of prepared references (``{gene}_WT.faa``).
        proteins_out: directory to write ``{accession}_{gene}.faa`` into.
        targets: gene names, or the :class:`~mutation_scan.references.ReferenceSet`
            returned by ``prepare_references``.
        manifest_csv: optional path to write the manifest to.
        threads: number of concurrent tblastn processes.
        genomes: explicit genome file list, bypassing discovery (for tests).

    Returns:
        The manifest: one row per attempted (genome, gene) pair.
    """
    refs_dir = Path(refs_dir)
    proteins_out = Path(proteins_out)
    proteins_out.mkdir(parents=True, exist_ok=True)

    gene_names = list(targets.genes) if isinstance(targets, ReferenceSet) else [
        str(gene).strip() for gene in targets if str(gene).strip()
    ]

    genome_files = list(genomes) if genomes is not None else discover_genomes(genomes_dir)

    if not genome_files:
        logger.warning("No genome assemblies found in %s", genomes_dir)
    if not gene_names:
        logger.warning("No target genes to extract")

    # Resolve each reference once, not once per genome.
    references = {gene: resolve_reference(refs_dir, gene) for gene in gene_names}
    for gene, path in references.items():
        if path is None:
            logger.error("No prepared reference for %s in %s", gene, refs_dir)

    tasks = [
        _Task(accession=genome.stem, genome=genome, gene=gene, reference=references[gene])
        for genome in genome_files
        for gene in gene_names
    ]

    if tasks:
        verify_tblastn(tblastn_binary)
        logger.info(
            "Extracting %d gene(s) from %d genome(s) = %d alignments on %d thread(s)",
            len(gene_names),
            len(genome_files),
            len(tasks),
            threads,
        )

    rows: list[dict] = []
    workers = max(1, int(threads))
    step = max(1, len(tasks) // 20)

    if workers == 1:
        for index, task in enumerate(tasks, start=1):
            rows.append(_run_task(task, proteins_out, tblastn_binary, evalue))
            if index % step == 0 or index == len(tasks):
                logger.info("  extraction %d/%d", index, len(tasks))
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [
                pool.submit(_run_task, task, proteins_out, tblastn_binary, evalue)
                for task in tasks
            ]
            for index, future in enumerate(futures, start=1):
                rows.append(future.result())
                if index % step == 0 or index == len(tasks):
                    logger.info("  extraction %d/%d", index, len(tasks))

    manifest = pd.DataFrame(rows, columns=list(MANIFEST_COLUMNS))
    if not manifest.empty:
        manifest = manifest.sort_values(["Accession", "Gene"], kind="stable").reset_index(drop=True)

    if manifest_csv is not None:
        manifest_csv = Path(manifest_csv)
        manifest_csv.parent.mkdir(parents=True, exist_ok=True)
        manifest.to_csv(manifest_csv, index=False)
        logger.info("Wrote extraction manifest: %s", manifest_csv)

    _log_summary(manifest, gene_names)
    return manifest


def _log_summary(manifest: pd.DataFrame, gene_names: Sequence[str]) -> None:
    if manifest.empty:
        logger.warning("Extraction produced no results")
        return

    counts = manifest["Status"].value_counts()
    logger.info(
        "Extraction complete: %s",
        ", ".join(f"{status}={count}" for status, count in counts.items()),
    )

    genomes = manifest["Accession"].nunique()
    for gene in gene_names:
        subset = manifest[manifest["Gene"] == gene]
        if subset.empty:
            continue
        extracted = int((subset["Status"] == STATUS_EXTRACTED).sum())
        logger.info(
            "  %s: extracted in %d/%d genomes (%.1f%%)",
            gene,
            extracted,
            genomes,
            extracted / genomes * 100 if genomes else 0.0,
        )


def read_manifest(manifest_csv: Path | str) -> pd.DataFrame:
    """Read an extraction manifest, keeping accessions and genes as strings."""
    manifest = pd.read_csv(manifest_csv, dtype={"Accession": str, "Gene": str, "Status": str})
    missing = [column for column in ("Accession", "Gene", "Status") if column not in manifest]
    if missing:
        raise ValueError(f"Manifest {manifest_csv} is missing column(s): {', '.join(missing)}")
    return manifest
