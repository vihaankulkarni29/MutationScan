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
import tempfile
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from functools import partial
from pathlib import Path

import pandas as pd

from .fasta import sanitize_protein_sequence, write_protein_record
from .parallel import WORKLOAD_THREAD, parallel_map
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
    "TBLASTN_THREADS",
    "TblastnError",
    "discover_genomes",
    "extract",
    "protein_path",
    "read_manifest",
    "run_tblastn",
    "run_tblastn_batch",
    "verify_tblastn",
]

#: Assembly file extensions treated as genomes. One file per genome.
GENOME_EXTENSIONS = frozenset({".fna", ".fa", ".fasta", ".fas", ".fsa"})

#: Standard homology cutoff. Loose enough for diverged orthologues, tight
#: enough to reject chance alignments; identity/coverage gating happens later.
DEFAULT_EVALUE = 1e-5

#: Per-alignment wall-clock ceiling, so one pathological genome cannot stall a
#: run. One call now aligns the whole gene panel, so this covers ~15x the work
#: the per-gene call did; a typical genome takes a couple of seconds, so the
#: headroom is still large.
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

#: Threads handed to tblastn itself. Our own pool already provides the
#: parallelism, and tblastn otherwise grabs every core it can find -- so N
#: concurrent tblastn processes each spawn N threads and thrash.
TBLASTN_THREADS = 1


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
        result = subprocess.run([binary, "-version"], capture_output=True, text=True, timeout=30)
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
        if path.is_file() and path.suffix.lower() in GENOME_EXTENSIONS and path.stat().st_size > 0
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


def _parse_row(line: str) -> tuple[str, float, TblastnHit | None]:
    """Parse one tabular row into (query id, bitscore, hit).

    The hit is ``None`` when the aligned segment sanitizes to nothing, which the
    single-query path reports as "no hit" and the batched path must treat the
    same way for that one query.
    """
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
    bitscore = float(fields[11])

    sequence = sanitize_protein_sequence(aligned_seq.replace("-", ""))
    if not sequence:
        return fields[0], bitscore, None

    covered = max(0, qend - qstart + 1)
    coverage = (covered / ref_length * 100.0) if ref_length else 0.0

    return (
        fields[0],
        bitscore,
        TblastnHit(
            sequence=sequence,
            pident=pident,
            aln_length=aln_length,
            ref_length=ref_length,
            ref_coverage_pct=round(coverage, 2),
            contig=contig,
            internal_stops=sequence.count("*"),
        ),
    )


def _parse_hit(line: str) -> TblastnHit | None:
    return _parse_row(line)[2]


def _tblastn_command(
    query: Path,
    subject: Path,
    binary: str,
    evalue: float,
) -> list[str]:
    return [
        binary,
        "-query",
        str(query),
        "-subject",
        str(subject),
        "-outfmt",
        f"6 {_OUTFMT_FIELDS}",
        "-evalue",
        repr(evalue),
        # Per query, not per call: with one multi-query FASTA this still yields
        # at most one target sequence per gene, which is what makes the batched
        # call equivalent to one call per gene.
        "-max_target_seqs",
        "1",
        "-max_hsps",
        "1",
        "-num_threads",
        str(TBLASTN_THREADS),
    ]


def _run_tblastn(command: list[str], timeout: int) -> str:
    """Run tblastn and return stdout, translating every failure mode."""
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
    return result.stdout or ""


def run_tblastn(
    ref_faa: Path,
    genome_fna: Path,
    binary: str = "tblastn",
    evalue: float = DEFAULT_EVALUE,
    timeout: int = TBLASTN_TIMEOUT,
) -> TblastnHit | None:
    """Align one reference protein against one genome; return the best HSP.

    Returns ``None`` when there is no hit above *evalue*.

    Kept as the single-query reference implementation: the batched path is
    checked against it, so it must stay honest rather than merely stay unused.

    Raises:
        TblastnError: tblastn failed, timed out, or produced unreadable output.
    """
    output = _run_tblastn(_tblastn_command(ref_faa, genome_fna, binary, evalue), timeout).strip()
    if not output:
        return None
    return _parse_hit(output.splitlines()[0])


def run_tblastn_batch(
    query_faa: Path,
    genome_fna: Path,
    genes: Sequence[str],
    binary: str = "tblastn",
    evalue: float = DEFAULT_EVALUE,
    timeout: int = TBLASTN_TIMEOUT,
) -> dict[str, TblastnHit]:
    """Align every gene in a multi-query FASTA against one genome, in one call.

    tblastn translates all six reading frames of the subject on every
    invocation. Running one call per gene therefore re-translates the whole
    genome once per gene, which is the bulk of the cost at cohort scale. One
    call with every gene as a query does that work once.

    Args:
        query_faa: FASTA whose record ids are the gene names.
        genome_fna: the subject assembly.
        genes: the gene names expected as query ids. Used to reject output for
            an id we did not ask for, rather than silently dropping it.

    Returns:
        The best hit per gene. A gene with no alignment is absent from the
        result, which is how the caller tells ``no_hit`` from a real error.

    Raises:
        TblastnError: tblastn failed, timed out, or produced unreadable output.
    """
    output = _run_tblastn(_tblastn_command(query_faa, genome_fna, binary, evalue), timeout)

    expected = set(genes)
    best: dict[str, tuple[float, TblastnHit]] = {}
    for line in output.splitlines():
        if not line.strip():
            continue
        query_id, bitscore, hit = _parse_row(line)
        if hit is None:
            continue
        if query_id not in expected:
            raise TblastnError(
                f"tblastn returned an unexpected query id {query_id!r}; "
                f"expected one of {sorted(expected)}"
            )
        # -max_hsps 1 should make this a no-op, but selecting on bitscore means
        # a future BLAST that reports more cannot quietly change which HSP wins.
        current = best.get(query_id)
        if current is None or bitscore > current[0]:
            best[query_id] = (bitscore, hit)

    return {gene: hit for gene, (_, hit) in best.items()}


@dataclass(frozen=True)
class _Task:
    """One (genome, gene) pair, aligned by its own tblastn call."""

    accession: str
    genome: Path
    gene: str
    reference: Path | None


@dataclass(frozen=True)
class _GenomeTask:
    """One genome against the whole gene panel."""

    accession: str
    genome: Path
    #: ``(gene, reference)`` in panel order. ``reference`` is ``None`` for a gene
    #: with no prepared reference, which still owes the manifest a row.
    genes: tuple[tuple[str, Path | None], ...]


def _build_query_fasta(references: dict[str, Path | None], destination: Path) -> list[str]:
    """Write every usable reference into one FASTA, id = gene name.

    Returns the gene names actually written, in panel order.

    Raises:
        TblastnError: a reference is missing, empty, or holds more than one
            record. Two records under one id would become two queries sharing an
            id, and the per-query mapping downstream would silently mis-assign.
    """
    written: list[str] = []
    with destination.open("w", encoding="utf-8") as handle:
        for gene, path in references.items():
            if path is None:
                continue
            if not gene or any(char.isspace() for char in gene):
                raise TblastnError(
                    f"Gene name {gene!r} cannot be used as a FASTA record id. "
                    "Rename the target or drop it from the panel."
                )
            try:
                text = path.read_text(encoding="utf-8")
            except OSError as exc:
                raise TblastnError(f"Could not read reference for {gene}: {exc}") from exc

            records = text.split(">")[1:]
            if len(records) != 1:
                raise TblastnError(
                    f"Reference for {gene} ({path.name}) holds {len(records)} records; "
                    "expected exactly 1. Split it into separate target genes."
                )
            sequence = "".join(records[0].splitlines()[1:]).strip()
            if not sequence:
                raise TblastnError(f"Reference for {gene} ({path.name}) is empty")

            handle.write(f">{gene}\n{sequence}\n")
            written.append(gene)

    return written


def _write_protein(accession: str, gene: str, proteins_dir: Path, hit: TblastnHit) -> str:
    """Write one protein. Returns an error detail, or an empty string."""
    try:
        write_protein_record(
            protein_path(proteins_dir, accession, gene),
            hit.sequence,
            seq_id=f"{accession}_{gene}",
            description=(
                f"gene={gene} accession={accession} contig={hit.contig} pident={hit.pident}"
            ),
        )
    except (OSError, ValueError) as exc:
        logger.error("%s/%s: could not write protein: %s", accession, gene, exc)
        return f"write failed: {exc}"
    return ""


def _run_genome_task(
    task: _GenomeTask,
    proteins_dir: Path,
    query_faa: Path | None,
    genes_in_query: tuple[str, ...],
    binary: str,
    evalue: float,
) -> list[dict]:
    """Align the whole panel against one genome; return one row per gene.

    Module level so the thread pool can call it directly.
    """
    no_ref = "no reference available"
    results: dict[str, dict] = {
        gene: _row(task.accession, gene, STATUS_NO_REFERENCE, no_ref)
        for gene, reference in task.genes
        if reference is None
    }

    wanted = [gene for gene, reference in task.genes if reference is not None]

    hits: dict[str, TblastnHit] = {}
    if wanted and query_faa is not None:
        try:
            hits = run_tblastn_batch(
                query_faa, task.genome, genes_in_query, binary=binary, evalue=evalue
            )
        except TblastnError as exc:
            # One call failed, so every gene in this genome failed. Reporting it
            # per row keeps "we could not look" from reading as "we looked and
            # found nothing" -- the distinction the manifest exists to preserve.
            logger.warning("%s: %s", task.accession, exc)
            return [
                _row(task.accession, gene, STATUS_ERROR, str(exc))
                if reference is not None
                else _row(task.accession, gene, STATUS_NO_REFERENCE, no_ref)
                for gene, reference in task.genes
            ]

    for gene in wanted:
        hit = hits.get(gene)
        if hit is None:
            results[gene] = _row(task.accession, gene, STATUS_NO_HIT, "no alignment above e-value")
            continue

        detail = _write_protein(task.accession, gene, proteins_dir, hit)
        if detail:
            results[gene] = _row(task.accession, gene, STATUS_ERROR, detail, hit)
            continue

        if hit.internal_stops:
            logger.info(
                "%s/%s: %d internal stop codon(s) in the extracted protein",
                task.accession,
                gene,
                hit.internal_stops,
            )
        results[gene] = _row(task.accession, gene, STATUS_EXTRACTED, "", hit)

    return [results[gene] for gene, _ in task.genes]


def _row(
    accession: str,
    gene: str,
    status: str,
    detail: str = "",
    hit: TblastnHit | None = None,
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
    """One (genome, gene) alignment via the single-query path.

    Retained so the batched path can be checked against it, rather than only
    against expectations.
    """
    if task.reference is None:
        return _row(task.accession, task.gene, STATUS_NO_REFERENCE, "no reference available")

    try:
        hit = run_tblastn(task.reference, task.genome, binary=binary, evalue=evalue)
    except TblastnError as exc:
        logger.warning("%s/%s: %s", task.accession, task.gene, exc)
        return _row(task.accession, task.gene, STATUS_ERROR, str(exc))

    if hit is None:
        return _row(task.accession, task.gene, STATUS_NO_HIT, "no alignment above e-value")

    detail = _write_protein(task.accession, task.gene, proteins_dir, hit)
    if detail:
        return _row(task.accession, task.gene, STATUS_ERROR, detail, hit)

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

    gene_names = (
        list(targets.genes)
        if isinstance(targets, ReferenceSet)
        else [str(gene).strip() for gene in targets if str(gene).strip()]
    )

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
        _GenomeTask(
            accession=genome.stem,
            genome=genome,
            genes=tuple((gene, references[gene]) for gene in gene_names),
        )
        for genome in genome_files
    ]

    rows: list[dict] = []
    if tasks:
        verify_tblastn(tblastn_binary)
        logger.info(
            "Extracting %d gene(s) from %d genome(s) = %d alignments on %d thread(s), "
            "batched to %d tblastn call(s)",
            len(gene_names),
            len(genome_files),
            len(tasks) * len(gene_names),
            threads,
            len(tasks),
        )

        with tempfile.TemporaryDirectory(prefix="mutationscan-") as tmp:
            query_faa = Path(tmp) / "panel.faa"
            genes_in_query = _build_query_fasta(references, query_faa)
            if len(genes_in_query) < len(gene_names):
                logger.warning(
                    "%d of %d gene(s) had no prepared reference and were not aligned",
                    len(gene_names) - len(genes_in_query),
                    len(gene_names),
                )
            logger.info(
                "Query panel: %d gene(s) in one FASTA, so each genome's six reading "
                "frames are translated once instead of once per gene",
                len(genes_in_query),
            )

            # Threads, not processes: tblastn is an external process, so the GIL
            # is released for essentially the whole task. A process pool would add
            # interpreter startup and pickle the task for no gain. tblastn is
            # pinned to one thread each, so the concurrency is ours to control.
            batches = list(
                parallel_map(
                    partial(
                        _run_genome_task,
                        proteins_dir=proteins_out,
                        query_faa=query_faa,
                        genes_in_query=tuple(genes_in_query),
                        binary=tblastn_binary,
                        evalue=evalue,
                    ),
                    tasks,
                    threads=threads,
                    workload=WORKLOAD_THREAD,
                    total=len(tasks),
                    label="extraction",
                )
            )
        rows = [row for batch in batches for row in batch]

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
