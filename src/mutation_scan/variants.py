"""Substitution calling by global alignment against the wild-type reference.

The algorithm, unchanged in substance from earlier MutationScan versions:

* global (Needleman-Wunsch) alignment of reference against query, BLOSUM62,
  gap open -10.0, gap extend -0.5 -- penalties chosen so that a point change is
  scored as a substitution rather than an indel pair;
* the **residue counter**: walk the alignment column by column and advance the
  position counter *only* when the reference column is not a gap. Positions are
  therefore reference coordinates -- position 83 is the 83rd residue of the
  reference protein, not the 83rd column of the alignment;
* a substitution is reported when reference and query differ and neither is a
  gap. Insertions and deletions are counted in the QC table but not reported as
  mutations.

Two gates decide whether a (genome, gene) pair is evaluated at all: identity to
the reference, and coverage of the reference. Coverage is what catches a
truncated or fragmented gene, which identity alone will happily call perfect.

Nothing here interprets a mutation. There is no resistance database, no drug
mapping, no severity score and no prediction: the output is what the alignment
says, and nothing more.
"""

from __future__ import annotations

import logging
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices

from .extract import STATUS_EXTRACTED, protein_path
from .fasta import read_protein_sequence
from .parallel import WORKLOAD_PROCESS, ParallelTaskError, parallel_map
from .references import resolve_reference

logger = logging.getLogger(__name__)

__all__ = [
    "GAP",
    "MUTATION_COLUMNS",
    "QC_COLUMNS",
    "VARIANT_STATUS_CALLED",
    "VARIANT_STATUS_ERROR",
    "VARIANT_STATUS_LOW_COVERAGE",
    "VARIANT_STATUS_LOW_IDENTITY",
    "VARIANT_STATUS_MISSING_PROTEIN",
    "VARIANT_STATUS_NO_REFERENCE",
    "Comparison",
    "VariantResult",
    "build_aligner",
    "call_variants",
    "compare_sequences",
    "read_mutations",
    "read_variant_qc",
]

GAP = "-"

GAP_OPEN_SCORE = -10.0
GAP_EXTEND_SCORE = -0.5

VARIANT_STATUS_CALLED = "called"
VARIANT_STATUS_LOW_IDENTITY = "low_identity"
VARIANT_STATUS_LOW_COVERAGE = "low_coverage"
VARIANT_STATUS_NO_REFERENCE = "no_reference"
VARIANT_STATUS_MISSING_PROTEIN = "missing_protein"
VARIANT_STATUS_ERROR = "error"

MUTATION_COLUMNS = (
    "Accession",
    "Gene",
    "Mutation",
    "Ref_AA",
    "Position",
    "Alt_AA",
    "Identity_pct",
    "Coverage_pct",
)

QC_COLUMNS = (
    "Accession",
    "Gene",
    "Status",
    "Identity_pct",
    "Coverage_pct",
    "Ref_Length",
    "Query_Length",
    "N_Mutations",
    "N_Deleted_Positions",
    "N_Inserted_Columns",
    "Detail",
)


def build_aligner() -> PairwiseAligner:
    """A global BLOSUM62 aligner with MutationScan's gap penalties."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = GAP_OPEN_SCORE
    aligner.extend_gap_score = GAP_EXTEND_SCORE
    return aligner


def _alignable_alphabet(aligner: PairwiseAligner) -> frozenset[str]:
    matrix = aligner.substitution_matrix
    return frozenset(matrix.alphabet) if matrix is not None else frozenset()


def _coerce_residues(sequence: str, alphabet: frozenset[str], label: str) -> str:
    """Replace residues the scoring matrix does not know with ``X``.

    Rare (``U``, ``O``, ``J`` and junk characters), but one odd residue must not
    abort a whole cohort.
    """
    if not alphabet:
        return sequence
    unknown = sorted({residue for residue in sequence if residue not in alphabet})
    if not unknown:
        return sequence
    logger.warning(
        "%s: residue(s) %s are not in the scoring matrix; treating them as X",
        label,
        "".join(unknown),
    )
    return "".join(residue if residue in alphabet else "X" for residue in sequence)


@dataclass(frozen=True)
class Comparison:
    """The result of aligning one query protein to one reference."""

    identity_pct: float
    coverage_pct: float
    ref_length: int
    query_length: int
    mutations: tuple[tuple[str, int, str], ...]
    deleted_positions: int
    inserted_columns: int


def compare_sequences(
    reference: str,
    query: str,
    aligner: PairwiseAligner | None = None,
) -> Comparison:
    """Align *query* to *reference* and return substitutions plus alignment stats.

    Identity and coverage are both expressed against the reference length, so a
    query covering half the reference perfectly scores ~50% identity and ~50%
    coverage -- which is how a truncated gene is distinguished from a diverged
    one.
    """
    aligner = aligner or build_aligner()
    alphabet = _alignable_alphabet(aligner)

    reference = _coerce_residues(reference.upper(), alphabet, "reference")
    query = _coerce_residues(query.upper(), alphabet, "query")

    if not reference:
        raise ValueError("reference sequence is empty")
    if not query:
        raise ValueError("query sequence is empty")

    # Order matters: alignment[0] is the reference, alignment[1] the query.
    # Index rather than list() -- a pair of long sequences can have an enormous
    # number of equally optimal alignments, and materializing them all is what
    # made the previous implementation liable to hang.
    alignments = aligner.align(reference, query)
    alignment = alignments[0]

    aligned_ref = str(alignment[0])
    aligned_query = str(alignment[1])
    if len(aligned_ref) != len(aligned_query):  # pragma: no cover - defensive
        raise ValueError("aligner returned rows of unequal length")

    ref_length = len(reference)
    identities = 0
    covered = 0
    deleted = 0
    inserted = 0
    mutations: list[tuple[str, int, str]] = []
    position = 0  # becomes 1 at the first non-gap reference column

    for ref_aa, query_aa in zip(aligned_ref, aligned_query, strict=False):
        if ref_aa != GAP:
            position += 1
            if query_aa != GAP:
                covered += 1
                if ref_aa == query_aa:
                    identities += 1
                else:
                    mutations.append((ref_aa, position, query_aa))
            else:
                deleted += 1
        elif query_aa != GAP:
            inserted += 1

    return Comparison(
        identity_pct=round(identities / ref_length * 100.0, 2) if ref_length else 0.0,
        coverage_pct=round(covered / ref_length * 100.0, 2) if ref_length else 0.0,
        ref_length=ref_length,
        query_length=len(query),
        mutations=tuple(mutations),
        deleted_positions=deleted,
        inserted_columns=inserted,
    )


@dataclass
class VariantResult:
    """Mutations found, plus one QC row per (genome, gene) pair considered."""

    mutations: pd.DataFrame
    qc: pd.DataFrame

    def __len__(self) -> int:
        return len(self.mutations)


def _pairs_from_manifest(manifest: pd.DataFrame) -> Iterator[tuple[str, str]]:
    missing = [
        column for column in ("Accession", "Gene", "Status") if column not in manifest.columns
    ]
    if missing:
        if manifest.empty:  # an empty run: nothing was attempted, nothing to call
            return
        raise ValueError(f"Manifest is missing column(s): {', '.join(missing)}")

    extracted = manifest[manifest["Status"] == STATUS_EXTRACTED]
    for row in extracted.itertuples():
        yield str(row.Accession), str(row.Gene)


def _pairs_from_directory(proteins_dir: Path) -> Iterator[tuple[str, str, Path]]:
    """Fallback discovery when no manifest is available.

    The gene name is read from the record description (``gene=...``) that
    extraction writes, falling back to splitting the filename on its last
    underscore. Reading the header first means a gene whose name contains an
    underscore is still handled correctly.
    """
    for path in sorted(proteins_dir.glob("*.faa")):
        gene: str | None = None
        accession: str | None = None
        try:
            record = read_protein_sequence(path)
        except ValueError:
            record = None
        if record is not None:
            tokens = dict(
                token.split("=", 1) for token in record.description.split() if token.count("=") == 1
            )
            gene = tokens.get("gene")
            accession = tokens.get("accession")

        if not gene or not accession:
            parts = path.stem.rsplit("_", 1)
            if len(parts) != 2:
                logger.warning(
                    "Skipping %s: cannot determine accession and gene "
                    "(expected {accession}_{gene}.faa)",
                    path.name,
                )
                continue
            accession, gene = parts

        yield accession, gene, path


def _qc_row(
    accession: str,
    gene: str,
    status: str,
    comparison: Comparison | None = None,
    n_mutations: int = 0,
    detail: str = "",
) -> dict:
    return {
        "Accession": accession,
        "Gene": gene,
        "Status": status,
        "Identity_pct": comparison.identity_pct if comparison else pd.NA,
        "Coverage_pct": comparison.coverage_pct if comparison else pd.NA,
        "Ref_Length": comparison.ref_length if comparison else pd.NA,
        "Query_Length": comparison.query_length if comparison else pd.NA,
        "N_Mutations": n_mutations if status == VARIANT_STATUS_CALLED else pd.NA,
        "N_Deleted_Positions": comparison.deleted_positions if comparison else pd.NA,
        "N_Inserted_Columns": comparison.inserted_columns if comparison else pd.NA,
        "Detail": detail,
    }


#: Per-process state, populated once by the pool initializer and then reused for
#: every task that worker handles. Building a PairwiseAligner means loading
#: BLOSUM62, and re-reading a reference means disk plus parsing; neither should
#: happen once per (genome, gene) pair when there are thousands of them.
_WORKER_ALIGNER: PairwiseAligner | None = None
_WORKER_REFERENCES: dict[str, str] = {}


def _init_variant_worker() -> None:
    """Pool initializer: build this worker's aligner once, before any task."""
    global _WORKER_ALIGNER
    _WORKER_ALIGNER = build_aligner()


def _worker_aligner() -> PairwiseAligner:
    global _WORKER_ALIGNER
    aligner = _WORKER_ALIGNER
    if aligner is None:  # serial path, or an initializer that did not run
        aligner = build_aligner()
        _WORKER_ALIGNER = aligner
    return aligner


def _worker_reference(gene: str, reference_path: str) -> str:
    """Read a reference protein, memoised for the life of this process."""
    cached = _WORKER_REFERENCES.get(reference_path)
    if cached is None:
        cached = str(read_protein_sequence(Path(reference_path), gene).seq)
        _WORKER_REFERENCES[reference_path] = cached
    return cached


@dataclass(frozen=True)
class _VariantTask:
    """One unit of work for the variant-calling stage.

    Carries *paths*, not sequences. Every task is small to pickle no matter how
    long the proteins are, and nothing is held in memory that has not been
    reached yet -- which is what lets the whole stage stream from a manifest of
    arbitrary size instead of materialising every query sequence up front.

    A task with a non-empty ``status`` has already been decided (a missing
    reference, a missing protein) and passes straight through the worker, so the
    QC table has exactly one row per pair no matter which path it took.
    """

    accession: str
    gene: str
    reference_path: str = ""
    protein_path: str = ""
    min_identity: float = 0.0
    min_coverage: float = 0.0
    status: str = ""
    detail: str = ""


def _process_variant_task(task: _VariantTask) -> dict:
    """Evaluate one (genome, gene) pair and apply the quality gates.

    Never raises for a data problem: a pair that cannot be evaluated comes back
    as an ``error`` result, so one unreadable protein cannot cost the run. A
    genuine defect still propagates, because a worker that dies must not be
    mistaken for a dataset quirk.
    """
    if task.status:
        return {
            "accession": task.accession,
            "gene": task.gene,
            "status": task.status,
            "comparison": None,
            "mutations": [],
            "detail": task.detail,
        }

    try:
        reference = _worker_reference(task.gene, task.reference_path)
        query = str(
            read_protein_sequence(Path(task.protein_path), f"{task.accession}_{task.gene}").seq
        )
        comparison = compare_sequences(reference, query, aligner=_worker_aligner())
    except (ValueError, OSError) as exc:
        return {
            "accession": task.accession,
            "gene": task.gene,
            "status": VARIANT_STATUS_ERROR,
            "comparison": None,
            "mutations": [],
            "detail": str(exc),
        }

    if comparison.identity_pct < task.min_identity:
        return {
            "accession": task.accession,
            "gene": task.gene,
            "status": VARIANT_STATUS_LOW_IDENTITY,
            "comparison": comparison,
            "mutations": [],
            "detail": f"identity {comparison.identity_pct:.2f}% < {task.min_identity:.2f}%",
        }

    if comparison.coverage_pct < task.min_coverage:
        return {
            "accession": task.accession,
            "gene": task.gene,
            "status": VARIANT_STATUS_LOW_COVERAGE,
            "comparison": comparison,
            "mutations": [],
            "detail": f"coverage {comparison.coverage_pct:.2f}% < {task.min_coverage:.2f}%",
        }

    mutations = [
        {
            "Accession": task.accession,
            "Gene": task.gene,
            "Mutation": f"{ref_aa}{position}{alt_aa}",
            "Ref_AA": ref_aa,
            "Position": position,
            "Alt_AA": alt_aa,
            "Identity_pct": comparison.identity_pct,
            "Coverage_pct": comparison.coverage_pct,
        }
        for ref_aa, position, alt_aa in comparison.mutations
    ]

    return {
        "accession": task.accession,
        "gene": task.gene,
        "status": VARIANT_STATUS_CALLED,
        "comparison": comparison,
        "mutations": mutations,
        "detail": "",
    }


def call_variants(
    proteins_dir: Path | str,
    refs_dir: Path | str,
    manifest: pd.DataFrame | None = None,
    min_identity_percent: float = 80.0,
    min_coverage_percent: float = 80.0,
    mutations_csv: Path | str | None = None,
    qc_csv: Path | str | None = None,
    threads: int = 1,
) -> VariantResult:
    """Call substitutions for every extracted protein.

    Args:
        proteins_dir: directory of ``{accession}_{gene}.faa`` files.
        refs_dir: prepared references, resolved the same way extraction does.
        manifest: extraction manifest. When given, it drives iteration, so the
            accession/gene of every protein is known exactly rather than parsed
            out of a filename. Falls back to scanning *proteins_dir*.
        min_identity_percent: reject the pair below this identity to the reference.
        min_coverage_percent: reject the pair below this coverage of the reference.
        mutations_csv / qc_csv: optional output paths.
        threads: number of parallel worker processes (default: 1).

    Returns:
        :class:`VariantResult` -- the mutation table, and a QC row for every pair
        considered (including those rejected by a gate, so a downstream stage can
        tell "no mutations" apart from "not evaluated").
    """
    proteins_dir = Path(proteins_dir)
    refs_dir = Path(refs_dir)

    if manifest is not None:
        pairs: Iterator[tuple[str, str, Path]] = (
            (accession, gene, protein_path(proteins_dir, accession, gene))
            for accession, gene in _pairs_from_manifest(manifest)
        )
        total: int | None = (
            int((manifest["Status"] == STATUS_EXTRACTED).sum())
            if not manifest.empty and "Status" in manifest.columns
            else None
        )
    else:
        discovered = list(_pairs_from_directory(proteins_dir))
        pairs = iter(discovered)
        total = len(discovered)

    logger.info(
        "Calling variants for %s (genome, gene) pair(s); "
        "gates: identity >= %.1f%%, coverage >= %.1f%% on %d thread(s)",
        total if total is not None else "all discovered",
        min_identity_percent,
        min_coverage_percent,
        threads,
    )

    mutation_rows: list[dict] = []
    qc_rows: list[dict] = []

    # Resolve each reference path once here, in the parent. Only the *path* is
    # handed to workers, which read and memoise the sequence themselves: a
    # reference is at most a few kB but a full cohort's worth of them would
    # otherwise be pickled once per task.
    reference_paths: dict[str, str | None] = {}

    def reference_path_for(gene: str) -> str | None:
        if gene not in reference_paths:
            path = resolve_reference(refs_dir, gene)
            if path is None:
                logger.error("No reference for %s in %s", gene, refs_dir)
                reference_paths[gene] = None
            else:
                # Fail here rather than in every worker, so a bad reference is
                # reported once instead of once per genome that lacks it.
                try:
                    read_protein_sequence(path, gene)
                except ValueError as exc:
                    logger.error("Unusable reference for %s: %s", gene, exc)
                    reference_paths[gene] = None
                else:
                    reference_paths[gene] = str(path)
        return reference_paths[gene]

    def _tasks() -> Iterator[_VariantTask]:
        """Stream tasks. Nothing is accumulated, so memory is flat in cohort size."""
        for accession, gene, faa_path in pairs:
            reference = reference_path_for(gene)
            if reference is None:
                yield _VariantTask(
                    accession=accession,
                    gene=gene,
                    status=VARIANT_STATUS_NO_REFERENCE,
                    detail="no usable reference",
                )
                continue
            if not faa_path.is_file():
                yield _VariantTask(
                    accession=accession,
                    gene=gene,
                    status=VARIANT_STATUS_MISSING_PROTEIN,
                    detail=f"missing {faa_path.name}",
                )
                continue
            yield _VariantTask(
                accession=accession,
                gene=gene,
                reference_path=reference,
                protein_path=str(faa_path),
                min_identity=min_identity_percent,
                min_coverage=min_coverage_percent,
            )

    try:
        for result in parallel_map(
            _process_variant_task,
            _tasks(),
            threads=threads,
            workload=WORKLOAD_PROCESS,
            total=total,
            initializer=_init_variant_worker,
            label="variant calling",
        ):
            _collect_result(result, qc_rows, mutation_rows)
    except ParallelTaskError as exc:
        # A worker dying is an infrastructure failure, not a data finding, and
        # the partial tables would be indistinguishable from a real result.
        logger.error("Variant calling aborted after %d pair(s): %s", exc.completed, exc)
        raise

    mutations = pd.DataFrame(mutation_rows, columns=list(MUTATION_COLUMNS))
    if not mutations.empty:
        mutations = mutations.sort_values(
            ["Accession", "Gene", "Position"], kind="stable"
        ).reset_index(drop=True)

    qc = pd.DataFrame(qc_rows, columns=list(QC_COLUMNS))
    if not qc.empty:
        qc = qc.sort_values(["Accession", "Gene"], kind="stable").reset_index(drop=True)

    for frame, path in ((mutations, mutations_csv), (qc, qc_csv)):
        if path is not None:
            path = Path(path)
            path.parent.mkdir(parents=True, exist_ok=True)
            frame.to_csv(path, index=False)
            logger.info("Wrote %s (%d rows)", path, len(frame))

    _log_summary(mutations, qc)
    return VariantResult(mutations=mutations, qc=qc)


def _collect_result(result: dict, qc_rows: list[dict], mutation_rows: list[dict]) -> None:
    """Collect results from a worker into the output lists."""
    accession = result["accession"]
    gene = result["gene"]
    status = result["status"]
    comparison = result["comparison"]
    mutations = result["mutations"]
    detail = result["detail"]

    if status == VARIANT_STATUS_CALLED:
        mutation_rows.extend(mutations)
        qc_rows.append(
            _qc_row(
                accession,
                gene,
                VARIANT_STATUS_CALLED,
                comparison,
                n_mutations=len(mutations),
            )
        )
    else:
        qc_rows.append(_qc_row(accession, gene, status, comparison, detail=detail))


def _log_summary(mutations: pd.DataFrame, qc: pd.DataFrame) -> None:
    if qc.empty:
        logger.warning("No (genome, gene) pairs were evaluated")
        return

    counts = qc["Status"].value_counts()
    logger.info(
        "Variant calling complete: %s",
        ", ".join(f"{status}={count}" for status, count in counts.items()),
    )
    evaluated = qc[qc["Status"] == VARIANT_STATUS_CALLED]
    logger.info(
        "  %d mutation(s) across %d evaluated pair(s) in %d genome(s)",
        len(mutations),
        len(evaluated),
        evaluated["Accession"].nunique() if not evaluated.empty else 0,
    )


def read_variant_qc(qc_csv: Path | str) -> pd.DataFrame:
    """Read a variant QC table, keeping accessions and genes as strings."""
    qc = pd.read_csv(qc_csv, dtype={"Accession": str, "Gene": str, "Status": str})
    missing = [column for column in ("Accession", "Gene", "Status") if column not in qc]
    if missing:
        raise ValueError(f"Variant QC {qc_csv} is missing column(s): {', '.join(missing)}")
    return qc


def read_mutations(mutations_csv: Path | str) -> pd.DataFrame:
    """Read a mutations table, keeping accessions and genes as strings."""
    mutations = pd.read_csv(mutations_csv, dtype={"Accession": str, "Gene": str, "Mutation": str})
    missing = [column for column in ("Accession", "Gene", "Mutation") if column not in mutations]
    if missing:
        raise ValueError(
            f"Mutations file {mutations_csv} is missing column(s): {', '.join(missing)}"
        )
    return mutations
