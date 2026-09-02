"""Target discovery and reference preparation.

This module is the reason MutationScan needs no gene list in its code: the
targets of a run are *whatever reference proteins the user supplies*. Drop
``geneA_WT.faa`` into the references directory and ``geneA`` becomes a target;
drop ``geneB.fasta`` in and so does ``geneB``. Gene names are only ever read from
filenames, never assumed.

References are normalized on the way into the run directory (one clean FASTA
record per gene, named ``{gene}_WT.faa``) so that BLAST+ and the aligner both
receive well-formed input, and so the run folder records exactly which
reference every call was made against.
"""

from __future__ import annotations

import logging
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional, Sequence

from .fasta import read_protein_sequence, sanitize_protein_sequence, write_protein_record

logger = logging.getLogger(__name__)

__all__ = [
    "REFERENCE_EXTENSIONS",
    "ReferenceSet",
    "canonical_reference_path",
    "discover_targets",
    "fetch_uniprot_reference",
    "gene_from_ref_stem",
    "load_reference_set",
    "missing_reference_genes",
    "prepare_references",
    "resolve_reference",
]

#: Filename extensions treated as protein references.
REFERENCE_EXTENSIONS = frozenset({".faa", ".fasta", ".fa", ".fas"})

_UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
_UNIPROT_TIMEOUT = 20
_UNIPROT_PAUSE = 0.5  # be a polite API citizen between fetches


def gene_from_ref_stem(stem: str) -> str:
    """Normalize a reference filename stem to a gene name (``geneA_WT`` -> ``geneA``)."""
    if stem.lower().endswith("_wt"):
        return stem[:-3]
    return stem


def _reference_candidates(directory: Path) -> list[Path]:
    """Non-empty reference files in *directory*, in a deterministic order."""
    if not directory.is_dir():
        return []
    candidates = []
    for path in sorted(directory.iterdir()):
        if path.suffix.lower() not in REFERENCE_EXTENSIONS:
            continue
        try:
            if not path.is_file() or path.stat().st_size == 0:
                continue
        except OSError:  # pragma: no cover - unreadable entry, skip it
            continue
        candidates.append(path)
    return candidates


def resolve_reference(refs_dir: Path | str, gene: str) -> Optional[Path]:
    """Find the reference file for *gene*, case-insensitively.

    Accepts ``{gene}_WT.faa``, ``{gene}.faa`` and the ``.fasta``/``.fa``/``.fas``
    equivalents. Returns ``None`` when the gene has no usable reference.
    """
    gene_key = str(gene).strip().lower()
    if not gene_key:
        return None

    for path in _reference_candidates(Path(refs_dir)):
        if gene_from_ref_stem(path.stem).lower() == gene_key:
            return path
    return None


def canonical_reference_path(refs_dir: Path | str, gene: str) -> Path:
    """Where a reference for *gene* is stored inside a run directory."""
    return Path(refs_dir) / f"{gene}_WT.faa"


def discover_targets(
    references_dir: Path | str,
    explicit_targets: Sequence[str] | None = None,
) -> list[str]:
    """Determine which genes a run targets.

    If *explicit_targets* is non-empty it is used verbatim (order preserved,
    de-duplicated case-insensitively) -- this is how a user analyses a subset of
    the references they keep on disk. Otherwise every gene with a reference file
    in *references_dir* is a target, sorted for reproducibility.
    """
    if explicit_targets:
        seen: set[str] = set()
        targets: list[str] = []
        for item in explicit_targets:
            name = str(item).strip()
            if not name or name.lower() in seen:
                continue
            seen.add(name.lower())
            targets.append(name)
        return targets

    discovered: dict[str, str] = {}
    for path in _reference_candidates(Path(references_dir)):
        gene = gene_from_ref_stem(path.stem)
        # First file wins; _reference_candidates is sorted, so this is stable.
        discovered.setdefault(gene.lower(), gene)
    return [discovered[key] for key in sorted(discovered)]


def missing_reference_genes(refs_dir: Path | str, genes: Iterable[str]) -> list[str]:
    """Genes with no resolvable reference in *refs_dir*."""
    return [gene for gene in genes if resolve_reference(refs_dir, gene) is None]


def fetch_uniprot_reference(gene: str, taxid: str) -> Optional[str]:
    """Fetch one reviewed UniProt protein FASTA for *gene* in *taxid*.

    Returns the raw FASTA text, or ``None`` if nothing matched or the request
    failed. Network problems are logged, never raised: a missing reference is
    reported as a missing target, which the caller surfaces to the user.
    """
    query = f"gene:{gene} AND taxonomy_id:{taxid} AND reviewed:true"
    url = (
        f"{_UNIPROT_SEARCH}?"
        + urllib.parse.urlencode({"query": query, "format": "fasta", "size": "1"})
    )

    try:
        request = urllib.request.Request(url, headers={"Accept": "text/plain"})
        with urllib.request.urlopen(request, timeout=_UNIPROT_TIMEOUT) as response:
            text = response.read().decode("utf-8").strip()
    except urllib.error.HTTPError as exc:
        logger.error("UniProt HTTP %s fetching %s (taxid %s): %s", exc.code, gene, taxid, exc.reason)
        return None
    except urllib.error.URLError as exc:
        logger.error("UniProt unreachable while fetching %s: %s", gene, exc.reason)
        return None
    except Exception as exc:  # noqa: BLE001 - network layers raise broadly
        logger.error("Failed to fetch UniProt reference for %s: %s", gene, exc)
        return None

    if not text:
        logger.warning("No reviewed UniProt entry for gene %s in taxid %s", gene, taxid)
        return None
    return text


@dataclass
class ReferenceSet:
    """The references a run will actually use."""

    refs_dir: Path
    genes: list[str] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)
    lengths: dict[str, int] = field(default_factory=dict)
    sources: dict[str, str] = field(default_factory=dict)

    def path(self, gene: str) -> Path:
        return canonical_reference_path(self.refs_dir, gene)

    def __len__(self) -> int:
        return len(self.genes)


def prepare_references(
    references_dir: Path | str,
    refs_dir: Path | str,
    targets: Sequence[str],
    uniprot_taxid: str = "",
) -> ReferenceSet:
    """Stage a normalized reference for every target into the run directory.

    For each target: use the local reference if present, otherwise fetch it from
    UniProt when *uniprot_taxid* is set. Whatever is found is rewritten as a
    single clean FASTA record at ``{refs_dir}/{gene}_WT.faa``.

    Genes that end up without a reference are returned in
    :attr:`ReferenceSet.missing` rather than raising, so the caller can report
    all of them at once.
    """
    references_dir = Path(references_dir)
    refs_dir = Path(refs_dir)
    refs_dir.mkdir(parents=True, exist_ok=True)

    result = ReferenceSet(refs_dir=refs_dir)
    fetched_any = False

    for gene in targets:
        source = resolve_reference(references_dir, gene)
        origin = f"local:{source.name}" if source is not None else ""
        fasta_text: Optional[str] = None

        if source is None and uniprot_taxid:
            if fetched_any:
                time.sleep(_UNIPROT_PAUSE)
            logger.info("No local reference for %s; querying UniProt (taxid %s)", gene, uniprot_taxid)
            fasta_text = fetch_uniprot_reference(gene, uniprot_taxid)
            fetched_any = True
            if fasta_text:
                origin = f"uniprot:{uniprot_taxid}"

        if source is None and fasta_text is None:
            logger.warning("No reference available for target %s", gene)
            result.missing.append(gene)
            continue

        try:
            if source is not None:
                record = read_protein_sequence(source, fallback_id=gene)
                sequence = str(record.seq)
                header_note = record.description or ""
                if header_note.startswith("headerless_source:"):
                    # Redundant with the local: origin we already record.
                    header_note = ""
            else:
                sequence, header_note = _parse_fetched_fasta(fasta_text or "")
        except ValueError as exc:
            logger.error("Unusable reference for %s: %s", gene, exc)
            result.missing.append(gene)
            continue

        if not sequence:
            logger.error("Reference for %s contained no residues", gene)
            result.missing.append(gene)
            continue

        destination = canonical_reference_path(refs_dir, gene)
        write_protein_record(
            destination,
            sequence,
            seq_id=f"{gene}_WT",
            description=f"{origin} {header_note}".strip(),
        )
        result.genes.append(gene)
        result.lengths[gene] = len(sequence)
        result.sources[gene] = origin
        logger.info("Reference ready: %s (%d aa, %s)", gene, len(sequence), origin)

    return result


#: Origins recorded in a prepared reference's FASTA header by prepare_references.
_ORIGIN_PREFIXES = ("local:", "uniprot:")


def _origin_from_description(description: str) -> str:
    """Recover the ``local:``/``uniprot:`` origin token from a header line."""
    for token in description.split():
        if token.startswith(_ORIGIN_PREFIXES):
            return token
    return ""


def load_reference_set(refs_dir: Path | str) -> ReferenceSet:
    """Rebuild a :class:`ReferenceSet` from a run's prepared references directory.

    :func:`prepare_references` normalizes every reference to ``{gene}_WT.faa`` and
    records where it came from in the header, so a later stage can report on what
    was actually used without re-resolving paths or re-fetching anything. That is
    what keeps the Snakemake wrappers free of logic.

    Nothing is reported as *missing* here: a gene with no file in *refs_dir* was
    never prepared, so it is simply absent.
    """
    refs_dir = Path(refs_dir)
    result = ReferenceSet(refs_dir=refs_dir)

    for path in _reference_candidates(refs_dir):
        gene = gene_from_ref_stem(path.stem)
        try:
            record = read_protein_sequence(path, fallback_id=gene)
        except ValueError as exc:
            logger.warning("Ignoring unreadable reference %s: %s", path.name, exc)
            continue
        result.genes.append(gene)
        result.lengths[gene] = len(record.seq)
        result.sources[gene] = _origin_from_description(record.description)

    return result


def _parse_fetched_fasta(text: str) -> tuple[str, str]:
    """Split downloaded FASTA text into (sequence, header description)."""
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        raise ValueError("empty FASTA response")

    header = ""
    body: list[str] = []
    for line in lines:
        if line.startswith(">"):
            if body:
                break  # only the first record
            header = line[1:].strip()
            continue
        body.append(line)

    sequence = sanitize_protein_sequence("".join(body))
    if not sequence:
        raise ValueError("FASTA response contained no residues")
    return sequence, header
