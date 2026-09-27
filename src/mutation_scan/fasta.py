"""Robust protein FASTA reading and writing.

Reference files in the wild are not always well-formed FASTA: they may be
headerless, wrapped at odd widths, or grouped into whitespace-separated blocks
of residues. BLAST+ and downstream alignment both want a clean single record,
so everything MutationScan reads goes through :func:`read_protein_sequence`,
and every reference it stores is rewritten with :func:`write_protein_record`.
"""

from __future__ import annotations

import io
import logging
import re
from pathlib import Path
from typing import cast

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

__all__ = [
    "looks_like_nucleotide",
    "read_first_record_sequence",
    "read_protein_sequence",
    "sanitize_protein_sequence",
    "write_protein_record",
]

#: Standard 20 residues plus ambiguity codes and the stop symbol. Anything else
#: is reported but tolerated -- we would rather analyse an odd file than refuse.
_PROTEIN_ALPHABET = re.compile(r"[ACDEFGHIKLMNPQRSTVWYBXZJUO*]+")

_WHITESPACE = re.compile(r"\s+")

#: Unambiguous DNA/RNA symbols. Deliberately excludes the IUPAC ambiguity codes
#: R/Y/S/W/K/M/B/D/H/V, every one of which is also an amino-acid letter --
#: counting them as nucleotide evidence is what makes a naive sniffer call a
#: protein a genome.
_NUCLEOTIDE_SYMBOLS = frozenset("ACGTUN")

#: Below this length the composition ratio is noise: a 20-residue peptide of
#: alanine and cysteine is indistinguishable from a short DNA fragment.
NUCLEOTIDE_SNIFF_MIN_LENGTH = 50

#: Real assemblies sit far above this; real proteins far below. The gap between
#: the two populations is wide enough that the exact cut-off does not matter.
NUCLEOTIDE_SNIFF_THRESHOLD = 0.9


def sanitize_protein_sequence(text: str) -> str:
    """Strip whitespace/digits from a residue block and upper-case it.

    Also drops a single trailing stop codon symbol, which reference downloads
    and translated ORFs both commonly carry.
    """
    sequence = _WHITESPACE.sub("", text).upper()
    sequence = re.sub(r"\d", "", sequence)
    if sequence.endswith("*"):
        sequence = sequence[:-1]
    return sequence


def looks_like_nucleotide(sequence: str) -> bool:
    """True when *sequence* is composed almost entirely of A/C/G/T/U/N.

    Used to catch the most common way of misusing this tool: putting protein
    FASTA in the genomes directory, nucleotide FASTA in the references
    directory, or swapping the two outright. Left to itself that mistake costs a
    full ``tblastn`` sweep and produces a run where every pair is ``no_hit``,
    with nothing in the output saying why.

    Sequences shorter than :data:`NUCLEOTIDE_SNIFF_MIN_LENGTH` return ``False``:
    the ratio is not informative at that length, and a false "this is a genome"
    is worse than no opinion.
    """
    residues = sanitize_protein_sequence(sequence)
    if len(residues) < NUCLEOTIDE_SNIFF_MIN_LENGTH:
        return False
    matches = sum(1 for residue in residues if residue in _NUCLEOTIDE_SYMBOLS)
    return matches / len(residues) >= NUCLEOTIDE_SNIFF_THRESHOLD


def read_first_record_sequence(path: Path | str, max_residues: int = 4000) -> str:
    """Residues of the first record in a FASTA file, up to *max_residues*.

    Deliberately not :func:`read_protein_sequence`, which insists on exactly one
    record: a genome assembly is many contigs, and validating one costs nothing
    only if we stop reading after the first. Streams line by line and stops at
    the second header or the residue cap, so cost is bounded on a 5 Mb assembly.

    Returns an empty string when the file has no residues. Raises nothing --
    callers decide whether an unreadable file is fatal.
    """
    path = Path(path)
    collected: list[str] = []
    total = 0
    try:
        with path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if collected:
                        break  # second record: the first one is complete
                    continue
                collected.append(line)
                total += len(line)
                if total >= max_residues:
                    break
    except OSError as exc:
        logger.debug("Could not read %s: %s", path, exc)
        return ""

    return sanitize_protein_sequence("".join(collected))[:max_residues]


def read_protein_sequence(path: Path | str, fallback_id: str = "sequence") -> SeqRecord:
    """Read a single protein sequence, tolerating imperfect FASTA.

    Headed files go through Biopython's FASTA dialects; headerless files are
    read as a bare residue block. Which branch applies is decided by looking for
    a ``>`` header rather than by trial and error, so a headerless reference does
    not provoke parser warnings about "comments before the first sequence".

    Raises:
        ValueError: the file is empty, or is FASTA that cannot be parsed.
    """
    path = Path(path)
    text = path.read_text(encoding="utf-8", errors="ignore")
    lines = [line.strip() for line in text.splitlines() if line.strip()]

    if not lines:
        raise ValueError(f"Sequence file is empty: {path}")

    if lines[0].startswith(">"):
        parse_errors: list[str] = []
        for fmt in ("fasta-pearson", "fasta-2line"):
            try:
                record = SeqIO.read(io.StringIO(text), fmt)
            except Exception as exc:  # noqa: BLE001 - dialect probing, keep going
                parse_errors.append(f"{fmt}: {exc}")
                continue
            cleaned = sanitize_protein_sequence(str(record.seq))
            if not cleaned:
                raise ValueError(f"FASTA record has no residues: {path}")
            record.seq = Seq(cleaned)
            return cast(SeqRecord, record)

        # It claims to be FASTA but no dialect accepted it -- most often two or
        # more records where exactly one is expected. Say so rather than guess.
        raise ValueError(
            f"Could not parse FASTA file {path} (expected exactly one record). "
            f"Attempts: {'; '.join(parse_errors)}"
        )

    sequence = sanitize_protein_sequence("".join(lines))
    if not sequence:
        raise ValueError(f"No amino-acid residues detected in {path}")

    if not _PROTEIN_ALPHABET.fullmatch(sequence):
        offending = sorted(set(sequence) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO*"))
        logger.warning(
            "Non-standard residue symbols %s in %s; proceeding with best effort",
            "".join(offending),
            path,
        )

    logger.debug("Loaded headerless sequence file: %s (%d aa)", path, len(sequence))
    return SeqRecord(
        Seq(sequence),
        id=fallback_id,
        description=f"headerless_source:{path.name}",
    )


def write_protein_record(
    path: Path | str,
    sequence: str,
    seq_id: str,
    description: str = "",
) -> Path:
    """Write one protein as wrapped, well-formed FASTA. Returns the path."""
    path = Path(path)
    cleaned = sanitize_protein_sequence(sequence)
    if not cleaned:
        raise ValueError(f"Refusing to write an empty sequence to {path}")

    path.parent.mkdir(parents=True, exist_ok=True)
    record = SeqRecord(Seq(cleaned), id=seq_id, description=description)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        SeqIO.write(record, handle, "fasta")
    return path
