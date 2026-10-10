"""dN/dS estimation from nucleotide CDS.

The pipeline extracts proteins with ``tblastn``, which reports the
*nucleotide* coordinates of each alignment.  This module uses those
coordinates to pull the actual coding sequence from the genome FASTA,
back-translates the *reference* protein with the standard genetic code
to obtain a reference CDS, and compares the two codon-by-codon to count
synonymous (dS) and non-synonymous (dN) nucleotide differences.

Because the reference is provided as a protein (no reference genome FASTA),
the reference CDS is a back-translation using the most common codon per
amino acid.  This is the standard approach when nucleotide references are
unavailable (e.g. SNAP, the ``ape`` ``code`` back-translation) and is
exact for the reference strain used to derive it.  Genomes from other
strains contribute genuine synonymous divergence (codon-usage differences)
to dS, which is the desired population-level signal.
"""

from __future__ import annotations

import functools
import logging
from collections import defaultdict
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

__all__ = [
    "DN_DS_COLUMNS",
    "GENETIC_CODE",
    "MOST_COMMON_CODON",
    "back_translate",
    "codon_sites",
    "compute_dn_ds",
    "extract_cds",
    "read_genome_contig",
    "reverse_complement",
]

#: Standard genetic code (NCBI translation table 1): codon -> amino acid.
GENETIC_CODE: dict[str, str] = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

#: Most common codon per amino acid (E. coli-biased for the default cohort).
MOST_COMMON_CODON: dict[str, str] = {
    "A": "GCT",
    "R": "CGT",
    "N": "AAT",
    "D": "GAT",
    "C": "TGT",
    "Q": "CAA",
    "E": "GAG",
    "G": "GGT",
    "H": "CAT",
    "I": "ATT",
    "L": "CTG",
    "K": "AAA",
    "M": "ATG",
    "F": "TTT",
    "P": "CCG",
    "S": "TCT",
    "T": "ACT",
    "W": "TGG",
    "Y": "TAT",
    "V": "GTT",
    "*": "TAA",
}

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcaN")
_STOP = "*"
_LOF_ALLELE = "LoF"


def back_translate(protein: str) -> str:
    """Convert a protein sequence to nucleotide using the most common codon."""
    return "".join(MOST_COMMON_CODON.get(aa, "GCT") for aa in protein.upper() if aa != _STOP)


def reverse_complement(nucleotide: str) -> str:
    """Reverse-complement a nucleotide sequence."""
    return nucleotide.translate(_COMPLEMENT)[::-1]


def codon_sites(codon: str) -> tuple[float, float]:
    """Count potential synonymous and non-synonymous sites in a codon.

    Uses the Nei-Gojobari method: each of the 9 possible point mutations
    (3 positions x 3 alternative nucleotides) is classified as synonymous
    (same amino acid) or non-synonymous (different amino acid).

    Returns:
        ``(non_syn_sites, syn_sites)``.  The sum is 9 for a standard codon.
    """
    codon = codon.upper()
    aa = GENETIC_CODE.get(codon)
    if aa is None:
        return (0.0, 0.0)

    non_syn = 0.0
    syn = 0.0
    for pos in range(3):
        for nt in "ACGT":
            if nt == codon[pos]:
                continue
            mutant = codon[:pos] + nt + codon[pos + 1 :]
            mutant_aa = GENETIC_CODE.get(mutant)
            if mutant_aa is None:
                non_syn += 1
            elif mutant_aa == aa:
                syn += 1
            else:
                non_syn += 1
    return (non_syn, syn)


@functools.lru_cache(maxsize=512)
def _read_genome(contig: str, genome_path: str) -> str | None:
    """Read one contig from a genome assembly, cached per genome file.

    ``genome_path`` is a string so the lru_cache key is stable.
    """
    path = Path(genome_path)
    if not path.is_file():
        return None
    from Bio import SeqIO

    target = contig.strip()
    for record in SeqIO.parse(path, "fasta"):
        if record.id == target or record.name == target:
            return str(record.seq).upper()
    return None


def read_genome_contig(genome_path: Path | str, contig: str) -> str:
    """Read a single contig from a multi-record genome FASTA.

    Returns an empty string when the contig is not found.
    """
    seq = _read_genome(contig, str(genome_path))
    return seq or ""


def extract_cds(
    genome_path: Path | str,
    contig: str,
    start: int,
    end: int,
) -> str:
    """Extract a nucleotide CDS from a genome using tblastn coordinates.

    tblastn reports 1-based, inclusive coordinates on the subject strand.
    When ``end < start`` the HSP is on the reverse strand and the slice is
    reverse-complemented.
    """
    seq = read_genome_contig(genome_path, contig)
    if not seq or start < 1 or end < 1:
        return ""
    if start <= end:
        return seq[start - 1 : end]
    return reverse_complement(seq[end - 1 : start])


def _translate(nucleotide: str) -> str:
    """Translate a nucleotide CDS to protein (no stop trimming)."""
    if len(nucleotide) % 3 != 0:
        return ""
    parts = []
    for i in range(0, len(nucleotide), 3):
        aa = GENETIC_CODE.get(nucleotide[i : i + 3].upper())
        if aa is None:
            parts.append("X")
        else:
            parts.append(aa)
    return "".join(parts)


def _compare_codons(ref_cds: str, query_cds: str) -> tuple[int, int]:
    """Count non-synonymous (dN) and synonymous (dS) nucleotide differences.

    Each nucleotide position where the two CDS differ is classified by
    whether the encoded amino acid also differs.  CDS sequences must be
    the same length and a multiple of three.
    """
    n_diffs = 0
    s_diffs = 0
    for i in range(0, min(len(ref_cds), len(query_cds)), 3):
        ref_codon = ref_cds[i : i + 3]
        query_codon = query_cds[i : i + 3]
        if len(ref_codon) < 3 or len(query_codon) < 3:
            continue
        ref_aa = GENETIC_CODE.get(ref_codon.upper())
        query_aa = GENETIC_CODE.get(query_codon.upper())
        if ref_aa is None or query_aa is None:
            continue
        diff = sum(1 for a, b in zip(ref_codon, query_codon, strict=True) if a != b)
        if diff == 0:
            continue
        if ref_aa == query_aa:
            s_diffs += diff
        else:
            n_diffs += diff
    return n_diffs, s_diffs


DN_DS_COLUMNS = (
    "Gene",
    "N_Sites",
    "S_Sites",
    "dN",
    "dS",
    "dN_dS",
    "N_Substitutions",
    "N_Genomes_Evaluated",
)


def _gene_dn_ds(
    reference_protein: str,
    genome_cdna: list[tuple[str, str, int]],
) -> dict:
    """Compute dN/dS for one gene.

    Args:
        reference_protein: the wild-type reference AA sequence.
        genome_cdna: list of ``(accession, nucleotide_cds,
            nucleotide_start)`` for genomes that carry this
            gene.  ``nucleotide_start`` is the 1-based
            subject coordinate tblastn reported, and is used
            to trim the reference CDS to the same region the
            genome CDS was extracted from.

    A tblastn hit does not always span the whole reference:
    low-complexity regions (a signal peptide, for instance)
    are masked, so the alignment -- and therefore the
    extracted CDS -- starts partway into the gene. Comparing
    a full-length reference against a partial CDS would
    misalign every codon, so the reference is trimmed to the
    hit's region before counting sites and differences. The
    site counts are accumulated per genome, because different
    genomes of the same gene can have different hit regions.
    """
    if not reference_protein:
        return _empty_dn_ds_result()

    ref_cdna = back_translate(reference_protein)
    if not ref_cdna or len(ref_cdna) % 3 != 0:
        return _empty_dn_ds_result()

    if genome_cdna is None or not genome_cdna:
        return _empty_dn_ds_result()

    dn_count = 0
    ds_count = 0
    substitutions = 0
    n_evaluated = 0
    total_n_sites = 0.0
    total_s_sites = 0.0

    for _accession, genome_cds, nucleotide_start in genome_cdna:
        if not genome_cds or len(genome_cds) % 3 != 0:
            continue
        # Trim the reference to the region the genome CDS was
        # extracted from. The subject coordinate is 1-based;
        # the reference CDS is the back-translation of the
        # whole protein, so the hit at subject position N
        # corresponds to reference CDS offset N-1.
        start = max(0, (nucleotide_start or 1) - 1)
        trimmed_ref = ref_cdna[start : start + len(genome_cds)]
        if len(trimmed_ref) != len(genome_cds):
            continue

        n_sites = 0.0
        s_sites = 0.0
        for i in range(0, len(trimmed_ref), 3):
            n, s = codon_sites(trimmed_ref[i : i + 3])
            n_sites += n
            s_sites += s
        total_n_sites += n_sites
        total_s_sites += s_sites

        n_diffs, s_diffs = _compare_codons(trimmed_ref, genome_cds)
        if n_diffs == 0 and s_diffs == 0:
            n_evaluated += 1
            continue
        dn_count += n_diffs
        ds_count += s_diffs
        substitutions += 1  # at least one codon differed
        n_evaluated += 1

    if substitutions == 0:
        return _empty_dn_ds_result(n_evaluated=n_evaluated)

    dN = dn_count / total_n_sites if total_n_sites > 0 else 0.0
    dS = ds_count / total_s_sites if total_s_sites > 0 else 0.0
    dn_ds = dN / dS if dS > 0 else float("inf")

    return {
        "N_Sites": round(total_n_sites, 1),
        "S_Sites": round(total_s_sites, 1),
        "dN": round(dN, 6),
        "dS": round(dS, 6),
        "dN_dS": dn_ds if dn_ds != float("inf") else float("inf"),
        "N_Substitutions": substitutions,
        "N_Genomes_Evaluated": n_evaluated,
    }


def _empty_dn_ds_result(n_evaluated: int = 0) -> dict:
    return {
        "N_Sites": 0.0,
        "S_Sites": 0.0,
        "dN": 0.0,
        "dS": 0.0,
        "dN_dS": 0.0,
        "N_Substitutions": 0,
        "N_Genomes_Evaluated": n_evaluated,
    }


def compute_dn_ds(
    mutations: pd.DataFrame | None,
    manifest: pd.DataFrame | None,
    refs_dir: Path | str,
    genomes_dir: Path | str,
    genes: list[str] | None = None,
) -> pd.DataFrame:
    """Compute per-gene dN/dS for every gene with at least one substitution.

    Args:
        mutations: the variant stage's mutation table. Must contain at
            least ``Gene`` and ``Mutation`` columns.  LoF alleles
            (``Mutation == "LoF"``) are excluded -- they are not amino-acid
            substitutions.
        manifest: the extraction manifest. Must have ``Accession``,
            ``Gene``, ``Status``, ``Contig``, ``Nucleotide_Start``,
            ``Nucleotide_End`` columns.
        refs_dir: directory of prepared reference proteins (``{gene}_WT.faa``).
        genomes_dir: directory of genome assemblies (one ``*.fna`` per genome).
        genes: explicit gene list to evaluate.  Defaults to the intersection of
            *manifest* and *refs_dir*.

    Returns:
        DataFrame with :data:`DN_DS_COLUMNS`, sorted by descending
        ``dN_dS`` (genes under the strongest positive selection first).
    """
    refs_dir = Path(refs_dir)
    genomes_dir = Path(genomes_dir)

    if mutations is None or mutations.empty or manifest is None or manifest.empty:
        return pd.DataFrame(columns=list(DN_DS_COLUMNS))

    required = ("Accession", "Gene", "Mutation")
    missing = [c for c in required if c not in mutations.columns]
    if missing:
        raise ValueError(f"Mutations table missing column(s): {', '.join(missing)}")

    # Only substitutions, not LoF alleles.
    substitutions = mutations[mutations["Mutation"].astype(str) != _LOF_ALLELE]

    if substitutions.empty:
        return pd.DataFrame(columns=list(DN_DS_COLUMNS))

    if genes is None:
        genes = sorted(substitutions["Gene"].astype(str).unique())

    # Index manifest entries by (accession, gene) for fast CDS extraction.
    extracted = manifest[
        (manifest["Status"].astype(str) == "extracted")
        & manifest["Contig"].notna()
        & manifest["Nucleotide_Start"].notna()
        & manifest["Nucleotide_End"].notna()
    ]

    # Group CDS by gene across genomes. Each entry carries
    # the subject start coordinate so _gene_dn_ds can trim
    # the reference CDS to the region that was actually
    # extracted -- a masked N-terminus (a signal peptide,
    # say) otherwise makes the extracted CDS shorter than
    # the reference and the gene is dropped entirely.
    by_gene: dict[str, list[tuple[str, str, int]]] = defaultdict(list)
    manifest_cols = {"Accession", "Gene", "Contig", "Nucleotide_Start", "Nucleotide_End"}
    if manifest_cols.issubset(extracted.columns):
        for row in extracted.itertuples(index=False):
            accession = str(row.Accession)
            gene = str(row.Gene)
            contig = str(getattr(row, "Contig", ""))
            ns = int(getattr(row, "Nucleotide_Start", 0) or 0)
            ne = int(getattr(row, "Nucleotide_End", 0) or 0)
            if not accession or not gene or not contig or ns == 0 or ne == 0:
                continue
            genome_path = genomes_dir / f"{accession}.fna"
            if not genome_path.is_file():
                genome_path = genomes_dir / f"{accession}.fa"
            if not genome_path.is_file():
                genome_path = genomes_dir / f"{accession}.fasta"
            if not genome_path.is_file():
                continue
            cdna = extract_cds(genome_path, contig, ns, ne)
            if cdna:
                by_gene[gene].append((accession, cdna, ns))

    from .fasta import read_protein_sequence

    rows: list[dict] = []
    for gene in genes:
        gene = str(gene)
        ref_path = refs_dir / f"{gene}_WT.faa"
        if not ref_path.is_file():
            continue

        try:
            record = read_protein_sequence(ref_path, fallback_id=gene)
        except ValueError:
            continue
        reference_protein = str(record.seq)

        result = _gene_dn_ds(reference_protein, by_gene.get(gene, []))
        result["Gene"] = gene
        rows.append(result)

    if not rows:
        return pd.DataFrame(columns=list(DN_DS_COLUMNS))

    frame = pd.DataFrame(rows, columns=list(DN_DS_COLUMNS))
    infinite = frame["dN_dS"] == float("inf")
    finite = frame[~infinite]
    frame = pd.concat(
        [finite.sort_values("dN_dS", ascending=False), frame[infinite]],
        ignore_index=True,
    )
    return frame
