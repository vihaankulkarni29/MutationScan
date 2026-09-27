"""Shared test helpers.

Fixtures deliberately use invented gene names and synthetic sequences: no test
may depend on a particular organism, gene or reference database, because the
tool itself must not.
"""

from __future__ import annotations

import shutil

import pytest

# One codon per residue -- enough to build a synthetic genome that tblastn can
# translate back to exactly the protein we started from.
_CODONS = {
    "A": "GCT",
    "R": "CGT",
    "N": "AAT",
    "D": "GAT",
    "C": "TGT",
    "Q": "CAA",
    "E": "GAA",
    "G": "GGT",
    "H": "CAT",
    "I": "ATT",
    "L": "CTT",
    "K": "AAA",
    "M": "ATG",
    "F": "TTT",
    "P": "CCT",
    "S": "TCT",
    "T": "ACT",
    "W": "TGG",
    "Y": "TAT",
    "V": "GTT",
    "*": "TAA",
}

#: Two unrelated synthetic proteins, long enough for tblastn to score reliably.
PROTEIN_ALPHA = (
    "MKTIALSYIRKPQWNDEFGHVLMACSTYNQPKRDEIGHFWVLASTMKNCQYPRDEIVGHFLWA"
    "STMKNCQYPRDEIVGHFLWASTMKNCQYPRDEIVGHFLWASTMKNCQYPRDEIVGHFLWAST"
)
PROTEIN_BETA = (
    "MQRSTLVNGGKKPPWWYYFFAACCDDEEHHIILLMMNNQQRRSSTTVVWWYYAAKKPPDDEE"
    "GGHHIILLMMNNQQRRSSTTVVWWYYFFAACCDDEEGGHHIILLMMNNQQRRSSTTVVWWYY"
)

requires_tblastn = pytest.mark.skipif(
    shutil.which("tblastn") is None,
    reason="BLAST+ tblastn is not installed",
)


def back_translate(protein: str) -> str:
    """Turn a protein into DNA using one fixed codon per residue."""
    return "".join(_CODONS[residue] for residue in protein.upper())


def synthetic_genome(
    proteins: dict[str, str],
    flank: str = "ATGCGCTTAAGGCCTTAAGGCATCGATCGGATCCTTAAGGCC" * 3,
) -> str:
    """Build a FASTA genome embedding each protein's coding sequence.

    Each protein goes on its own contig, padded with non-coding flanks so the
    gene is not at a sequence boundary.
    """
    contigs = []
    for name, protein in proteins.items():
        dna = f"{flank}{back_translate(protein)}TAA{flank}"
        contigs.append(f">contig_{name}\n{dna}\n")
    return "".join(contigs)


def mutate(protein: str, *changes: tuple[int, str]) -> str:
    """Apply 1-based (position, new residue) substitutions to a protein."""
    residues = list(protein)
    for position, new in changes:
        assert 1 <= position <= len(residues), f"position {position} outside protein"
        residues[position - 1] = new
    return "".join(residues)


@pytest.fixture
def references_dir(tmp_path):
    """A references directory defining two targets: alphaX and betaQ."""
    directory = tmp_path / "references"
    directory.mkdir()
    (directory / "alphaX_WT.faa").write_text(f">alphaX\n{PROTEIN_ALPHA}\n", encoding="utf-8")
    (directory / "betaQ_WT.faa").write_text(f">betaQ\n{PROTEIN_BETA}\n", encoding="utf-8")
    return directory
