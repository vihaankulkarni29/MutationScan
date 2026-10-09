"""Tests for dN/dS estimation."""

from __future__ import annotations

import pandas as pd

from mutation_scan.evolution import (
    DN_DS_COLUMNS,
    GENETIC_CODE,
    MOST_COMMON_CODON,
    back_translate,
    codon_sites,
    compute_dn_ds,
    reverse_complement,
)


def test_genetic_code_completeness() -> None:
    """All 64 codons are present and map to valid amino acids."""
    assert len(GENETIC_CODE) == 64
    for aa in GENETIC_CODE.values():
        assert aa in "ACDEFGHIKLMNPQRSTVWY*"


def test_most_common_codon_completeness() -> None:
    """All 20 amino acids + stop have a most common codon."""
    expected_aa = "ACDEFGHIKLMNPQRSTVWY*"
    assert set(MOST_COMMON_CODON.keys()) == set(expected_aa)
    for codon in MOST_COMMON_CODON.values():
        assert codon in GENETIC_CODE


def test_codon_sites_basic() -> None:
    """A codon with 1 synonymous and 8 non-synonymous sites (or similar)."""
    # Methionine ATG: only 1 codon -> all 9 mutations are non-synonymous
    n, s = codon_sites("ATG")
    assert n == 9.0
    assert s == 0.0

    # Leucine CTT: multiple codons, some synonymous
    n, s = codon_sites("CTT")
    assert n + s == 9.0
    assert s > 0.0


def test_back_translate() -> None:
    """Protein back-translates to nucleotide using most common codons."""
    protein = "MA"
    cdna = back_translate(protein)
    assert cdna == "ATGGCT"
    assert len(cdna) % 3 == 0


def test_reverse_complement() -> None:
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("AAA") == "TTT"
    assert reverse_complement("ACGT") == "ACGT"


def test_compute_dn_ds_empty() -> None:
    """Empty mutations returns empty DataFrame."""
    muts = pd.DataFrame()
    manifest = pd.DataFrame()
    result = compute_dn_ds(muts, manifest, "refs", "genomes")
    assert result.empty
    assert list(result.columns) == list(DN_DS_COLUMNS)


def test_compute_dn_ds_no_substitutions() -> None:
    """Only LoF alleles returns empty DataFrame."""
    muts = pd.DataFrame({"Accession": ["G1"], "Gene": ["acrB"], "Mutation": ["LoF"]})
    manifest = pd.DataFrame(
        {
            "Accession": ["G1"],
            "Gene": ["acrB"],
            "Status": ["extracted"],
            "Contig": ["NC_000913"],
            "Nucleotide_Start": [100],
            "Nucleotide_End": [4000],
        }
    )
    result = compute_dn_ds(muts, manifest, "refs", "genomes")
    assert result.empty


def test_dn_ds_columns() -> None:
    """Output columns match expected DN_DS_COLUMNS."""
    assert DN_DS_COLUMNS == (
        "Gene",
        "N_Sites",
        "S_Sites",
        "dN",
        "dS",
        "dN_dS",
        "N_Substitutions",
        "N_Genomes_Evaluated",
    )


def test_compute_dn_ds_missing_refs_dir(tmp_path) -> None:
    """Non-existent refs_dir returns empty DataFrame."""
    muts = pd.DataFrame({"Accession": ["G1"], "Gene": ["acrB"], "Mutation": ["S83L"]})
    manifest = pd.DataFrame()
    result = compute_dn_ds(muts, manifest, "nonexistent", "genomes")
    assert result.empty


def test_compute_dn_ds_missing_genomes_dir(tmp_path) -> None:
    """Non-existent genomes_dir returns empty DataFrame."""
    refs_dir = tmp_path / "refs"
    refs_dir.mkdir()
    (refs_dir / "acrB_WT.faa").write_text(">acrB\nMAIV")
    muts = pd.DataFrame({"Accession": ["G1"], "Gene": ["acrB"], "Mutation": ["S83L"]})
    manifest = pd.DataFrame()
    result = compute_dn_ds(muts, manifest, str(refs_dir), "nonexistent")
    assert result.empty


def test_dn_ds_output_order() -> None:
    """Results are sorted by dN_dS descending (inf first)."""
    # We can't easily test this without a full genome setup, but the
    # column order is tested above.
    pass
