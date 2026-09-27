"""Tests for descriptive co-occurrence counting.

The denominator is the part that can be silently wrong, so most of these tests
are about eligibility: a gene that was never evaluated in a genome must not
appear in that genome's counts, in either direction.
"""

from __future__ import annotations

import pandas as pd
import pytest

from mutation_scan.cooccurrence import (
    GENE_COLUMNS,
    VARIANT_COLUMNS,
    cooccurrence,
    eligibility_from_qc,
    read_cooccurrence,
)
from mutation_scan.variants import (
    VARIANT_STATUS_CALLED,
    VARIANT_STATUS_LOW_COVERAGE,
    VARIANT_STATUS_LOW_IDENTITY,
)


def mutations_frame(rows):
    """rows: (accession, gene, mutation) triples."""
    return pd.DataFrame(
        [
            {"Accession": accession, "Gene": gene, "Mutation": mutation}
            for accession, gene, mutation in rows
        ],
        columns=["Accession", "Gene", "Mutation"],
    )


def qc_frame(rows, status=VARIANT_STATUS_CALLED):
    """rows: (accession, gene) pairs that were evaluated, or (accession, gene, status)."""
    records = []
    for row in rows:
        accession, gene = row[0], row[1]
        records.append(
            {
                "Accession": accession,
                "Gene": gene,
                "Status": row[2] if len(row) > 2 else status,
            }
        )
    return pd.DataFrame(records, columns=["Accession", "Gene", "Status"])


# -- eligibility ------------------------------------------------------------


def test_eligibility_counts_only_called_pairs():
    qc = qc_frame(
        [
            ("g1", "alphaX"),
            ("g1", "betaQ", VARIANT_STATUS_LOW_IDENTITY),
            ("g2", "alphaX"),
            ("g2", "betaQ"),
        ]
    )

    assert eligibility_from_qc(qc) == {"alphaX": {"g1", "g2"}, "betaQ": {"g2"}}


def test_eligibility_of_an_empty_or_absent_table():
    assert eligibility_from_qc(pd.DataFrame()) == {}
    assert eligibility_from_qc(None) == {}


def test_eligibility_rejects_a_foreign_table():
    with pytest.raises(ValueError, match="missing column"):
        eligibility_from_qc(pd.DataFrame([{"a": 1}]))


# -- gene-level counting ----------------------------------------------------


def test_counts_are_plain_intersections():
    #        alphaX  betaQ
    # g1       x      x
    # g2       x      -
    # g3       -      x
    # g4       -      -
    mutations = mutations_frame(
        [
            ("g1", "alphaX", "A1V"),
            ("g1", "betaQ", "C2D"),
            ("g2", "alphaX", "A1V"),
            ("g3", "betaQ", "C2D"),
        ]
    )
    qc = qc_frame([(g, gene) for g in ("g1", "g2", "g3", "g4") for gene in ("alphaX", "betaQ")])

    result = cooccurrence(mutations, qc)

    assert list(result.genes.columns) == list(GENE_COLUMNS)
    row = result.genes.iloc[0]
    assert (row["Gene_A"], row["Gene_B"]) == ("alphaX", "betaQ")
    assert int(row["N_Both"]) == 1
    assert int(row["N_A"]) == 2
    assert int(row["N_B"]) == 2
    assert int(row["N_Eligible"]) == 4
    assert row["Frac_B_given_A"] == 0.5
    assert row["Frac_A_given_B"] == 0.5


def test_a_gene_never_evaluated_in_a_genome_is_excluded_from_that_denominator():
    """The core honesty guarantee: betaQ is only evaluated in g1."""
    mutations = mutations_frame(
        [("g1", "alphaX", "A1V"), ("g1", "betaQ", "C2D"), ("g2", "alphaX", "A1V")]
    )
    qc = qc_frame([("g1", "alphaX"), ("g1", "betaQ"), ("g2", "alphaX")])

    row = cooccurrence(mutations, qc).genes.iloc[0]

    # g2 had no betaQ result at all, so it is not in the pair's universe.
    assert int(row["N_Eligible"]) == 1
    assert int(row["N_A"]) == 1  # alphaX mutated in g1 only, within that universe
    assert int(row["N_Both"]) == 1
    assert row["Frac_B_given_A"] == 1.0


def test_a_gate_rejection_removes_a_genome_from_the_denominator():
    mutations = mutations_frame([("g1", "alphaX", "A1V"), ("g1", "betaQ", "C2D")])
    qc = qc_frame(
        [
            ("g1", "alphaX"),
            ("g1", "betaQ"),
            ("g2", "alphaX"),
            ("g2", "betaQ", VARIANT_STATUS_LOW_COVERAGE),
        ]
    )

    row = cooccurrence(mutations, qc).genes.iloc[0]

    assert int(row["N_Eligible"]) == 1  # only g1 had both genes evaluated


def test_a_pair_with_no_shared_genome_is_omitted_entirely():
    """Not zero -- absent. A fraction over an empty universe is meaningless."""
    mutations = mutations_frame([("g1", "alphaX", "A1V"), ("g2", "betaQ", "C2D")])
    qc = qc_frame([("g1", "alphaX"), ("g2", "betaQ")])

    result = cooccurrence(mutations, qc, min_count=0)

    assert result.genes.empty
    assert result.genes_matrix.loc["alphaX", "betaQ"] == 0


def test_unmutated_but_evaluated_genes_still_appear_in_the_matrix():
    mutations = mutations_frame([("g1", "alphaX", "A1V")])
    qc = qc_frame([("g1", "alphaX"), ("g1", "betaQ"), ("g2", "alphaX"), ("g2", "betaQ")])

    result = cooccurrence(mutations, qc, min_count=0)

    assert list(result.genes_matrix.index) == ["alphaX", "betaQ"]
    assert result.gene_totals == {"alphaX": 1, "betaQ": 0}
    assert result.eligible_per_gene == {"alphaX": 2, "betaQ": 2}
    row = result.genes.iloc[0]
    assert int(row["N_Both"]) == 0
    assert row["Frac_A_given_B"] == 0.0  # not NaN, not a division by zero


def test_multiple_mutations_in_one_gene_count_the_genome_once():
    mutations = mutations_frame(
        [
            ("g1", "alphaX", "A1V"),
            ("g1", "alphaX", "K9R"),
            ("g1", "alphaX", "T20P"),
            ("g1", "betaQ", "C2D"),
        ]
    )
    qc = qc_frame([("g1", "alphaX"), ("g1", "betaQ")])

    row = cooccurrence(mutations, qc).genes.iloc[0]

    assert int(row["N_A"]) == 1
    assert int(row["N_Both"]) == 1


def test_min_count_filters_the_long_table_but_not_the_matrix():
    mutations = mutations_frame(
        [("g1", "alphaX", "A1V"), ("g1", "betaQ", "C2D"), ("g2", "alphaX", "A1V")]
    )
    qc = qc_frame([(g, gene) for g in ("g1", "g2") for gene in ("alphaX", "betaQ")])

    assert len(cooccurrence(mutations, qc, min_count=1).genes) == 1
    strict = cooccurrence(mutations, qc, min_count=2)
    assert strict.genes.empty
    assert strict.genes_matrix.loc["alphaX", "betaQ"] == 1


def test_pairs_are_unordered_and_sorted_by_count():
    mutations = mutations_frame(
        [
            ("g1", "alphaX", "A1V"),
            ("g1", "betaQ", "C2D"),
            ("g1", "gammaZ", "E3K"),
            ("g2", "alphaX", "A1V"),
            ("g2", "gammaZ", "E3K"),
            ("g3", "alphaX", "A1V"),
            ("g3", "gammaZ", "E3K"),
        ]
    )
    qc = qc_frame([(g, gene) for g in ("g1", "g2", "g3") for gene in ("alphaX", "betaQ", "gammaZ")])

    genes = cooccurrence(mutations, qc).genes

    assert len(genes) == 3  # 3 choose 2, each unordered
    assert genes[["Gene_A", "Gene_B"]].values.tolist()[0] == ["alphaX", "gammaZ"]
    assert genes["N_Both"].tolist() == sorted(genes["N_Both"].tolist(), reverse=True)


def test_matrix_is_symmetric_with_gene_totals_on_the_diagonal():
    mutations = mutations_frame(
        [("g1", "alphaX", "A1V"), ("g1", "betaQ", "C2D"), ("g2", "alphaX", "A1V")]
    )
    qc = qc_frame([(g, gene) for g in ("g1", "g2") for gene in ("alphaX", "betaQ")])

    matrix = cooccurrence(mutations, qc).genes_matrix

    assert matrix.loc["alphaX", "betaQ"] == matrix.loc["betaQ", "alphaX"] == 1
    assert matrix.loc["alphaX", "alphaX"] == 2
    assert matrix.loc["betaQ", "betaQ"] == 1


# -- variant level ----------------------------------------------------------


def test_variant_level_is_off_unless_requested():
    mutations = mutations_frame([("g1", "alphaX", "A1V"), ("g1", "betaQ", "C2D")])
    qc = qc_frame([("g1", "alphaX"), ("g1", "betaQ")])

    assert cooccurrence(mutations, qc).variants is None
    assert cooccurrence(mutations, qc, levels=("gene", "variant")).variants is not None


def test_variant_pairs_are_labelled_gene_colon_mutation():
    mutations = mutations_frame(
        [
            ("g1", "alphaX", "A1V"),
            ("g1", "betaQ", "C2D"),
            ("g2", "alphaX", "A1V"),
            ("g2", "betaQ", "C2D"),
            ("g3", "alphaX", "K9R"),
        ]
    )
    qc = qc_frame([(g, gene) for g in ("g1", "g2", "g3") for gene in ("alphaX", "betaQ")])

    variants = cooccurrence(mutations, qc, levels=("gene", "variant")).variants

    assert list(variants.columns) == list(VARIANT_COLUMNS)
    top = variants.iloc[0]
    assert (top["Variant_A"], top["Variant_B"]) == ("alphaX:A1V", "betaQ:C2D")
    assert int(top["N_Both"]) == 2
    assert int(top["N_Eligible"]) == 3


def test_two_variants_of_the_same_gene_use_that_genes_own_denominator():
    mutations = mutations_frame(
        [("g1", "alphaX", "A1V"), ("g1", "alphaX", "K9R"), ("g2", "alphaX", "A1V")]
    )
    qc = qc_frame([("g1", "alphaX"), ("g2", "alphaX"), ("g3", "alphaX")])

    variants = cooccurrence(mutations, qc, levels=("variant",)).variants
    row = variants.iloc[0]

    assert (row["Gene_A"], row["Gene_B"]) == ("alphaX", "alphaX")
    assert int(row["N_Eligible"]) == 3
    assert int(row["N_Both"]) == 1


# -- edges and I/O ----------------------------------------------------------


def test_no_mutations_gives_empty_but_well_formed_tables(tmp_path):
    qc = qc_frame([("g1", "alphaX"), ("g1", "betaQ")])

    result = cooccurrence(
        mutations_frame([]),
        qc,
        levels=("gene", "variant"),
        min_count=0,
        genes_csv=tmp_path / "genes.csv",
        matrix_csv=tmp_path / "matrix.csv",
        variants_csv=tmp_path / "variants.csv",
    )

    assert list(result.genes.columns) == list(GENE_COLUMNS)
    assert result.variants.empty
    assert int(result.genes.iloc[0]["N_Both"]) == 0
    assert (result.genes_matrix.values == 0).all()
    assert list(read_cooccurrence(tmp_path / "genes.csv").columns) == list(GENE_COLUMNS)


def test_a_single_gene_produces_no_pairs():
    mutations = mutations_frame([("g1", "alphaX", "A1V")])
    qc = qc_frame([("g1", "alphaX")])

    result = cooccurrence(mutations, qc, min_count=0)

    assert result.genes.empty
    assert result.genes_matrix.shape == (1, 1)


def test_missing_qc_falls_back_and_warns(caplog):
    mutations = mutations_frame(
        [("g1", "alphaX", "A1V"), ("g1", "betaQ", "C2D"), ("g2", "alphaX", "A1V")]
    )

    with caplog.at_level("WARNING"):
        result = cooccurrence(mutations)

    assert "overstated" in caplog.text
    # Without QC every mutated genome looks eligible for every gene.
    assert int(result.genes.iloc[0]["N_Eligible"]) == 2


def test_unknown_level_is_rejected():
    with pytest.raises(ValueError, match="Unknown co-occurrence level"):
        cooccurrence(mutations_frame([]), qc_frame([]), levels=("gene", "epistasis"))


def test_foreign_mutations_table_is_rejected():
    with pytest.raises(ValueError, match="missing column"):
        cooccurrence(pd.DataFrame([{"a": 1}]), qc_frame([("g1", "alphaX")]))


def test_numeric_accessions_are_not_coerced():
    mutations = mutations_frame([("9999.12345", "alphaX", "A1V"), ("9999.12345", "betaQ", "C2D")])
    qc = qc_frame([("9999.12345", "alphaX"), ("9999.12345", "betaQ")])

    assert int(cooccurrence(mutations, qc).genes.iloc[0]["N_Both"]) == 1


def test_output_round_trips_through_csv(tmp_path):
    mutations = mutations_frame(
        [("g1", "alphaX", "A1V"), ("g1", "betaQ", "C2D"), ("g2", "alphaX", "A1V")]
    )
    qc = qc_frame([(g, gene) for g in ("g1", "g2") for gene in ("alphaX", "betaQ")])

    result = cooccurrence(
        mutations,
        qc,
        genes_csv=tmp_path / "genes.csv",
        matrix_csv=tmp_path / "matrix.csv",
    )

    reloaded = read_cooccurrence(tmp_path / "genes.csv")
    pd.testing.assert_frame_equal(reloaded, result.genes)
    matrix = pd.read_csv(tmp_path / "matrix.csv", index_col=0)
    assert matrix.loc["alphaX", "betaQ"] == 1


def test_no_inference_columns_are_emitted():
    """A regression guard on the design constraint, not on the arithmetic."""
    mutations = mutations_frame([("g1", "alphaX", "A1V"), ("g1", "betaQ", "C2D")])
    qc = qc_frame([("g1", "alphaX"), ("g1", "betaQ")])

    result = cooccurrence(mutations, qc, levels=("gene", "variant"))
    columns = set(result.genes.columns) | set(result.variants.columns)

    forbidden = ("p_value", "pvalue", "odds", "enrich", "severity", "score", "significan")
    assert not [column for column in columns if any(word in column.lower() for word in forbidden)]
