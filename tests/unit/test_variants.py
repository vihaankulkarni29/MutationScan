"""Tests for substitution calling.

The residue-counter numbering and the two quality gates are the scientific core
of the tool, so they are tested directly on synthetic sequences where the right
answer is known by construction.
"""

from __future__ import annotations

import pandas as pd
import pytest

from mutation_scan.extract import STATUS_EXTRACTED, STATUS_NO_HIT, protein_path
from mutation_scan.fasta import write_protein_record
from mutation_scan.variants import (
    MUTATION_COLUMNS,
    QC_COLUMNS,
    VARIANT_STATUS_CALLED,
    VARIANT_STATUS_ERROR,
    VARIANT_STATUS_LOW_COVERAGE,
    VARIANT_STATUS_LOW_IDENTITY,
    VARIANT_STATUS_MISSING_PROTEIN,
    VARIANT_STATUS_NO_REFERENCE,
    build_aligner,
    call_variants,
    compare_sequences,
    read_mutations,
    read_variant_qc,
)

from ..conftest import PROTEIN_ALPHA, PROTEIN_BETA, mutate


# -- the residue counter ----------------------------------------------------


def test_single_substitution_is_numbered_in_reference_coordinates():
    comparison = compare_sequences("MKTIALSY", "MKTIALSL")

    assert comparison.mutations == (("Y", 8, "L"),)
    assert comparison.identity_pct == 87.5
    assert comparison.coverage_pct == 100.0


def test_multiple_substitutions_are_all_reported():
    comparison = compare_sequences("MKTIALSY", "MRTVALSY")

    assert comparison.mutations == (("K", 2, "R"), ("I", 4, "V"))
    assert comparison.identity_pct == 75.0


def test_identical_sequences_yield_no_mutations():
    comparison = compare_sequences(PROTEIN_ALPHA, PROTEIN_ALPHA)

    assert comparison.mutations == ()
    assert comparison.identity_pct == 100.0
    assert comparison.coverage_pct == 100.0


def test_positions_skip_reference_gaps():
    """An insertion in the query must not shift downstream numbering.

    Reference  M K T - - I A L S Y   (positions 1..8)
    Query      M K T V V I A L S L
    The inserted VV occupies no reference position, so the final change is
    reported at position 8, not 10.
    """
    reference = PROTEIN_ALPHA
    query = PROTEIN_ALPHA[:20] + "WW" + mutate(PROTEIN_ALPHA[20:], (5, "P"))

    comparison = compare_sequences(reference, query)

    assert comparison.inserted_columns == 2
    positions = [position for _ref, position, _alt in comparison.mutations]
    assert 25 in positions  # 20 + 5, in reference coordinates
    assert max(positions) <= len(reference)


def test_a_deletion_is_counted_but_not_called_as_a_mutation():
    reference = PROTEIN_ALPHA
    query = PROTEIN_ALPHA[:30] + PROTEIN_ALPHA[35:]  # 5 residues removed

    comparison = compare_sequences(reference, query)

    assert comparison.deleted_positions == 5
    assert comparison.mutations == ()
    assert comparison.coverage_pct < 100.0


def test_terminal_truncation_lowers_coverage_not_identity_of_the_aligned_part():
    """Coverage is the metric that catches a truncated gene."""
    reference = PROTEIN_ALPHA
    query = PROTEIN_ALPHA[: len(PROTEIN_ALPHA) // 2]

    comparison = compare_sequences(reference, query)

    assert comparison.mutations == ()
    assert comparison.coverage_pct == pytest.approx(50.0, abs=1.0)
    # Identity is measured over the reference, so a half-length perfect match
    # scores ~50% -- it does not look like a perfect gene.
    assert comparison.identity_pct == pytest.approx(50.0, abs=1.0)


def test_internal_stop_is_reported_as_a_substitution_to_star():
    reference = "MKTIALSY"
    query = "MKT*ALSY"

    comparison = compare_sequences(reference, query)

    assert comparison.mutations == (("I", 4, "*"),)


def test_unknown_residues_do_not_abort_the_comparison():
    comparison = compare_sequences("MKTIALSY", "MKTUALSY")

    # U is not in BLOSUM62; it is treated as X rather than crashing the run.
    assert comparison.mutations == (("I", 4, "X"),)


def test_empty_sequences_are_rejected():
    with pytest.raises(ValueError, match="reference sequence is empty"):
        compare_sequences("", "MKT")
    with pytest.raises(ValueError, match="query sequence is empty"):
        compare_sequences("MKT", "")


def test_a_shared_aligner_gives_the_same_answer_as_a_fresh_one():
    aligner = build_aligner()
    reference, query = PROTEIN_ALPHA, mutate(PROTEIN_ALPHA, (10, "P"), (40, "W"))

    assert compare_sequences(reference, query, aligner=aligner) == compare_sequences(
        reference, query
    )


# -- call_variants over a directory -----------------------------------------


@pytest.fixture
def run_dirs(tmp_path):
    refs = tmp_path / "refs"
    proteins = tmp_path / "proteins"
    write_protein_record(refs / "alphaX_WT.faa", PROTEIN_ALPHA, "alphaX_WT")
    write_protein_record(refs / "betaQ_WT.faa", PROTEIN_BETA, "betaQ_WT")
    proteins.mkdir(parents=True, exist_ok=True)
    return refs, proteins


def put_protein(proteins_dir, accession, gene, sequence):
    return write_protein_record(
        protein_path(proteins_dir, accession, gene),
        sequence,
        seq_id=f"{accession}_{gene}",
        description=f"gene={gene} accession={accession}",
    )


def manifest_for(pairs, status=STATUS_EXTRACTED):
    return pd.DataFrame(
        [{"Accession": accession, "Gene": gene, "Status": status} for accession, gene in pairs],
        columns=["Accession", "Gene", "Status"],
    )


def test_call_variants_is_driven_by_the_manifest(run_dirs):
    refs, proteins = run_dirs
    put_protein(proteins, "g1", "alphaX", mutate(PROTEIN_ALPHA, (12, "W")))
    put_protein(proteins, "g1", "betaQ", PROTEIN_BETA)

    result = call_variants(
        proteins, refs, manifest=manifest_for([("g1", "alphaX"), ("g1", "betaQ")])
    )

    assert list(result.mutations.columns) == list(MUTATION_COLUMNS)
    assert list(result.qc.columns) == list(QC_COLUMNS)
    assert result.mutations["Mutation"].tolist() == [
        f"{PROTEIN_ALPHA[11]}12W"
    ]
    assert set(result.qc["Status"]) == {VARIANT_STATUS_CALLED}
    # betaQ was evaluated and found unmutated -- a row exists saying so.
    beta = result.qc[result.qc["Gene"] == "betaQ"].iloc[0]
    assert beta["Status"] == VARIANT_STATUS_CALLED
    assert int(beta["N_Mutations"]) == 0


def test_manifest_rows_that_were_not_extracted_are_not_evaluated(run_dirs):
    refs, proteins = run_dirs
    put_protein(proteins, "g1", "alphaX", PROTEIN_ALPHA)
    manifest = pd.DataFrame(
        [
            {"Accession": "g1", "Gene": "alphaX", "Status": STATUS_EXTRACTED},
            {"Accession": "g1", "Gene": "betaQ", "Status": STATUS_NO_HIT},
        ]
    )

    result = call_variants(proteins, refs, manifest=manifest)

    assert result.qc["Gene"].tolist() == ["alphaX"]


def test_identity_gate_rejects_a_diverged_protein(run_dirs):
    refs, proteins = run_dirs
    put_protein(proteins, "g1", "alphaX", PROTEIN_BETA)  # unrelated sequence

    result = call_variants(
        proteins, refs, manifest=manifest_for([("g1", "alphaX")]), min_coverage_percent=0
    )

    assert result.mutations.empty
    row = result.qc.iloc[0]
    assert row["Status"] == VARIANT_STATUS_LOW_IDENTITY
    assert "identity" in row["Detail"]
    assert pd.isna(row["N_Mutations"])


def test_coverage_gate_rejects_a_truncated_protein(run_dirs):
    refs, proteins = run_dirs
    put_protein(proteins, "g1", "alphaX", PROTEIN_ALPHA[:40])

    result = call_variants(
        proteins,
        refs,
        manifest=manifest_for([("g1", "alphaX")]),
        min_identity_percent=0,  # identity would not save us here either
        min_coverage_percent=80,
    )

    row = result.qc.iloc[0]
    assert row["Status"] == VARIANT_STATUS_LOW_COVERAGE
    assert "coverage" in row["Detail"]


def test_lowering_the_gates_lets_a_borderline_pair_through(run_dirs):
    """The old integration test failed because 75% identity met an 80% gate."""
    refs, proteins = run_dirs
    put_protein(proteins, "g1", "alphaX", mutate(PROTEIN_ALPHA, (3, "P"), (9, "P"), (60, "P")))

    strict = call_variants(
        proteins, refs, manifest=manifest_for([("g1", "alphaX")]), min_identity_percent=99.9
    )
    lenient = call_variants(
        proteins, refs, manifest=manifest_for([("g1", "alphaX")]), min_identity_percent=90
    )

    assert strict.qc.iloc[0]["Status"] == VARIANT_STATUS_LOW_IDENTITY
    assert lenient.qc.iloc[0]["Status"] == VARIANT_STATUS_CALLED
    assert len(lenient.mutations) == 3


def test_missing_protein_file_and_missing_reference_are_reported(run_dirs):
    refs, proteins = run_dirs
    put_protein(proteins, "g1", "alphaX", PROTEIN_ALPHA)

    result = call_variants(
        proteins,
        refs,
        manifest=manifest_for([("g1", "alphaX"), ("g1", "betaQ"), ("g1", "gammaZ")]),
    )
    status = dict(zip(result.qc["Gene"], result.qc["Status"]))

    assert status["alphaX"] == VARIANT_STATUS_CALLED
    assert status["betaQ"] == VARIANT_STATUS_MISSING_PROTEIN  # in manifest, not on disk
    assert status["gammaZ"] == VARIANT_STATUS_NO_REFERENCE


def test_an_unreadable_protein_is_recorded_as_an_error(run_dirs):
    refs, proteins = run_dirs
    put_protein(proteins, "g1", "alphaX", PROTEIN_ALPHA)
    bad = protein_path(proteins, "g2", "alphaX")
    bad.write_text(">a\nSEQ\n>b\nSEQ\n", encoding="utf-8")  # two records

    result = call_variants(
        proteins, refs, manifest=manifest_for([("g1", "alphaX"), ("g2", "alphaX")])
    )
    status = dict(zip(result.qc["Accession"], result.qc["Status"]))

    assert status["g1"] == VARIANT_STATUS_CALLED
    assert status["g2"] == VARIANT_STATUS_ERROR


def test_directory_scan_fallback_without_a_manifest(run_dirs):
    refs, proteins = run_dirs
    put_protein(proteins, "562.112393", "alphaX", mutate(PROTEIN_ALPHA, (7, "W")))

    result = call_variants(proteins, refs)

    assert result.mutations["Accession"].tolist() == ["562.112393"]
    assert result.mutations["Gene"].tolist() == ["alphaX"]


def test_directory_scan_uses_the_header_for_genes_containing_underscores(run_dirs):
    """Filename splitting alone would mis-parse `g1_some_gene.faa`."""
    refs, proteins = run_dirs
    write_protein_record(refs / "some_gene_WT.faa", PROTEIN_ALPHA, "some_gene_WT")
    write_protein_record(
        proteins / "g1_some_gene.faa",
        mutate(PROTEIN_ALPHA, (4, "W")),
        seq_id="g1_some_gene",
        description="gene=some_gene accession=g1",
    )

    result = call_variants(proteins, refs)

    assert result.mutations["Gene"].tolist() == ["some_gene"]
    assert result.mutations["Accession"].tolist() == ["g1"]


def test_empty_inputs_produce_well_formed_empty_tables(run_dirs, tmp_path):
    refs, proteins = run_dirs

    result = call_variants(
        proteins,
        refs,
        manifest=manifest_for([]),
        mutations_csv=tmp_path / "mutations.csv",
        qc_csv=tmp_path / "variant_qc.csv",
    )

    assert result.mutations.empty and result.qc.empty
    assert list(read_mutations(tmp_path / "mutations.csv").columns) == list(MUTATION_COLUMNS)
    assert list(read_variant_qc(tmp_path / "variant_qc.csv").columns) == list(QC_COLUMNS)


def test_a_manifest_missing_required_columns_is_rejected(run_dirs):
    refs, proteins = run_dirs

    with pytest.raises(ValueError, match="missing column"):
        call_variants(proteins, refs, manifest=pd.DataFrame([{"a": 1}]))


def test_output_is_sorted_and_round_trips_through_csv(run_dirs, tmp_path):
    refs, proteins = run_dirs
    put_protein(proteins, "g2", "alphaX", mutate(PROTEIN_ALPHA, (50, "W")))
    put_protein(proteins, "g1", "alphaX", mutate(PROTEIN_ALPHA, (30, "W"), (10, "P")))
    mutations_csv = tmp_path / "mutations.csv"

    result = call_variants(
        proteins,
        refs,
        manifest=manifest_for([("g2", "alphaX"), ("g1", "alphaX")]),
        mutations_csv=mutations_csv,
    )

    assert result.mutations["Accession"].tolist() == ["g1", "g1", "g2"]
    assert result.mutations["Position"].tolist() == [10, 30, 50]
    reloaded = read_mutations(mutations_csv)
    assert reloaded["Mutation"].tolist() == result.mutations["Mutation"].tolist()


def test_numeric_accessions_survive_the_csv_round_trip(run_dirs, tmp_path):
    refs, proteins = run_dirs
    put_protein(proteins, "562.112393", "alphaX", mutate(PROTEIN_ALPHA, (5, "W")))
    mutations_csv = tmp_path / "mutations.csv"

    call_variants(
        proteins,
        refs,
        manifest=manifest_for([("562.112393", "alphaX")]),
        mutations_csv=mutations_csv,
    )

    assert read_mutations(mutations_csv)["Accession"].iloc[0] == "562.112393"
