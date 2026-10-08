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
    LOF_ALLELE,
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
    loss_of_function_alleles,
    read_mutations,
    read_variant_qc,
    with_loss_of_function,
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


def test_internal_stop_is_reported_as_a_truncation():
    """A stop codon is a truncation event, not a substitution to ``*``."""
    reference = "MKTIALSY"
    query = "MKT*ALSY"

    comparison = compare_sequences(reference, query)

    assert comparison.mutations == ()
    assert comparison.truncations == ((4, "I", 0.5),)
    assert comparison.functional_knockouts == (4,)


def test_a_terminal_stop_loses_nothing_downstream():
    reference = "MKTIALSY"
    query = "MKTIALS*"

    comparison = compare_sequences(reference, query)

    assert comparison.truncations == ((8, "Y", 0.0),)
    assert comparison.functional_knockouts == ()


def test_a_stop_that_matches_the_reference_is_no_truncation():
    """A reference carrying its own stop cannot be truncated by a match."""
    reference = "MKT*ALSY"
    query = "MKT*ALSY"

    comparison = compare_sequences(reference, query)

    assert comparison.mutations == ()
    assert comparison.truncations == ()


def test_a_truncation_losing_exactly_a_fifth_is_not_a_knockout():
    reference = "MKTAQPLYWN"  # 10 residues; a stop at 8 loses exactly 20%
    query = "MKTAQPL*WN"

    comparison = compare_sequences(reference, query)

    assert comparison.truncations == ((8, "Y", 0.2),)
    assert comparison.functional_knockouts == ()


def test_a_truncation_losing_more_than_a_fifth_is_a_knockout():
    reference = "MKTAQPLYWN"
    query = "MKTAQ*YWN"  # stop at 7: 30% of the coding sequence is lost

    comparison = compare_sequences(reference, query)

    assert comparison.truncations == ((7, "L", 0.3),)
    assert comparison.functional_knockouts == (7,)


def test_the_knockout_threshold_is_a_setting():
    """The same truncation crosses the bar only when the bar is low enough."""
    reference = "MKTAQPLYWN"
    query = "MKTAQ*YWN"  # stop at 7: 30% of the coding sequence is lost

    assert compare_sequences(reference, query).functional_knockouts == (7,)
    # A stricter bar -- past 30% -- leaves it a plain truncation.
    strict = compare_sequences(reference, query, knockout_lost_fraction=0.5)
    assert strict.functional_knockouts == ()
    assert strict.truncations == ((7, "L", 0.3),)


def test_call_variants_takes_the_knockout_threshold(run_dirs):
    refs, proteins = run_dirs
    query = PROTEIN_ALPHA[:40] + PROTEIN_ALPHA[41:]
    query = query[:40] + "*" + query[41:]
    put_protein(proteins, "g1", "alphaX", query)

    lenient = call_variants(
        proteins,
        refs,
        manifest=manifest_for([("g1", "alphaX")]),
        min_identity_percent=0,
        min_coverage_percent=0,
        knockout_lost_percent=10,
    )
    assert lenient.qc.iloc[0]["N_Functional_Knockouts"] == 1

    strict = call_variants(
        proteins,
        refs,
        manifest=manifest_for([("g1", "alphaX")]),
        min_identity_percent=0,
        min_coverage_percent=0,
        knockout_lost_percent=90,
    )
    assert strict.qc.iloc[0]["N_Functional_Knockouts"] == 0
    # The truncation itself is reported either way; the threshold
    # only decides whether it counts as a knockout.
    assert strict.qc.iloc[0]["N_Truncations"] == 1


# A reference with no repeated residues, so the aligner can only
# place an indel where it actually is.
_DISTINCT = "MKTAQPLYWNDCH"


def test_an_out_of_frame_deletion_is_a_frameshift():
    reference = _DISTINCT
    query = _DISTINCT[:6] + _DISTINCT[7:]  # residue 7 deleted

    comparison = compare_sequences(reference, query)

    assert comparison.deleted_positions == 1
    assert comparison.frameshifts == ((7, -1),)


def test_an_in_frame_deletion_is_not_a_frameshift():
    reference = _DISTINCT
    query = _DISTINCT[:6] + _DISTINCT[9:]  # residues 7-9: a whole codon

    comparison = compare_sequences(reference, query)

    assert comparison.deleted_positions == 3
    assert comparison.frameshifts == ()


def test_an_out_of_frame_insertion_is_a_frameshift():
    reference = _DISTINCT
    query = _DISTINCT[:5] + "GG" + _DISTINCT[5:]  # GG inserted after position 5

    comparison = compare_sequences(reference, query)

    assert comparison.inserted_columns == 2
    assert comparison.frameshifts == ((5, 2),)


def test_an_in_frame_insertion_is_not_a_frameshift():
    reference = _DISTINCT
    query = _DISTINCT[:5] + "GGG" + _DISTINCT[5:]

    comparison = compare_sequences(reference, query)

    assert comparison.inserted_columns == 3
    assert comparison.frameshifts == ()


def test_a_terminal_deletion_is_not_a_frameshift():
    """A run at the alignment end is a coverage overhang, not an indel."""
    reference = _DISTINCT
    query = _DISTINCT[:-1]  # last residue missing

    comparison = compare_sequences(reference, query)

    assert comparison.deleted_positions == 1
    assert comparison.frameshifts == ()


def test_a_terminal_insertion_is_not_a_frameshift():
    reference = _DISTINCT
    query = _DISTINCT + "GG"  # extends past the reference end

    comparison = compare_sequences(reference, query)

    assert comparison.inserted_columns == 2
    assert comparison.frameshifts == ()


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
    assert result.mutations["Mutation"].tolist() == [f"{PROTEIN_ALPHA[11]}12W"]
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


def test_called_pairs_report_frameshift_and_truncation_events(run_dirs):
    """Events are part of the QC row for a pair that passes the gates."""
    refs, proteins = run_dirs
    # One residue deleted (out of frame) and a stop codon near the
    # middle of the reference, past the knockout threshold.
    query = PROTEIN_ALPHA[:40] + PROTEIN_ALPHA[41:]
    query = query[:40] + "*" + query[41:]
    put_protein(proteins, "g1", "alphaX", query)

    result = call_variants(
        proteins,
        refs,
        manifest=manifest_for([("g1", "alphaX")]),
        min_identity_percent=0,
        min_coverage_percent=0,
    )

    row = result.qc.iloc[0]
    assert row["Status"] == VARIANT_STATUS_CALLED
    assert row["N_Frameshifts"] == 1
    assert row["N_Truncations"] == 1
    assert row["N_Functional_Knockouts"] == 1
    assert "frameshift fs" in row["Detail"]
    assert "truncation" in row["Detail"]
    assert "knockout" in row["Detail"]


def test_events_are_not_substitutions_and_leave_the_mutation_table(run_dirs):
    refs, proteins = run_dirs
    query = PROTEIN_ALPHA[:40] + PROTEIN_ALPHA[41:]
    query = query[:40] + "*" + query[41:]
    put_protein(proteins, "g1", "alphaX", query)

    result = call_variants(
        proteins,
        refs,
        manifest=manifest_for([("g1", "alphaX")]),
        min_identity_percent=0,
        min_coverage_percent=0,
    )

    assert result.mutations.empty
    assert "truncation" not in result.mutations.columns


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
    status = dict(zip(result.qc["Gene"], result.qc["Status"], strict=False))

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
    status = dict(zip(result.qc["Accession"], result.qc["Status"], strict=False))

    assert status["g1"] == VARIANT_STATUS_CALLED
    assert status["g2"] == VARIANT_STATUS_ERROR


def test_directory_scan_fallback_without_a_manifest(run_dirs):
    refs, proteins = run_dirs
    put_protein(proteins, "9999.12345", "alphaX", mutate(PROTEIN_ALPHA, (7, "W")))

    result = call_variants(proteins, refs)

    assert result.mutations["Accession"].tolist() == ["9999.12345"]
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
    put_protein(proteins, "9999.12345", "alphaX", mutate(PROTEIN_ALPHA, (5, "W")))
    mutations_csv = tmp_path / "mutations.csv"

    call_variants(
        proteins,
        refs,
        manifest=manifest_for([("9999.12345", "alphaX")]),
        mutations_csv=mutations_csv,
    )

    assert read_mutations(mutations_csv)["Accession"].iloc[0] == "9999.12345"


# -- parallel variant calling ------------------------------------------------


def _mixed_cohort(proteins_dir):
    """Enough (genome, gene) pairs, of mixed quality, to exercise chunking.

    Deliberately includes a pair that must be rejected by the identity gate and
    one rejected by the coverage gate, so a parallel run cannot agree with a
    sequential one merely by getting every pair right.
    """
    pairs = []
    for index in range(24):
        accession = f"g{index:02d}"
        if index % 7 == 0:
            # Diverged past the identity gate.
            sequence = "".join("W" if i % 3 else aa for i, aa in enumerate(PROTEIN_ALPHA))
            put_protein(proteins_dir, accession, "alphaX", sequence)
            pairs.append((accession, "alphaX"))
            continue
        if index % 7 == 1:
            # Truncated past the coverage gate.
            put_protein(proteins_dir, accession, "betaQ", PROTEIN_BETA[: len(PROTEIN_BETA) // 3])
            pairs.append((accession, "betaQ"))
            continue
        put_protein(proteins_dir, accession, "alphaX", mutate(PROTEIN_ALPHA, (12, "W")))
        put_protein(proteins_dir, accession, "betaQ", PROTEIN_BETA)
        pairs.extend([(accession, "alphaX"), (accession, "betaQ")])
    return pairs


def test_parallel_variant_calling_matches_the_sequential_result(run_dirs):
    """The parallel path must be a scheduling change and nothing else.

    Identical inputs have to produce identical mutations *and* identical QC
    rows, including which pairs the gates reject. Anything less means the
    parallel path is silently changing results, not just speeding them up.
    """
    refs, proteins = run_dirs
    pairs = _mixed_cohort(proteins)
    manifest = manifest_for(pairs)

    sequential = call_variants(proteins, refs, manifest=manifest, threads=1)
    parallel = call_variants(proteins, refs, manifest=manifest, threads=3)

    assert len(parallel.qc) == len(pairs)
    pd.testing.assert_frame_equal(sequential.mutations, parallel.mutations)
    pd.testing.assert_frame_equal(sequential.qc, parallel.qc)


def test_parallel_variant_calling_reports_every_gate_outcome(run_dirs):
    """Chunking must not drop or merge the rejected pairs.

    A worker crash previously turned every pair in a chunk into an ``error``
    row that looked like a data problem; this pins the real status mix so a
    regression there is visible.
    """
    refs, proteins = run_dirs
    pairs = _mixed_cohort(proteins)

    result = call_variants(proteins, refs, manifest=manifest_for(pairs), threads=3)

    statuses = set(result.qc["Status"])
    assert VARIANT_STATUS_ERROR not in statuses
    assert VARIANT_STATUS_CALLED in statuses
    assert VARIANT_STATUS_LOW_IDENTITY in statuses
    assert result.qc["Accession"].nunique() == len({a for a, _ in pairs})


def test_variant_threads_below_one_is_clamped_not_rejected(run_dirs):
    """``threads=0`` must not silently disable the stage."""
    refs, proteins = run_dirs
    put_protein(proteins, "g1", "alphaX", mutate(PROTEIN_ALPHA, (12, "W")))

    result = call_variants(proteins, refs, manifest=manifest_for([("g1", "alphaX")]), threads=0)

    assert result.mutations["Mutation"].tolist() == [f"{PROTEIN_ALPHA[11]}12W"]


# -- loss-of-function collapsing ---------------------------------------------


def _qc_row(accession: str, gene: str, status: str, **events: int) -> dict:
    return {
        "Accession": accession,
        "Gene": gene,
        "Species": "species",
        "Status": status,
        "Identity_pct": 95.0,
        "Coverage_pct": 95.0,
        "Ref_Length": 100,
        "Query_Length": 100,
        "N_Mutations": 0,
        "N_Deleted_Positions": 0,
        "N_Inserted_Columns": 0,
        "N_Frameshifts": events.get("N_Frameshifts", 0),
        "N_Truncations": events.get("N_Truncations", 0),
        "N_Functional_Knockouts": events.get("N_Functional_Knockouts", 0),
        "Detail": events.get("Detail", ""),
    }


def test_a_functional_knockout_collapses_to_one_lof_allele():
    qc = pd.DataFrame(
        [
            _qc_row("g1", "acrB", VARIANT_STATUS_CALLED, N_Functional_Knockouts=1),
        ]
    )
    alleles = loss_of_function_alleles(qc)

    assert len(alleles) == 1
    row = alleles.iloc[0]
    assert row["Accession"] == "g1"
    assert row["Gene"] == "acrB"
    assert row["Mutation"] == LOF_ALLELE


def test_a_frameshift_collapses_to_one_lof_allele():
    qc = pd.DataFrame([_qc_row("g1", "acrB", VARIANT_STATUS_CALLED, N_Frameshifts=1)])
    alleles = loss_of_function_alleles(qc)

    assert len(alleles) == 1
    assert alleles.iloc[0]["Mutation"] == LOF_ALLELE


def test_a_truncation_below_the_knockout_bar_is_not_a_lof_allele():
    """Most of the protein survives, so it is a QC event only, not a LoF."""
    qc = pd.DataFrame(
        [
            _qc_row(
                "g1",
                "acrB",
                VARIANT_STATUS_CALLED,
                N_Truncations=1,
                Detail="truncation at 90 position 90/100 (90.0% lost)",
            )
        ]
    )
    alleles = loss_of_function_alleles(qc)

    assert alleles.empty


def test_one_allele_per_gene_even_when_multiple_events():
    """Frameshift and knockout on the same pair collapse to a single allele."""
    qc = pd.DataFrame(
        [
            _qc_row(
                "g1",
                "acrB",
                VARIANT_STATUS_CALLED,
                N_Frameshifts=1,
                N_Functional_Knockouts=1,
                Detail="frameshift fs@42; truncation at 7 (30.0% lost)",
            )
        ]
    )
    alleles = loss_of_function_alleles(qc)

    assert len(alleles) == 1
    row = alleles.iloc[0]
    assert row["Mutation"] == LOF_ALLELE
    assert "frameshift" in row["LoF_Events"]
    assert "truncation" in row["LoF_Events"]


def test_rejected_pairs_do_not_contribute_lof_alleles():
    qc = pd.DataFrame(
        [
            _qc_row("g1", "acrB", "low_coverage", N_Frameshifts=1),
            _qc_row("g2", "acrB", VARIANT_STATUS_CALLED, N_Functional_Knockouts=1),
        ]
    )
    alleles = loss_of_function_alleles(qc)

    assert len(alleles) == 1
    assert alleles.iloc[0]["Accession"] == "g2"


def test_qc_without_event_columns_yields_no_alleles():
    qc = pd.DataFrame(
        [
            {"Accession": "g1", "Gene": "acrB", "Status": VARIANT_STATUS_CALLED},
        ]
    )
    assert loss_of_function_alleles(qc).empty


def test_with_loss_of_function_is_idempotent():
    mutations = pd.DataFrame(
        [{"Accession": "g1", "Gene": "acrB", "Mutation": "H596N", "Position": 596}]
    )
    qc = pd.DataFrame([_qc_row("g1", "acrB", VARIANT_STATUS_CALLED, N_Functional_Knockouts=1)])

    once = with_loss_of_function(mutations, qc)
    twice = with_loss_of_function(once, qc)

    assert pd.testing.assert_frame_equal(once, twice) is None
    assert (once["Mutation"] == LOF_ALLELE).sum() == 1


def test_lof_alleles_carry_the_mutation_schema_for_concat():
    qc = pd.DataFrame([_qc_row("g1", "alphaX", VARIANT_STATUS_CALLED, N_Frameshifts=1)])
    alleles = loss_of_function_alleles(qc)

    assert list(alleles.columns) == ["Accession", "Gene", "Mutation", "LoF_Events"]
    merged = with_loss_of_function(pd.DataFrame(columns=list(MUTATION_COLUMNS)), qc)
    assert LOF_ALLELE in merged["Mutation"].astype(str).values


def test_lof_event_drops_through_call_variants_into_qc(run_dirs):
    """End-to-end: a frameshift in a passing pair appears in the QC table."""
    refs, proteins = run_dirs
    query = PROTEIN_ALPHA[:40] + PROTEIN_ALPHA[41:]  # one residue deleted (out of frame)
    put_protein(proteins, "g1", "alphaX", query)

    result = call_variants(
        proteins,
        refs,
        manifest=manifest_for([("g1", "alphaX")]),
        min_identity_percent=0,
        min_coverage_percent=0,
    )

    assert result.qc.iloc[0]["Status"] == VARIANT_STATUS_CALLED
    assert result.qc.iloc[0]["N_Frameshifts"] == 1
    assert not loss_of_function_alleles(result.qc).empty


def test_lof_carriers_appear_in_the_genome_map(run_dirs):
    """A genome whose only mutation is a LoF still shows one row in the map."""
    refs, proteins = run_dirs
    query = PROTEIN_ALPHA[:40] + PROTEIN_ALPHA[41:] + "*" + PROTEIN_ALPHA[42:]
    put_protein(proteins, "g1", "alphaX", query)

    result = call_variants(
        proteins,
        refs,
        manifest=manifest_for([("g1", "alphaX"), ("g1", "betaQ")]),
        min_identity_percent=0,
        min_coverage_percent=0,
    )

    from mutation_scan.epistasis import genome_mutation_map

    frame = genome_mutation_map(result.mutations, result.qc)
    g1 = frame[frame["Accession"] == "g1"].iloc[0]
    assert g1["N_LoF"] >= 1
    assert g1["N_Mutated_Genes"] >= 1
