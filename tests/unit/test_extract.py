"""Tests for tblastn-based extraction and the extraction manifest."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from mutation_scan import extract as extract_module
from mutation_scan.extract import (
    MANIFEST_COLUMNS,
    STATUS_ERROR,
    STATUS_EXTRACTED,
    STATUS_NO_HIT,
    STATUS_NO_REFERENCE,
    TblastnError,
    TblastnHit,
    discover_genomes,
    extract,
    protein_path,
    read_manifest,
    verify_tblastn,
)
from mutation_scan.references import prepare_references

from ..conftest import (
    PROTEIN_ALPHA,
    PROTEIN_BETA,
    mutate,
    requires_tblastn,
    synthetic_genome,
)


def write_genome(directory, accession, proteins, suffix=".fna"):
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{accession}{suffix}"
    path.write_text(synthetic_genome(proteins), encoding="utf-8")
    return path


# -- genome discovery -------------------------------------------------------


def test_discover_genomes_accepts_extensions_and_sorts(tmp_path):
    for name in ("b.fna", "a.fasta", "c.fa", "d.fas"):
        (tmp_path / name).write_text(">x\nACGT\n", encoding="utf-8")
    (tmp_path / "notes.txt").write_text("ignore", encoding="utf-8")
    (tmp_path / "empty.fna").write_text("", encoding="utf-8")

    found = [path.name for path in discover_genomes(tmp_path)]

    assert found == ["a.fasta", "b.fna", "c.fa", "d.fas"]


def test_discover_genomes_on_missing_dir(tmp_path):
    assert discover_genomes(tmp_path / "absent") == []


def test_accession_keeps_dotted_names(tmp_path):
    """Genome identifiers often contain dots and must survive intact."""
    path = write_genome(tmp_path, "9999.12345", {"alphaX": PROTEIN_ALPHA})

    assert discover_genomes(tmp_path) == [path]
    assert path.stem == "9999.12345"
    assert protein_path(tmp_path, "9999.12345", "alphaX").name == "9999.12345_alphaX.faa"


# -- binary verification ----------------------------------------------------


def test_verify_tblastn_rejects_a_missing_binary():
    with pytest.raises(TblastnError, match="not found"):
        verify_tblastn("definitely_not_a_real_binary_xyz")


def test_verify_tblastn_rejects_a_nonzero_exit(monkeypatch, tmp_path):
    """The old implementation returned success on a non-zero exit."""
    import subprocess

    fake = tmp_path / "tblastn"
    fake.write_text("", encoding="utf-8")

    monkeypatch.setattr(extract_module.shutil, "which", lambda _name: str(fake))
    monkeypatch.setattr(
        extract_module.subprocess,
        "run",
        lambda *a, **k: subprocess.CompletedProcess(a[0], 1, "", "libfoo.so missing"),
    )

    with pytest.raises(TblastnError, match="exited 1"):
        verify_tblastn(str(fake))


# -- manifest bookkeeping (no BLAST needed) ---------------------------------


@pytest.fixture
def fake_tblastn(monkeypatch):
    """Replace tblastn with a scripted stub keyed by (gene, accession).

    Patches both the single-query and the batched entry point, so these tests
    exercise manifest bookkeeping without needing BLAST installed. A scripted
    *Exception* now fails the whole genome rather than one pair -- see
    ``test_a_failed_call_errors_every_gene_in_that_genome`` for why.
    """
    behaviour: dict = {}

    def _hit(gene: str, sequence: str) -> TblastnHit:
        return TblastnHit(
            sequence=sequence,
            pident=99.0,
            aln_length=len(sequence),
            ref_length=len(sequence),
            ref_coverage_pct=100.0,
            contig=f"contig_{gene}",
            internal_stops=sequence.count("*"),
        )

    def fake_run(ref_faa, genome_fna, binary="tblastn", evalue=None, timeout=None):
        gene = ref_faa.stem.replace("_WT", "")
        outcome = behaviour.get((genome_fna.stem, gene), behaviour.get(gene, "hit"))
        if isinstance(outcome, Exception):
            raise outcome
        if outcome is None:
            return None
        return _hit(gene, outcome if isinstance(outcome, str) else PROTEIN_ALPHA)

    def fake_run_batch(query_faa, genome_fna, genes, binary="tblastn", evalue=None, timeout=None):
        hits = {}
        for gene in genes:
            outcome = behaviour.get((genome_fna.stem, gene), behaviour.get(gene, "hit"))
            if isinstance(outcome, Exception):
                raise outcome
            if outcome is None:
                continue
            hits[gene] = _hit(gene, outcome if isinstance(outcome, str) else PROTEIN_ALPHA)
        return hits

    monkeypatch.setattr(extract_module, "run_tblastn", fake_run)
    monkeypatch.setattr(extract_module, "run_tblastn_batch", fake_run_batch)
    monkeypatch.setattr(extract_module, "verify_tblastn", lambda binary="tblastn": "stub")
    return behaviour


def test_manifest_has_one_row_per_genome_gene_pair(tmp_path, references_dir, fake_tblastn):
    refs = prepare_references(references_dir, tmp_path / "refs", ["alphaX", "betaQ"])
    for accession in ("g1", "g2", "g3"):
        write_genome(tmp_path / "genomes", accession, {"alphaX": PROTEIN_ALPHA})

    manifest = extract(
        tmp_path / "genomes",
        refs.refs_dir,
        tmp_path / "proteins",
        refs,
        manifest_csv=tmp_path / "extraction_manifest.csv",
    )

    assert list(manifest.columns) == list(MANIFEST_COLUMNS)
    assert len(manifest) == 6  # 3 genomes x 2 genes
    assert set(manifest["Accession"]) == {"g1", "g2", "g3"}
    assert set(manifest["Gene"]) == {"alphaX", "betaQ"}
    # Deterministic ordering.
    assert manifest[["Accession", "Gene"]].values.tolist() == [
        ["g1", "alphaX"],
        ["g1", "betaQ"],
        ["g2", "alphaX"],
        ["g2", "betaQ"],
        ["g3", "alphaX"],
        ["g3", "betaQ"],
    ]


def test_no_hit_and_error_are_distinguished_from_success(tmp_path, references_dir, fake_tblastn):
    refs = prepare_references(references_dir, tmp_path / "refs", ["alphaX", "betaQ"])
    write_genome(tmp_path / "genomes", "g1", {"alphaX": PROTEIN_ALPHA})
    write_genome(tmp_path / "genomes", "g2", {"alphaX": PROTEIN_ALPHA})
    fake_tblastn[("g1", "betaQ")] = None  # absent gene
    fake_tblastn[("g2", "betaQ")] = TblastnError("tblastn timed out after 180s")

    manifest = extract(tmp_path / "genomes", refs.refs_dir, tmp_path / "proteins", refs)
    status = {(row.Accession, row.Gene): row.Status for row in manifest.itertuples()}

    assert status[("g1", "alphaX")] == STATUS_EXTRACTED
    assert status[("g1", "betaQ")] == STATUS_NO_HIT
    # Batching changed this deliberately: one tblastn call covers the whole gene
    # panel, so a call that fails takes every gene in that genome with it. What
    # matters is that it is reported as an error rather than silently degrading
    # to no_hit, which would read as "we looked and found nothing".
    assert status[("g2", "alphaX")] == STATUS_ERROR
    assert status[("g2", "betaQ")] == STATUS_ERROR
    detail = manifest.loc[
        (manifest["Accession"] == "g2") & (manifest["Gene"] == "betaQ"), "Detail"
    ].iloc[0]
    assert "timed out" in detail


def test_a_failed_call_errors_every_gene_in_that_genome(tmp_path, references_dir, fake_tblastn):
    """A genome-level failure must not look like a clean empty result."""
    refs = prepare_references(references_dir, tmp_path / "refs", ["alphaX", "betaQ"])
    write_genome(tmp_path / "genomes", "bad", {"alphaX": PROTEIN_ALPHA})
    write_genome(tmp_path / "genomes", "good", {"alphaX": PROTEIN_ALPHA})
    fake_tblastn[("bad", "alphaX")] = TblastnError("tblastn exited 1: bad subject")

    manifest = extract(tmp_path / "genomes", refs.refs_dir, tmp_path / "proteins", refs)
    status = {(row.Accession, row.Gene): row.Status for row in manifest.itertuples()}

    assert status[("bad", "alphaX")] == STATUS_ERROR
    assert status[("bad", "betaQ")] == STATUS_ERROR
    # The healthy genome is untouched, so the failure is scoped, not global.
    assert status[("good", "alphaX")] == STATUS_EXTRACTED
    assert status[("good", "betaQ")] == STATUS_EXTRACTED
    assert not (tmp_path / "proteins" / "bad_alphaX.faa").exists()


def test_target_without_a_reference_is_recorded_not_skipped(tmp_path, references_dir, fake_tblastn):
    refs = prepare_references(references_dir, tmp_path / "refs", ["alphaX"])
    write_genome(tmp_path / "genomes", "g1", {"alphaX": PROTEIN_ALPHA})

    manifest = extract(
        tmp_path / "genomes", refs.refs_dir, tmp_path / "proteins", ["alphaX", "gammaZ"]
    )

    row = manifest[manifest["Gene"] == "gammaZ"].iloc[0]
    assert row["Status"] == STATUS_NO_REFERENCE
    assert len(manifest) == 2


def test_internal_stops_are_counted(tmp_path, references_dir, fake_tblastn):
    refs = prepare_references(references_dir, tmp_path / "refs", ["alphaX"])
    write_genome(tmp_path / "genomes", "g1", {"alphaX": PROTEIN_ALPHA})
    fake_tblastn[("g1", "alphaX")] = "MKTI*ALSY*IRK"

    manifest = extract(tmp_path / "genomes", refs.refs_dir, tmp_path / "proteins", refs)

    assert int(manifest.iloc[0]["Internal_Stops"]) == 2


def test_empty_inputs_yield_an_empty_but_well_formed_manifest(
    tmp_path, references_dir, fake_tblastn
):
    refs = prepare_references(references_dir, tmp_path / "refs", ["alphaX"])
    (tmp_path / "genomes").mkdir()

    manifest = extract(
        tmp_path / "genomes",
        refs.refs_dir,
        tmp_path / "proteins",
        refs,
        manifest_csv=tmp_path / "manifest.csv",
    )

    assert manifest.empty
    assert list(manifest.columns) == list(MANIFEST_COLUMNS)
    assert list(read_manifest(tmp_path / "manifest.csv").columns) == list(MANIFEST_COLUMNS)


def test_threaded_and_serial_extraction_agree(tmp_path, references_dir, fake_tblastn):
    refs = prepare_references(references_dir, tmp_path / "refs", ["alphaX", "betaQ"])
    for index in range(6):
        write_genome(tmp_path / "genomes", f"g{index}", {"alphaX": PROTEIN_ALPHA})

    serial = extract(tmp_path / "genomes", refs.refs_dir, tmp_path / "p1", refs, threads=1)
    threaded = extract(tmp_path / "genomes", refs.refs_dir, tmp_path / "p2", refs, threads=4)

    pd.testing.assert_frame_equal(serial, threaded)


def test_read_manifest_rejects_a_foreign_csv(tmp_path):
    path = tmp_path / "wrong.csv"
    path.write_text("a,b\n1,2\n", encoding="utf-8")

    with pytest.raises(ValueError, match="missing column"):
        read_manifest(path)


def test_read_manifest_keeps_numeric_accessions_as_strings(tmp_path, references_dir, fake_tblastn):
    refs = prepare_references(references_dir, tmp_path / "refs", ["alphaX"])
    write_genome(tmp_path / "genomes", "9999.12345", {"alphaX": PROTEIN_ALPHA})
    manifest_csv = tmp_path / "manifest.csv"

    extract(
        tmp_path / "genomes",
        refs.refs_dir,
        tmp_path / "proteins",
        refs,
        manifest_csv=manifest_csv,
    )

    assert read_manifest(manifest_csv)["Accession"].iloc[0] == "9999.12345"


# -- real tblastn -----------------------------------------------------------


@requires_tblastn
def test_real_tblastn_recovers_the_embedded_protein(tmp_path, references_dir):
    """End-to-end extraction: a protein encoded in synthetic DNA comes back exactly."""
    refs = prepare_references(references_dir, tmp_path / "refs", ["alphaX", "betaQ"])
    variant = mutate(PROTEIN_ALPHA, (5, "W"))
    write_genome(tmp_path / "genomes", "g1", {"alphaX": variant, "betaQ": PROTEIN_BETA})

    manifest = extract(
        tmp_path / "genomes",
        refs.refs_dir,
        tmp_path / "proteins",
        refs,
        manifest_csv=tmp_path / "manifest.csv",
    )

    assert set(manifest["Status"]) == {STATUS_EXTRACTED}

    from mutation_scan.fasta import read_protein_sequence

    extracted = str(read_protein_sequence(protein_path(tmp_path / "proteins", "g1", "alphaX")).seq)
    assert extracted == variant

    row = manifest[manifest["Gene"] == "alphaX"].iloc[0]
    assert row["Ref_Length"] == len(PROTEIN_ALPHA)
    assert row["Ref_Coverage_pct"] == pytest.approx(100.0, abs=0.5)
    assert 90.0 <= float(row["Pident"]) < 100.0  # one substitution
    assert row["Contig"] == "contig_alphaX"


@requires_tblastn
def test_real_tblastn_reports_no_hit_for_an_absent_gene(tmp_path, references_dir):
    refs = prepare_references(references_dir, tmp_path / "refs", ["alphaX", "betaQ"])
    write_genome(tmp_path / "genomes", "g1", {"alphaX": PROTEIN_ALPHA})  # no betaQ

    manifest = extract(tmp_path / "genomes", refs.refs_dir, tmp_path / "proteins", refs)
    status = dict(zip(manifest["Gene"], manifest["Status"], strict=False))

    assert status["alphaX"] == STATUS_EXTRACTED
    assert status["betaQ"] == STATUS_NO_HIT
    assert not protein_path(tmp_path / "proteins", "g1", "betaQ").exists()


@requires_tblastn
def test_real_tblastn_measures_partial_coverage(tmp_path, references_dir):
    """A truncated gene must be reported as low coverage, not silently accepted."""
    refs = prepare_references(references_dir, tmp_path / "refs", ["alphaX"])
    truncated = PROTEIN_ALPHA[: len(PROTEIN_ALPHA) // 2]
    write_genome(tmp_path / "genomes", "g1", {"alphaX": truncated})

    manifest = extract(tmp_path / "genomes", refs.refs_dir, tmp_path / "proteins", refs)
    row = manifest.iloc[0]

    assert row["Status"] == STATUS_EXTRACTED
    assert float(row["Ref_Coverage_pct"]) < 60.0
    # Identity alone would look perfect -- coverage is what catches this.
    assert float(row["Pident"]) > 95.0


# -- query panel construction ------------------------------------------------


def test_build_query_fasta_writes_one_record_per_gene(tmp_path):
    refs = tmp_path / "refs"
    refs.mkdir()
    (refs / "alphaX_WT.faa").write_text(f">orig_id\n{PROTEIN_ALPHA}\n", encoding="utf-8")
    (refs / "betaQ_WT.faa").write_text(f">orig_id\n{PROTEIN_BETA}\n", encoding="utf-8")

    written = extract_module._build_query_fasta(
        {"alphaX": refs / "alphaX_WT.faa", "betaQ": refs / "betaQ_WT.faa", "gone": None},
        tmp_path / "panel.faa",
    )

    assert written == ["alphaX", "betaQ"]
    text = (tmp_path / "panel.faa").read_text(encoding="utf-8")
    # The record id must be the gene, not whatever the reference called itself:
    # run_tblastn_batch maps results back by this id.
    assert text.count(">") == 2
    assert ">alphaX" in text and ">betaQ" in text
    assert "orig_id" not in text


def test_build_query_fasta_rejects_a_multi_record_reference(tmp_path):
    """Two records under one id would become two queries sharing an id."""
    refs = tmp_path / "refs"
    refs.mkdir()
    (refs / "alphaX_WT.faa").write_text(
        f">one\n{PROTEIN_ALPHA}\n>two\n{PROTEIN_BETA}\n", encoding="utf-8"
    )

    with pytest.raises(TblastnError, match="2 records"):
        extract_module._build_query_fasta({"alphaX": refs / "alphaX_WT.faa"}, tmp_path / "p.faa")


def test_build_query_fasta_rejects_an_empty_reference(tmp_path):
    refs = tmp_path / "refs"
    refs.mkdir()
    (refs / "alphaX_WT.faa").write_text(">alphaX\n\n", encoding="utf-8")

    with pytest.raises(TblastnError, match="empty"):
        extract_module._build_query_fasta({"alphaX": refs / "alphaX_WT.faa"}, tmp_path / "p.faa")


def test_build_query_fasta_rejects_a_gene_name_that_cannot_be_a_record_id(tmp_path):
    refs = tmp_path / "refs"
    refs.mkdir()
    (refs / "alphaX_WT.faa").write_text(f">x\n{PROTEIN_ALPHA}\n", encoding="utf-8")

    with pytest.raises(TblastnError, match="FASTA record id"):
        extract_module._build_query_fasta({"bad name": refs / "alphaX_WT.faa"}, tmp_path / "p.faa")


def test_a_multi_record_reference_stops_the_run_early(tmp_path, references_dir, fake_tblastn):
    """The failure must surface before any genome is aligned."""
    refs = prepare_references(references_dir, tmp_path / "refs", ["alphaX"])
    write_genome(tmp_path / "genomes", "g1", {"alphaX": PROTEIN_ALPHA})
    doubled = refs.refs_dir / "alphaX_WT.faa"
    doubled.write_text(doubled.read_text(encoding="utf-8") * 2, encoding="utf-8")

    with pytest.raises(TblastnError, match="2 records"):
        extract(tmp_path / "genomes", refs.refs_dir, tmp_path / "proteins", refs)


# -- batched tblastn --------------------------------------------------------


@requires_tblastn
def test_batched_call_maps_each_query_to_its_own_gene(tmp_path, references_dir):
    """The core risk of batching: results must not be mis-assigned between genes."""
    refs = prepare_references(references_dir, tmp_path / "refs", ["alphaX", "betaQ"])
    variant_a = mutate(PROTEIN_ALPHA, (5, "W"))
    write_genome(tmp_path / "genomes", "g1", {"alphaX": variant_a, "betaQ": PROTEIN_BETA})

    panel = tmp_path / "panel.faa"
    extract_module._build_query_fasta(
        {"alphaX": refs.refs_dir / "alphaX_WT.faa", "betaQ": refs.refs_dir / "betaQ_WT.faa"},
        panel,
    )
    hits = extract_module.run_tblastn_batch(
        panel, tmp_path / "genomes" / "g1.fna", ["alphaX", "betaQ"]
    )

    assert set(hits) == {"alphaX", "betaQ"}
    assert hits["alphaX"].sequence == variant_a
    assert hits["betaQ"].sequence == PROTEIN_BETA
    assert hits["alphaX"].contig == "contig_alphaX"
    assert hits["betaQ"].contig == "contig_betaQ"


@requires_tblastn
def test_max_target_seqs_is_applied_per_query_not_per_call(tmp_path, references_dir):
    """A regression guard for the worst batching failure mode.

    If ``-max_target_seqs 1`` were interpreted per call rather than per query,
    a single call would report one gene and the other 14 would be reported as
    ``no_hit`` -- an almost-empty manifest that still looks like a clean run.
    Placing the two genes on separate contigs makes the two interpretations
    distinguishable.
    """
    refs = prepare_references(references_dir, tmp_path / "refs", ["alphaX", "betaQ"])
    genomes = tmp_path / "genomes"
    write_genome(genomes, "split", {"alphaX": PROTEIN_ALPHA})
    with (genomes / "split.fna").open("a", encoding="utf-8") as handle:
        handle.write(synthetic_genome({"betaQ": PROTEIN_BETA}))

    panel = tmp_path / "panel.faa"
    extract_module._build_query_fasta(
        {"alphaX": refs.refs_dir / "alphaX_WT.faa", "betaQ": refs.refs_dir / "betaQ_WT.faa"},
        panel,
    )
    hits = extract_module.run_tblastn_batch(panel, genomes / "split.fna", ["alphaX", "betaQ"])

    assert set(hits) == {"alphaX", "betaQ"}


@requires_tblastn
def test_batched_and_single_query_agree_on_every_gene(tmp_path, references_dir):
    """Parity against the one-call-per-gene path, on real BLAST output."""
    refs = prepare_references(references_dir, tmp_path / "refs", ["alphaX", "betaQ"])
    variant_a = mutate(PROTEIN_ALPHA, (5, "W"))
    genome = write_genome(tmp_path / "genomes", "g1", {"alphaX": variant_a, "betaQ": PROTEIN_BETA})

    panel = tmp_path / "panel.faa"
    extract_module._build_query_fasta(
        {"alphaX": refs.refs_dir / "alphaX_WT.faa", "betaQ": refs.refs_dir / "betaQ_WT.faa"},
        panel,
    )
    batched = extract_module.run_tblastn_batch(panel, genome, ["alphaX", "betaQ"])

    for gene in ("alphaX", "betaQ"):
        single = extract_module.run_tblastn(refs.refs_dir / f"{gene}_WT.faa", genome)
        assert single is not None
        assert batched[gene] == single, f"{gene} differs between batched and single-query"


@requires_tblastn
def test_a_genome_missing_every_gene_yields_all_no_hit(tmp_path, references_dir):
    refs = prepare_references(references_dir, tmp_path / "refs", ["alphaX", "betaQ"])
    genomes = tmp_path / "genomes"
    genomes.mkdir(parents=True)
    (genomes / "empty.fna").write_text(">c1\nACGTACGTACGTACGTACGT\n", encoding="utf-8")

    manifest = extract(tmp_path / "genomes", refs.refs_dir, tmp_path / "proteins", refs)

    assert len(manifest) == 2
    assert set(manifest["Status"]) == {STATUS_NO_HIT}


def test_run_tblastn_batch_rejects_an_unexpected_query_id(monkeypatch):
    """A query id we never asked for means the mapping assumption is broken."""
    line = "\t".join(
        [
            "ghost",
            "contig",
            "99.0",
            "10",
            "0",
            "0",
            "1",
            "10",
            "1",
            "30",
            "1e-05",
            "50.0",
            "10",
            "20",
            "M" * 10,
        ]
    )
    monkeypatch.setattr(extract_module, "_run_tblastn", lambda command, timeout: line + "\n")

    with pytest.raises(TblastnError, match="unexpected query id"):
        extract_module.run_tblastn_batch(Path("panel.faa"), Path("g.fna"), ["alphaX"])


def test_run_tblastn_batch_keeps_the_highest_scoring_hsp_per_query(monkeypatch):
    """Belt and braces: -max_hsps 1 should make this moot, but not silently lose."""

    def row(query: str, bitscore: float, sequence: str) -> str:
        return "\t".join(
            [
                query,
                "contig",
                "99.0",
                str(len(sequence)),
                "0",
                "0",
                "1",
                str(len(sequence)),
                "1",
                "30",
                "1e-05",
                str(bitscore),
                str(len(sequence)),
                "20",
                sequence,
            ]
        )

    monkeypatch.setattr(
        extract_module,
        "_run_tblastn",
        lambda command, timeout: "\n".join(
            [row("alphaX", 10.0, "MAAAAAAAAA"), row("alphaX", 99.0, "MCCCCCCCCC")]
        ),
    )

    hits = extract_module.run_tblastn_batch(Path("panel.faa"), Path("g.fna"), ["alphaX"])

    assert hits["alphaX"].sequence == "MCCCCCCCCC"
