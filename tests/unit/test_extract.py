"""Tests for tblastn-based extraction and the extraction manifest."""

from __future__ import annotations

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
    """Replace tblastn with a scripted stub keyed by (gene, accession)."""
    behaviour: dict = {}

    def fake_run(ref_faa, genome_fna, binary="tblastn", evalue=None, timeout=None):
        gene = ref_faa.stem.replace("_WT", "")
        key = (genome_fna.stem, gene)
        outcome = behaviour.get(key, behaviour.get(gene, "hit"))
        if isinstance(outcome, Exception):
            raise outcome
        if outcome is None:
            return None
        sequence = outcome if isinstance(outcome, str) else PROTEIN_ALPHA
        return TblastnHit(
            sequence=sequence,
            pident=99.0,
            aln_length=len(sequence),
            ref_length=len(sequence),
            ref_coverage_pct=100.0,
            contig=f"contig_{gene}",
            internal_stops=sequence.count("*"),
        )

    monkeypatch.setattr(extract_module, "run_tblastn", fake_run)
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
    assert status[("g2", "betaQ")] == STATUS_ERROR
    detail = manifest.loc[
        (manifest["Accession"] == "g2") & (manifest["Gene"] == "betaQ"), "Detail"
    ].iloc[0]
    assert "timed out" in detail


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
