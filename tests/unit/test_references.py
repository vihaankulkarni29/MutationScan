"""Tests for target discovery and reference preparation.

These tests use invented gene names on purpose: nothing in MutationScan may
depend on a particular gene, organism or naming convention.
"""

from __future__ import annotations

import pytest

from mutation_scan.references import (
    discover_targets,
    gene_from_ref_stem,
    missing_reference_genes,
    prepare_references,
    resolve_reference,
)

SEQ_A = "MKTIALSYIRKPQWNDEFGH"
SEQ_B = "MQRSTLVNGGKKPPWWYYFF"


def write_ref(directory, filename, sequence, header=None):
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / filename
    body = f">{header}\n{sequence}\n" if header else f"{sequence}\n"
    path.write_text(body, encoding="utf-8")
    return path


# -- gene naming ------------------------------------------------------------


@pytest.mark.parametrize(
    ("stem", "expected"),
    [
        ("alphaX", "alphaX"),
        ("alphaX_WT", "alphaX"),
        ("alphaX_wt", "alphaX"),
        ("some_gene_WT", "some_gene"),
    ],
)
def test_gene_from_ref_stem(stem, expected):
    assert gene_from_ref_stem(stem) == expected


# -- resolution -------------------------------------------------------------


@pytest.mark.parametrize(
    "filename",
    ["zetaK_WT.faa", "zetaK.faa", "zetaK_WT.fasta", "zetaK.fa", "zetaK.fas"],
)
def test_resolve_reference_accepts_layouts_and_extensions(tmp_path, filename):
    write_ref(tmp_path, filename, SEQ_A)

    assert resolve_reference(tmp_path, "zetaK") is not None


def test_resolve_reference_is_case_insensitive(tmp_path):
    write_ref(tmp_path, "ZetaK_WT.faa", SEQ_A)

    assert resolve_reference(tmp_path, "zetak") is not None
    assert resolve_reference(tmp_path, "ZETAK") is not None


def test_resolve_reference_ignores_empty_and_foreign_files(tmp_path):
    (tmp_path / "zetaK_WT.faa").write_text("", encoding="utf-8")
    write_ref(tmp_path, "notes.txt", SEQ_A)
    write_ref(tmp_path, "sidecar.json", SEQ_A)

    assert resolve_reference(tmp_path, "zetaK") is None
    assert resolve_reference(tmp_path, "notes") is None


def test_resolve_reference_handles_missing_dir_and_blank_gene(tmp_path):
    assert resolve_reference(tmp_path / "absent", "zetaK") is None
    assert resolve_reference(tmp_path, "  ") is None


# -- discovery --------------------------------------------------------------


def test_discover_targets_reads_the_directory(tmp_path):
    write_ref(tmp_path, "betaQ_WT.faa", SEQ_A)
    write_ref(tmp_path, "alphaX.fasta", SEQ_B)
    write_ref(tmp_path, "README.md", "ignore me")

    assert discover_targets(tmp_path) == ["alphaX", "betaQ"]


def test_discover_targets_deduplicates_gene_variants(tmp_path):
    write_ref(tmp_path, "alphaX.faa", SEQ_A)
    write_ref(tmp_path, "alphaX_WT.faa", SEQ_A)

    assert discover_targets(tmp_path) == ["alphaX"]


def test_explicit_targets_win_and_keep_their_order(tmp_path):
    write_ref(tmp_path, "alphaX.faa", SEQ_A)
    write_ref(tmp_path, "betaQ.faa", SEQ_B)

    assert discover_targets(tmp_path, ["betaQ", "alphaX", "betaq", ""]) == ["betaQ", "alphaX"]


def test_discover_targets_on_empty_dir(tmp_path):
    assert discover_targets(tmp_path) == []


# -- preparation ------------------------------------------------------------


def test_prepare_references_normalizes_into_the_run_dir(tmp_path):
    source = tmp_path / "references"
    # Headerless, whitespace-grouped: the legacy layout found in real refs dirs.
    write_ref(source, "alphaX_WT.faa", "MKTIALSY IRKPQWND\nEFGH")
    refs_out = tmp_path / "run" / "refs"

    result = prepare_references(source, refs_out, ["alphaX"])

    assert result.genes == ["alphaX"]
    assert result.missing == []
    staged = refs_out / "alphaX_WT.faa"
    assert staged.is_file()
    text = staged.read_text(encoding="utf-8")
    assert text.startswith(">alphaX_WT")
    assert "".join(text.splitlines()[1:]) == SEQ_A
    assert result.lengths["alphaX"] == len(SEQ_A)
    assert result.sources["alphaX"].startswith("local:")


def test_prepare_references_canonicalizes_alternate_filenames(tmp_path):
    source = tmp_path / "references"
    write_ref(source, "alphaX.fasta", SEQ_A, header="sp|P0|ALPHAX some description")
    refs_out = tmp_path / "refs"

    result = prepare_references(source, refs_out, ["alphaX"])

    # Regardless of the input name/extension, the run stores {gene}_WT.faa.
    assert (refs_out / "alphaX_WT.faa").is_file()
    assert resolve_reference(refs_out, "alphaX") == refs_out / "alphaX_WT.faa"
    assert result.lengths["alphaX"] == len(SEQ_A)


def test_prepare_references_reports_all_missing_targets(tmp_path):
    source = tmp_path / "references"
    write_ref(source, "alphaX.faa", SEQ_A)
    refs_out = tmp_path / "refs"

    result = prepare_references(source, refs_out, ["alphaX", "betaQ", "gammaZ"])

    assert result.genes == ["alphaX"]
    assert result.missing == ["betaQ", "gammaZ"]
    assert missing_reference_genes(refs_out, ["alphaX", "betaQ"]) == ["betaQ"]


def test_prepare_references_never_touches_the_network_without_a_taxid(tmp_path, monkeypatch):
    def explode(*args, **kwargs):  # pragma: no cover - must not be called
        raise AssertionError("UniProt was contacted without a configured taxid")

    monkeypatch.setattr("mutation_scan.references.fetch_uniprot_reference", explode)

    result = prepare_references(tmp_path / "references", tmp_path / "refs", ["alphaX"])

    assert result.missing == ["alphaX"]


def test_prepare_references_fetches_missing_targets_when_taxid_is_set(tmp_path, monkeypatch):
    calls = []

    def fake_fetch(gene, taxid):
        calls.append((gene, taxid))
        return f">sp|X1|{gene} fetched protein\n{SEQ_B}\n"

    monkeypatch.setattr("mutation_scan.references.fetch_uniprot_reference", fake_fetch)
    monkeypatch.setattr("mutation_scan.references.time.sleep", lambda _seconds: None)
    refs_out = tmp_path / "refs"

    result = prepare_references(tmp_path / "references", refs_out, ["betaQ"], uniprot_taxid="9999")

    assert calls == [("betaQ", "9999")]
    assert result.genes == ["betaQ"]
    assert result.sources["betaQ"] == "uniprot:9999"
    assert "".join((refs_out / "betaQ_WT.faa").read_text(encoding="utf-8").splitlines()[1:]) == SEQ_B


def test_prepare_references_prefers_local_over_fetching(tmp_path, monkeypatch):
    source = tmp_path / "references"
    write_ref(source, "alphaX_WT.faa", SEQ_A)

    def explode(*args, **kwargs):  # pragma: no cover - must not be called
        raise AssertionError("fetched a reference that was already local")

    monkeypatch.setattr("mutation_scan.references.fetch_uniprot_reference", explode)

    result = prepare_references(source, tmp_path / "refs", ["alphaX"], uniprot_taxid="9999")

    assert result.sources["alphaX"].startswith("local:")


def test_prepare_references_marks_a_failed_fetch_as_missing(tmp_path, monkeypatch):
    monkeypatch.setattr("mutation_scan.references.fetch_uniprot_reference", lambda g, t: None)

    result = prepare_references(tmp_path / "references", tmp_path / "refs", ["betaQ"], uniprot_taxid="1")

    assert result.genes == []
    assert result.missing == ["betaQ"]


def test_prepare_references_survives_an_unparsable_reference(tmp_path):
    source = tmp_path / "references"
    source.mkdir(parents=True)
    # Two records where one is expected: unusable, but must not abort the run.
    (source / "alphaX.faa").write_text(f">one\n{SEQ_A}\n>two\n{SEQ_B}\n", encoding="utf-8")
    write_ref(source, "betaQ.faa", SEQ_B)

    result = prepare_references(source, tmp_path / "refs", ["alphaX", "betaQ"])

    assert result.missing == ["alphaX"]
    assert result.genes == ["betaQ"]
