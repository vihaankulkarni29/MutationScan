"""Tests for tolerant protein FASTA reading/writing."""

from __future__ import annotations

import pytest

from mutation_scan.fasta import (
    read_protein_sequence,
    sanitize_protein_sequence,
    write_protein_record,
)

SEQ = "MKTIALSYIRKPQWNDEFGH"


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("MKT IAL\nSY", "MKTIALSY"),
        ("mktialsy", "MKTIALSY"),
        ("MKTIALSY*", "MKTIALSY"),
        ("  10 MKTIAL 20 SY  ", "MKTIALSY"),  # numbered/wrapped legacy layout
        ("", ""),
    ],
)
def test_sanitize_protein_sequence(raw, expected):
    assert sanitize_protein_sequence(raw) == expected


def test_reads_standard_fasta(tmp_path):
    path = tmp_path / "ref.faa"
    path.write_text(f">sp|P00000|ALPHAX description here\n{SEQ}\n", encoding="utf-8")

    record = read_protein_sequence(path)

    assert str(record.seq) == SEQ
    assert record.id == "sp|P00000|ALPHAX"


def test_reads_wrapped_fasta(tmp_path):
    path = tmp_path / "ref.faa"
    path.write_text(">alphaX\nMKTIALSY\nIRKPQWND\nEFGH\n", encoding="utf-8")

    assert str(read_protein_sequence(path).seq) == SEQ


def test_reads_headerless_whitespace_grouped_sequence(tmp_path):
    """The legacy layout real reference directories actually contain."""
    path = tmp_path / "ref.faa"
    path.write_text("MKTIALSY IRKPQWND\nEFGH\n", encoding="utf-8")

    record = read_protein_sequence(path, fallback_id="alphaX")

    assert str(record.seq) == SEQ
    assert record.id == "alphaX"


def test_headerless_read_emits_no_parser_warnings(tmp_path, recwarn):
    path = tmp_path / "ref.faa"
    path.write_text("MKTIALSY IRKPQWND\nEFGH\n", encoding="utf-8")

    read_protein_sequence(path)

    assert [str(w.message) for w in recwarn.list] == []


def test_empty_file_is_rejected(tmp_path):
    path = tmp_path / "ref.faa"
    path.write_text("   \n\n", encoding="utf-8")

    with pytest.raises(ValueError, match="empty"):
        read_protein_sequence(path)


def test_multi_record_fasta_is_rejected_with_a_clear_message(tmp_path):
    path = tmp_path / "ref.faa"
    path.write_text(f">one\n{SEQ}\n>two\n{SEQ}\n", encoding="utf-8")

    with pytest.raises(ValueError, match="exactly one record"):
        read_protein_sequence(path)


def test_header_only_fasta_is_rejected(tmp_path):
    path = tmp_path / "ref.faa"
    path.write_text(">alphaX\n", encoding="utf-8")

    with pytest.raises(ValueError):
        read_protein_sequence(path)


def test_write_then_read_round_trips(tmp_path):
    path = write_protein_record(
        tmp_path / "nested" / "out.faa", "mkt ialsy\nirkpqwndefgh", "alphaX_WT", "local:src.faa"
    )

    text = path.read_text(encoding="utf-8")
    assert text.startswith(">alphaX_WT local:src.faa")
    assert str(read_protein_sequence(path).seq) == SEQ


def test_write_refuses_empty_sequences(tmp_path):
    with pytest.raises(ValueError, match="empty"):
        write_protein_record(tmp_path / "out.faa", "   ", "alphaX")
