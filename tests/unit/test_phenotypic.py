"""Tests for the phenotypic concordance module.

The module makes three claims, and each needs its own tests:

1. **Classification is honest** -- a MIC with no breakpoint, or in a
   band the standard leaves undefined, is ``?`` rather than a
   default. A test suite that only ever feeds well-formed inputs will
   pass while the code happily invents susceptible phenotypes.
2. **Loaders tolerate real-world tables** -- headers differ between a
   BV-BRC export and a hand-written one, MICs arrive as ``<=0.25`` and
   ``>64``, and both TSV and CSV dialects have to parse.
3. **Metrics handle degenerate cohorts** -- one genome, no carriers,
   no resistant genomes, no breakpoints. Each of these must return a
   usable result rather than raising or reporting a misleading zero.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from benchmarks.phenotypic import (
    INTERMEDIATE,
    RESISTANT,
    SUSCEPTIBLE,
    UNKNOWN,
    Breakpoint,
    GeneDrugAssociation,
    MicRecord,
    PhenotypeRules,
    burden_correlation,
    classify_mic,
    concordance,
    gene_drug_associations,
    read_breakpoints,
    read_mic_records,
    read_rules,
    write_report,
)


def _mic(
    accession: str,
    drug: str = "ciprofloxacin",
    mic: float = 0.25,
    organism: str = "e. coli",
) -> MicRecord:
    """One MIC record.

    The drug defaults to the ``ciprofloxacin`` fixture's name so a test
    that forgets it still exercises a breakpoint that exists, rather
    than silently testing a ``?`` phenotype.
    """
    assert drug, "an empty drug can never match a breakpoint"
    return MicRecord(accession=accession, drug=drug, mic=mic, organism=organism)


def _table(*breakpoints: Breakpoint) -> dict[tuple[str, str], Breakpoint]:
    return {(bp.organism.lower(), bp.drug.lower()): bp for bp in breakpoints}


@pytest.fixture
def ciprofloxacin() -> Breakpoint:
    """EUCAST E. coli ciprofloxacin: S <= 0.25, R > 0.25.

    No intermediate band, so anything above ``s_max`` is resistant.
    """
    return Breakpoint(organism="e. coli", drug="ciprofloxacin", s_max=0.25, standard="EUCAST")


@pytest.fixture
def gentamicin() -> Breakpoint:
    """A breakpoint with an explicit intermediate band."""
    return Breakpoint(
        organism="e. coli",
        drug="gentamicin",
        s_max=2.0,
        i_max=4.0,
        r_min=8.0,
        standard="CLSI",
    )


# ---------------------------------------------------------------------------
# classify_mic
# ---------------------------------------------------------------------------


class TestClassifyMic:
    def test_no_band_above_s_max_is_resistant(self, ciprofloxacin) -> None:
        assert classify_mic(0.125, ciprofloxacin) == SUSCEPTIBLE
        assert classify_mic(0.25, ciprofloxacin) == SUSCEPTIBLE
        assert classify_mic(0.5, ciprofloxacin) == RESISTANT
        assert classify_mic(64.0, ciprofloxacin) == RESISTANT

    def test_intermediate_band(self, gentamicin) -> None:
        assert classify_mic(1.0, gentamicin) == SUSCEPTIBLE
        assert classify_mic(2.0, gentamicin) == SUSCEPTIBLE
        assert classify_mic(4.0, gentamicin) == INTERMEDIATE
        assert classify_mic(8.0, gentamicin) == RESISTANT
        assert classify_mic(16.0, gentamicin) == RESISTANT

    def test_no_breakpoint_is_unknown_not_susceptible(self) -> None:
        """An unknown phenotype must never read as a negative one.

        If this returned ``S``, a cohort with no breakpoints at all
        would score perfect agreement against any caller.
        """
        assert classify_mic(0.01, None) == UNKNOWN

    def test_gap_between_bands_without_i_is_unknown(self) -> None:
        """A gap the standard does not define stays unknown."""
        bp = Breakpoint(organism="e coli", drug="ampicillin", s_max=8.0, r_min=32.0)
        assert classify_mic(4.0, bp) == SUSCEPTIBLE
        assert classify_mic(8.0, bp) == SUSCEPTIBLE
        # Between 8 and 32 with no i_max: EUCAST would say R, but the
        # table has not said so, and guessing invents a phenotype.
        assert classify_mic(16.0, bp) == UNKNOWN
        assert classify_mic(32.0, bp) == RESISTANT
        assert classify_mic(64.0, bp) == RESISTANT

    def test_only_r_min_means_below_is_susceptible(self) -> None:
        bp = Breakpoint(organism="e coli", drug="chloramphenicol", r_min=32.0)
        assert classify_mic(0.5, bp) == SUSCEPTIBLE
        assert classify_mic(31.0, bp) == SUSCEPTIBLE
        assert classify_mic(32.0, bp) == RESISTANT

    def test_only_s_max_means_above_is_resistant(self, ciprofloxacin) -> None:
        assert ciprofloxacin.r_min is None
        assert classify_mic(0.25, ciprofloxacin) == SUSCEPTIBLE
        assert classify_mic(0.5, ciprofloxacin) == RESISTANT

    def test_non_positive_mic_is_unknown(self, ciprofloxacin) -> None:
        assert classify_mic(0.0, ciprofloxacin) == UNKNOWN
        assert classify_mic(-1.0, ciprofloxacin) == UNKNOWN

    def test_non_finite_mic_is_unknown(self, ciprofloxacin) -> None:
        assert classify_mic(float("nan"), ciprofloxacin) == UNKNOWN
        assert classify_mic(float("inf"), ciprofloxacin) == UNKNOWN


# ---------------------------------------------------------------------------
# MicRecord
# ---------------------------------------------------------------------------


class TestMicRecord:
    def test_log2_mic(self) -> None:
        assert _mic("A", "ciprofloxacin", 0.25).log2_mic == pytest.approx(-2.0)
        assert _mic("A", "ciprofloxacin", 1.0).log2_mic == pytest.approx(0.0)
        assert _mic("A", "ciprofloxacin", 64.0).log2_mic == pytest.approx(6.0)

    def test_rejects_non_positive_mic(self) -> None:
        with pytest.raises(ValueError, match="positive finite"):
            MicRecord(accession="A", drug="ciprofloxacin", mic=0.0)
        with pytest.raises(ValueError, match="positive finite"):
            MicRecord(accession="A", drug="ciprofloxacin", mic=-1.0)
        with pytest.raises(ValueError, match="positive finite"):
            MicRecord(accession="A", drug="ciprofloxacin", mic=float("nan"))

    def test_key(self) -> None:
        assert _mic("A", "ciprofloxacin", 1.0).key == ("A", "ciprofloxacin")


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------


class TestLoaders:
    def test_read_mic_records_tsv(self, tmp_path: Path) -> None:
        path = tmp_path / "mic.tsv"
        path.write_text(
            "Accession\tDrug\tMIC\tOrganism\n"
            "GCF_0001\tciprofloxacin\t0.25\tEscherichia coli\n"
            "GCF_0002\tciprofloxacin\t4\tEscherichia coli\n",
            encoding="utf-8",
        )
        records = read_mic_records(path)
        assert [(r.accession, r.drug, r.mic) for r in records] == [
            ("GCF_0001", "ciprofloxacin", 0.25),
            ("GCF_0002", "ciprofloxacin", 4.0),
        ]

    def test_read_mic_records_csv_with_spacey_headers(self, tmp_path: Path) -> None:
        path = tmp_path / "mic.csv"
        path.write_text(
            "Genome ID,Antibiotic,MIC Value\nGCF_0001,ciprofloxacin,0.5\n",
            encoding="utf-8",
        )
        records = read_mic_records(path)
        assert len(records) == 1
        assert records[0].accession == "GCF_0001"
        assert records[0].drug == "ciprofloxacin"

    def test_read_mic_records_censored_and_off_dilution(self, tmp_path: Path) -> None:
        path = tmp_path / "mic.tsv"
        path.write_text(
            "accession\tdrug\tmic\nA\tcip\t<=0.25\nB\tcip\t>64\nC\tcip\t0.5\n",
            encoding="utf-8",
        )
        records = read_mic_records(path)
        assert [(r.accession, r.mic) for r in records] == [
            ("A", 0.25),
            ("B", 64.0),
            ("C", 0.5),
        ]

    def test_read_mic_records_skips_unparseable(self, tmp_path: Path) -> None:
        """An unreadable MIC is skipped, never coerced to zero.

        Coercing would manufacture a susceptible genome and inflate
        agreement.
        """
        path = tmp_path / "mic.tsv"
        path.write_text(
            "accession\tdrug\tmic\nA\tcip\t0.5\nB\tcip\tnot-a-number\nC\tcip\t\nD\tcip\t0.25\n",
            encoding="utf-8",
        )
        records = read_mic_records(path)
        assert [r.accession for r in records] == ["A", "D"]

    def test_read_mic_records_missing_columns(self, tmp_path: Path) -> None:
        path = tmp_path / "bad.tsv"
        path.write_text("accession\tdrug\nA\tcip\n", encoding="utf-8")
        with pytest.raises(SystemExit, match="needs"):
            read_mic_records(path)

    def test_read_breakpoints_requires_a_threshold(self, tmp_path: Path) -> None:
        """A row with no threshold at all cannot classify anything."""
        path = tmp_path / "bp.tsv"
        path.write_text("organism\tdrug\nEscherichia coli\tciprofloxacin\n", encoding="utf-8")
        with pytest.raises(SystemExit, match="at least one of"):
            read_breakpoints(path)

    def test_read_breakpoints_keys_lowercased(self, tmp_path: Path) -> None:
        path = tmp_path / "bp.tsv"
        path.write_text(
            "organism\tdrug\ts_max\tr_min\tstandard\n"
            "Escherichia Coli\tCiprofloxacin\t0.25\t\tEUCAST\n",
            encoding="utf-8",
        )
        table = read_breakpoints(path)
        assert ("escherichia coli", "ciprofloxacin") in table
        bp = table[("escherichia coli", "ciprofloxacin")]
        assert bp.s_max == 0.25
        assert bp.r_min is None
        assert bp.standard == "EUCAST"

    def test_read_breakpoints_drops_thresholdless_rows(self, tmp_path: Path) -> None:
        path = tmp_path / "bp.tsv"
        path.write_text(
            "organism\tdrug\ts_max\nE coli\tcip\t0.25\nE coli\t\t\nE coli\tcolistin\t\n",
            encoding="utf-8",
        )
        table = read_breakpoints(path)
        assert set(table) == {("e coli", "cip")}

    def test_read_rules(self, tmp_path: Path) -> None:
        path = tmp_path / "rules.tsv"
        path.write_text(
            "gene\tdrug\tallele\n"
            "gyrA\tCiprofloxacin\tS83L\n"
            "gyrA\tciprofloxacin\tS83L\n"
            "parC\tCiprofloxacin\tS80I\n"
            "acrB\tcefazolin\t\n",
            encoding="utf-8",
        )
        rules = read_rules(path)
        # gyrA's two rows are the same gene-drug pair once lowercased,
        # so they collapse to one relationship rather than two.
        assert rules.predicts_resistance("gyrA", "ciprofloxacin")
        assert rules.predicts_resistance("gyrA", "Ciprofloxacin")
        assert rules.predicts_resistance("parC", "CIPROFLOXACIN")
        assert rules.predicts_resistance("acrB", "cefazolin")
        assert not rules.predicts_resistance("acrB", "ciprofloxacin")
        assert not rules.predicts_resistance("gyrA", "gentamicin")

    def test_read_rules_drugs_for(self, tmp_path: Path) -> None:
        path = tmp_path / "rules.tsv"
        path.write_text("gene\tdrug\ngyrA\tciprofloxacin\nacrB\tcefazolin\n", encoding="utf-8")
        rules = read_rules(path)
        assert rules.drugs_for(["gyrA", "acrB"]) == {"ciprofloxacin", "cefazolin"}
        assert rules.drugs_for([]) == set()
        assert rules.drugs_for(["unknown_gene"]) == set()

    def test_dialect_picks_tab_over_comma(self, tmp_path: Path) -> None:
        """A table with both separators is decided by the header count.

        ``csv.Sniffer`` raises "Could not determine delimiter" on small
        tables, which is exactly the case that matters: a curated
        breakpoint file with one data row and an empty optional column.
        """
        from benchmarks.phenotypic import _dialect_for

        path = tmp_path / "bp.tsv"
        path.write_text("organism\tdrug\ts_max\nE coli\tcip\t0.25\n", encoding="utf-8")
        assert _dialect_for(path) is csv.excel_tab

    def test_dialect_picks_comma_when_dominant(self, tmp_path: Path) -> None:
        from benchmarks.phenotypic import _dialect_for

        path = tmp_path / "bp.csv"
        path.write_text("organism,drug,s_max\nE coli,cip,0.25\n", encoding="utf-8")
        assert _dialect_for(path) is csv.excel

    def test_dialect_survives_empty_optional_column(self, tmp_path: Path) -> None:
        """The row the sniffer chokes on: a blank intermediate cell."""
        from benchmarks.phenotypic import _dialect_for

        path = tmp_path / "bp.tsv"
        path.write_text(
            "organism\tdrug\ts_max\ti_max\tr_min\tstandard\n"
            "Escherichia coli\tgentamicin\t2\t4\t8\tCLSI\n",
            encoding="utf-8",
        )
        assert _dialect_for(path) is csv.excel_tab

    def test_comments_and_blank_lines_do_not_become_the_header(self, tmp_path: Path) -> None:
        """A documented table's prose is not its header row.

        A table that carries its own provenance and caveats is more
        useful than one that does not; the loader has to look past
        them. A blank line between the prose and the column names is
        especially easy to trip over, because it is neither a comment
        nor a field.
        """
        path = tmp_path / "bp.tsv"
        path.write_text(
            "# Breakpoints for the drugs this project bears on.\n"
            "#\n"
            "# Two conventions, and the difference matters.\n"
            "\n"
            "organism\tdrug\ts_max\tr_min\tstandard\n"
            "Escherichia coli\tciprofloxacin\t0.25\t\tEUCAST\n",
            encoding="utf-8",
        )
        table = read_breakpoints(path)
        assert set(table) == {("escherichia coli", "ciprofloxacin")}
        assert table[("escherichia coli", "ciprofloxacin")].s_max == 0.25

    def test_comment_lines_are_skipped_but_rows_are_not(self, tmp_path: Path) -> None:
        """A '#' *inside* a data field never appears at line start."""
        path = tmp_path / "bp.tsv"
        path.write_text(
            "organism\tdrug\ts_max\tsource\n"
            "Escherichia coli\tciprofloxacin\t0.25\ttable #4\n"
            "# a mid-file comment\n"
            "Escherichia coli\tgentamicin\t2\tCLSI\n",
            encoding="utf-8",
        )
        table = read_breakpoints(path)
        assert set(table) == {
            ("escherichia coli", "ciprofloxacin"),
            ("escherichia coli", "gentamicin"),
        }


# ---------------------------------------------------------------------------
# Breakpoint lookup
# ---------------------------------------------------------------------------


class TestBreakpointLookup:
    def test_organism_beats_empty_fallback(self, ciprofloxacin) -> None:
        """A species-specific row must win over the generic one."""
        generic = Breakpoint("", "ciprofloxacin", s_max=0.06)
        table = _table(generic, ciprofloxacin)
        record = _mic("A", "ciprofloxacin", 0.125, organism="e. coli")
        from benchmarks.phenotypic import _lookup_breakpoint

        assert _lookup_breakpoint(table, record) is ciprofloxacin

    def test_empty_fallback_when_no_species_row(self) -> None:
        table = {("", "gentamicin"): Breakpoint("", "gentamicin", s_max=2.0, r_min=8.0)}
        record = _mic("A", "gentamicin", 0.5, organism="k. pneumoniae")
        from benchmarks.phenotypic import _lookup_breakpoint

        assert _lookup_breakpoint(table, record) is table[("", "gentamicin")]

    def test_missing_lookup_gives_unknown(self) -> None:
        from benchmarks.phenotypic import _lookup_breakpoint

        assert _lookup_breakpoint({}, _mic("A", "gentamicin", 0.5)) is None


# ---------------------------------------------------------------------------
# Gene-drug association
# ---------------------------------------------------------------------------


class TestGeneDrugAssociation:
    def test_perfect_association_is_significant(self, ciprofloxacin) -> None:
        """A gene carried by every resistant genome and no other."""
        records = [_mic(f"G{i}", "ciprofloxacin", 0.5 if i < 6 else 0.125) for i in range(12)]
        carriers = {
            **{f"G{i}": {"gyrA": True} for i in range(6)},
            **{f"G{i}": {"gyrA": False} for i in range(6, 12)},
        }
        associations = gene_drug_associations(carriers, records, _table(ciprofloxacin))

        assert len(associations) == 1
        item = associations[0]
        assert item.gene == "gyrA"
        assert item.n_carrier == 6
        assert item.n_resistant == 6
        assert item.n_both == 6
        assert item.n_evaluable == 12
        assert item.p_value < 0.01
        assert item.significant

    def test_no_carriers_means_no_test(self, ciprofloxacin) -> None:
        """A gene nobody carries must not produce a test with no data."""
        records = [_mic("A", "ciprofloxacin", 0.5), _mic("B", "ciprofloxacin", 0.125)]
        carriers: dict[str, dict[str, bool]] = {"A": {}, "B": {}}
        assert gene_drug_associations(carriers, records, _table(ciprofloxacin)) == []

    def test_no_resistant_genomes_means_no_test(self, ciprofloxacin) -> None:
        records = [_mic("A", "ciprofloxacin", 0.125), _mic("B", "ciprofloxacin", 0.06)]
        carriers = {"A": {"gyrA": True}, "B": {"gyrA": True}}
        assert gene_drug_associations(carriers, records, _table(ciprofloxacin)) == []

    def test_missing_breakpoint_skips_drug(self, ciprofloxacin) -> None:
        """A drug with no breakpoint has no testable phenotype."""
        records = [_mic("A", "ciprofloxacin", 0.5), _mic("A", "gentamicin", 32.0)]
        carriers = {"A": {"gyrA": True, "aac": True}}
        associations = gene_drug_associations(carriers, records, _table(ciprofloxacin))
        assert {item.drug for item in associations} == {"ciprofloxacin"}

    def test_unknown_phenotype_is_excluded_from_the_table(self, gentamicin) -> None:
        """MICs in an undefined band are not counted as susceptible."""
        records = [
            _mic("A", "gentamicin", 0.5),
            _mic("B", "gentamicin", 16.0),
            # 16 with s_max=2, r_min=8 is resistant; use the gap instead.
            _mic("C", "gentamicin", 0.5),
        ]
        carriers = {"A": {"aac": True}, "B": {"aac": True}, "C": {"aac": False}}
        associations = gene_drug_associations(carriers, records, _table(gentamicin))
        assert all(item.n_evaluable == 3 for item in associations)

    def test_genome_without_mic_is_not_a_negative(self, ciprofloxacin) -> None:
        """A genome whose MIC is unknown contributes nothing."""
        records = [_mic("A", "ciprofloxacin", 0.5), _mic("B", "ciprofloxacin", 0.125)]
        # "C" carries the gene but has no MIC for this drug.
        carriers = {"A": {"gyrA": True}, "B": {"gyrA": False}, "C": {"gyrA": True}}
        associations = gene_drug_associations(carriers, records, _table(ciprofloxacin))
        item = associations[0]
        assert item.n_carrier == 1
        assert item.n_evaluable == 2

    def test_bh_alpha_is_honoured(self) -> None:
        """`significant` reads the threshold the caller asked for."""
        items = [
            GeneDrugAssociation(gene="g", drug="d", p_value=0.01, q_value=0.01),
            GeneDrugAssociation(gene="h", drug="d", p_value=0.6, q_value=0.6),
        ]
        assert items[0].significant  # default alpha 0.05
        assert not items[1].significant
        items[0].alpha = 0.001
        assert not items[0].significant

    def test_empty_inputs_return_no_tests(self, ciprofloxacin) -> None:
        assert gene_drug_associations({}, [], _table(ciprofloxacin)) == []

    def test_as_row_shape(self, ciprofloxacin) -> None:
        records = [_mic("A", "ciprofloxacin", 0.5), _mic("B", "ciprofloxacin", 0.125)]
        carriers = {"A": {"gyrA": True}, "B": {"gyrA": False}}
        (item,) = gene_drug_associations(carriers, records, _table(ciprofloxacin))
        row = item.as_row()
        assert set(row) == {
            "gene",
            "drug",
            "n_carrier",
            "n_resistant",
            "n_both",
            "n_evaluable",
            "odds_ratio",
            "p_value",
            "q_value",
            "significant",
        }


# ---------------------------------------------------------------------------
# Burden correlation
# ---------------------------------------------------------------------------


class TestBurdenCorrelation:
    def test_burden_tracks_mic(self) -> None:
        """More mutated genes, higher MIC: rho should be positive."""
        records = [
            _mic("A", "ciprofloxacin", 0.06),
            _mic("B", "ciprofloxacin", 0.125),
            _mic("C", "ciprofloxacin", 0.25),
            _mic("D", "ciprofloxacin", 0.5),
            _mic("E", "ciprofloxacin", 1.0),
        ]
        carriers = {
            "A": {},
            "B": {"g1": True},
            "C": {"g1": True, "g2": True},
            "D": {"g1": True, "g2": True, "g3": True},
            "E": {"g1": True, "g2": True, "g3": True},
        }
        table = {("", "ciprofloxacin"): Breakpoint("", "ciprofloxacin", s_max=0.25)}
        results = burden_correlation(carriers, records, table)
        (item,) = results
        assert item.n == 5
        assert item.rho > 0.9
        assert item.n_carriers == 4

    def test_constant_burden_is_not_a_test(self) -> None:
        """No variance in burden means no correlation to compute."""
        records = [_mic("A", "ciprofloxacin", 1.0), _mic("B", "ciprofloxacin", 2.0)]
        carriers = {"A": {"g1": True}, "B": {"g1": True}}
        table = {("", "ciprofloxacin"): Breakpoint("", "ciprofloxacin", s_max=0.25)}
        (item,) = burden_correlation(carriers, records, table)
        assert item.rho == 0.0
        assert item.p_value == 1.0
        assert not item.significant

    def test_constant_mic_is_not_a_test(self) -> None:
        records = [_mic("A", "ciprofloxacin", 1.0), _mic("B", "ciprofloxacin", 1.0)]
        carriers = {"A": {}, "B": {"g1": True}}
        table = {("", "ciprofloxacin"): Breakpoint("", "ciprofloxacin", s_max=0.25)}
        (item,) = burden_correlation(carriers, records, table)
        assert not item.significant

    def test_drug_without_breakpoint_is_skipped(self) -> None:
        """Breakpoints are needed only to know the drug is evaluated."""
        records = [_mic("A", "ciprofloxacin", 1.0)]
        assert burden_correlation({"A": {}}, records, {}) == []

    def test_single_genome_is_not_a_test(self) -> None:
        records = [_mic("A", "ciprofloxacin", 1.0)]
        table = {("", "ciprofloxacin"): Breakpoint("", "ciprofloxacin", s_max=0.25)}
        (item,) = burden_correlation({"A": {"g1": True}}, records, table)
        assert item.n == 1
        assert item.rho == 0.0

    def test_no_carriers_gives_zero_burden(self) -> None:
        records = [_mic("A", "ciprofloxacin", 0.5), _mic("B", "ciprofloxacin", 1.0)]
        table = {("", "ciprofloxacin"): Breakpoint("", "ciprofloxacin", s_max=0.25)}
        (item,) = burden_correlation({"A": {}, "B": {}}, records, table)
        assert item.n_carriers == 0
        assert not item.significant

    def test_gene_filter_narrows_the_count(self) -> None:
        """Only phenotypically relevant genes add to the burden."""
        records = [_mic("A", "ciprofloxacin", 0.5), _mic("B", "ciprofloxacin", 1.0)]
        carriers = {"A": {"g1": True, "noise": True}, "B": {"noise": True}}
        table = {("", "ciprofloxacin"): Breakpoint("", "ciprofloxacin", s_max=0.25)}
        (item,) = burden_correlation(
            carriers, records, table, resistance_genes={"g1": frozenset({"ciprofloxacin"})}
        )
        assert item.n_carriers == 1

    def test_empty_filter_is_not_a_filter(self) -> None:
        """An empty gene filter must not filter every gene out.

        The CLI derives its filter from the significant associations;
        when nothing is significant the filter is empty, and passing
        ``{}`` down would mean "no gene is phenotypically relevant" --
        reporting zero carriers for a cohort that plainly has some.
        """
        records = [_mic("A", "ciprofloxacin", 0.5), _mic("B", "ciprofloxacin", 1.0)]
        carriers = {"A": {"gyrA": True}, "B": {}}
        table = {("", "ciprofloxacin"): Breakpoint("", "ciprofloxacin", s_max=0.25)}
        (item,) = burden_correlation(carriers, records, table, resistance_genes={})
        assert item.n_carriers == 1


# ---------------------------------------------------------------------------
# Concordance
# ---------------------------------------------------------------------------


class TestConcordance:
    def test_perfect_agreement(self, ciprofloxacin) -> None:
        records = [_mic("A", "ciprofloxacin", 0.5), _mic("B", "ciprofloxacin", 0.125)]
        carriers = {"A": {"gyrA": True}, "B": {}}
        rules = PhenotypeRules({"gyrA": frozenset({"ciprofloxacin"})})
        result = concordance(carriers, records, _table(ciprofloxacin), rules)
        assert result.n_evaluable == 2
        assert result.agreement == 1.0
        assert result.very_major_error == 0.0
        assert result.major_error == 0.0

    def test_very_major_error_is_the_dangerous_direction(self, ciprofloxacin) -> None:
        """Reported susceptible, measured resistant."""
        records = [_mic("A", "ciprofloxacin", 0.5)]  # resistant
        carriers = {"A": {}}  # no resistance gene called
        rules = PhenotypeRules({"gyrA": frozenset({"ciprofloxacin"})})
        result = concordance(carriers, records, _table(ciprofloxacin), rules)
        assert result.n_very_major == 1
        assert result.n_major == 0
        assert result.very_major_error == 1.0
        assert result.major_error == 0.0

    def test_major_error_is_the_other_direction(self, ciprofloxacin) -> None:
        """Reported resistant, measured susceptible."""
        records = [_mic("A", "ciprofloxacin", 0.125)]  # susceptible
        carriers = {"A": {"gyrA": True}}
        rules = PhenotypeRules({"gyrA": frozenset({"ciprofloxacin"})})
        result = concordance(carriers, records, _table(ciprofloxacin), rules)
        assert result.n_very_major == 0
        assert result.n_major == 1
        assert result.major_error == 1.0

    def test_unmapped_drug_counts_as_not_predicted(self, ciprofloxacin) -> None:
        records = [_mic("A", "gentamicin", 32.0)]
        carriers = {"A": {"gyrA": True}}
        rules = PhenotypeRules({"gyrA": frozenset({"ciprofloxacin"})})
        result = concordance(carriers, records, _table(ciprofloxacin), rules)
        assert result.n_evaluable == 0

    def test_unknown_phenotype_is_excluded_not_counted(self, ciprofloxacin) -> None:
        """A MIC with no breakpoint does not silently score as correct."""
        records = [_mic("A", "gentamicin", 32.0), _mic("B", "ciprofloxacin", 0.125)]
        carriers = {"A": {"gyrA": True}, "B": {}}
        rules = PhenotypeRules({"gyrA": frozenset({"ciprofloxacin"})})
        result = concordance(carriers, records, _table(ciprofloxacin), rules)
        assert result.n_unknown_phenotype == 1
        assert result.n_evaluable == 1
        assert result.agreement == 1.0

    def test_intermediate_phenotype_is_minor_error(self, gentamicin) -> None:
        records = [_mic("A", "gentamicin", 4.0)]  # I
        carriers = {"A": {}}
        rules = PhenotypeRules({"aac": frozenset({"gentamicin"})})
        result = concordance(carriers, records, _table(gentamicin), rules)
        assert result.n_evaluable == 1
        assert result.n_minor == 1
        assert result.n_very_major == 0
        assert result.n_major == 0
        assert result.agreement == 0.0
        assert result.minor_error == 1.0

    def test_genome_without_carriers_is_susceptible(self, ciprofloxacin) -> None:
        records = [_mic("A", "ciprofloxacin", 0.125)]
        result = concordance({"A": {}}, records, _table(ciprofloxacin), PhenotypeRules())
        assert result.agreement == 1.0

    def test_empty_everything(self, ciprofloxacin) -> None:
        result = concordance({}, [], _table(ciprofloxacin), PhenotypeRules())
        assert result.n_evaluable == 0
        assert result.agreement == 0.0
        assert result.n_very_major == 0
        assert result.per_drug == {}

    def test_per_drug_breakdown(self, ciprofloxacin, gentamicin) -> None:
        records = [
            _mic("A", "ciprofloxacin", 0.5),
            _mic("A", "gentamicin", 32.0),
            _mic("B", "ciprofloxacin", 0.125),
            _mic("B", "gentamicin", 0.5),
        ]
        carriers = {"A": {"gyrA": True, "aac": True}, "B": {}}
        rules = PhenotypeRules(
            {
                "gyrA": frozenset({"ciprofloxacin"}),
                "aac": frozenset({"gentamicin"}),
            }
        )
        result = concordance(carriers, records, _table(ciprofloxacin, gentamicin), rules)
        assert set(result.per_drug) == {"ciprofloxacin", "gentamicin"}
        assert result.per_drug["ciprofloxacin"]["n"] == 2
        assert result.per_drug["ciprofloxacin"]["agreement"] == pytest.approx(1.0)
        assert result.per_drug["gentamicin"]["n"] == 2
        assert result.per_drug["gentamicin"]["agreement"] == pytest.approx(1.0)

    def test_as_row_shape(self, ciprofloxacin) -> None:
        records = [_mic("A", "ciprofloxacin", 0.125), _mic("B", "ciprofloxacin", 0.5)]
        carriers = {"A": {}, "B": {"gyrA": True}}
        rules = PhenotypeRules({"gyrA": frozenset({"ciprofloxacin"})})
        row = concordance(carriers, records, _table(ciprofloxacin), rules).as_row()
        assert set(row) == {
            "n_evaluable",
            "n_unknown_phenotype",
            "agreement",
            "very_major_error",
            "major_error",
            "minor_error",
            "n_very_major",
            "n_major",
            "n_minor",
            "per_drug",
        }


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


class TestWriteReport:
    def test_writes_both_formats(self, tmp_path: Path) -> None:
        paths = write_report([], [], None, tmp_path)
        assert Path(paths["markdown"]).is_file()
        assert Path(paths["json"]).is_file()

    def test_json_payload_without_concordance(self, tmp_path: Path) -> None:
        """The supervised section is optional and reports as null."""
        paths = write_report([], [], None, tmp_path)
        payload = json.loads(Path(paths["json"]).read_text(encoding="utf-8"))
        assert payload["gene_drug_associations"] == []
        assert payload["burden_correlations"] == []
        assert payload["concordance"] is None

    def test_report_mentions_data_sources(self, tmp_path: Path) -> None:
        paths = write_report([], [], None, tmp_path)
        text = Path(paths["markdown"]).read_text(encoding="utf-8")
        assert "BV-BRC" in text
        assert "EUCAST" in text
        assert "CLSI" in text

    def test_empty_sections_are_labelled_not_blank(self, tmp_path: Path) -> None:
        paths = write_report([], [], None, tmp_path)
        text = Path(paths["markdown"]).read_text(encoding="utf-8")
        assert "_none_" in text
        assert "0 of 0 tests significant" in text
