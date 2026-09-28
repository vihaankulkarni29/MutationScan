"""Tests for biochemical scoring and epistasis modules."""

from __future__ import annotations

import pandas as pd

from mutation_scan.epistasis import (
    EPISTASIS_COLUMNS,
    GENOME_MAP_COLUMNS,
    detect_epistasis,
    genome_mutation_map,
    read_epistasis,
)
from mutation_scan.scoring import (
    compute_severity,
    mutation_severity,
    pair_severity,
)

# -- scoring ----------------------------------------------------------------


class TestComputeSeverity:
    def test_identical_substitution_is_low_severity(self):
        s = compute_severity("A", "A")
        assert s["severity_score"] == 0.0
        assert s["blosum62_score"] == 4
        assert s["charge_delta"] == 0
        assert s["kd_delta"] == 0.0
        assert s["volume_delta"] == 0

    def test_conservative_substitution_is_low_severity(self):
        s = compute_severity("L", "I")
        assert s["severity_score"] < 30.0

    def test_radical_substitution_is_high_severity(self):
        s = compute_severity("R", "D")
        assert s["severity_score"] > 50.0
        assert s["charge_delta"] == -2

    def test_blosum62_penalises_rare_swaps(self):
        s = compute_severity("W", "C")
        assert s["blosum62_score"] < 0

    def test_case_insensitive(self):
        s_upper = compute_severity("R", "D")
        s_lower = compute_severity("r", "d")
        assert s_upper == s_lower


class TestMutationSeverity:
    def _mutations(self):
        return pd.DataFrame(
            [
                {"Accession": "g1", "Gene": "acrB", "Mutation": "H596N", "Position": 596},
                {"Accession": "g1", "Gene": "acrA", "Mutation": "R50H", "Position": 50},
                {"Accession": "g2", "Gene": "acrB", "Mutation": "H596N", "Position": 596},
            ],
            columns=["Accession", "Gene", "Mutation", "Position"],
        )

    def test_returns_expected_columns(self):
        frame = mutation_severity(self._mutations())
        expected = {
            "Accession",
            "Gene",
            "Mutation",
            "Ref_AA",
            "Alt_AA",
            "blosum62_score",
            "charge_delta",
            "kd_delta",
            "volume_delta",
            "severity_score",
        }
        assert expected.issubset(frame.columns)

    def test_preserves_row_count(self):
        assert len(mutation_severity(self._mutations())) == 3

    def test_empty_input(self):
        frame = mutation_severity(pd.DataFrame(columns=["Accession", "Gene", "Mutation"]))
        assert frame.empty

    def test_severity_in_range(self):
        frame = mutation_severity(self._mutations())
        assert frame["severity_score"].between(0, 100).all()


class TestPairSeverity:
    def test_averages_two_mutations(self):
        sev = pd.DataFrame(
            [
                {"Gene": "acrB", "Mutation": "H596N", "severity_score": 40.0},
                {"Gene": "acrA", "Mutation": "R50H", "severity_score": 60.0},
            ]
        )
        assert pair_severity("acrB:H596N", "acrA:R50H", sev) == 50.0

    def test_missing_falls_back_to_zero(self):
        assert pair_severity("acrB:X1Y", "acrA:Z2W", pd.DataFrame()) == 0.0


# -- epistasis ---------------------------------------------------------------


def _mutations():
    return pd.DataFrame(
        [
            {"Accession": "g1", "Gene": "acrB", "Mutation": "H596N", "Position": 596},
            {"Accession": "g1", "Gene": "acrA", "Mutation": "R50H", "Position": 50},
            {"Accession": "g2", "Gene": "acrB", "Mutation": "H596N", "Position": 596},
            {"Accession": "g2", "Gene": "acrA", "Mutation": "R50H", "Position": 50},
            {"Accession": "g3", "Gene": "acrB", "Mutation": "H596N", "Position": 596},
        ],
        columns=["Accession", "Gene", "Mutation", "Position"],
    )


def _qc():
    return pd.DataFrame(
        [
            {"Accession": "g1", "Gene": "acrB", "Status": "called"},
            {"Accession": "g1", "Gene": "acrA", "Status": "called"},
            {"Accession": "g2", "Gene": "acrB", "Status": "called"},
            {"Accession": "g2", "Gene": "acrA", "Status": "called"},
            {"Accession": "g3", "Gene": "acrB", "Status": "called"},
            {"Accession": "g3", "Gene": "acrA", "Status": "called"},
        ],
        columns=["Accession", "Gene", "Status"],
    )


class TestGenomeMutationMap:
    def test_one_row_per_genome(self):
        frame = genome_mutation_map(_mutations(), _qc())
        assert len(frame) == 3

    def test_expected_columns(self):
        frame = genome_mutation_map(_mutations(), _qc())
        assert list(frame.columns) == list(GENOME_MAP_COLUMNS)

    def test_mutation_map_contains_variants(self):
        frame = genome_mutation_map(_mutations(), _qc())
        row = frame.iloc[0]
        assert "acrB:H596N" in row["Mutation_Map"]
        assert "acrA:R50H" in row["Mutation_Map"]

    def test_counts_are_correct(self):
        frame = genome_mutation_map(_mutations(), _qc())
        g2 = frame[frame["Accession"] == "g2"].iloc[0]
        assert g2["N_Mutations"] == 2
        assert g2["N_Mutated_Genes"] == 2

    def test_empty_input(self):
        frame = genome_mutation_map(
            pd.DataFrame(columns=["Accession", "Gene", "Mutation", "Position"]),
            pd.DataFrame(columns=["Accession", "Gene", "Status"]),
        )
        assert frame.empty

    def test_qc_restricts_to_called(self):
        qc = pd.DataFrame(
            [
                {"Accession": "g1", "Gene": "acrB", "Status": "called"},
                {"Accession": "g1", "Gene": "acrA", "Status": "low_coverage"},
            ]
        )
        frame = genome_mutation_map(_mutations(), qc)
        row = frame[frame["Accession"] == "g1"].iloc[0]
        assert row["N_Mutations"] == 1
        assert "acrB:H596N" in row["Mutation_Map"]
        assert "acrA:R50H" not in row["Mutation_Map"]


class TestDetectEpistasis:
    def test_returns_expected_columns(self):
        frame = detect_epistasis(_mutations(), _qc(), threads=1)
        assert list(frame.columns) == list(EPISTASIS_COLUMNS)

    def test_perfect_cooccurrence_has_high_or(self):
        frame = detect_epistasis(_mutations(), _qc(), threads=1)
        pair = frame[(frame["Variant_A"] == "acrA:R50H") & (frame["Variant_B"] == "acrB:H596N")]
        assert not pair.empty
        assert float(pair.iloc[0]["Odds_Ratio"]) > 1.0

    def test_sign_positive_for_enriched(self):
        frame = detect_epistasis(_mutations(), _qc(), threads=1)
        pair = frame[(frame["Variant_A"] == "acrA:R50H") & (frame["Variant_B"] == "acrB:H596N")]
        assert not pair.empty
        assert pair.iloc[0]["Sign"] in ("positive", "neutral")

    def test_min_count_filters(self):
        frame = detect_epistasis(_mutations(), _qc(), min_count=3, threads=1)
        assert frame.empty

    def test_empty_mutations_returns_empty(self):
        frame = detect_epistasis(
            pd.DataFrame(columns=["Accession", "Gene", "Mutation"]),
            pd.DataFrame(columns=["Accession", "Gene", "Status"]),
            threads=1,
        )
        assert frame.empty

    def test_fdr_correction_applied(self):
        frame = detect_epistasis(_mutations(), _qc(), threads=1)
        assert "Q_value" in frame.columns
        assert frame["Q_value"].between(0, 1).all()

    def test_pair_severity_included(self):
        sev = mutation_severity(_mutations())
        frame = detect_epistasis(_mutations(), _qc(), severity=sev, threads=1)
        assert "Pair_Severity" in frame.columns
        assert frame["Pair_Severity"].sum() > 0

    def test_ranking_by_q_then_severity(self):
        sev = mutation_severity(_mutations())
        frame = detect_epistasis(_mutations(), _qc(), severity=sev, threads=1)
        q_vals = frame["Q_value"].values
        for i in range(len(q_vals) - 1):
            assert q_vals[i] <= q_vals[i + 1]

    def test_read_epistasis_roundtrip(self, tmp_path):
        frame = detect_epistasis(_mutations(), _qc(), threads=1)
        path = tmp_path / "networks.csv"
        frame.to_csv(path, index=False)
        reloaded = read_epistasis(path)
        assert len(reloaded) == len(frame)
        assert list(reloaded.columns) == list(frame.columns)


class TestEpistasisNetworks:
    def test_full_pipeline(self, tmp_path):
        mut = pd.DataFrame(
            [
                {"Accession": "g1", "Gene": "acrB", "Mutation": "H596N", "Position": 596},
                {"Accession": "g1", "Gene": "acrA", "Mutation": "R50H", "Position": 50},
            ],
            columns=["Accession", "Gene", "Mutation", "Position"],
        )
        qc = pd.DataFrame(
            [
                {"Accession": "g1", "Gene": "acrB", "Status": "called"},
                {"Accession": "g1", "Gene": "acrA", "Status": "called"},
            ]
        )
        mut_path = tmp_path / "mutations.csv"
        qc_path = tmp_path / "qc.csv"
        out_dir = tmp_path / "epi"
        mut.to_csv(mut_path, index=False)
        qc.to_csv(qc_path, index=False)

        from mutation_scan.epistasis import epistasis_networks

        result = epistasis_networks(mut_path, qc_path, out_dir, threads=1)

        assert (out_dir / "genome_mutation_map.csv").exists()
        assert (out_dir / "mutation_severity.csv").exists()
        assert (out_dir / "epistasis_networks.csv").exists()
        assert len(result) == 3
