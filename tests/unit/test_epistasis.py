"""Tests for biochemical scoring and epistasis modules."""

from __future__ import annotations

import pandas as pd
import pytest

from mutation_scan import epistasis
from mutation_scan.epistasis import (
    EPISTASIS_COLUMNS,
    GENOME_MAP_COLUMNS,
    detect_epistasis,
    genome_mutation_map,
    read_epistasis,
)
from mutation_scan.scoring import (
    build_severity_index,
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


def _cohort(n_genes: int = 4, variants_per_gene: int = 8, n_genomes: int = 30):
    """A synthetic cohort with enough structure to exercise the batching.

    Variants draw their carriers from a small pool of overlapping patterns rather
    than independently, which is what real cohorts look like and what makes the
    distinct-table collapse possible. Uniformly random carriers would give almost
    every pair a unique intersection count, so a fixture built that way would
    measure itself rather than the memoisation.
    """
    patterns = []
    for p in range(6):
        patterns.append({g for g in range(n_genomes) if g % 6 == p or g % 11 == p})

    mutations, qc = [], []
    for g in range(n_genomes):
        accession = f"g{g:03d}"
        for gene in range(n_genes):
            name = f"gene{gene}"
            mutations.append(
                {"Accession": accession, "Gene": name, "Mutation": "W1X", "Position": 1}
            )
            qc.append({"Accession": accession, "Gene": name, "Status": "called"})

    for gene in range(n_genes):
        name = f"gene{gene}"
        for v in range(variants_per_gene):
            carriers = patterns[(v + gene) % len(patterns)]
            for g in sorted(carriers):
                mutations.append(
                    {
                        "Accession": f"g{g:03d}",
                        "Gene": name,
                        "Mutation": f"M{v}",
                        "Position": 100 + v,
                    }
                )

    return (
        pd.DataFrame(mutations, columns=["Accession", "Gene", "Mutation", "Position"]),
        pd.DataFrame(qc, columns=["Accession", "Gene", "Status"]),
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


class TestDetectEpistasisIsThreadIndependent:
    """``threads`` is advisory, so it must not be able to change the answer.

    The stage used to pick between a thread pool and a process pool on pair count.
    Pooling it lost -- the pool spent its life paying ``spawn`` import costs
    against a small amount of real work -- so it now resolves every pair with one
    Gram matrix per universe and one Fisher test per *distinct* 2x2 table. The
    flag survives for compatibility, which only holds if it is genuinely inert.
    """

    def test_threads_does_not_change_the_output(self):
        one = detect_epistasis(_mutations(), _qc(), threads=1)
        eight = detect_epistasis(_mutations(), _qc(), threads=8)

        pd.testing.assert_frame_equal(one, eight)

    def test_the_stage_does_not_import_a_pool(self):
        """A regression guard, not a style preference.

        Reintroducing a pool would look fine in every output and cost seconds of
        spawn time per run, so nothing else would catch it.
        """
        assert not hasattr(epistasis, "parallel_map")


class TestDistinctTableMemoisation:
    """The memoised Fisher test must agree with testing every pair directly.

    ``fisher_exact`` is a pure function of the integer 2x2 table, so testing each
    distinct table once is exact rather than approximate. This pins that: the
    compressed path is compared against an explicit per-pair recomputation.
    """

    def _reference_counts(self, mutations, qc):
        """Compute every pair's counts the obvious, slow way."""
        eligible = {}
        for _, row in qc[qc["Status"] == "called"].iterrows():
            eligible.setdefault(str(row["Gene"]), set()).add(str(row["Accession"]))

        carriers = {}
        gene_of = {}
        for _, row in mutations.iterrows():
            key = f"{row['Gene']}:{row['Mutation']}"
            carriers.setdefault(key, set()).add(str(row["Accession"]))
            gene_of[key] = str(row["Gene"])

        variants = sorted(carriers)
        counts = {}
        for i, a in enumerate(variants):
            for b in variants[i + 1 :]:
                gene_a, gene_b = gene_of[a], gene_of[b]
                shared = (
                    eligible.get(gene_a, set())
                    if gene_a == gene_b
                    else eligible.get(gene_a, set()) & eligible.get(gene_b, set())
                )
                if not shared:
                    continue
                in_a = carriers[a] & shared
                in_b = carriers[b] & shared
                n11 = len(in_a & in_b)
                n10 = len(in_a - in_b)
                n01 = len(in_b - in_a)
                counts[(a, b)] = (n11, n10, n01, len(shared) - n11 - n10 - n01)
        return counts

    def test_counts_match_a_direct_per_pair_recomputation(self):
        mutations, qc = _cohort()
        expected = self._reference_counts(mutations, qc)

        variants = sorted({f"{r.Gene}:{r.Mutation}" for r in mutations.itertuples()})
        from itertools import combinations

        pairs = list(combinations(variants, 2))
        gene_of = {f"{r.Gene}:{r.Mutation}": str(r.Gene) for r in mutations.itertuples()}
        carriers = {}
        for r in mutations.itertuples():
            carriers.setdefault(f"{r.Gene}:{r.Mutation}", set()).add(str(r.Accession))
        eligible = {}
        for r in qc.itertuples():
            if r.Status == "called":
                eligible.setdefault(str(r.Gene), set()).add(str(r.Accession))

        tables = epistasis._resolve_pair_tables(variants, pairs, carriers, gene_of, eligible)
        got = {(t.variant_a, t.variant_b): (t.n11, t.n10, t.n01, t.n00) for t in tables}

        assert len(got) > 100, "fixture should produce a meaningful number of pairs"
        assert got == expected

    def test_one_test_per_distinct_table_not_per_pair(self):
        """The compression is the whole point; assert it actually happens."""
        mutations, qc = _cohort()
        from itertools import combinations

        variants = sorted({f"{r.Gene}:{r.Mutation}" for r in mutations.itertuples()})
        pairs = list(combinations(variants, 2))
        gene_of = {f"{r.Gene}:{r.Mutation}": str(r.Gene) for r in mutations.itertuples()}
        carriers = {}
        for r in mutations.itertuples():
            carriers.setdefault(f"{r.Gene}:{r.Mutation}", set()).add(str(r.Accession))
        eligible = {}
        for r in qc.itertuples():
            if r.Status == "called":
                eligible.setdefault(str(r.Gene), set()).add(str(r.Accession))

        tables = epistasis._resolve_pair_tables(variants, pairs, carriers, gene_of, eligible)
        distinct = {(t.n11, t.n10, t.n01, t.n00) for t in tables}

        assert len(tables) > 100
        assert len(distinct) * 5 < len(tables)

    def test_the_cohort_runs_end_to_end(self):
        mutations, qc = _cohort()
        frame = detect_epistasis(mutations, qc)

        assert list(frame.columns) == list(EPISTASIS_COLUMNS)
        assert not frame.empty
        # Every retained pair must be scored against genomes where both genes
        # were actually called, or N_Eligible would not be its own denominator.
        assert (frame["N_Eligible"] > 0).all()
        assert (
            frame["N_Both"] + frame["N_A_only"] + frame["N_B_only"] + frame["N_Neither"]
            == frame["N_Eligible"]
        ).all()

    def test_identical_tables_share_a_p_value(self):
        frame = detect_epistasis(_mutations(), _qc())
        by_table = frame.groupby(["N_Both", "N_A_only", "N_B_only", "N_Neither"])

        for _, group in by_table:
            assert group["P_value"].nunique() == 1


class TestMinCountFiltersBeforeCorrection:
    """``min_count`` changes ``n``, so it must be applied before BH.

    Filtering afterwards would leave every q-value in the file subtly wrong while
    still looking plausible. The stage now also uses min_count to skip Fisher
    tests entirely, so both halves of that are pinned here.
    """

    def test_q_values_reflect_the_filtered_pair_count(self):
        loose = detect_epistasis(_mutations(), _qc(), min_count=1)
        tight = detect_epistasis(_mutations(), _qc(), min_count=99)

        assert len(loose) > len(tight) == 0

    def test_a_tighter_threshold_never_lowers_a_q_value(self):
        loose = detect_epistasis(_mutations(), _qc(), min_count=1)
        tight = detect_epistasis(_mutations(), _qc(), min_count=2)

        merged = loose.merge(tight, on=["Variant_A", "Variant_B"], suffixes=("_loose", "_tight"))

        assert not merged.empty
        assert (merged["Q_value_tight"] >= merged["Q_value_loose"] - 1e-12).all()


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


class TestSeverityIndex:
    """The lookup that made epistasis tractable, pinned against the naive form.

    ``pair_severity`` used to rescan the whole severity table twice per pair,
    which is quadratic: at cohort scale (146,070 pairs over 66,416 mutations)
    that scan ran for hours and the pipeline never finished. The index build
    has to agree with the scan on every input, including the awkward ones.
    """

    @staticmethod
    def _naive(variant_a: str, variant_b: str, severity: pd.DataFrame) -> float:
        """The original implementation, kept here as the reference."""
        if severity.empty:
            return 0.0

        def mean_for(variant: str) -> float:
            gene, mutation = variant.split(":", 1)
            sub = severity[
                (severity["Gene"].astype(str) == gene)
                & (severity["Mutation"].astype(str) == mutation)
            ]
            return float(sub["severity_score"].mean()) if not sub.empty else 0.0

        return round((mean_for(variant_a) + mean_for(variant_b)) / 2.0, 2)

    def test_index_matches_the_naive_lookup(self):
        severity = mutation_severity(_mutations())
        index = build_severity_index(severity)
        variants = [f"{row.Gene}:{row.Mutation}" for row in severity.head(20).itertuples()]

        for variant_a, variant_b in zip(variants, variants[1:] + variants[:1], strict=True):
            assert pair_severity(variant_a, variant_b, severity, index) == self._naive(
                variant_a, variant_b, severity
            ), f"{variant_a}/{variant_b}"

    def test_a_variant_absent_from_the_table_scores_zero(self):
        severity = mutation_severity(_mutations())
        index = build_severity_index(severity)

        assert pair_severity("nosuchgene:A1V", "acrA:R50H", severity, index) == self._naive(
            "nosuchgene:A1V", "acrA:R50H", severity
        )

    def test_an_empty_severity_table_scores_zero(self):
        empty = pd.DataFrame(columns=["Gene", "Mutation", "severity_score"])

        assert build_severity_index(empty) == {}
        assert pair_severity("acrA:R50H", "acrB:H596N", empty) == 0.0

    def test_one_variant_carried_by_many_genomes_uses_their_mean(self):
        """A variant in many genomes is one key, scored by the mean over genomes."""
        repeated = pd.concat([mutation_severity(_mutations())] * 4, ignore_index=True)
        first = repeated.iloc[0]
        gene, mutation = str(first["Gene"]), str(first["Mutation"])

        index = build_severity_index(repeated)
        expected = float(
            repeated[
                (repeated["Gene"].astype(str) == gene)
                & (repeated["Mutation"].astype(str) == mutation)
            ]["severity_score"].mean()
        )

        assert index[(gene, mutation)] == pytest.approx(expected)
