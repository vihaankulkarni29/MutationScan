"""Species-aware reference routing.

The regression these guard against is not a crash or a changed column count. It
is a run that completes cleanly and reports every difference between two species'
proteins as an acquired mutation. On the 570-genome three-species cohort that
produced 61,858 Klebsiella mutations from E. coli reference genes while the
equivalent per-species runs produced 414.

Nothing here asserts on a real organism: the fixtures are invented genes and
synthetic sequences, because the tool itself must not know any biology.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from mutation_scan.extract import (
    STATUS_EXTRACTED,
    STATUS_NO_SPECIES_REFERENCE,
    extract,
)
from mutation_scan.references import prepare_species_references
from mutation_scan.species import (
    SPECIES_UNASSIGNED,
    SpeciesConfigError,
    SpeciesPanel,
    parse_species_panel,
)

from ..conftest import (
    PROTEIN_ALPHA,
    PROTEIN_BETA,
    requires_tblastn,
    synthetic_genome,
)

#: Same length, unrelated composition: distinct enough that aligning one against
#: the other cannot score well, which is what makes a cross-species mix-up
#: detectable at all.
PROTEIN_GAMMA = PROTEIN_ALPHA[::-1]


def write_genome(directory: Path, accession: str, proteins: dict[str, str]) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{accession}.fna"
    path.write_text(synthetic_genome(proteins), encoding="utf-8")
    return path


def write_reference(directory: Path, gene: str, protein: str) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / f"{gene}_WT.faa").write_text(
        f">{gene}_WT synthetic\n{protein}\n", encoding="utf-8"
    )


def two_species_panel(root: Path) -> tuple[SpeciesPanel, dict]:
    """A panel with two species, each owning a different gene."""
    panel = parse_species_panel(
        {
            "alpha_species": {
                "genome_pattern": "Alpha_*",
                "references": "alpha_refs",
            },
            "beta_species": {
                "genome_pattern": "Beta_*",
                "references": "beta_refs",
            },
        },
        root,
    )
    write_reference(root / "alpha_refs", "geneA", PROTEIN_ALPHA)
    write_reference(root / "beta_refs", "geneB", PROTEIN_BETA)
    staged = prepare_species_references(root, root / "run_refs", panel)
    return panel, staged


class TestPanelParsing:
    def test_an_absent_block_disables_species_routing(self, tmp_path):
        panel = parse_species_panel(None, tmp_path)

        assert not panel.enabled
        assert panel.profiles == ()

    def test_an_empty_block_disables_species_routing(self, tmp_path):
        assert not parse_species_panel({}, tmp_path).enabled

    def test_a_relative_reference_path_resolves_against_the_root(self, tmp_path):
        panel = parse_species_panel({"s": {"genome_pattern": "X_*", "references": "sub"}}, tmp_path)

        assert panel.profiles[0].references_dir == tmp_path / "sub"

    def test_an_absolute_reference_path_is_left_alone(self, tmp_path):
        elsewhere = tmp_path / "elsewhere"
        panel = parse_species_panel(
            {"s": {"genome_pattern": "X_*", "references": str(elsewhere)}}, tmp_path
        )

        assert panel.profiles[0].references_dir == elsewhere

    @pytest.mark.parametrize(
        "value,message",
        [
            (["not", "a", "mapping"], "must be a mapping"),
            ({"s": "not a mapping"}, "must be a mapping"),
            ({"s": {"references": "r"}}, "genome_pattern"),
            ({"s": {"genome_pattern": "X_*"}}, "references"),
            ({"s": {"genome_pattern": "", "references": "r"}}, "genome_pattern"),
            ({"s": {"genome_pattern": "X_*", "references": ""}}, "references"),
            ({"s": {"genome_pattern": "X_*", "references": "r", "extra": 1}}, "Unknown key"),
        ],
    )
    def test_a_malformed_block_is_rejected(self, tmp_path, value, message):
        with pytest.raises(SpeciesConfigError, match=message):
            parse_species_panel(value, tmp_path)

    def test_the_reserved_unassigned_name_is_refused(self, tmp_path):
        with pytest.raises(SpeciesConfigError, match="reserved"):
            parse_species_panel(
                {SPECIES_UNASSIGNED: {"genome_pattern": "X_*", "references": "r"}}, tmp_path
            )


class TestGenomeAssignment:
    @pytest.fixture
    def panel(self, tmp_path):
        return parse_species_panel(
            {
                "alpha": {"genome_pattern": "Alpha_*", "references": "a"},
                "beta": {"genome_pattern": "Beta_*", "references": "b"},
            },
            tmp_path,
        )

    def test_a_matching_genome_gets_its_species(self, panel):
        assert panel.species_for(Path("Alpha_GCA_1.fna")) == "alpha"

    def test_the_pattern_matches_the_stem_not_the_suffix(self, panel):
        assert panel.species_for(Path("/deep/path/Alpha_x.fna.gz")) == "alpha"

    def test_an_unmatched_genome_is_reported_not_guessed(self, panel):
        assert panel.species_for(Path("Gamma_x.fna")) == SPECIES_UNASSIGNED

    def test_an_overlapping_pattern_resolves_to_the_more_specific_one(self, tmp_path):
        panel = parse_species_panel(
            {
                "broad": {"genome_pattern": "*", "references": "b"},
                "narrow": {"genome_pattern": "Alpha_*", "references": "n"},
            },
            tmp_path,
        )

        assert panel.species_for(Path("Alpha_x.fna")) == "narrow"
        assert panel.species_for(Path("Beta_x.fna")) == "broad"

    def test_grouping_separates_species_and_drops_the_rest(self, panel):
        genomes = [Path("Alpha_1.fna"), Path("Beta_1.fna"), Path("Beta_2.fna")]

        grouped = panel.grouped(genomes)

        assert [p.name for p in grouped["beta"]] == ["Beta_1.fna", "Beta_2.fna"]
        assert [p.name for p in grouped["alpha"]] == ["Alpha_1.fna"]
        assert SPECIES_UNASSIGNED not in grouped

    def test_unassigned_lists_only_the_unmatched(self, panel):
        unassigned = panel.unassigned([Path("Alpha_1.fna"), Path("Gamma_1.fna")])

        assert [p.name for p in unassigned] == ["Gamma_1.fna"]


@requires_tblastn
class TestExtractionRoutesEachGenomeToItsOwnSpecies:
    """The core guard: no genome is aligned against another species' protein."""

    def test_each_genome_gets_its_own_species_panel(self, tmp_path):
        genomes_dir = tmp_path / "genomes"
        write_genome(genomes_dir, "Alpha_one", {"geneA": PROTEIN_ALPHA})
        write_genome(genomes_dir, "Beta_one", {"geneB": PROTEIN_BETA})
        panel, staged = two_species_panel(tmp_path)

        targets = sorted({g for s in staged.values() for g in s.genes})
        manifest = extract(
            genomes_dir,
            tmp_path / "run_refs",
            tmp_path / "proteins",
            targets,
            genomes=sorted(genomes_dir.glob("*.fna")),
            panel=panel,
            species_references=staged,
        )

        by_pair = manifest.set_index(["Accession", "Gene"])
        assert by_pair.loc[("Alpha_one", "geneA"), "Status"] == STATUS_EXTRACTED
        assert by_pair.loc[("Beta_one", "geneB"), "Status"] == STATUS_EXTRACTED
        # And each row records which species it was compared within.
        assert by_pair.loc[("Alpha_one", "geneA"), "Species"] == "alpha_species"
        assert by_pair.loc[("Beta_one", "geneB"), "Species"] == "beta_species"

    def test_a_gene_absent_from_a_species_is_never_borrowed(self, tmp_path):
        """The failure the module exists to prevent.

        ``geneB`` belongs to the beta species only. The alpha genome must be told
        it cannot be compared, not handed the beta protein -- which would align
        well enough to pass an identity gate and then report every difference as
        an acquired mutation.
        """
        genomes_dir = tmp_path / "genomes"
        write_genome(genomes_dir, "Alpha_one", {"geneA": PROTEIN_ALPHA})
        panel, staged = two_species_panel(tmp_path)

        targets = sorted({g for s in staged.values() for g in s.genes})
        manifest = extract(
            genomes_dir,
            tmp_path / "run_refs",
            tmp_path / "proteins",
            targets,
            genomes=sorted(genomes_dir.glob("*.fna")),
            panel=panel,
            species_references=staged,
        )

        borrowed = manifest[(manifest["Accession"] == "Alpha_one") & (manifest["Gene"] == "geneB")]
        assert len(borrowed) == 1
        assert borrowed.iloc[0]["Status"] == STATUS_NO_SPECIES_REFERENCE
        assert "alpha_species" in str(borrowed.iloc[0]["Detail"])
        # Crucially, no protein was written for it.
        assert not (tmp_path / "proteins" / "Alpha_one_geneB.faa").exists()

    def test_a_genome_matching_no_species_is_not_compared_at_all(self, tmp_path):
        genomes_dir = tmp_path / "genomes"
        write_genome(genomes_dir, "Gamma_one", {"geneA": PROTEIN_ALPHA})
        panel, staged = two_species_panel(tmp_path)

        targets = sorted({g for s in staged.values() for g in s.genes})
        manifest = extract(
            genomes_dir,
            tmp_path / "run_refs",
            tmp_path / "proteins",
            targets,
            genomes=sorted(genomes_dir.glob("*.fna")),
            panel=panel,
            species_references=staged,
        )

        assert set(manifest["Status"]) == {STATUS_NO_SPECIES_REFERENCE}
        assert not list((tmp_path / "proteins").glob("*.faa"))

    def test_the_manifest_records_the_species_for_every_row(self, tmp_path):
        genomes_dir = tmp_path / "genomes"
        write_genome(genomes_dir, "Alpha_one", {"geneA": PROTEIN_ALPHA})
        write_genome(genomes_dir, "Beta_one", {"geneB": PROTEIN_BETA})
        panel, staged = two_species_panel(tmp_path)

        targets = sorted({g for s in staged.values() for g in s.genes})
        manifest = extract(
            genomes_dir,
            tmp_path / "run_refs",
            tmp_path / "proteins",
            targets,
            genomes=sorted(genomes_dir.glob("*.fna")),
            panel=panel,
            species_references=staged,
        )

        assert "Species" in manifest.columns
        assert set(manifest["Species"]) <= {"alpha_species", "beta_species"}
        # One row per (genome, gene) regardless of outcome.
        assert len(manifest) == 2 * len(targets)


@requires_tblastn
class TestVariantCallingUsesTheSpeciesReference:
    def test_a_pair_is_compared_against_its_own_species_protein(self, tmp_path):
        from mutation_scan.variants import call_variants

        genomes_dir = tmp_path / "genomes"
        write_genome(genomes_dir, "Alpha_one", {"geneA": PROTEIN_ALPHA})
        write_genome(genomes_dir, "Beta_one", {"geneB": PROTEIN_BETA})
        panel, staged = two_species_panel(tmp_path)

        targets = sorted({g for s in staged.values() for g in s.genes})
        manifest = extract(
            genomes_dir,
            tmp_path / "run_refs",
            tmp_path / "proteins",
            targets,
            genomes=sorted(genomes_dir.glob("*.fna")),
            panel=panel,
            species_references=staged,
        )

        result = call_variants(
            tmp_path / "proteins",
            tmp_path / "run_refs",
            manifest=manifest,
            min_identity_percent=80.0,
            min_coverage_percent=80.0,
            species_refs_dirs={k: v.refs_dir for k, v in staged.items()},
        )

        called = result.qc[result.qc["Status"] == "called"]
        assert len(called) == 2
        # Identical proteins back to their own references: no mutations at all.
        assert result.mutations.empty
        assert set(called["Identity_pct"]) == {100.0}
        assert set(called["Species"]) == {"alpha_species", "beta_species"}

    def test_a_mismatched_species_reference_is_not_used_as_a_substitute(self, tmp_path):
        """Swapping the reference must change the answer, not be silently absorbed."""
        from mutation_scan.variants import call_variants

        genomes_dir = tmp_path / "genomes"
        write_genome(genomes_dir, "Alpha_one", {"geneA": PROTEIN_ALPHA})
        panel, staged = two_species_panel(tmp_path)

        manifest = pd.DataFrame(
            [
                {
                    "Accession": "Alpha_one",
                    "Gene": "geneA",
                    "Status": STATUS_EXTRACTED,
                    "Species": "alpha_species",
                }
            ]
        )
        (tmp_path / "proteins").mkdir(exist_ok=True)
        (tmp_path / "proteins" / "Alpha_one_geneA.faa").write_text(
            f">Alpha_one_geneA\n{PROTEIN_GAMMA}\n", encoding="utf-8"
        )

        result = call_variants(
            tmp_path / "proteins",
            tmp_path / "run_refs",
            manifest=manifest,
            min_identity_percent=80.0,
            min_coverage_percent=80.0,
            species_refs_dirs={k: v.refs_dir for k, v in staged.items()},
        )

        # The reversed protein is not the alpha reference, so it must not be
        # reported as a near-identical match.
        identity = result.qc["Identity_pct"].iloc[0]
        assert identity != 100.0


class TestMultiRecordBundlesAreNotTargets:
    def test_a_combined_fasta_is_skipped_rather_than_failing_the_run(self, tmp_path):
        """Bundles sit next to real references and are not genes."""
        from mutation_scan.references import discover_targets

        refs = tmp_path / "refs"
        write_reference(refs, "geneA", PROTEIN_ALPHA)
        write_reference(refs, "geneB", PROTEIN_BETA)
        (refs / "combined_refs.faa").write_text(
            f">geneA\n{PROTEIN_ALPHA}\n>geneB\n{PROTEIN_BETA}\n", encoding="utf-8"
        )

        assert discover_targets(refs) == ["geneA", "geneB"]

    def test_single_record_references_are_still_discovered(self, tmp_path):
        from mutation_scan.references import discover_targets

        refs = tmp_path / "refs"
        write_reference(refs, "geneA", PROTEIN_ALPHA)

        assert discover_targets(refs) == ["geneA"]
