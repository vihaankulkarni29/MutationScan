"""End-to-end tests: the whole tool, on a cohort we planted ourselves.

The cohort is built from invented genes (``alphaX``, ``betaQ``) and synthetic
sequences, so these tests prove the pipeline works on *whatever* the user
supplies rather than on one organism's data. Every number asserted below is
derivable by hand from the four genomes constructed in :func:`cohort`:

===========  ================  ===============  ==================================
genome       alphaX            betaQ            note
===========  ================  ===============  ==================================
``GEN001``   wild type         wild type        evaluated, nothing mutated
``GEN002``   ``Y8L``           wild type        one gene mutated
``GEN003``   ``Y8L``           ``T5A``          both genes mutated -- the one
                                                genome behind ``N_Both``
``GEN004``   absent            absent           no target genes at all
===========  ================  ===============  ==================================

So for the pair (alphaX, betaQ): ``N_Both=1``, ``N_A=2``, ``N_B=1`` and
``N_Eligible=3`` -- GEN004 contributes to no denominator because neither gene was
ever evaluated there.
"""

from __future__ import annotations

import json

import pytest
import yaml

from mutation_scan import cli
from mutation_scan.config import Config
from mutation_scan.cooccurrence import read_cooccurrence
from mutation_scan.extract import STATUS_EXTRACTED, STATUS_NO_HIT, read_manifest
from mutation_scan.variants import (
    VARIANT_STATUS_CALLED,
    read_mutations,
    read_variant_qc,
)

from ..conftest import (
    PROTEIN_ALPHA,
    PROTEIN_BETA,
    mutate,
    requires_tblastn,
    synthetic_genome,
)

#: The substitutions planted in the cohort, in reference coordinates.
ALPHA_MUTATION = "Y8L"  # PROTEIN_ALPHA position 8 is Y
BETA_MUTATION = "T5A"  # PROTEIN_BETA position 5 is T


def _write_cohort(root):
    """Lay out references and genomes; return (references_dir, genomes_dir)."""
    references = root / "references"
    genomes = root / "genomes"
    references.mkdir(parents=True)
    genomes.mkdir(parents=True)

    (references / "alphaX_WT.faa").write_text(f">alphaX\n{PROTEIN_ALPHA}\n", encoding="utf-8")
    (references / "betaQ_WT.faa").write_text(f">betaQ\n{PROTEIN_BETA}\n", encoding="utf-8")

    alpha_mutant = mutate(PROTEIN_ALPHA, (8, "L"))
    beta_mutant = mutate(PROTEIN_BETA, (5, "A"))

    (genomes / "GEN001.fna").write_text(
        synthetic_genome({"alphaX": PROTEIN_ALPHA, "betaQ": PROTEIN_BETA}),
        encoding="utf-8",
    )
    (genomes / "GEN002.fna").write_text(
        synthetic_genome({"alphaX": alpha_mutant, "betaQ": PROTEIN_BETA}),
        encoding="utf-8",
    )
    (genomes / "GEN003.fna").write_text(
        synthetic_genome({"alphaX": alpha_mutant, "betaQ": beta_mutant}),
        encoding="utf-8",
    )
    # A genome carrying neither target: it must be recorded as searched-and-empty,
    # never silently dropped, and must not land in anyone's denominator.
    (genomes / "GEN004.fna").write_text(
        ">contig_empty\n" + "ATGCGCTTAAGGCCTTAAGGCATCGATCGGATCC" * 20 + "\n",
        encoding="utf-8",
    )
    return references, genomes


@pytest.fixture(scope="module")
def completed_run(tmp_path_factory):
    """Run the whole tool once, through the real command line, and return its Config.

    Module-scoped because tblastn is the slow part: every assertion below reads
    the outputs of this single run.
    """
    root = tmp_path_factory.mktemp("cohort")
    references, genomes = _write_cohort(root)
    out = root / "output"

    exit_code = cli.main(
        [
            "run",
            "--genomes",
            str(genomes),
            "--references",
            str(references),
            "--out",
            str(out),
            "--job-name",
            "integration",
            "--cooccurrence-level",
            "gene,variant",
            "--threads",
            "2",
        ]
    )
    assert exit_code == 0

    return Config.from_mapping(
        {
            "job_name": "integration",
            "genomes_dir": str(genomes),
            "references_dir": str(references),
            "output_root": str(out),
            "cooccurrence": {"levels": ["gene", "variant"], "min_count": 1},
        }
    )


@requires_tblastn
def test_every_declared_output_exists_and_reloads(completed_run):
    """config.outputs() is a promise; each file must be there and be readable."""
    for path in completed_run.outputs():
        assert path.is_file(), f"missing declared output: {path}"

    read_manifest(completed_run.manifest_csv)
    read_mutations(completed_run.mutations_csv)
    read_variant_qc(completed_run.variant_qc_csv)
    read_cooccurrence(completed_run.cooccurrence_genes_csv)
    read_cooccurrence(completed_run.cooccurrence_variants_csv)
    json.loads(completed_run.run_summary_json.read_text(encoding="utf-8"))


@requires_tblastn
def test_the_targets_came_from_the_reference_filenames(completed_run):
    manifest = read_manifest(completed_run.manifest_csv)
    assert set(manifest["Gene"]) == {"alphaX", "betaQ"}
    # Nothing was configured: two files on disk produced two targets.
    assert completed_run.targets == ()


@requires_tblastn
def test_a_genome_without_the_targets_is_recorded_not_dropped(completed_run):
    manifest = read_manifest(completed_run.manifest_csv)
    gen004 = manifest[manifest["Accession"] == "GEN004"]

    assert len(gen004) == 2, "both targets should have been attempted"
    assert set(gen004["Status"]) == {STATUS_NO_HIT}

    # Attempted-and-empty is not the same as evaluated: it must not reach QC.
    qc = read_variant_qc(completed_run.variant_qc_csv)
    assert "GEN004" not in set(qc["Accession"])


@requires_tblastn
def test_planted_substitutions_are_called_at_their_reference_positions(completed_run):
    mutations = read_mutations(completed_run.mutations_csv)
    found = {(row.Accession, row.Gene, row.Mutation) for row in mutations.itertuples()}

    assert ("GEN002", "alphaX", ALPHA_MUTATION) in found
    assert ("GEN003", "alphaX", ALPHA_MUTATION) in found
    assert ("GEN003", "betaQ", BETA_MUTATION) in found

    # The wild-type genome is evaluated and clean -- absent here, present in QC.
    assert "GEN001" not in set(mutations["Accession"])
    qc = read_variant_qc(completed_run.variant_qc_csv)
    gen001 = qc[qc["Accession"] == "GEN001"]
    assert set(gen001["Status"]) == {VARIANT_STATUS_CALLED}


@requires_tblastn
def test_wild_type_genomes_pass_both_gates_at_full_identity(completed_run):
    qc = read_variant_qc(completed_run.variant_qc_csv)
    gen001 = qc[qc["Accession"] == "GEN001"]

    assert (gen001["Identity_pct"] == 100.0).all()
    assert (gen001["Coverage_pct"] == 100.0).all()
    assert (gen001["N_Mutations"] == 0).all()


@requires_tblastn
def test_gene_pair_counts_match_the_planted_cohort(completed_run):
    """The whole point of the tool, checked against numbers we can count by hand."""
    genes = read_cooccurrence(completed_run.cooccurrence_genes_csv)
    assert len(genes) == 1  # two targets -> exactly one unordered pair

    row = genes.iloc[0]
    assert {row["Gene_A"], row["Gene_B"]} == {"alphaX", "betaQ"}
    assert int(row["N_Eligible"]) == 3  # GEN001-3 evaluated for both; GEN004 not
    assert int(row["N_Both"]) == 1  # only GEN003 carries both

    # Orientation-independent: alphaX is mutated twice, betaQ once.
    counts = {row["Gene_A"]: int(row["N_A"]), row["Gene_B"]: int(row["N_B"])}
    assert counts == {"alphaX": 2, "betaQ": 1}


@requires_tblastn
def test_the_matrix_diagonal_counts_mutated_genomes(completed_run):
    matrix = read_cooccurrence(completed_run.cooccurrence_matrix_csv)
    matrix = matrix.set_index(matrix.columns[0])

    assert int(matrix.loc["alphaX", "alphaX"]) == 2
    assert int(matrix.loc["betaQ", "betaQ"]) == 1
    assert int(matrix.loc["alphaX", "betaQ"]) == 1
    assert int(matrix.loc["betaQ", "alphaX"]) == 1


@requires_tblastn
def test_variant_level_pairs_use_the_gene_level_denominator(completed_run):
    variants = read_cooccurrence(completed_run.cooccurrence_variants_csv)
    labels = {frozenset((row.Variant_A, row.Variant_B)): row for row in variants.itertuples()}

    pair = labels[frozenset((f"alphaX:{ALPHA_MUTATION}", f"betaQ:{BETA_MUTATION}"))]
    assert int(pair.N_Both) == 1
    assert int(pair.N_Eligible) == 3  # from the genes, not from the mutations


@requires_tblastn
def test_the_run_summary_records_provenance_without_inference(completed_run):
    summary = json.loads(completed_run.run_summary_json.read_text(encoding="utf-8"))

    assert summary["mutation_scan_version"]
    assert summary["environment"]["tblastn"]
    assert summary["config"]["job_name"] == "integration"
    assert sorted(summary["references"]["genes"]) == ["alphaX", "betaQ"]
    assert summary["references"]["missing"] == []
    assert summary["references"]["lengths"]["alphaX"] == len(PROTEIN_ALPHA)

    assert summary["extraction"]["genomes"] == 4
    assert summary["extraction"]["status_counts"][STATUS_EXTRACTED] == 6
    assert summary["extraction"]["status_counts"][STATUS_NO_HIT] == 2
    assert summary["variants"]["pairs_called"] == 6
    assert summary["cooccurrence"]["mutated_per_gene"] == {"alphaX": 2, "betaQ": 1}
    assert summary["cooccurrence"]["eligible_per_gene"] == {"alphaX": 3, "betaQ": 3}

    # Co-occurrence is descriptive only: no p-values, odds, enrichment scores
    # or severities in the cooccurrence section of the provenance record.
    cooccurrence_text = json.dumps(summary["cooccurrence"]).lower()
    for forbidden in ("p_value", "pvalue", "odds", "enrich", "severity", "significan"):
        assert forbidden not in cooccurrence_text


@requires_tblastn
def test_the_run_is_reproducible_when_repeated(completed_run, tmp_path):
    """Same inputs, same settings, byte-identical tables."""
    second = tmp_path / "again"
    exit_code = cli.main(
        [
            "run",
            "--genomes",
            str(completed_run.genomes_dir),
            "--references",
            str(completed_run.references_dir),
            "--out",
            str(second),
            "--job-name",
            "integration",
            "--cooccurrence-level",
            "gene,variant",
        ]
    )
    assert exit_code == 0

    for name in (
        "extraction_manifest.csv",
        "mutations.csv",
        "variant_qc.csv",
        "cooccurrence_genes.csv",
        "cooccurrence_genes_matrix.csv",
        "cooccurrence_variants.csv",
    ):
        first_text = (completed_run.run_dir / name).read_text(encoding="utf-8")
        second_text = (second / "integration" / name).read_text(encoding="utf-8")
        assert first_text == second_text, f"{name} differs between identical runs"


# -- tests that need no BLAST+ ----------------------------------------------


def _config_file(tmp_path, **overrides):
    references, genomes = _write_cohort(tmp_path)
    settings = {
        "job_name": "layout_check",
        "genomes_dir": str(genomes),
        "references_dir": str(references),
        "output_root": str(tmp_path / "output"),
    }
    settings.update(overrides)
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(settings), encoding="utf-8")
    return path


def test_config_check_prints_the_resolved_layout(tmp_path, capsys):
    path = _config_file(tmp_path)

    assert cli.main(["config-check", "--config", str(path), "--json"]) == 0

    layout = json.loads(capsys.readouterr().out)
    assert layout["inputs_present"] == {"genomes_dir": True, "references_dir": True}
    assert layout["run_dir"].endswith("layout_check")
    assert any(name.endswith("run_summary.json") for name in layout["expected_outputs"])
    assert not any(
        name.endswith("cooccurrence_variants.csv") for name in layout["expected_outputs"]
    )


def test_config_check_names_the_targets_it_discovered(tmp_path, capsys):
    path = _config_file(tmp_path)

    assert cli.main(["config-check", "--config", str(path)]) == 0

    out = capsys.readouterr().out
    assert "targets (2)" in out
    assert "alphaX" in out and "betaQ" in out


def test_config_check_fails_on_a_missing_input_directory(tmp_path):
    path = _config_file(tmp_path, genomes_dir=str(tmp_path / "nowhere"))

    assert cli.main(["config-check", "--config", str(path)]) == 1


def test_variant_level_output_is_declared_only_when_requested(tmp_path):
    gene_only = Config.from_mapping({"cooccurrence": {"levels": ["gene"]}})
    both = Config.from_mapping({"cooccurrence": {"levels": ["gene", "variant"]}})

    assert gene_only.cooccurrence_variants_csv not in gene_only.outputs()
    assert both.cooccurrence_variants_csv in both.outputs()


def test_a_references_directory_with_no_proteins_is_a_clean_failure(tmp_path):
    (tmp_path / "empty_refs").mkdir()
    (tmp_path / "genomes").mkdir()

    exit_code = cli.main(
        [
            "run",
            "--genomes",
            str(tmp_path / "genomes"),
            "--references",
            str(tmp_path / "empty_refs"),
            "--out",
            str(tmp_path / "output"),
        ]
    )
    assert exit_code == 2  # told what to do, not a traceback


def test_a_missing_genomes_directory_is_a_clean_failure(tmp_path):
    (tmp_path / "refs").mkdir()

    exit_code = cli.main(
        [
            "run",
            "--genomes",
            str(tmp_path / "nowhere"),
            "--references",
            str(tmp_path / "refs"),
            "--out",
            str(tmp_path / "output"),
        ]
    )
    assert exit_code == 2
