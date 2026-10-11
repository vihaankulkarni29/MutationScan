"""Tests for run provenance (run_summary.json).

The summary is the only artifact a run leaves behind, so the LoF-integrity fix
must be visible here: a genome whose only call is a LoF must still show up in
the mutated-genome count, while older runs (QC without event columns) are
unchanged.
"""

from __future__ import annotations

import pandas as pd
import pytest

from mutation_scan.config import Config
from mutation_scan.summary import build_run_summary
from mutation_scan.variants import MUTATION_COLUMNS, VARIANT_STATUS_CALLED


def _qc_row(
    accession: str,
    gene: str,
    status: str,
    n_frameshifts: int = 0,
    n_knockouts: int = 0,
) -> dict:
    return {
        "Accession": accession,
        "Gene": gene,
        "Status": status,
        "N_Frameshifts": n_frameshifts,
        "N_Functional_Knockouts": n_knockouts,
    }


def _mutation_row(accession: str, gene: str, mutation: str) -> dict:
    return {
        "Accession": accession,
        "Gene": gene,
        "Mutation": mutation,
        "Ref_AA": "X",
        "Position": 1,
        "Alt_AA": "Y",
        "Identity_pct": 99.0,
        "Coverage_pct": 100.0,
    }


@pytest.fixture
def config(tmp_path):
    return Config.from_mapping({}, base_dir=tmp_path)


def _variants(config, qc=None, mutations=None):
    return build_run_summary(config, qc=qc, mutations=mutations)["variants"]


def test_lof_only_genome_is_counted_as_mutated(config):
    """A genome whose only call is a LoF still counts as mutated."""
    mutations = pd.DataFrame([_mutation_row("g1", "alphaX", "H596N")])
    qc = pd.DataFrame(
        [
            _qc_row("g1", "alphaX", VARIANT_STATUS_CALLED),
            _qc_row("g2", "alphaX", VARIANT_STATUS_CALLED, n_knockouts=1),
        ]
    )
    v = _variants(config, qc, mutations)

    assert v["mutations"] == 1
    assert v["loss_of_function_alleles"] == 1
    assert v["genomes_with_loss_of_function"] == 1
    assert v["genomes_with_any_mutation"] == 2


def test_qc_without_event_columns_is_a_noop(config):
    """Older runs (QC without event columns) are unchanged by the LoF layer."""
    mutations = pd.DataFrame([_mutation_row("g1", "alphaX", "H596N")])
    qc = pd.DataFrame(
        [
            {"Accession": "g1", "Gene": "alphaX", "Status": VARIANT_STATUS_CALLED},
        ]
    )
    v = _variants(config, qc, mutations)

    assert v["loss_of_function_alleles"] == 0
    assert v["genomes_with_loss_of_function"] == 0
    assert v["genomes_with_any_mutation"] == 1


def test_empty_mutations_with_lof_carriers_counts_genomes(config):
    """A genome whose only mutation is a LoF is counted even with no substitutions."""
    mutations = pd.DataFrame(columns=list(MUTATION_COLUMNS))
    qc = pd.DataFrame([_qc_row("g2", "alphaX", VARIANT_STATUS_CALLED, n_knockouts=1)])
    v = _variants(config, qc, mutations)

    assert v["mutations"] == 0
    assert v["loss_of_function_alleles"] == 1
    assert v["genomes_with_loss_of_function"] == 1
    assert v["genomes_with_any_mutation"] == 1


def test_none_mutations_reports_no_lof(config):
    """A dry/partial run with no mutations stream reports 0 LoF and 0 mutated."""
    qc = pd.DataFrame([_qc_row("g2", "alphaX", VARIANT_STATUS_CALLED, n_knockouts=1)])
    v = _variants(config, qc, None)

    assert v["loss_of_function_alleles"] == 0
    assert v["genomes_with_loss_of_function"] == 0
    assert v["genomes_with_any_mutation"] == 0


def test_substitutions_without_lof_count_normally(config):
    """Pure substitution runs report LoF zeros and the usual genome count."""
    mutations = pd.DataFrame(
        [
            _mutation_row("g1", "alphaX", "H596N"),
            _mutation_row("g1", "betaQ", "K12R"),
        ]
    )
    qc = pd.DataFrame(
        [
            _qc_row("g1", "alphaX", VARIANT_STATUS_CALLED),
            _qc_row("g1", "betaQ", VARIANT_STATUS_CALLED),
        ]
    )
    v = _variants(config, qc, mutations)

    assert v["mutations"] == 2
    assert v["loss_of_function_alleles"] == 0
    assert v["genomes_with_loss_of_function"] == 0
    assert v["genomes_with_any_mutation"] == 1
