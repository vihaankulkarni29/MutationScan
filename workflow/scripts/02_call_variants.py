"""Snakemake wrapper: call substitutions against the wild-type reference.

Thin by design -- unpack what Snakemake handed us and call the library. The
alignment and the gates live in ``mutation_scan.variants``.
"""

import logging

from mutation_scan.extract import read_manifest
from mutation_scan.variants import call_variants

logging.basicConfig(level=logging.INFO, format="%(levelname)-7s %(message)s")

call_variants(
    snakemake.input.proteins_dir,
    snakemake.input.refs_dir,
    manifest=read_manifest(snakemake.input.manifest),
    min_identity_percent=snakemake.params.min_identity_percent,
    min_coverage_percent=snakemake.params.min_coverage_percent,
    mutations_csv=snakemake.output.mutations,
    qc_csv=snakemake.output.qc,
)
