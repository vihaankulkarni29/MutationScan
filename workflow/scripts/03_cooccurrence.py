"""Snakemake wrapper: count how often pairs of genes are mutated together.

Thin by design -- unpack what Snakemake handed us and call the library. The QC
table is passed in because it, not the manifest, defines each pair's denominator.
"""

import logging

from mutation_scan.cooccurrence import cooccurrence
from mutation_scan.variants import read_mutations, read_variant_qc

logging.basicConfig(level=logging.INFO, format="%(levelname)-7s %(message)s")

cooccurrence(
    read_mutations(snakemake.input.mutations),
    read_variant_qc(snakemake.input.qc),
    levels=snakemake.params.levels,
    min_count=snakemake.params.min_count,
    genes_csv=snakemake.output.genes,
    matrix_csv=snakemake.output.matrix,
    variants_csv=snakemake.output.get("variants"),
)
