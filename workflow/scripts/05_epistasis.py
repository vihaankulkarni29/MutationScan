"""Snakemake wrapper: epistasis analysis.

Thin by design -- unpack what Snakemake handed us and call the library. Reads
the mutations table and QC from earlier stages, writes the genome map,
per-mutation severity, and ranked epistasis networks.
"""

import logging

from mutation_scan.epistasis import epistasis_networks

logging.basicConfig(level=logging.INFO, format="%(levelname)-7s %(message)s")

epistasis_networks(
    mutations_csv=snakemake.input.mutations,
    qc_csv=snakemake.input.qc,
    output_dir=snakemake.output.dir,
    fdr_threshold=snakemake.params.fdr_threshold,
    min_count=snakemake.params.min_count,
    threads=snakemake.threads,
)
