"""Snakemake wrapper: write run_summary.json.

Thin by design -- everything here is a read of a file an earlier rule already
wrote, handed straight to ``mutation_scan.summary.write_run_summary``. Both entry
points call that one function, so a run is described identically whether it was
launched with ``mutationscan run`` or with Snakemake.
"""

import logging

from mutation_scan.config import Config
from mutation_scan.cooccurrence import load_result
from mutation_scan.extract import read_manifest, verify_tblastn
from mutation_scan.references import load_reference_set
from mutation_scan.summary import write_run_summary
from mutation_scan.variants import read_mutations, read_variant_qc

logging.basicConfig(level=logging.INFO, format="%(levelname)-7s %(message)s")

config = Config.from_mapping(snakemake.config)
qc = read_variant_qc(snakemake.input.qc)

write_run_summary(
    config,
    path=snakemake.output.summary,
    references=load_reference_set(snakemake.input.refs_dir),
    manifest=read_manifest(snakemake.input.manifest),
    qc=qc,
    mutations=read_mutations(snakemake.input.mutations),
    cooccurrence=load_result(
        snakemake.input.genes,
        snakemake.input.matrix,
        snakemake.params.variants or None,
        qc=qc,
    ),
    tblastn_version=verify_tblastn(snakemake.params.tblastn_binary),
)