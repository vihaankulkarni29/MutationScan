"""Snakemake wrapper: extract target proteins from every genome.

Thin by design -- unpack what Snakemake handed us and call the library. All the
logic lives in ``mutation_scan.references`` and ``mutation_scan.extract``.
"""

import logging

from mutation_scan.extract import extract
from mutation_scan.references import discover_targets, prepare_references

logging.basicConfig(level=logging.INFO, format="%(levelname)-7s %(message)s")

targets = discover_targets(snakemake.input.references_dir, snakemake.params.targets)
if not targets:
    raise SystemExit(
        f"No reference proteins found in {snakemake.input.references_dir}. "
        "Add <gene>_WT.faa files; each filename becomes a target."
    )

references = prepare_references(
    snakemake.input.references_dir,
    snakemake.output.refs_dir,
    targets,
    snakemake.params.uniprot_taxid,
)
if references.missing:
    raise SystemExit(
        "No reference available for: "
        + ", ".join(references.missing)
        + ". Supply the missing <gene>_WT.faa file, or set uniprot_taxid in the config."
    )

extract(
    snakemake.input.genomes_dir,
    references.refs_dir,
    snakemake.output.proteins_dir,
    references,
    manifest_csv=snakemake.output.manifest,
    tblastn_binary=snakemake.params.tblastn_binary,
    threads=snakemake.threads,
)
