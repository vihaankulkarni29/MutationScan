# MutationScan workflow.
#
#     snakemake --cores 4                       # uses config/config.yaml
#     snakemake --cores 4 --config job_name=x   # override any config key
#     snakemake -n                              # dry run: check the DAG
#
# Three stages, in order: extract proteins from genomes, call substitutions
# against the reference, count co-occurrence. Every path below is derived from a
# single validated Config object, so this file and `mutationscan run` cannot
# drift apart -- they call the same functions with the same layout.
#
# Nothing here names a gene, organism or database. Targets are whatever
# reference files are present in references_dir.

from mutation_scan.config import Config

configfile: "config/config.yaml"

# Snakemake runs with the workdir as cwd, so relative config paths resolve there.
CONFIG = Config.from_mapping(config)

WANT_VARIANT_LEVEL = "variant" in CONFIG.cooccurrence.levels


rule all:
    input:
        [str(path) for path in CONFIG.outputs()],


rule extract_proteins:
    """tblastn every reference protein against every genome."""
    input:
        genomes_dir=str(CONFIG.genomes_dir),
        references_dir=str(CONFIG.references_dir),
    output:
        manifest=str(CONFIG.manifest_csv),
        proteins_dir=directory(str(CONFIG.proteins_dir)),
        refs_dir=directory(str(CONFIG.refs_dir)),
    params:
        targets=list(CONFIG.targets),
        uniprot_taxid=CONFIG.uniprot_taxid,
        tblastn_binary=CONFIG.tblastn_binary,
    threads: CONFIG.threads
    script:
        "workflow/scripts/01_extract.py"


rule call_variants:
    """Global-align each extracted protein to its reference and call substitutions."""
    input:
        manifest=str(CONFIG.manifest_csv),
        proteins_dir=str(CONFIG.proteins_dir),
        refs_dir=str(CONFIG.refs_dir),
    output:
        mutations=str(CONFIG.mutations_csv),
        qc=str(CONFIG.variant_qc_csv),
    params:
        min_identity_percent=CONFIG.min_identity_percent,
        min_coverage_percent=CONFIG.min_coverage_percent,
    script:
        "workflow/scripts/02_call_variants.py"


rule cooccurrence:
    """Count how often pairs of genes are mutated in the same genome."""
    input:
        mutations=str(CONFIG.mutations_csv),
        qc=str(CONFIG.variant_qc_csv),
    output:
        genes=str(CONFIG.cooccurrence_genes_csv),
        matrix=str(CONFIG.cooccurrence_matrix_csv),
        **(
            {"variants": str(CONFIG.cooccurrence_variants_csv)}
            if WANT_VARIANT_LEVEL
            else {}
        ),
    params:
        levels=list(CONFIG.cooccurrence.levels),
        min_count=CONFIG.cooccurrence.min_count,
    script:
        "workflow/scripts/03_cooccurrence.py"


rule run_summary:
    """Record what was run, on what, with which settings, and what came out."""
    input:
        manifest=str(CONFIG.manifest_csv),
        mutations=str(CONFIG.mutations_csv),
        qc=str(CONFIG.variant_qc_csv),
        genes=str(CONFIG.cooccurrence_genes_csv),
        matrix=str(CONFIG.cooccurrence_matrix_csv),
        refs_dir=str(CONFIG.refs_dir),
    output:
        summary=str(CONFIG.run_summary_json),
    params:
        variants=(str(CONFIG.cooccurrence_variants_csv) if WANT_VARIANT_LEVEL else ""),
        tblastn_binary=CONFIG.tblastn_binary,
    script:
        "workflow/scripts/04_summary.py"
