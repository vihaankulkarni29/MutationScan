# Snakefile - MutationScan Master Workflow
configfile: "config/config.yaml"

import os

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------
GENOMES_DIR = config.get("local_genomes", "data/local_genomes")
TARGETS_FILE = config.get("targets_file", "config/targets.txt")
JOB_NAME = config.get("job_name", "default_run")
SKIP_EXTRACTION = str(config.get("skip_extraction", False)).lower() in ("1", "true", "yes", "y")
LEGACY_RESULTS_DIR = config.get("legacy_results_dir", "data/results")

OUT_DIR = f"data/output/{JOB_NAME}"

PROTEINS_INPUT_DIR = config.get(
    "proteins_dir",
    f"{LEGACY_RESULTS_DIR}/proteins" if SKIP_EXTRACTION else f"{OUT_DIR}/proteins"
)
REFS_INPUT_DIR = config.get(
    "refs_dir",
    f"{LEGACY_RESULTS_DIR}/refs" if SKIP_EXTRACTION else f"{OUT_DIR}/refs"
)

# ---------------------------------------------------------------------------
# MASTER RULE
# ---------------------------------------------------------------------------
rule all:
    input:
        f"{OUT_DIR}/1_genomics_report.csv",
        f"{OUT_DIR}/2_epistasis_networks.csv",
        f"{OUT_DIR}/ControlScan_Networks",

# ---------------------------------------------------------------------------
# PHASE 1A: PROTEIN EXTRACTION
# ---------------------------------------------------------------------------
rule extract_proteins:
    input:
        genomes_dir=GENOMES_DIR,
        targets_file=TARGETS_FILE
    output:
        proteins_dir=directory(f"{OUT_DIR}/proteins"),
        refs_dir=directory(f"{OUT_DIR}/refs"),
        marker=f"{OUT_DIR}/proteins/.proteins_extracted"
    params:
        uniprot_taxid=config.get("uniprot_taxid", ""),
        reference_seed_dir=config.get("reference_seed_dir", ""),
        out_dir=OUT_DIR,
        skip_extraction=config.get("skip_extraction", False)
    script:
        "src/scripts/02a_extract_proteins.py"

# ---------------------------------------------------------------------------
# PHASE 1B: VARIANT CALLING
# ---------------------------------------------------------------------------
rule call_variants:
    input:
        proteins_dir=PROTEINS_INPUT_DIR,
        refs_dir=REFS_INPUT_DIR
    output:
        report=f"{OUT_DIR}/1_genomics_report.csv",
        marker=f"{OUT_DIR}/.variants_called"
    params:
        out_dir=OUT_DIR
    script:
        "src/scripts/02b_call_variants.py"

# ---------------------------------------------------------------------------
# PHASE 2: BIOCHEMICAL EPISTASIS
# ---------------------------------------------------------------------------
rule biochemical_epistasis:
    input:
        report=f"{OUT_DIR}/1_genomics_report.csv"
    output:
        networks=f"{OUT_DIR}/2_epistasis_networks.csv",
        plots_dir=directory(f"{OUT_DIR}/ControlScan_Networks")
    params:
        out_dir=OUT_DIR
    script:
        "src/scripts/03_biochemical_epistasis.py"
