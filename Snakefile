# Snakefile - MutationScan Master Workflow
configfile: "config/config.yaml"

# ---------------------------------------------------------
# DYNAMIC ROUTING LOGIC
# ---------------------------------------------------------
# If the user provides local genomes, bypass acquisition.
GENOMES_DIR = config.get("local_genomes", "")
if not GENOMES_DIR:
    GENOMES_DIR = "data/results/genomes"

# ---------------------------------------------------------
# MASTER TARGET RULE
# ---------------------------------------------------------
rule all:
    input:
        "data/results/1_genomics_report.csv"

# ---------------------------------------------------------
# PHASE 1a: ACQUISITION (Will be skipped if GENOMES_DIR is local)
# ---------------------------------------------------------
rule acquire_genomes:
    input:
        raw_csv=config["input_csv"]
    output:
        curated_csv="data/results/curated_metadata.csv",
        genomes_dir=directory("data/results/genomes")
    params:
        email=config["email"],
        api_key=config.get("api_key", "")
    script:
        "src/scripts/01_acquire_genomes.py"

# ---------------------------------------------------------
# PHASE 1b: EXTRACTION & CALLING
# ---------------------------------------------------------
rule extract_and_call:
    input:
        genomes_dir=GENOMES_DIR,
        targets_file=config["targets_file"]
    output:
        proteins_dir=directory("data/results/proteins"),
        refs_dir=directory("data/results/refs"),
        mutations_csv="data/results/1_genomics_report.csv"
    params:
        uniprot_taxid=config.get("uniprot_taxid", "")
    script:
        "src/scripts/02_extract_and_call.py"
