# Snakefile - MutationScan Master Workflow
configfile: "config/config.yaml"

import os

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------
# The entry point is now strictly a local directory of provided genomes.
GENOMES_DIR = config.get("local_genomes", "data/local_genomes")
TARGETS_FILE = config.get("targets_file", "config/acr_targets.txt")
DEFAULT_PDB = config.get("default_pdb", "data/5o66.pdb")

# ---------------------------------------------------------------------------
# MASTER RULE
# ---------------------------------------------------------------------------
rule all:
    input:
        "data/results/1_genomics_report.csv",
        "data/results/2_epistasis_networks.csv",
        "data/results/ControlScan_Networks",
        "data/results/3_biophysics_docking.csv"

# ---------------------------------------------------------------------------
# PHASE 1: EXTRACTION & VARIANT CALLING (The New Entry Point)
# ---------------------------------------------------------------------------
rule extract_and_call:
    input:
        genomes=GENOMES_DIR,
        targets=TARGETS_FILE
    output:
        proteins=directory("data/results/proteins"),
        report="data/results/1_genomics_report.csv"
    shell:
        """
        python src/scripts/02_extract_and_call.py \
            --genomes {input.genomes} \
            --targets {input.targets} \
            --output-dir data/results
        """

# ---------------------------------------------------------------------------
# PHASE 2: BIOCHEMICAL EPISTASIS
# ---------------------------------------------------------------------------
rule biochemical_epistasis:
    input:
        report="data/results/1_genomics_report.csv"
    output:
        networks="data/results/2_epistasis_networks.csv",
        plots=directory("data/results/ControlScan_Networks")
    shell:
        """
        python src/scripts/03_biochemical_epistasis.py \
            --report {input.report} \
            --output-dir data/results
        """

# ---------------------------------------------------------------------------
# PHASE 3: OPENMM DYNAMICS & HTVS DOCKING
# ---------------------------------------------------------------------------
rule htvs_biophysics:
    input:
        networks="data/results/2_epistasis_networks.csv",
        proteins="data/results/proteins",
        pdb=DEFAULT_PDB
    output:
        docking_report="data/results/3_biophysics_docking.csv",
        mutated_pdbs=directory("data/results/Mutated_Structures")
    shell:
        """
        python src/scripts/04_htvs_biophysics.py \
            --networks {input.networks} \
            --proteins {input.proteins} \
            --reference-pdb {input.pdb} \
            --output-dir data/results
        """
