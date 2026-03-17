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
JOB_NAME = config.get("job_name", "default_run")
OUT_DIR = f"data/output/{JOB_NAME}"

# ---------------------------------------------------------------------------
# MASTER RULE
# ---------------------------------------------------------------------------
rule all:
    input:
        f"{OUT_DIR}/1_genomics_report.csv",
        f"{OUT_DIR}/2_epistasis_networks.csv",
        f"{OUT_DIR}/ControlScan_Networks",
        f"{OUT_DIR}/3_biophysics_docking.csv"

# ---------------------------------------------------------------------------
# PHASE 1: EXTRACTION & VARIANT CALLING (The New Entry Point)
# ---------------------------------------------------------------------------
rule extract_and_call:
    input:
        genomes_dir=GENOMES_DIR,
        targets_file=TARGETS_FILE
    output:
        proteins_dir=directory(f"{OUT_DIR}/proteins"),
        report=f"{OUT_DIR}/1_genomics_report.csv"
    params:
        refs_dir=f"{OUT_DIR}/refs",
        uniprot_taxid=config.get("uniprot_taxid", ""),
        out_dir=OUT_DIR
    script:
        "src/scripts/02_extract_and_call.py"

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

# ---------------------------------------------------------------------------
# PHASE 3: OPENMM DYNAMICS & HTVS DOCKING
# ---------------------------------------------------------------------------
rule htvs_biophysics:
    input:
        networks=f"{OUT_DIR}/2_epistasis_networks.csv",
        proteins_dir=f"{OUT_DIR}/proteins",
        pdb_file=DEFAULT_PDB
    output:
        docking_report=f"{OUT_DIR}/3_biophysics_docking.csv",
        mutated_pdbs=directory(f"{OUT_DIR}/Mutated_Structures"),
        readme=f"{OUT_DIR}/README_Biophysics.txt"
    params:
        pdb=DEFAULT_PDB,
        chain_map=config.get("chain_map", ""),
        ligand=config.get("ligand", ""),
        center_x=config.get("center_x", 0.0),
        center_y=config.get("center_y", 0.0),
        center_z=config.get("center_z", 0.0),
        stiffness=config.get("md_stiffness", 500.0),
        out_dir=OUT_DIR
    script:
        "src/scripts/04_htvs_biophysics.py"
