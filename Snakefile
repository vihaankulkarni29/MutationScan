# Snakefile - MutationScan Master Workflow
configfile: "config/config.yaml"

# ---------------------------------------------------------
# DYNAMIC ROUTING & EXECUTION LOGIC
# ---------------------------------------------------------
# If the user provides local genomes, bypass acquisition.
GENOMES_DIR = config.get("local_genomes", "")
if not GENOMES_DIR:
    GENOMES_DIR = "data/results/genomes"

# Check if user wants to trigger Phase 4
RUN_BIOPHYSICS = bool(config.get("default_pdb", ""))

# ---------------------------------------------------------
# MASTER TARGET RULE
# ---------------------------------------------------------
TARGETS = [
    "data/results/1_genomics_report.csv",
    "data/results/2_epistasis_networks.csv",
    "data/results/ControlScan_Networks",
]

if RUN_BIOPHYSICS:
    TARGETS.extend([
        "data/results/3_biophysics_docking.csv",
        "data/results/Mutated_Structures",
        "data/results/README_Biophysics.txt",
    ])


rule all:
    input: TARGETS

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
        api_key=config.get("api_key", ""),
        filter_country=config.get("filter_country", "india"),
        filter_min_year=config.get("filter_min_year", 2015)
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
        mutations_csv="data/results/1_genomics_report.csv"
    params:
        refs_dir="data/results/refs",
        uniprot_taxid=config.get("uniprot_taxid", "")
    script:
        "src/scripts/02_extract_and_call.py"

# ---------------------------------------------------------
# PHASE 2 & 3: CONTROLSCAN & EPISTASIS
# ---------------------------------------------------------
rule biochemical_epistasis:
    input:
        mutations_csv="data/results/1_genomics_report.csv"
    output:
        epistasis_csv="data/results/2_epistasis_networks.csv",
        networks_dir=directory("data/results/ControlScan_Networks")
    script:
        "src/scripts/03_biochemical_epistasis.py"


# ---------------------------------------------------------
# PHASE 4: HTVS BIOPHYSICS (OPTIONAL)
# ---------------------------------------------------------
rule htvs_biophysics:
    input:
        epistasis_csv="data/results/2_epistasis_networks.csv"
    output:
        docking_csv="data/results/3_biophysics_docking.csv",
        mutated_dir=directory("data/results/Mutated_Structures"),
        readme="data/results/README_Biophysics.txt"
    params:
        pdb=config.get("default_pdb", ""),
        chain_map=config.get("chain_map", ""),
        ligand=config.get("ligand", ""),
        center_x=config.get("center_x", 0.0),
        center_y=config.get("center_y", 0.0),
        center_z=config.get("center_z", 0.0),
        stiffness=config.get("md_stiffness", 500.0)
    script:
        "src/scripts/04_htvs_biophysics.py"
