# MutationScan

MutationScan is a Snakemake-orchestrated AMR analytics pipeline that transforms local bacterial genome assemblies into:

1. Mutation call reports
2. Biochemical epistasis network rankings
3. Optional structure-guided docking deltas (WT vs mutant)

The repository is structured for production use with deterministic workflow steps, job-scoped output directories, and strict separation of source code vs runtime state.

## What MutationScan Does

MutationScan executes a staged workflow:

1. Sequence extraction and variant calling from local `.fna` genomes
2. Biochemical scoring and co-occurrence epistasis network generation
3. Optional biophysics docking against a provided protein structure

Current design principle:

- Local genomes are the input source (no built-in metadata download stage in the production DAG).
- Every run is namespaced by `job_name` and writes to `data/output/{job_name}/`.

## Recent Updates (March 2026)

The current main branch includes several pipeline correctness and quality upgrades:

- Variant-calling identity filter to suppress weak-homology mutation inflation.
- MVBM docking refinements with fixed-pocket targeting and flexible-residue mutant docking.
- Fast steric quality control with explicit `FAILED_QC` status for non-physical mutant models.
- Confidence and interpretation annotations in biophysics outputs for easier triage.

These updates are now the documented baseline behavior for new runs.

## Production Workflow (Snakemake)

The active workflow in [Snakefile](Snakefile) calls exactly these scripts:

1. [src/scripts/02a_extract_proteins.py](src/scripts/02a_extract_proteins.py)
2. [src/scripts/02b_call_variants.py](src/scripts/02b_call_variants.py)
3. [src/scripts/03_biochemical_epistasis.py](src/scripts/03_biochemical_epistasis.py)
4. [src/scripts/04_htvs_biophysics.py](src/scripts/04_htvs_biophysics.py)

No legacy acquisition script is used in the current production DAG.

## Inputs and Outputs

Required inputs:

- Local genomes directory (default `data/local_genomes`)
- Target gene list (default `config/acr_targets.txt`)
- Optional reference PDB for Phase 3 biophysics (default `data/5o66.pdb` in config)
- Optional ligand path from config (`ligand`)

Primary outputs for a run:

- `data/output/{job_name}/1_genomics_report.csv`
- `data/output/{job_name}/2_epistasis_networks.csv`
- `data/output/{job_name}/ControlScan_Networks/`
- `data/output/{job_name}/3_biophysics_docking.csv`
- `data/output/{job_name}/Mutated_Structures/`
- `data/output/{job_name}/README_Biophysics.txt`

## Configuration

Edit [config/config.yaml](config/config.yaml) to control run behavior.

Minimum important keys:

- `job_name`: output namespace for this run
- `local_genomes`: folder containing `.fna` files
- `targets_file`: target genes list
- `variant_min_identity_percent`: minimum alignment identity threshold for variant emission (default `80`)
- `default_pdb`: structure file for biophysics stage
- `ligand`: optional ligand file path for docking
- `pocket_center_x`/`pocket_center_y`/`pocket_center_z`: optional override for docking pocket center (default AcrB center)
- `exhaustiveness`: docking search exhaustiveness (default `16`)

Example:

```yaml
job_name: "trial_001"
local_genomes: "data/local_genomes"
targets_file: "config/acr_targets.txt"
variant_min_identity_percent: 80
default_pdb: "data/5o66.pdb"
ligand: "data/ligands/ligand.sdf"
exhaustiveness: 16
```

Identity filtering note:

- Alignments below `variant_min_identity_percent` are skipped before mutation emission.
- If you want a broader but noisier search, reduce to `75`; for stricter calls, keep `80` or raise it.

## Quick Start

### Option A: Local environment

Use the project Conda environment definition:

```bash
conda env create -f environment.yml
conda activate mutationscan
pip install -e .
```

Dry-run the DAG:

```bash
python -m snakemake -n --cores 1 --config job_name="smoke_test"
```

Run the workflow:

```bash
python -m snakemake --cores 4 --config job_name="run_2026_03_18"
```

### Option B: Docker

```bash
docker compose build
docker compose run --rm mutationscan python -m snakemake -n --cores 1 --config job_name="docker_smoke"
docker compose run --rm mutationscan python -m snakemake --cores 4 --config job_name="docker_run"
```

## CI/CD Notes

Repository CI validates:

- Unit tests
- Snakemake DAG buildability

Runtime data/state folders are intentionally quarantined via ignore rules, and `.snakemake/` is not tracked.

## Scientific and Operational Disclaimers

This pipeline is intended for research and engineering triage workflows.

1. Not a clinical diagnostic device.
2. Mutation-to-phenotype inference is model- and rule-dependent, not ground truth.
3. Docking outputs are best-effort relative estimates, not absolute binding free-energy truth.
4. Fast local docking does not fully model large conformational changes, explicit solvent, long-timescale dynamics, or complete thermodynamic integration.
5. For high-confidence mechanistic conclusions, use full molecular dynamics and dedicated free-energy methods.

## Repository Hygiene Policy

Tracked assets should remain source/config/documentation only.

Not shipped as production code or tracked outputs:

- `.snakemake/` runtime state
- Generated output under `data/output/*`
- Downloaded genome payloads under `data/local_genomes/*`
- Ad hoc local experiment files

Keep placeholders only (`.gitkeep`) in runtime data folders.

## Troubleshooting

Common causes of failed runs:

1. Missing `.fna` files in `local_genomes`
2. Missing/incorrect target genes file
3. Missing PDB when biophysics stage is enabled
4. Missing external binaries in local environment (`tblastn`, docking dependencies)

Recommended first check:

```bash
python -m snakemake -n --cores 1 --config job_name="debug_run"
```

## License

See [LICENSE](LICENSE).
