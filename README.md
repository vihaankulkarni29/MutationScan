# MutationScan v2.1 ??
**The Deterministic Genotype-to-Phenotype-to-Biophysics Engine for AMR Triage**

MutationScan is an industry-grade, Snakemake-orchestrated bioinformatics pipeline designed to analyze Antimicrobial Resistance (AMR). It translates clinical genomic assemblies ($1D$) into structural thermodynamic data ($3D$) by identifying mutations, scoring their evolutionary and biochemical severity, mapping their epistatic networks, and performing High-Throughput Virtual Screening (HTVS) biophysics to calculate drug binding affinity loss ($\Delta\Delta G$).

Unlike black-box machine learning predictors, MutationScan is a **purely deterministic mathematical and physical engine**.

---

## ?? Scientific Rigor & Methodology

MutationScan operates on four core scientific pillars:

1. **Frameshift-Proof Extraction:** Bacterial assemblies are notoriously messy. MutationScan autonomously fetches wild-type references from UniProt and uses dynamic `tblastn` translation alignments to extract perfectly in-frame clinical target proteins.
2. **Deterministic Biochemical Scoring (ControlScan):** Mutations are mathematically penalized using a composite score that normalizes **Grantham physicochemical distances** against **BLOSUM62 evolutionary probabilities**, avoiding reliance on noisy/missing transcriptomic data.
3. **Composite Epistasis:** The pipeline maps co-occurring mutation networks using Jaccard-style frequencies weighted by the average ControlScan biochemical severity of the pair, identifying the most dangerous evolutionary pathways.
4. **HTVS Soft-Core Dynamics:** To solve the "steric clash" problem of in-silico mutagenesis without requiring months of supercomputer Molecular Dynamics, Phase 4 utilizes **OpenMM L-BFGS Energy Minimization**. This applies a *Rigid-Backbone, Flexible-Sidechain Approximation* to relax the mutated pocket into a local energy minimum before using **SMINA** to dock the ligand. Results are rigorously reported as $\Delta\Delta G$ (Mutant vs. Wild-type) to isolate the exact thermodynamic cost of the mutation.

---

## ??? Pipeline Architecture (Snakemake DAG)

The workflow is completely modular. Users can start from raw CSV metadata (triggering autonomous downloads) or provide their own local `.fna` assemblies.

* **Phase 1a: Acquisition (Optional)** - Resolves IDs, applies geographic/temporal filters, and idempotently downloads NCBI/BV-BRC assemblies.
* **Phase 1b: Extraction & Calling** - Autonomously fetches UniProt wild-types, extracts sequences via `tblastn`, and calls variants.
* **Phase 2 & 3: Biochemical Epistasis** - Applies ControlScan mathematics and generates publication-ready NetworkX co-occurrence plots.
* **Phase 4: HTVS Biophysics (Optional)** - If a base PDB is provided, dynamically mutates the top 5 epistatic networks, minimizes energy, and runs thermodynamic docking.

---

## ?? Quick Start Guide (Dockerized)

Because MutationScan relies on heavy cross-platform bioinformatics binaries (`tblastn`, `openmm`, `smina`), the entire engine runs inside an isolated Docker container. **It will run identically on Windows, Mac, or Linux.**

### 1. Prerequisites & Setup
Ensure you have [Docker](https://docs.docker.com/get-docker/) and `docker-compose` installed.

```bash
git clone [https://github.com/vihaankulkarni29/MutationScan.git](https://github.com/vihaankulkarni29/MutationScan.git)
cd MutationScan

# Build the bioinformatics container
docker-compose build
```
