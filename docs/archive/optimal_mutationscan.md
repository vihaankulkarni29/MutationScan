# Optimal MutationScan

**Vision:** A high-throughput, laptop-friendly AMR sequence analysis tool that democratizes bioinformatics — fast on ordinary hardware, grounded in hard biochemistry and mathematics, with no structural docking stage.

---

## Design Principles

1. **Sequence-only** — no PDB, no docking, no OpenMM/Smina/Meeko
2. **Throughput-first** — batch and parallelize across genomes; deduplicate identical haplotypes
3. **Biochemistry over heuristics** — every score traceable to matrices, scales, or statistics
4. **Honest epistasis** — co-occurrence plus significance testing, not just raw counts
5. **Single-command, single output folder** — `data/output/{job_name}/`

---

## What We Remove

| Removed | Reason |
|---------|--------|
| Phase 3 biophysics (PDBFixer, Smina, ligand prep) | Heavy deps, slow, laptop-unfriendly |
| OpenMM, pdbfixer, smina, meeko, rdkit (docking) | Not needed for sequence pipeline |
| Arbitrary top-5 network cap | Loses signal; rank all, filter by significance |
| Runtime binary downloads | Non-reproducible, fragile |

---

## Optimal Pipeline (4 stages)

```
Genomes (.fna)
    │
    ▼
┌─────────────────────────────────────────────────────────┐
│  STAGE 1 — INGEST & EXTRACT          [parallel tblastn] │
│  QC genomes → extract target proteins → dedupe haplotypes │
└─────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────┐
│  STAGE 2 — VARIANT CALLING           [BLOSUM62 align]   │
│  Query vs WT → substitutions only → identity gate       │
└─────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────┐
│  STAGE 3 — BIOCHEMICAL SEVERITY      [ControlScan++]    │
│  Per-mutation: evolutionary + chemical + property delta │
└─────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────┐
│  STAGE 4 — EPISTASIS & COHORT STATS  [math + networks]  │
│  Frequencies → enrichment → composite epistatic scores    │
└─────────────────────────────────────────────────────────┘
    │
    ▼
Outputs: genomics report, severity table, epistasis networks, cohort summary
```

---

## Stage Details

### Stage 1 — Ingest & Extract (HT)

- **Input:** `data/local_genomes/*.fna`, `config/acr_targets.txt`, seeded `refs/*_WT.faa`
- **Parallel tblastn** across genomes (multiprocessing pool)
- **QC:** genome size, N-count, extraction success rate per gene
- **Haplotype deduplication:** identical protein sequences across accessions counted once for alignment; frequency tracked in metadata
- **Output:** `proteins/`, `refs/`, `extraction_qc.csv`

### Stage 2 — Variant Calling

- Global alignment (BLOSUM62, gap-aware residue counter) — keep current algorithm
- Configurable `variant_min_identity_percent` (default 80)
- Substitutions only; indels flagged separately (optional column)
- **Output:** `1_genomics_report.csv` — Accession, Gene, Mutation, Frequency_context, …

### Stage 3 — Biochemical Severity (ControlScan++)

Per mutation, compute a **transparent severity vector** (not a black box):

| Metric | Source | Meaning |
|--------|--------|---------|
| BLOSUM62 score | Evolutionary substitution matrix | How rare is this AA change? |
| Grantham distance | Composition, polarity, volume | Chemical dissimilarity |
| Property deltas | Kyte-Doolittle, charge at pH 7, volume | Hydrophobicity / charge / size shift |
| Disorder delta | Optional IUPred-style propensity | Local flexibility change proxy |
| **Composite Severity (0–100)** | Weighted combination | Single triage score |

**Conformation proxy (no structure):** large Grantham + charge reversal + proline/glycine introduction → high "fold-disrupting potential" flag. This replaces docking-based "binding change" with sequence-derived biophysics.

- **Output:** `2_mutation_severity.csv` — merged into genomics report or sidecar

### Stage 4 — Epistasis & Cohort Statistics

For each mutation pair (A, B) co-occurring in ≥1 genome:

| Stat | Formula / method |
|------|------------------|
| Observed count | Genomes with both A and B |
| Expected count | `P(A) × P(B) × N` under independence |
| Enrichment ratio | observed / expected |
| Fisher's exact p-value | Significance of co-occurrence |
| Pair severity | `(severity_A + severity_B) / 2` |
| **Epistatic score** | `enrichment × pair_severity × −log10(p)` |

- Build full ranked network table (no arbitrary top-5 cut)
- Optional: NetworkX export for top significant pairs only
- **Output:** `3_epistasis_networks.csv`, `ControlScan_Networks/` (plots for significant pairs)

### Final Outputs

```
data/output/{job_name}/
├── 1_genomics_report.csv       # mutations per genome
├── 2_mutation_severity.csv     # biochemical scores per mutation
├── 3_epistasis_networks.csv    # ranked co-occurrence + statistics
├── cohort_summary.json         # N genomes, genes, top mutations
├── extraction_qc.csv             # per-genome extraction stats
└── ControlScan_Networks/         # optional plots (significant pairs only)
```

---

## Dependencies (Laptop-Friendly)

**Required:** Python 3.10+, pandas, biopython, networkx, matplotlib, pyyaml, snakemake  
**External binary:** `tblastn` (BLAST+) only  
**Removed:** openmm, pdbfixer, smina, meeko, rdkit (docking use)

---

## Snakemake Rules (Optimal)

```
rule all:
    → 1_genomics_report.csv
    → 2_mutation_severity.csv
    → 3_epistasis_networks.csv
    → cohort_summary.json

rule extract_proteins:     # Stage 1
rule call_variants:        # Stage 2
rule biochemical_scoring:  # Stage 3 (split from old combined script)
rule epistasis_networks:   # Stage 4 (statistics + plots)
```

No `htvs_biophysics` rule.

---

## Performance Targets

| Cohort size | Target (laptop, 4 cores) |
|-------------|--------------------------|
| 100 genomes × 5 genes | < 10 min |
| 1,000 genomes × 5 genes | < 2 hr (with parallel tblastn) |
| 10,000 genomes | Batch mode; haplotype dedup critical |

---

## Migration from Current MutationScan

| Current | Optimal |
|---------|---------|
| `03_biochemical_epistasis.py` (combined) | Split into `03_biochemical_scoring.py` + `04_epistasis_networks.py` |
| `04_htvs_biophysics.py` | **Delete** |
| Snakefile Phase 3 rule | **Remove** |
| `environment.yml` docking deps | **Remove** openmm, pdbfixer, smina, meeko |
| ControlScan `MutationScorer` | **Extend** with property deltas + fold-disrupt flags |
| Epistasis `freq × severity` only | **Add** Fisher's exact, enrichment, −log10(p) |
