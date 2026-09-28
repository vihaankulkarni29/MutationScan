# MutationScan

Find amino-acid substitutions in your proteins of interest across a cohort of
bacterial genomes, and count how often pairs of those proteins are mutated in the
same genome.

**No gene, organism, drug or database is hardcoded anywhere in this tool.** You
supply reference proteins and genomes; the reference filenames *are* the target
list. Point it at different proteins and it studies different proteins — there is
nothing to reconfigure and nothing in the code to edit.

---

## What it does

Five stages, in order.

**1. Extract.** For every (genome, gene) pair, `tblastn` aligns the reference
protein against the six-frame translation of the genome and the best hit is
written out as a protein sequence. This finds the gene without needing the genome
to be annotated.

**2. Call substitutions.** Each extracted protein is globally aligned to its
reference (Needleman–Wunsch, BLOSUM62). Walking the alignment, the position
counter advances only on non-gap reference columns — so `Y83F` means *the 83rd
residue of your reference protein*, not the 83rd column of some alignment. Two
gates decide whether a pair is called at all: **identity** to the reference and
**coverage** of it. Coverage is what catches a truncated or fragmented gene,
which identity alone will happily score as near-perfect.

**3. Count co-occurrence.** For each pair of genes, count the genomes where both
are mutated. The denominator is only those genomes where **both genes were
actually evaluated**, so "gene B is not mutated here" is never confused with
"gene B was never observed here". A pair with no genome in common is omitted
entirely rather than reported as zero. Counts are descriptive only — no p-values
or enrichment statistics at this stage.

**4. Epistasis.** For every pair of mutations that co-occur in at least one
genome, Fisher's exact test probes whether they appear together more (or less)
often than chance. Benjamini–Hochberg FDR correction is applied across all
tested pairs, and significant pairs are classified as *positive* (enriched
together), *negative* (mutually exclusive), or *neutral*. Pairs are ranked by
FDR-adjusted p-value, then by composite biochemical severity
(BLOSUM62 + charge + hydropathy + volume), then by co-occurrence count. A
per-genome **mutation map** is also written so you can see every mutation in
every genome at a glance.

---

## Install

```bash
conda env create -f environment.yml
conda activate mutationscan
pip install -e .
```

The only external binary is `tblastn` (BLAST+), which the environment provides.
Check it:

```bash
tblastn -version
```

## Run it

Put one protein FASTA per target gene in a references folder, and your genome
assemblies in another:

```
data/references/geneA_WT.faa      data/genomes/SAMPLE_001.fna
data/references/geneB_WT.faa      data/genomes/SAMPLE_002.fna
data/references/geneC_WT.faa      data/genomes/SAMPLE_003.fna
```

Then:

```bash
mutationscan run --genomes data/genomes --references data/references --out data/output --job-name my_run --threads 4
```

That is the whole interface. Three genes went in because three files were in the
folder. Filenames become gene names (`geneA_WT.faa` → `geneA`; the `_WT` is
optional), and genome filenames become the accessions in every output table.

Check a config without running anything:

```bash
mutationscan config-check --config config/config.yaml
```

### Useful flags

| Flag | Effect |
|---|---|
| `--targets geneA,geneB` | analyse a subset of the references on disk |
| `--min-identity 80` | reject a protein below this % identity to the reference |
| `--min-coverage 80` | reject a protein below this % coverage of the reference |
| `--cooccurrence-level gene,variant` | also count pairs at `gene:mutation` resolution |
| `--min-count 3` | only report pairs co-mutated in ≥ 3 genomes |
| `--taxid 1234` | fetch a missing reference from UniProt for that organism |
| `--threads 8` | concurrent `tblastn` processes |
| `--config FILE` | start from a YAML config; flags override it |

The epistasis stage runs automatically at the end of `mutationscan run`. To
run it separately on an existing mutations table:

```bash
mutationscan epistasis --mutations data/output/my_run/mutations.csv \
                      --qc data/output/my_run/variant_qc.csv \
                      --out data/output/my_run \
                      --fdr-threshold 0.05 --min-count 1 --threads 4
```

## Configuration

Everything above can live in [config/config.yaml](config/config.yaml) instead.
The whole schema:

```yaml
job_name: default_run

genomes_dir: data/genomes           # nucleotide assemblies, one *.fna per genome
references_dir: data/references     # wild-type proteins, {gene}_WT.faa or {gene}.faa
output_root: data/output            # results land in {output_root}/{job_name}/

targets: []                         # [] means "every reference present"
uniprot_taxid: ""                   # "" disables all network access

min_identity_percent: 80
min_coverage_percent: 80

tblastn_binary: tblastn
threads: 4

cooccurrence:
  levels: [gene]                    # add "variant" for {gene}:{mutation} pairs
  min_count: 1

epistasis:
  fdr_threshold: 0.05             # Benjamini-Hochberg q-value cutoff
  min_count: 1                    # require >= this many co-occurring genomes
  threads: 1                       # parallel pairwise Fisher tests
```

Unknown keys are rejected with a message naming the valid ones, so a typo fails
before the run starts rather than silently doing the wrong thing.

## Outputs

All under `data/output/{job_name}/`:

| File | What it holds |
|---|---|
| `extraction_manifest.csv` | One row per attempted (genome, gene) pair, with the tblastn hit stats and a status: `extracted`, `no_hit`, `no_reference` or `error`. A genome where nothing was found is recorded here, not dropped. |
| `mutations.csv` | One row per substitution: genome, gene, `{ref}{position}{alt}`, and the identity/coverage of the alignment it came from. |
| `variant_qc.csv` | One row per pair *considered*, including the ones a gate rejected (`low_identity`, `low_coverage`). This is what distinguishes "not mutated" from "not evaluated", and it supplies the denominators below. |
| `cooccurrence_genes.csv` | One row per gene pair: `N_Both`, `N_A`, `N_B`, `N_Eligible`, and the two plain fractions. |
| `cooccurrence_genes_matrix.csv` | The same counts as a gene × gene matrix; the diagonal is the genomes mutating that gene. |
| `cooccurrence_variants.csv` | The same at `gene:mutation` resolution. Only written when you ask for the `variant` level. |
| `genome_mutation_map.csv` | One row per genome with all its mutations listed (`gene:mut | gene:mut …`). Written by the epistasis stage. |
| `mutation_severity.csv` | Per-mutation biochemical severity (BLOSUM62, charge, hydropathy, volume, composite score). |
| `epistasis_networks.csv` | Ranked epistatic pairs: Fisher's exact p-value, BH-adjusted q-value, odds ratio, sign (positive/negative/neutral), and composite pair severity. |
| `run_summary.json` | Versions, settings, per-stage status counts, and the numerator and denominator behind every count. Attach this to a manuscript. |
| `refs/`, `proteins/` | The exact reference used for every call, and every protein extracted. |

## Epistasis output

`epistasis_networks.csv` is the main statistical output. Each row is an unordered
pair of mutations. Read the columns together:

- `N_Both`, `N_A_only`, `N_B_only`, `N_Neither`, `N_Eligible` — the 2×2 table and
  its denominator (only genomes where **both** mutations were evaluated).
- `Observed` vs `Expected` — how many co-occurrences you saw vs the null model.
- `Odds_Ratio` — enrichment direction. `>1` means they co-occur more than expected,
  `<1` means they exclude each other, `inf` means a structural zero (handled
  gracefully).
- `P_value` — Fisher's exact test on the 2×2 table.
- `Q_value` — Benjamini–Hochberg FDR-adjusted p-value.
- `Sign` — `positive` (enriched, q < threshold), `negative` (excluded, q < threshold),
  or `neutral` (not significant).
- `Pair_Severity` — average biochemical severity of the two substitutions, used
  only for ranking ties.

Pairs are sorted by `Q_value` ascending (most significant first), then severity
descending, then `N_Both` descending.

**Read `N_Eligible` before comparing two pairs.** Different pairs can have
different denominators — that is the honest behaviour, not a bug, and it is why
the column exists.

## What this tool does not do

Co-occurrence here is **a descriptive count and nothing more**. When one protein
is mutated, was another protein also mutated in that same genome? That is the
whole question, and the answer is a number.

It also does not predict phenotype, map mutations to drugs, model structure or
estimate binding. It reports what the alignment says and counts what it sees.

## Snakemake

For large cohorts, the same stages run under Snakemake with checkpointing and
cluster support. It reads the same `config/config.yaml` and calls the same
library functions, so the two entry points cannot produce different results:

```bash
snakemake -n --cores 1
snakemake --cores 8
```

## Docker

```bash
docker compose build
docker compose run --rm mutationscan run --genomes data/genomes --references data/references --out data/output
```

Your `./data` and `./config` are mounted from the host, so results appear on your
machine.

## Tests

```bash
python -m pytest
```

Unit tests need no external tools; the integration tests build a small synthetic
cohort with substitutions planted at known positions and are skipped
automatically if BLAST+ is absent.

## Troubleshooting

| Symptom | Cause |
|---|---|
| `No reference proteins found` | The references folder is empty or holds no `.faa`/`.fasta`/`.fa`/`.fas` files. |
| `tblastn not found` | BLAST+ is not on `PATH`; set `tblastn_binary` to its full path. |
| Everything is `no_hit` | The references are from too distant an organism, or the genome files are not nucleotide assemblies. |
| Many pairs are `low_coverage` | The gene is fragmented across contigs in those assemblies. Lower `--min-coverage` deliberately, or accept the exclusion. |
| A pair is missing from the output | No genome evaluated both of its genes. A fraction over an empty denominator is undefined, so the pair is omitted rather than reported as zero. |

<img width="6107" height="1751" alt="Genome Mutation-2026-09-28-151342" src="https://github.com/user-attachments/assets/91246564-7e77-4342-a45d-e77e305ea506" />

## License

See [LICENSE](LICENSE).
