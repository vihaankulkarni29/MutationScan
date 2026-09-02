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

Three stages, in order.

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
entirely rather than reported as zero.

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
| `run_summary.json` | Versions, settings, per-stage status counts, and the numerator and denominator behind every count. Attach this to a manuscript. |
| `refs/`, `proteins/` | The exact reference used for every call, and every protein extracted. |

**Read `N_Eligible` before comparing two pairs.** Different pairs can have
different denominators — that is the honest behaviour, not a bug, and it is why
the column exists.

## What this tool does not do

Co-occurrence here is **a descriptive count and nothing more**. When one protein
is mutated, was another protein also mutated in that same genome? That is the
whole question, and the answer is a number.

There is no p-value, no enrichment test, no odds ratio, no severity score, no
network ranking and no causal claim anywhere in the output. Two genes appearing
together often may reflect shared ancestry, sampling bias in how the cohort was
assembled, population structure, or nothing at all. Distinguishing those is your
work, not the tool's.

It also does not predict phenotype, map mutations to drugs, model structure or
estimate binding. It reports what the alignment says.

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

## License

See [LICENSE](LICENSE).
