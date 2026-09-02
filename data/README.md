# `data/`

Everything in here is **yours**. Nothing in this folder is part of the tool, and
git tracks only the empty scaffold — the `.gitkeep` files and this README.

You supply two things; MutationScan produces the third.

## What you put in

### `references/` — the reference proteins, which *are* the target list

One protein FASTA per gene you want to study, named after the gene:

```
references/geneA_WT.faa
references/geneB_WT.faa
```

There is no separate list of targets anywhere. **The filenames are the target
list.** Add `geneC_WT.faa` and `geneC` is analysed on the next run; remove a file
and it is not. `{gene}.faa`, `.fasta`, `.fa` and `.fas` work too, and the `_WT`
suffix is optional.

Each file should hold one protein sequence. Headerless files and odd line
wrapping are read fine; two records in one file are rejected with an error that
says so.

If a reference is missing you can have it fetched instead: set `uniprot_taxid` in
the config to your organism's NCBI taxonomy ID, and any absent target is pulled
from reviewed UniProt entries for that organism. Leave it `""` to disable all
network access.

### `genomes/` — the cohort

One nucleotide assembly per genome, as `.fna`, `.fasta`, `.fa` or `.fna.gz`:

```
genomes/SAMPLE_001.fna
genomes/SAMPLE_002.fna
```

**The filename becomes the genome's accession** in every output table, so name
the files however you want to see them in your results. Dots are fine.

Point the tool at a folder anywhere on disk with `--genomes`; this one is just
the default.

## What comes out

### `output/{job_name}/`

One folder per run, holding the manifest, the mutation and QC tables, the
co-occurrence counts and `run_summary.json`. See the main
[README](../README.md) for what each file contains.

Re-running the same `job_name` overwrites that folder. Use a different
`--job-name` to keep runs side by side.

## The rest

`databases/` and `results/` are scaffold folders kept for convenience. Any other
folder you create under `data/` is ignored by git automatically, so private
cohorts and scratch directories stay out of version control without you having to
edit `.gitignore`.
