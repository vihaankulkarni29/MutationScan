# Benchmarks

Measured performance for MutationScan, with the commands to reproduce each
number. Every figure here was taken on the machine described below; your own
numbers will differ, and the *shape* of the results matters more than the
absolute values.

## Test machine

| | |
|---|---|
| CPU | 12 logical cores |
| OS | Windows (workload uses the `spawn` start method) |
| Python | 3.13 |
| tblastn | 2.17.0 |
| Cohort | 570 genomes (`Beta Barrel Project`), 15 genes |

**Read run-to-run variance before trusting a single figure.** This machine
swung by up to 2x on identical work — the same full pipeline run reported
148.9s and 165.2s, and one repeated micro-benchmark measured `fisher_exact` at
787µs on one pass and 2581µs on another. Where a number moved, it is reported as
best-of-N and labelled as such. Prefer running `benchmarks/parallel_bench.py`
yourself before drawing conclusions.

## Full pipeline

570 genomes, `--threads 12 --variant-threads 12`, wall clock per stage:

| stage | single reference panel | with `species:` routing |
|---|---|---|
| preflight | 2.7s | 0.9s |
| extraction | 87.8s | 68.5s |
| variants | 49.1s | 11.3s |
| co-occurrence | 0.1s | 0.1s |
| epistasis | 9.3s | 3.4s |
| **total** | **~149–165s** | **84.1s** |

Species routing roughly halves extraction and variant calling because each genome
is then aligned only against the genes its own species has, rather than the whole
panel. See `config/config.yaml` for the `species:` block.

Extraction is now the dominant cost at ~59% of wall clock. That floor is
`tblastn` itself: it translates all six reading frames of the subject on every
invocation, once per genome, regardless of how many queries it is given.
Batching to one multi-query call per genome already removed the per-gene
multiplier (570 calls instead of 8,550).

## Executor scaling

`python benchmarks/parallel_bench.py --verify-thresholds` — CPU-bound process
workload, 4 workers, serial as the baseline:

| n | serial | parallel | speedup |
|---|---|---|---|
| 8 | 0.20s | 0.31s | 0.66x |
| 32 | 2.18s | 0.54s | 4.04x |
| 64 | 4.49s | 0.65s | 6.93x |
| 128 | 8.87s | 0.99s | 8.98x |
| 256 | 12.53s | 2.05s | 6.10x |
| 512 | 33.77s | 5.31s | 6.36x |

All sizes produced output identical to the serial run.

Measured pool start-up: **382ms** for 4 workers, ~95ms per worker. The n=8
regression is the expected shape — below roughly the start-up cost, a pool costs
more than it saves. Chunk sizing keeps this from stranding workers at scale, and
the threshold check fails if that regresses.

## Extraction scaling

60-genome subset, 900 alignments (15 genes × 60 genomes), best of 2:

| threads | time | speedup |
|---|---|---|
| 1 | 96.0s | 1.00x |
| 4 | 26.0s | 3.69x |
| 8 | 15.6s | 6.16x |

Near-linear, as expected: each `tblastn` call is ~1.6s of external work and pays
process start-up once per *genome*, not once per alignment.

## Epistasis stage

Before the restructure (146,070 variant pairs, one Fisher test each, in a process
pool): **115.1s** at 12 cores, with the pool spending most of its life paying
`spawn` import costs rather than doing work.

After, best of 3, same data:

| component | time |
|---|---|
| `mutation_severity` | 0.44s |
| `genome_mutation_map` | 3.43s |
| `detect_epistasis` | 5.84s |
| **stage total** | **9.71s** |

12x on the stage. The win is algorithmic, not parallel: pair denominators depend
only on the two genes involved, which yields ~27 distinct universes rather than
146,070, and Fisher's exact test is a pure function of the integer 2×2 table, so
tables are memoised and tested once each — 146,070 tests become ~1,700. Output
was verified byte-identical.

## Mutation counts

Not a speed metric, but the largest behavioural change measured, and the one
most likely to invalidate a comparison against older output.

| configuration | mutations |
|---|---|
| one reference panel across a 3-species cohort | 66,416 |
| `species:` routing | **4,157** |

Per-species mean alignment identity after routing: 97.03% / 99.74% / 99.66%.
The single-panel figure counts interspecies sequence differences as acquired
mutations; see the `Species` section of the README for the mechanism.

## Resource scaling

`python -m benchmarks.scaling --references data/references` measures
wall-clock and peak RSS across cohort sizes and thread counts, writing
`scaling.csv`, `scaling.json` and a two-panel figure. Synthetic cohorts
are regenerated at each size, so the measurement reflects a real run end
to end -- generation is excluded from both time and memory, which is the
honest cost of the pipeline itself.

**Panel A -- cohort size vs. wall-clock, at 2 threads:**

| genomes | time | genomes/hour |
|---|---|---|
| 10 | 6.8s | 5,294 |
| 20 | 11.9s | 6,050 |
| 40 | 17.7s | 13,559 |
| 80 | 33.3s | 8,649 |

Linear `O(N)`: doubling the cohort roughly doubles the time, as expected
when the dominant cost is one `tblastn` call per genome. The 10→20 step
is the noisy end -- the smallest cohort pays fixed start-up twice for
little work to amortise it.

**Panel B -- CPU cores vs. peak RSS, at 10 genomes:**

| cores | peak RSS |
|---|---|
| 2 | 168 MiB |
| 4 | 200 MiB |
| 8 | 264 MiB |

Flat under 2 GiB. The parallelism is process-level (`tblastn`), so each
worker loads its own copy of the reference panel into memory rather than
sharing it; the curve is the honest footprint, not a shared-memory
optimisation.

Peak RSS is now measured, not reported as 0: `Tool.run` spawns the child
with `Popen` and samples the summed RSS of the process tree on a
background thread at 50 ms, started *before* `communicate()`. The old
`peak_memory_mb` returned 0 on every platform except Linux, where it
reached for `/usr/bin/time -v`; the benchmark now records a real number
wherever `psutil` is installed. See `notes/developer-notes.md` for the
2026-10-10 entry.

## AMR tool comparison

`python -m benchmarks.run_benchmark --dataset synthetic` scores MutationScan
against AMRFinderPlus, Mykrobe, Kleborate and StarAMR on a shared ground
truth. Synthetic spike-in (default) is the only source with exact ground
truth for sensitivity, specificity and dN/dS; PATRIC and NCTC 3000 are
stubs pending manual download. 20 genomes, 5 genes (`data/references`):

| metric | MutationScan |
|---|---|
| sensitivity | 1.000 |
| precision | 1.000 |
| dN/dS Pearson | 0.991 |
| lineage concordance | 1.000 |

The synthetic generator back-translates each reference to a CDS, introduces
a known mix of non-synonymous substitutions (the mutations under test) and
silent changes (so dS is non-zero), and records every change. Because the
genome *is* the CDS, tblastn finds it at 1..N and the pipeline's own dN/dS
stage reads back exactly what was written, so a disagreement is a bug rather
than a coordinate artefact.

**A finding this harness caught.** `tolC` carries a signal peptide
(`MKKLLPILIGLSLSGFSSLS`, a hydrophobic N-terminus) that tblastn's
low-complexity filter masks, so the alignment starts at residue 21 and the
extracted CDS is 1419 nt against a 1479 nt reference. The old
`compute_dn_ds` skipped any genome whose CDS length differed from the
reference, which dropped *every* genome of `tolC` -- the gene was simply
absent from `dn_ds.csv`. The fix trims the reference CDS to the hit's
subject region before counting sites and differences. See
`notes/developer-notes.md` for the 2026-10-10 entry.