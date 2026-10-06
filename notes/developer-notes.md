# Developer notes

Dated record of substantive changes: what moved, why, and what is still open.
Newest first. Measurements are in [`benchmarks.md`](benchmarks.md).

For historical v2 → v3 rebuild documents see
[`docs/archive/`](../docs/archive/README.md); that code was removed in `d5f0bb0`.

---

## 2026-10-05 — Species-aware reference routing (`e3d8bc0`)

**Problem.** The cohort spans three species (353 *Klebsiella pneumoniae*, 163
*E. coli*, 54 *Pseudomonas aeruginosa*) but a single reference panel was applied
to all 570 genomes. E. coli K-12 reference proteins align to *Klebsiella*
genomes at ~78% mean identity — near the 80% gate — so ordinary interspecies
differences were being reported as acquired mutations.

**Effect.** 66,416 mutations → 4,157. Per-species mean identity is now 97.03% /
99.74% / 99.66%. Mutations from alignment pairs below 95% identity fell from 97.7%
to 67.0%, and the remainder is within-*E. coli* diversity rather than
interspecies noise.

**Note for review.** Reciprocal best hit does not address this and was
considered before the routing approach. *E. coli* acrB and *K. pneumoniae* acrB
are true orthologs, so RBH confirms them and passes them through — it separates
paralog from ortholog, whereas the issue is an ortholog from the wrong species. An
RBH pass against the available reference directory would also have been inert,
since it contains only the 15 target proteins and therefore no competing
paralogs to disambiguate against.

**Design.** New `src/mutation_scan/species.py` owns panel parsing and
genome→species assignment (`fnmatch` on the filename stem, with a literal-length
specificity tie-break). `references.py` stages per-species references under
`refs/{species}/`. `extract.py` builds one query panel per species, so genomes of
a species still share a single `tblastn` call. `variants.py` reads the genome's
species from the extraction manifest rather than re-deriving it, so the two
stages cannot disagree about which reference applied. A gene a species has no
reference for is reported `no_reference_for_species` and is never compared against
another species' protein; genomes matching no pattern abort the run.

**Also fixed.** `discover_targets` treated multi-record `combined_*.faa` bundles
as genes and aborted the run on an unreadable "reference" it never needed.

**Compatibility.** Opt-in. With no `species:` block the flat path is unchanged
and verified byte-identical to the pre-feature baseline.

## 2026-10-05 — Epistasis restructure (`3e33d49`)

**Problem.** The stage ran one Fisher exact test per variant pair in a process
pool. Pooling was a net loss: it gained ~2x at four workers and ran *slower* at
eight, because the pool spent most of its life paying `spawn` import costs
(~3.4s per worker) against only ~50s of real test work. Measured serial I/O
floor was 0.54s, which rules out an Amdahl explanation.

**Change.** Two exact rearrangements replace the per-pair loop:

- A pair's denominator depends only on its two genes, so pairs group into ~27
  universes instead of 146,070; each universe's co-occurrence counts come from one
  Gram matrix (`G.T @ G`) rather than one set intersection per pair. float64 is
  exact at these counts, so results are bit-identical regardless of BLAS thread
  count.
- Fisher's exact is a pure function of the integer 2×2 table, so tables are
  memoised and tested once each: 146,070 tests → ~1,700.

Row order is preserved because the final sort is stable and ties fall back to
enumeration order.

**Also in this commit.** `mutation_severity` memoises `compute_severity` on the
residue pair (~600 distinct values) instead of recomputing per row through
`iterrows`; `genome_mutation_map` drops a row-wise `apply`. Stage total went
115.1s → 9.7s. `run_summary.json` gained `cores_available` so requested-vs-effective
is interpretable from the artifact alone.

**Verification.** Byte-identical against the pre-refactor baseline for
`mutations.csv`, `cooccurrence_genes.csv`, `cooccurrence_genes_matrix.csv`,
`genome_mutation_map.csv`, `mutation_severity.csv` and `epistasis_networks.csv`;
`variant_qc.csv` and `extraction_manifest.csv` identical on all pre-existing
columns.

`epistasis.threads` is retained and recorded as advisory rather than removed, so
existing scripts and config files keep working.

## 2026-09-30 — Parallelisation and extraction batching (`d965bb6`)

Every stage parallelised. `tblastn` batched to one multi-query call per genome
(570 calls instead of 8,550 — it re-translates all six reading frames on every
invocation, so per-gene calls re-translated each genome 15 times over). Per-stage
timings and requested-vs-effective parallelism added to the run summary. Full
570-genome run at this point: 423.7s.

## 2026-09-28 — PyPI release preparation (`afa3d64`, `57a84c8`, `2eb7ac9`)

Project URLs, release workflow, version 1.0.0, PEP 639 licence fix, and release
workflow moved to API-token auth.

## 2026-09-28 — CI hardening (`3da2962`, `e5c0048`)

Snakemake DAG validation and `config-check` added to the integration workflow;
`AGENTS.md` introduced; `config-check` made tolerant of absent input data, which
is the expected state in CI.

## 2026-09-27 — Epistasis integration and dependency fix (`d264950`, `67c369e`, `504b4a6`)

Epistasis module (genome mutation map, Fisher exact, BH FDR, severity scoring)
wired into the pipeline and documented. scipy moved to runtime dependencies.

## 2026-09-01 → 2026-09-02 — v3 rebuild (`91e19df` … `c66d9a4`)

v2 pipeline deleted and rebuilt in stages: config and tolerant FASTA layer, extract
and variants, co-occurrence, orchestration (CLI, run summary, Snakemake workflow),
then removal of every trace of the study the tool was built from.

---

## Open items

Ordered by how much they affect published results.

1. **Frameshift / premature-stop calling.** Internal indels are counted in QC
   (`N_Deleted_Positions`, `N_Inserted_Columns`) but never reported as events. A
   1-bp insertion currently presents as one deleted plus one inserted residue,
   indistinguishable from two substitutions.
2. **Within-*E. coli* identity.** 67% of the remaining 4,157 mutations come from
   pairs below 95% identity — genuine phylogroup diversity against K-12, but
   `tolC` averaged 80.12% identity, meaning the identity gate is doing real work
   it should not have to. Per-phylogroup references would resolve this.
3. **Incomplete per-species reference coverage.** 7 *E. coli* genes, 2
   *Klebsiella*, 6 *Pseudomonas*. Routing reports the gap honestly instead of
   corrupting it, but the *Klebsiella* arm is currently 2 genes rather than 9.
   Acquiring the missing references is data work, not code.
4. **Promoter / non-coding pass.** Absent. Blocked on reference promoter/operator
   sequences, which this dataset does not contain — it needs a new input concept
   rather than a flag.
5. **Reference proteome for orthology checking.** The reference directory holds
   only the 15 target proteins. A meaningful reciprocal-best-hit check would need
   a full reference-strain proteome to supply competing paralogs.
6. **Extraction below ~88s.** Now ~59% of wall clock and `tblastn`-bound. Would
   need contig splitting or letting `tblastn` use its own threads; gains are
   modest and there is a risk of altering results silently.