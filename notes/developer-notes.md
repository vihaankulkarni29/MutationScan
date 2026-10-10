# Developer notes

Dated record of substantive changes: what moved, why, and what is still open.
Newest first. Measurements are in [`benchmarks.md`](benchmarks.md).

For historical v2 → v3 rebuild documents see
[`docs/archive/`](../docs/archive/README.md); that code was removed in `d5f0bb0`.

---

## 2026-10-10 — AMR benchmarking harness & signal-peptide dN/dS fix

**Problem.** MutationScan had no way to answer "how does it compare to the
tools a reader would reach for?" -- AMRFinderPlus, Mykrobe, Kleborate,
StarAMR. Each speaks its own dialect, and a sensitivity number means
something different in each. Separately, the dN/dS stage silently dropped
every gene whose tblastn hit did not span the whole reference.

**Fix (harness).** A new `benchmarks/` package, run as
`python -m benchmarks.run_benchmark --dataset synthetic`:

- `schemas.py` -- canonical mutation/epistasis/dN/dS/lineage records, so a
  match is a string comparison and every tool is scored by the same code.
- `tools/` -- one wrapper per external tool, each normalising its native
  output (AMRFinderPlus TSV, Mykrobe JSON, Kleborate per-genome TSV,
  StarAMR ResFinder/PointFinder TSVs). A wrapper's `available` check skips
  a tool whose binary is absent, so a missing tool is never scored as
  "found nothing".
- `amr_benchmark.py` -- metrics: sensitivity/precision/specificity/novel
  recall at the mutation level, precision/recall/F1 at the epistasis level
  (matched on *gene* pairs, so a tool reporting variant pairs is comparable
  to gene-level ground truth), Pearson/Spearman/MAE at the dN/dS level,
  lineage concordance, and performance (wall time, peak RSS, throughput).
- `download_data.py` -- a synthetic spike-in generator with exact ground
  truth: it back-translates each reference to a CDS, introduces a known
  mix of non-synonymous substitutions (the mutations under test) and
  silent changes (so dS is non-zero), and records every change. Because the
  genome *is* the CDS, tblastn finds it at 1..N and the pipeline's own
  dN/dS stage reads back exactly what was written -- so a disagreement is a
  bug, not a coordinate artefact.
- `report.py` -- Markdown + HTML comparison tables.
- `.github/workflows/amr-benchmark.yml` -- weekly + push CI.

**Fix (dN/dS).** `compute_dn_ds` used to skip any genome whose extracted
CDS length differed from the reference, which dropped *every* genome of a
gene whose tblastn hit was partial. The benchmark surfaced this: `tolC`
carries a signal peptide (`MKKLLPILIGLSLSGFSSLS`, hydrophobic N-terminus)
that tblastn's low-complexity filter masks, so the alignment starts at
residue 21 and the extracted CDS is 1419 nt against a 1479 nt reference.
The fix trims the reference CDS to the hit's subject region (using
`Nucleotide_Start`) before counting sites and differences, and accumulates
site counts per genome so different hit regions cannot collide. Genes are
no longer silently absent from `dn_ds.csv`.

**Verified by.** Self-benchmark on a 20-genome synthetic cohort: sensitivity
1.0, precision 1.0, dN/dS Pearson 0.991 (MAE 0.031), lineage concordance
1.0. `tests/unit/test_evolution.py` gained three partial-hit cases. All 378
tests (359 unit + 19 integration) pass; `ruff` + `mypy` clean.

**Known limitation.** The synthetic ground-truth dN/dS is computed over
the full CDS while the pipeline computes over the tblastn-aligned region,
so for a gene with a masked N-terminus the absolute values differ slightly
(the correlation is what the metric reports, and it is 0.991). The
epistasis precision is conservative: the ground truth records only the
*designed* co-mutating pairs, while the pipeline also reports genuinely
significant random co-associations in a small cohort.

## 2026-10-09 — Lineage tagging & dN/dS selection pressure

**Problem.** The reviewer warned that population structure can confound epistasis:
if a pandemic clone (e.g. ST131) dominates the cohort, shared mutations by
descent inflate Fisher-exact p-values. A companion requirement was evidence
of positive selection on mutated genes (dN/dS > 1).

**Fix (lineage).** A new `lineage_manifest` config option (TSV/CSV mapping
Accession → Lineage/ST) wires into the pipeline:
- Genome mutation map gains a `Lineage` column.
- Epistasis networks gain `Lineage_Distribution` (e.g. `ST131:3,ST11:2`)
  and `Clonality` columns (`Clonal` = single lineage, `Convergent` =
  ≥2 lineages). Convergent pairs are the gold standard for true
  biological epistasis.

**Fix (dN/dS).** A new `evolution` module computes per-gene dN/dS:
- Extracts actual CDS from genome FASTA using the tblastn coordinates now
  recorded in the extraction manifest (`Nucleotide_Start`, `Nucleotide_End`).
- Back-translates the reference protein with the standard genetic code to
  obtain a reference CDS.
- Compares codon-by-codon to count synonymous (dS) and non-synonymous (dN)
  nucleotide differences.
- Outputs `dn_ds.csv` with per-gene `N_Sites`, `S_Sites`, `dN`, `dS`,
  `dN_dS`, `N_Substitutions`, `N_Genomes_Evaluated`. Genes under positive
  selection (dN/dS > 1) appear first.

**Verified by.** `tests/unit/test_lineage.py` (14 cases: manifest parsing,
distribution counting, clonal/convergent classification, genome map and
epistasis annotation), `tests/unit/test_evolution.py` (14 cases: genetic
code completeness, codon site counting, back-translation, empty/edge cases).
All 375 tests pass (356 unit + 19 integration). `ruff` + `mypy` clean.

## 2026-10-08 — Summary provenance for loss-of-function-only genomes

**Problem.** `build_run_summary` recorded `genomes_with_any_mutation` from
`mutations.csv` alone. A genome whose only called allele is a LoF knockout or
frameshift was counted in the co-occurrence genome map and the epistasis
network (both LoF-aware from `a447db3`) but dropped from the run summary's
mutated-genome count, so the provenance record under-reported mutated genomes
and offered no LoF fields at all.

**Fix.** `summary.build_run_summary` now merges LoF alleles through the same
`variants.with_loss_of_function(mutations, qc)` entry point the downstream
stages already use, so `genomes_with_any_mutation` sees LoF-only genomes too.
The merge is a no-op on a `None` mutations stream (a dry/partial run) and on a
QC lacking the event columns (an older run), so pre-existing summaries are
unchanged. Three fields were touched: `loss_of_function_alleles`,
`genomes_with_loss_of_function`, and the now-LoF-augmented
`genomes_with_any_mutation`. `mutations` itself stays substitution-only
(`int(len(mutations))`) so the two counts remain independently auditable.

**Verified by.** `tests/unit/test_summary.py` (5 cases: LoF-only genome counted,
older QC without event columns is a no-op, empty substitutions + LoF carrier,
`None` mutations reports zero, pure substitution run unaffected). `mypy` and
`ruff` clean. All 328 unit + 19 integration tests pass.

## 2026-10-07 — Frameshift and premature-stop calling

**Problem.** The variant caller reported only substitutions, deletions and
insertions. A premature stop codon in the extracted protein was reported as a
substitution to `*` (e.g. `E67*`), and an out-of-frame indel was invisible:
it either aligned as a long stretch of mismatches or, when tblastn's best HSP
ended at the disruption, as a plain coverage drop. Neither was distinguishable
from ordinary divergence.

**Fix.** `compare_sequences` walks the alignment once and calls two events
alongside substitutions:

- A query stop codon aligned against a reference residue is a **truncation**,
  recorded as `(position, ref_aa, fraction_lost)` -- the fraction lost is the
  share of the reference protein downstream of the stop -- and rendered
  `Q30*`. A stop that matches the reference's own `*` is a match, not an
  event.
- An **internal** indel run whose length is not a multiple of three is a
  **frameshift**, recorded as `(position, signed_length)` and rendered
  `fs181 (-4)`. A run touching either end of the alignment is a coverage
  overhang, not an indel, and is never reported: when a frame shift occurs,
  tblastn's best HSP is the clean shifted-frame fragment, so the disruption
  lands at an alignment end and the honest signal is the coverage drop.

A truncation removing more than a fifth of the reference is a **functional
knockout** (`knockout_lost_percent`, default 20, settable in `config.yaml`
or `--knockout-lost-percent`; the Level 4 biophysical disruption tier).
Stops are no longer substitutions: `mutations.csv` carries
only residue-to-residue changes, and the events surface in `variant_qc.csv`
as three new columns (`N_Frameshifts`, `N_Truncations`,
`N_Functional_Knockouts`) with the event text in `Detail` for called pairs,
e.g. `truncation E67* (68.8% lost) [functional knockout]`.

**Notation.** Event positions are reference coordinates, matching the
mutation table. The signed length in a frameshift event is measured in
**amino-acid residues** -- the unit of the alignment -- not nucleotides;
`fs181 (-4)` means a four-residue deletion starting at reference residue
181. Indel runs are judged independently, so a gene carrying several
out-of-frame runs reports each of them separately (and the QC counts
them separately), which is the compound case. A run touching either end
of the alignment is a coverage overhang, not an indel, and is never
reported.

**Verified by.** 303 unit tests (up from 285) and 19 tblastn
integration tests. An end-to-end run with a planted TAA stop codon reports
`truncation Q30* (76.0% lost) [functional knockout]` at 99.2% identity; a
planted 1-bp insertion shows the coverage signal (84%) with no spurious
event. On the 570-genome baseline: all 66,414 non-stop mutation rows and all
nine pre-existing QC columns are byte-identical, and the only mutation rows
lost are the two `acrR E67*` stop substitutions (a GCA/GCF mirror pair), now
truncation events. The cohort carries 9 truncation pairs (8 knockouts) and
5,152 frameshift pairs, 555 of them called.

## 2026-10-07 — Loss-of-function alleles in the epistasis stage (LoEPI)

**Problem.** The frameshift, truncation and functional-knockout events live in
`variant_qc.csv` only. Co-occurrence, the genome map and epistasis still saw
substitutions alone, so a knockout such as `acrR E67*` was invisible to every
stage after variant calling: `acrR:LoF` never entered a Fisher table against
`acrB:V139F`, and the severity scorer on `mutations.csv` could never see it.

**Change.** Loss-of-function is a collapsed allele class that rides in the same
stream as substitutions, so the downstream stages answer "is this gene mutated
in this genome?" with one consistent denominator.

- `variants.loss_of_function_alleles(qc)` emits one `LoF` row per called
  `(genome, gene)` carrying a functional knockout or a frameshift; the original
  event text rides along in `LoF_Events` so the collapsed class can be
  stratified by event type later. Non-knockout truncations (the protein is
  mostly intact) stay in QC only and do **not** cross the LoF bar.
- `variants.with_loss_of_function(mutations, qc)` merges that stream in and is
  idempotent: a frame that already carries `LoF` rows is returned unchanged, so
  a stage cannot double-count its own output.
- `cooccurrence`, `genome_mutation_map` and `detect_epistasis` each merge once on
  entry; `epistasis_networks` merges for `mutation_severity` and lets the other
  two no-op on the already-merged frame. `GENOME_MAP_COLUMNS` gained `N_LoF`.
- `scoring`: `compute_severity`/`mutation_severity` now guard on
  `Mutation == "LoF"` before parsing the residue string. Without the guard the
  literal "LoF" parses as L→F and emits a bogus BLOSUM62 score; a
  loss-of-function allele is not a substitution, so it is scored categorically
  at `LOF_SEVERITY = 100.0` (maximal functional impact) with empty
  `Ref_AA`/`Alt_AA` and zeroed biochemical deltas. `pair_severity` then averages
  that 100 with the partner's score.

**Design decisions (agreed with review).** Frameshifts and nonsense/truncations
collapse into a single `LoF` class for statistical power; the event-type records
are retained in `LoF_Events` for the stratification escape hatch. A gene with a
truncation in one genome and a frameshift in another tests as one functional
class.

**Caveat — frameshift composition.** The detector runs on the *protein*
alignment, so a flagged "frameshift" is a residue gap run whose length is not a
multiple of three. On this cohort that turns out to be **in-frame nucleotide
indels** (a residue gap of L is 3L nucleotides, always in-frame): 732 flagged
events are all lengths {1, 2, 4, 8, 64} residues (483x1, 194x4, 27x2, 24x64),
none a multiple of three. True nucleotide frameshifts (1–2 nt) scramble the
translation and surface instead as alignment breaks → coverage drops,
which the detector treats as coverage overhangs and never tags as frameshifts.
So the `LoF` class here leans on (mostly short) in-frame deletions plus rare
true premature-stop knockouts; `LoF_Events` preserves the distinction so the
collapsed signal can be re-split if it proves noisy.

**Verified by.** `tests/unit/test_variants.py` (LoF extraction, knockout bar,
rejected-pair exile, idempotence, end-to-end through `call_variants`),
`tests/unit/test_epistasis.py` (LoF in the genome map and pair tables, severity
guard, idempotent merge), `tests/unit/test_cooccurrence.py` (LoF carriers count
as mutated, rejected LoF excluded). 353 unit + 19 integration tests green; ruff
and mypy clean.

**Regression on the 570-genome baseline** (re-run `call_variants` on the cohort
proteins + references; no `tblastn` required for this stage):

- Substitution stream byte-identical (66,414 rows) to the committed variant
  layer; the only rows that layer moves to events are the two `acrR E67*` stop
  substitutions.
- 557 called pairs carry a `LoF` allele (555 frameshifts + 2 functional
  knockouts; no pair carries both). 484 genomes carry ≥1; 557 severity rows
  emitted at `severity_score = 100.0`.
- 960 `gene:LoF` pairs enter the epistasis network; every one co-occurs
  (`N_Both > 0`).
- **Non-regression:** the 49,052 pre-existing (non-LoF) variant pairs keep
  byte-identical `N_Both`/`N_A`/`N_B`/`N_Eligible`/`Odds_Ratio`/`P_value`;
  gene-level co-occurrence is byte-identical (0 genes gain a carrier set);
  non-carrier genome-map rows are identical. The only drift on common pairs is
  `Q_value` (max 0.0019, the FDR family grows by 960 tests) with **0**
  significance-class (`Sign`) flips.

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

1. **Loss-of-function alleles in the epistasis stage.** — resolved, see the "2026-10-07 — Loss-of-function alleles in the epistasis stage" entry above. The three open questions were answered in design review: `acrR:E67*` tests against the same gene-pair space as `acrB:V139F`; frameshifts + nonsense collapse into one `LoF` class with the event type retained in `LoF_Events` for stratification; a truncation in one genome and a frameshift in another are one functional class.
2. **Promoter / intergenic frameshifts.** The frameshift detector works
   on protein alignments, so it captures coding-sequence events only.
   Regulatory frameshifts (e.g. in a promoter region of `ramA`) or
   intergenic indels affecting transcriptional attenuation are
   invisible. Absent as a stage: blocked on reference
   promoter/operator sequences, which this dataset does not contain --
   it needs a new input concept (nucleotide-aware alignments) rather
   than a flag.
3. **Within-*E. coli* identity.** 67% of the remaining 4,157 mutations come from
   pairs below 95% identity — genuine phylogroup diversity against K-12, but
   `tolC` averaged 80.12% identity, meaning the identity gate is doing real work
   it should not have to. Per-phylogroup references would resolve this.
4. **Incomplete per-species reference coverage.** 7 *E. coli* genes, 2
   *Klebsiella*, 6 *Pseudomonas*. Routing reports the gap honestly instead of
   corrupting it, but the *Klebsiella* arm is currently 2 genes rather than 9.
   Acquiring the missing references is data work, not code.
5. **Reference proteome for orthology checking.** The reference directory holds
   only the 15 target proteins. A meaningful reciprocal-best-hit check would need
   a full reference-strain proteome to supply competing paralogs.
6. **Extraction below ~88s.** Now ~59% of wall clock and `tblastn`-bound. Would
   need contig splitting or letting `tblastn` use its own threads; gains are
   modest and there is a risk of altering results silently.