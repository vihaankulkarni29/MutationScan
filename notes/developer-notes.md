# Developer notes

Dated record of substantive changes: what moved, why, and what is still open.
Newest first. Measurements are in [`benchmarks.md`](benchmarks.md).

For historical v2 → v3 rebuild documents see
[`docs/archive/`](../docs/archive/README.md); that code was removed in `d5f0bb0`.

---

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

1. **Loss-of-function alleles in the epistasis stage.** The new
   truncation and frameshift events live in variant calling and QC only.
   Open questions before wiring them in: should `acrR:E67*` be tested
   for epistasis against `acrB:V139F`? Should frameshifts be compared
   only to other frameshifts, or to every variant in the same gene?
   And should a gene with a truncation in one genome and a frameshift
   in another collapse into one loss-of-function allele class, so the
   three frameshift carriers in a cohort test as one functional class
   rather than three distinct variants? Note that the biochemical
   severity scorer never sees these events today: it runs on
   `mutations.csv`, which carries residue-to-residue changes only, so
   BLOSUM62/charge/hydropathy scoring on a stop or frameshift allele
   is undefined by construction and needs a separate severity model
   (a knockout is already flagged in QC, which is the honest answer
   until that model exists).
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