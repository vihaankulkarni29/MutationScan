# MutationScan — Code & Design Audit

**Date:** 2026-08-30
**Scope:** Full pipeline (genomics extraction → variant calling → epistasis; biophysics phase reviewed but out of DAG).
**Method:** Every finding below was verified against the actual source and the on-disk state (1148 genomes in `data/local_genomes/`, references in `data/results/refs/`, test execution, and `snakemake -n`). This document supersedes `review.md`, which was written against a different version of the code and is substantially inaccurate (see "Corrections to `review.md`").

---

## Bottom line

- `snakemake -n` (dry-run) **passes**. There is no hard-coded `5o66.pdb` in the DAG and no biophysics rule in the `Snakefile` — the previous review's headline is wrong.
- The DAG is three rules only: `extract_proteins → call_variants → biochemical_epistasis` (`Snakefile:29-81`).
- A fresh **real** run fails in Phase 1a for two config reasons (Tier 1).
- If those are fixed, the pipeline runs but the genomics report is **wrong on every row** (Tier 2, bug #3) and carries **no resistance knowledge at all** (bug #4).

---

## Tier 1 — Fresh-run blockers

**1. `config.yaml` points at an empty targets file.**
`config/config.yaml:4` → `targets_file: "config/targets.txt"`, but `config/targets.txt` is comments only. `load_target_genes()` returns `[]`, and `src/scripts/02a_extract_proteins.py:60` exits. The real gene list is in `config/acr_targets.txt`, which nothing references.

**2. References are never seeded on a default run.**
`reference_seed_dir: ""` and `uniprot_taxid: ""` (`config/config.yaml:7,10`). With the seed dir empty, `src/scripts/02a_extract_proteins.py:88` falls back to `Path(".")` (repo root, no `.faa`), seeds 0, then `:96-104` computes `missing_reference_genes` and exits. The WT refs in `data/results/refs/` are only used when `skip_extraction=True` (`Snakefile:22-24`).

---

## Tier 2 — Correctness bugs

**3. Swapped/mislabeled tuple in the no-ML path (headline bug).**
Production uses `enable_ml=False` (`src/scripts/02b_call_variants.py:62`). Every mutation reaches `src/mutation_scan/analysis/variant_caller.py:606-610`:

```python
if not self.enable_ml:
    drug = self.drug_mapping.get(gene_name.lower())
    phenotype = f"{drug} resistance (VUS)" if drug else "Variant of Unknown Significance"
    return (phenotype, "N/A", "N/A", None, "Clinical DB")
```

Contract is `(status, phenotype, pdb, score, source)` (cf. correct hit at `:582-588`). Result in the CSV: **Status** = phenotype text, **Phenotype** = `"N/A"`, **Source** = `"Clinical DB"` (false — nothing was matched). Affects every row.

**4. No resistance database in production.**
`variant_caller.py:96` looks for `resistance_db.json` inside `refs_dir` (the freshly created output refs dir, which only holds `.faa`). There is no `config/resistance_db.json`; the only copy is a test fixture. `_load_resistance_db()` returns `{}` (`:162-164`), so the "Resistant" branch never fires.

**5. VUS counters are always zero.**
`variant_caller.py:310` and `:769` test `df['Status'] == 'VUS'`, but Status is never literally `'VUS'` on any path.

**6. Reference resolution is inconsistent.**
Extraction resolves refs flexibly and case-insensitively (`resolve_reference`, `tblastn_extractor.py:42-64`); variant calling hard-codes `refs_dir / f"{gene_name}_WT.faa"` (`variant_caller.py:333`). A ref seeded as `acrB.faa` extracts fine but yields zero variants (silent `return []`).

**7. Case-sensitive resistance-DB lookup.**
`variant_caller.py:565` uses `gene_name not in self.resistance_db` with no normalization, while adjacent drug-mapping lookups use `.lower()`. Epistasis lowercases gene names (`03_biochemical_epistasis.py:51`) — casing will diverge once a DB is populated.

**8. Integration test encodes a filter failure.**
`tests/integration/test_variant_caller.py::test_multiple_mutations` uses sequences at 75% identity against the 80% default filter → the expected mutations are rejected. Four tests also share one output dir with no isolation.

**9. `_verify_tblastn_available` passes a broken binary.**
`tblastn_extractor.py:172-192`: a non-zero `-version` exit takes no branch and returns `None`; only `FileNotFoundError` is caught.

---

## Tier 3 — Scientific / methodological lacunae

- **"Epistasis" is co-occurrence, not epistasis.** `03_biochemical_epistasis.py:83-89` links mutations sharing a genome, with no background rate, association test, or correction for clonal population structure. On a largely clonal cohort this surfaces lineage markers, not interactions.
- **Severity fallback silently equals `1.0`** (`03_biochemical_epistasis.py:67,69`) — reads as a real low score, not a sentinel.
- **Networks truncated to 5 edges** (`:120`, `.head(5)`), unlogged; plots are 2-node graphs.
- **No frameshift/partial QC beyond percent identity.** Single best HSP (`-max_target_seqs 1 -max_hsps 1`, `tblastn_extractor.py:334-335`), column-13 `sseq` taken verbatim; split genes truncate, internal stops pass through.
- Positional numbering itself is sound (reference-coordinate residue counter after re-alignment).

---

## Tier 4 — Dead code, orphaned phase, packaging, CI, hygiene

- **Biophysics phase is orphaned.** `src/scripts/04_htvs_biophysics.py` reads `snakemake.input/output` keys that do not exist in the `Snakefile` → unreachable. Internally: hard-coded `5o66`/pocket centers, `infer_protein_type` always returns `"efflux_pump"` (`:656`), `pocket_center_from_mutations` never called (`:148`), naive ΔΔG `mut − wt`, asymmetric rigid-WT vs flexible-mutant docking, unreachable `return` after `raise` (`:326`).
- **Placeholder modules exported as real:** `alignment.SequenceAligner.{global,local}_alignment` return hard-coded `{"score":0,…}` (`alignment.py:44,66`); `variant_detector.detect_variants` → `[]`; `translator.find_orfs` → `[]`; `coordinate_parser` → `[]`; `biophysics/__init__.py` empty.
- **`publish.yml` broken:** installs on Python 3.9 while `pyproject.toml` requires `>=3.10`.
- **CI tests little:** `tests.yml` runs only `tests/unit/` and dry-runs with `local_genomes="data/genomes"` (nonexistent), never exercising the shipped config or the Tier-1 blockers.
- **Version drift:** `__init__.py` = `2.0.0`, `pyproject.toml` = `2.1.0`; deps omit `snakemake`, include unused `requests`.
- **Dead/bloat:** `config/settings.yaml` (97 unused lines); 35 MB `utility scripts/bvbrc-cli` committed; `drug_mapping.json` maps gyrA/parC/marR/acrR but **not acrA/acrB/tolC** (3 of 5 real targets).

---

## Corrections to `review.md`

| `review.md` claim | Reality |
|---|---|
| Dry-run crashes on hard-coded `data/5o66.pdb` | No such reference in the DAG; `snakemake -n` **passes** |
| `config.yaml` has a `default_pdb` key | It does not |
| `variant_min_identity_percent` missing from config | Present (`config.yaml:12`), read at `02b:57` |
| Snakefile has a Phase-3 biophysics rule (`md_stiffness`, `params.ligand`, `center_x`) | No biophysics rule exists |
| `gyrA` sentinel/skip bug in extraction | No such code; skip logic is `resolve_reference`-based (`02a:111`) |
| `.fasta` vs `.faa` mismatch breaks refs | Both accepted (`_REF_EXTENSIONS`); UniProt fetch saves canonical `.faa` |
| No reference seeding exists | Full seeding path exists (`seed_references_from_dir`, `_ensure_references_exist`) |

Correctly identified by `review.md`: the `.head(5)` truncation, shared-directory integration tests, and the orphaned-biophysics smell (though attributed to the wrong mechanisms).

---

## Prioritized fixes

1. Fix the tuple order and source label at `variant_caller.py:610`; stop labeling non-hits `"Clinical DB"`.
2. Ship/load a real resistance DB from `config/`, or make an empty DB a loud warning.
3. Point config at real inputs (`acr_targets.txt` + `reference_seed_dir: data/results/refs`).
4. Unify reference resolution (use `resolve_reference` in `_call_variants_single`).
5. Fix VUS counters and the case-sensitive DB lookup.
6. Fix the integration test identity/threshold; isolate test dirs.
7. Pin CI to 3.10+, run integration tests, dry-run the shipped config.
8. Decide the biophysics phase's fate; reconcile versions; drop `settings.yaml` and the vendored CLI.
9. Replace co-occurrence with a real association test + phylogeny correction, or reframe it explicitly as descriptive-only; make the severity fallback a true sentinel; remove/annotate the `.head(5)` cap.
