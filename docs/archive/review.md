# MutationScan Pipeline Audit

**Date:** 2026-08-30  
**Branch:** main  
**Scope:** Full pipeline audit — lacunae, dead code, hardcoded noise, inconsistencies

---

## Executive Summary

The pipeline **does not run end-to-end on this machine right now**:

- Snakemake dry-run fails: missing `data/5o66.pdb`
- Phase 1 has a **reference bootstrap gap** on fresh runs
- Tests: **18/19 pass**, 1 integration test fails (cross-test pollution + identity filter)
- CI only runs unit tests and a partial dry-run — it would **not catch** most of these

---

## 1. Pipeline Blockers (won't run or will silently fail)

| Issue | Where | Detail |
|-------|--------|--------|
| **Missing default PDB** | `config/config.yaml`, dry-run | `default_pdb: data/5o66.pdb` — file does not exist. `data/AF-AcrB.pdb` exists but isn't wired in. |
| **Reference bootstrap chicken-and-egg** | `02a_extract_proteins.py`, `tblastn_extractor.py` | Fresh runs write refs to `data/output/{job}/refs/` which starts **empty**. Legacy refs live in `data/results/refs/` but **nothing copies them**. Without `uniprot_taxid`, extraction fails at pre-flight. |
| **`.fasta` vs `.faa` mismatch** | `tblastn_extractor.py` | UniProt auto-fetch saves `{gene}.fasta`. Pre-flight accepts `.fasta` OR `{gene}_WT.faa`. But tblastn only loads `*.faa` — fetched UniProt refs are **invisible** to tblastn. |
| **Empty resistance DB** | `variant_caller.py` | `resistance_db.json` expected in per-run `refs/` — **not in repo**. Every mutation falls through to fallback logic; "Resistant" status is effectively unreachable in production. |
| **`local_genomes` has data but no refs path** | Snakefile | ~1000+ `.fna` files exist, but without seeded refs the pipeline stops at Phase 1a. |

---

## 2. Real Bugs (wrong output, not just missing polish)

### Variant caller — swapped Status/Phenotype when ML is off

`src/mutation_scan/analysis/variant_caller.py` (lines 606–610):

```python
if not self.enable_ml:
    drug = self.drug_mapping.get(gene_name.lower())
    phenotype = f"{drug} resistance (VUS)" if drug else "Variant of Unknown Significance"
    return (phenotype, "N/A", "N/A", None, "Clinical DB")
```

Should return `(status, phenotype, pdb, score, source)` but returns **phenotype in the Status column** and `"N/A"` as Phenotype. Production sets `enable_ml=False` in `02b_call_variants.py`, so this path is always hit for unknown mutations.

### `skip_extraction` checks wrong gene

`src/scripts/02a_extract_proteins.py` (line 66):

```python
if protein_files and (refs_dir / "gyrA_WT.faa").exists():
```

Pipeline targets are `acrA/acrB/tolC/acrR/marR`. Skip logic checks legacy **gyrA** — wrong sentinel for this project.

### Gene name case sensitivity in resistance DB lookup

`src/mutation_scan/analysis/variant_caller.py` (line 565):

```python
if gene_name not in self.resistance_db:
```

No `.lower()` normalization. Filename parsing yields `acrB`; DB keys would need exact match.

### Integration test failure (confirmed)

`test_multiple_mutations` fails because:

1. Prior tests leave extra `.faa` files in shared `data/test_variant_caller/proteins/`
2. 75% identity < 80% threshold → mutations silently dropped

---

## 3. Config / Docs / Snakefile Mismatches

| Documented (README) | Actual code | Impact |
|---------------------|-------------|--------|
| `pocket_center_x/y/z` | `center_x/y/z` in `04_htvs_biophysics.py` | User overrides silently ignored |
| `exhaustiveness: 16` | `biophysics_exhaustiveness` | README key does nothing |
| `variant_min_identity_percent` | Works, but **not in `config.yaml`** | Only via CLI override |
| `config/settings.yaml` | 97 lines of NCBI/ABRicate/BLAST settings | **Completely unused** by Snakefile |
| Snakefile `params.ligand`, `center_x`, `md_stiffness` | Script reads `snakemake.config`, not `snakemake.params` | Snakefile params are dead |
| CI dry-run uses `local_genomes="data/genomes"` | Production default is `data/local_genomes` | CI validates a different path |
| `pyproject.toml` v2.1.0 | `__init__.py` v2.0.0 | Version drift |

---

## 4. Phase-by-Phase Lacunae

### Phase 1a — Protein Extraction

| Issue | Detail |
|-------|--------|
| **No ref seeding** | `data/results/refs/*.faa` exist but aren't copied to `OUT_DIR/refs` |
| **tblastn takes only best HSP** | `-max_target_seqs 1 -max_hsps 1` — paralogs/duplicates ignored |
| **Silent per-gene failure** | Failed genes logged but run continues; downstream gets partial data |
| **Gene filter is substring match** | `tg in f.stem.lower()` — `acr` would match `acrA`, `acrB`, `acrR` loosely |
| **sys.path hack** | All scripts use `sys.path.append(...)` instead of proper package install |

### Phase 1b — Variant Calling

| Issue | Detail |
|-------|--------|
| **No resistance DB shipped** | All mutations → VUS/unknown path with buggy tuple |
| **ML module referenced but missing** | `mutation_scan.ml_predictor` — import always fails, logs warning noise |
| **80% identity filter** | Drops alignments aggressively; may miss real partial genes (test proves this) |
| **Drug mapping missing pipeline targets** | `acrA`, `acrB`, `tolC` not in `drug_mapping.json` — only `acrR`, `marR` mapped |
| **Global alignment on partial proteins** | Can produce spurious substitutions at alignment boundaries |
| **No indel calling** | Only substitutions; insertions/deletions ignored by design but never documented in output |

### Phase 2 — Epistasis / ControlScan

| Issue | Detail |
|-------|--------|
| **Arbitrary top-5 cap** | `.head(5)` — biophysics only ever sees 5 networks, rest discarded |
| **Silent score fallback** | Parse failure → `ControlScan_Score = 1.0` — masks bad mutation strings |
| **No statistical test** | Co-occurrence = raw count; no Fisher's exact, no FDR, no cohort normalization |
| **Composite score = freq × severity** | Favors common mild pairs over rare severe ones; no tunable weighting |
| **Phase label confusion** | Script logs "Phase 2" and "Phase 3" but it's one Snakemake rule |
| **External package fallback** | `from controlscan.scorer import MutationScorer` — package doesn't exist in repo |
| **Plots are 2-node graphs only** | NetworkX spring layout on single edges — presentation artifact, not analysis |
| **Gene lowercasing** | `df['Gene'] = df['Gene'].str.lower()` — may break PDB/chain lookups downstream |

### Phase 3 — Biophysics / Docking

| Issue | Detail |
|-------|--------|
| **"OpenMM Dynamics" in Snakefile rule name** | No OpenMM MD runs. `biophysics_deep_relaxed_md` only logs a warning |
| **`proteins_dir` input declared, never used** | Dead Snakefile input |
| **`pocket_center_from_mutations()` defined, never called** | Dead code; fixed pocket used instead |
| **Three conflicting pocket centers** | AcrB preset `(18,-24,5)`, H596 CA fallback, generic `(18.5,55.2,20.1)` |
| **Default ligand auto-downloads ciprofloxacin** | From PubChem at runtime if missing — may not match user's antibiotic study |
| **Smina binary auto-download** | From SourceForge at runtime — fragile, non-reproducible |
| **`infer_protein_type()` always returns `efflux_pump`** | Interpretation labels are cosmetic, not target-aware |
| **ΔΔG = mut_affinity − wt_affinity** | Sign convention undocumented; Vina scores aren't true ΔG |
| **Rigid receptor for WT, flex sidechains for mutant** | Asymmetric protocol biases ΔΔG |
| **`md_stiffness` in Snakefile** | Never read anywhere |
| **`stiffness` param passed** | Dead |
| **QC clash threshold hardcoded** | `2.0 Å`, max 20 pairs — not configurable |
| **README says MVBM/CDSS** | `biophysics/cdss_engine.py` and package namespace are **empty/missing** |

---

## 5. Dead Code / Stub Modules (exported but unused)

These are in the package, exported in `__init__.py`, but **not called by the production DAG**:

| Module | State |
|--------|-------|
| `analysis/variant_detector.py` | Stub — `return []` |
| `analysis/alignment.py` | Stub — `return {"score": 0, ...}` |
| `core/coordinate_parser.py` | Stub |
| `core/translator.py` | `find_orfs()` stub; only `translate()` tested |
| `core/gene_finder.py` | Full ABRicate/BLAST impl — legacy, unused |
| `core/sequence_extractor.py` | Coordinate-based extraction — replaced by tblastn |
| `core/reference_builder.py` | NCBI ref fetcher — unused |
| `core/genome_processor.py` | Unused |
| `utils/config_parser.py` | Tested but not used by pipeline scripts |
| `biophysics/__init__.py` | Empty namespace |

Roughly **half the `mutation_scan` package is legacy scaffolding**.

---

## 6. Hardcoded / Magic Values

| Value | Location | Concern |
|-------|----------|---------|
| `80.0` % identity | variant caller default | Drops 75% alignments silently |
| `.head(5)` | epistasis script | Arbitrary network cap |
| `1.0` ControlScan fallback | epistasis script | Hides parse errors |
| `seed=42` | spring_layout + docking | Fixed seed, not configurable except `biophysics_seed` |
| `BOX_SIZE = 25.0` | biophysics | Hardcoded docking box |
| `EXHAUSTIVENESS = 16` | biophysics | Not wired to README's `exhaustiveness` key |
| `DOCKING_SEED_DEFAULT = 42` | biophysics | |
| `QC_CLASH_DISTANCE = 2.0 Å` | biophysics | |
| `CIPRO_URL` PubChem CID 2764 | biophysics | Default ligand |
| `find_ca_coord(pdb, "A", 596)` | biophysics | AcrB H596 hardcoded fallback |
| `TARGET_POCKET_CENTERS["acrb"]` | biophysics | Only one target defined |
| `default_chain = "A"` | biophysics | No validation against PDB |
| `antibiotic = "Ciprofloxacin"` | variant caller default | Irrelevant to Acr targets |
| `gyrA_WT.faa` | skip_extraction check | Wrong gene |
| `tblastn` timeout 60s | extractor | May fail on large genomes |
| `time.sleep(0.5)` | UniProt fetch | Rate limit band-aid |

---

## 7. CI / Test Gaps

| Gap | Impact |
|-----|--------|
| CI runs **only** `tests/unit/` | Integration test failure not caught |
| CI dry-run doesn't override `default_pdb` | Would also fail on missing PDB in clean checkout |
| No end-to-end smoke test with fixture genomes | Bootstrap ref problem undetected |
| No test for epistasis or biophysics stages | Largest scientific risk unvalidated |
| `pyproject.toml` missing `snakemake` dep | Only works if conda env used |
| No `[tool.setuptools.packages.find]` | Package discovery relies on implicit behavior |

---

## 8. Repo / Hygiene Noise

- **`utility scripts/bvbrc-cli/`** — huge third-party bundle, untracked, unrelated to production DAG
- **`config/settings.yaml`** — 90+ lines of dead configuration
- **`data/results/`** — legacy run artifacts mixed with placeholder structure
- **Package `__init__.py` docstring** still claims "Genome downloading" — contradicts local-only design
- **`TOOLKIT_USAGE.md`** references runs/paths (`Ciproflaxcin_Run` typo) not in production config

---

## 9. Scientific / Interpretation Lacunae

These aren't code bugs but **methodological gaps** that produce noisy or misleading outputs:

1. **No linkage between mutation severity and clinical resistance** for Acr targets (no DB entries)
2. **Epistasis pairs are co-occurrence, not epistasis** — no interaction term, no background rate
3. **Docking ΔΔG on efflux pump** — substrate binding site mutations ≠ efflux mechanism; interpretation labels (`REDUCED_SUBSTRATE_AFFINITY`) are heuristic
4. **tblastn translated subject sequence** may include frameshift artifacts at contig edges — no QC on translation quality beyond identity %
5. **Top-5 networks → biophysics** — highest composite score may miss biologically critical singleton mutations
6. **No sample metadata** in pipeline — can't stratify by lineage, geography, or resistance phenotype

---

## Priority Fix List

### P0 — pipeline actually runs

1. Ship or copy seed refs (`*_WT.faa`) into `OUT_DIR/refs` on first run
2. Fix `.fasta`/`.faa` inconsistency in tblastn extractor
3. Point `default_pdb` to an existing file (`AF-AcrB.pdb`) or ship `5o66.pdb`
4. Fix `_fallback_to_ml` tuple swap bug

### P1 — output correctness

5. Add `resistance_db.json` for Acr targets (or remove dead "Resistant" path)
6. Add `acrA/acrB/tolC` to `drug_mapping.json`
7. Align README config keys with code (`center_x`, `biophysics_exhaustiveness`)
8. Remove or wire `skip_extraction` gyrA check to actual targets

### P2 — noise reduction

9. Delete or quarantine stub modules
10. Remove dead Snakefile params / unused inputs
11. Make epistasis cap configurable (or pass all networks to biophysics)
12. Fix integration test isolation; add to CI
