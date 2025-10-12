## MutationScan Test Plan and Acceptance Criteria

This document defines how we will validate the entire MutationScan pipeline end-to-end and per-domino, including structure, dependencies, integration, correctness, and user experience.

### Goals
- Ensure all dependencies install cleanly on supported environments.
- Verify repository structure and required files/manifests exist at each domino stage.
- Confirm each domino tool runs standalone (help + smoke) and full with sample inputs.
- Validate logic produces correct, non-empty, consistent outputs and handoffs.
- Prove tools integrate in the orchestrated pipeline with proper error handling.
- Audit user experience (messages, progress, troubleshooting) for clarity.

---

## Scope
- Platforms:
  - Ubuntu (native) – primary supported environment
  - Windows 10/11 with WSL Ubuntu – recommended for ABRicate + native tools
  - Windows native – allowed for mock-only/demo runs (ABRicate not supported natively)
- Python versions: 3.9 – 3.11
- Pipeline dominos: Harvester → Annotator → Extractor → Aligner → Analyzer → CoOccurrence → Reporter

Out of scope: Large-scale benchmarking, non-x86 architectures, bespoke HPC clusters.

---

## Environments and Prerequisites
- Python 3.9+ virtual environment or conda env
- Git, pip, build tools (pip>=23, setuptools>=61, wheel)
- Orchestrator and tools from this repository (up to date on main)
- External tools/packages installed in the same environment:
  - ncbi-genome-extractor (GitHub)
  - ABRicate (Linux/WSL) + databases setup (abricate --setupdb)
  - abricate-automator (GitHub)
  - FastaAAExtractor (GitHub) – CLI-capable
  - wildtype-aligner (GitHub)
  - subscan-analyzer (GitHub)

Note: On Windows native runs, use the mock extractor and skip ABRicate-dependent tests; full runs should use Ubuntu/WSL.

---

## Acceptance Criteria (Overview)
1. Dependencies
   - Editable install of `subscan` succeeds (`pip install -e subscan/`).
   - All Python imports succeed for core libs and installed GitHub tools.
   - ABRicate available on Ubuntu/WSL (`abricate --version` and `--setupdb` done).

2. Structure
   - Required scripts present: `subscan/tools/run_*.py` for all dominos.
   - Orchestrator present and executable; manifests named as per contract.

3. Domino functionality (standalone)
   - Each `run_*.py` supports `--help` and a minimal smoke run, exiting 0 on success.
   - For Extractor: `fasta_aa_extractor --help` or `python -m fasta_aa_extractor --help` works.

4. Integration and handoff
   - Each domino produces its manifest/file with expected schema/keys.
   - Next domino accepts previous manifest and completes.

5. Correctness checks
   - Non-empty outputs where expected, file counts within reasonable bounds.
   - Basic schema validation and value sanity (e.g., FASTA headers, TSV columns).

6. UX quality
   - Clear progress indicators per domino.
   - Actionable error messages with remediation hints.
   - README/docs contain accurate instructions and troubleshooting.

7. End-to-end
   - Dry-run completes and produces all placeholder manifests and final HTML shell.
   - Small full-run completes and produces a valid `final_report.html`.

---

## Test Matrix

- Smoke tests (fast):
  - Tools: `--help` runs for each domino and external CLI
  - Extractor: `subscan/tools/smoke_test_extractor.py` PASS
  - Orchestrator dry-run produces all placeholder manifests

- Full tests (sample data):
  - Run Harvester with 1–3 accessions; downstream dominos generate real outputs
  - Validate manifests and outputs at each stage
  - Produce final report

---

## Per-Domino Contracts and Checks

For each domino, we define inputs, outputs, CLI, and checks.

### Domino 1: Harvester
- Input: `--accessions <file>`, `--email <string>`
- Output: `01_harvester_results/genome_manifest.json`, FASTAs, metadata CSV
- Smoke: `python subscan/tools/run_harvester.py --help`
- Full: small accession list (1–3 entries)
- Checks:
  - Manifest JSON exists and includes genomes list with FASTA paths
  - Metadata CSV present and non-empty

### Domino 2: Annotator
- Input: `--manifest genome_manifest.json`
- Output: `02_annotator_results/annotation_manifest.json`, ABRicate TSVs
- Smoke: `python subscan/tools/run_annotator.py --help`
- Full: requires ABRicate on Ubuntu/WSL
- Checks:
  - Manifest lists genomes with `amr_card_results` TSVs
  - TSVs have header lines and at least one data row for known genes

### Domino 3: Extractor
- Input: `--manifest annotation_manifest.json`, `--gene-list gene_list.txt`
- Output: `03_extractor_results/protein_manifest.json`, `*_proteins.faa`
- Smoke: `python subscan/tools/smoke_test_extractor.py`
- Full: real FastaAAExtractor installed (or mock for Windows-only demo)
- Checks:
  - Protein FASTA files created, contain `>` headers
  - Manifest enumerates protein files with counts and sizes

### Domino 4: Aligner
- Input: `--manifest protein_manifest.json` (+ optional `--sepi-species`)
- Output: `04_aligner_results/alignment_manifest.json` (+ alignments)
- Smoke: `python subscan/tools/run_aligner.py --help`
- Full: requires wildtype-aligner package
- Checks:
  - Alignment files exist and non-empty; manifest lists them

### Domino 5: Analyzer
- Input: `--manifest alignment_manifest.json`
- Output: `05_analyzer_results/analysis_manifest.json`, metrics/tables
- Smoke: `python subscan/tools/run_analyzer.py --help`
- Full: requires subscan-analyzer package
- Checks:
  - Manifest present; metrics/tables non-empty with reasonable values

### Domino 6: CoOccurrence
- Input: `--manifest analysis_manifest.json`
- Output: `06_cooccurrence_results/cooccurrence_manifest.json`, matrices/plots
- Smoke: `python subscan/tools/run_cooccurrence_analyzer.py --help`
- Full: uses outputs from Analyzer
- Checks:
  - Manifest present; co-occurrence data non-empty

### Domino 7: Reporter
- Input: `--manifest cooccurrence_manifest.json`
- Output: `07_reporter_results/`, `final_report.html`
- Smoke: `python subscan/tools/run_reporter.py --help`
- Full: generates final HTML report
- Checks:
  - `final_report.html` exists, opens without errors (basic HTML contains expected sections)

---

## Orchestrator Tests

### Dry-run (fast):
```bash
python subscan/tools/run_pipeline.py \
  --accessions data_input/accession_list.txt \
  --gene-list data_input/gene_list.txt \
  --email you@example.com \
  --output-dir /tmp/ms_dryrun \
  --dry-run
```
Acceptance:
- Creates run directory with 01..06 result folders
- Writes placeholder manifests for each domino
- Prints progress for all 7 stages and finishes cleanly

### Full-run (sample data, Ubuntu/WSL):
```bash
python subscan/tools/run_pipeline.py \
  --accessions subscan/sample_data/sample_accessions.txt \
  --gene-list subscan/sample_data/gene_list.txt \
  --email you@example.com \
  --output-dir /tmp/ms_full \
  --threads 4
```
Acceptance:
- Each domino completes; required manifests/files exist
- `final_report.html` generated in the run directory

---

## Dependency and External Tool Checks

Quick checks to embed in validation scripts:
- Python imports: `requests`, `tqdm`, `pandas`, `Bio`, `jinja2`, `plotly`
- Console/module availability:
  - `fasta_aa_extractor --help` or `python -m fasta_aa_extractor --help`
  - `abricate --version` and `abricate --setupdb` (Ubuntu/WSL)
  - `python -c "import importlib; importlib.import_module('ncbi_genome_extractor')"`
  - `... 'wildtype_aligner'`, `'subscan'`, `'subscan_analyzer'`

Acceptance:
- All core imports succeed; external CLIs respond with 0 exit code

---

## Output Verification Harness (Principles)
- Presence checks: files exist at documented paths
- Schema checks: JSON required keys, TSV columns, FASTA headers
- Sanity checks: counts > 0 where expected, sizes reasonable, no empty critical files
- Handoff checks: previous manifest path referenced in next manifest metadata

---

## User Experience Audit
- Messages are ASCII-safe and clear
- Progress shows domino index and description
- Errors include actionable next steps (what to install/configure)
- README/docs match current CLI flags and workflows

Acceptance:
- No ambiguous errors; missing tool errors list exact remediation

---

## How to Run All Checks (Operator Quick Start)

1) Environment and install
```bash
python -m venv ~/mutationscan_env && source ~/mutationscan_env/bin/activate
python -m pip install -U pip setuptools wheel
cd subscan && pip install -e . && cd ..
```

2) Smoke checks
```bash
python subscan/tools/run_harvester.py --help
python subscan/tools/run_annotator.py --help
python subscan/tools/smoke_test_extractor.py
python subscan/tools/run_aligner.py --help
python subscan/tools/run_analyzer.py --help
python subscan/tools/run_cooccurrence_analyzer.py --help
python subscan/tools/run_reporter.py --help
```

3) Orchestrator dry-run
```bash
python subscan/tools/run_pipeline.py \
  --accessions data_input/accession_list.txt \
  --gene-list data_input/gene_list.txt \
  --email you@example.com \
  --output-dir /tmp/ms_dryrun \
  --dry-run
```

4) Full run (Ubuntu/WSL)
```bash
python subscan/tools/run_pipeline.py \
  --accessions subscan/sample_data/sample_accessions.txt \
  --gene-list subscan/sample_data/gene_list.txt \
  --email you@example.com \
  --output-dir /tmp/ms_full \
  --threads 4
```

---

## Reporting
- For each run, collect:
  - Console logs
  - Manifest JSONs
  - Any error stack traces
- Summarize PASS/FAIL per domino and for end-to-end

---

## Next Steps
- Implement automated validators and smoke tests per this plan
- Wire CI to run smoke tests on push/PR
- Iterate on UX based on audit findings
