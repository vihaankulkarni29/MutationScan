# AGENTS.md

## Lint + typecheck
```bash
ruff check src/ tests/ benchmarks/
ruff format --check src/ tests/ benchmarks/
mypy src/mutation_scan/
```

## Tests
```bash
python -m pytest tests/unit/ -q
python -m pytest tests/integration/ -v
```

`tests/integration/` needs `tblastn` on PATH -- it runs the real binary rather
than skipping. In `tests/unit/`, the `@requires_tblastn` tests skip when BLAST+
is absent, so a green unit run proves nothing about extraction unless you check
the skip count.

## Benchmarks
```bash
python benchmarks/parallel_bench.py
```
Startup cost, scaling, serial/parallel parity, and the chunk-size and
serial-fallback thresholds. Not part of CI; run it after touching
`src/mutation_scan/parallel.py`.

## Key paths
- Source: `src/mutation_scan/`
- Tests: `tests/unit/`, `tests/integration/`
- Workflow: `Snakefile`, `config/config.yaml`
- Workflow scripts: `workflow/scripts/0*.py`
- Parallelism: `src/mutation_scan/parallel.py` (executor),
  `src/mutation_scan/parallel_plan.py` (what a run will actually use -- the
  run summary and `config-check` both read this, so change it in one place)
- Species routing: `src/mutation_scan/species.py` (panel parsing + genome
  assignment), and the `species:` block in `config/config.yaml`. A multi-species
  cohort MUST configure this -- one reference panel across species silently calls
  every interspecies difference a mutation, and RBH does not catch it.
