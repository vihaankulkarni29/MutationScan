# AGENTS.md

## Lint + typecheck
```bash
ruff check src/ tests/
ruff format --check src/ tests/
mypy src/mutation_scan/
```

## Tests
```bash
python -m pytest tests/unit/ -q
python -m pytest tests/integration/ -v
```

## Key paths
- Source: `src/mutation_scan/`
- Tests: `tests/unit/`, `tests/integration/`
- Workflow: `Snakefile`, `config/config.yaml`
- Workflow scripts: `workflow/scripts/0*.py`
