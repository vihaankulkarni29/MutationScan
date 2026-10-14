# AMR benchmark results

## Summary

| tool | version | sensitivity | precision | F1 | epistasis F1 | dN/dS r | lineage | runtime | genomes/hr |
|------|---------|-------------|-----------|-----|--------------|---------|---------|---------|-----------|
| MutationScan | mutationscan 1.0.0 | 0.0% | 0.0% | 0.000 | 0.000 | 0.000 | 0.0% | 34.1s | 1266 |

## Mutation calling

| tool | TP | FP | FN | TN | novel recall |
|------|----|----|----|----|--------------|
| MutationScan | 0 | 39 | 0 | -- | -- |

## Epistasis

| tool | TP | FP | FN | precision | recall | F1 |
|------|----|----|----|-----------|--------|-----|
| MutationScan | 0 | 12 | 0 | 0.0% | 0.0% | 0.000 |

## dN/dS

| tool | Pearson | Spearman | MAE | genes |
|------|---------|----------|-----|-------|
| MutationScan | 0.000 | 0.000 | 0.000 | 0 |

## Performance

| tool | runtime | peak RSS | genomes | genomes/hr |
|------|---------|----------|---------|------------|
| MutationScan | 34.1s | 272 MiB | 12 | 1266 |
