# AMR benchmark results

## Summary

| tool | version | sensitivity | precision | F1 | epistasis F1 | dN/dS r | lineage | runtime | genomes/hr |
|------|---------|-------------|-----------|-----|--------------|---------|---------|---------|-----------|
| Kleborate | Kleborate v3.3.0b1 | 0.0% | 0.0% | 0.000 | 0.000 | 0.000 | 0.0% | 120.3s | 359 |

## Mutation calling

| tool | TP | FP | FN | TN | novel recall |
|------|----|----|----|----|--------------|
| Kleborate | 0 | 91 | 0 | -- | -- |

## Epistasis

| tool | TP | FP | FN | precision | recall | F1 |
|------|----|----|----|-----------|--------|-----|
| Kleborate | 0 | 0 | 0 | 0.0% | 0.0% | 0.000 |

## dN/dS

| tool | Pearson | Spearman | MAE | genes |
|------|---------|----------|-----|-------|
| Kleborate | 0.000 | 0.000 | 0.000 | 0 |

## Performance

| tool | runtime | peak RSS | genomes | genomes/hr |
|------|---------|----------|---------|------------|
| Kleborate | 120.3s | 363 MiB | 12 | 359 |
