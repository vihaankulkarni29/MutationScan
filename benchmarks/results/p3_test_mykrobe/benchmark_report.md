# AMR benchmark results

## Summary

| tool | version | sensitivity | precision | F1 | epistasis F1 | dN/dS r | lineage | runtime | genomes/hr |
|------|---------|-------------|-----------|-----|--------------|---------|---------|---------|-----------|
| Mykrobe | mykrobe v0.13.0 | 0.0% | 0.0% | 0.000 | 0.000 | 0.000 | 0.0% | 22.5s | 1917 |

## Mutation calling

| tool | TP | FP | FN | TN | novel recall |
|------|----|----|----|----|--------------|
| Mykrobe | 0 | 0 | 0 | -- | -- |

## Epistasis

| tool | TP | FP | FN | precision | recall | F1 |
|------|----|----|----|-----------|--------|-----|
| Mykrobe | 0 | 0 | 0 | 0.0% | 0.0% | 0.000 |

## dN/dS

| tool | Pearson | Spearman | MAE | genes |
|------|---------|----------|-----|-------|
| Mykrobe | 0.000 | 0.000 | 0.000 | 0 |

## Performance

| tool | runtime | peak RSS | genomes | genomes/hr |
|------|---------|----------|---------|------------|
| Mykrobe | 22.5s | 82 MiB | 12 | 1917 |
