# AMR benchmark results

## Summary

| tool | version | sensitivity | precision | F1 | epistasis F1 | dN/dS r | lineage | runtime | genomes/hr |
|------|---------|-------------|-----------|-----|--------------|---------|---------|---------|-----------|
| AMRFinderPlus | 4.2.7 | 0.0% | 0.0% | 0.000 | 0.000 | 0.000 | 0.0% | 247.2s | 175 |

## Mutation calling

| tool | TP | FP | FN | TN | novel recall |
|------|----|----|----|----|--------------|
| AMRFinderPlus | 0 | 0 | 0 | -- | -- |

## Epistasis

| tool | TP | FP | FN | precision | recall | F1 |
|------|----|----|----|-----------|--------|-----|
| AMRFinderPlus | 0 | 0 | 0 | 0.0% | 0.0% | 0.000 |

## dN/dS

| tool | Pearson | Spearman | MAE | genes |
|------|---------|----------|-----|-------|
| AMRFinderPlus | 0.000 | 0.000 | 0.000 | 0 |

## Performance

| tool | runtime | peak RSS | genomes | genomes/hr |
|------|---------|----------|---------|------------|
| AMRFinderPlus | 247.2s | 481 MiB | 12 | 175 |
