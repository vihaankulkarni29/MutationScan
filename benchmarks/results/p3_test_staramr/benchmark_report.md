# AMR benchmark results

## Summary

| tool | version | sensitivity | precision | F1 | epistasis F1 | dN/dS r | lineage | runtime | genomes/hr |
|------|---------|-------------|-----------|-----|--------------|---------|---------|---------|-----------|
| StarAMR | staramr 0.12.3 | 0.0% | 0.0% | 0.000 | 0.000 | 0.000 | 0.0% | 36.7s | 1178 |

## Mutation calling

| tool | TP | FP | FN | TN | novel recall |
|------|----|----|----|----|--------------|
| StarAMR | 0 | 0 | 0 | -- | -- |

## Epistasis

| tool | TP | FP | FN | precision | recall | F1 |
|------|----|----|----|-----------|--------|-----|
| StarAMR | 0 | 0 | 0 | 0.0% | 0.0% | 0.000 |

## dN/dS

| tool | Pearson | Spearman | MAE | genes |
|------|---------|----------|-----|-------|
| StarAMR | 0.000 | 0.000 | 0.000 | 0 |

## Performance

| tool | runtime | peak RSS | genomes | genomes/hr |
|------|---------|----------|---------|------------|
| StarAMR | 36.7s | -- | 12 | 1178 |
