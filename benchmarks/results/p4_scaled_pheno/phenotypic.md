# Phenotypic concordance

Does a called mutation predict the phenotype the laboratory measures? Three measurements, in increasing order of how much they assume about the biology.

## Gene-drug association (unsupervised)

Carrying a mutation in a gene vs. being phenotypically resistant, Fisher exact two-sided, Benjamini-Hochberg across the gene x drug family. No curated knowledge used.

| gene | drug | carriers | resistant | both | n | odds ratio | p | q | sig |
|------|------|----------|-----------|------|---|------------|---|---|-----|
| _none_ | | | | | | | | | |

**0 of 0 tests significant** at FDR < 0.05.

## Burden vs. MIC (unsupervised)

Spearman between the number of mutated genes a genome carries and its MIC on the log2 scale. Uses the concentration, not a dichotomised threshold, so it does not depend on the breakpoint being right.

| drug | n | carriers | rho | p | sig |
|------|---|----------|-----|---|-----|
| _none_ | | | | | |

## Concordance vs. curated rules

VME = reported susceptible, measured resistant. ME = reported resistant, measured susceptible. They are not interchangeable.

| metric | value |
|--------|-------|
| evaluable | 0 |
| unknown phenotype | 17 |
| agreement | 0.000 |
| very major error | 0.000 (0) |
| major error | 0.000 (0) |
| minor error | 0.000 (0) |

## Data sources

- MIC: BV-BRC (formerly PATRIC) AMR phenotype service, or a curated table. `benchmarks/phenotypic.py` reads accession, drug, MIC and organism.
- Breakpoints: EUCAST `clinical_breakpoints` (open access) or CLSI M100 (paywalled). Both parse the same way.
- Rules: a curated gene -> drug mapping, one row per resistance relationship.
