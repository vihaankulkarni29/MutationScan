# Phenotypic concordance

Does a called mutation predict the phenotype the laboratory measures? Three measurements, in increasing order of how much they assume about the biology.

## Gene-drug association (unsupervised)

Carrying a mutation in a gene vs. being phenotypically resistant, Fisher exact two-sided, Benjamini-Hochberg across the gene x drug family. No curated knowledge used.

| gene | drug | carriers | resistant | both | n | odds ratio | p | q | sig |
|------|------|----------|-----------|------|---|------------|---|---|-----|
| acrA | ceftriaxone | 1 | 3 | 1 | 3 | -- | 1.0000 | 1.0000 | no |
| marR | ceftriaxone | 1 | 3 | 1 | 3 | -- | 1.0000 | 1.0000 | no |
| tolC | ceftriaxone | 1 | 3 | 1 | 3 | -- | 1.0000 | 1.0000 | no |
| acrA | ciprofloxacin | 2 | 9 | 2 | 9 | -- | 1.0000 | 1.0000 | no |
| acrB | ciprofloxacin | 3 | 9 | 3 | 9 | -- | 1.0000 | 1.0000 | no |
| acrR | ciprofloxacin | 1 | 9 | 1 | 9 | -- | 1.0000 | 1.0000 | no |
| marR | ciprofloxacin | 9 | 9 | 9 | 9 | -- | 1.0000 | 1.0000 | no |
| tolC | ciprofloxacin | 1 | 9 | 1 | 9 | -- | 1.0000 | 1.0000 | no |
| marR | gentamicin | 3 | 3 | 3 | 3 | -- | 1.0000 | 1.0000 | no |

**0 of 9 tests significant** at FDR < 0.05.

## Burden vs. MIC (unsupervised)

Spearman between the number of mutated genes a genome carries and its MIC on the log2 scale. Uses the concentration, not a dichotomised threshold, so it does not depend on the breakpoint being right.

| drug | n | carriers | rho | p | sig |
|------|---|----------|-----|---|-----|
| cefepime | 2 | 0 | 0.000 | 1.0000 | no |
| ceftriaxone | 3 | 1 | 0.000 | 1.0000 | no |
| ciprofloxacin | 9 | 9 | 0.779 | 0.0134 | **yes** |
| gentamicin | 3 | 3 | 0.000 | 1.0000 | no |

## Concordance vs. curated rules

VME = reported susceptible, measured resistant. ME = reported resistant, measured susceptible. They are not interchangeable.

| metric | value |
|--------|-------|
| evaluable | 17 |
| unknown phenotype | 0 |
| agreement | 0.765 |
| very major error | 0.235 (4) |
| major error | 0.000 (0) |
| minor error | 0.000 (0) |

| drug | n | agreement | VME | ME | minor |
|------|---|-----------|-----|----|-------|
| cefepime | 2 | 0.000 | 1.000 | 0.000 | 0.000 |
| ceftriaxone | 3 | 0.333 | 0.667 | 0.000 | 0.000 |
| ciprofloxacin | 9 | 1.000 | 0.000 | 0.000 | 0.000 |
| gentamicin | 3 | 1.000 | 0.000 | 0.000 | 0.000 |

## Data sources

- MIC: BV-BRC (formerly PATRIC) AMR phenotype service, or a curated table. `benchmarks/phenotypic.py` reads accession, drug, MIC and organism.
- Breakpoints: EUCAST `clinical_breakpoints` (open access) or CLSI M100 (paywalled). Both parse the same way.
- Rules: a curated gene -> drug mapping, one row per resistance relationship.
