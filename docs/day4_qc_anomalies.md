# Day 4 — QC anomalies and observations

This document records descriptive QC findings for the primary PD1 cohort.

No sample was removed, no clinical value was imputed, and no predictive inference was performed.

## Expression matrix integrity

- Matrix dimensions: 17002 genes × 36 PD1 samples.
- Missing count values: 0.
- Negative count values: 0.
- Non-integer-like count values: 0.

## Potential sample-level QC outliers

### PD1_9

- Library size: 4269738.
- Genes detected: 16069.
- Zero-count percentage: 5.488%.
- Library-size flag: FALSE.
- Genes-detected flag: TRUE.
- Percent-zero flag: TRUE.

This is a descriptive QC flag only. The sample is retained pending further assessment.

Genes detected and zero-count percentage are mathematically related because both use the same count > 0 detection criterion. They should therefore not be interpreted as independent QC signals.

## Clinical missingness

- Age at treatment: 0/36 missing (0%).
- Gender: 0/36 missing (0%).
- BRAF mutation: 6/36 missing (16.7%).
- LDH high: 5/36 missing (13.9%).
- Biopsy site: 0/36 missing (0%).
- Biopsy type: 0/36 missing (0%).

Missing BRAF and LDH values are retained as source-data missingness and are not treated as pipeline failures.

## Clinical category consistency

### Biopsy type

- visceral / Visceral

### Biopsy site

- small bowel / Small bowel

Raw clinical categories are preserved at this stage.

Capitalization, punctuation or abbreviation differences are documented rather than silently harmonized.

Any future harmonization must be explicit, reproducible and traceable.

## Interpretation boundary

All Day 4 analyses are descriptive quality-control analyses.

No association testing, predictive modelling or causal interpretation is performed at this stage.
