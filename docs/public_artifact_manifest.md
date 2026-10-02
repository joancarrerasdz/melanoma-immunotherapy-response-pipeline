# Public artifact manifest

## Purpose

This document defines which repository artifacts represent the current strict
methodological rebuild and which artifacts are retained only for exploratory,
audit or historical provenance purposes.

The classification does not alter, recompute or remove any analytical result.

---

## 1. Canonical public figures

The following figures represent the primary current results of the strict
rebuild.

### Week 3 — internal nested cross-validation

- `figures/week3_nested_cv_roc_strict.png`

  Primary visual summary of strict internal nested-cross-validation
  discrimination.

- `figures/week3_gene_selection_stability_strict.png`

  Primary visual summary of feature-selection stability across the five
  outer training folds.

### Week 4 — locked external validation

- `figures/week4_external_validation_roc_strict.png`

  Primary visual summary of the locked external evaluation.

These figures may be displayed directly in the public README.

---

## 2. Supporting strict analytical artifacts

Files carrying the `_strict` suffix are part of the reconstructed analytical
and audit trail.

Examples include:

- Week 3 outer-fold validation outputs;
- Week 3 gene-selection stability outputs;
- Week 4 expression-representation audits;
- frozen deployment-model specifications;
- blind external predictions;
- locked external-validation metrics and predictions;
- freeze manifests and provenance records.

These files support reproducibility and methodological auditability but do not
all need to be displayed directly in the public README.

---

## 3. QA artifacts

Files prefixed with `qa_` document formal QA gates.

The principal Week 3–4 remediation gates are:

- preprocessing contract;
- reproducible environment;
- confounding / cohort-heterogeneity assessment.

Their consolidated state is:

- 74 PASS;
- 0 FAIL.

QA artifacts are evidence of methodological verification rather than primary
scientific results.

---

## 4. Exploratory artifacts

Exploratory outputs include Week 1–2 QC and descriptive analyses such as:

- `figures/week2_pca_exploratory.png`;
- `figures/week2_sample_correlation_exploratory.png`;
- Week 2 differential-expression exploratory outputs;
- early cohort-QC figures such as `day4_*`.

These artifacts remain valid descriptive records but are not estimates of
predictive performance and are not part of the primary Week 3–4 result set.

---

## 5. Historical thesis artifacts

Artifacts associated with the original master's-thesis workflow are retained
for provenance and comparison.

This includes files whose names contain patterns such as:

- `final_100`;
- `signature_size`;
- historical 100-gene external-validation outputs;
- historical Random Forest importance outputs.

Examples include:

- `figures/ROC_GSE160638_random_forest_final_100_no_leakage.png`;
- `figures/ROC_GSE91061_external_validation_final_100.png`;
- `figures/ROC_GSE78220_external_validation_final_100.png`;
- `figures/ROC_external_validation_combined_final_100.png`;
- `figures/external_validation_auc_summary_with_thresholds_final_100.png`;
- `figures/external_validation_balanced_accuracy_summary_with_thresholds_final_100.png`;
- `figures/Top20_gene_importance_RF_final_100_exploratory.png`;
- figures beginning with `signature_size_`.

These files must not be interpreted as the primary results of the current
strict reconstructed pipeline.

Historical numerical results are summarized separately in:

- `docs/historical_thesis_results.md`

The original academic repository state is preserved at:

- `tfm-original`

The original thesis document:

- `docs/Memoria.pdf`

is also retained as historical academic provenance.

---

## 6. Public interpretation rule

When an artifact from the strict rebuild and an artifact from the historical
workflow address the same analytical question, the strict rebuild artifact is
the canonical current result.

Historical artifacts remain available for transparency, provenance and
methodological comparison.

No historical output may replace or be combined with a strict result in a way
that obscures the distinction between the two analytical workflows.

---

## 7. Current public result set

The minimal public-facing result set is therefore:

1. strict Week 3 nested-cross-validation performance;
2. strict Week 3 feature-selection stability;
3. strict Week 4 cohort-specific external validation;
4. methodological QA status;
5. reproducibility and provenance documentation.

This classification is presentation-only and does not modify the frozen Week 3
or Week 4 analyses.
