# Historical thesis results — reference only

## Status

This document preserves selected numerical outputs from the original master's
thesis workflow for provenance and comparison.

These results are **historical reference outputs**. They are not the primary
results of the reconstructed strict pipeline.

The original academic repository state is preserved at tag:

`tfm-original`

The current methodological conclusions are reported in the repository
`README.md`, together with the strict Week 3 nested-cross-validation and
Week 4 locked external-validation results.

---

## Historical internal validation — final 100-gene model

| Metric | Value |
|---|---:|
| Accuracy | 0.753425 |
| Kappa | 0.506006 |
| Sensitivity | 0.810811 |
| Specificity | 0.694444 |
| Balanced Accuracy | 0.752628 |
| AUC | 0.827703 |
| Mean best mtry | 5.200000 |
| Unique genes selected | 292 |
| Genes selected in >=4 folds | 29 |
| Genes selected in 5/5 folds | 13 |

Historical source:

`results/rf_cv_metrics_summary_final_100_no_leakage.csv`

---

## Historical signature-size comparison

| Number of genes | AUC | Accuracy | Balanced Accuracy |
|---:|---:|---:|---:|
| 50 | 0.800300 | 0.712329 | 0.712087 |
| 100 | 0.832958 | 0.767123 | 0.766517 |
| 250 | 0.831456 | 0.780822 | 0.780405 |
| 500 | 0.827703 | 0.780822 | 0.780030 |

Historical source:

`results_reference/signature_size_comparison_metrics.csv`

---

## Historical external validation — final 100-gene model

| Dataset | AUC | Accuracy | Sensitivity | Specificity | Balanced Accuracy |
|---|---:|---:|---:|---:|---:|
| GSE91061 | 0.670513 | 0.612245 | 0.600000 | 0.615385 | 0.607692 |
| GSE78220 | 0.611111 | 0.444444 | 0.000000 | 1.000000 | 0.500000 |

Historical source:

`results/external_validation_metrics_by_dataset_final_100.csv`

---

## Historical gene-stability distribution

| Folds selected | Number of genes |
|---:|---:|
| 1 | 180 |
| 2 | 58 |
| 3 | 25 |
| 4 | 16 |
| 5 | 13 |

Historical source:

`results/gene_stability_summary_final_100.csv`

---

## Relationship to the strict rebuild

The reconstructed pipeline identified 12 genes selected in all five strict
Week 3 outer folds.

The historical workflow reported 13 genes selected in all five folds.

The strict reproducibility audit found no overlap between those two 5/5 sets.

This observation documents lack of reproduction of the historical
feature-selection stability under the stricter nested-cross-validation
procedure. It does not establish biological irrelevance of the historical
genes.

See:

`docs/week3_historical_vs_strict_stability_audit.md`
