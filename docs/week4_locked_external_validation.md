# Week 4 — Locked external validation

## Objective

Week 4 evaluates whether the predictive signal reconstructed during Week 3
transfers to independent melanoma anti-PD-1 cohorts.

The external-validation procedure is explicitly locked before external
response labels are consulted.

The two external cohorts are:

- GSE91061
- GSE78220

## Frozen predictive specification

The deployment model is based on the 12 genes selected in all 5/5 strict
Week 3 outer training folds:

- DENND2A
- FCRL1
- FCRL2
- MS4A1
- PLCG2
- RASL10B
- SCARF2
- SDR42E1
- TICRR
- TLR1
- TRIM71
- VARS2

The frozen deployment specification is:

- training dataset: GSE160638;
- training samples: 36;
- genes: 12;
- expression representation: `samplewise_rank`;
- ranking scope: within each sample across the frozen 12-gene core;
- Random Forest `mtry`: 1;
- Random Forest trees: 500;
- classification threshold: 0.48.

The threshold was selected from internal out-of-fold predictions before
external response labels were opened.

## External-input contract

### GSE91061

- 49 retained pre-treatment samples;
- 39 NonResponders;
- 10 Responders;
- all 49 samples matched to expression data;
- all 12 frozen genes available;
- no duplicated sample identifiers.

### GSE78220

- 27 retained pre-treatment samples;
- 12 NonResponders;
- 15 Responders;
- all 27 samples matched to expression data;
- all 12 frozen genes available;
- no duplicated sample identifiers.

Across both cohorts, 76 external samples are evaluated.

## Pre-unblinding freeze

The deployment model and all 76 external predictions were frozen before
external response labels were opened.

The blind-prediction freeze is preserved at commit:

`3e427be1b885b9514cb8acc857e3eb3ca90e4ceb`

A manifest stores MD5 checksums for the frozen deployment artifacts.

The final Week 4 reproduction gate verifies that these checksums remain
unchanged.

No model component is regenerated or modified after unblinding.

## Locked external-validation results

Dataset-specific results are the primary external-validation results.

| Dataset | Samples | AUC | Accuracy | Sensitivity | Specificity | Balanced Accuracy |
|---|---:|---:|---:|---:|---:|---:|
| GSE91061 | 49 | 0.6641 | 0.3265 | 0.9000 | 0.1795 | 0.5397 |
| GSE78220 | 27 | 0.6444 | 0.6667 | 0.8000 | 0.5000 | 0.6500 |

The pooled result is reported only as a secondary summary:

| Scope | Samples | AUC | Accuracy | Sensitivity | Specificity | Balanced Accuracy |
|---|---:|---:|---:|---:|---:|---:|
| Pooled | 76 | 0.5718 | 0.4474 | 0.8400 | 0.2549 | 0.5475 |

The pooled estimate is secondary because the two external cohorts differ in
response prevalence and cohort composition.

## Interpretation

The locked external evaluation shows limited and cohort-dependent
transportability.

GSE91061 retains high sensitivity but has low specificity, indicating a
strong tendency to predict the Responder class in that cohort.

GSE78220 shows a more balanced classification profile, although
discrimination remains moderate.

These results are reported exactly as observed.

No external cohort was used to:

- reselect genes;
- change the expression representation;
- retune `mtry`;
- refit the deployment model;
- modify the classification threshold;
- optimize the threshold using external outcomes.

The external results therefore assess the frozen Week 4 deployment procedure
rather than a model adapted to the external cohorts.

## Methodological boundary

The Week 3 AUC of 0.8896 is an internal nested-cross-validation estimate and
must not be presented as external-validation performance.

The Week 4 dataset-specific external results are the primary assessment of
transportability.

The pooled external result is a secondary descriptive summary.

The 12 genes are stability-derived consensus candidates. Week 4 evaluates a
deployment model built from those candidates but does not establish a
clinically validated molecular signature.

## Reproduction

From the repository root:

```bash
Rscript --vanilla scripts/run_week4.R
```

The Week 4 runner:

1. reconstructs and validates the strict external-input contract;
2. reproduces the expression-representation audit;
3. verifies the frozen deployment artifacts by MD5;
4. verifies the 76 blind external predictions;
5. reproduces the locked external evaluation;
6. validates the expected external metrics and methodological safeguards.

The runner deliberately does not execute
`scripts/16_freeze_deployment_model_strict.R`.

The model and blind predictions must remain those frozen before external
response labels were opened.

A successful reproduction terminates with:

```text
WEEK 4 REPRODUCTION: PASS
```
