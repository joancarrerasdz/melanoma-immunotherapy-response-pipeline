# P1 preprocessing contract

## Purpose

This document formalizes the preprocessing contract used by the
melanoma-immunotherapy-response-pipeline methodological rebuild.

It addresses portfolio QA criteria U08 (preprocessing) and P1-04
(transcriptomic preprocessing).

This document does not introduce a new analytical procedure and does not
modify any frozen Week 3 or Week 4 result.

## 1. Validated modelling foundation

The predictive rebuild starts from the validated Week 1 inputs:

- data_processed/GSE160638_raw_counts_PD1_aligned.rds
- data_processed/GSE160638_TableS2A_PD1_clinical.csv

The validated internal cohort contains:

- 36 anti-PD-1 samples;
- 22 Responders;
- 14 NonResponders;
- 17,002 genes in the validated raw-count input.

Clinical/sample alignment is established before predictive modelling.

No exploratory Week 2 result is used to redefine the cohort or endpoint.

## 2. Week 2 exploratory preprocessing

Week 2 provides an exploratory expression-analysis layer.

Its outputs include descriptive or exploratory analyses such as:

- expression filtering summaries;
- TMM/logCPM-based exploratory expression matrices;
- PCA;
- sample correlations;
- exploratory differential-expression analysis.

These analyses are useful for data understanding and quality assessment.

They are not globally fitted predictive preprocessing inputs and are not
used as a globally preprocessed feature matrix for Week 3 modelling.

Week 2 exploratory transformations do not replace the
training-partition-specific preprocessing required during predictive
validation.

## 3. Week 3 predictive preprocessing

Week 3 implements predictive modelling using strict nested
cross-validation.

The central leakage-control rule is that cohort-level or supervised
operations relevant to prediction are estimated from training data only.

Accordingly:

- outer-test samples are excluded from supervised preprocessing estimation;
- expression filtering is estimated within training partitions;
- TMM normalization used for differential-expression feature ranking is
  estimated within the corresponding training partition;
- predictive expression transformation uses sample-wise logCPM calculated
  independently for each sample and does not estimate parameters from
  held-out samples;
- feature selection is performed inside the resampling procedure;
- hyperparameter tuning is restricted to training data;
- the held-out outer fold is used only for out-of-fold evaluation.

The strict Week 3 workflow uses:

- 5 outer folds;
- 3 inner folds;
- 100 selected genes per training partition.

Validation performance is calculated exclusively from pooled outer
out-of-fold predictions.

Gene-frequency analysis after nested cross-validation is a
post-validation stability analysis and is not independent validation.

## 4. Week 3 stability-derived candidate set

Across the five Week 3 outer training partitions:

- 312 genes were selected in at least 1/5 folds;
- 100 genes were selected in at least 2/5 folds;
- 51 genes were selected in at least 3/5 folds;
- 25 genes were selected in at least 4/5 folds;
- 12 genes were selected in all 5/5 folds.

The 12 genes selected in all five outer folds are treated as
stability-derived consensus candidates.

They are not described as an independently or clinically validated
molecular signature.

The historical 13-gene result is retained only as a reproducibility
reference and does not define the strict Week 3 candidate set.

## 5. Week 4 external-expression representation

External cohorts use different expression scales or platform-derived
representations and therefore are not transferred directly onto the
internal numerical expression scale.

Before opening external response labels, Week 4 audits candidate
cross-dataset representations using expression information only.

The selected representation is samplewise_rank.

Ranking is calculated independently within each sample across the frozen
12-gene core.

This representation:

- uses the same frozen 12-gene space for all cohorts;
- does not require external response labels;
- does not fit a transformation using pooled external outcomes;
- does not use external performance to choose the representation.

The representation decision is documented in:

- docs/week4_expression_representation_decision.md

## 6. Frozen external-validation contract

Before external outcome labels were opened, Week 4 froze:

- the 12-gene candidate space;
- the samplewise_rank representation;
- the Random Forest model;
- mtry = 1;
- 500 trees;
- classification threshold = 0.48;
- 76 blind external predictions.

The external cohorts are:

- GSE91061: 49 pre-treatment samples;
- GSE78220: 27 pre-treatment samples.

After unblinding, no external cohort is permitted to drive:

- gene reselection;
- representation changes;
- model refitting;
- hyperparameter retuning;
- classification-threshold changes;
- threshold optimization using external outcomes.

This separation is documented in:

- docs/week4_locked_external_validation.md

## 7. End-to-end preprocessing boundary

The methodological sequence is:

    Validated raw counts + clinical endpoint
                    |
                    v
    Week 2 exploratory preprocessing
    (descriptive / exploratory only)
                    |
                    v
    Week 3 strict nested CV
                    |
                    +--> filtering estimated from training data
                    |
                    +--> TMM-based DE ranking fitted in training data
                    |
                    +--> predictive sample-wise logCPM transformation
                    |    with no parameters estimated from held-out samples
                    |
                    +--> feature selection inside resampling
                    |
                    +--> hyperparameter tuning inside training data
                    |
                    v
    Outer out-of-fold validation
                    |
                    v
    Post-validation gene-stability analysis
                    |
                    v
    Frozen 12-gene consensus candidate set
                    |
                    v
    Week 4 label-free cross-dataset representation audit
                    |
                    v
    Frozen samplewise_rank representation
                    |
                    v
    Frozen deployment model + blind predictions
                    |
                    v
    External outcome unblinding and locked evaluation

## 8. Leakage safeguards

The following operations are explicitly prohibited:

- global supervised feature selection before cross-validation;
- use of outer-test samples to estimate predictive preprocessing;
- use of external response labels to select the Week 4 representation;
- model refitting after external-label unblinding;
- gene reselection after external-label unblinding;
- hyperparameter tuning using external outcomes;
- classification-threshold optimization using external outcomes.

These safeguards separate:

- exploratory analysis;
- internal model development;
- internal validation;
- stability analysis;
- external transportability assessment.

## 9. Reproducibility

The principal reproduction command for Week 3 is:

    Rscript --vanilla scripts/run_week3.R

The principal reproduction command for Week 4 is:

    Rscript --vanilla scripts/run_week4.R

The Week 4 runner verifies the frozen deployment artifacts and does not
rebuild the frozen model after external outcomes have been opened.

## 10. QA interpretation

This contract provides formal methodological evidence relevant to:

- U08 — preprocessing;
- P1-04 — transcriptomic preprocessing.

The underlying analytical results are unchanged.

The purpose of this remediation is to document the already implemented
preprocessing boundaries, ordering, leakage safeguards and reproducibility
contract without retrospectively modifying the frozen Week 3 or Week 4
analysis.
