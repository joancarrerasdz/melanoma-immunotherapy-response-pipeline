# Week 3 — Historical nested-CV audit

## Scope

This document records the methodological audit of the historical
random-forest cross-validation workflow before the Week 3 rebuild.

The objective is to distinguish components that are already suitable
for reuse from components that require reconstruction to guarantee
strict separation between model development and evaluation data.

## Historical outer cross-validation

The historical workflow uses:

- 5 outer cross-validation folds;
- explicit outer training and test partitions;
- feature selection restricted to the corresponding outer-training set;
- model fitting on the outer-training set;
- prediction on the untouched outer-test samples;
- storage of out-of-fold class predictions and Responder probabilities;
- aggregation of all outer-fold predictions for the final confusion
  matrix and ROC/AUC calculation.

No evidence was identified that outer-test response labels are used for
feature selection, hyperparameter tuning, or model fitting.

The out-of-fold prediction aggregation strategy is therefore retained
as a valid design principle for the rebuild.

## Inner cross-validation limitation

Within each outer fold, the historical workflow performs feature
selection before the inner cross-validation procedure.

Conceptually, the sequence is:

1. define the complete outer-training set;
2. perform filtering and differential-expression-based feature
   selection using the complete outer-training set;
3. generate the outer-training expression matrix;
4. perform inner cross-validation only for random-forest tuning.

Consequently, samples subsequently used as inner-validation samples
have already contributed to preprocessing and feature selection.

This does not expose the outer-test labels to the model, but it means
that the inner cross-validation layer is not fully isolated.

The historical implementation therefore does not satisfy the strict
nested-CV standard adopted for the rebuilt pipeline.

## Normalization limitation

The historical workflow calculates normalization factors separately for
the complete outer-training matrix and the complete outer-test matrix.

Although this operation does not use outer-test response labels, the
transformation of an outer-test sample may depend on the other samples
present in the same test fold.

The Week 3 rebuild will avoid test-fold-dependent fitting and will use a
strictly training-derived preprocessing strategy for predictive
evaluation.

## Full-dataset model

The historical script also fits a model to the complete dataset after
cross-validation and calculates variable importance.

This step is acceptable only as a post-validation exploratory analysis.

It must not contribute to:

- cross-validation predictions;
- reported validation performance;
- hyperparameter selection;
- claims of predictive generalization.

## Week 3 rebuild decision

The predictive workflow will be reconstructed so that every
data-dependent preprocessing and feature-selection operation used for
model development is fitted within the appropriate training partition.

For every outer fold:

1. the outer-test samples remain isolated;
2. the outer-training set is divided into inner folds;
3. preprocessing is fitted using inner-training samples only;
4. feature selection is performed using inner-training samples only;
5. candidate models are evaluated on inner-validation samples;
6. hyperparameters are selected using inner-validation performance;
7. the complete predictive pipeline is refitted using the complete
   outer-training set;
8. the fitted pipeline is applied once to the outer-test samples;
9. outer-test predictions are stored as out-of-fold predictions.

After all outer folds have completed, global performance is calculated
exclusively from the pooled out-of-fold predictions.

Any model subsequently fitted to the complete dataset is considered
exploratory and is kept separate from validation metrics.

## Audit conclusion

Historical outer-test separation: PASS

Historical out-of-fold aggregation: PASS

Strict inner-CV preprocessing isolation: REBUILD REQUIRED

Strict inner-CV feature-selection isolation: REBUILD REQUIRED

Strict test-independent preprocessing: REBUILD REQUIRED

Week 3 implementation status: READY FOR RECONSTRUCTION
